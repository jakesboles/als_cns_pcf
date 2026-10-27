import os
import sys
import glob
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

from pcf_utils import log, to_polygon, find_files, load_cell_geoms

sys.stdout.reconfigure(line_buffering=True)


def load_features(paths):
  geoms = [to_polygon(g) for p in paths for g in gpd.read_file(p).geometry]
  kept = np.array([g for g in geoms if g is not None], dtype=object)
  return kept, len(geoms) - len(kept)


def pixel_size(base, override):
  if override:
    return override
  ok = base["n_pixels"] > 0
  px = np.sqrt(base.loc[ok, "area_um2"] / base.loc[ok, "n_pixels"])
  return float(np.median(px)) if len(px) and np.isfinite(px).all() else None


def distances(code, args, cell_files, feature_files):
  base = pd.read_parquet(os.path.join(args.in_dir, f"{code}.parquet"),
                         columns=["cell_id", "cell_type", "n_pixels", "area_um2"])
  px_um = pixel_size(base, args.px_um)
  if px_um is None:
    raise RuntimeError(f"{code}: can't derive pixel size from 01 output; pass --px_um")
  cells = base[["cell_id", "cell_type"]].merge(
    load_cell_geoms(cell_files.get(code, {}), set(base["cell_id"]), code), on="cell_id", how="left")
  has_geom = cells["geometry"].notna().to_numpy()
  if (~has_geom).any():
    log(f"WARNING {code}: {(~has_geom).sum()} cells in 01 output have no polygon; distances NA")
  geoms = cells["geometry"].to_numpy()
  idx = np.flatnonzero(has_geom)

  out = pd.DataFrame({"code": code, "cell_id": cells["cell_id"]})
  summary = {}
  for label in sorted(feature_files):
    col = f"dist_{label}_um"
    paths = feature_files[label].get(code)
    feats, n_drop = load_features(paths) if paths else (np.array([], dtype=object), 0)
    if not paths or not len(feats):
      log(f"{code}: no {label} features, {col} = NA")
      out[col] = pd.array([pd.NA] * len(out), dtype="Float64")
      continue
    log(f"{code}: {label} {len(feats)} features ({n_drop} dropped)")
    d = np.full(len(out), np.nan)
    (qi, _), dist = shapely.STRtree(feats).query_nearest(geoms[idx], return_distance=True, all_matches=False)
    d[idx[qi]] = dist * px_um
    out[col] = pd.array(np.where(np.isnan(d), None, d), dtype="Float64")
    summary[label] = round(float(np.nanmedian(d)), 1)

  out.to_parquet(os.path.join(args.data_dir, f"{code}.parquet"), index=False)
  dist_cols = [c for c in out.columns if c.startswith("dist_")]
  out[dist_cols].astype(float).groupby(cells["cell_type"].values).median().to_csv(
    os.path.join(args.results_dir, f"{code}_median_dist_um.csv"))
  return code, px_um, summary


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--in_dir", default="data/01_quantify_cells/")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--vessel_dir", default="annotations/vasculature")
  p.add_argument("--pathology_dir", default="annotations/pathology")
  p.add_argument("--data_dir", default="data/04_feature_distances/")
  p.add_argument("--results_dir", default="results/04_feature_distances/")
  p.add_argument("--px_um", type=float)
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  os.makedirs(args.data_dir, exist_ok=True)
  os.makedirs(args.results_dir, exist_ok=True)

  codes = args.codes or sorted(os.path.basename(f).removesuffix(".parquet")
                               for f in glob.glob(os.path.join(args.in_dir, "*.parquet")))
  if not codes:
    sys.exit(f"No 01 outputs in {args.in_dir}")

  feature_files = find_files(args.pathology_dir)
  vessels = find_files(args.vessel_dir, by_code_only=True)
  if vessels:
    feature_files["vessel"] = vessels
  if not feature_files:
    sys.exit(f"No feature geoJSONs under {args.vessel_dir} or {args.pathology_dir}/<feature>/")

  cell_files = defaultdict(dict)
  for cell_type, by_code in find_files(args.cells_dir).items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths
  log(f"{len(codes)} images; features: {', '.join(sorted(feature_files))}")

  jobs = [(c, args, dict(cell_files), feature_files) for c in codes]
  with mp.get_context("fork").Pool(min(args.workers, len(codes))) as pool:
    for code, px_um, summary in pool.starmap(distances, jobs):
      log(f"{code}: {px_um:.4f} um/px; median distance (um) {summary}")
  log("Done")


if __name__ == "__main__":
  main()
