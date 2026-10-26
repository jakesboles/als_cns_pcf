import os
import sys
import glob
import json
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

from pcf_utils import log, get_code, to_polygon

sys.stdout.reconfigure(line_buffering=True)

QC_COLORS = {"pos": [220, 30, 30], "touch": [30, 120, 220], "neg": [120, 120, 120], "feature": [255, 200, 0]}


def find_files(root):
  files = defaultdict(lambda: defaultdict(list))
  for f in sorted(glob.glob(os.path.join(root, "*", "*.geojson"))):
    code = get_code(f)
    if code is None:
      log(f"WARNING: can't parse code from {f}, skipping")
      continue
    files[os.path.basename(os.path.dirname(f))][code].append(f)
  return {k: dict(v) for k, v in files.items()}


def load_cells(cell_files, keep_ids, code):
  ids, geoms = [], []
  for cell_type, paths in sorted(cell_files.items()):
    for p in paths:
      for k, g in enumerate(gpd.read_file(p).geometry, start=1):
        ids.append(f"{code}_{cell_type}_{k}")
        geoms.append(g)
  df = pd.DataFrame({"cell_id": ids, "geometry": geoms})
  df = df[df["cell_id"].isin(keep_ids)]
  df["geometry"] = [to_polygon(g) for g in df["geometry"]]
  return df


def load_features(paths):
  geoms = [to_polygon(g) for p in paths for g in gpd.read_file(p).geometry]
  kept = np.array([g for g in geoms if g is not None], dtype=object)
  return kept, len(geoms) - len(kept)


def feature(geom, props):
  return {"type": "Feature", "geometry": json.loads(shapely.to_geojson(geom)), "properties": props}


def write_qc_geojson(path, label, min_frac, cell_ids, cell_geoms, pos, touch, context, feats):
  out = []
  for k in np.flatnonzero(pos | touch | context):
    kind = "pos" if pos[k] else "touch" if touch[k] else "neg"
    cls = {"pos": f"{label}+", "touch": f"{label} touch <{min_frac:g}", "neg": f"{label}-"}[kind]
    out.append(feature(cell_geoms[k], {
      "objectType": "detection", "name": cell_ids[k],
      "classification": {"name": cls, "color": QC_COLORS[kind]}}))
  for g in feats:
    out.append(feature(g, {
      "objectType": "annotation", "isLocked": True,
      "classification": {"name": f"{label} feature", "color": QC_COLORS["feature"]}}))
  with open(path, "w") as fh:
    json.dump({"type": "FeatureCollection", "features": out}, fh)


def annotate(code, args, cell_files, path_files):
  base = pd.read_parquet(os.path.join(args.in_dir, f"{code}.parquet"), columns=["cell_id", "cell_type"])
  cells = load_cells(cell_files.get(code, {}), set(base["cell_id"]), code)
  cells = base.merge(cells, on="cell_id", how="left")
  missing = cells["geometry"].isna().sum()
  if missing:
    log(f"WARNING {code}: {missing} cells in 01 output have no polygon in annotations/cells; their calls will be NA")
  has_geom = cells["geometry"].notna().to_numpy()
  geoms = cells["geometry"].to_numpy()
  shapely.prepare(geoms[has_geom])

  out = pd.DataFrame({"code": code, "cell_id": cells["cell_id"]})
  summary = {}
  for label in sorted(path_files):
    paths = path_files[label].get(code)
    if not paths:
      log(f"{code}: no {label} annotation, has_{label} = NA")
      for col, dt in [(f"has_{label}", "boolean"), (f"n_{label}", "Int32"),
                      (f"{label}_overlap_px", "Float64"), (f"{label}_max_frac", "Float64")]:
        out[col] = pd.array([pd.NA] * len(out), dtype=dt)
      continue
    feats, n_drop = load_features(paths)
    log(f"{code}: {label} {len(feats)} features ({n_drop} dropped)")

    n_hit = np.zeros(len(out), dtype=np.int32)
    ov_area = np.zeros(len(out))
    max_frac = np.zeros(len(out))
    context = np.zeros(len(out), dtype=bool)
    touch = np.zeros(len(out), dtype=bool)
    if len(feats):
      tree = shapely.STRtree(feats)
      idx = np.flatnonzero(has_geom)
      ci, fi = tree.query(geoms[idx], predicate="intersects")
      ci = idx[ci]
      inter = shapely.area(shapely.intersection(geoms[ci], feats[fi]))
      frac = inter / shapely.area(feats[fi])
      np.add.at(n_hit, ci, frac >= args.min_frac)
      np.add.at(ov_area, ci, inter)
      np.maximum.at(max_frac, ci, frac)
      touch[ci] = True
      cj, _ = tree.query(geoms[idx], predicate="dwithin", distance=args.context_px)
      context[idx[cj]] = True

    pos = n_hit > 0
    na = ~has_geom
    out[f"has_{label}"] = pd.array(np.where(na, None, pos), dtype="boolean")
    out[f"n_{label}"] = pd.array(np.where(na, None, n_hit), dtype="Int32")
    out[f"{label}_overlap_px"] = pd.array(np.where(na, None, ov_area), dtype="Float64")
    out[f"{label}_max_frac"] = pd.array(np.where(na, None, max_frac), dtype="Float64")
    summary[label] = int(pos.sum())

    pd.crosstab(cells["cell_type"], out[f"has_{label}"].astype(str), margins=True).to_csv(
      os.path.join(args.results_dir, f"{code}_{label}_counts.csv"))
    if not args.no_qc_geojson:
      write_qc_geojson(os.path.join(args.results_dir, f"{code}_{label}_qc.geojson"), label,
                       args.min_frac, cells["cell_id"].to_numpy(), geoms, pos & has_geom,
                       touch & ~pos & has_geom, context & has_geom, feats)

  out.to_parquet(os.path.join(args.data_dir, f"{code}.parquet"), index=False)
  return code, summary


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--in_dir", default="data/01_quantify_cells/")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--pathology_dir", default="annotations/pathology")
  p.add_argument("--data_dir", default="data/03_annotate_pathology/")
  p.add_argument("--results_dir", default="results/03_annotate_pathology/")
  p.add_argument("--min_frac", type=float, default=0.05)
  p.add_argument("--context_px", type=float, default=50)
  p.add_argument("--no_qc_geojson", action="store_true")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  os.makedirs(args.data_dir, exist_ok=True)
  os.makedirs(args.results_dir, exist_ok=True)

  codes = args.codes or sorted(os.path.basename(f).removesuffix(".parquet")
                               for f in glob.glob(os.path.join(args.in_dir, "*.parquet")))
  if not codes:
    sys.exit(f"No 01 outputs in {args.in_dir}")
  path_files = find_files(args.pathology_dir)
  if not path_files:
    sys.exit(f"No pathology geoJSONs under {args.pathology_dir}/<feature>/")
  by_type = find_files(args.cells_dir)
  cell_files = defaultdict(dict)
  for cell_type, by_code in by_type.items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths
  log(f"{len(codes)} images; pathology: {', '.join(sorted(path_files))}; min_frac {args.min_frac}")

  jobs = [(c, args, dict(cell_files), path_files) for c in codes]
  with mp.get_context("fork").Pool(min(args.workers, len(codes))) as pool:
    for code, summary in pool.starmap(annotate, jobs):
      log(f"{code}: positive cells {summary}")
  log("Done")


if __name__ == "__main__":
  main()
