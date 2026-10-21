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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pcf_utils import log, get_code, to_polygon

sys.stdout.reconfigure(line_buffering=True)

PLOT_MAX_CELLS = 200000
COLORS = {"unassigned": "lightgray", "multiple": "red"}


def find_region_files(region_dir):
  files = defaultdict(lambda: defaultdict(list))
  for f in sorted(glob.glob(os.path.join(region_dir, "*", "*.geojson"))):
    code = get_code(f)
    if code is None:
      log(f"WARNING: can't parse code from {f}, skipping")
      continue
    label = os.path.basename(os.path.dirname(f)).removesuffix("_annotations").lower()
    files[label][code].append(f)
  return {label: dict(by_code) for label, by_code in files.items()}


def load_region(paths):
  geoms, n_drop = [], 0
  for p in paths:
    for g in gpd.read_file(p).geometry:
      poly = to_polygon(g)
      if poly is None:
        n_drop += 1
      else:
        geoms.append(poly)
  return (shapely.union_all(geoms) if geoms else None), len(geoms), n_drop


def annotate(code, in_dir, data_dir, results_dir, region_files):
  cells = pd.read_parquet(os.path.join(in_dir, f"{code}.parquet"), columns=["cell_id", "cell_type", "x", "y"])
  x, y = cells["x"].to_numpy(), cells["y"].to_numpy()
  out = pd.DataFrame({"code": code, "cell_id": cells["cell_id"]})
  regions = {}

  for label in sorted(region_files):
    col = f"in_{label}"
    paths = region_files[label].get(code)
    if not paths:
      log(f"{code}: no {label} annotation, {col} = NA")
      out[col] = pd.array([pd.NA] * len(cells), dtype="boolean")
      continue
    if len(paths) > 1:
      log(f"{code}: {len(paths)} {label} files, using their union")
    geom, n_poly, n_drop = load_region(paths)
    log(f"{code}: {label} {n_poly} polygons ({n_drop} dropped)")
    if geom is None:
      out[col] = pd.array([False] * len(cells), dtype="boolean")
      continue
    shapely.prepare(geom)
    out[col] = pd.array(shapely.contains_xy(geom, x, y), dtype="boolean")
    regions[label] = geom

  labels = sorted(region_files)
  hits = np.column_stack([out[f"in_{l}"].fillna(False).to_numpy(bool) for l in labels])
  n_hits = hits.sum(axis=1)
  first = np.array([l.upper() for l in labels])[hits.argmax(axis=1)]
  out["region"] = np.where(n_hits == 1, first, np.where(n_hits > 1, "multiple", "unassigned"))

  out.to_parquet(os.path.join(data_dir, f"{code}.parquet"), index=False)
  pd.crosstab(cells["cell_type"], out["region"], margins=True).to_csv(
    os.path.join(results_dir, f"{code}_region_counts.csv"))
  plot_regions(code, regions, x, y, out["region"].to_numpy(), os.path.join(results_dir, f"{code}_regions.png"))
  return code, out["region"].value_counts().to_dict()


def plot_regions(code, regions, x, y, region, out_png):
  idx = np.arange(len(x))
  if len(idx) > PLOT_MAX_CELLS:
    idx = np.random.default_rng(1).choice(idx, PLOT_MAX_CELLS, replace=False)
  cmap = plt.get_cmap("tab10")
  colors = {lab.upper(): cmap(k) for k, lab in enumerate(sorted(regions))} | COLORS
  fig, ax = plt.subplots(figsize=(10, 10))
  for label, geom in regions.items():
    gpd.GeoSeries([geom]).plot(ax=ax, color=colors[label.upper()], alpha=0.15, edgecolor="black", linewidth=0.3)
  for r in sorted(set(region)):
    sel = idx[region[idx] == r]
    ax.scatter(x[sel], y[sel], s=0.2, color=colors.get(r, "black"), label=f"{r} ({(region == r).sum()})",
               rasterized=True)
  ax.invert_yaxis()
  ax.set_aspect("equal")
  ax.axis("off")
  ax.legend(markerscale=20, loc="upper right")
  ax.set_title(code)
  fig.savefig(out_png, dpi=150, bbox_inches="tight")
  plt.close(fig)


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--in_dir", default="data/01_quantify_cells/")
  p.add_argument("--region_dir", default="annotations/region")
  p.add_argument("--data_dir", default="data/02_annotate_regions/")
  p.add_argument("--results_dir", default="results/02_annotate_regions/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  os.makedirs(args.data_dir, exist_ok=True)
  os.makedirs(args.results_dir, exist_ok=True)

  codes = args.codes or sorted(os.path.basename(f).removesuffix(".parquet")
                               for f in glob.glob(os.path.join(args.in_dir, "*.parquet")))
  if not codes:
    sys.exit(f"No 01 outputs in {args.in_dir}")
  region_files = find_region_files(args.region_dir)
  if not region_files:
    sys.exit(f"No region geoJSONs under {args.region_dir}/<region>_annotations/")
  log(f"{len(codes)} images; regions: {', '.join(sorted(region_files))}")

  jobs = [(c, args.in_dir, args.data_dir, args.results_dir, region_files) for c in codes]
  with mp.get_context("fork").Pool(min(args.workers, len(codes))) as pool:
    for code, counts in pool.starmap(annotate, jobs):
      log(f"{code}: {counts}")
  log("Done")


if __name__ == "__main__":
  main()
