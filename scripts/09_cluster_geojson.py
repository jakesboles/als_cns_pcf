import os
import sys
import json
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely

from pcf_utils import log, find_files, load_cell_geoms

sys.stdout.reconfigure(line_buffering=True)

COLORS = {"astrocytes": [40, 170, 40], "cd4_tcell": [0, 190, 255], "cd8_tcell": [30, 60, 220],
          "cortical_neurons": [150, 60, 200], "myeloid": [220, 30, 30], "oligodendrocyte": [255, 170, 0]}
OTHER = [120, 120, 120]
MARKERS = ["NeuN", "GFAP", "Vimentin", "Iba-1", "Olig-2", "CD3e", "CD4", "CD8"]


def write_code(code, args, cells, cell_files):
  q = pd.read_parquet(os.path.join(args.quant_dir, f"{code}.parquet"), columns=["cell_id"] + args.markers)
  cells = cells.merge(q, on="cell_id", how="left")
  cells = cells.merge(load_cell_geoms(cell_files.get(code, {}), set(cells["cell_id"]), code),
                      on="cell_id", how="left")
  missing = int(cells["geometry"].isna().sum())
  if missing:
    log(f"WARNING {code}: {missing} clustered cells have no polygon; skipped")
  cells = cells[cells["geometry"].notna()]

  out_dir = os.path.join(args.out_dir, code)
  os.makedirs(out_dir, exist_ok=True)
  meas_cols = ["area_um2"] + args.markers
  for k, sub in cells.groupby("cluster", sort=True, observed=True):
    if args.clusters and k not in args.clusters:
      continue
    geoms = [json.loads(g) for g in shapely.to_geojson(sub["geometry"].to_numpy())]
    meas = sub[meas_cols].astype(float).to_numpy()
    feats = [{"type": "Feature", "geometry": g,
              "properties": {"objectType": "detection", "name": cid,
                             "classification": {"name": ct, "color": COLORS.get(ct, OTHER)},
                             "measurements": {c: (None if np.isnan(v) else v) for c, v in zip(meas_cols, m)}}}
             for g, cid, ct, m in zip(geoms, sub["cell_id"], sub["cell_type"], meas)]
    with open(os.path.join(out_dir, f"{code}_cluster{k}.geojson"), "w") as fh:
      json.dump({"type": "FeatureCollection", "features": feats}, fh)

  counts = pd.crosstab(cells["cell_type"], cells["cluster"], margins=True)
  counts.to_csv(os.path.join(out_dir, f"{code}_cluster_counts.csv"))
  return code, len(cells)


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--resolution", default="0.3")
  p.add_argument("--clusters", nargs="*")
  p.add_argument("--codes", nargs="*")
  p.add_argument("--markers", nargs="*", default=MARKERS)
  p.add_argument("--clusters_file", default="data/08_loadings_clustering/clusters.parquet")
  p.add_argument("--flags", default="data/05_object_creation/cell_flags.parquet")
  p.add_argument("--quant_dir", default="data/01_quantify_cells/")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--results_dir", default="results/09_cluster_geojson/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  col = f"leiden_{args.resolution}"
  cl = pd.read_parquet(args.clusters_file)
  if col not in cl.columns:
    sys.exit(f"{col} not in {args.clusters_file}; available: {', '.join(c for c in cl.columns if c != 'cell_id')}")
  cl = cl[["cell_id", col]].rename(columns={col: "cluster"})
  cl["cluster"] = cl["cluster"].astype(str)
  cl["cluster"] = pd.Categorical(cl["cluster"], categories=sorted(cl["cluster"].unique(), key=lambda s: (len(s), s)))

  flags = pd.read_parquet(args.flags, columns=["code", "cell_id", "cell_type", "area_um2"])
  cells = cl.merge(flags, on="cell_id", how="left")
  if cells["code"].isna().any():
    sys.exit(f"{int(cells['code'].isna().sum())} clustered cells not in {args.flags}; rerun 05/07/08 consistently")

  if args.clusters:
    unknown = set(args.clusters) - set(cells["cluster"].astype(str))
    if unknown:
      sys.exit(f"Clusters not in {col}: {', '.join(sorted(unknown))}")
    args.clusters = set(args.clusters)
  codes = args.codes or sorted(cells["code"].unique())
  unknown = set(codes) - set(cells["code"])
  if unknown:
    sys.exit(f"Codes with no clustered cells: {', '.join(sorted(unknown))}")

  first = pq.read_schema(os.path.join(args.quant_dir, f"{codes[0]}.parquet")).names
  absent = [m for m in args.markers if m not in first]
  if absent:
    sys.exit(f"Markers not in 01 output: {', '.join(absent)}")

  cell_files = defaultdict(dict)
  for cell_type, by_code in find_files(args.cells_dir).items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths

  args.out_dir = os.path.join(args.results_dir, col)
  os.makedirs(args.out_dir, exist_ok=True)
  pd.crosstab(cells["cell_type"], cells["cluster"], margins=True).to_csv(
    os.path.join(args.out_dir, "cluster_counts_all.csv"))
  log(f"{col}: {cells['cluster'].nunique()} clusters, {len(codes)} images")

  jobs = [(c, args, cells[cells["code"] == c], dict(cell_files)) for c in codes]
  with mp.get_context("fork").Pool(min(args.workers, len(codes))) as pool:
    for code, n in pool.starmap(write_code, jobs):
      log(f"{code}: {n} cells written")
  log("Done")


if __name__ == "__main__":
  main()
