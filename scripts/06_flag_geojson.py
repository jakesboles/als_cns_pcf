import os
import sys
import json
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
import shapely

from pcf_utils import log, find_files, load_cell_geoms

sys.stdout.reconfigure(line_buffering=True)

COLORS = {"kept": [40, 170, 40], "size": [240, 150, 0], "region": [150, 60, 200],
          "lineage": [220, 30, 30], "multiple": [255, 0, 255], "NA": [120, 120, 120]}


def label_of(keep, reason):
  if pd.isna(keep):
    return "NA"
  if keep:
    return "kept"
  return reason if ";" not in reason else f"multiple ({reason})"


def color_of(label):
  return COLORS["multiple"] if label.startswith("multiple") else COLORS.get(label, COLORS["NA"])


def write_geojson(code, args, flags, cell_files):
  cells = flags.merge(load_cell_geoms(cell_files.get(code, {}), set(flags["cell_id"]), code),
                      on="cell_id", how="left")
  missing = cells["geometry"].isna().sum()
  if missing:
    log(f"WARNING {code}: {missing} flagged cells have no polygon; skipped")
  cells = cells[cells["geometry"].notna()]
  if args.dropped_only:
    cells = cells[cells["keep"] != True]

  labels = [label_of(k, r) for k, r in zip(cells["keep"], cells["flag_reason"])]
  geoms = [json.loads(g) for g in shapely.to_geojson(cells["geometry"].to_numpy())]
  feats = [{"type": "Feature", "geometry": g,
            "properties": {"objectType": "detection", "name": cid,
                           "classification": {"name": lab, "color": color_of(lab)},
                           "measurements": {"area_um2": float(a)}}}
           for g, cid, lab, a in zip(geoms, cells["cell_id"], labels, cells["area_um2"])]
  with open(os.path.join(args.results_dir, f"{code}_flags.geojson"), "w") as fh:
    json.dump({"type": "FeatureCollection", "features": feats}, fh)

  counts = pd.crosstab(cells["cell_type"], pd.Series(labels, index=cells.index, name="status"), margins=True)
  counts.to_csv(os.path.join(args.results_dir, f"{code}_flag_counts.csv"))
  return code, len(feats), pd.Series(labels).value_counts().to_dict()


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--flags", default="data/05_object_creation/cell_flags.parquet")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--results_dir", default="results/06_flag_geojson/")
  p.add_argument("--dropped_only", action="store_true")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  os.makedirs(args.results_dir, exist_ok=True)
  if not os.path.exists(args.flags):
    sys.exit(f"{args.flags} not found; run 05 (its last section exports it)")
  flags = pd.read_parquet(args.flags)
  codes = args.codes or sorted(flags["code"].unique())
  unknown = set(codes) - set(flags["code"])
  if unknown:
    sys.exit(f"Codes not in {args.flags}: {', '.join(sorted(unknown))}")

  cell_files = defaultdict(dict)
  for cell_type, by_code in find_files(args.cells_dir).items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths
  log(f"{len(codes)} images")

  jobs = [(c, args, flags[flags["code"] == c], dict(cell_files)) for c in codes]
  with mp.get_context("fork").Pool(min(args.workers, len(codes))) as pool:
    for code, n, summary in pool.starmap(write_geojson, jobs):
      log(f"{code}: {n} cells written {summary}")
  log("Done")


if __name__ == "__main__":
  main()
