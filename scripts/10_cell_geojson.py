import os
import re
import sys
import glob
import json
import argparse
import multiprocessing as mp
from collections import defaultdict

import pandas as pd
import shapely

from pcf_utils import log, find_files, load_cell_geoms

sys.stdout.reconfigure(line_buffering=True)

COLORS = {"astrocyte": [40, 170, 40], "cortical_neuron": [150, 60, 200], "oligodendrocyte": [255, 170, 0],
          "Microglia": [220, 30, 30], "Macrophage": [255, 0, 255], "cd4_tcell": [0, 190, 255],
          "cd8_tcell": [30, 60, 220]}
OTHER = [120, 120, 120]
TCELLS = ["cd4_tcell", "cd8_tcell"]
IMAGE_RE = re.compile(r"(JSB(\d+-\d+)_Scan\d+\.er)\.qptiff$")


def as_bool(s):
  return s.astype(str).str.upper().isin(["TRUE", "T", "1"])


def layers(cells):
  # (folder, subfolder, class suffix, cells)
  out = []
  for ct, d in cells.groupby("cell_type2"):
    out.append(("cell_type", ct, "", d))
    out.append(("ptdp_positive", ct, " pTDP+", d[as_bool(d["has_ptdp"])]))
    out.append(("dpr_positive", ct, " DPR+", d[as_bool(d["has_dpr"])]))
    if ct in TCELLS:
      out.append(("tcell_not_touching_vessel", ct, " off vessel", d[pd.to_numeric(d["dist_vessel_um"], errors="coerce") > 0]))
  return out


def write_code(code, stem, args, cells, cell_files):
  cells = cells.merge(load_cell_geoms(cell_files, set(cells["cell_id"]), code), on="cell_id", how="left")
  missing = cells["geometry"].isna().sum()
  if missing:
    log(f"WARNING {code}: {missing} cells have no polygon; skipped")
  cells = cells[cells["geometry"].notna()]
  counts = []
  for folder, ct, suffix, d in layers(cells):
    counts.append({"code": code, "layer": folder, "cell_type2": ct, "n_cells": len(d)})
    if d.empty:
      continue
    geoms = [json.loads(g) for g in shapely.to_geojson(d["geometry"].to_numpy())]
    cls = {"name": ct + suffix, "color": COLORS.get(ct, OTHER)}
    feats = [{"type": "Feature", "geometry": g,
              "properties": {"objectType": "detection", "name": cid, "classification": cls,
                             "measurements": {"area_um2": float(a), **({} if pd.isna(v) else {"dist_vessel_um": float(v)})}}}
             for g, cid, a, v in zip(geoms, d["cell_id"], d["area_um2"],
                                     pd.to_numeric(d["dist_vessel_um"], errors="coerce"))]
    out_dir = os.path.join(args.results_dir, folder, ct)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, f"{stem}.geojson"), "w") as fh:
      json.dump({"type": "FeatureCollection", "features": feats}, fh)
  return counts


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--cells", default="data/09_qc_dimred/filtered_data.csv")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--image_dir", default="pcf_images")
  p.add_argument("--results_dir", default="results/10_cell_geojson/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  os.makedirs(args.results_dir, exist_ok=True)
  if not os.path.exists(args.cells):
    sys.exit(f"{args.cells} not found; run 09 (its first array task writes it)")
  cells = pd.read_csv(args.cells, usecols=["code", "cell_id", "cell_type2", "area_um2", "has_ptdp", "has_dpr",
                                           "dist_vessel_um"], dtype={"code": str})
  codes = args.codes or sorted(cells["code"].unique())
  unknown = set(codes) - set(cells["code"])
  if unknown:
    sys.exit(f"Codes not in {args.cells}: {', '.join(sorted(unknown))}")

  stems = {}
  for f in glob.glob(os.path.join(args.image_dir, "*.qptiff")):
    m = IMAGE_RE.search(os.path.basename(f))
    if m:
      stems[m.group(2)] = m.group(1)
  no_image = [c for c in codes if c not in stems]
  if no_image:
    sys.exit(f"No JSB<code>_Scan<n>.er.qptiff in {args.image_dir} for: {', '.join(no_image)}")

  cell_files = defaultdict(dict)
  for cell_type, by_code in find_files(args.cells_dir).items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths
  log(f"{len(codes)} images, {len(cells)} cells")

  jobs = [(c, stems[c], args, cells[cells["code"] == c], cell_files.get(c, {})) for c in codes]
  with mp.get_context("fork").Pool(max(1, min(args.workers, len(codes)))) as pool:
    counts = [r for res in pool.starmap(write_code, jobs) for r in res]
  counts = pd.DataFrame(counts)
  counts.to_csv(os.path.join(args.results_dir, "cell_counts.csv"), index=False)
  log(counts.pivot_table(index=["layer", "cell_type2"], values="n_cells", aggfunc="sum").to_string())
  log("Done")


if __name__ == "__main__":
  main()
