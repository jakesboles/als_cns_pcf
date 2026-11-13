import os
import sys
import argparse
import multiprocessing as mp

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

from pcf_utils import log, to_polygon, find_files

sys.stdout.reconfigure(line_buffering=True)


def as_bool(s):
  return s.astype(str).str.upper().isin(["TRUE", "T", "1"])


def load_polys(paths):
  geoms = [to_polygon(g) for p in paths for g in gpd.read_file(p).geometry]
  return np.array([g for g in geoms if g is not None], dtype=object)


def band_labels(radii):
  return ["0"] + [f"{a}-{b}" for a, b in zip([0] + radii[:-1], radii)] + [f">{radii[-1]}"]


def to_band(d_um, radii):
  return pd.cut(d_um, [-np.inf, 0] + radii + [np.inf], labels=band_labels(radii), right=True)


def regions_of(code, region_files):
  # inner region first (sc: GM inside WM); later regions exclude earlier ones
  regs, taken = {}, None
  for label in sorted(region_files, key=lambda l: (l != "gm", l)):
    paths = region_files[label].get(code)
    if not paths:
      continue
    polys = load_polys(paths)
    if not len(polys):
      continue
    r = shapely.union_all(polys)
    if taken is not None:
      r = shapely.difference(r, taken)
    taken = r if taken is None else shapely.union(taken, r)
    regs[label.removesuffix("_annotations").upper()] = r
  return regs


def section(code, cells, feature_files, region_files, px_um, radii):
  regs = regions_of(code, region_files)
  if not regs:
    log(f"WARNING {code}: no region annotations; skipped")
    return [], []
  cells = cells[cells["region"].isin(list(regs))].reset_index(drop=True)
  pts = shapely.points(cells["x"].to_numpy(), cells["y"].to_numpy())
  pts_tree = shapely.STRtree(pts)
  bands_out, incl_out = [], []
  reg_area = {r: shapely.area(g) * px_um ** 2 for r, g in regs.items()}
  reg_n = cells.groupby(["region", "cell_type2"]).size().rename("n_type_region").reset_index()
  reg_tot = cells.groupby("region").size().rename("n_region")

  for label, by_code in feature_files.items():
    paths = by_code.get(code)
    feats = load_polys(paths) if paths else np.array([], dtype=object)
    if not len(feats):
      log(f"{code}: no {label} features")
      continue
    (qi, fi), d = shapely.STRtree(feats).query_nearest(pts, return_distance=True, all_matches=False)
    dist = np.full(len(cells), np.inf)
    dist[qi] = d * px_um
    pos = as_bool(cells[f"has_{label}"]) if f"has_{label}" in cells else pd.Series(False, index=cells.index)

    # band areas: nested buffers of all features, clipped to each region
    unions = [shapely.union_all(feats)] + [shapely.union_all(shapely.buffer(feats, r / px_um, quad_segs=32)) for r in radii]
    labs = band_labels(radii)
    area_rows = []
    for rname, rg in regs.items():
      a = [shapely.area(shapely.intersection(u, rg)) * px_um ** 2 for u in unions]
      areas = [a[0]] + [a[i + 1] - a[i] for i in range(len(radii))] + [reg_area[rname] - a[-1]]
      area_rows += [{"region": rname, "band": b, "band_area_um2": ar} for b, ar in zip(labs, areas)]
    areas = pd.DataFrame(area_rows)

    df = cells.assign(band=to_band(dist, radii).astype(str), positive=pos.to_numpy())
    cnt = (df.groupby(["region", "band", "cell_type2"]).agg(n_cells=("cell_id", "size"), n_positive=("positive", "sum"))
             .reset_index())
    grid = (areas.merge(reg_n, on="region")[["region", "band", "cell_type2"]])
    b = (grid.merge(cnt, on=["region", "band", "cell_type2"], how="left").fillna({"n_cells": 0, "n_positive": 0})
             .merge(areas, on=["region", "band"]).merge(reg_n, on=["region", "cell_type2"])
             .merge(reg_tot, left_on="region", right_index=True))
    b["n_band"] = b.groupby(["region", "band"])["n_cells"].transform("sum")
    b["region_area_um2"] = b["region"].map(reg_area)
    b["exp_label"] = b["n_band"] * b["n_type_region"] / b["n_region"]
    b["exp_area"] = b["n_type_region"] * b["band_area_um2"] / b["region_area_um2"]
    b["oe_label"] = b["n_cells"] / b["exp_label"]
    b["oe_area"] = b["n_cells"] / b["exp_area"]
    b.insert(0, "label", label)
    b.insert(0, "code", code)
    bands_out.append(b)

    # per inclusion: every cell within the largest radius of each feature (not exclusive)
    # polygonal buffers sit slightly inside the true distance; pad, then filter on exact distance
    fidx, pidx = pts_tree.query(shapely.buffer(feats, radii[-1] * 1.02 / px_um), predicate="intersects")
    dd = shapely.distance(feats[fidx], pts[pidx]) * px_um
    keep = dd <= radii[-1]
    near = pd.DataFrame({"inclusion": fidx[keep], "ring": to_band(dd[keep], radii).astype(str),
                         "cell_type2": cells["cell_type2"].to_numpy()[pidx[keep]],
                         "positive": pos.to_numpy()[pidx[keep]]})
    rc = (near.groupby(["inclusion", "ring", "cell_type2"]).agg(n_cells=("positive", "size"), n_positive=("positive", "sum"))
              .reset_index())
    cen = shapely.centroid(feats)
    inc_region = np.full(len(feats), "outside", dtype=object)
    for rname, rg in regs.items():
      inc_region[shapely.contains_xy(rg, shapely.get_x(cen), shapely.get_y(cen))] = rname
    meta = pd.DataFrame({"inclusion": np.arange(len(feats)), "inclusion_area_um2": shapely.area(feats) * px_um ** 2,
                         "inclusion_region": inc_region, "x": shapely.get_x(cen), "y": shapely.get_y(cen)})
    rc = meta.merge(rc, on="inclusion", how="left")
    rc.insert(0, "label", label)
    rc.insert(0, "code", code)
    incl_out.append(rc)
    log(f"{code}: {label} {len(feats)} features, {int((dist <= radii[-1]).sum())} cells within {radii[-1]} um")
  return bands_out, incl_out


def worker(code, args, cells, feature_files, region_files, px_um):
  try:
    return section(code, cells, feature_files, region_files, px_um, args.radii)
  except Exception as e:
    log(f"ERROR {code}: {e}")
    raise


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--codes", nargs="*")
  p.add_argument("--cells", default="data/09_qc_dimred/filtered_data.csv")
  p.add_argument("--quant_dir", default="data/01_quantify_cells/")
  p.add_argument("--pathology_dir", default="annotations/pathology")
  p.add_argument("--region_dir", default="annotations/region")
  p.add_argument("--radii", type=float, nargs="+", default=[10, 25, 50, 100])
  p.add_argument("--px_um", type=float)
  p.add_argument("--data_dir", default="data/13_inclusion_neighborhood/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()
  args.radii = sorted(int(r) if float(r).is_integer() else r for r in args.radii)

  os.makedirs(args.data_dir, exist_ok=True)
  if not os.path.exists(args.cells):
    sys.exit(f"{args.cells} not found; run 09 (its first array task writes it)")
  feature_files = find_files(args.pathology_dir)
  if not feature_files:
    sys.exit(f"No feature geoJSONs under {args.pathology_dir}/<label>/")
  region_files = find_files(args.region_dir)
  head = pd.read_csv(args.cells, nrows=0).columns
  use = ["code", "cell_id", "cell_type2", "region", "x", "y", "group", "sample_id"] + \
        [c for c in head if c in [f"has_{l}" for l in feature_files]]
  cells = pd.read_csv(args.cells, usecols=use, dtype={"code": str})
  codes = args.codes or sorted(cells["code"].unique())
  unknown = set(codes) - set(cells["code"])
  if unknown:
    sys.exit(f"Codes not in {args.cells}: {', '.join(sorted(unknown))}")
  log(f"{len(codes)} sections, {len(cells)} cells, labels {', '.join(sorted(feature_files))}, radii {args.radii} um")

  px = {}
  for c in codes:
    if args.px_um:
      px[c] = args.px_um
      continue
    q = pd.read_parquet(os.path.join(args.quant_dir, f"{c}.parquet"), columns=["n_pixels", "area_um2"])
    q = q[q["n_pixels"] > 0]
    px[c] = float(np.median(np.sqrt(q["area_um2"] / q["n_pixels"])))

  jobs = [(c, args, cells[cells["code"] == c], feature_files, region_files, px[c]) for c in codes]
  with mp.get_context("fork").Pool(max(1, min(args.workers, len(codes)))) as pool:
    res = pool.starmap(worker, jobs)
  meta = cells.groupby("code")[["group", "sample_id"]].first().reset_index()
  bands = pd.concat([b for r in res for b in r[0]], ignore_index=True).merge(meta, on="code")
  incl = pd.concat([i for r in res for i in r[1]], ignore_index=True).merge(meta, on="code")
  bands.to_csv(os.path.join(args.data_dir, "band_counts.csv"), index=False)
  incl.to_parquet(os.path.join(args.data_dir, "inclusion_ring_counts.parquet"), index=False)
  log(f"wrote {len(bands)} band rows, {incl[['code', 'label', 'inclusion']].drop_duplicates().shape[0]} inclusions")
  log("Done")


if __name__ == "__main__":
  main()
