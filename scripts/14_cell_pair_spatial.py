import os
import sys
import zlib
import warnings
import argparse
import multiprocessing as mp

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from pcf_utils import log

sys.stdout.reconfigure(line_buffering=True)


def nearest_medians(xy, t, n_types):
  # median distance from each cell of type a to the nearest other cell of type b -> (a, b)
  out = np.full((n_types, n_types), np.nan)
  order = np.argsort(t, kind="stable")
  bounds = np.searchsorted(t[order], np.arange(n_types + 1))
  for b in range(n_types):
    tb = np.flatnonzero(t == b)
    if len(tb) == 0:
      continue
    k = 2 if len(tb) > 1 else 1
    d, _ = cKDTree(xy[tb]).query(xy, k=k, workers=1)
    if k == 2:
      d = np.where(t == b, d[:, 1], d[:, 0])
    elif len(tb) == 1:
      d = np.where(t == b, np.nan, d)
    d = d[order]
    for a in range(n_types):
      seg = d[bounds[a]:bounds[a + 1]]
      if len(seg) and not np.all(np.isnan(seg)):
        out[a, b] = np.nanmedian(seg)
  return out


def pair_counts(pairs, t, n_types):
  # number of b cells within r of each a cell, summed over a cells -> (a, b); symmetric
  c = np.bincount(t[pairs[:, 0]] * n_types + t[pairs[:, 1]], minlength=n_types ** 2).reshape(n_types, n_types)
  return c + c.T


def summarise(obs, null):
  # pairs with a single cell (e.g. one T cell in a section) are all NaN
  with warnings.catch_warnings(), np.errstate(invalid="ignore", divide="ignore"):
    warnings.simplefilter("ignore", RuntimeWarning)
    mu = np.nanmean(null, axis=0)
    sd = np.nanstd(null, axis=0, ddof=1)
    ext = np.sum(np.abs(null - mu) >= np.abs(obs - mu), axis=0)
    p = (1 + ext) / (1 + np.sum(~np.isnan(null), axis=0))
    z = (obs - mu) / sd
  return mu, sd, z, np.where(np.isnan(obs), np.nan, p)


def section(code, cells, types, px_um, args):
  n_types = len(types)
  xy = cells[["x", "y"]].to_numpy(float) * px_um
  t = pd.Categorical(cells["cell_type2"], categories=types).codes.astype(np.int64)
  n = np.bincount(t, minlength=n_types)
  rng = np.random.default_rng([args.seed, zlib.crc32(code.encode())])
  pairs = {r: cKDTree(xy).query_pairs(r, output_type="ndarray") for r in args.radii}

  obs_d = nearest_medians(xy, t, n_types)
  obs_c = {r: pair_counts(p, t, n_types) for r, p in pairs.items()}
  null_d = np.empty((args.n_perm, n_types, n_types))
  null_c = {r: np.empty((args.n_perm, n_types, n_types)) for r in args.radii}
  for i in range(args.n_perm):
    tp = rng.permutation(t)
    null_d[i] = nearest_medians(xy, tp, n_types)
    for r, p in pairs.items():
      null_c[r][i] = pair_counts(p, tp, n_types)

  a_idx, b_idx = np.meshgrid(np.arange(n_types), np.arange(n_types), indexing="ij")
  base = pd.DataFrame({"code": code, "from_type": np.array(types)[a_idx.ravel()], "to_type": np.array(types)[b_idx.ravel()],
                       "n_from": n[a_idx.ravel()], "n_to": n[b_idx.ravel()]})

  mu, sd, z, p = summarise(obs_d, null_d)
  dist = base.assign(obs_median_um=obs_d.ravel(), null_median_um=mu.ravel(), null_sd_um=sd.ravel(),
                     log2_ratio=np.log2(obs_d / mu).ravel(), z=z.ravel(), p_perm=p.ravel())

  counts = []
  for r in args.radii:
    mu, sd, z, p = summarise(obs_c[r].astype(float), null_c[r])
    with np.errstate(invalid="ignore", divide="ignore"):
      nb = obs_c[r] / n[:, None]
      log2_oe = np.log2((obs_c[r] + 0.5) / (mu + 0.5))
    counts.append(base.assign(radius_um=r, obs_pairs=obs_c[r].ravel(), obs_mean_nbrs=nb.ravel(), null_pairs=mu.ravel(),
                              null_sd=sd.ravel(), log2_oe=log2_oe.ravel(), z=z.ravel(), p_perm=p.ravel()))
  counts = pd.concat(counts, ignore_index=True)
  keep = lambda d: d[(d["n_from"] > 0) & (d["n_to"] > 0)]
  log(f"{code}: {len(cells)} cells, {args.n_perm} permutations, {', '.join(f'{len(p)} pairs <= {r} um' for r, p in pairs.items())}")
  return keep(dist), keep(counts)


def worker(code, cells, types, px_um, args):
  try:
    return section(code, cells, types, px_um, args)
  except Exception as e:
    log(f"ERROR {code}: {e}")
    raise


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--region", required=True)
  p.add_argument("--codes", nargs="*")
  p.add_argument("--cells", default="data/09_qc_dimred/filtered_data.csv")
  p.add_argument("--quant_dir", default="data/01_quantify_cells/")
  p.add_argument("--radii", type=float, nargs="+", default=[30])
  p.add_argument("--n_perm", type=int, default=1000)
  p.add_argument("--seed", type=int, default=1)
  p.add_argument("--px_um", type=float)
  p.add_argument("--data_dir", default="data/14_cell_pair_spatial/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()
  args.radii = sorted(int(r) if float(r).is_integer() else r for r in args.radii)
  args.region = args.region.upper()
  out_dir = os.path.join(args.data_dir, args.region)
  os.makedirs(out_dir, exist_ok=True)

  if not os.path.exists(args.cells):
    sys.exit(f"{args.cells} not found; run 09 (its first array task writes it)")
  cells = pd.read_csv(args.cells, usecols=["code", "cell_id", "cell_type2", "region", "x", "y", "group", "sample_id"],
                      dtype={"code": str})
  if args.region not in set(cells["region"]):
    sys.exit(f"No cells with region {args.region}; regions present: {', '.join(sorted(cells['region'].dropna().unique()))}")
  meta = cells.groupby("code")[["group", "sample_id"]].first().reset_index()
  cells = cells[cells["region"] == args.region]
  types = sorted(cells["cell_type2"].unique())
  codes = args.codes or sorted(cells["code"].unique())
  unknown = set(codes) - set(cells["code"])
  if unknown:
    sys.exit(f"Codes with no {args.region} cells in {args.cells}: {', '.join(sorted(unknown))}")
  log(f"{args.region}: {len(codes)} sections, {len(cells)} cells, types {', '.join(types)}, radii {args.radii} um, {args.n_perm} permutations")

  px = {}
  for c in codes:
    if args.px_um:
      px[c] = args.px_um
      continue
    q = pd.read_parquet(os.path.join(args.quant_dir, f"{c}.parquet"), columns=["n_pixels", "area_um2"])
    q = q[q["n_pixels"] > 0]
    px[c] = float(np.median(np.sqrt(q["area_um2"] / q["n_pixels"])))

  jobs = [(c, cells[cells["code"] == c].reset_index(drop=True), types, px[c], args) for c in codes]
  with mp.get_context("fork").Pool(max(1, min(args.workers, len(codes)))) as pool:
    res = pool.starmap(worker, jobs)
  for k, name in enumerate(["nearest_distance.csv", "neighbor_counts.csv"]):
    out = pd.concat([r[k] for r in res], ignore_index=True).merge(meta, on="code")
    out.insert(1, "region", args.region)
    out.to_csv(os.path.join(out_dir, name), index=False)
    log(f"wrote {len(out)} rows to {os.path.join(out_dir, name)}")
  log("Done")


if __name__ == "__main__":
  main()
