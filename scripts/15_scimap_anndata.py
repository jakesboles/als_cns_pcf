import os
import sys
import argparse
import multiprocessing as mp

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import anndata as ad

sys.stdout.reconfigure(line_buffering=True)

RENAME = {"Mac2/Galectin-3": "Galectin-3"}


def log(msg):
  import time
  print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def load_section(code, quant_dir, markers, ids, raw):
  path = os.path.join(quant_dir, f"{code}.parquet")
  cols = {c: RENAME.get(c, c) for c in pq.read_schema(path).names if RENAME.get(c, c) in markers}
  missing = set(markers) - set(cols.values())
  if missing:
    raise ValueError(f"{code}: channels missing from 01 output: {', '.join(sorted(missing))}")
  d = pd.read_parquet(path, columns=["cell_id", "n_pixels", "area_um2"] + list(cols)).rename(columns=cols)
  q = d[d["n_pixels"] > 0]
  px_um = float(np.median(np.sqrt(q["area_um2"] / q["n_pixels"])))
  d = d.set_index("cell_id").reindex(ids)
  if d[markers[0]].isna().any():
    raise ValueError(f"{code}: {int(d[markers[0]].isna().sum())} cells not in 01 output")
  return code, px_um, (d[markers].to_numpy(np.float32) if raw else None)


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--cells", default="data/09_qc_dimred/filtered_data.csv")
  p.add_argument("--quant_dir", default="data/01_quantify_cells/")
  p.add_argument("--markers", default="results/05_object_creation/marker_cofactors_thresholds.csv")
  p.add_argument("--z_cap", type=float, default=10)
  p.add_argument("--out", default="data/15_scimap_anndata/cells.h5ad")
  p.add_argument("--no_intensities", action="store_true")
  p.add_argument("--raw_layer", action="store_true")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  for f in [args.cells, args.markers]:
    if not os.path.exists(f):
      sys.exit(f"{f} not found")
  obs = pd.read_csv(args.cells, dtype={"code": str}, low_memory=False)
  # 05's kept channels
  markers = pd.read_csv(args.markers)["marker"].tolist()
  codes = sorted(obs["code"].unique())
  log(f"{len(obs)} cells, {len(codes)} sections, {len(markers)} markers")

  jobs = [(c, args.quant_dir, markers, obs.loc[obs["code"] == c, "cell_id"].to_numpy(), not args.no_intensities)
          for c in codes]
  with mp.get_context("fork").Pool(max(1, min(args.workers, len(codes)))) as pool:
    res = {c: (px, m) for c, px, m in pool.starmap(load_section, jobs)}

  obs = pd.concat([obs[obs["code"] == c] for c in codes], ignore_index=True)
  px = obs["code"].map({c: res[c][0] for c in codes}).to_numpy()
  obs["X_centroid"] = obs["x"] * px
  obs["Y_centroid"] = obs["y"] * px
  obs["imageid"] = obs["code"]
  obs["phenotype"] = obs["cell_type2"]
  for c in obs.columns:
    if c.startswith("has_"):
      obs[c] = obs[c].map({True: True, False: False, "TRUE": True, "FALSE": False}).astype("boolean")
    elif c != "cell_id" and (obs[c].dtype == object or pd.api.types.is_string_dtype(obs[c])):
      obs[c] = obs[c].astype("category")
  obs.index = obs["cell_id"].astype(str).to_numpy()

  if args.no_intensities:
    X, layers, var = np.zeros((len(obs), 0), np.float32), {}, None
  else:
    raw = np.vstack([res[c][1] for c in codes])
    X = np.log1p(np.maximum(raw, 0))
    mu, sd = X.mean(axis=0, dtype=np.float64), X.std(axis=0, dtype=np.float64)
    sd[~(sd > 0)] = 1
    z = np.clip((X - mu.astype(np.float32)) / sd.astype(np.float32), -args.z_cap, args.z_cap)
    layers = {"log_dataset_z": z, **({"raw": raw} if args.raw_layer else {})}
    var = pd.DataFrame({"log1p_mean": mu, "log1p_sd": sd}, index=markers)
  adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
  adata.obsm["spatial"] = obs[["X_centroid", "Y_centroid"]].to_numpy()
  adata.uns["units"] = "X_centroid/Y_centroid and obsm['spatial'] in um; x/y in px; X = log1p(max(raw, 0)); layers['log_dataset_z'] = X z-scored per marker over all cells, capped at +-z_cap"
  os.makedirs(os.path.dirname(args.out), exist_ok=True)
  adata.write_h5ad(args.out)
  log(f"wrote {args.out}: {adata.n_obs} cells x {adata.n_vars} markers")


if __name__ == "__main__":
  main()
