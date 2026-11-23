import os
import sys
import time
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm
from scipy.stats import norm

sys.stdout.reconfigure(line_buffering=True)


def log(msg):
  print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def save_png(fn, path, **kw):
  d, f = os.path.split(path)
  fn(saveDir=d, fileName=f, **kw)
  plt.close("all")


def nearest_medians(adata, meta_cols):
  d = adata.uns["spatial_distance"].join(adata.obs[meta_cols])
  return (d.melt(id_vars=meta_cols, var_name="to_type", value_name="dist_um")
           .groupby(meta_cols + ["to_type"], observed=True)["dist_um"]
           .median().rename("median_dist_um").reset_index()
           .rename(columns={"imageid": "code", "phenotype": "from_type"}))


def interaction_table(si, meta):
  codes = [c for c in si.columns if c not in ("phenotype", "neighbour_phenotype") and not c.startswith("pvalue_")]
  rows = []
  for c in codes:
    t = si[["phenotype", "neighbour_phenotype", c, "pvalue_" + c]].copy()
    t.columns = ["from_type", "to_type", "score", "p"]
    t["code"] = c
    rows.append(t)
  out = pd.concat(rows, ignore_index=True)
  out["z"] = np.sign(out["score"]) * norm.isf(out["p"].clip(lower=1e-300) / 2)
  return out.merge(meta, on="code")


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--region", required=True)
  p.add_argument("--codes", nargs="*")
  p.add_argument("--h5ad", default="data/15_scimap_anndata/cells.h5ad")
  p.add_argument("--radii", type=float, nargs="+", default=[30])
  p.add_argument("--n_perm", type=int, default=1000)
  p.add_argument("--k", type=int, default=10)
  p.add_argument("--nbhd_radius", type=float)
  p.add_argument("--no_neighborhoods", action="store_true")
  p.add_argument("--results_dir", default="results/16_scimap_spatial/")
  p.add_argument("--data_dir", default="data/16_scimap_spatial/")
  args = p.parse_args()
  args.region = args.region.upper()
  fmt = lambda r: f"{r:g}"
  res_dir = os.path.join(args.results_dir, args.region)
  dat_dir = os.path.join(args.data_dir, args.region)
  os.makedirs(res_dir, exist_ok=True)
  os.makedirs(dat_dir, exist_ok=True)

  if not os.path.exists(args.h5ad):
    sys.exit(f"{args.h5ad} not found; run 15 first")
  backed = ad.read_h5ad(args.h5ad, backed="r")
  keep = (backed.obs["region"] == args.region).to_numpy()
  if args.codes:
    unknown = set(args.codes) - set(backed.obs["code"].astype(str))
    if unknown:
      sys.exit(f"Codes not in {args.h5ad}: {', '.join(sorted(unknown))}")
    keep &= backed.obs["code"].astype(str).isin(args.codes).to_numpy()
  if not keep.any():
    sys.exit(f"No cells with region {args.region}")
  adata = backed[keep].to_memory()
  backed.file.close()
  for c in adata.obs.columns:
    if isinstance(adata.obs[c].dtype, pd.CategoricalDtype):
      adata.obs[c] = adata.obs[c].cat.remove_unused_categories()
  meta = (adata.obs[["imageid", "group", "sample_id"]].drop_duplicates()
          .rename(columns={"imageid": "code"}).astype(str))
  log(f"{args.region}: {adata.n_obs} cells, {adata.obs['imageid'].nunique()} sections, radii {args.radii} um, "
      f"{args.n_perm} permutations")
  pd.DataFrame({"param": ["region", "radii_um", "n_perm", "k", "nbhd_radius_um", "n_cells", "n_sections"],
                "value": [args.region, " ".join(map(fmt, args.radii)), args.n_perm, args.k,
                          fmt(args.nbhd_radius or args.radii[0]), adata.n_obs, adata.obs["imageid"].nunique()]}
               ).to_csv(os.path.join(res_dir, "run_params.csv"), index=False)

  save_png(sm.pl.stacked_barplot, os.path.join(res_dir, "composition_by_group.png"),
           adata=adata, x_axis="group", y_axis="phenotype", method="percent")
  save_png(sm.pl.stacked_barplot, os.path.join(res_dir, "composition_by_section.png"),
           adata=adata, x_axis="imageid", y_axis="phenotype", method="percent")

  # A: nearest distance
  t = time.time()
  adata = sm.tl.spatial_distance(adata, phenotype="phenotype", imageid="imageid", label="spatial_distance",
                                 verbose=False)
  log(f"spatial_distance: {time.time() - t:.0f} s")
  nearest_medians(adata, ["imageid", "group", "sample_id", "phenotype"]).to_csv(
    os.path.join(res_dir, "nearest_distance_medians.csv"), index=False)
  save_png(sm.pl.spatial_distance, os.path.join(res_dir, "spatial_distance_heatmap.png"),
           adata=adata, heatmap_summarize=True, log=True)
  for ct in adata.obs["phenotype"].cat.categories:
    save_png(sm.pl.spatial_distance, os.path.join(res_dir, f"distance_from_{ct}_by_group.png"),
             adata=adata, method="numeric", distance_from=ct, imageid="group")

  # B: interactions within each radius
  for r in args.radii:
    t = time.time()
    label = f"spatial_interaction_r{fmt(r)}"
    adata = sm.tl.spatial_interaction(adata, phenotype="phenotype", imageid="imageid", method="radius", radius=r,
                                      permutation=args.n_perm, label=label, verbose=False)
    log(f"spatial_interaction r = {fmt(r)} um: {time.time() - t:.0f} s")
    interaction_table(adata.uns[label], meta).to_csv(os.path.join(res_dir, f"interaction_r{fmt(r)}.csv"),
                                                     index=False)
    # summarised: grey where the mean p over sections > 0.05; _all: unmasked; _by_section: pair x section
    save_png(sm.pl.spatial_interaction, os.path.join(res_dir, f"interaction_r{fmt(r)}.png"),
             adata=adata, spatial_interaction=label, p_val=0.05)
    save_png(sm.pl.spatial_interaction, os.path.join(res_dir, f"interaction_r{fmt(r)}_all.png"),
             adata=adata, spatial_interaction=label, p_val=1)
    if adata.obs["imageid"].nunique() > 1:
      save_png(sm.pl.spatial_interaction, os.path.join(res_dir, f"interaction_r{fmt(r)}_by_section.png"),
               adata=adata, spatial_interaction=label, p_val=0.05, summarize_plot=False)

  # neighbourhoods
  if not args.no_neighborhoods:
    r = args.nbhd_radius or args.radii[0]
    adata = sm.tl.spatial_count(adata, phenotype="phenotype", imageid="imageid", method="radius", radius=r,
                                label="spatial_count", verbose=False)
    adata = sm.tl.spatial_cluster(adata, df_name="spatial_count", method="kmeans", k=args.k, label="neighborhood",
                                  verbose=False)
    adata.obs["neighborhood"] = adata.obs["neighborhood"].astype(str).astype("category")
    for x, y, f in [("neighborhood", "phenotype", "neighborhood_composition.png"),
                    ("imageid", "neighborhood", "neighborhood_by_section.png"),
                    ("group", "neighborhood", "neighborhood_by_group.png")]:
      save_png(sm.pl.stacked_barplot, os.path.join(res_dir, f), adata=adata, x_axis=x, y_axis=y, method="percent")
    (adata.obs.groupby(["imageid", "neighborhood"], observed=False).size().rename("n_cells").reset_index()
       .rename(columns={"imageid": "code"}).astype({"code": str}).merge(meta, on="code")
       .to_csv(os.path.join(res_dir, "neighborhood_counts.csv"), index=False))
    log(f"neighbourhoods: k = {args.k}, radius {fmt(r)} um")

  adata.write_h5ad(os.path.join(dat_dir, "cells_scimap.h5ad"))
  log(f"wrote {os.path.join(dat_dir, 'cells_scimap.h5ad')}")
  log("Done")


if __name__ == "__main__":
  main()
