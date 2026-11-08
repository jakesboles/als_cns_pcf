import os
import sys
import glob
import html
import argparse
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
import tifffile
from PIL import Image, ImageDraw

from pcf_utils import log, get_code, find_files, load_cell_geoms

sys.stdout.reconfigure(line_buffering=True)

RENAME = {"Mac2/Galectin-3": "Galectin-3"}
TYPE_COLORS = {"astrocyte": (40, 170, 40), "cd4_tcell": (0, 190, 255), "cd8_tcell": (30, 60, 220),
               "cortical_neuron": (150, 60, 200), "myeloid": (220, 30, 30), "oligodendrocyte": (255, 170, 0)}
REQUIRED = {"myeloid": ["CD11c", "Iba-1", "TMEM119", "CD163"],
            "cortical_neuron": ["MAP-2", "NeuN"],
            "tcell": ["CD3e", "CD4", "CD8"],
            "astrocyte": ["Vimentin", "GFAP", "Collagen IV", "Claudin-5"],
            "oligodendrocyte": ["Olig-2"]}
OTHER = (150, 150, 150)
TARGET = (255, 0, 255)
GAP = 4


def read_region(page, fh, y0, y1, x0, x1):
  H, W = page.shape[-2:]
  out = np.zeros((y1 - y0, x1 - x0), dtype=page.dtype)
  cy0, cy1, cx0, cx1 = max(y0, 0), min(y1, H), max(x0, 0), min(x1, W)
  if cy0 >= cy1 or cx0 >= cx1:
    return out
  th, tw = (page.tilelength, page.tilewidth) if page.is_tiled else (min(page.rowsperstrip or H, H), W)
  ntx = -(-W // tw)
  for ty in range(cy0 // th, (cy1 - 1) // th + 1):
    for tx in range(cx0 // tw, (cx1 - 1) // tw + 1):
      i = ty * ntx + tx
      if not page.databytecounts[i]:
        continue
      fh.seek(page.dataoffsets[i])
      seg, (_, _, sy, sx, _), _ = page.decode(fh.read(page.databytecounts[i]), i, jpegtables=page.jpegtables)
      seg = seg[0, :, :, 0]
      ay0, ay1 = max(cy0, sy), min(cy1, sy + seg.shape[0])
      ax0, ax1 = max(cx0, sx), min(cx1, sx + seg.shape[1])
      if ay0 < ay1 and ax0 < ax1:
        out[ay0 - y0:ay1 - y0, ax0 - x0:ax1 - x0] = seg[ay0 - sy:ay1 - sy, ax0 - sx:ax1 - sx]
  return out


def page(pages, i):
  p = pages[i]
  return p if isinstance(p, tifffile.TiffPage) else p.aspage()


def renamed_columns(path, markers):
  cols = {}
  for c in pq.read_schema(path).names:
    m = RENAME.get(c, c)
    if m in markers:
      cols[c] = m
  return cols


def cluster_sums(code, args, cofactor, members):
  path = os.path.join(args.quant_dir, f"{code}.parquet")
  cols = renamed_columns(path, cofactor)
  d = pd.read_parquet(path, columns=["cell_id"] + list(cols)).rename(columns=cols)
  out = {}
  for group, m in members.items():
    sub = d.merge(m, on="cell_id", how="inner")
    if sub.empty:
      continue
    markers = [c for c in cofactor if c in sub.columns]
    a = np.arcsinh(sub[markers].to_numpy(np.float64) / np.array([cofactor[c] for c in markers]))
    s = pd.DataFrame(a, columns=markers).groupby(sub["cluster"].to_numpy()).sum()
    n = sub.groupby("cluster").size()
    out[group] = (s, n)
  return code, out


def enrichment(sums, counts, top):
  total, N = sums.sum(), counts.sum()
  means = sums.div(counts, axis=0)
  rest_n = (N - counts).replace(0, np.nan)
  rest = (total.to_numpy()[None, :] - sums).div(rest_n, axis=0)
  diff = means - rest.fillna(0)
  tops = {k: list(diff.loc[k].sort_values(ascending=False).index[:top]) for k in diff.index}
  channels = list(dict.fromkeys(m for k in diff.index for m in tops[k]))
  return diff, tops, channels


def find_image(code, image_dir):
  hits = [f for f in glob.glob(os.path.join(image_dir, "*.qptiff")) if get_code(f) == code]
  return hits[0] if len(hits) == 1 else None


def crop_code(code, args, targets, needed, cell_files):
  image = find_image(code, args.image_dir)
  if image is None:
    log(f"WARNING {code}: no unique qptiff in {args.image_dir}; {len(targets)} cells skipped")
    return code, {}
  ch = pd.read_csv(os.path.join(args.channels_dir, f"{code}_channels.csv"))
  idx = {RENAME.get(c, c): int(i) for i, c in zip(ch["index"], ch["channel"])}
  missing = needed - set(idx)
  if missing:
    log(f"WARNING {code}: channels not in image: {', '.join(sorted(missing))}")

  q = pd.read_parquet(os.path.join(args.quant_dir, f"{code}.parquet"),
                      columns=["cell_id", "cell_type", "x", "y", "n_pixels", "area_um2"])
  px_um = float(np.median(np.sqrt(q["area_um2"] / q["n_pixels"])))
  half = int(round(args.window_um / px_um / 2))
  geoms = load_cell_geoms(cell_files, set(q["cell_id"]), code).merge(q[["cell_id", "cell_type"]], on="cell_id")
  geoms = geoms[geoms["geometry"].notna()].reset_index(drop=True)
  tree = shapely.STRtree(geoms["geometry"].to_numpy())
  xy = q.set_index("cell_id")[["x", "y"]]

  out = {}
  with tifffile.TiffFile(image) as tif:
    pages, fh = tif.series[0].pages, tif.filehandle
    for t in targets.itertuples():
      cx, cy = xy.loc[t.cell_id]
      x0, y0 = int(round(cx)) - half, int(round(cy)) - half
      x1, y1 = x0 + 2 * half, y0 + 2 * half
      crops = {c: read_region(page(pages, idx[c]), fh, y0, y1, x0, x1) for c in t.channels if c in idx}
      outlines = []
      for j in tree.query(shapely.box(x0, y0, x1, y1)):
        g = geoms.at[j, "geometry"]
        for part in getattr(g, "geoms", [g]):
          ring = np.asarray(part.exterior.coords) - (x0, y0)
          outlines.append((ring, geoms.at[j, "cell_type"], geoms.at[j, "cell_id"] == t.cell_id))
      out[t.cell_id] = {"crops": crops, "outlines": outlines}
  log(f"{code}: {len(out)} cells cropped, window {2 * half} px ({px_um:.3f} um/px)")
  return code, out


def render(crops, outlines, channels, ranges, scale, quality, path):
  size = next(iter(crops.values())).shape[0] * scale if crops else 1
  strip = Image.new("RGB", (len(channels) * (size + GAP) - GAP, size), (0, 0, 0))
  for k, c in enumerate(channels):
    if c not in crops:
      continue
    lo, hi = ranges[c]
    a = np.clip((crops[c].astype(np.float32) - lo) / max(hi - lo, 1e-9), 0, 1)
    panel = Image.fromarray((a * 255).astype(np.uint8)).resize((size, size), Image.NEAREST).convert("RGB")
    d = ImageDraw.Draw(panel)
    for ring, ct, is_target in sorted(outlines, key=lambda o: o[2]):
      pts = [tuple(p) for p in (ring * scale)]
      d.line(pts, fill=TARGET if is_target else TYPE_COLORS.get(ct, OTHER), width=3 if is_target else 1)
    strip.paste(panel, (k * (size + GAP), 0))
  strip.save(path, quality=quality)


def write_html(path, group, res, channels, required, tops, diff, rows, counts, size):
  ratio = GAP / size
  legend = "".join(f'<span class="sw" style="border-color:rgb{c}"></span>{html.escape(t)} ' for t, c in TYPE_COLORS.items())
  nav = " ".join(f'<a href="#c{k}">{k}</a>' for k in tops)
  parts = [f"""<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(group)} leiden {res}</title>
<style>
:root{{--h:160px}}
body{{font-family:sans-serif;margin:12px;background:#111;color:#ddd}}
a{{color:#8cf}}
.sw{{display:inline-block;width:14px;height:10px;border:3px solid;margin:0 3px 0 8px;vertical-align:middle}}
.top{{position:sticky;top:0;background:#111;padding:6px 0;z-index:2}}
.head{{display:flex;margin-left:260px}}
.h{{width:var(--h);margin-right:calc(var(--h)*{ratio:.5f});font-size:11px;overflow:hidden;white-space:nowrap;text-overflow:ellipsis}}
.h.req{{color:#8cf}}
.h.hi{{color:#ff0;font-weight:bold}}
.row{{display:flex;align-items:center;margin:2px 0}}
.lab{{width:256px;flex:none;font-size:11px;padding-right:4px;word-break:break-all}}
.row img{{height:var(--h);image-rendering:auto}}
h2{{margin:18px 0 4px}}
</style></head><body>
<div class="top"><b>{html.escape(group)}</b>, Leiden {res} &nbsp; panel size <input type="range" min="60" max="400" value="160" oninput="document.documentElement.style.setProperty('--h',this.value+'px')">
&nbsp; headers: <span style="color:#8cf">always shown</span>, <span style="color:#ff0;font-weight:bold">cluster's top markers</span>
&nbsp; outlines: <span class="sw" style="border-color:rgb{TARGET}"></span>sampled cell {legend}<br>clusters: {nav}</div>
"""]
  for k in tops:
    hi = set(tops[k])
    hc = "".join(f'<div class="h{" hi" if c in hi else " req" if c in required else ""}">{html.escape(c)}</div>' for c in channels)
    top_txt = ", ".join(f"{html.escape(m)} (+{diff.loc[k, m]:.2f})" for m in tops[k])
    parts.append(f'<h2 id="c{k}">Cluster {k}</h2><div>{counts[k]} cells; top markers vs rest of group (mean asinh difference): {top_txt}</div>'
                 f'<div class="head">{hc}</div>')
    for r in rows.get(k, []):
      lab = (f'{html.escape(r["cell_id"])}<br>{html.escape(r["code"])} · {html.escape(r["cell_type"])}'
             f'<br>{html.escape(r["flag_status"])} · {r["area_um2"]:.0f} µm²')
      img = f'<img loading="lazy" src="{html.escape(r["img"])}">' if r["img"] else "<i>no image</i>"
      parts.append(f'<div class="row"><div class="lab">{lab}</div>{img}</div>')
  parts.append("</body></html>")
  with open(path, "w") as fh:
    fh.write("\n".join(parts))


def main():
  p = argparse.ArgumentParser()
  p.add_argument("--groups", nargs="*")
  p.add_argument("--resolution", default="0.3")
  p.add_argument("--n_cells", type=int, default=20)
  p.add_argument("--top", type=int, default=4)
  p.add_argument("--window_um", type=float, default=64)
  p.add_argument("--scale", type=int, default=2)
  p.add_argument("--clip", type=float, nargs=2, default=[1, 99.5])
  p.add_argument("--quality", type=int, default=90)
  p.add_argument("--seed", type=int, default=1)
  p.add_argument("--clusters_dir", default="data/07_cell_clustering")
  p.add_argument("--flags", default="data/05_object_creation/cell_flags.parquet")
  p.add_argument("--cofactors", default="results/05_object_creation/marker_cofactors_thresholds.csv")
  p.add_argument("--quant_dir", default="data/01_quantify_cells/")
  p.add_argument("--channels_dir", default="results/01_quantify_cells/")
  p.add_argument("--image_dir", default="pcf_images")
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--results_dir", default="results/08_cluster_montage/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  col = f"leiden_{args.resolution}"
  groups = args.groups or sorted(os.path.basename(os.path.dirname(f))
                                 for f in glob.glob(os.path.join(args.clusters_dir, "*", "clusters.parquet")))
  if not groups:
    sys.exit(f"No clusters.parquet under {args.clusters_dir}/<group>/; run 07")
  cof = pd.read_csv(args.cofactors)
  cofactor = dict(zip(cof["marker"], cof["cofactor"]))
  flags = pd.read_parquet(args.flags, columns=["code", "cell_id", "cell_type", "area_um2", "keep", "flag_reason"])
  flags["flag_status"] = np.where(flags["keep"].isna(), "NA",
                                  np.where(flags["keep"].fillna(False).astype(bool), "kept", flags["flag_reason"].astype(str)))

  members = {}
  for g in groups:
    f = os.path.join(args.clusters_dir, g, "clusters.parquet")
    if not os.path.exists(f):
      sys.exit(f"{f} not found")
    cl = pd.read_parquet(f)
    if col not in cl.columns:
      sys.exit(f"{col} not in {f}; available: {', '.join(c for c in cl.columns if c.startswith('leiden_'))}")
    cl = cl[["cell_id", col]].rename(columns={col: "cluster"})
    cl["cluster"] = cl["cluster"].astype(str)
    cl = cl.merge(flags[["cell_id", "code"]], on="cell_id", how="left")
    if cl["code"].isna().any():
      sys.exit(f"{int(cl['code'].isna().sum())} {g} cells not in {args.flags}; rerun 05/07 consistently")
    members[g] = cl
  bad = {g: [c for c in REQUIRED.get(g, []) if c not in cofactor] for g in groups}
  bad = {g: v for g, v in bad.items() if v}
  if bad:
    sys.exit("REQUIRED markers not among 05's channels: " + "; ".join(f"{g}: {', '.join(v)}" for g, v in bad.items()))
  codes = sorted(set().union(*(set(m["code"]) for m in members.values())))
  log(f"{len(groups)} groups ({', '.join(groups)}), {col}, {len(codes)} sections")

  pool_n = max(1, min(args.workers, len(codes)))
  jobs = [(c, args, cofactor, {g: m.loc[m["code"] == c, ["cell_id", "cluster"]] for g, m in members.items()}) for c in codes]
  sums, counts = defaultdict(list), defaultdict(list)
  with mp.get_context("fork").Pool(pool_n) as pool:
    for code, res in pool.starmap(cluster_sums, jobs):
      for g, (s, n) in res.items():
        sums[g].append(s)
        counts[g].append(n)

  rng = np.random.default_rng(args.seed)
  plan = {}
  targets = []
  for g in groups:
    s = pd.concat(sums[g]).groupby(level=0).sum()
    n = pd.concat(counts[g]).groupby(level=0).sum()
    order = sorted(n.index, key=lambda x: (len(x), x))
    s, n = s.loc[order], n.loc[order]
    diff, tops, enriched = enrichment(s, n, args.top)
    required = REQUIRED.get(g, [])
    channels = required + [c for c in enriched if c not in required]
    out_dir = os.path.join(args.results_dir, g, col)
    os.makedirs(os.path.join(out_dir, "img"), exist_ok=True)
    diff.to_csv(os.path.join(out_dir, "marker_enrichment.csv"))
    pd.DataFrame([(k, r + 1, m, diff.loc[k, m]) for k in tops for r, m in enumerate(tops[k])],
                 columns=["cluster", "rank", "marker", "mean_asinh_diff"]).to_csv(
      os.path.join(out_dir, "top_markers.csv"), index=False)
    plan[g] = {"channels": channels, "required": required, "tops": tops, "diff": diff, "counts": n, "out_dir": out_dir}
    log(f"{g}: {len(order)} clusters, {len(channels)} channels ({len(required)} always shown): {', '.join(channels)}")
    m = members[g]
    for k in order:
      ids = m.loc[m["cluster"] == k, ["cell_id", "code"]]
      pick = ids.iloc[np.sort(rng.choice(len(ids), size=min(args.n_cells, len(ids)), replace=False))]
      for i, r in enumerate(pick.itertuples()):
        targets.append({"group": g, "cluster": k, "i": i, "cell_id": r.cell_id, "code": r.code, "channels": channels})
  targets = pd.DataFrame(targets)

  cell_files = defaultdict(dict)
  for cell_type, by_code in find_files(args.cells_dir).items():
    for code, paths in by_code.items():
      cell_files[code][cell_type] = paths
  jobs = []
  for code, t in targets.groupby("code"):
    needed = set().union(*map(set, t["channels"]))
    jobs.append((code, args, t, needed, dict(cell_files.get(code, {}))))
  crops = {}
  with mp.get_context("fork").Pool(max(1, min(args.workers, len(jobs)))) as pool:
    for code, res in pool.starmap(crop_code, jobs):
      crops.update(res)

  meta = flags.set_index("cell_id")
  for g in groups:
    pl = plan[g]
    t = targets[targets["group"] == g]
    ranges = {}
    for c in pl["channels"]:
      px = [crops[cid]["crops"][c].ravel() for cid in t["cell_id"] if cid in crops and c in crops[cid]["crops"]]
      ranges[c] = tuple(np.percentile(np.concatenate(px), args.clip)) if px else (0, 1)
    pd.DataFrame([(c, *ranges[c]) for c in pl["channels"]], columns=["channel", "display_min", "display_max"]).to_csv(
      os.path.join(pl["out_dir"], "display_ranges.csv"), index=False)
    rows = defaultdict(list)
    size = None
    for r in t.itertuples():
      rel = os.path.join("img", f"c{r.cluster}_{r.i:02d}.webp")
      ok = r.cell_id in crops and crops[r.cell_id]["crops"]
      if ok:
        cr = crops[r.cell_id]
        render(cr["crops"], cr["outlines"], pl["channels"], ranges, args.scale, args.quality, os.path.join(pl["out_dir"], rel))
        size = size or next(iter(cr["crops"].values())).shape[0] * args.scale
      m = meta.loc[r.cell_id]
      rows[r.cluster].append({"cell_id": r.cell_id, "code": r.code, "cell_type": m["cell_type"],
                              "flag_status": m["flag_status"], "area_um2": float(m["area_um2"]), "img": rel if ok else None})
    write_html(os.path.join(pl["out_dir"], "index.html"), g, args.resolution, pl["channels"], pl["required"], pl["tops"],
               pl["diff"], rows, pl["counts"], size or 256)
    log(f"{g}: wrote {os.path.join(pl['out_dir'], 'index.html')}")

  links = "".join(f'<li><a href="{html.escape(g)}/{col}/index.html">{html.escape(g)}</a></li>' for g in groups)
  with open(os.path.join(args.results_dir, f"index_{col}.html"), "w") as fh:
    fh.write(f'<!doctype html><html><head><meta charset="utf-8"><title>Cluster montages {col}</title></head>'
             f'<body style="font-family:sans-serif"><h1>Cluster montages, {col}</h1><ul>{links}</ul></body></html>')
  log("Done")


if __name__ == "__main__":
  main()
