import os
import re
import sys
import glob
import time
import argparse
import multiprocessing as mp

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from rasterio.features import rasterize
from rasterio.transform import Affine
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(line_buffering=True)

CODE_RE = re.compile(r"JSB(\d+-\d+)_")
BAND_ROWS = 2048

LABELS = None
OVERLAP_MASKS = None
N_CELLS = None
IMAGE_PATH = None


def log(msg):
  print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def get_code(path):
  m = CODE_RE.search(os.path.basename(path))
  return m.group(1) if m else None


def to_polygon(g):
  if g is None or g.is_empty:
    return None
  if g.geom_type == "LineString":
    coords = list(dict.fromkeys(g.coords))
    if len(coords) < 3:
      return None
    g = Polygon(coords)
  if not g.is_valid:
    g = shapely.make_valid(g)
  if g.geom_type == "GeometryCollection":
    polys = [p for p in g.geoms if p.geom_type in ("Polygon", "MultiPolygon")]
    g = shapely.union_all(polys) if polys else None
  if g is None or g.geom_type not in ("Polygon", "MultiPolygon") or g.area == 0:
    return None
  return g


def load_cells(cell_files, code):
  frames = []
  for cell_type, path in cell_files:
    gdf = gpd.read_file(path)
    gdf = gdf[["geometry"]].set_crs(None, allow_override=True)
    gdf["feature_index"] = np.arange(1, len(gdf) + 1)
    gdf["cell_type"] = cell_type
    gdf["geometry"] = gdf.geometry.apply(to_polygon)
    n_bad = gdf.geometry.isna().sum()
    if n_bad:
      log(f"  {cell_type}: dropped {n_bad} empty/degenerate polygons")
    frames.append(gdf[gdf.geometry.notna()])
    log(f"  {cell_type}: {len(frames[-1])} cells")
  cells = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), geometry="geometry")
  cells["cell_id"] = code + "_" + cells["cell_type"] + "_" + cells["feature_index"].astype(str)
  cells["label"] = np.arange(1, len(cells) + 1, dtype=np.int32)
  return cells


def find_overlaps(geoms):
  tree = shapely.STRtree(geoms)
  i, j = tree.query(geoms, predicate="intersects")
  keep = i < j
  i, j = i[keep], j[keep]
  inter = shapely.area(shapely.intersection(geoms[i], geoms[j])) > 0
  return np.unique(np.concatenate([i[inter], j[inter]]))


def cell_mask(geom, shape):
  xmin, ymin, xmax, ymax = geom.bounds
  c0, r0 = max(int(np.floor(xmin)), 0), max(int(np.floor(ymin)), 0)
  c1, r1 = min(int(np.ceil(xmax)) + 1, shape[1]), min(int(np.ceil(ymax)) + 1, shape[0])
  if c1 <= c0 or r1 <= r0:
    return None
  m = rasterize([(geom, 1)], out_shape=(r1 - r0, c1 - c0),
                transform=Affine.translation(c0, r0), fill=0, dtype="uint8").astype(bool)
  return (r0, r1, c0, c1, m)


def channel_names(series):
  names = []
  for k, page in enumerate(series.pages):
    desc = page.description or ""
    m = re.search(r"<Biomarker>(.*?)</Biomarker>", desc) or re.search(r"<Name>(.*?)</Name>", desc)
    names.append(m.group(1).strip() if m else f"ch{k}")
  seen = {}
  out = []
  for nm in names:
    seen[nm] = seen.get(nm, 0) + 1
    out.append(nm if seen[nm] == 1 else f"{nm}_{seen[nm]}")
  return out


def pixel_size_um(page):
  try:
    num, den = page.tags["XResolution"].value
    unit = page.tags["ResolutionUnit"].value
    per_unit = num / den
    if unit == 3:
      return 1e4 / per_unit
    if unit == 2:
      return 25400 / per_unit
  except (KeyError, ZeroDivisionError, TypeError, ValueError):
    pass
  return np.nan


def quantify_channel(c):
  with tifffile.TiffFile(IMAGE_PATH) as tif:
    plane = tif.series[0].pages[c].asarray()
  sums = np.zeros(N_CELLS + 1)
  for r0 in range(0, plane.shape[0], BAND_ROWS):
    r1 = r0 + BAND_ROWS
    sums += np.bincount(LABELS[r0:r1].ravel(), weights=plane[r0:r1].ravel(), minlength=N_CELLS + 1)
  over = {lab: plane[r0:r1, c0:c1][m].mean() for lab, (r0, r1, c0, c1, m) in OVERLAP_MASKS.items()}
  return c, sums, over


def qc_overlay(path, dapi_idx, cells, shape, out_png, half=1000):
  cx, cy = np.median(cells["x"]), np.median(cells["y"])
  c0, r0 = int(max(cx - half, 0)), int(max(cy - half, 0))
  c1, r1 = int(min(cx + half, shape[1])), int(min(cy + half, shape[0]))
  with tifffile.TiffFile(path) as tif:
    crop = tif.series[0].pages[dapi_idx].asarray()[r0:r1, c0:c1]
  sub = cells.cx[c0:c1, r0:r1]
  fig, ax = plt.subplots(figsize=(10, 10))
  ax.imshow(crop, cmap="gray", vmax=np.percentile(crop, 99.5), extent=[c0, c1, r1, r0])
  sub.boundary.plot(ax=ax, linewidth=0.4, color="yellow")
  ax.set_xlim(c0, c1)
  ax.set_ylim(r1, r0)
  ax.set_title(f"{os.path.basename(path)}: {len(sub)} cells in crop")
  fig.savefig(out_png, dpi=150, bbox_inches="tight")
  plt.close(fig)


def main():
  global LABELS, OVERLAP_MASKS, N_CELLS, IMAGE_PATH

  p = argparse.ArgumentParser()
  p.add_argument("--image", required=True)
  p.add_argument("--cells_dir", default="annotations/cells")
  p.add_argument("--data_dir", default="data/01_quantify_cells/")
  p.add_argument("--results_dir", default="results/01_quantify_cells/")
  p.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
  args = p.parse_args()

  IMAGE_PATH = args.image
  if not os.path.isfile(IMAGE_PATH):
    sys.exit(f"Image not found: {IMAGE_PATH!r}")
  code = get_code(IMAGE_PATH)
  if code is None:
    sys.exit(f"Can't parse JSB<code>_ from {IMAGE_PATH}")
  os.makedirs(args.data_dir, exist_ok=True)
  os.makedirs(args.results_dir, exist_ok=True)
  log(f"{code}: {IMAGE_PATH}")

  cell_files = sorted(
    (os.path.basename(os.path.dirname(f)), f)
    for f in glob.glob(os.path.join(args.cells_dir, "*", "*.geojson"))
    if get_code(f) == code
  )
  if not cell_files:
    sys.exit(f"No cell geoJSONs for {code} under {args.cells_dir}/<cell_type>/")

  with tifffile.TiffFile(IMAGE_PATH) as tif:
    series = tif.series[0]
    if series.axes[-2:] != "YX" or len(series.shape) != 3:
      sys.exit(f"Unexpected series axes {series.axes} / shape {series.shape}")
    shape = series.shape[1:]
    channels = channel_names(series)
    px_um = pixel_size_um(series.pages[0])
  log(f"Image {shape[0]} x {shape[1]} px, {len(channels)} channels, {px_um:.4f} um/px")

  cells = load_cells(cell_files, code)
  N_CELLS = len(cells)
  geoms = cells.geometry.values

  xmin, ymin, xmax, ymax = cells.total_bounds
  if xmin < 0 or ymin < 0 or xmax > shape[1] or ymax > shape[0]:
    log(f"WARNING: cell bounds {cells.total_bounds} exceed image {shape[1]} x {shape[0]}")

  overlap_idx = find_overlaps(geoms)
  cells["overlap"] = False
  cells.iloc[overlap_idx, cells.columns.get_loc("overlap")] = True
  log(f"{N_CELLS} cells, {len(overlap_idx)} overlapping another cell")

  free = ~cells["overlap"].values
  if free.any():
    LABELS = rasterize(zip(geoms[free], cells["label"].values[free]), out_shape=shape,
                       transform=Affine.identity(), fill=0, dtype="int32")
  else:
    LABELS = np.zeros(shape, dtype=np.int32)
  log("Rasterized label image")

  n_px = np.zeros(N_CELLS + 1, dtype=np.int64)
  for r0 in range(0, shape[0], BAND_ROWS):
    n_px += np.bincount(LABELS[r0:r0 + BAND_ROWS].ravel(), minlength=N_CELLS + 1)
  OVERLAP_MASKS = {}
  for k in overlap_idx:
    cm = cell_mask(geoms[k], shape)
    if cm is not None and cm[4].any():
      lab = int(cells["label"].iat[k])
      OVERLAP_MASKS[lab] = cm
      n_px[lab] = cm[4].sum()

  cent = shapely.centroid(geoms)
  cells["x"] = shapely.get_x(cent)
  cells["y"] = shapely.get_y(cent)
  cells["n_pixels"] = n_px[cells["label"].values]
  cells["area_um2"] = cells["n_pixels"] * px_um ** 2
  n_zero = int((cells["n_pixels"] == 0).sum())
  if n_zero:
    log(f"{n_zero} cells cover no pixel centers; dropped from the parquet, listed in {code}_dropped_cells.csv")

  log(f"Quantifying {len(channels)} channels with {args.workers} workers")
  means = np.full((N_CELLS, len(channels)), np.nan, dtype=np.float32)
  lab = cells["label"].values
  with np.errstate(invalid="ignore", divide="ignore"):
    with mp.get_context("fork").Pool(args.workers) as pool:
      for c, sums, over in pool.imap_unordered(quantify_channel, range(len(channels))):
        col = sums / n_px
        for l, v in over.items():
          col[l] = v
        means[:, c] = col[lab]
        log(f"  done {channels[c]}")

  out = pd.DataFrame({
    "code": code,
    "cell_id": cells["cell_id"].values,
    "cell_type": cells["cell_type"].values,
    "feature_index": cells["feature_index"].values,
    "x": cells["x"].values,
    "y": cells["y"].values,
    "n_pixels": cells["n_pixels"].values,
    "area_um2": cells["area_um2"].values,
    "overlap": cells["overlap"].values,
  })
  out = pd.concat([out, pd.DataFrame(means, columns=channels)], axis=1)
  pd.DataFrame({"index": range(len(channels)), "channel": channels}).to_csv(
    os.path.join(args.results_dir, f"{code}_channels.csv"), index=False)
  out.groupby("cell_type").agg(
    n_cells=("cell_id", "size"),
    n_overlap=("overlap", "sum"),
    n_zero_px=("n_pixels", lambda s: int((s == 0).sum())),
    median_area_um2=("area_um2", "median"),
  ).to_csv(os.path.join(args.results_dir, f"{code}_qc.csv"))

  zero = out["n_pixels"].values == 0
  bx0, by0, bx1, by1 = shapely.bounds(geoms).T
  out.loc[zero, ["cell_id", "cell_type", "feature_index", "x", "y"]].assign(
    poly_area_px=shapely.area(geoms[zero]),
    inside_image=(bx0[zero] >= 0) & (by0[zero] >= 0) & (bx1[zero] <= shape[1]) & (by1[zero] <= shape[0]),
  ).to_csv(os.path.join(args.results_dir, f"{code}_dropped_cells.csv"), index=False)
  out[~zero].to_parquet(os.path.join(args.data_dir, f"{code}.parquet"), index=False)

  dapi = [k for k, nm in enumerate(channels) if "DAPI" in nm.upper()]
  qc_overlay(IMAGE_PATH, dapi[0] if dapi else 0, cells, shape,
             os.path.join(args.results_dir, f"{code}_overlay.png"))
  log("Done")


if __name__ == "__main__":
  main()
