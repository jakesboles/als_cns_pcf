import os
import re
import time

import shapely
from shapely.geometry import Polygon

CODE_RE = re.compile(r"JSB(\d+-\d+)_")


def log(msg):
  print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def get_code(path):
  m = CODE_RE.search(os.path.basename(path))
  return m.group(1) if m else None


def to_polygon(g):
  if g is None or g.is_empty:
    return None
  if g.geom_type == "LineString":
    coords = list(g.coords)
    if len(set(coords)) < 3:
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


def find_files(root, by_code_only=False):
  import glob
  from collections import defaultdict
  pattern = "*.geojson" if by_code_only else os.path.join("*", "*.geojson")
  files = defaultdict(lambda: defaultdict(list))
  for f in sorted(glob.glob(os.path.join(root, pattern))):
    code = get_code(f)
    if code is None:
      log(f"WARNING: can't parse code from {f}, skipping")
      continue
    label = "" if by_code_only else os.path.basename(os.path.dirname(f))
    files[label][code].append(f)
  out = {k: dict(v) for k, v in files.items()}
  return out.get("", {}) if by_code_only else out


def load_cell_geoms(cell_files, keep_ids, code):
  import geopandas as gpd
  import pandas as pd
  ids, geoms = [], []
  for cell_type, paths in sorted(cell_files.items()):
    for p in paths:
      for k, g in enumerate(gpd.read_file(p).geometry, start=1):
        ids.append(f"{code}_{cell_type}_{k}")
        geoms.append(g)
  df = pd.DataFrame({"cell_id": ids, "geometry": geoms})
  df = df[df["cell_id"].isin(keep_ids)]
  df["geometry"] = [to_polygon(g) for g in df["geometry"]]
  return df
