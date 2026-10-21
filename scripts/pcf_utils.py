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
