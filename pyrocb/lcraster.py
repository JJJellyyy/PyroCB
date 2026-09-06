"""Land cover as paletted PNG overlays, one per fire.

The polygon route hit a wall (D-39): true 30 m vector land cover over every pyroCb
perimeter is ~14 MB compressed for 2023 alone, and even the 1 km tier costs thousands
of canvas paths to pan and zoom. A raster does the same job better on every axis:

* **Size.** An 8-bit paletted PNG of a 15-class map is dominated by large uniform
  regions and compresses ~50:1. At ~120 m ground resolution 2023 comes to ~1.9 MB for
  all 69 fires - roughly 8x finer than the old 1 km tier for less than 2x the bytes.
* **Speed.** `L.imageOverlay` is one `<img>` per fire, composited by the GPU. Panning
  and zooming stay smooth no matter how detailed the land cover is; the old polygon
  layer had to re-stroke every path.
* **Build time.** No `rasterio.features.shapes`, no dissolve, no simplify - the slow,
  topology-exception-prone part of the old path disappears.

The image is warped to **EPSG:3857** so Leaflet's linear placement of an image over a
lat/lon rectangle is exact rather than stretched (`L.imageOverlay` maps the image
linearly in *display* projection space).

A 3x3 majority filter is applied before encoding: at map scale it is invisible, and it
removes the isolated-pixel speckle that PNG spends ~35% of its bytes on.
"""
from __future__ import annotations

import base64
import io as _io

import numpy as np
import rasterio
from PIL import Image
from rasterio.warp import Resampling, calculate_default_transform, reproject
from rasterio.windows import from_bounds
from scipy import ndimage

from . import config as C
from .landcover import LC_CLASSES, _tif


def _palette() -> bytes:
    pal = np.zeros((256, 3), "uint8")
    for code, (_label, hexc) in LC_CLASSES.items():
        pal[code] = [int(hexc[k:k + 2], 16) for k in (1, 3, 5)]
    return pal.tobytes()


def _despeckle(a: np.ndarray) -> np.ndarray:
    """3x3 majority filter. Invisible at map scale, ~35% off the PNG."""
    out = np.zeros_like(a)
    best = np.zeros(a.shape, "float32")
    for v in np.unique(a):
        if v == 0:
            continue
        c = ndimage.uniform_filter((a == v).astype("float32"), 3)
        m = c > best
        best = np.where(m, c, best)
        out = np.where(m, v, out)
    return out.astype("uint8")


def _bounds_wgs84(transform, w, h):
    """(south, west), (north, east) of a 3857 raster, for L.imageOverlay."""
    import pyproj
    x0, y0 = transform * (0, 0)
    x1, y1 = transform * (w, h)
    tf = pyproj.Transformer.from_crs(3857, C.WGS84, always_xy=True)
    lon0, lat0 = tf.transform(x0, y0)
    lon1, lat1 = tf.transform(x1, y1)
    return [[min(lat0, lat1), min(lon0, lon1)], [max(lat0, lat1), max(lon0, lon1)]]


def fire_overlays(year: int, mergids, geoms_3979) -> dict:
    """{str(mergid): {"b": [[s,w],[n,e]], "d": "data:image/png;base64,..."}}"""
    pal = _palette()
    out: dict[str, dict] = {}
    with rasterio.open(_tif(year)) as r:
        for mid, g in zip(mergids, geoms_3979):
            if g is None or g.is_empty:
                continue
            bg = g.buffer(C.LC_BUFFER_M)
            minx, miny, maxx, maxy = bg.bounds
            ext = max(maxx - minx, maxy - miny)
            cell = max(30.0, C.LC_PNG_MIN_M, ext / C.LC_PNG_CAP_PX)
            ow = max(1, int((maxx - minx) / cell))
            oh = max(1, int((maxy - miny) / cell))
            a = r.read(1, window=from_bounds(minx, miny, maxx, maxy, r.transform),
                       out_shape=(oh, ow), resampling=Resampling.mode,
                       boundless=True, fill_value=0).astype("uint8")
            a[~np.isin(a, list(LC_CLASSES))] = 0
            if not a.any():
                continue                       # outside the Canadian raster (Siberia)
            a = _despeckle(a)

            src_tr = rasterio.transform.from_bounds(minx, miny, maxx, maxy, ow, oh)
            dtr, dw, dh = calculate_default_transform(3979, 3857, ow, oh,
                                                      minx, miny, maxx, maxy)
            dst = np.zeros((dh, dw), "uint8")
            reproject(a, dst, src_transform=src_tr, src_crs=3979,
                      dst_transform=dtr, dst_crs=3857, resampling=Resampling.nearest)
            if not dst.any():
                continue

            im = Image.fromarray(dst)
            im.putpalette(pal)
            buf = _io.BytesIO()
            im.save(buf, "PNG", optimize=True, transparency=0, bits=8)
            out[str(int(mid))] = {
                "b": _bounds_wgs84(dtr, dw, dh),
                "d": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(),
            }
    return out
