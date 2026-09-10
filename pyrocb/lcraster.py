"""Land cover as paletted PNG overlays, one per fire.

The polygon route hit a wall (D-39): true 30 m vector land cover over every pyroCb
perimeter is ~14 MB compressed for 2023 alone, and even the 1 km tier costs thousands
of canvas paths to pan and zoom. A raster does the same job better on every axis:

* **Size.** An 8-bit paletted PNG of a 15-class map is dominated by large uniform
  regions and compresses ~50:1.
* **Speed.** `L.imageOverlay` is one `<img>` per fire, composited by the GPU. Panning
  and zooming stay smooth no matter how detailed the land cover is; the old polygon
  layer had to re-stroke every path.
* **Build time.** No `rasterio.features.shapes`, no dissolve, no simplify - the slow,
  topology-exception-prone part of the old path disappears.

**Grid alignment (D-64).** `LC_PNG_CELL_M` is deliberately equal to `BURN_CELL_M` (375 m,
one VIIRS pixel), and the output grid is snapped to the same convention
`burngrid.py` uses: cell edges fall on exact multiples of the cell size measured from
the EPSG:3857 Mercator origin, not from each fire's own bounding box. Two rasters built
this way - land cover here, fire activity in `burngrid.py` - always share cell
boundaries at a given fire, because both grids are ultimately the same infinite tiling
of Mercator space, just cropped to different windows. The old version resampled in the
source CRS at each fire's own bbox-relative origin, then reprojected with an
auto-computed transform - two fires' land-cover grids never lined up with anything,
including the fire-activity layer, even at matching resolutions.

The image is warped to **EPSG:3857** so Leaflet's linear placement of an image over a
lat/lon rectangle is exact rather than stretched (`L.imageOverlay` maps the image
linearly in *display* projection space) - already true before D-64, kept the same.

A 3x3 majority filter is applied before encoding: at map scale it is invisible, and it
removes the isolated-pixel speckle that PNG spends bytes on.
"""
from __future__ import annotations

import base64
import io as _io
import math

import numpy as np
import rasterio
from PIL import Image
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform, transform_bounds
from rasterio.windows import from_bounds
from scipy import ndimage

from . import config as C
from .landcover import LC_CLASSES, _tif

_SRC_CRS = 3979
_DST_CRS = 3857


def _palette() -> bytes:
    pal = np.zeros((256, 3), "uint8")
    for code, (_label, hexc) in LC_CLASSES.items():
        pal[code] = [int(hexc[k:k + 2], 16) for k in (1, 3, 5)]
    return pal.tobytes()


def _despeckle(a: np.ndarray) -> np.ndarray:
    """3x3 majority filter. Invisible at map scale, trims PNG bytes."""
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


def _bounds_wgs84(tr, w, h):
    """(south, west), (north, east) of a 3857 raster, for L.imageOverlay."""
    x0, y0 = tr * (0, 0)
    x1, y1 = tr * (w, h)
    lons, lats = transform(_DST_CRS, C.WGS84, [x0, x1], [y0, y1])
    return [[min(lats), min(lons)], [max(lats), max(lons)]]


def _mercator_cell_m(lat_deg: float, ground_m: float) -> float:
    """Ground metres -> EPSG:3857 metres at this latitude - the same Mercator
    inflation `burngrid.encode` applies, so a cell here is `ground_m` on the actual
    ground, not `ground_m` of already-distorted map distance."""
    return ground_m / max(0.15, math.cos(math.radians(lat_deg)))


def _snapped_transform(bx0: float, by0: float, bx1: float, by1: float, cell: float):
    """A 3857 transform whose cell edges are exact multiples of `cell` from the
    Mercator origin (0, 0) - `burngrid.encode`'s own `floor(mx / cell)` convention.
    Any two rasters gridded this way with the same `cell` share cell boundaries
    exactly, independent of which window of the world either one covers."""
    i0, i1 = math.floor(bx0 / cell), math.ceil(bx1 / cell)
    j0, j1 = math.floor(by0 / cell), math.ceil(by1 / cell)
    w, h = max(1, i1 - i0), max(1, j1 - j0)
    return rasterio.transform.from_origin(i0 * cell, j1 * cell, cell, cell), w, h


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

            cen = bg.centroid
            _, clat = transform(_SRC_CRS, C.WGS84, [cen.x], [cen.y])
            bx0, by0, bx1, by1 = transform_bounds(_SRC_CRS, _DST_CRS, minx, miny, maxx, maxy)

            cell = _mercator_cell_m(clat[0], C.LC_PNG_CELL_M)
            span = max(bx1 - bx0, by1 - by0)
            if span / cell > C.LC_PNG_CAP_PX:      # continent-sized fire: widen the
                cell = span / C.LC_PNG_CAP_PX       # cell, keep the same grid origin
            dst_transform, dw, dh = _snapped_transform(bx0, by0, bx1, by1, cell)

            # Source read: decimated to roughly 3x the destination grid, not to native
            # 30 m - the destination is already capped (LC_PNG_CAP_PX), so this stays
            # cheap even for a huge fire, while still giving the majority resample
            # below several source samples per destination cell.
            oversample = 3
            rw, rh = max(1, dw * oversample), max(1, dh * oversample)
            src = r.read(1, window=from_bounds(minx, miny, maxx, maxy, r.transform),
                         out_shape=(rh, rw), resampling=Resampling.mode,
                         boundless=True, fill_value=0).astype("uint8")
            if not src.any():
                continue                           # outside the Canadian raster (Siberia)
            src_transform = rasterio.transform.from_bounds(minx, miny, maxx, maxy, rw, rh)

            dst = np.zeros((dh, dw), "uint8")
            reproject(src, dst, src_transform=src_transform, src_crs=_SRC_CRS,
                     dst_transform=dst_transform, dst_crs=_DST_CRS,
                     resampling=Resampling.mode)
            dst[~np.isin(dst, list(LC_CLASSES))] = 0
            if not dst.any():
                continue
            dst = _despeckle(dst)

            im = Image.fromarray(dst)
            im.putpalette(pal)
            buf = _io.BytesIO()
            im.save(buf, "PNG", optimize=True, transparency=0, bits=8)
            out[str(int(mid))] = {
                "b": _bounds_wgs84(dst_transform, dw, dh),
                "d": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(),
            }
    return out
