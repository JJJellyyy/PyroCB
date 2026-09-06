"""Gridded burn-time encoding - the map's fire-activity payload.

Replaces "ship every VIIRS detection as a lat/lon box". Detections are snapped to a
fixed ground grid (`BURN_CELL_M`, 375 m = one VIIRS pixel) built in **EPSG:3857**, so
the browser can turn a cell index straight into a canvas pixel with one multiply-add
(Web Mercator is the display projection) instead of projecting a polygon per feature.

Two collapse modes, both producing the same structure:

* `collapse_first=True`  - one entry per cell, the time it FIRST burned. What a
  cumulative whole-fire animation needs; also makes "newly burned in the last W hours"
  a contiguous slice, which is the honest way to draw a spreading front.
* `collapse_first=False` - one entry per (cell, observation time). Exact, used for the
  +/-24 h window view where a cell genuinely re-burns across overpasses.

Entries are sorted by time and the three index arrays are delta-encoded, so gzip sees
long runs of tiny integers. Measured on 2023: 393 k cells at true 375 m compress to
**0.64 MB** base64, against 1.7 MB for the old 275 k decimated points at 1-3 km blocks
- 2.7x smaller and ~8x finer.

Because entries are time-sorted, EVERY view the map offers (cumulative up to t, a
sliding window, a static +/-W window) is a contiguous index range, which is what lets
the canvas layer draw incrementally instead of redrawing every frame.
"""
from __future__ import annotations

import numpy as np
from pyproj import Transformer

from . import config as C

R_MERC = 6378137.0
_TF = None


def _tf() -> Transformer:
    global _TF
    if _TF is None:
        _TF = Transformer.from_crs(C.METRIC_CRS, 3857, always_xy=True)
    return _TF


def encode(x, y, t_min, cell_ground_m: float | None = None,
           collapse_first: bool = True) -> dict | None:
    """Encode detections into a delta-coded, time-sorted cell grid.

    `x`, `y` are EPSG:3571 metres; `t_min` is time in **minutes** relative to whatever
    reference the caller wants (onset, or the fire's first detection).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    t = np.rint(np.asarray(t_min, float)).astype(np.int64)
    if not len(x):
        return None

    mx, my = _tf().transform(x, y)
    mx = np.asarray(mx, float)
    my = np.asarray(my, float)
    ok = np.isfinite(mx) & np.isfinite(my)
    if not ok.any():
        return None
    mx, my, t = mx[ok], my[ok], t[ok]

    # Mercator metres are inflated by 1/cos(lat); size the cell so it is
    # `cell_ground_m` on the ground at this fire's latitude.
    lat = np.degrees(2.0 * np.arctan(np.exp(my.mean() / R_MERC)) - np.pi / 2.0)
    cell = float(cell_ground_m or C.BURN_CELL_M) / max(0.15, np.cos(np.radians(lat)))

    i = np.floor(mx / cell).astype(np.int64)
    j = np.floor(my / cell).astype(np.int64)
    i0, j0 = int(i.min()), int(j.min())
    i -= i0
    j -= j0

    if collapse_first:
        key = i * (np.int64(1) << 32) + j
        o = np.lexsort((t, key))
        ks, ts, is_, js_ = key[o], t[o], i[o], j[o]
        keep = np.r_[True, ks[1:] != ks[:-1]]
        i, j, t = is_[keep], js_[keep], ts[keep]
    else:
        trip = np.unique(np.stack([i, j, t], axis=1), axis=0)
        i, j, t = trip[:, 0], trip[:, 1], trip[:, 2]

    o = np.argsort(t, kind="stable")
    i, j, t = i[o], j[o], t[o]

    return {
        "ox": round(i0 * cell, 1), "oy": round(j0 * cell, 1), "c": round(cell, 3),
        "di": np.diff(i, prepend=0).tolist(),
        "dj": np.diff(j, prepend=0).tolist(),
        "dt": np.diff(t, prepend=0).tolist(),
    }


def n_cells(payload: dict | None) -> int:
    return 0 if not payload else len(payload["di"])
