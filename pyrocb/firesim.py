"""Whole-fire evolution series for the map's "animate the entire event" mode.

Every NFP detection of a pyroCb fire (the full merge complex, *not* localised - the
point is to watch the fronts grow and merge), from first detection to last, snapped to
a 375 m grid and reduced to **the time each cell first burned** (`burngrid.encode`,
`collapse_first=True`). That is exactly what a cumulative growth animation needs, and
it makes "newly burned in the last W hours" a contiguous slice of the same array.

This replaced a tiered decimation that kept detections near a pyroCb onset at native
resolution and thinned the quiet spans to 1 km / 3 km, inflating each surviving
detection to its dedup cell. On a small fire a single 3 km block covered most of the
perimeter. The grid is both finer (375 m everywhere) and smaller: 0.64 MB vs 1.7 MB
base64 for 2023.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import burngrid
from . import config as C
from . import io


def fire_series(year: int, onsets: dict) -> dict:
    """{str(mergid): payload} for every mergid in `onsets`.

    `onsets[mergid]` is a list of {"t": UTC, "lat": ..., "lon": ..., "id": pcb_id}.
    """
    nfp = io.read_nfp_year(year)
    by = {k: v for k, v in nfp.groupby("mergid")}
    out = {}

    for mergid, ons in onsets.items():
        g = by.get(int(mergid))
        if g is None or not len(g):
            continue
        t = pd.to_datetime(g["datetime_utc"]).to_numpy()
        t0 = t.min()
        mins = (t - t0) / np.timedelta64(1, "m")

        grid = burngrid.encode(g["x"].to_numpy(float), g["y"].to_numpy(float), mins,
                               collapse_first=True)
        if grid is None:
            continue

        ons = sorted(ons, key=lambda o: pd.Timestamp(o["t"]))
        ons_h = [round(float((np.datetime64(pd.Timestamp(o["t"])) - t0)
                             / np.timedelta64(1, "m")) / 60.0, 2) for o in ons]
        out[str(int(mergid))] = {
            "t0": pd.Timestamp(t0).strftime("%Y-%m-%d %H:%M"),
            "dur_h": round(float(mins.max()) / 60.0, 1),
            "pcb_h": ons_h,
            "pcb_ll": [[round(float(o["lat"]), 4), round(float(o["lon"]), 4)] for o in ons],
            "pcb_id": [str(o["id"]) for o in ons],
            "g": grid,
            "n_raw": int(len(g)),
            "n_cell": burngrid.n_cells(grid),
        }
    return out
