"""Per-fire spread-rate time series, from the 12-hourly Snapshot attributes.

Each Snapshot step carries `avgspread` and `spread95` per fire (FireAtlas'
mean / 95th-percentile fire-line advance for that step, metres). They ship as
strings and are 0 on inactive steps (F7). This module reads them once per year,
aggregates the merge components per `mergid` the same way `io.read_snapshot`
does (mean of avgspread, max of spread95), and returns a compact per-fire series
the map colours the whole-fire animation by.

Step time is the **nominal end of the 12 h interval** in UTC: AM -> that day
12:00, PM -> next 00:00. Monotonic and 12 h-spaced, which is exactly the real
observation cadence the animation scrubber snaps to.
"""
from __future__ import annotations

import os

import geopandas as gpd
import pandas as pd

from . import config as C
from . import io


def _step_end_utc(step: str) -> pd.Timestamp:
    d = pd.Timestamp(f"{step[:4]}-{step[4:6]}-{step[6:8]}")
    return d + pd.Timedelta(hours=12 if step[8:10] == "AM" else 24)


def fire_spread_steps(year: int, mergids) -> dict:
    """{str(mergid): {"t": [iso...], "av": [float...], "s95": [float...]}}.

    Only steps where the fire is present are kept; the series is sorted by time.
    """
    want = {int(m) for m in mergids}
    acc: dict[int, list] = {}
    for step in io.list_steps(year):
        fp = os.path.join(C.SNAP_DIR.format(y=year), f"{step}.gpkg")
        try:
            g = gpd.read_file(fp, columns=["mergid", "avgspread", "spread95"],
                              ignore_geometry=True)
        except Exception:
            continue
        g = g[g["mergid"].isin(want)]
        if not len(g):
            continue
        g["avgspread"] = pd.to_numeric(g["avgspread"], errors="coerce").fillna(0.0)
        g["spread95"] = pd.to_numeric(g["spread95"], errors="coerce").fillna(0.0)
        a = g.groupby("mergid").agg(av=("avgspread", "mean"), s95=("spread95", "max"))
        iso = _step_end_utc(step).strftime("%Y-%m-%dT%H:%M")
        for mid, row in a.iterrows():
            av = 0.0 if pd.isna(row["av"]) else round(float(row["av"]), 1)
            s95 = 0.0 if pd.isna(row["s95"]) else round(float(row["s95"]), 1)
            acc.setdefault(int(mid), []).append((iso, av, s95))

    out = {}
    for mid, rows in acc.items():
        rows.sort()
        out[str(mid)] = {"t": [r[0] for r in rows],
                         "av": [r[1] for r in rows],
                         "s95": [r[2] for r in rows]}
    return out


def clamp(spread_by_year: dict) -> dict:
    """Global 2nd..98th percentile of the non-zero values, for the colour ramp."""
    import numpy as np
    av, s95 = [], []
    for ydict in spread_by_year.values():
        for s in ydict.values():
            av += [v for v in s["av"] if v > 0]
            s95 += [v for v in s["s95"] if v > 0]
    def rng(a):
        a = np.asarray(a, float)
        return [0.0, 1.0] if not len(a) else [float(np.percentile(a, 2)),
                                              float(np.percentile(a, 98))]
    return {"av": rng(av), "s95": rng(s95)}
