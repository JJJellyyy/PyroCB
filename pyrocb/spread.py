"""Per-fire step-attribute time series, from the 12-hourly Snapshot attributes.

Each Snapshot step carries `avgspread`/`spread95` (FireAtlas' mean / 95th-percentile
fire-line advance for that step, metres) and `meanFRP`/`FRP95` (fire radiative power:
area-normalised mean in W/m2, 95th-percentile in MW - FireAtlas ships these on
different bases, which is a known asymmetry in the product, not a bug here). They
ship as strings and are 0 on inactive steps (F7). This module reads them once per
year, aggregates the merge components per `mergid` the same way `io.read_snapshot`
does (mean for the area-normalised/average fields, max for the percentile fields),
and returns a compact per-fire series the map colours the whole-fire animation by -
by spread rate or by radiative power, the viewer's choice.

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


_COLS = ["avgspread", "spread95", "meanFRP", "FRP95"]
_KEYS = {"avgspread": "av", "spread95": "s95", "meanFRP": "frp", "FRP95": "frp95"}
_AGG  = {"av": "mean", "s95": "max", "frp": "mean", "frp95": "max"}


def fire_spread_steps(year: int, mergids) -> dict:
    """{str(mergid): {"t": [iso...], "av": [...], "s95": [...], "frp": [...], "frp95": [...]}}.

    Only steps where the fire is present are kept; the series is sorted by time.
    """
    want = {int(m) for m in mergids}
    acc: dict[int, list] = {}
    for step in io.list_steps(year):
        fp = os.path.join(C.SNAP_DIR.format(y=year), f"{step}.gpkg")
        try:
            g = gpd.read_file(fp, columns=["mergid"] + _COLS, ignore_geometry=True)
        except Exception:
            continue
        g = g[g["mergid"].isin(want)]
        if not len(g):
            continue
        for c in _COLS:
            g[c] = pd.to_numeric(g[c], errors="coerce").fillna(0.0)
        a = g.groupby("mergid").agg(**{_KEYS[c]: (c, _AGG[_KEYS[c]]) for c in _COLS})
        iso = _step_end_utc(step).strftime("%Y-%m-%dT%H:%M")
        for mid, row in a.iterrows():
            vals = tuple(0.0 if pd.isna(row[k]) else round(float(row[k]), 1)
                        for k in ("av", "s95", "frp", "frp95"))
            acc.setdefault(int(mid), []).append((iso, *vals))

    out = {}
    for mid, rows in acc.items():
        rows.sort()
        out[str(mid)] = {"t":    [r[0] for r in rows],
                         "av":   [r[1] for r in rows],
                         "s95":  [r[2] for r in rows],
                         "frp":  [r[3] for r in rows],
                         "frp95":[r[4] for r in rows]}
    return out


def clamp(spread_by_year: dict) -> dict:
    """Global 2nd..98th percentile of the non-zero values, for the colour ramp."""
    import numpy as np
    acc = {"av": [], "s95": [], "frp": [], "frp95": []}
    for ydict in spread_by_year.values():
        for s in ydict.values():
            for k in acc:
                acc[k] += [v for v in s[k] if v > 0]
    def rng(a):
        a = np.asarray(a, float)
        return [0.0, 1.0] if not len(a) else [float(np.percentile(a, 2)),
                                              float(np.percentile(a, 98))]
    return {k: rng(v) for k, v in acc.items()}
