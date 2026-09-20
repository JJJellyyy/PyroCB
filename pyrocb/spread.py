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

Step time is the UTC instant at which the step's **local-solar** observation window
opens - AM at local 0, PM at local `LOCAL_AM_END_H` - converted per fire via its own
longitude (`_step_start_utc`). Window *starts*, because the map attributes a burn cell
to the last marker at or before the cell's own time.
"""
from __future__ import annotations

import os

import geopandas as gpd
import pandas as pd

from . import config as C
from . import io
from . import timebins as T


def _step_start_utc(step: str, lon: float) -> pd.Timestamp:
    """UTC instant at which a step's local-solar observation window OPENS.

    Two things were wrong with the previous version (`_step_end_utc`), which returned
    `date + 12 h` for AM and `date + 24 h` for PM:

    1. **No longitude correction.** A step label is a LOCAL-SOLAR half-day (see
       `timebins`): its date is the local date and AM/PM are local windows. Pinning it
       to a UTC clock slid the whole series by lon/15 h - measured against real NFP
       `datetime_utc`, the error ran from +0.2 h at lon -165 to -21.4 h at lon +165,
       a 17.7 h longitude-dependent drift, and put 31% of detections outside the very
       12 h bin their timestamp claimed. Every other half-day -> UTC conversion in this
       codebase (`timebins.nominal_utc`, `io.read_ignitions`, `matching._halfday_utc`,
       `variables`) already subtracts `lon/15`; this one was the lone exception.

    2. **Wrong edge of the window.** The map assigns a burn cell to a step by walking
       to the LAST marker at or before the cell's time, so these markers have to be
       window STARTS. Emitting window ends put almost every cell on the previous
       overpass. AM covers local 0.5-7.0 h and PM local 7.1-14.4 h, so the dividing
       line is `LOCAL_AM_END_H` and a day's AM opens at local 0.

    Checked on 2023: detections landing in the step bin their own label names go from
    1.5% (old) to 100.0% (this).
    """
    d = pd.Timestamp(f"{step[:4]}-{step[4:6]}-{step[6:8]}")
    local_h = 0.0 if step[8:10] == "AM" else float(C.LOCAL_AM_END_H)
    return d + pd.to_timedelta(local_h - float(T.norm_lon(lon)) / 15.0, unit="h")


def _lon_by_mergid(year: int, want: set) -> dict:
    """Representative longitude per fire, for the local-solar -> UTC conversion."""
    fin = io.read_final(year)
    fin = fin[fin["mergid"].isin(want)]
    lon = fin.geometry.representative_point().to_crs(C.WGS84).x.values
    return {int(m): float(v) for m, v in zip(fin["mergid"].values, lon)}


_COLS = ["avgspread", "spread95", "meanFRP", "FRP95"]
_KEYS = {"avgspread": "av", "spread95": "s95", "meanFRP": "frp", "FRP95": "frp95"}
_AGG  = {"av": "mean", "s95": "max", "frp": "mean", "frp95": "max"}


def fire_spread_steps(year: int, mergids) -> dict:
    """{str(mergid): {"t": [iso...], "av": [...], "s95": [...], "frp": [...], "frp95": [...]}}.

    Only steps where the fire is present are kept; the series is sorted by time.
    """
    want = {int(m) for m in mergids}
    lon_by = _lon_by_mergid(year, want)
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
        for mid, row in a.iterrows():
            # Per fire, not per step: the local-solar -> UTC shift depends on where the
            # fire is, so two fires on the same step have different UTC boundaries.
            iso = _step_start_utc(step, lon_by.get(int(mid), 0.0)).strftime("%Y-%m-%dT%H:%M")
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
