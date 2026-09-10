"""Half-day step handling and the 'which snapshot precedes onset' rule.

FEDS bins observations into two VIIRS overpass groups per day, labelled AM and PM.
These are *local solar* groups, not calendar halves: measured against real NFP
`datetime_utc`, AM pixels fall at local 0.5-7.0 h and PM at local 7.1-14.4 h, with
NO observation between local 14.4 h and 00.5 h the next day.

pyroCb onset peaks at local ~15.0 h - inside that blind window - so picking "the
snapshot before onset" needs care.  See `step_before`.
"""
from __future__ import annotations

import datetime as dt
import re

import numpy as np
import pandas as pd

from . import config as C

STEP_RE = re.compile(r"^(\d{4})(\d{2})(\d{2})(AM|PM)$")

# Nominal local-solar hour of each overpass group, used only by the Tier B fallback.
NOMINAL_LOCAL_H = {"AM": 2.0, "PM": 13.0}


def parse_step(step: str) -> tuple[dt.date, str]:
    m = STEP_RE.match(step)
    if not m:
        raise ValueError(f"bad step {step!r}")
    y, mo, d, ap = m.groups()
    return dt.date(int(y), int(mo), int(d)), ap


def step_sort_key(step: str) -> tuple:
    d, ap = parse_step(step)
    return (d, 0 if ap == "AM" else 1)


def make_step(d: dt.date, ampm: str) -> str:
    return f"{d:%Y%m%d}{ampm}"


def prev_step(step: str) -> str:
    d, ap = parse_step(step)
    return make_step(d, "AM") if ap == "PM" else make_step(d - dt.timedelta(days=1), "PM")


def next_step(step: str) -> str:
    d, ap = parse_step(step)
    return make_step(d, "PM") if ap == "AM" else make_step(d + dt.timedelta(days=1), "AM")


def norm_lon(lon):
    """Normalise longitude to (-180, 180]."""
    lon = np.asarray(lon, dtype=float)
    out = ((lon + 180.0) % 360.0) - 180.0
    return np.where(out == -180.0, 180.0, out)


def local_solar_hour(utc: pd.Series | pd.Timestamp, lon) -> np.ndarray:
    """Local solar hour (0-24) = UTC hour + lon/15."""
    t = pd.to_datetime(utc)
    hours = (t.dt.hour + t.dt.minute / 60 + t.dt.second / 3600
             if isinstance(t, pd.Series)
             else t.hour + t.minute / 60 + t.second / 3600)
    return np.asarray(hours + norm_lon(lon) / 15.0) % 24.0


def local_date(utc, lon):
    """Calendar date in local solar time."""
    t = pd.to_datetime(utc)
    off = pd.to_timedelta(norm_lon(lon) / 15.0, unit="h")
    return (t + off).date() if not isinstance(t, pd.Series) else (t + off).dt.date


def nominal_utc(step: str, lon: float) -> pd.Timestamp:
    """Approximate UTC instant of a step's overpass at a given longitude."""
    d, ap = parse_step(step)
    return (pd.Timestamp(d) + pd.Timedelta(hours=NOMINAL_LOCAL_H[ap])
            - pd.Timedelta(hours=float(norm_lon(lon)) / 15.0))


def step_before_local_clock(onset_utc: pd.Timestamp, lon: float) -> str:
    """Tier B: local-solar rule, used when a fire has no NFP observations.

    h < 7.0            -> previous local date, PM
    7.0 <= h < 14.4    -> local date, AM   (that day's PM may straddle onset)
    h >= 14.4          -> local date, PM
    """
    h = float(local_solar_hour(pd.Series([onset_utc]), [lon])[0])
    ld = local_date(onset_utc, lon)
    if h < C.LOCAL_AM_END_H:
        return make_step(ld - dt.timedelta(days=1), "PM")
    if h < C.LOCAL_PM_END_H:
        return make_step(ld, "AM")
    return make_step(ld, "PM")


def step_containing(onset_utc: pd.Timestamp, lon: float) -> str:
    """The overpass group that best represents the fire's state AT onset.

    Used for *matching* (which fire produced this pyroCb), where we want the fire as
    it existed that day - unlike `step_before`, which is used for the pre-band and
    deliberately excludes the pyroCb-day run.

    h < 7.0  -> that day's AM (the overpass just gone)
    h >= 7.0 -> that day's PM.  Above local 14.4 h we are in the observation blind
                window, so the same day's PM is still the most recent look at the fire.
    """
    h = float(local_solar_hour(pd.Series([onset_utc]), [lon])[0])
    ld = local_date(onset_utc, lon)
    return make_step(ld, "AM" if h < C.LOCAL_AM_END_H else "PM")


def snapshot_offsets_h(onset_utc: pd.Timestamp, lon: float) -> tuple[float, float, float, float]:
    """(b1_h, b2_h, a1_h, a2_h): nominal hours from onset to the 1st/2nd overpass
    step before onset and the 1st/2nd strictly after, on the fixed AM/PM schedule -
    independent of whether that step's file or detections actually exist (contrast
    `pyrocb/nearby.py`'s `_walk_available`, which is this same walk gated on data
    availability). b1/b2 are typically negative, a1/a2 positive, but see below.

    Anchor is `step_containing`, and b1 IS the anchor (not `prev_step` of it) - the
    same "anchor counts as the nearest before" convention `nearby.py`'s "before"
    tiers use, so a map feature built on this and v2's `before`-method both mean the
    same overpass by "1st before" for the same event. Most onsets fall after both
    that day's AM and PM windows close (module docstring: onset peaks around local
    15.0 h, past the 14.4 h blind-window boundary), so b1's NOMINAL time usually
    precedes onset - but for an onset early in the PM window (local 7.1-13.0 h) b1's
    nominal representative time (13.0 h) can fall slightly AFTER the actual onset.
    That is the same coarse half-day approximation `nominal_utc`/`_dt_hours` already
    carry everywhere else in this codebase, not a new one introduced here.
    """
    anchor = step_containing(onset_utc, lon)
    b1, b2 = anchor, prev_step(anchor)
    a1 = next_step(anchor)
    a2 = next_step(a1)
    to_h = lambda s: (nominal_utc(s, lon) - onset_utc).total_seconds() / 3600.0
    return to_h(b1), to_h(b2), to_h(a1), to_h(a2)


def step_before(onset_utc: pd.Timestamp, obs: pd.DataFrame | None, lon: float,
                available: set[str] | None = None) -> tuple[str | None, str, float]:
    """Return (step, tier, obs_lag_h) for the last observation before `onset_utc`.

    Tier A (preferred): the latest step whose MEDIAN observation time for this fire
    precedes onset.  Uses real timestamps, so no local-time model is needed, and one
    stray late pixel cannot disqualify an overpass whose bulk clearly precedes onset.
    `obs` must have columns step, t_med (for one mergid).

    Tier B: local-solar clock, when the fire has no usable observations.
    Tier C: walk backwards if the chosen step has no snapshot file.
    """
    if obs is not None and len(obs):
        ok = obs[obs["t_med"] < onset_utc]
        if len(ok):
            row = ok.loc[ok["t_med"].idxmax()]
            step = str(row["step"])
            lag = (onset_utc - row["t_med"]).total_seconds() / 3600.0
            step, back = _walk_back(step, available)
            return step, ("A" if not back else "A+C"), lag

    step = step_before_local_clock(onset_utc, lon)
    step, back = _walk_back(step, available)
    lag = (onset_utc - nominal_utc(step, lon)).total_seconds() / 3600.0 if step else np.nan
    return step, ("B" if not back else "B+C"), lag


def _walk_back(step: str | None, available: set[str] | None) -> tuple[str | None, bool]:
    """Tier C: step back until a step that actually exists on disk. Gaps are real."""
    if step is None or available is None or step in available:
        return step, False
    for _ in range(C.MAX_BACKWALK_STEPS):
        step = prev_step(step)
        if step in available:
            return step, True
    return None, True
