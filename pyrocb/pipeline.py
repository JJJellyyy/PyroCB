"""Orchestration: run the per-year stages with caching, then export.

Every stage is per-year and independent, so switching SNAPSHOT_YEARS from [2023] to
None (all years) reuses 2023 untouched and only computes the rest.
"""
from __future__ import annotations

import glob
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from . import active as A
from . import bands as B
from . import config as C
from . import io
from . import matching as M
from . import variables as V


def _inputs(year: int) -> list:
    """Files whose mtime/size invalidate a year's cached stages."""
    snaps = sorted(glob.glob(os.path.join(C.SNAP_DIR.format(y=year), "*.gpkg")))
    return [C.XLSX_FP, C.PERIMS_TPL.format(y=year), C.IGN_TPL.format(y=year)] + snaps[:1] + snaps[-1:]


def match_cached(year: int) -> pd.DataFrame:
    params = {"cap": C.MATCH_CAP_KM, "crs": C.METRIC_CRS, "am": C.LOCAL_AM_END_H}
    df = C.cached("match", year, params, _inputs(year),
                  lambda: _norm(M.match_year(year)))
    return _denorm(df)


def bands_cached(year: int, matches: pd.DataFrame):
    """Geometry and per-event bookkeeping come from one build; memoise so a cold
    cache does not run the (expensive) band construction twice."""
    params = {"prec": C.PRECISION_M, "sliver": C.SLIVER_KM2, "cap": C.MATCH_CAP_KM}
    inputs = _inputs(year)
    held = {}

    def both():
        if "v" not in held:
            held["v"] = B.build_bands(year, matches)
        return held["v"]

    g = C.cached("bands", year, params, inputs, lambda: both()[0], geo=True)
    e = C.cached("bandev", year, params, inputs, lambda: both()[1])
    return g, e


def events_cached(year: int, matches, band_ev) -> pd.DataFrame:
    params = {"win": list(C.NFP_WINDOWS_H), "cap": C.MATCH_CAP_KM, "prec": C.PRECISION_M}
    df = C.cached("events", year, params, _inputs(year),
                  lambda: _norm(V.build_events(year, matches, band_ev)))
    return _denorm(df)


def _norm(df: pd.DataFrame) -> pd.DataFrame:
    """Make a frame parquet-safe: nullable ints, string objects."""
    df = df.copy()
    if "mergid" in df:
        df["mergid"] = pd.to_numeric(df["mergid"], errors="coerce").astype("Int64")
    if "nearest_mergid" in df:
        df["nearest_mergid"] = pd.to_numeric(df["nearest_mergid"], errors="coerce").astype("Int64")
    for c in df.columns:
        if df[c].dtype != object:
            continue
        vals = df[c].dropna()
        # bool-valued object columns must not become the strings "True"/"False"
        if len(vals) and vals.map(lambda v: isinstance(v, (bool, np.bool_))).all():
            df[c] = df[c].astype("boolean")
        else:
            df[c] = df[c].astype("string")
    return df


def _denorm(df: pd.DataFrame) -> pd.DataFrame:
    return df


def active_cached(year: int, events: pd.DataFrame):
    params = {"win": list(C.ACTIVE_WINDOWS_H), "cap": C.MATCH_CAP_KM, "prec": C.PRECISION_M}
    return C.cached("active", year, params, _inputs(year),
                    lambda: A.active_polygons(year, events), geo=True)


def run_year(year: int, verbose: bool = True):
    if verbose:
        print(f"[{year}]", flush=True)
    m = match_cached(year)
    g, be = bands_cached(year, m)
    ev = events_cached(year, m, be)
    act = active_cached(year, ev)
    # fold the active-burn areas back onto the event rows
    s = A.summarise(act)
    if len(s):
        ev = ev.merge(s, on="pcb_id", how="left")
    return m, g, ev, act


def run(years: list[int], verbose: bool = True):
    """Run every requested year; returns (events, bands, active) concatenated."""
    evs, gs, acts = [], [], []
    for y in years:
        _, g, ev, act = run_year(y, verbose)
        evs.append(ev)
        if len(g):
            gs.append(g)
        if len(act):
            acts.append(act)
    events = pd.concat(evs, ignore_index=True) if evs else pd.DataFrame()
    bands = pd.concat(gs, ignore_index=True) if gs else B._empty_bands()
    active = pd.concat(acts, ignore_index=True) if acts else A._empty()
    bands = gpd.GeoDataFrame(bands, geometry="geometry", crs=C.METRIC_CRS) if len(bands) else bands
    active = gpd.GeoDataFrame(active, geometry="geometry", crs=C.METRIC_CRS) if len(active) else active
    return events, bands, active


# ------------------------------------------------------------------ export
DICT = [
    ("pcb_id", "-", "Stable event id PCB_YYYY_MMDD_HHMM (UTC); a/b/c suffix where a timestamp carries several events", "PCBAtlas.xlsx"),
    ("mergid", "-", "FireAtlas post-merge fire id; the join key across all products", "snapshot/final/ignitions/NFP"),
    ("match_method", "-", "inside | nearest | inside_dormant | nearest_dormant | unmatched", "derived"),
    ("match_dist_km", "km", "Distance from the pyroCb point to the matched fire; 0 when inside", "derived"),
    ("match_step", "-", "Half-day snapshot used for matching", "derived"),
    ("step_before", "-", "Last half-day step whose median observation precedes onset", "derived"),
    ("step_tier", "-", "A = NFP observation times, B = local-solar clock, +C = walked back to an existing file", "derived"),
    ("obs_lag_h", "h", "Hours from that step's median observation to onset. LARGE VALUES MEAN A STALE PRE-BAND", "derived"),
    ("area_pre_km2", "km2", "Fire area already burned as of step_before (lower bound; see obs_lag_h)", "snapshot"),
    ("area_post_km2", "km2", "Final area minus area_pre", "derived"),
    ("area_final_km2", "km2", "Area of the cumulative union S_final", "final perimeter"),
    ("frac_burned_before", "-", "area_pre / area_final, in [0,1] by construction. BIASED LOW, correlated with time of day", "derived"),
    ("band_id", "-", "Index of the growth band whose inner boundary is this event's step_before", "derived"),
    ("n_events_in_band", "-", "Events sharing this band; >1 means step-derived columns are identical for them", "derived"),
    ("dedup_collapsed", "bool", "True when this event shares its half-day step with another on the same fire", "derived"),
    ("pcb_index", "-", "1-based rank of this event among its fire's pyroCb", "derived"),
    ("n_pcb_total", "-", "pyroCb events matched to this fire over the season", "derived"),
    ("n_pcb_before", "-", "pyroCb on this fire strictly before this one", "derived"),
    ("n_pcb_after", "-", "pyroCb on this fire strictly after this one", "derived"),
    ("hours_since_prev_pcb", "h", "Gap to the previous pyroCb on the same fire", "derived"),
    ("hours_to_next_pcb", "h", "Gap to the next pyroCb on the same fire", "derived"),
    ("n_ign_total", "-", "Ignitions in this fire's merge complex (mergid group)", "ignitions"),
    ("n_ign_before", "-", "Of those, ignited before onset", "ignitions"),
    ("n_ign_after", "-", "Of those, ignited after onset", "ignitions"),
    ("n_ign_same_halfday", "-", "Ignitions stamped within +/-7 h of onset, i.e. unresolvable", "ignitions"),
    ("ign_merge_span_days", "d", "Span from first to last ignition of the complex. LARGE => the concentric band model is questionable", "ignitions"),
    ("n_pixels_before", "-", "VIIRS detections before onset (exact datetime_utc)", "NFP"),
    ("n_pixels_after", "-", "VIIRS detections at/after onset", "NFP"),
    ("n_pixels_total", "-", "Season total; equals the final perimeter n_pixels (asserted)", "NFP"),
    ("frp_p95_prev24h", "MW", "95th percentile FRP in the 24 h before onset. HEADLINE intensity variable", "NFP"),
    ("n_obs_prev24h", "-", "Distinct observation times in that window; 0 means no overpass, not no fire", "NFP"),
    ("qc_stale_obs", "bool", "No overpass within 14 h of onset, so every 'before' quantity is stale", "derived"),
    ("qc_final_smaller", "bool", "Final perimeter smaller than the last snapshot (expected: never)", "derived"),
    ("snap_farea", "km2", "Fire area in the snapshot at match_step (cumulative burned, not actively burning)", "snapshot"),
    ("snap_flinelen", "km", "Active fire-line length at match_step", "snapshot"),
    ("snap_meanFRP", "W/m2", "Area-normalised mean FRP at match_step", "snapshot"),
    ("snap_n_components", "-", "Snapshot rows sharing this mergid (un-executed merge components)", "snapshot"),
    ("active_area_km2_12h", "km2", "ACTIVE BURN AREA: union of VIIRS pixel footprints within +/-12 h of onset. The sampling geometry for land cover / ERA5", "NFP"),
    ("active_area_km2_6h / _24h", "km2", "Same at +/-6 h and +/-24 h, for sensitivity", "NFP"),
    ("active_n_pixels_12h", "-", "Detections inside that window", "NFP"),
    ("active_n_obs_12h", "-", "Distinct overpass times in the window; 1 means a single look", "NFP"),
    ("active_n_before_onset_12h", "-", "Of those detections, how many precede onset", "NFP"),
    ("active_frp_p95_12h", "MW", "95th percentile FRP of the active-burn pixels", "NFP"),
    ("match_fwd_steps", "-", "Overpasses searched AHEAD to find a fire that blew up before it was first observed; 0 normally", "derived"),
    ("fire_farea_km2", "km2", "Season-final fire area", "final perimeter"),
    ("fire_duration_days", "d", "Fire duration", "final perimeter"),
    ("fire_lcc_final", "-", "Dominant land cover of the final perimeter", "final perimeter"),
    ("lat / lon", "deg", "APPROXIMATE FIRE location (not the plume), stored to 0.1 deg ~ 11 km", "PCBAtlas.xlsx"),
    ("max_inject_alt", "m", "Maximum smoke injection altitude", "PCBAtlas.xlsx"),
    ("trop_relative_alt", "m", "Injection altitude minus MERRA-2 tropopause; negative = below", "PCBAtlas.xlsx"),
]


def export(events: pd.DataFrame, active=None, verbose: bool = True):
    df = events.copy()
    for c in df.columns:
        if pd.api.types.is_datetime64_any_dtype(df[c]):
            df[c] = df[c].dt.strftime("%Y-%m-%d %H:%M")
    df.to_csv(C.OUT_TSV, sep="\t", index=False, encoding="utf-8")
    try:
        df.to_excel(C.OUT_XLSX, index=False)
    except Exception as e:  # openpyxl hiccup shouldn't kill the run
        print(f"  (xlsx skipped: {e})", flush=True)
    pd.DataFrame(DICT, columns=["column", "units", "definition", "source"]).to_csv(
        C.OUT_DICT, sep="\t", index=False, encoding="utf-8")
    if verbose:
        print(f"  wrote {C.OUT_TSV.name} ({len(df)} rows x {df.shape[1]} cols), "
              f"{C.OUT_XLSX.name}, {C.OUT_DICT.name}", flush=True)
    return df
