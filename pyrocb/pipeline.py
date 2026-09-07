"""Orchestration: run the per-year stages with caching, then export.

Every stage is per-year and independent, so switching SNAPSHOT_YEARS from [2023] to
None (all years) reuses 2023 untouched and only computes the rest.
"""
from __future__ import annotations

import glob
import json
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from . import active as A
from . import bands as B
from . import config as C
from . import firesim as FS
from . import io
from . import landcover as LC
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


_SCHEMA_DEFAULTS = {"match_fwd_steps": 0, "nearest_mergid": pd.NA,
                    "qc_nfp_repaired": False, "qc_nfp_pixel_exact": True}


def _norm(df: pd.DataFrame) -> pd.DataFrame:
    """Make a frame parquet-safe: nullable ints, string objects, stable schema.

    Columns only created when a year actually has (say) a forward-matched event would
    otherwise make the per-year frames un-concatenable.
    """
    df = df.copy()
    for c, v in _SCHEMA_DEFAULTS.items():
        if c not in df.columns:
            df[c] = v
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
    params = {"win": [list(w) for w in C.ACTIVE_TABLE_WINDOWS], "cap": C.MATCH_CAP_KM,
              "prec": C.PRECISION_M, "bridge": C.ACTIVE_BRIDGE_M,
              "attach": C.ACTIVE_ATTACH_CAP_KM}
    return C.cached("active", year, params, _inputs(year),
                    lambda: A.active_polygons(year, events), geo=True)


def landcover_cached(year: int, mergids) -> "gpd.GeoDataFrame":
    """Coarse (10 km) land-cover polygon backdrop around the year's pyroCb fires.
    Everything closer in is a PNG overlay - see `lcpng_cached`."""
    mergids = sorted(int(m) for m in mergids)
    params = {"buf": C.LC_BUFFER_M, "bg_cell": C.LC_BG_CELL_M,
              "bg_buf": C.LC_BG_BUFFER_M, "bg_simp": C.LC_BG_SIMPLIFY_M,
              "split": C.LC_YEAR_SPLIT, "mergids": mergids, "v": 2}
    tif = LC._tif(year)

    def build():
        per = io.read_final(year).to_crs(3979)
        g = per[per["mergid"].isin(mergids)].geometry
        return LC.landcover_layer(year, list(g))

    return C.cached("landcover", year, params, [tif], build, geo=True)


def lcpng_cached(year: int, mergids) -> dict:
    """{str(mergid): {b, d}} paletted-PNG land-cover overlays, one per fire.

    A dict of base64 blobs, so it uses the same small JSON cache as `firesim_cached`
    rather than parquet.
    """
    mergids = sorted(int(m) for m in mergids)
    params = {"buf": C.LC_BUFFER_M, "min_m": C.LC_PNG_MIN_M, "cap": C.LC_PNG_CAP_PX,
              "split": C.LC_YEAR_SPLIT, "mergids": mergids}
    return _json_cached("lcpng", year, params, [LC._tif(year)],
                        lambda: _build_lcpng(year, mergids))


def _build_lcpng(year: int, mergids) -> dict:
    from . import lcraster
    per = io.read_final(year).to_crs(3979)
    sub = per[per["mergid"].isin(mergids)]
    return lcraster.fire_overlays(year, sub["mergid"].tolist(), list(sub.geometry))


def spread_steps_cached(year: int, mergids) -> dict:
    """{str(mergid): {t, av, s95}} 12-hourly spread-rate series - see pyrocb/spread.py."""
    from . import spread as SP
    mergids = sorted(int(m) for m in mergids)
    return _json_cached("spread", year, {"mergids": mergids, "v": 2},
                        _inputs(year), lambda: SP.fire_spread_steps(year, mergids))


def _json_cached(stage: str, year: int, params: dict, inputs, builder) -> dict:
    """Cache a plain dict (not a frame) as JSON, keyed like `config.cached`."""
    key = C.params_hash({"p": params, "inputs": C._stamp(inputs), "code": C.CODE_VERSION})
    C.CACHE_DIR.mkdir(exist_ok=True)
    fp = C.CACHE_DIR / f"{stage}_{year}_{key}.json"
    if fp.exists():
        try:
            return json.loads(fp.read_text())
        except Exception:
            pass
    for old in C.CACHE_DIR.glob(f"{stage}_{year}_*.json"):
        old.unlink()
    d = builder()
    fp.write_text(json.dumps(d))
    return d


def landcover_fractions_cached(year: int, active_gdf, events: pd.DataFrame) -> pd.DataFrame:
    """Land-cover top-3 composition of each event's `LC_ACTIVE_WINDOW` (6 h pre-onset)
    active-burn polygon, at the raster's native 30 m resolution.

    Blow-up events have an empty (6, 0) polygon; for those, fall back to the first
    `LC_BLOWUP_SPAN_H` hours of detected burning (`active.first_burn_polys`)."""
    pre, post = C.LC_ACTIVE_WINDOW
    sub = active_gdf[(active_gdf["pre_h"] == pre) & (active_gdf["post_h"] == post)]
    sub = sub[~sub.geometry.is_empty & sub.geometry.notna()][["pcb_id", "geometry"]]
    have = set(sub["pcb_id"])

    ev = events[(events["year"] == year) & events["mergid"].notna()]
    need = set(ev["pcb_id"]) - have
    if need:
        fb = A.first_burn_polys(year, events)
        if len(fb):
            fb = fb[fb["pcb_id"].isin(need)]
            sub = pd.concat([sub, fb[["pcb_id", "geometry"]]], ignore_index=True)
    sub = gpd.GeoDataFrame(sub, geometry="geometry", crs=C.METRIC_CRS)
    ids = sorted(sub["pcb_id"].tolist())
    params = {"win": [pre, post], "ids": ids, "split": C.LC_YEAR_SPLIT,
              "fallback": ["first", C.LC_BLOWUP_SPAN_H]}

    def build():
        if not len(sub):
            return pd.DataFrame(columns=["pcb_id", "lc_n_px"])
        g = sub.to_crs(3979)
        return LC.active_fractions(year, g["pcb_id"].tolist(), list(g.geometry))

    return C.cached("lcfrac", year, params, [LC._tif(year)], build)


def firesim_cached(year: int, events: pd.DataFrame, linked) -> dict:
    """{str(mergid): whole-fire evolution payload} for the year's pyroCb fires.

    Returns a plain dict (not a frame), so it uses a small JSON cache keyed on the
    firesim params, the pyroCb onsets and the NFP file stamps.
    """
    linked = {int(m) for m in linked}
    ev = events[(events["year"] == year) & events["mergid"].notna()].copy()
    ev = ev[ev["mergid"].astype("int64").isin(linked)]
    onsets = {}
    for m, g in ev.groupby("mergid"):
        onsets[int(m)] = [
            {"t": str(pd.Timestamp(r.datetime)), "lat": float(r.lat),
             "lon": float(r.lon), "id": str(r.pcb_id),
             "ix": int(r.pcb_index) if pd.notna(r.pcb_index) else None}
            for r in g.sort_values("datetime").itertuples()]
    return _json_cached("firesim", year, {"cell": C.BURN_CELL_M, "onsets": onsets},
                        _inputs(year), lambda: FS.fire_series(year, onsets))


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
    ("active_area_km2_pre12_post12", "km2", "ACTIVE BURN AREA: union of VIIRS pixel footprints from 12 h before to 12 h after onset, localised to the pyroCb's own fire cluster. The sampling geometry for land cover / ERA5", "NFP"),
    ("active_area_km2_pre{A}_post{B}", "km2", "Same for other windows: pre6_post6, pre24_post24, pre0_post24 (pure post, for blow-ups), pre24_post0, pre6_post0 (the 6 h-before-onset burn, sampled for land cover)", "NFP"),
    ("lc_top{1,2,3}_cls / _frac", "- / frac", "Top 3 Land Cover of Canada classes and their fraction of the pre6_post0 active-burn polygon, at native 30 m; lc_n_px = 30 m cells sampled", "Land Cover of Canada 2015/2020"),
    ("active_n_pixels_pre12_post12", "-", "Detections in the localised window", "NFP"),
    ("active_n_obs_pre12_post12", "-", "Distinct overpass times in the window; 1 means a single look", "NFP"),
    ("active_n_before_onset_pre12_post12", "-", "Of those detections, how many precede onset (0 => blow-up, see active_is_blowup)", "NFP"),
    ("active_n_after_onset_pre12_post12", "-", "Of those detections, how many follow onset", "NFP"),
    ("active_frp_p95_pre12_post12", "MW", "95th percentile FRP of the active-burn pixels", "NFP"),
    ("active_detach_km_pre12_post12", "km", "Distance from the pyroCb location to the kept fire cluster; 0 if it sits inside one", "NFP"),
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

    if active is not None and len(active):
        a = active.copy()
        # GPKG cannot store an empty geometry; a null geom marks "no detection in window".
        a.loc[a.geometry.is_empty, "geometry"] = None
        for c in a.columns:
            if pd.api.types.is_datetime64_any_dtype(a[c]):
                a[c] = a[c].dt.strftime("%Y-%m-%d %H:%M:%S")
        if a.crs is None:
            a = a.set_crs(C.METRIC_CRS)
        a.to_file(C.OUT_ACTIVE, driver="GPKG", layer="active_burn")
        if verbose:
            nn = int(a.geometry.notna().sum())
            nw = a.groupby(["pre_h", "post_h"]).ngroups
            print(f"  wrote {C.OUT_ACTIVE.name} ({len(a)} rows, {nn} with geometry, "
                  f"{nw} windows)", flush=True)
    return df
