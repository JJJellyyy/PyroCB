"""Matched fire, v2: a nearby-detection centroid that ignores fire identity.

`matching.py` answers "which *fire* (mergid) produced this pyroCb?" and is careful
about causality: it only looks at fires as they stood at-or-before onset, and it
keeps the pyroCb attached to one tracked fire even across merges. That carefulness
is also its main assumption - if the atlas coordinate is a few km off, or two fires
are burning shoulder-to-shoulder, the "one correct mergid" framing can point the
matched location at the wrong lobe of ground.

This module is a second opinion that makes no such commitment. For each pyroCb event,
independently at every radius in `config.NEARBY_RADII_KM` and every strategy in
`config.NEARBY_METHODS`:

    1. Take the reported (lat, lon, datetime) at face value.
    2. Pick which overpass(es) to look at - see "Methods" below.
    3. Keep every NFP detection across those overpass(es) whose centre is within the
       radius of the reported point, regardless of mergid.
    4. The v2 position is the UNWEIGHTED mean of those detection centres, and the
       polygon behind it is the union of their VIIRS pixel footprints.

**Why NFP detections and not the Snapshot polygon.** The snapshot geometry is an
alpha-shape hull over *cumulative* burned area: it contains interior ground that is
not burning, and it counts everything burned so far rather than what was alight at
onset. Its area centroid is therefore the middle of a hull, not the middle of the
fire. NFP is incremental - consecutive steps share no coordinates - so one step's
detections are exactly the pixels newly seen burning at that overpass. That is the
thing worth taking a mean of.

Unweighted is deliberate: every detection counts once, so this is the centre of the
*detected burning area*, not of radiative power. (FRP is carried per event anyway,
via `active.py`.)

**Methods** (`config.NEARBY_METHODS`) - which overpass(es) feed step 3, in order of
preference; a method's later tiers are only tried if an earlier tier's widest radius
finds nothing at all (a regional swath gap, not a real absence - D-57):
    "closest"          The single nearest-in-time overpass, before OR after onset,
                       whichever is nearer. No causal claim, just "closest look".
    "before"           The nearest overpass AT OR BEFORE onset. Falls through to
                       AFTER-onset overpasses only once every before-onset tier is
                       exhausted - i.e. for a blow-up fire (F4: no detections exist
                       before the atlas onset time at all), not as a first resort.
                       This is structural, not a separate blow-up flag: "before" tiers
                       are tried first because they are listed first, "after" tiers
                       are what's left once they've all come up empty.
    "before_after"     The nearest before-onset AND nearest after-onset overpass,
                       POOLED into one search - captures the fire's state straddling
                       onset in a single answer, at roughly double the sample size of
                       a single-overpass method.
    "before2_after2"   The nearest TWO before-onset and nearest TWO after-onset
                       overpasses, pooled - the widest net, ~2 days of overpasses
                       either side at the normal 12 h cadence.
No method is "the" answer; they are intentionally different ways of drawing the same
line, shipped side by side so a reviewer can see which one performs best. See D-66.

Two entry points, mirroring `active.py`'s split:
    `nearby_polygons`  one row per (pcb_id, radius_km, method) with the pixel-
                       footprint geometry - what the map draws as "active fire used
                       for centering", and what a reviewer would want to inspect
                       directly.
    `summarise`        pivots that into wide per-event columns
                       (`fire2_r{R}_{method}_lat/lon/...`) for the attribute table.
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree

from . import config as C
from . import io
from . import timebins as T


def nearby_polygons(year: int, radii_km=None, methods=None, verbose: bool = True) -> gpd.GeoDataFrame:
    """One row per (pcb_id, radius_km, method) for every pyroCb event of `year`.

    Columns: pcb_id, radius_km, method, step, dt_h, det_dt_h, lat, lon, area_km2,
    n_pixels, n_fires, mergids, geometry (METRIC_CRS; empty Polygon where there was
    no v2 position).
    """
    radii = tuple(sorted(radii_km or C.NEARBY_RADII_KM))
    methods = tuple(methods or C.NEARBY_METHODS)
    pcb = io.read_pcb()
    p = pcb[pcb["year"] == year].copy()
    if not len(p):
        return _empty()

    # The overpasses that actually carry detections - not the snapshot file list,
    # since a step with a snapshot but no NFP rows is useless to this method.
    nfp = io.read_nfp_year(year)
    avail = frozenset(nfp["step"].unique())
    by_step = {k: v for k, v in nfp.groupby("step")}
    # One KD-tree per overpass, built on first use and shared across every event and
    # method that touches that step - a step's tree does not depend on who's asking.
    trees: dict[str, cKDTree] = {}

    rows = []
    for r in p.itertuples():
        anchor = T.step_containing(r.datetime, r.lon)
        before = _walk_available(anchor, avail, "back", C.NEARBY_MAX_STEPS)
        after = _walk_available(T.next_step(anchor), avail, "fwd", C.NEARBY_MAX_STEPS)
        px, py = _to_metric([r.lon], [r.lat])
        px, py = float(px[0]), float(py[0])

        for method in methods:
            tiers = _method_tiers(before, after, method, r.datetime, r.lon)
            rows += _match_event(r, tiers, by_step, trees, px, py, radii, method)

    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=C.METRIC_CRS)
    if verbose:
        _report(year, out, radii, methods)
    return out


def summarise(gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Wide per-event table: `fire2_r{R}_{method}_{lat,lon,step,dt_h,det_dt_h,
    area_km2,n_pixels,n_fires}` for every (radius, method) present in `gdf`."""
    if not len(gdf):
        return pd.DataFrame(columns=["pcb_id"])
    out = gdf[["pcb_id"]].drop_duplicates().set_index("pcb_id")
    for (radius, method), sub in gdf.groupby(["radius_km", "method"]):
        pre = f"fire2_{_key(radius, method)}"
        sub = sub.set_index("pcb_id")
        for c in ("lat", "lon", "step", "dt_h", "det_dt_h", "area_km2",
                  "n_pixels", "n_fires"):
            out[f"{pre}_{c}"] = sub[c]
    return out.reset_index()


# ---------------------------------------------------------------- method -> tiers
def _method_tiers(before: list[str], after: list[str], method: str,
                  t: pd.Timestamp, lon: float) -> list[list[str]]:
    """Ordered list of step-groups to try for `method`, most-preferred first. Each
    group's detections are pooled into ONE search; `_match_event` stops at the first
    group that finds anything within the widest radius."""
    if method == "closest":
        merged = sorted(before + after, key=lambda s: abs(_dt_hours(s, lon, t)))
        return [[s] for s in merged]
    if method == "before":
        return [[s] for s in before] + [[s] for s in after]
    if method == "before_after":
        n = max(len(before), len(after))
        groups = []
        for i in range(n):
            g = ([before[i]] if i < len(before) else []) + ([after[i]] if i < len(after) else [])
            if g:
                groups.append(g)
        return groups
    if method == "before2_after2":
        g = before[:2] + after[:2]
        return [g] if g else []
    raise ValueError(f"unknown NEARBY_METHODS entry: {method!r}")


def _walk_available(start: str, avail: frozenset, direction: str, n: int) -> list[str]:
    """Up to `n` available steps, nearest first, walking `direction` ('back' or
    'fwd') from `start` INCLUSIVE. For "before" (`start=anchor`, direction='back')
    this matches `timebins._walk_back`'s at-or-before contract; for "after"
    (`start=next_step(anchor)`, direction='fwd') it is strictly-after by construction
    of the caller's start point - the D-57 fix that made "after" mean something
    different from "before" whenever the anchor itself has data."""
    step_fn = T.prev_step if direction == "back" else T.next_step
    out, s, tries = [], start, 0
    max_tries = C.MAX_BACKWALK_STEPS * max(1, n)
    while len(out) < n and tries < max_tries:
        if s in avail:
            out.append(s)
        s = step_fn(s)
        tries += 1
    return out


# ---------------------------------------------------------------- matching
def _query_step(step: str, by_step: dict, trees: dict, px: float, py: float,
                radius_m: float) -> np.ndarray:
    d = by_step.get(step)
    if d is None or not len(d):
        return np.empty(0, dtype=int)
    if step not in trees:
        trees[step] = cKDTree(np.column_stack([d["x"].to_numpy(float), d["y"].to_numpy(float)]))
    return np.asarray(trees[step].query_ball_point([px, py], radius_m), dtype=int)


def _pool(group: list[str], by_step: dict, trees: dict, px: float, py: float,
         radius_m: float) -> pd.DataFrame | None:
    """Every detection from every step in `group` within `radius_m` of (px, py),
    concatenated with a `_step` tag recording which overpass each row came from."""
    parts = []
    for step in group:
        idx = _query_step(step, by_step, trees, px, py, radius_m)
        if len(idx):
            sub = by_step[step].iloc[idx].copy()
            sub["_step"] = step
            parts.append(sub)
    return pd.concat(parts, ignore_index=True) if parts else None


def _match_event(r, tiers: list[list[str]], by_step: dict, trees: dict,
                 px: float, py: float, radii: tuple[float, ...], method: str) -> list[dict]:
    widest = max(radii)
    for group in tiers:
        wide_pool = _pool(group, by_step, trees, px, py, widest * 1000.0)
        if wide_pool is None or not len(wide_pool):
            continue                               # this tier saw nothing - try the next
        rows = []
        for radius in radii:
            pool = wide_pool if radius == widest else \
                _pool(group, by_step, trees, px, py, radius * 1000.0)
            rows.append(_one_radius(r, group, pool, radius, method))
        return rows
    # No tier found anything at any radius.
    first_step = tiers[0][0] if tiers and tiers[0] else None
    return [_null_row(r, radius, first_step, method) for radius in radii]


def _one_radius(r, group: list[str], pool: pd.DataFrame | None,
                radius_km: float, method: str) -> dict:
    step_dt = min((_dt_hours(s, r.lon, r.datetime) for s in group), key=abs) if group else np.nan
    base = dict(pcb_id=r.pcb_id, radius_km=float(radius_km), method=method,
                step=",".join(group), dt_h=round(step_dt, 2) if step_dt == step_dt else np.nan)
    if pool is None or not len(pool):
        return {**base, **_null_fields(), "geometry": shapely.Polygon()}

    x = pool["x"].to_numpy(float)
    y = pool["y"].to_numpy(float)
    # THE v2 POSITION: unweighted mean of the (possibly multi-overpass) detection centres.
    clat, clon = _to_wgs84(float(x.mean()), float(y.mean()))

    # Each detection's real DS x DT footprint (km -> half-width in metres), unioned.
    hx = pool["DS"].to_numpy(float) * 500.0
    hy = pool["DT"].to_numpy(float) * 500.0
    geom = io.ensure_valid([shapely.union_all(
        shapely.box(x - hx, y - hy, x + hx, y + hy), grid_size=C.PRECISION_M)])[0]

    off_h = ((pd.to_datetime(pool["datetime_utc"]) - pd.Timestamp(r.datetime))
             .dt.total_seconds().to_numpy() / 3600.0)
    mergids = sorted({int(m) for m in pool["mergid"]})
    return {
        **base,
        "det_dt_h": round(float(np.median(off_h)), 2),
        "lat": round(clat, 5), "lon": round(clon, 5),
        "area_km2": round(geom.area / 1e6, 4),
        "n_pixels": int(len(pool)), "n_fires": len(mergids),
        "mergids": ",".join(str(m) for m in mergids),
        "geometry": geom,
    }


def _dt_hours(step: str, lon: float, t: pd.Timestamp) -> float:
    """Signed hours from onset to a step's nominal time: + after onset, - before."""
    return (T.nominal_utc(step, lon) - t).total_seconds() / 3600.0


# ---------------------------------------------------------------- geometry helpers
def _to_metric(lon, lat):
    g = gpd.GeoSeries(shapely.points(np.asarray(lon, float), np.asarray(lat, float)),
                      crs=C.WGS84).to_crs(C.METRIC_CRS)
    return g.x.values, g.y.values


def _to_wgs84(x: float, y: float) -> tuple[float, float]:
    """(lat, lon) for one METRIC_CRS point."""
    g = gpd.GeoSeries(shapely.points([x], [y]), crs=C.METRIC_CRS).to_crs(C.WGS84)
    return float(g.y.iloc[0]), float(g.x.iloc[0])


def _key(radius_km: float, method: str) -> str:
    """(10.0, "before_after") -> 'r10_before_after' (whole-km radii are the expected case)."""
    r = int(radius_km) if float(radius_km).is_integer() else radius_km
    return f"r{r}_{method}"


# ---------------------------------------------------------------- rows / reporting
def _null_fields() -> dict:
    return {"lat": np.nan, "lon": np.nan, "det_dt_h": np.nan, "area_km2": 0.0,
            "n_pixels": 0, "n_fires": 0, "mergids": ""}


def _null_row(r, radius_km: float, step: str | None, method: str) -> dict:
    dt = _dt_hours(step, r.lon, r.datetime) if step else np.nan
    return {"pcb_id": r.pcb_id, "radius_km": float(radius_km), "method": method, "step": step,
            "dt_h": round(dt, 2) if dt == dt else np.nan,
            **_null_fields(), "geometry": shapely.Polygon()}


def _empty() -> gpd.GeoDataFrame:
    cols = ["pcb_id", "radius_km", "method", "step", "dt_h", "det_dt_h", "lat", "lon",
            "area_km2", "n_pixels", "n_fires", "mergids"]
    return gpd.GeoDataFrame({c: [] for c in cols},
                            geometry=gpd.GeoSeries([], crs=C.METRIC_CRS))


def _report(year: int, out: gpd.GeoDataFrame, radii: tuple[float, ...],
           methods: tuple[str, ...]) -> None:
    n_events = out["pcb_id"].nunique()
    print(f"  {year}: nearby v2 ({n_events} events)", flush=True)
    for method in methods:
        parts = []
        for radius in radii:
            sub = out[(out["radius_km"] == radius) & (out["method"] == method)]
            parts.append(f"{radius:g}km {int(sub['lat'].notna().sum())}/{n_events}")
        print(f"    {method:<16} " + ", ".join(parts), flush=True)
