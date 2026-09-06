"""Active burn area around pyroCb onset.

The polygon of ground actually detected burning around the time of a pyroCb - the
sampling geometry for fractional land cover, ERA5 and every other covariate to come.

Built from NFP pixel footprints, not the published polygons, because:

* NFP carries a real `datetime_utc`, so the window resolves to the hour. Anything
  from the 12-hourly snapshots is stuck at half-day resolution and cannot separate
  two pyroCb on the same afternoon.
* it is genuine *detected* burning, not a hull fitted around it.

**Window** = (hours before onset, hours after onset), stepped over `ACTIVE_WINDOWS`.
Asymmetric on purpose: a blow-up fire (F4) is first seen by the overpass *after*
onset, so a pure-post window is the honest choice there.

**Localisation.** NFP `mergid` is the season-*final* id (F5), so a narrow time filter
on `by[mergid]` still pulls in pixels from every fire that eventually merged into this
one - which for a big fire can be several fronts 100+ km apart, only one of which made
the pyroCb. Measured on the matched set, half the in-window pixels were of that kind.
So after the time filter the detections are clustered (KD-tree connected components at
`ACTIVE_BRIDGE_M` + a pixel diagonal, which bridges swath gaps) and only the cluster
holding the pyroCb's own location is kept - or the nearest cluster within
`ACTIVE_ATTACH_CAP_KM` if the (0.1 deg-rounded) location lands in a gap.
`n_clusters` and `detach_km` record what happened, so the choice stays auditable.

Each VIIRS detection is its DS x DT footprint (along-scan x along-track, km, both
scan-angle dependent) as an axis-aligned rectangle in EPSG:3571. Orientation is not
in NFP so the rectangle is not rotated to the scan direction - a small approximation
at 375-780 m per side.
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import config as C
from . import io


def active_polygons(year: int, events: pd.DataFrame,
                    windows=None, verbose: bool = True) -> gpd.GeoDataFrame:
    """One polygon per (pcb_id, pre_h, post_h) for every matched event of `year`.

    Localisation (the expensive part) is done ONCE per event on the widest window;
    each narrower window is a time-subset of that localised pixel set.
    """
    windows = windows or C.ACTIVE_TABLE_WINDOWS
    rows = []
    for e, core in _cores(year, events):
        for pre, post in windows:
            rows.append(_window_row(e, core, pre, post))
    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=C.METRIC_CRS)
    if verbose:
        _report(year, out, windows)
    return out


def _grid_areas(core, cell_m=100.0) -> dict:
    """{"pre_post": area_km2} for every window in ACTIVE_WINDOWS, rasterised.

    A real `union_all` per window is ~70 ms, which over 24 windows x 400 events is half
    an hour; at 100 m against 375 m detections the rasterised area is under a percent
    off, and the authoritative polygons (ACTIVE_TABLE_WINDOWS) stay exact unions.
    """
    x, y, dx, dy, tr = core["x"], core["y"], core["dx"], core["dy"], core["t_rel_h"]
    x0, y0 = (x - dx).min(), (y - dy).min()
    nx = int((x + dx).max() - x0) // int(cell_m) + 2
    ny = int((y + dy).max() - y0) // int(cell_m) + 2
    i0 = ((x - dx - x0) / cell_m).astype(int); i1 = ((x + dx - x0) / cell_m).astype(int) + 1
    j0 = ((y - dy - y0) / cell_m).astype(int); j1 = ((y + dy - y0) / cell_m).astype(int) + 1
    g = {}
    for pre, post in C.ACTIVE_WINDOWS:
        m = np.flatnonzero((tr >= -pre) & (tr <= post))
        if not len(m):
            g[f"{pre}_{post}"] = 0.0
            continue
        grid = np.zeros((nx, ny), bool)
        for k in m:
            grid[i0[k]:i1[k], j0[k]:j1[k]] = True
        g[f"{pre}_{post}"] = round(float(grid.sum()) * (cell_m ** 2) / 1e6, 2)
    return g


def slider_data(year: int, events: pd.DataFrame):
    """One localisation pass -> (area_grid, cell_grids) for the map's window slider /
    animation. `area_grid` is {pcb_id: {"pre_post": km2}}; `cell_grids` is
    {pcb_id: burngrid payload}. Sharing `_cores` avoids re-grouping the NFP file.

    `collapse_first=False`: over a +/-24 h window a cell genuinely re-burns across
    successive overpasses, and the static window view has to be exact, so entries are
    unique (cell, observation time) rather than first-burn-only.
    """
    from . import burngrid
    grid, cells = {}, {}
    for e, core in _cores(year, events):
        if core is None or core.get("empty"):
            continue
        grid[e["pcb_id"]] = _grid_areas(core)
        enc = burngrid.encode(core["x"], core["y"], core["t_rel_h"] * 60.0,
                              collapse_first=False)
        if enc is not None:
            cells[e["pcb_id"]] = enc
    return grid, cells


def first_burn_polys(year: int, events: pd.DataFrame, span_h: float = None) -> gpd.GeoDataFrame:
    """Land-cover fallback geometry for blow-up fires (F4), which have no detected
    burning in the 6 h *before* onset. Union of the localised VIIRS pixel footprints
    in the first `span_h` hours *at or after* onset; if the fire is not even seen by
    then, the first `span_h` hours of whatever detections exist.

    Same construction as the (6, 0) window - actual detected-burning pixels only,
    localised to the pyroCb's own cluster - just a different time slice.
    """
    span_h = C.LC_BLOWUP_SPAN_H if span_h is None else span_h
    rows = []
    for e, core in _cores(year, events):
        if core is None or core.get("empty"):
            continue
        t = core["t_rel_h"]
        m = (t >= 0.0) & (t <= span_h)
        if not m.any():
            m = t <= float(t.min()) + span_h
        if not m.any():
            continue
        x, y, dx, dy = core["x"][m], core["y"][m], core["dx"][m], core["dy"][m]
        geom = io.ensure_valid([shapely.union_all(
            shapely.box(x - dx, y - dy, x + dx, y + dy), grid_size=C.PRECISION_M)])[0]
        rows.append({"pcb_id": e["pcb_id"], "geometry": geom})
    if not rows:
        return gpd.GeoDataFrame({"pcb_id": []},
                                geometry=gpd.GeoSeries([], crs=C.METRIC_CRS))
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=C.METRIC_CRS)


def _cores(year: int, events: pd.DataFrame):
    """Yield (event_row, core) for every matched event. `core` is a dict with the
    widest-window pixels already localised to the event's own fire cluster, or None
    if there were no detections at all in the widest window."""
    wide = max(C.ACTIVE_STEPS_H)
    ev = events[(events["year"] == year) & events["mergid"].notna()].copy()
    if not len(ev):
        return
    ev["mergid"] = ev["mergid"].astype("int64")
    pts = gpd.GeoSeries(gpd.points_from_xy(ev["lon"], ev["lat"]), crs=C.WGS84) \
        .to_crs(C.METRIC_CRS).values
    by = {k: v for k, v in io.read_nfp_year(year).groupby("mergid")}

    for (_, e), pt in zip(ev.iterrows(), pts):
        g = by.get(e["mergid"])
        t = pd.Timestamp(e["datetime"])
        if g is None or not len(g):
            yield e, None
            continue
        tt = g["datetime_utc"].values
        sel = (tt >= np.datetime64(t - pd.Timedelta(hours=wide))) & \
              (tt <= np.datetime64(t + pd.Timedelta(hours=wide)))
        if not sel.any():
            yield e, None
            continue
        s = g[sel]
        x = s["x"].to_numpy(float); y = s["y"].to_numpy(float)
        dx = s["DS"].to_numpy(float) * 500.0; dy = s["DT"].to_numpy(float) * 500.0
        mask, n_clusters, detach_km = _localise(x, y, dx, dy, pt)
        if not mask.any():
            yield e, {"empty": True, "n_clusters": n_clusters, "detach_km": detach_km}
            continue
        s = s[mask]
        yield e, {
            "onset": t,
            "t_rel_h": (pd.to_datetime(s["datetime_utc"]) - t).dt.total_seconds().to_numpy() / 3600.0,
            "dt_utc": s["datetime_utc"].to_numpy(),
            "x": x[mask], "y": y[mask], "dx": dx[mask], "dy": dy[mask],
            "frp": s["frp"].to_numpy(float),
            "n_clusters": n_clusters, "detach_km": detach_km,
        }


def _localise(x, y, dx, dy, pt):
    """Cluster detections and keep only the one at the pyroCb.

    Returns (mask, n_clusters, detach_km). Clustering is a KD-tree connected-component
    join at `ACTIVE_BRIDGE_M` plus a pixel diagonal - equivalent to buffering the boxes
    together, but O(n log n) instead of a shapely buffer/unbuffer over thousands of
    parts, which dominated the runtime.
    """
    n = len(x)
    if n == 0:
        return np.zeros(0, bool), 0, np.nan
    eps = C.ACTIVE_BRIDGE_M + 2.0 * max(float(np.median(dx)), float(np.median(dy)))
    pairs = cKDTree(np.column_stack([x, y])).query_pairs(eps, output_type="ndarray")
    if len(pairs):
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    else:
        g = coo_matrix((n, n))
    ncomp, lab = connected_components(g, directed=False)

    px, py = pt.x, pt.y
    inside = (np.abs(x - px) <= dx) & (np.abs(y - py) <= dy)
    if inside.any():
        pick, detach = lab[int(np.argmax(inside))], 0.0
    else:
        d = np.maximum(np.hypot(x - px, y - py) - np.hypot(dx, dy), 0.0)
        j = int(np.argmin(d))
        if d[j] > C.ACTIVE_ATTACH_CAP_KM * 1000.0:
            return np.zeros(n, bool), ncomp, round(float(d[j]) / 1000.0, 2)
        pick, detach = lab[j], round(float(d[j]) / 1000.0, 2)
    return lab == pick, ncomp, detach


def _window_row(e, core, pre, post):
    base = dict(pcb_id=e["pcb_id"], year=int(e["year"]), mergid=int(e["mergid"]),
                pre_h=int(pre), post_h=int(post), datetime=pd.Timestamp(e["datetime"]))
    if core is None:
        return {**base, **_null(), "geometry": shapely.Polygon()}
    if core.get("empty"):
        return {**base, **_null(), "n_clusters": core["n_clusters"],
                "detach_km": core["detach_km"], "geometry": shapely.Polygon()}

    m = (core["t_rel_h"] >= -pre) & (core["t_rel_h"] <= post)
    if not m.any():
        return {**base, **_null(), "n_clusters": core["n_clusters"],
                "detach_km": core["detach_km"], "geometry": shapely.Polygon()}

    x, y, dx, dy = core["x"][m], core["y"][m], core["dx"][m], core["dy"][m]
    geom = io.ensure_valid([shapely.union_all(
        shapely.box(x - dx, y - dy, x + dx, y + dy), grid_size=C.PRECISION_M)])[0]
    frp = core["frp"][m]
    dtu = core["dt_utc"][m]
    ts = pd.to_datetime(dtu)
    t = core["onset"]
    return {**base,
            "n_pixels": int(m.sum()),
            "n_obs": int(ts.nunique()),
            "area_km2": round(io.km2(geom), 4),
            "footprint_km2": round(float((dx * dy * 4).sum()) / 1e6, 4),
            "frp_sum": round(float(frp.sum()), 2),
            "frp_max": round(float(frp.max()), 2),
            "frp_p95": round(float(np.percentile(frp, 95)), 2),
            "t_first": ts.min(), "t_last": ts.max(),
            "n_before_onset": int((dtu < np.datetime64(t)).sum()),
            "n_after_onset": int((dtu >= np.datetime64(t)).sum()),
            "n_clusters": int(core["n_clusters"]), "detach_km": core["detach_km"],
            "geometry": geom}


def _null():
    return dict(n_pixels=0, n_obs=0, area_km2=0.0, footprint_km2=0.0,
                frp_sum=np.nan, frp_max=np.nan, frp_p95=np.nan,
                t_first=pd.NaT, t_last=pd.NaT, n_before_onset=0, n_after_onset=0,
                n_clusters=0, detach_km=np.nan)


def _empty():
    cols = (["pcb_id", "year", "mergid", "pre_h", "post_h", "datetime"] + list(_null()))
    return gpd.GeoDataFrame({c: [] for c in cols},
                            geometry=gpd.GeoSeries([], crs=C.METRIC_CRS))


def _tag(pre, post):
    return f"pre{int(pre)}_post{int(post)}"


def summarise(active: gpd.GeoDataFrame) -> pd.DataFrame:
    """Wide per-event columns to merge into the attribute table.

    The headline window (`ACTIVE_HEADLINE`) gets the full metric set; the other
    `ACTIVE_TABLE_WINDOWS` get area + pixel count only.
    """
    if not len(active):
        return pd.DataFrame(columns=["pcb_id"])
    out = active[["pcb_id"]].drop_duplicates().set_index("pcb_id")
    full = ["area_km2", "n_pixels", "n_obs", "frp_p95",
            "n_before_onset", "n_after_onset", "detach_km"]
    for pre, post in C.ACTIVE_TABLE_WINDOWS:
        sub = active[(active.pre_h == pre) & (active.post_h == post)].set_index("pcb_id")
        cols = full if (pre, post) == C.ACTIVE_HEADLINE else ["area_km2", "n_pixels"]
        for c in cols:
            out[f"active_{c}_{_tag(pre, post)}"] = sub[c]
    return out.reset_index()


def _report(year, out, windows):
    hp, ho = C.ACTIVE_HEADLINE
    h = out[(out.pre_h == hp) & (out.post_h == ho)]
    det = int((h.detach_km.fillna(0) > 0).sum())
    print(f"  {year}: active-burn polygons for {out.pcb_id.nunique()} events "
          f"x {len(windows)} windows "
          f"| headline (±{hp}h): median {h.area_km2.median():.1f} km2, "
          f"max {h.area_km2.max():.1f}, {det} attached to nearest cluster", flush=True)
