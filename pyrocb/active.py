"""Active burn area at pyroCb onset.

The polygon of ground actually detected burning around the time of a pyroCb - the
sampling geometry for fractional land cover, ERA5 and every other covariate to come.

Built from NFP pixel footprints rather than the published polygons, for two reasons:

* NFP carries a real `datetime_utc`, so the window resolves to the hour. Everything
  derived from the 12-hourly snapshots is stuck at half-day resolution and cannot
  separate two pyroCb on the same afternoon.
* The area is genuine *detected* burning, not a hull fitted around it.

The window is **symmetric** around onset. A purely-backward window would usually be
empty or badly stale: there is no observation between local 14.4 h and 00.5 h and
pyroCb onset peaks at local ~15.0 h (FINDINGS F3), so the run that fed the plume is
typically first seen by the overpass *after* onset.

Each VIIRS detection is represented by its own DS x DT footprint (along-scan x
along-track, km, both scan-angle dependent) as an axis-aligned rectangle in EPSG:3571.
Orientation is not recorded in NFP, so the rectangle is not rotated to the scan
direction - an approximation worth stating, though at 375-780 m per side it is small
against a fire front.
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from . import config as C
from . import io


def active_polygons(year: int, events: pd.DataFrame,
                    windows=None, verbose: bool = True) -> gpd.GeoDataFrame:
    """One polygon per (pcb_id, window_h) for every matched event of `year`."""
    windows = windows or C.ACTIVE_WINDOWS_H
    ev = events[(events["year"] == year) & events["mergid"].notna()].copy()
    if not len(ev):
        return _empty()
    ev["mergid"] = ev["mergid"].astype("int64")

    nfp = io.read_nfp_year(year)
    by = {k: v for k, v in nfp.groupby("mergid")}

    rows = []
    for _, e in ev.iterrows():
        g = by.get(e["mergid"])
        t = pd.Timestamp(e["datetime"])
        for w in windows:
            rows.append(_one(e, g, t, w))
    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=C.METRIC_CRS)
    if verbose:
        _report(year, out, windows)
    return out


def _one(e, g, t, w):
    base = dict(pcb_id=e["pcb_id"], year=int(e["year"]), mergid=int(e["mergid"]),
                window_h=w, datetime=t)
    if g is None or not len(g):
        return {**base, **_null(), "geometry": shapely.Polygon()}

    lo, hi = np.datetime64(t - pd.Timedelta(hours=w)), np.datetime64(t + pd.Timedelta(hours=w))
    tt = g["datetime_utc"].values
    sel = (tt >= lo) & (tt <= hi)
    if not sel.any():
        return {**base, **_null(), "geometry": shapely.Polygon()}

    s = g[sel]
    x, y = s["x"].to_numpy(float), s["y"].to_numpy(float)
    # DS/DT are the along-scan / along-track pixel dimensions in km.
    dx, dy = s["DS"].to_numpy(float) * 500.0, s["DT"].to_numpy(float) * 500.0
    boxes = shapely.box(x - dx, y - dy, x + dx, y + dy)
    geom = shapely.union_all(boxes, grid_size=C.PRECISION_M)
    geom = io.ensure_valid([geom])[0]

    frp = s["frp"].to_numpy(float)
    ts = pd.to_datetime(s["datetime_utc"])
    return {**base,
            "n_pixels": int(sel.sum()),
            "n_obs": int(ts.nunique()),
            "area_km2": round(io.km2(geom), 4),
            "footprint_km2": round(float((dx * dy * 4).sum()) / 1e6, 4),
            "frp_sum": round(float(frp.sum()), 2),
            "frp_max": round(float(frp.max()), 2),
            "frp_p95": round(float(np.percentile(frp, 95)), 2),
            "t_first": ts.min(), "t_last": ts.max(),
            "n_before_onset": int((s["datetime_utc"] < t).sum()),
            "n_after_onset": int((s["datetime_utc"] >= t).sum()),
            "geometry": geom}


def _null():
    return dict(n_pixels=0, n_obs=0, area_km2=0.0, footprint_km2=0.0,
                frp_sum=np.nan, frp_max=np.nan, frp_p95=np.nan,
                t_first=pd.NaT, t_last=pd.NaT, n_before_onset=0, n_after_onset=0)


def _empty():
    cols = ["pcb_id", "year", "mergid", "window_h", "datetime"] + list(_null())
    return gpd.GeoDataFrame({c: [] for c in cols},
                            geometry=gpd.GeoSeries([], crs=C.METRIC_CRS))


def summarise(active: gpd.GeoDataFrame) -> pd.DataFrame:
    """Wide per-event columns to merge into the attribute table."""
    if not len(active):
        return pd.DataFrame(columns=["pcb_id"])
    keep = ["area_km2", "n_pixels", "n_obs", "frp_p95", "n_before_onset", "n_after_onset"]
    w = active.pivot_table(index="pcb_id", columns="window_h", values=keep, aggfunc="first")
    w.columns = [f"active_{a}_{int(b)}h" for a, b in w.columns]
    return w.reset_index()


def _report(year, out, windows):
    print(f"  {year}: active-burn polygons for {out.pcb_id.nunique()} events "
          f"x {len(windows)} windows", flush=True)
    for w in windows:
        s = out[out.window_h == w]
        empty = int((s.area_km2 == 0).sum())
        print(f"        +/-{w:>2}h: median {s.area_km2.median():8.2f} km2  "
              f"max {s.area_km2.max():9.2f}  empty {empty}/{len(s)}", flush=True)
