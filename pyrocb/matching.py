"""Associate each pyroCb event with the fire (`mergid`) that produced it.

Rule: containment first, then nearest ACTIVE fire within a distance cap.

Matching against the *snapshot* rather than the season-final perimeter matters. The
final footprint of a big fire can be 100 km across; a point falling inside it says
nothing about whether that fire had reached the spot on the day of the pyroCb, and a
fire that only ignited weeks later would still swallow the point. The snapshot at
onset is the fire as it actually was.

The floor on match quality is the atlas itself: lat/lon are stored to 0.1 deg, a
worst-case displacement of ~7.8 km, so `match_dist_km` above ~5 km is provisional.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import shapely

from . import config as C
from . import io
from . import timebins as T


def match_year(year: int, verbose: bool = True) -> pd.DataFrame:
    """One row per pyroCb event of `year`, with its matched fire and provenance."""
    pcb = io.read_pcb()
    p = pcb[pcb["year"] == year].copy()
    if not len(p):
        return _empty()

    avail = io.steps_set(year)
    p["match_step"] = [
        T._walk_back(T.step_containing(t, lon), avail)[0]
        for t, lon in zip(p["datetime"], p["lon"])
    ]

    rows = []
    # Group by step so each snapshot is read once (279 reads across all years, not 484).
    for step, grp in p.groupby("match_step", dropna=False):
        if step is None:
            rows += [_unmatched(r) for _, r in grp.iterrows()]
            continue
        snap = io.read_snapshot(year, step)
        pts = shapely.points(
            *np.array([[x, y] for x, y in zip(*_to_metric(grp["lon"], grp["lat"]))]).T
        )
        tree = shapely.STRtree(snap.geometry.values)
        areas = snap.geometry.area.values

        for (_, r), pt in zip(grp.iterrows(), np.atleast_1d(pts)):
            inside = tree.query(pt, predicate="within")
            if len(inside):
                # Overlapping candidates are rare; take the largest fire.
                k = inside[int(np.argmax(areas[inside]))]
                rows.append(_row(r, step, snap, k, 0.0, "inside"))
                continue
            k = tree.nearest(pt)
            d_km = float(shapely.distance(pt, snap.geometry.values[k])) / 1000.0
            if d_km <= C.MATCH_CAP_KM:
                rows.append(_row(r, step, snap, k, d_km, "nearest"))
            else:
                rows.append(_unmatched(r, step, d_km, int(snap["mergid"].values[k])))

    out = pd.DataFrame(rows).sort_values("datetime").reset_index(drop=True)
    out = _rescue_new_fire(year, out)
    out = _rescue_dormant(year, out)
    out = _apply_pixel_matching(out)
    if verbose:
        _report(year, out)
    return out


def _rescue_new_fire(year: int, out: pd.DataFrame) -> pd.DataFrame:
    """Look FORWARD for fires that blew up before they were first observed.

    There is a ~10 h nightly gap between overpass groups (no observation between local
    14.4 h and 00.5 h) and pyroCb onset peaks at local ~15.0 h - right inside it. A fire
    that ignites and blows up in that window does not exist in any snapshot at or before
    onset; its first look is the NEXT overpass.

    2023-R90877 is the type case: pyroCb at 2023-06-25 00:00 UTC, fire 11918 first
    detected 10.3 h later already carrying 117 pixels, so it was plainly burning at
    onset. Searching forward recovers it; searching only backwards cannot.

    This affects matching only. The pre-band still uses the last observation strictly
    BEFORE onset, so such an event correctly gets band0_undefined / area_pre = 0.
    """
    avail = io.steps_set(year)
    for k in range(C.MATCH_FORWARD_STEPS):
        miss = out["match_method"] == "unmatched"
        if not miss.any():
            break
        by_step = {}
        for i, r in out.loc[miss].iterrows():
            s = r["match_step"]
            if not isinstance(s, str):
                continue
            for _ in range(k + 1):
                s = T.next_step(s)
            if s in avail:
                by_step.setdefault(s, []).append(i)
        for step, idx in by_step.items():
            snap = io.read_snapshot(year, step)
            tree = shapely.STRtree(snap.geometry.values)
            areas = snap.geometry.area.values
            xs, ys = _to_metric(out.loc[idx, "lon"], out.loc[idx, "lat"])
            for i, x, y in zip(idx, xs, ys):
                pt = shapely.points(x, y)
                inside = tree.query(pt, predicate="within")
                if len(inside):
                    j = inside[int(np.argmax(areas[inside]))]
                    d_km, method = 0.0, "inside_newfire"
                else:
                    j = tree.nearest(pt)
                    d_km = float(shapely.distance(pt, snap.geometry.values[j])) / 1000.0
                    if d_km > C.MATCH_CAP_KM:
                        continue
                    method = "nearest_newfire"
                out.loc[i, ["mergid", "match_dist_km", "match_method",
                            "nearest_dist_km", "match_fwd_steps"]] = [
                    int(snap["mergid"].values[j]), round(d_km, 3), method, round(d_km, 3), k + 1]
    return out


def _apply_pixel_matching(out: pd.DataFrame) -> pd.DataFrame:
    """Overlay the dataset's pixel-based matching (pcbcube.matching, with the mass-based credit fix).

    Where it credits a fire, that fire's `mergid` wins - it was chosen from the detections actually
    seen within 12 h of onset, whereas the snapshot rules above can pick a neighbour (4 of 396 events
    differ) or lose the event (1 rescued). The snapshot rules stay the fallback for events the pixel
    matcher leaves unmatched. `match_method` becomes "pixels" for overridden rows.
    """
    fp = C.DATA.parent / "OUTPUT_DATA" / "datacube" / "pcb_events.parquet"
    if not fp.exists() or out.empty:
        return out
    ev = pd.read_parquet(fp, columns=["pcb_id", "matched", "mergid"])
    ev = ev[ev["matched"].astype(bool)].set_index("pcb_id")["mergid"]
    out = out.copy()
    out["mergid"] = pd.array(out["mergid"], dtype="Int64")
    for i, r in out.iterrows():
        new = ev.get(r["pcb_id"])
        if new is None or pd.isna(new):
            continue
        if pd.isna(r["mergid"]) or int(r["mergid"]) != int(new):
            out.loc[i, ["mergid", "match_method"]] = [int(new), "pixels"]
    return out


def _rescue_dormant(year: int, out: pd.DataFrame) -> pd.DataFrame:
    """Second tier for events unmatched against the snapshot.

    A snapshot lists only fires that were *actively detected* that overpass, so a fire
    obscured by cloud or missed by the swath vanishes from the file even though it is
    burning. Retry those events against the season-final perimeters, restricted to
    fires whose [start, end] interval brackets onset - which keeps the timing
    constraint that makes snapshot matching worth doing in the first place.
    """
    miss = out["match_method"] == "unmatched"
    if not miss.any():
        return out

    fin = io.read_final(year).copy()
    fin["t0"] = _halfday_utc(fin, "tst")
    fin["t1"] = _halfday_utc(fin, "ted")
    tree = shapely.STRtree(fin.geometry.values)
    areas = fin.geometry.area.values

    xs, ys = _to_metric(out.loc[miss, "lon"], out.loc[miss, "lat"])
    for (i, r), x, y in zip(out.loc[miss].iterrows(), xs, ys):
        pt = shapely.points(x, y)
        t = r["datetime"]
        cand = tree.query(pt, predicate="dwithin", distance=C.MATCH_CAP_KM * 1000.0)
        cand = [k for k in cand if fin["t0"].values[k] <= t <= fin["t1"].values[k]]
        if not cand:
            continue
        inside = [k for k in cand if shapely.within(pt, fin.geometry.values[k])]
        if inside:
            k = inside[int(np.argmax(areas[inside]))]
            d_km, method = 0.0, "inside_dormant"
        else:
            d = [float(shapely.distance(pt, fin.geometry.values[k])) / 1000.0 for k in cand]
            j = int(np.argmin(d))
            k, d_km, method = cand[j], d[j], "nearest_dormant"
        out.loc[i, ["mergid", "match_dist_km", "match_method", "nearest_dist_km"]] = [
            int(fin["mergid"].values[k]), round(d_km, 3), method, round(d_km, 3)]
    return out


def _halfday_utc(df, pre):
    """Nominal UTC instant of a tst_*/ted_* half-day stamp, via the fire's own centroid."""
    lon = df.geometry.centroid.to_crs(C.WGS84).x.values
    d = pd.to_datetime(dict(year=df[f"{pre}_year"], month=df[f"{pre}_month"],
                            day=df[f"{pre}_day"]), errors="coerce")
    h = df[f"{pre}_ampm"].map(T.NOMINAL_LOCAL_H).fillna(7.0).values
    return d + pd.to_timedelta(h - T.norm_lon(lon) / 15.0, unit="h")


def _to_metric(lon, lat):
    g = (__import__("geopandas")
         .GeoSeries(shapely.points(np.asarray(lon, float), np.asarray(lat, float)),
                    crs=C.WGS84).to_crs(C.METRIC_CRS))
    return g.x.values, g.y.values


_SNAP_KEEP = ["farea", "n_pixels", "n_newpixels", "flinelen", "pixden", "meanFRP",
              "FRP95", "t_inactive", "duration", "lake_area", "avgspread", "spread95",
              "n_components"]


def _row(r, step, snap, k, d_km, method):
    d = _base(r)
    d.update(mergid=int(snap["mergid"].values[k]), match_step=step,
             match_dist_km=round(d_km, 3), match_method=method,
             nearest_dist_km=round(d_km, 3))
    for c in _SNAP_KEEP:
        if c in snap.columns:
            v = snap[c].values[k]
            d["snap_" + c] = None if pd.isna(v) else float(v)
    return d


def _unmatched(r, step=None, d_km=np.nan, nearest=None):
    d = _base(r)
    d.update(mergid=pd.NA, match_step=step, match_dist_km=np.nan,
             match_method="unmatched", nearest_dist_km=round(d_km, 3) if d_km == d_km else np.nan,
             nearest_mergid=nearest, match_fwd_steps=0)
    return d


_PASS = ["pcb_id", "year", "datetime", "doy", "fire_name", "region", "status",
         "lat", "lon", "pcb_min_bt", "trop_alt", "max_inject_alt",
         "trop_relative_alt", "min_smoke_mass", "max_smoke_mass", "plume_code",
         "pulse_code", "notes"]


def _base(r):
    return {c: (r[c] if c in r.index else None) for c in _PASS}


def _empty():
    return pd.DataFrame(columns=_PASS + ["mergid", "match_step", "match_dist_km",
                                         "match_method", "nearest_dist_km"])


def _report(year, out):
    n = len(out)
    vc = out["match_method"].value_counts()
    print(f"  {year}: {n} events -> " +
          ", ".join(f"{k} {v}" for k, v in vc.items()) +
          f" | {out['mergid'].dropna().nunique()} distinct fires", flush=True)
    d = out.loc[out["match_method"] == "nearest", "match_dist_km"]
    if len(d):
        print(f"        nearest-match dist km: median {d.median():.2f} p90 {d.quantile(.9):.2f} max {d.max():.2f}",
              flush=True)
    # Cap sensitivity: the nearest fire never changes, only whether we accept it.
    nd = out["nearest_dist_km"]
    caps = ", ".join(f"{c:g}km={int((nd <= c).sum())}" for c in C.MATCH_CAP_SENSITIVITY)
    print(f"        matched at cap: {caps}  (of {n})", flush=True)
    bad = out.loc[out["match_method"] == "unmatched", "region"].value_counts()
    if len(bad):
        print("        unmatched by region: " + ", ".join(f"{k} {v}" for k, v in bad.items()),
              flush=True)
