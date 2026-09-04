"""Concentric growth bands: a fire's footprint split at every pyroCb it produced.

For a fire with N distinct pyroCb half-day steps:

    S_0     = V(snap(t_0))                  # burned before the first pyroCb
    S_i     = S_{i-1} u V(snap(t_i))        # cumulative, monotone BY CONSTRUCTION
    S_final = S_{N-1} u V(final_perim)
    B_i     = S_i \\ S_{i-1}

Building a running union rather than differencing raw snapshots is what makes the
bands disjoint, non-negative and exactly tiling S_final. Consecutive snapshots are not
perfectly nested (re-polygonisation loses ~0.17% of the earlier polygon on average),
so plain differencing would produce negative-area bands and rings that punch holes in
each other on the map. The union costs nothing extra and removes the whole class of
problem - but we record what it papered over in the qc_* columns.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

from . import config as C
from . import io
from . import timebins as T


def build_bands(year: int, matches: pd.DataFrame, verbose: bool = True):
    """Return (bands GeoDataFrame in METRIC_CRS, per-event step/QC DataFrame)."""
    m = matches[matches["mergid"].notna()].copy()
    if not len(m):
        return _empty_bands(), pd.DataFrame()
    m["mergid"] = m["mergid"].astype("int64")

    obs = io.nfp_obs_index(year)
    obs_by_fire = {k: v.sort_values("t_med") for k, v in obs.groupby("mergid")}
    avail = io.steps_set(year)

    # ---- 1. per event: which step precedes onset -------------------------
    recs = []
    for _, r in m.iterrows():
        o = obs_by_fire.get(r["mergid"])
        step, tier, lag = T.step_before(r["datetime"], o, r["lon"], avail)
        recs.append({"pcb_id": r["pcb_id"], "mergid": r["mergid"],
                     "datetime": r["datetime"], "step_before": step,
                     "step_tier": tier, "obs_lag_h": round(lag, 2) if lag == lag else np.nan})
    ev = pd.DataFrame(recs).sort_values(["mergid", "datetime"]).reset_index(drop=True)

    # ---- 2. fetch each fire's geometry at each needed step ---------------
    need = ev.dropna(subset=["step_before"])[["mergid", "step_before"]].drop_duplicates()
    geom = _fetch(year, need, avail, verbose)

    # ---- 3. cumulative union chain + differences -------------------------
    fin = io.read_final(year).set_index("mergid")
    band_rows, ev_rows = [], []

    for mergid, g in ev.groupby("mergid"):
        g = g.sort_values("datetime")
        steps = [s for s in dict.fromkeys(g["step_before"].dropna())]   # ordered, deduped
        steps.sort(key=T.step_sort_key)
        band_of_step = {s: i for i, s in enumerate(steps)}

        cum, cum_areas, dropped = None, [], 0.0
        rows = []
        for i, s in enumerate(steps):
            piece = geom.get((mergid, s))
            new = piece if cum is None else (shapely.union_all([cum, piece]) if piece is not None else cum)
            regress = 0.0 if cum is None or piece is None else io.km2(shapely.difference(cum, piece))
            band = new if cum is None else shapely.difference(new, cum)
            band, drp = io.drop_slivers(band)
            dropped += drp
            rows.append(dict(band_i=i, geometry=band, area_km2=io.km2(band),
                             step_lo="" if i == 0 else steps[i - 1], step_hi=s,
                             qc_regress_km2=round(regress, 4)))
            cum = new
            cum_areas.append(io.km2(cum))

        # final band: everything after the last pyroCb
        f = fin.geometry.get(mergid)
        if f is not None:
            f = io.snap([f])[0]
            s_final = f if cum is None else shapely.union_all([cum, f])
            shrink = 0.0 if cum is None else io.km2(shapely.difference(cum, f))
            band = s_final if cum is None else shapely.difference(s_final, cum)
            band, drp = io.drop_slivers(band)
            dropped += drp
            rows.append(dict(band_i=len(steps), geometry=band, area_km2=io.km2(band),
                             step_lo=steps[-1] if steps else "", step_hi="final",
                             qc_regress_km2=round(shrink, 4)))
            cum = s_final
        area_final = io.km2(cum)

        n_bands = len(rows)
        for r in rows:
            r.update(mergid=mergid, year=year, n_bands=n_bands,
                     is_post_final_band=(r["band_i"] == n_bands - 1),
                     cum_area_km2=round(area_final, 4),
                     qc_sliver_dropped_km2=round(dropped, 4))
            band_rows.append(r)

        # ---- 4. per-event area / band bookkeeping ------------------------
        for _, e in g.iterrows():
            s = e["step_before"]
            bi = band_of_step.get(s)
            pre = cum_areas[bi] if bi is not None and bi < len(cum_areas) else 0.0
            ev_rows.append(dict(
                pcb_id=e["pcb_id"], mergid=mergid, step_before=s, step_tier=e["step_tier"],
                obs_lag_h=e["obs_lag_h"], band_id=bi,
                band0_undefined=bool(s is None),
                area_pre_km2=round(pre, 4),
                area_final_km2=round(area_final, 4),
                area_post_km2=round(area_final - pre, 4),
                frac_burned_before=round(pre / area_final, 5) if area_final > 0 else np.nan,
                qc_final_smaller=bool(rows and rows[-1]["qc_regress_km2"] > 0.10 * area_final),
            ))

    bands = gpd.GeoDataFrame(band_rows, geometry="geometry", crs=C.METRIC_CRS)
    evdf = pd.DataFrame(ev_rows)
    evdf = _add_dedup(evdf)
    if verbose:
        _report(year, bands, evdf)
    return bands, evdf


def _fetch(year, need, avail, verbose):
    """Geometry of each (mergid, step), reading every snapshot at most once.

    A fire absent from its snapshot is walked backwards: snapshots list only actively
    detected fires, so a quiet fire's cumulative footprint is unreported, not shrunk -
    using its last observed footprint is correct, not an approximation.
    """
    out, back = {}, {}
    # (mergid, originally requested step, step currently being tried)
    pending = [(int(a), b, b) for a, b in need.itertuples(index=False)]
    for depth in range(C.MAX_BACKWALK_STEPS + 1):
        if not pending:
            break
        by_step, nxt = {}, []
        for mid, orig, cur in pending:
            by_step.setdefault(cur, []).append((mid, orig))
        for s, items in sorted(by_step.items(), key=lambda kv: T.step_sort_key(kv[0])):
            snap = io.read_snapshot(year, s).set_index("mergid")
            for mid, orig in items:
                if mid in snap.index:
                    out[(mid, orig)] = io.snap([snap.geometry.loc[mid]])[0]
                    back[(mid, orig)] = depth
                else:
                    p = T.prev_step(s)
                    if p in avail:
                        nxt.append((mid, orig, p))
        pending = nxt
    if verbose and back:
        w = sum(1 for v in back.values() if v)
        print(f"        geometry: {len(out)}/{len(need)} (mergid,step) resolved, "
              f"{w} needed a back-walk", flush=True)
    return out


def _add_dedup(ev: pd.DataFrame) -> pd.DataFrame:
    """Events sharing a fire and a step_before share a band - and share every
    step-derived number. Flag it so identical values do not look like a bug."""
    if not len(ev):
        return ev
    n = ev.groupby(["mergid", "step_before"], dropna=False)["pcb_id"].transform("size")
    ev["n_events_in_band"] = n.values
    ev["dedup_collapsed"] = ev["n_events_in_band"] > 1
    return ev


def _empty_bands():
    return gpd.GeoDataFrame(
        {c: [] for c in ["mergid", "year", "band_i", "n_bands", "area_km2", "step_lo",
                         "step_hi", "cum_area_km2", "is_post_final_band",
                         "qc_regress_km2", "qc_sliver_dropped_km2"]},
        geometry=gpd.GeoSeries([], crs=C.METRIC_CRS))


def _report(year, bands, ev):
    print(f"  {year}: {len(bands)} bands over {bands['mergid'].nunique()} fires", flush=True)
    nb = bands.groupby("mergid")["n_bands"].first().value_counts().sort_index()
    print("        bands per fire: " + ", ".join(f"{k}:{v}" for k, v in nb.items()), flush=True)
    print(f"        sliver area dropped {bands.groupby('mergid')['qc_sliver_dropped_km2'].first().sum():.3f} km2"
          f" | qc_final_smaller fires: {int(ev['qc_final_smaller'].sum())}", flush=True)
    if "obs_lag_h" in ev:
        l = ev["obs_lag_h"].dropna()
        if len(l):
            print(f"        obs_lag_h: median {l.median():.1f} p90 {l.quantile(.9):.1f} max {l.max():.1f}"
                  f" | events with lag<6h: {int((l < 6).sum())}/{len(l)}", flush=True)
    print(f"        tier: " + ", ".join(f"{k} {v}" for k, v in ev["step_tier"].value_counts().items()),
          flush=True)
