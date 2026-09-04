"""Assemble the per-pyroCb-event attribute table.

One row per pyroCb event, every event kept - unmatched ones carry null fire columns
and `match_method == "unmatched"` rather than being dropped, because the unmatched set
is non-random (Cordillera terrain/cloud, unmapped steppe fuel) and silently discarding
it would bias the sample toward well-detected boreal forest fires.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import io
from . import timebins as T


# ------------------------------------------------------------------ pyroCb sequence
def pcb_sequence(ev: pd.DataFrame) -> pd.DataFrame:
    """Rank each event within its fire: does one pyroCb beget more?"""
    d = ev.sort_values(["mergid", "datetime"]).copy()
    g = d.groupby("mergid", dropna=False)
    d["pcb_index"] = g.cumcount() + 1
    d["n_pcb_total"] = g["pcb_id"].transform("size")
    d["n_pcb_before"] = d["pcb_index"] - 1
    d["n_pcb_after"] = d["n_pcb_total"] - d["pcb_index"]
    prev = g["datetime"].shift(1)
    nxt = g["datetime"].shift(-1)
    d["hours_since_prev_pcb"] = (d["datetime"] - prev).dt.total_seconds() / 3600
    d["hours_to_next_pcb"] = (nxt - d["datetime"]).dt.total_seconds() / 3600
    assert (d["n_pcb_before"] + d["n_pcb_after"] + 1 == d["n_pcb_total"]).all()
    return d.round({"hours_since_prev_pcb": 2, "hours_to_next_pcb": 2})


# ------------------------------------------------------------------ ignitions
def ignition_vars(year: int, ev: pd.DataFrame) -> pd.DataFrame:
    """Ignition counts before/after onset, by `mergid`.

    `mergid` groups the separately-ignited fires that later grew together (1,444 of
    12,151 rows in 2023; one 2023 mergid absorbs 115 ignitions), so n_ign_total is the
    size of the merge complex - and the variable that tells you whether the concentric
    band model applies to this fire at all.
    """
    # ign_t / ign_index now live on read_ignitions so the map and this table agree.
    ign = io.read_ignitions(year).assign(ign_halfday=lambda d: d["tst_ampm"])
    by = {k: v for k, v in ign.groupby("mergid")}

    rows = []
    for _, e in ev.iterrows():
        mid, t = e["mergid"], e["datetime"]
        if pd.isna(mid) or mid not in by:
            rows.append({})
            continue
        g = by[mid]
        before = int((g["ign_t"] < t).sum())
        # ignitions stamped in the same half-day as onset are genuinely unresolvable
        same = int((abs(g["ign_t"] - t) < pd.Timedelta(hours=7)).sum())
        span = (g["ign_t"].max() - g["ign_t"].min()).total_seconds() / 86400
        rows.append(dict(
            n_ign_total=len(g), n_ign_before=before, n_ign_after=len(g) - before,
            n_ign_same_halfday=same,
            first_ign_date=g["ign_t"].min(),
            days_first_ign_to_pcb=round((t - g["ign_t"].min()).total_seconds() / 86400, 3),
            ign_merge_span_days=round(span, 3),
        ))
    return pd.DataFrame(rows, index=ev.index)


# ------------------------------------------------------------------ NFP
def nfp_vars(year: int, ev: pd.DataFrame) -> pd.DataFrame:
    """Exact pixel/FRP splits around onset, from real `datetime_utc`.

    These are the cleanest variables in the product - unlike everything derived from
    the 12-hourly steps, they separate two pyroCb that fall in the same half-day.

    FRP caveat: not comparable across scan angles or satellites and season sums mix
    wildly uneven sample counts, so `frp_p95_prev24h` is the headline intensity
    variable and `n_obs_prev{W}h` must always be carried, so "zero because no
    overpass" stays distinguishable from "zero because no fire".
    """
    nfp = io.read_nfp_year(year)
    by = {k: v.sort_values("datetime_utc") for k, v in nfp.groupby("mergid")}

    rows = []
    for _, e in ev.iterrows():
        mid, t = e["mergid"], e["datetime"]
        if pd.isna(mid) or mid not in by:
            rows.append({})
            continue
        g = by[mid]
        tt, frp = g["datetime_utc"].values, g["frp"].values
        pre, post = tt < np.datetime64(t), tt >= np.datetime64(t)
        r = dict(
            n_pixels_total=len(g),
            n_pixels_before=int(pre.sum()), n_pixels_after=int(post.sum()),
            frp_sum_before=_f(frp[pre].sum()), frp_sum_after=_f(frp[post].sum()),
            frp_max_before=_f(frp[pre].max() if pre.any() else np.nan),
            frp_p95_before=_f(np.percentile(frp[pre], 95) if pre.any() else np.nan),
            hours_since_first_detection=round(
                (t - g["datetime_utc"].iloc[0]).total_seconds() / 3600, 2),
            hours_since_last_detection=round(
                (t - g.loc[pre, "datetime_utc"].max()).total_seconds() / 3600, 2)
            if pre.any() else np.nan,
        )
        for w in C.NFP_WINDOWS_H:
            lo = np.datetime64(t - pd.Timedelta(hours=w))
            hi = np.datetime64(t + pd.Timedelta(hours=w))
            a, b = (tt >= lo) & (tt < np.datetime64(t)), (tt >= np.datetime64(t)) & (tt < hi)
            r[f"n_pixels_prev{w}h"] = int(a.sum())
            r[f"n_pixels_next{w}h"] = int(b.sum())
            r[f"n_obs_prev{w}h"] = int(pd.unique(tt[a]).size)
            r[f"frp_p95_prev{w}h"] = _f(np.percentile(frp[a], 95) if a.any() else np.nan)
            r[f"frp_max_prev{w}h"] = _f(frp[a].max() if a.any() else np.nan)
            r[f"frp_p95_next{w}h"] = _f(np.percentile(frp[b], 95) if b.any() else np.nan)
        r["newpix_rate_prev24h"] = round(r["n_pixels_prev24h"] / 24.0, 3)
        r["qc_stale_obs"] = bool(r["hours_since_last_detection"] > C.STALE_OBS_H) \
            if r["hours_since_last_detection"] == r["hours_since_last_detection"] else True
        rows.append(r)
    return pd.DataFrame(rows, index=ev.index)


def _f(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 3)


# ------------------------------------------------------------------ fire attributes
def fire_vars(year: int, ev: pd.DataFrame) -> pd.DataFrame:
    fin = io.read_final(year).set_index("mergid")
    lon = fin.geometry.centroid.to_crs(C.WGS84).x
    keep = fin[["farea", "fperim", "duration", "n_pixels", "lcc_final"]].rename(
        columns={"farea": "fire_farea_km2", "fperim": "fire_fperim",
                 "duration": "fire_duration_days", "n_pixels": "fire_n_pixels",
                 "lcc_final": "fire_lcc_final"})
    keep = keep.assign(
        fire_tst=_halfday(fin, "tst", lon), fire_ted=_halfday(fin, "ted", lon),
        fire_centroid_lon=lon.round(4),
        fire_centroid_lat=fin.geometry.centroid.to_crs(C.WGS84).y.round(4))
    out = ev[["mergid"]].join(keep, on="mergid")
    out["days_to_fire_end"] = (out["fire_ted"] - ev["datetime"]).dt.total_seconds() / 86400
    return out.drop(columns="mergid").round({"days_to_fire_end": 3})


def _halfday(fin, pre, lon):
    d = pd.to_datetime(dict(year=fin[f"{pre}_year"], month=fin[f"{pre}_month"],
                            day=fin[f"{pre}_day"]), errors="coerce")
    h = fin[f"{pre}_ampm"].map(T.NOMINAL_LOCAL_H).fillna(7.0).values
    return d + pd.to_timedelta(h - T.norm_lon(lon.values) / 15.0, unit="h")


# ------------------------------------------------------------------ assemble
def build_events(year: int, matches: pd.DataFrame, band_ev: pd.DataFrame,
                 verbose: bool = True) -> pd.DataFrame:
    ev = pcb_sequence(matches)
    ev = ev.merge(band_ev, on=["pcb_id", "mergid"], how="left")
    ev = pd.concat([ev.reset_index(drop=True),
                    ignition_vars(year, ev).reset_index(drop=True),
                    nfp_vars(year, ev).reset_index(drop=True),
                    fire_vars(year, ev).reset_index(drop=True)], axis=1)
    ev = ev.sort_values("datetime").reset_index(drop=True)

    # The strongest end-to-end check on the mergid join chain: NFP pixels summed per
    # fire must equal the final perimeter's n_pixels. Asserted for years whose NFP ids
    # are sound; for a repaired year (2021 ships lon/lat coordinates and stale ids, so
    # pixels are re-attributed by geometry) exactness is not recoverable, so it becomes
    # a per-row QC flag instead of a hard failure.
    repaired = io.nfp_repaired(year)
    ev["qc_nfp_repaired"] = repaired
    ok = ev["n_pixels_total"].notna() & ev["fire_n_pixels"].notna()
    bad = int((ev.loc[ok, "n_pixels_total"] != ev.loc[ok, "fire_n_pixels"]).sum()) if ok.any() else 0
    ev["qc_nfp_pixel_exact"] = ev["n_pixels_total"] == ev["fire_n_pixels"]
    if not repaired:
        assert bad == 0, f"{year}: NFP pixel total != final n_pixels for {bad} events"
    if verbose:
        note = ("NFP ids REPAIRED by geometry - NFP-derived columns are approximate; "
                f"{len(ev) - bad}/{len(ev)} events still exact"
                if repaired else "NFP/final pixel identity verified")
        print(f"  {year}: {len(ev)} event rows x {ev.shape[1]} columns  ({note})", flush=True)
    return ev
