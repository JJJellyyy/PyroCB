"""Whole-fire evolution series for the map's "animate the entire event" mode.

Every NFP detection of a pyroCb fire (the full merge complex, *not* localised - the
point is to watch the fronts grow and merge), from first detection to last, decimated
adaptively by how close in time it is to one of the fire's pyroCb onsets:

    | proximity to a pyroCb onset | tier | spatial dedup | time bin |
    |-----------------------------|------|---------------|----------|
    | <= FIRESIM_PCB_NEAR_H[0] h  |  A   | native        | native   |
    | <= FIRESIM_PCB_NEAR_H[1] h  |  B   | 1 km grid     | 6 h      |
    | otherwise                   |  C   | 3 km grid     | 24 h     |

So the ±24 h around every pyroCb stays near full resolution and the long quiet spans
between events are coarse. If a fire is still over FIRESIM_MAX_PTS after that, the B/C
grids are widened until it fits.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd

from . import config as C
from . import io


def _decimate(sec, x, y, tier, grid_m, bin_h):
    """Boolean keep-mask: tier A untouched; tiers B/C keep the first detection in each
    (tier, spatial cell, time bin)."""
    keep = np.ones(len(sec), bool)
    coarse = np.flatnonzero(tier > 0)
    if not len(coarse):
        return keep
    gm = np.where(tier[coarse] == 1, grid_m[1], grid_m[2])
    bh = np.where(tier[coarse] == 1, bin_h[1], bin_h[2])
    key = np.stack([
        tier[coarse].astype(np.int64),
        np.floor(x[coarse] / gm).astype(np.int64),
        np.floor(y[coarse] / gm).astype(np.int64),
        np.floor(sec[coarse] / (bh * 3600.0)).astype(np.int64),
    ], axis=1)
    _, first = np.unique(key, axis=0, return_index=True)
    drop = np.ones(len(coarse), bool)
    drop[first] = False
    keep[coarse[drop]] = False
    return keep


def fire_series(year: int, onsets: dict) -> dict:
    """{str(mergid): payload} for every mergid in `onsets`.

    `onsets[mergid]` is a list of {"t": UTC, "lat": ..., "lon": ..., "id": pcb_id}.
    """
    nfp = io.read_nfp_year(year)
    by = {k: v for k, v in nfp.groupby("mergid")}
    grid = np.array(C.FIRESIM_GRID_M, float)
    binh = np.array(C.FIRESIM_BIN_H, float)
    nearA, nearB = C.FIRESIM_PCB_NEAR_H
    out = {}

    for mergid, ons in onsets.items():
        g = by.get(int(mergid))
        if g is None or not len(g):
            continue
        g = g.sort_values("datetime_utc")
        t = pd.to_datetime(g["datetime_utc"]).to_numpy()
        t0 = t.min()
        sec = (t - t0) / np.timedelta64(1, "s")
        x = g["x"].to_numpy(float); y = g["y"].to_numpy(float)
        dx = g["DS"].to_numpy(float) * 500.0
        dy = g["DT"].to_numpy(float) * 500.0

        ons = sorted(ons, key=lambda o: pd.Timestamp(o["t"]))
        ons_s = np.array([(np.datetime64(pd.Timestamp(o["t"])) - t0) / np.timedelta64(1, "s")
                          for o in ons], float)
        d_h = np.abs(sec[:, None] - ons_s[None, :]).min(axis=1) / 3600.0
        tier = np.where(d_h <= nearA, 0, np.where(d_h <= nearB, 1, 2)).astype(np.int8)

        gm, bh = grid.copy(), binh.copy()
        for _ in range(6):
            keep = _decimate(sec, x, y, tier, gm, bh)
            if keep.sum() <= C.FIRESIM_MAX_PTS:
                break
            gm[1:] *= 1.6; bh[1:] *= 1.5
        k = np.flatnonzero(keep)

        # coarse tiers get a footprint the size of their dedup cell so the blocks tile
        hw = np.where(tier[k] == 0, dx[k], np.where(tier[k] == 1, gm[1] / 2, gm[2] / 2))
        hh = np.where(tier[k] == 0, dy[k], np.where(tier[k] == 1, gm[1] / 2, gm[2] / 2))

        ll = gpd.GeoSeries(gpd.points_from_xy(x[k], y[k]), crs=C.METRIC_CRS).to_crs(C.WGS84)
        olat = round(float(ll.y.mean()), 4); olon = round(float(ll.x.mean()), 4)
        out[str(int(mergid))] = {
            "o": [olat, olon],
            "t0": pd.Timestamp(t0).strftime("%Y-%m-%d %H:%M"),
            "dur_h": round(float(sec.max()) / 3600.0, 1),
            "pcb_h": [round(float(s) / 3600.0, 1) for s in ons_s],
            "pcb_ll": [[round(float(o["lat"]), 4), round(float(o["lon"]), 4)] for o in ons],
            "pcb_id": [str(o["id"]) for o in ons],
            "dlat": ((ll.y.to_numpy() - olat) * 1e5).round().astype(int).tolist(),
            "dlon": ((ll.x.to_numpy() - olon) * 1e5).round().astype(int).tolist(),
            "hw": np.round(hw).astype(int).tolist(),
            "hh": np.round(hh).astype(int).tolist(),
            "t":  np.round(sec[k] / 3600.0, 2).tolist(),
            "tier": tier[k].tolist(),
            "n_raw": int(len(g)),
        }
    return out
