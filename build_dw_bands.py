"""Growth bands of DISTANCE-WEIGHTED burning complexes (test) for the HTML map (band option "Distance-weighted complexes (test)").

Complexes exactly as now (pcbcube.complexes: everything burning in the last 72 h with decay weights, <= 2 undetected pixels between
footprints), plus a DISTANCE decay - symmetric, the same for every fire: each complex's hotspot = the new pixel where the Gaussian-
smoothed (SIGMA_KM) new-pixel power of the complex is highest; every pixel gets weight w_age x exp(-d^2 / 2 SIGMA^2), d = distance to
the hotspot. Far ends of large complexes count for (almost) nothing, without a hard cut.
PyroCb complex = the current pick (strongest pull, weighted power / max(d, 1 km)^2, pixels up to onset).
One feature per complex x overpass x distance class (core d <= SIGMA, middle SIGMA < d <= 2 SIGMA, outer > 2 SIGMA) of the complex's
NEW pixels (closed over the 2-pixel gap): opacity by class on the map. Properties: mergid, step, doy, band_i / n_bands, dclass,
n_new, area_new_km2 / area_new_d_km2 (distance-weighted), frp_new_mw / frp_new_d_mw, area_w_km2 / area_wd_km2 (whole complex,
time-weighted / time x distance-weighted), hot_d_km (distance of the pyroCb point to the hotspot, pyroCb complexes), pcb, pcb_id, region, rg.
Output dw_bands/dw_{year}.json next to the map.

    python build_dw_bands.py [years ...]          (default 2023)
"""
from __future__ import annotations

import base64
import gzip
import json
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pcbcube import config as C                                    # noqa: E402
from pcbcube import prepulse                                       # noqa: E402
from pcbcube.atlas import StepIndex                                # noqa: E402
from pcbcube.complexes import D_FLOOR_KM, GAP_FACTOR, WEIGHTS, components  # noqa: E402
from pcbcube.geometry import pixel_footprints                      # noqa: E402

DATA = HERE / "map_src" / "data.json"
OUTDIR = HERE / "dw_bands"
SIGMA_KM = 5.0
GRID_M = 1000.0
SIMPLIFY_M = 60.0
MARGIN_M = 80e3


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def hotspot(x, y, p, sigma_km=SIGMA_KM, grid_m=GRID_M):
    """Index of the point where the Gaussian-smoothed power (grid) is highest."""
    from scipy.ndimage import gaussian_filter
    if len(x) == 1 or p.sum() <= 0:
        return int(np.argmax(p)) if len(p) else 0
    x0, y0 = x.min(), y.min()
    ix = ((x - x0) // grid_m).astype(int); iy = ((y - y0) // grid_m).astype(int)
    g = np.zeros((iy.max() + 1, ix.max() + 1)); np.add.at(g, (iy, ix), p)
    g = gaussian_filter(g, sigma_km * 1000 / grid_m, mode="constant")
    return int(np.argmax(g[iy, ix]))


def build_year(year, d):
    t0 = time.time()
    st = StepIndex()
    region, rgt = {}, {}
    for f in unpack(d["bnd"])["features"]:
        p = f["properties"]; region[int(p["mergid"])] = p.get("region"); rgt[int(p["mergid"])] = p.get("rg", 0)
    keys = np.array(sorted(year * 10_000_000 + m for m in region), dtype=np.int64)
    P = pq.read_table(C.OUT / "pixels.parquet", columns=["fire_key", "si", "x", "y", "DS", "DT", "frp", "pix_area_km2", "datetime_utc"],
                      filters=[("year", "=", year)]).to_pandas()
    P["frp_mw"] = P.frp * P.pix_area_km2; P["hx"] = P.DS * 500; P["hy"] = P.DT * 500
    P["mine"] = np.isin(P.fire_key.to_numpy(np.int64), keys)
    by_si = {s: g for s, g in P.groupby("si")}
    ev = prepulse.events(prepulse.GAP_DAYS)
    ev = ev[ev.fire_key.notna() & (ev.onset_utc.dt.year == year)]
    ev_at = {int(s): g for s, g in ev.groupby("si_grid")}
    rows = []
    steps = sorted(P.loc[P.mine, "si"].unique())
    s2 = 2 * (SIGMA_KM * 1e3) ** 2
    for n_s, s in enumerate(steps):
        A = pd.concat([by_si[s - a].assign(age=a, w=WEIGHTS[a]) for a in range(len(WEIGHTS)) if (s - a) in by_si], ignore_index=True)
        M = A[A.mine & (A.age == 0)]
        if M.empty:
            continue
        A = A[(A.x >= M.x.min() - MARGIN_M) & (A.x <= M.x.max() + MARGIN_M) & (A.y >= M.y.min() - MARGIN_M) & (A.y <= M.y.max() + MARGIN_M)].reset_index(drop=True)
        lab = components(A.x.to_numpy(), A.y.to_numpy(), A.hx.to_numpy(), A.hy.to_numpy())
        pick = {}
        for _, e in ev_at.get(int(s), pd.DataFrame()).iterrows():
            ok = (A.datetime_utc <= e.onset_utc).to_numpy()
            dd = np.hypot(A.x - e.x, A.y - e.y).to_numpy() / 1e3
            cg = pd.DataFrame({"c": lab[ok], "d": dd[ok], "fw": (A.frp_mw * A.w).to_numpy()[ok]}).groupby("c").agg(d=("d", "min"), fw=("fw", "sum"))
            if len(cg):
                k = int((cg.fw / np.maximum(cg.d, D_FLOOR_KM) ** 2).idxmax())
                pick.setdefault(k, []).append((e.pcb_id, float(e.x), float(e.y)))
        new = (A.age == 0).to_numpy()
        for c in np.unique(lab[new & A.mine.to_numpy()]):
            m = lab == c
            C_ = A[m]; nw = (C_.age == 0).to_numpy()
            h = hotspot(C_.x.to_numpy()[nw], C_.y.to_numpy()[nw], C_.frp_mw.to_numpy()[nw])
            hx0, hy0 = C_.x.to_numpy()[nw][h], C_.y.to_numpy()[nw][h]
            dist = np.hypot(C_.x - hx0, C_.y - hy0).to_numpy()
            kd = np.exp(-dist ** 2 / s2)
            aw = (C_.pix_area_km2 * C_.w).to_numpy()
            tot = dict(area_w_km2=round(float(aw.sum()), 2), area_wd_km2=round(float((aw * kd).sum()), 2))
            pc = pick.get(int(c))
            hot_d = round(float(min(np.hypot(px - hx0, py - hy0) for _, px, py in pc) / 1e3), 1) if pc else None
            Q = C_[nw]; dq = dist[nw]; kq = kd[nw]
            cls = np.where(dq <= SIGMA_KM * 1e3, "core", np.where(dq <= 2 * SIGMA_KM * 1e3, "middle", "outer"))
            r = GAP_FACTOR * float(np.median(np.maximum(Q.hx, Q.hy)))
            for dc in ("outer", "middle", "core"):
                mm = cls == dc
                if not mm.any():
                    continue
                q = Q[mm]
                g = pixel_footprints(q.x, q.y, q.DS, q.DT)
                g = shapely.simplify(g.buffer(r, quad_segs=3).buffer(-r, quad_segs=3), SIMPLIFY_M, preserve_topology=True)
                if g.is_empty:
                    continue
                for fk in np.unique(Q.fire_key[Q.mine]):
                    mid = int(fk - year * 10_000_000)
                    rows.append(dict(mergid=mid, si=int(s), step=str(st.steps[s]), dclass=dc, n_new=int(nw.sum()),
                                     area_new_km2=round(float(Q.pix_area_km2.sum()), 2), area_new_d_km2=round(float((Q.pix_area_km2 * kq).sum()), 2),
                                     frp_new_mw=round(float(Q.frp_mw.sum()), 1), frp_new_d_mw=round(float((Q.frp_mw * kq).sum()), 1), **tot,
                                     hot_d_km=hot_d, pcb=int(pc is not None), pcb_id=";".join(p for p, _, _ in pc) if pc else None,
                                     region=region.get(mid), rg=rgt.get(mid, 0), geometry=g))
        if n_s % 50 == 0:
            print(f"    {year}: step {n_s + 1}/{len(steps)} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    G = gpd.GeoDataFrame(rows, geometry="geometry", crs=3571).to_crs(4326)
    G["geometry"] = shapely.transform(G.geometry.values, lambda a: np.round(a, 4))
    G = G[~G.geometry.is_empty]
    G["doy"] = pd.to_datetime(G["step"].str[:8], format="%Y%m%d").dt.dayofyear
    G["band_i"] = G.groupby("mergid")["si"].rank(method="dense").astype(int) - 1
    G["n_bands"] = G.groupby("mergid")["band_i"].transform("max") + 1
    order = {"outer": 0, "middle": 1, "core": 2}
    G = G.assign(_o=G.dclass.map(order)).sort_values(["mergid", "si", "_o"]).drop(columns="_o")
    pc = G[(G.pcb == 1)].drop_duplicates(["pcb_id"])
    print(f"  {year}: {len(G):,} features over {G.mergid.nunique()} fires; pyroCb complexes {len(pc)}; pyroCb point to hotspot median "
          f"{pc.hot_d_km.median():.1f} km; distance-weighted / unweighted area of pyroCb complexes median "
          f"{(pc.area_wd_km2 / pc.area_w_km2).median():.2f} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return json.loads(G.drop(columns=["si"]).to_json(drop_id=True))


def main(years=None):
    data = json.loads(DATA.read_text(encoding="utf-8"))
    OUTDIR.mkdir(exist_ok=True)
    for y in years or [2023]:
        gj = build_year(int(y), data[str(y)])
        raw = json.dumps(gj, separators=(",", ":")).encode()
        (OUTDIR / f"dw_{y}.json").write_bytes(raw)
        print(f"       dw_bands/dw_{y}.json: {len(raw) / 1e6:.1f} MB (~{len(gzip.compress(raw, 6)) / 1e6:.1f} MB gzipped)", flush=True)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
