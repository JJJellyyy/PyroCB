"""Growth bands split by SPREAD-DIRECTION FRONT (pcbcube.fronts, test method v2) for the HTML map (band option "Spread fronts (test)").

Same set-up as build_complex_bands.py, but the unit is a front: per overpass the burning pixels around the pyroCb-linked fires (last
72 h, decay weights) are split by pcbcube.fronts.segment (direction smoothed over 5 km, fronts within 60 deg of their mean direction,
small fronts merged). One feature per front x overpass holding new pixels of a map fire: the area of the front's NEW pixels, closed
over the 2-pixel gap and simplified. PyroCb fronts: at each matched pyroCb's pre-onset overpass, pcbcube.fronts.assign_pcb on the
pixels detected up to onset -> the front with the strongest pull (role "main") and every front with >= 25 % of that pull ("co-feeding").
Properties: mergid, step, doy, band_i / n_bands, n_new, area_new_km2, frp_new_mw, area_w_km2 (front incl. attached decaying pixels),
dir_deg (power-weighted mean spread direction, compass, 0 = north), n_fronts (fronts of the fire at that overpass), pcb (0/1),
pcb_role, pcb_id, pcb_kept, n_cofeed, region, rg.
Output front_bands/fb_{year}.json next to the map (fetched per year on demand).

    python build_front_bands.py [years ...]          (default 2023)
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
from pcbcube.complexes import GAP_FACTOR, WEIGHTS                  # noqa: E402
from pcbcube.fronts import assign_pcb, segment                     # noqa: E402
from pcbcube.geometry import pixel_footprints                      # noqa: E402

DATA = HERE / "map_src" / "data.json"
OUTDIR = HERE / "front_bands"
SIMPLIFY_M = 60.0
MARGIN_M = 80e3


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def build_year(year, d):
    t0 = time.time()
    st = StepIndex()
    region, rgt = {}, {}
    for f in unpack(d["bnd"])["features"]:
        p = f["properties"]; region[int(p["mergid"])] = p.get("region"); rgt[int(p["mergid"])] = p.get("rg", 0)
    keys = np.array(sorted(year * 10_000_000 + m for m in region), dtype=np.int64)
    P = pq.read_table(C.OUT / "pixels.parquet", columns=["fire_key", "si", "x", "y", "DS", "DT", "frp", "pix_area_km2", "datetime_utc",
                                                          "spread_ux", "spread_uy"], filters=[("year", "=", year)]).to_pandas()
    P["frp_mw"] = P.frp * P.pix_area_km2; P["hx"] = P.DS * 500; P["hy"] = P.DT * 500
    P["mine"] = np.isin(P.fire_key.to_numpy(np.int64), keys)
    by_si = {s: g for s, g in P.groupby("si")}
    ev = prepulse.events(prepulse.GAP_DAYS)
    ev = ev[ev.fire_key.notna() & (ev.onset_utc.dt.year == year)]
    ev_at = {int(s): g for s, g in ev.groupby("si_grid")}
    rows = []
    steps = sorted(P.loc[P.mine, "si"].unique())
    for n_s, s in enumerate(steps):
        parts = [by_si[s - a].assign(age=a, w=WEIGHTS[a]) for a in range(len(WEIGHTS)) if (s - a) in by_si]
        A = pd.concat(parts, ignore_index=True)
        M = A[A.mine & (A.age == 0)]
        if M.empty:
            continue
        A = A[(A.x >= M.x.min() - MARGIN_M) & (A.x <= M.x.max() + MARGIN_M) & (A.y >= M.y.min() - MARGIN_M) & (A.y <= M.y.max() + MARGIN_M)].reset_index(drop=True)
        lab = segment(A)
        role = {}
        for _, e in ev_at.get(int(s), pd.DataFrame()).iterrows():
            ok = (A.datetime_utc <= e.onset_utc).to_numpy()
            ids, _ = assign_pcb(A[ok].reset_index(drop=True), lab[ok], float(e.x), float(e.y))
            for i, f in enumerate(ids):
                r = role.setdefault(f, {"ids": [], "role": "co-feeding", "kept": False, "n": len(ids)})
                r["ids"].append(e.pcb_id); r["kept"] |= bool(e.kept); r["n"] = max(r["n"], len(ids))
                if i == 0:
                    r["role"] = "main"
        new = (A.age == 0).to_numpy()
        hit = np.unique(lab[new & A.mine.to_numpy() & (lab >= 0)])
        nfr = {}
        for f in hit:
            for fk in np.unique(A.fire_key[(lab == f) & new & A.mine.to_numpy()]):
                nfr[fk] = nfr.get(fk, 0) + 1
        for f in hit:
            m = lab == f; mn = m & new
            q = A[mn]
            r = GAP_FACTOR * float(np.median(np.maximum(q.hx, q.hy)))
            g = pixel_footprints(q.x, q.y, q.DS, q.DT)
            g = shapely.simplify(g.buffer(r, quad_segs=3).buffer(-r, quad_segs=3), SIMPLIFY_M, preserve_topology=True)
            if g.is_empty:
                continue
            rr = A[m & (A.age <= 1).to_numpy() & A.spread_ux.notna().to_numpy()]
            if len(rr):
                wv = np.maximum(rr.frp_mw, 1e-3); nn = np.hypot(rr.spread_ux, rr.spread_uy).replace(0, np.nan)
                dx = np.nansum(wv * rr.spread_ux / nn); dy = np.nansum(wv * rr.spread_uy / nn)
                dir_deg = float((np.degrees(np.arctan2(dx, dy)) + 360) % 360)      # x = east, y = north in EPSG:3571 (approx.)
            else:
                dir_deg = np.nan
            info = role.get(int(f))
            for fk in np.unique(q.fire_key[q.mine]):
                mid = int(fk - year * 10_000_000)
                rows.append(dict(mergid=mid, si=int(s), step=str(st.steps[s]), n_new=int(mn.sum()), area_new_km2=round(float(q.pix_area_km2.sum()), 2),
                                 frp_new_mw=round(float(q.frp_mw.sum()), 1), area_w_km2=round(float((A.pix_area_km2[m] * A.w[m]).sum()), 2),
                                 dir_deg=None if np.isnan(dir_deg) else round(dir_deg), n_fronts=int(nfr.get(fk, 1)),
                                 pcb=int(info is not None), pcb_role=info["role"] if info else None, pcb_id=";".join(info["ids"]) if info else None,
                                 pcb_kept=bool(info["kept"]) if info else None, n_cofeed=int(info["n"]) if info else None,
                                 region=region.get(mid), rg=rgt.get(mid, 0), geometry=g))
        if n_s % 50 == 0:
            print(f"    {year}: step {n_s + 1}/{len(steps)} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    G = gpd.GeoDataFrame(rows, geometry="geometry", crs=3571).to_crs(4326)
    G["geometry"] = shapely.transform(G.geometry.values, lambda a: np.round(a, 4))
    G = G[~G.geometry.is_empty]
    G["doy"] = pd.to_datetime(G["step"].str[:8], format="%Y%m%d").dt.dayofyear
    G["band_i"] = G.groupby("mergid")["si"].rank(method="dense").astype(int) - 1
    G["n_bands"] = G.groupby("mergid")["band_i"].transform("max") + 1
    G = G.sort_values(["mergid", "si", "pcb"])
    pc = G[G.pcb == 1]
    print(f"  {year}: {len(G):,} front bands over {G.mergid.nunique()} fires; pyroCb fronts {len(pc)} "
          f"(main {int((pc.pcb_role == 'main').sum())}, co-feeding {int((pc.pcb_role == 'co-feeding').sum())}) in {(time.time() - t0) / 60:.1f} min", flush=True)
    return json.loads(G.drop(columns=["si"]).to_json(drop_id=True))


def main(years=None):
    data = json.loads(DATA.read_text(encoding="utf-8"))
    OUTDIR.mkdir(exist_ok=True)
    for y in years or [2023]:
        gj = build_year(int(y), data[str(y)])
        raw = json.dumps(gj, separators=(",", ":")).encode()
        (OUTDIR / f"fb_{y}.json").write_bytes(raw)
        print(f"       front_bands/fb_{y}.json: {len(raw) / 1e6:.1f} MB (~{len(gzip.compress(raw, 6)) / 1e6:.1f} MB gzipped)", flush=True)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
