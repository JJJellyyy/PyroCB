"""Growth bands split by BURNING COMPLEX for the HTML map (band option "Burning complexes").

For every pyroCb-linked fire in the map payload and every overpass: the burning complexes of that overpass (pcbcube.complexes rules:
pixels of the last 72 h with decay weights, linked when at most GAP_FACTOR undetected pixels lie between their footprints) that hold
new pixels of the fire. One feature per complex x overpass: the area of the complex's NEW pixels (detected at that overpass),
closed over the same gap (smoke) and simplified - so successive overpasses do not cover each other. Properties: mergid (the map
fire), step, doy, band_i / n_bands (overpass order within the fire), n_new, area_new_km2, frp_new_mw, area_w_km2 (decay-weighted
area of the whole complex), rg_km (spread), n_fire_ids, pcb (1 = the complex picked for a kept pyroCb at its pre-onset overpass,
matched to complex_panel by step and centroid), pcb_id, region, rg (pyroCb-region tag of the fire).
Output complex_bands/cb_{year}.json next to the map (fetched per year on demand).

    python build_complex_bands.py [years ...]
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
from pcbcube.atlas import StepIndex                                # noqa: E402
from pcbcube.complexes import GAP_FACTOR, WEIGHTS, components      # noqa: E402
from pcbcube.complex_panel import PANEL_FP                         # noqa: E402
from pcbcube.geometry import pixel_footprints                      # noqa: E402

DATA = HERE / "map_src" / "data.json"
OUTDIR = HERE / "complex_bands"
SIMPLIFY_M = 60.0


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def build_year(year, d, panel):
    t0 = time.time()
    st = StepIndex()
    bnd = unpack(d["bnd"])
    region, rgt = {}, {}
    for f in bnd["features"]:
        p = f["properties"]; region[int(p["mergid"])] = p.get("region"); rgt[int(p["mergid"])] = p.get("rg", 0)
    keys = {year * 10_000_000 + m for m in region}
    P = pq.read_table(C.OUT / "pixels.parquet", columns=["fire_key", "si", "x", "y", "DS", "DT", "frp", "pix_area_km2"],
                      filters=[("year", "=", year)]).to_pandas().sort_values("si", kind="stable").reset_index(drop=True)
    si = P.si.to_numpy(); x_all = P.x.to_numpy(float); y_all = P.y.to_numpy(float)
    hx = P.DS.to_numpy(float) * 500; hy = P.DT.to_numpy(float) * 500; fk = P.fire_key.to_numpy(np.int64)
    area = P.pix_area_km2.to_numpy(float); fmw = P.frp.to_numpy(float) * area
    mine = np.isin(fk, list(keys))
    pos = panel[panel.year == year]
    rows = []
    for s in np.unique(si[mine]):
        lo, hi = np.searchsorted(si, s - len(WEIGHTS) + 1, "left"), np.searchsorted(si, s, "right")
        ix = np.arange(lo, hi); age = s - si[ix]; new = age == 0
        lab = components(x_all[ix], y_all[ix], hx[ix], hy[ix], GAP_FACTOR)
        w = np.asarray(WEIGHTS)[age]
        hit = np.unique(lab[new & mine[ix]])
        ps = pos[pos.si == s]
        for k in hit:
            m = lab == k; mn = m & new
            q = ix[mn]
            g = pixel_footprints(x_all[q], y_all[q], hx[q] / 500, hy[q] / 500)
            r = GAP_FACTOR * float(np.median(np.maximum(hx[q], hy[q])))
            g = shapely.simplify(g.buffer(r, quad_segs=3).buffer(-r, quad_segs=3), SIMPLIFY_M, preserve_topology=True)
            if g.is_empty:
                continue
            aw = area[ix][m] * w[m]
            cxw = fmw[ix][m] * w[m]; cxw = cxw if cxw.sum() > 0 else aw
            cx, cy = np.average(x_all[ix][m], weights=cxw), np.average(y_all[ix][m], weights=cxw)
            rg = np.sqrt(np.average((x_all[ix][m] - cx) ** 2 + (y_all[ix][m] - cy) ** 2, weights=aw)) / 1e3
            fks = fk[q]; mids = [int(f - year * 10_000_000) for f in np.unique(fks) if f in keys]
            pcb, pid = 0, None
            if len(ps):
                dd = np.hypot(ps.cx.to_numpy() - cx, ps.cy.to_numpy() - cy)
                j = int(np.argmin(dd))
                if dd[j] < 2000:
                    pcb, pid = 1, ps.pcb_id.iloc[j]
            for mid in mids:                                   # one copy per map fire it belongs to
                rows.append(dict(mergid=mid, si=int(s), step=str(st.steps[s]), n_new=int(mn.sum()), area_new_km2=round(float(area[q].sum()), 2),
                                 frp_new_mw=round(float(fmw[q].sum()), 1), area_w_km2=round(float(aw.sum()), 2), rg_km=round(float(rg), 2),
                                 n_fire_ids=int(len(np.unique(fk[ix][m]))), pcb=pcb, pcb_id=pid, region=region.get(mid), rg=rgt.get(mid, 0),
                                 geometry=g))
    if not rows:
        return {"type": "FeatureCollection", "features": []}
    G = gpd.GeoDataFrame(rows, geometry="geometry", crs=3571).to_crs(4326)
    G["geometry"] = shapely.transform(G.geometry.values, lambda a: np.round(a, 4))
    G = G[~G.geometry.is_empty]
    G["doy"] = pd.to_datetime(G["step"].str[:8], format="%Y%m%d").dt.dayofyear
    G["band_i"] = G.groupby("mergid")["si"].rank(method="dense").astype(int) - 1
    G["n_bands"] = G.groupby("mergid")["band_i"].transform("max") + 1
    G = G.sort_values(["mergid", "si", "pcb"])
    print(f"  {year}: {len(G):,} complex bands over {G.mergid.nunique()} fires ({G.pcb.sum()} pyroCb complexes) in {time.time() - t0:.0f}s", flush=True)
    return json.loads(G.drop(columns=["si"]).to_json(drop_id=True))


def main(years=None):
    data = json.loads(DATA.read_text(encoding="utf-8"))
    panel = pd.read_parquet(PANEL_FP, columns=["year", "si", "cx", "cy", "cls", "pcb_id"])
    panel = panel[panel.cls == 1]
    OUTDIR.mkdir(exist_ok=True)
    for y in years or sorted(data, key=int):
        gj = build_year(int(y), data[str(y)], panel)
        raw = json.dumps(gj, separators=(",", ":")).encode()
        (OUTDIR / f"cb_{y}.json").write_bytes(raw)
        print(f"       complex_bands/cb_{y}.json: {len(raw) / 1e6:.1f} MB (~{len(gzip.compress(raw, 6)) / 1e6:.1f} MB gzipped)", flush=True)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
