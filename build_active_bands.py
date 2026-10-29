"""12 h ACTIVE-AREA bands for the HTML map: one band per fire per overpass = the 12 h growth ring (pcbcube.rings) united with
the footprints (DS x DT) of every new VIIRS pixel the fire had at that overpass - inside or outside the ring.

The plain 12 h polygon bands (build_snapshot_bands.py) only show how the atlas perimeter grew (and merge quiet steps); about half
of each overpass's new pixels lie inside the old perimeter (edge margin of the alpha-shape hull, interior islands burning later,
late first detections). These bands show where the fire was actually active at each overpass: no merging, and bands of later
overpasses may overlap earlier ones (drawn on top). Same fires as the other band options (the pyroCb-linked fires in the payload).

Display geometry: the footprint union is closed by CLOSE_M (+/- buffer: far fewer vertices than raw pixel boxes), united with
the ring and simplified; the areas in the properties are exact (before smoothing). Writes active_bands/ab_{year}.json next to
the map (GeoJSON, WGS84), fetched on demand when the "12 h active area" band option is chosen - not in data.json, so the map
file does not grow. Properties: mergid, band_i / n_bands (time order), step, doy, area_km2 (ring u footprints), ring_km2,
fp_km2 (dissolved footprints), n_pix, frp_mw, cum_area_km2, region, rg. No map_tool join needed; deploy the folder with the map.

    python build_active_bands.py            # all years
    python build_active_bands.py 2023       # one year (size check)
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
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
from pcbcube.geometry import pixel_footprints                      # noqa: E402
from pcbcube.rings import rings_path                               # noqa: E402

DATA = HERE / "map_src" / "data.json"
MIN_TOL_M, MAX_TOL_M = 90.0, 200.0
DECIMALS = 4                   # ~11 m
CLOSE_M = 150.0                # display smoothing of the footprint union
OUTDIR = HERE / "active_bands"


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def pack(o):
    raw = json.dumps(o, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode(), len(raw)


def _fire(args):
    """All active-area bands of one fire: list of dicts with WKB geometry."""
    mid, R, P = args                                     # R: rings (si, wkb, ring_km2, cum_km2); P: pixels (si, x, y, DS, DT, frp_mw)
    out = []
    for si in sorted(set(R) | set(P)):
        parts, disp = [], []
        ring_km2 = cum = 0.0
        if si in R:
            g, ring_km2, cum = R[si]
            parts.append(shapely.from_wkb(g)); disp.append(parts[-1])
        fp_km2, n, frp = 0.0, 0, 0.0
        if si in P:
            x, y, ds, dt, f = P[si]
            fp = pixel_footprints(x, y, ds, dt)
            if fp is not None and not fp.is_empty:
                parts.append(fp); fp_km2 = fp.area / 1e6
                disp.append(fp.buffer(CLOSE_M, quad_segs=2).buffer(-CLOSE_M, quad_segs=2))
            n, frp = len(x), float(np.nansum(f))
        if not parts:
            continue
        geom = shapely.union_all(parts)
        if geom.is_empty:
            continue
        shown = shapely.union_all(disp)
        out.append(dict(si=int(si), wkb=shapely.to_wkb(geom if shown.is_empty else shown), area_km2=round(geom.area / 1e6, 3),
                        ring_km2=round(ring_km2, 3),
                        fp_km2=round(fp_km2, 3), n_pix=n, frp_mw=round(frp, 1), cum_km2=cum))
    return mid, out


def build_year(year: int, d: dict, workers: int):
    t0 = time.time()
    st = StepIndex()
    bnd = unpack(d["bnd"])
    region, rg = {}, {}
    for f in bnd["features"]:
        p = f["properties"]; region[int(p["mergid"])] = p.get("region"); rg[int(p["mergid"])] = p.get("rg")
    if "sb" in d:
        for f in unpack(d["sb"])["features"]:
            p = f["properties"]
            if p.get("rg") is not None:
                rg.setdefault(int(p["mergid"]), p["rg"])
    ids = sorted(region)
    if not ids:
        return {"type": "FeatureCollection", "features": []}
    keys = [year * 10_000_000 + m for m in ids]
    G = gpd.read_parquet(rings_path(year))
    G = G[G["fire_key"].isin(keys) & (G["ring_area_km2"] > 0)]
    Pp = pq.read_table(C.OUT / "pixels.parquet", columns=["fire_key", "si", "x", "y", "DS", "DT", "frp", "pix_area_km2"],
                       filters=[("year", "=", year)]).to_pandas()
    Pp = Pp[Pp["fire_key"].isin(keys)]
    Pp["frp_mw"] = Pp["frp"] * Pp["pix_area_km2"]
    jobs = []
    for k in keys:
        R = {int(r.si): (shapely.to_wkb(r.geometry), float(r.ring_area_km2), float(r.cum_area_km2)) for r in G[G.fire_key == k].itertuples()}
        P = {int(si): (g.x.to_numpy(), g.y.to_numpy(), g.DS.to_numpy(), g.DT.to_numpy(), g.frp_mw.to_numpy())
             for si, g in Pp[Pp.fire_key == k].groupby("si")}
        if R or P:
            jobs.append((k - year * 10_000_000, R, P))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_fire, jobs, chunksize=4))
    rows = []
    for mid, bands in res:
        if not bands:
            continue
        cum_max = max([b["cum_km2"] for b in bands] + [0.0])
        tol = float(np.clip(0.0012 * np.sqrt(max(cum_max, 1.0) * 1e6), MIN_TOL_M, MAX_TOL_M))
        cum = 0.0
        for j, b in enumerate(bands):
            cum = max(cum, b["cum_km2"])
            geom = shapely.simplify(shapely.from_wkb(b["wkb"]), tol, preserve_topology=True)
            if geom.is_empty:
                continue
            step = str(st.steps[b["si"]])
            rows.append(dict(mergid=mid, band_i=j, n_bands=len(bands), step=step, step_hi=step, area_km2=b["area_km2"],
                             ring_km2=b["ring_km2"], fp_km2=b["fp_km2"], n_pix=b["n_pix"], frp_mw=b["frp_mw"],
                             cum_area_km2=round(cum, 2), region=region.get(mid), rg=rg.get(mid, 0), geometry=geom))
    if not rows:
        return {"type": "FeatureCollection", "features": []}
    g = gpd.GeoDataFrame(rows, geometry="geometry", crs=3571).to_crs(4326)
    g["geometry"] = shapely.transform(g.geometry.values, lambda a: np.round(a, DECIMALS))
    g = g[~g.geometry.is_empty]
    g["doy"] = pd.to_datetime(g["step"].str[:8], format="%Y%m%d", errors="coerce").dt.dayofyear
    g = g.sort_values(["mergid", "band_i"])
    g["band_i"] = g.groupby("mergid").cumcount()
    g["n_bands"] = g.groupby("mergid")["band_i"].transform("size")
    print(f"  {year}: {len(g):,} active-area bands over {g.mergid.nunique()} fires (median {int(g.groupby('mergid').size().median())}/fire; "
          f"{(g.fp_km2 > 0).mean():.0%} with pixels, {(g.ring_km2 > 0).mean():.0%} with perimeter growth) in {time.time() - t0:.0f}s", flush=True)
    return json.loads(g.drop(columns=["step_hi"]).to_json(drop_id=True))


def main(years=None):
    data = json.loads(DATA.read_text(encoding="utf-8"))
    years = years or sorted(data, key=int)
    workers = max(1, (os.cpu_count() or 2) - 2)
    OUTDIR.mkdir(exist_ok=True)
    tot = gz = 0.0
    for y in years:
        gj = build_year(int(y), data[str(y)], workers)
        raw = json.dumps(gj, separators=(",", ":")).encode()
        (OUTDIR / f"ab_{y}.json").write_bytes(raw)
        z = len(gzip.compress(raw, 6)) / 1e6
        tot += len(raw) / 1e6; gz += z
        print(f"       active_bands/ab_{y}.json: {len(raw) / 1e6:.1f} MB (~{z:.1f} MB gzipped over http)", flush=True)
    print(f"wrote {len(years)} files: {tot:.0f} MB raw, ~{gz:.0f} MB transferred")


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
