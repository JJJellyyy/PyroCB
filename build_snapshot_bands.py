"""Growth bands from the 12 h snapshots as REAL polygon rings (not 375 m burn-age cells).

For every fire that has a pyroCb band set in the map payload, walk its Snapshot polygons in time:

    S_k    = S_{k-1} u V(snap_k)          running union (monotone, disjoint rings by construction,
                                           same construction as pyrocb/bands.py)
    band_j = S_k \\ S_last_kept            one ring per kept snapshot

A snapshot only becomes its own band when the fire grew by at least
max(MIN_GROW_KM2, FRAC_OF_FINAL * final area) since the last kept band - quiet 12 h steps are
merged into the next band so a three-month fire does not become 150 hairline rings. The first
band is the footprint at the fire's first snapshot; the last band ends at its last snapshot.

Writes DATA[year]["sb"] (gzip+base64 GeoJSON, WGS84, same property names as the pyroCb bands)
into map_src/data.json. Then:  python map_tool.py join

    python build_snapshot_bands.py            # all years
    python build_snapshot_bands.py 2023       # one year (size check)
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

import numpy as np
import pandas as pd
import geopandas as gpd
import pyogrio
import shapely

from pyrocb import config as C
from pyrocb import io

DATA = "map_src/data.json"
MIN_GROW_KM2 = 0.5
FRAC_OF_FINAL = 0.012
COORD_DECIMALS = 5          # ~1 m


def _read_step(args):
    year, step, ids = args
    fp = os.path.join(C.SNAP_DIR.format(y=year), f"{step}.gpkg")
    try:
        g = pyogrio.read_dataframe(fp, columns=["mergid"], where=f"mergid IN ({','.join(map(str, ids))})")
    except Exception:                                    # noqa: BLE001  (empty / unreadable step)
        return step, {}
    if g.empty:
        return step, {}
    if g.crs is None or g.crs.to_epsg() != C.METRIC_CRS:
        g = g.to_crs(C.METRIC_CRS)
    geoms = io.ensure_valid(g.geometry.values)
    out = {}
    for mid, idx in pd.Series(range(len(g))).groupby(g["mergid"].to_numpy()):
        out[int(mid)] = shapely.to_wkb(shapely.union_all(geoms[idx.to_numpy()]))
    return step, out


def _build_fire(args):
    mid, seq = args                                      # seq: [(step, wkb)] sorted by time
    steps = [s for s, _ in seq]
    geoms = io.snap([shapely.from_wkb(w) for _, w in seq])
    final_km2 = max(io.km2(g) for g in geoms)
    thresh = max(MIN_GROW_KM2, FRAC_OF_FINAL * final_km2)
    cum = kept = None
    last_step, rows = "", []
    for i, (s, g) in enumerate(zip(steps, geoms)):
        new = g if cum is None else shapely.union_all([cum, g])
        grow = io.km2(new) if kept is None else io.km2(shapely.difference(new, kept))
        if kept is None or grow >= thresh or i == len(steps) - 1:
            band = new if kept is None else shapely.difference(new, kept)
            band, _ = io.drop_slivers(band)
            if band is not None and not band.is_empty:
                rows.append({"step_lo": last_step, "step_hi": s, "i_hi": i,
                             "area_km2": round(io.km2(band), 4), "wkb": shapely.to_wkb(band)})
                last_step = s
            kept = new
        cum = new
    return mid, rows, io.km2(cum), final_km2


def _round(g, nd=COORD_DECIMALS):
    return shapely.transform(g, lambda a: np.round(a, nd))


def build_year(year, d, workers):
    t0 = time.time()
    bnd = json.loads(gzip.decompress(base64.b64decode(d["bnd"])))
    region = {}
    for f in bnd["features"]:
        region[int(f["properties"]["mergid"])] = f["properties"].get("region")
    ids = sorted(region)
    if not ids:
        return {"type": "FeatureCollection", "features": []}, 0
    steps = list(io.list_steps(year))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        per_step = list(ex.map(_read_step, [(year, s, ids) for s in steps], chunksize=8))
    seqs = {}
    for s, got in per_step:                                  # steps are in time order
        for mid, w in got.items():
            seqs.setdefault(mid, []).append((s, w))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_build_fire, sorted(seqs.items()), chunksize=1))

    rows = []
    for mid, bands, cum_km2, final_km2 in res:
        n = len(bands)
        tol = float(np.clip(0.0012 * np.sqrt(final_km2 * 1e6), C.BAND_SIMPLIFY_M, 150.0))
        for j, b in enumerate(bands):
            geom = shapely.simplify(shapely.from_wkb(b["wkb"]), tol, preserve_topology=True)
            if geom.is_empty:
                continue
            rows.append(dict(mergid=mid, year=year, band_i=j, n_bands=n, area_km2=b["area_km2"],
                             cum_area_km2=round(cum_km2, 4), step_lo=b["step_lo"], step_hi=b["step_hi"],
                             is_post_final_band=False, qc_regress_km2=0.0,
                             region=region.get(mid), geometry=geom))
    if not rows:
        return {"type": "FeatureCollection", "features": []}, 0
    g = gpd.GeoDataFrame(rows, geometry="geometry", crs=C.METRIC_CRS).to_crs(C.WGS84)
    g["geometry"] = [_round(x) for x in g.geometry.values]
    g = g[~g.geometry.is_empty]
    g["doy"] = pd.to_datetime(g["step_hi"].str[:8], format="%Y%m%d", errors="coerce").dt.dayofyear
    # band_i / n_bands only count bands that survived simplification
    g = g.sort_values(["mergid", "band_i"])
    g["band_i"] = g.groupby("mergid").cumcount()
    g["n_bands"] = g.groupby("mergid")["band_i"].transform("size")
    print(f"  {year}: {len(g)} bands over {g.mergid.nunique()} fires "
          f"(median {int(g.groupby('mergid').size().median())}/fire) in {time.time() - t0:.0f}s", flush=True)
    return json.loads(g.to_json(drop_id=True)), len(g)


def pack(gj):
    raw = json.dumps(gj, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode(), len(raw)


def main(years=None):
    with open(DATA, encoding="utf-8") as f:
        data = json.load(f)
    years = years or sorted(data)
    workers = max(1, (os.cpu_count() or 2) - 2)
    total_mb = 0.0
    for y in years:
        gj, n = build_year(int(y), data[str(y)], workers)
        packed, raw_len = pack(gj)
        data[str(y)]["sb"] = packed
        total_mb += len(packed) / 1e6
        print(f"       payload sb[{y}]: {len(packed) / 1e6:.2f} MB (gz+b64), raw {raw_len / 1e6:.1f} MB", flush=True)
    bak = DATA + ".pre_sb"
    if not os.path.exists(bak):
        shutil.copyfile(DATA, bak)
    with open(DATA, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, separators=(",", ":"))
    print(f"added sb payload for {len(years)} years: {total_mb:.1f} MB; data.json now "
          f"{os.path.getsize(DATA) / 1e6:.1f} MB (backup: {bak})")


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
