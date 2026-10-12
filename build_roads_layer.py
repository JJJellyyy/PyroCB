"""Roads for the HTML map: OpenStreetMap roads and tracks within BUFFER_KM of each year's pyroCb-linked fire perimeters.

Reads the extracted lines (INPUT_DATA/ROADS/roads_{region}.parquet, from `python -m pcbcube.roads extract`) and the
perimeters of DATA[y].per in map_src/data.json; writes one lazily loaded file per year:

    roads/roads_{year}.json     FeatureCollection of 3 MultiLineStrings (lon/lat, 4 decimals ~ 8 m), simplified to 40 m:
                                  cls 2 = major roads (motorway ... secondary), 1 = other public roads, 0 = tracks

Only pyroCb-linked fires get roads on the map (every fire of every year would be hundreds of MB); the distance to the
nearest road exists in the dataset for every pixel and can be used for colouring (variable "Distance to nearest road").

    python build_roads_layer.py [year ...] [--all]      (--all: roads near every fire -> roads/roads_all_{year}.json)
"""
from __future__ import annotations

import base64
import glob
import gzip
import json
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
from shapely import STRtree

HERE = Path(__file__).resolve().parent
ROADS = HERE.parent / "INPUT_DATA" / "ROADS"
DATA = HERE / "map_src" / "data.json"
OUT = HERE / "roads"
BUFFER_KM = 10.0
SIMPLIFY_M = 40.0
BUFFER_KM_ALL = 5.0                    # all-fires mode: tighter corridor and a little more simplification (file size)
SIMPLIFY_M_ALL = 60.0
METRIC = 3571
NAMES = {2: "major", 1: "public", 0: "track"}


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def main(years=None, all_fires=False):
    """all_fires: roads within BUFFER_KM_ALL of EVERY fire (DATA[y].per + .perAll) -> roads/roads_all_{year}.json."""
    t0 = time.time()
    data = json.loads(DATA.read_text(encoding="utf-8"))
    years = [str(y) for y in (years or sorted(int(k) for k in data))]
    buf_km, simp_m, stem = (BUFFER_KM_ALL, SIMPLIFY_M_ALL, "roads_all") if all_fires else (BUFFER_KM, SIMPLIFY_M, "roads")
    # one buffered area per year, in metres
    areas = {}
    for y in years:
        feats = unpack(data[y]["per"])["features"] + (unpack(data[y]["perAll"])["features"] if all_fires and "perAll" in data[y] else [])
        per = gpd.GeoDataFrame.from_features(feats, crs=4326).to_crs(METRIC)
        areas[y] = shapely.union_all(shapely.buffer(shapely.make_valid(per.geometry.to_numpy()), buf_km * 1000.0, quad_segs=4))
    print(f"buffers ready ({time.time() - t0:.0f}s)", flush=True)
    found = {y: {0: [], 1: [], 2: []} for y in years}
    for fp in sorted(glob.glob(str(ROADS / "roads_*.parquet"))):
        region = Path(fp).stem.replace("roads_", "")
        g = gpd.read_parquet(fp)
        geoms = g.geometry.to_numpy()
        lonlat_bounds = g.geometry.bounds.to_numpy()
        for y in years:
            a = gpd.GeoSeries([areas[y]], crs=METRIC).to_crs(4326).iloc[0]
            minx, miny, maxx, maxy = a.bounds
            cand = np.flatnonzero((lonlat_bounds[:, 2] >= minx) & (lonlat_bounds[:, 0] <= maxx)
                                  & (lonlat_bounds[:, 3] >= miny) & (lonlat_bounds[:, 1] <= maxy))
            if not len(cand):
                continue
            shapely.prepare(a)
            sel = shapely.intersects(geoms[cand], a)
            idx = cand[sel]
            if len(idx):
                for c in (0, 1, 2):
                    found[y][c].append(geoms[idx[g["cls"].to_numpy()[idx] == c]])
        print(f"  {region}: {len(g):,} lines scanned ({time.time() - t0:.0f}s)", flush=True)
        del g, geoms
    OUT.mkdir(exist_ok=True)
    for y in years:
        feats = []
        n = 0
        for c in (2, 1, 0):
            parts = [p for p in found[y][c] if len(p)]
            if not parts:
                continue
            ls = np.concatenate(parts)
            n += len(ls)
            m = gpd.GeoSeries(ls, crs=4326).to_crs(METRIC)
            m = shapely.simplify(m.to_numpy(), simp_m)
            ll = gpd.GeoSeries(m, crs=METRIC).to_crs(4326).to_numpy()
            parts = [p for geom in ll for p in shapely.get_parts(geom)]
            ll = shapely.transform(shapely.multilinestrings(parts), lambda xy: np.round(xy, 4))
            feats.append({"type": "Feature", "properties": {"cls": c, "name": NAMES[c]},
                          "geometry": json.loads(shapely.to_geojson(ll))})
        (OUT / f"{stem}_{y}.json").write_text(json.dumps({"type": "FeatureCollection", "features": feats}, separators=(",", ":")),
                                             encoding="utf-8")
        print(f"  {y}: {n:,} road segments -> {stem}_{y}.json {(OUT / f'{stem}_{y}.json').stat().st_size / 1e6:.1f} MB", flush=True)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main([int(a) for a in args] or None, all_fires="--all" in sys.argv)
