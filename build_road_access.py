"""Road access of every fire for the HTML map (map only, not the dataset).

road_km = the minimum distance (km, EPSG:3571) between a fire's final perimeter (DATA[y].per and .perAll in
map_src/data.json) and the nearest OpenStreetMap road of a public class OR a service road:
    motorway, trunk, primary, secondary, tertiary (+ _link), unclassified, residential, living_street, road, service
(no tracks, paths, footways). 0 when a road touches, crosses or lies inside the perimeter. Capped at CAP_KM.
OSM is today's road network, not the network of the fire year.

Method: road lines and perimeter boundaries densified to <= DENSIFY_M in EPSG:3571; nearest road vertex to any boundary
vertex (scipy cKDTree) = boundary distance (error < ~DENSIFY_M / 2); fires whose polygon intersects a road line
(shapely STRtree, exact) get 0. Done region by region (one OSM extract at a time) and the minimum taken over regions.

Inputs   INPUT_DATA/ROADS/roads_{region}.parquet          public roads (cls >= 1) from `python -m pcbcube.roads extract`
         INPUT_DATA/ROADS/roads_service_{region}.parquet  highway=service, written by the `extract` step here
Writes   cache/road_access.tsv                            year, mergid, road_km (all fires)
         map_src/data.json                                road_km on per / perAll, propagated by mergid to pcb / ign / ignAll

    python build_road_access.py extract [region ...]     (service roads from the .osm.pbf, ~25 min for all regions)
    python build_road_access.py compute                  (distances -> cache/road_access.tsv)
    python build_road_access.py apply                    (tsv -> data.json)        then   python map_tool.py join
"""
from __future__ import annotations

import base64
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROADS = HERE.parent / "INPUT_DATA" / "ROADS"
DATA = HERE / "map_src" / "data.json"
TSV = HERE / "cache" / "road_access.tsv"
REGIONS = ["alaska", "canada", "russia", "mongolia", "kazakhstan", "china", "finland", "sweden", "norway"]
METRIC = 3571
DENSIFY_M = 100.0
CAP_KM = 100.0


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def pack(o):
    return base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def load_fires():
    """All final perimeters of all years as one GeoDataFrame (year, mergid, key, geometry lon/lat)."""
    import geopandas as gpd
    data = json.loads(DATA.read_text(encoding="utf-8"))
    parts = []
    for y in sorted(data, key=int):
        for key in ("per", "perAll"):
            if key in data[y]:
                feats = unpack(data[y][key])["features"]
                if feats:
                    g = gpd.GeoDataFrame.from_features(feats, crs=4326)[["mergid", "geometry"]]
                    g["year"], g["key"] = int(y), key
                    parts.append(g)
    return pd.concat(parts, ignore_index=True)


# ---- extract -----------------------------------------------------------------------------------------------
def corridor_mask(fires, cell=0.5, lon_pad=4, lat_pad=2):
    """Cells (lon x lat, `cell` deg) within ~CAP_KM of any fire perimeter bbox (same padding as pcbcube.roads)."""
    from scipy.ndimage import binary_dilation
    nx, ny = int(360 / cell), int(180 / cell)
    occ = np.zeros((nx, ny), bool)
    b = fires.geometry.bounds.to_numpy()
    i0 = np.clip(((b[:, 0] + 180) / cell).astype(int), 0, nx - 1)
    i1 = np.clip(((b[:, 2] + 180) / cell).astype(int), 0, nx - 1)
    j0 = np.clip(((b[:, 1] + 90) / cell).astype(int), 0, ny - 1)
    j1 = np.clip(((b[:, 3] + 90) / cell).astype(int), 0, ny - 1)
    for a, bb, c, d in zip(i0, i1, j0, j1):
        if bb - a > nx / 2:                              # wraps the antimeridian: mark both ends only
            occ[a, c:d + 1] = occ[bb, c:d + 1] = True
        else:
            occ[a:bb + 1, c:d + 1] = True
    return binary_dilation(occ, structure=np.ones((2 * lon_pad + 1, 2 * lat_pad + 1), bool)), cell


def extract(regions=None):
    import pyogrio
    mask, cell = corridor_mask(load_fires())
    for r in regions or REGIONS:
        out = ROADS / f"roads_service_{r}.parquet"
        pbf = ROADS / f"{r}.osm.pbf"
        if out.exists() or not pbf.exists():
            print(f"  {r}: {'have' if out.exists() else 'no pbf'}", flush=True)
            continue
        t0 = time.time()
        g = pyogrio.read_dataframe(pbf, layer="lines", where="highway = 'service'", columns=["highway"], use_arrow=False)
        b = g.geometry.bounds.to_numpy()
        i = np.clip(((b[:, [0, 2]] + 180) / cell).astype(int), 0, mask.shape[0] - 1)
        j = np.clip(((b[:, [1, 3]] + 90) / cell).astype(int), 0, mask.shape[1] - 1)
        keep = mask[i[:, 0], j[:, 0]] | mask[i[:, 1], j[:, 1]] | mask[i[:, 0], j[:, 1]] | mask[i[:, 1], j[:, 0]]
        n = len(g)
        g = g[keep][["geometry"]].copy()
        tmp = out.with_suffix(".part")
        g.to_parquet(tmp)
        tmp.replace(out)
        print(f"  {r}: {len(g):,} service ways in the corridor (of {n:,} read) ({(time.time() - t0) / 60:.1f} min)", flush=True)
        del g


# ---- compute -----------------------------------------------------------------------------------------------
def region_roads(r):
    """Public + service road lines of one region in EPSG:3571 (numpy array of geometries), or None."""
    import geopandas as gpd
    import shapely
    geoms = []
    fp = ROADS / f"roads_{r}.parquet"
    if fp.exists():
        g = gpd.read_parquet(fp)
        geoms.append(g.loc[g["cls"] >= 1].to_crs(METRIC).geometry.to_numpy())
    fs = ROADS / f"roads_service_{r}.parquet"
    if fs.exists():
        geoms.append(gpd.read_parquet(fs).to_crs(METRIC).geometry.to_numpy())
    else:
        print(f"  WARNING: no service roads for {r}", flush=True)
    if not geoms:
        return None
    a = np.concatenate(geoms)
    return a[~shapely.is_empty(a)]


def compute():
    import shapely
    from scipy.spatial import cKDTree
    from shapely import STRtree
    t0 = time.time()
    fires = load_fires()
    polys = shapely.make_valid(fires.to_crs(METRIC).geometry.to_numpy())
    bxy, bidx = shapely.get_coordinates(shapely.segmentize(shapely.boundary(polys), DENSIFY_M), return_index=True)
    print(f"{len(fires):,} fires, {len(bxy):,} boundary vertices ({time.time() - t0:.0f}s)", flush=True)
    cap = CAP_KM * 1000.0
    best = np.full(len(fires), cap)
    inside = np.zeros(len(fires), bool)
    for r in REGIONS:
        t1 = time.time()
        lines = region_roads(r)
        if lines is None or not len(lines):
            continue
        # exact: a road line touching / crossing / inside the perimeter
        fi, _ = STRtree(lines).query(polys, predicate="intersects")
        inside[np.unique(fi)] = True
        # nearest road vertex to the densified boundary
        rxy = shapely.get_coordinates(shapely.segmentize(lines, DENSIFY_M))
        del lines
        tree = cKDTree(rxy)
        d, _ = tree.query(bxy, k=1, distance_upper_bound=cap, workers=-1)
        d = np.minimum(d, cap)
        m = np.full(len(fires), cap)
        np.minimum.at(m, bidx, d)
        best = np.minimum(best, m)
        print(f"  {r}: {len(rxy):,} road vertices, {len(np.unique(fi)):,} fires touched; "
              f"{(time.time() - t1):.0f}s", flush=True)
        del rxy, tree, d
    km = np.where(inside, 0.0, np.round(best / 1000.0, 1))
    out = pd.DataFrame({"year": fires["year"], "mergid": fires["mergid"], "key": fires["key"], "road_km": km})
    TSV.parent.mkdir(exist_ok=True)
    out.to_csv(TSV, sep="\t", index=False)
    print(f"{TSV.name}: {len(out):,} fires ({time.time() - t0:.0f}s)")
    stats(out)


def stats(out):
    bins = [("0 (touch/inside)", out.road_km == 0), ("< 1 km", out.road_km < 1), ("< 5 km", out.road_km < 5),
            ("< 20 km", out.road_km < 20), ("< 50 km", out.road_km < 50), (">= 100 km (cap)", out.road_km >= CAP_KM)]
    for key, lab in (("per", "pyroCb-linked"), ("perAll", "other fires")):
        s = out[out.key == key]
        msk = {n: m[s.index] for n, m in bins}
        print(f"  {lab}: n={len(s):,}, median {s.road_km.median():.1f} km; " +
              ", ".join(f"{n} {msk[n].mean():.1%}" for n, _ in bins))


# ---- apply ---------------------------------------------------------------------------------------------------
def apply():
    t0 = time.time()
    tab = pd.read_csv(TSV, sep="\t")
    data = json.loads(DATA.read_text(encoding="utf-8"))
    for y in sorted(data, key=int):
        d = data[y]
        lut = {int(m): float(v) for m, v in tab.loc[tab.year == int(y), ["mergid", "road_km"]].itertuples(index=False)}
        n = {}
        for key in ("per", "perAll", "pcb", "ign", "ignAll"):
            if key not in d:
                continue
            fc = unpack(d[key])
            n[key] = 0
            for f in fc["features"]:
                p = f["properties"]
                p.pop("road_km", None)
                v = lut.get(p.get("mergid"))
                if v is not None:
                    p["road_km"] = v
                    n[key] += 1
            d[key] = pack(fc)
        print(f"  {y}: " + ", ".join(f"{k} {v:,}" for k, v in n.items()), flush=True)
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
    print(f"data.json {DATA.stat().st_size / 1e6:.1f} MB ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else ""
    if step == "extract":
        extract(sys.argv[2:] or None)
    elif step == "compute":
        compute()
    elif step == "apply":
        apply()
    elif step == "stats":
        stats(pd.read_csv(TSV, sep="\t"))
    else:
        print(__doc__)
