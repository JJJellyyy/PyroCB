"""Permafrost for the HTML map (map only, not the dataset).

Source: Obu et al. (2019), Northern Hemisphere permafrost map based on TTOP modelling for 2000-2016 at 1 km2 scale,
Earth-Science Reviews 193 - the PERZONES product (PANGAEA 888600, CC-BY 3.0): permafrost zones from the modelled permafrost
probability - Continuous (> 90 %), Discontinuous (50-90 %), Sporadic (10-50 %), Isolated (< 10 %).

Writes
  permafrost.json            the four zones as simplified polygons (lon/lat), loaded lazily by the map's "Permafrost" layer
  map_src/data.json          every fire perimeter (DATA[y].per and .perAll) gets the share of its area on each zone, in
                             whole percent: pf_c (continuous), pf_d (discontinuous), pf_s (sporadic), pf_i (isolated);
                             a zone with 0 % is omitted. Areas are computed in the equal-area EPSG:3571 from the
                             un-simplified zone polygons.

    python build_permafrost.py [year ...] [--propagate-only]     then   python map_tool.py join
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
SHP = HERE.parent / "INPUT_DATA" / "PERMAFROST"
DATA = HERE / "map_src" / "data.json"
OUT = HERE / "permafrost.json"
CODE = {4: ("c", "Continuous"), 3: ("d", "Discontinuous"), 2: ("s", "Sporadic"), 1: ("i", "Isolated")}
MIN_LAT = 40.0


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def pack(o):
    return base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def load_zones() -> gpd.GeoDataFrame:
    f = glob.glob(str(SHP / "**" / "*PERZONES*.shp"), recursive=True)[0]
    g = gpd.read_file(f).set_crs(3995, allow_override=True)          # Arctic polar stereographic (the .prj axis tags confuse GDAL)
    g["geometry"] = shapely.make_valid(g.geometry.to_numpy())
    return g


def fix_antimeridian(geom):
    """Polygons that cross the 180th meridian come out of the reprojection with vertices jumping from +180 to -180, which a
    web map draws as a bar across the whole map. Shift such polygons into 0..360, then cut them at 180 and move the east
    part back to -180..."""
    parts = []
    for poly in shapely.get_parts(geom):
        minx, _, maxx, _ = shapely.bounds(poly)
        if maxx - minx > 180:
            poly = shapely.transform(poly, lambda c: np.column_stack([np.where(c[:, 0] < 0, c[:, 0] + 360.0, c[:, 0]), c[:, 1]]))
            west = shapely.clip_by_rect(poly, 0, -90, 180, 90)
            east = shapely.transform(shapely.clip_by_rect(poly, 180, -90, 360, 90), lambda c: np.column_stack([c[:, 0] - 360.0, c[:, 1]]))
            parts += [west, east]
        else:
            parts.append(poly)
    flat = [q for p in parts for q in shapely.get_parts(p) if shapely.get_type_id(q) == 3 and not shapely.is_empty(q)]   # polygons only
    return shapely.multipolygons(flat)


def write_zones(g: gpd.GeoDataFrame):
    simp = g.copy()
    simp["geometry"] = shapely.simplify(simp.geometry.to_numpy(), 500.0)        # 1 km cells: 0.5 km keeps the shapes
    feats = []
    for code, (short, name) in CODE.items():
        part = simp[simp["GRIDCODE"] == code]
        u = shapely.make_valid(shapely.union_all(shapely.make_valid(part.geometry.to_numpy()), grid_size=1.0))
        # cut along the 180th meridian BEFORE reprojecting: in polar stereographic (central meridian 0) it is the +y axis.
        # A 100 m wide slit means no piece can straddle +-180 and wrap around the map after reprojection.
        u = shapely.make_valid(shapely.difference(u, shapely.box(-50.0, 0.0, 50.0, 2.0e7)))
        pieces = [q for q in shapely.get_parts(u) if shapely.get_type_id(q) == 3]
        gg = gpd.GeoDataFrame({"geometry": pieces}, crs=3995).to_crs(4326)
        ok = gg.geometry.bounds.eval("maxx - minx") < 180                     # anything still wrapping is dropped (none expected)
        geom = shapely.multipolygons([q for g_ in gg.geometry[ok] for q in shapely.get_parts(shapely.make_valid(g_)) if shapely.get_type_id(q) == 3])
        geom = shapely.clip_by_rect(geom, -180, MIN_LAT, 180, 90)
        print(f"  {name}: {len(pieces):,} pieces, {int((~ok).sum())} still wrapping (dropped)", flush=True)
        feats.append({"type": "Feature", "properties": {"zone": short, "name": name},
                      "geometry": json.loads(shapely.to_geojson(shapely.transform(geom, lambda c: np.round(c, 3))))})
    OUT.write_text(json.dumps({"type": "FeatureCollection", "features": feats,
                               "meta": {"source": "Obu et al. 2019, Earth-Science Reviews 193 (PANGAEA 888600, CC-BY 3.0), 2000-2016 TTOP model, 1 km",
                                        "zones": {v[0]: v[1] for v in CODE.values()}}}, separators=(",", ":")), encoding="utf-8")
    print(f"permafrost.json: {OUT.stat().st_size / 1e6:.1f} MB")


def overlaps(fc: dict, zones3571, tree, codes) -> int:
    """Add pf_* percentages to the features of a perimeter FeatureCollection (modifies fc); returns features touched."""
    feats = fc["features"]
    if not feats:
        return 0
    g = gpd.GeoDataFrame.from_features(feats, crs=4326).to_crs(3571)
    fg = shapely.make_valid(g.geometry.to_numpy())
    area = shapely.area(fg)
    fi, zi = tree.query(fg, predicate="intersects")
    n = 0
    if len(fi):
        inter = shapely.area(shapely.intersection(fg[fi], zones3571[zi]))
        for code, (short, _) in CODE.items():
            m = codes[zi] == code
            tot = np.bincount(fi[m], weights=inter[m], minlength=len(fg))
            pct = np.where(area > 0, np.clip(np.rint(100.0 * tot / np.where(area > 0, area, 1)), 0, 100), 0).astype(int)
            for k in np.flatnonzero(pct):
                feats[k]["properties"][f"pf_{short}"] = int(pct[k])
        n = int(np.unique(fi).size)
    for f in feats:                                      # idempotent: clear old values of fires that no longer overlap
        pass
    return n


def propagate(d: dict) -> int:
    """Copy each fire's pf_* shares (from DATA[y].per / .perAll, keyed by mergid within the year) onto the features of
    DATA[y].pcb / .ign / .ignAll, so popups and the info panel need no lookup."""
    pf = {}
    for key in ("per", "perAll"):
        if key in d:
            for f in unpack(d[key])["features"]:
                p = f["properties"]
                pf[p["mergid"]] = {k: v for k, v in p.items() if k.startswith("pf_")}
    n = 0
    for key in ("pcb", "ign", "ignAll"):
        if key not in d:
            continue
        fc = unpack(d[key])
        for f in fc["features"]:
            p = f["properties"]
            for k in ("pf_c", "pf_d", "pf_s", "pf_i"):
                p.pop(k, None)
            if p.get("mergid") in pf:
                p.update(pf[p["mergid"]])
                n += bool(pf[p["mergid"]])
        d[key] = pack(fc)
    return n


def main(years=None, propagate_only=False):
    t0 = time.time()
    if propagate_only:
        data = json.loads(DATA.read_text(encoding="utf-8"))
        for y in years or sorted(data, key=int):
            print(f"  {y}: {propagate(data[str(y)]):,} pyroCb / ignition features carry permafrost shares", flush=True)
        DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
        return
    g = load_zones()
    if "--zones-only" in sys.argv:
        write_zones(g)
        return
    if not OUT.exists() or not years:
        write_zones(g)
    z = g.to_crs(3571)
    zones3571, codes = z.geometry.to_numpy(), z["GRIDCODE"].to_numpy()
    tree = STRtree(zones3571)
    print(f"zones ready ({len(zones3571):,} polygons, {time.time() - t0:.0f}s)", flush=True)
    data = json.loads(DATA.read_text(encoding="utf-8"))
    for y in years or sorted(data, key=int):
        d = data[str(y)]
        for key in ("per", "perAll"):
            if key not in d:
                continue
            fc = unpack(d[key])
            for f in fc["features"]:                     # start clean so reruns do not keep stale values
                for s in ("c", "d", "s", "i"):
                    f["properties"].pop(f"pf_{s}", None)
            n = overlaps(fc, zones3571, tree, codes)
            d[key] = pack(fc)
            share = np.mean([any(f["properties"].get(f"pf_{s}", 0) for s in ("c", "d", "s")) for f in fc["features"]]) if fc["features"] else 0
            print(f"  {y} {key}: {len(fc['features']):>6} fires, {n:>6} touch a zone, {share:.0%} on continuous/discontinuous/sporadic", flush=True)
        propagate(d)
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
    print(f"data.json {DATA.stat().st_size / 1e6:.1f} MB ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]   # --zones-only: rewrite permafrost.json only
    main([int(a) for a in args] or None, propagate_only="--propagate-only" in sys.argv)
