"""Difference layers for compare.html: which AREA one atlas version adds or loses relative to another.
Computed in EPSG:3571 (equal-area) from the union of all final perimeters per atlas and year; pieces < MIN_KM2 dropped,
then simplified (SIMPLIFY_M) and written as compare/<layer>_<year>.json.gz with properties kind, area_km2 (and sats).

Layers (key -> kinds)
  new_vs_old     gained = in the joint new atlas, not in the old (published) atlas; lost = in old, not in joint
  modis_vs_old   modis_not_in_old = all MCD64A1 burned cells (in/add/near/omit) outside the old perimeters;
                 modis_add_not_in_old = only the 'add' cells outside the old perimeters
  multisat_gain  joint_only = in joint but in none of the available single-satellite versions;
                 singles_only = in the union of the single-satellite versions but not in joint (sats = versions used)

    python build_compare_diff.py [--years ...] [--layers new_vs_old modis_vs_old multisat_gain]
"""
from __future__ import annotations

import argparse
import glob
import json
import time

import psutil

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from build_compare_map import (OUT, REL, ROOT, SIMPLIFY_M, SIN_CELL, SIN_X0, SIN_Y0, SINU, C, _write)

MIN_KM2 = 0.05
SINGLES = ("snpp", "n20", "n21")
CRS = 3571


def union_perims(layer, year):
    """Union (EPSG:3571) of all final perimeters of one atlas version and year, or None."""
    if layer == "old":
        fp = glob.glob(str(C.INPUT / "VIIRS" / "**" / f"final_viirs{year}.gpkg"), recursive=True)
        if not fp or year > 2023:
            return None
        fp = fp[0]
    else:
        fp = REL[layer] / "final_perims" / f"final_viirs{year}.gpkg"
        if not fp.exists():
            return None
    G = gpd.read_file(fp, columns=[]).to_crs(CRS)
    g = shapely.make_valid(G.geometry.values)
    return _union(g)


def modis_union(year, classes):
    """Union (EPSG:3571) of the MCD64A1 cells of the given classes, polygonised on the sinusoidal grid."""
    import rasterio.features
    from affine import Affine
    fp = ROOT / "OUTPUT_DATA" / "mcd64" / f"cells_{year}.parquet"
    if not fp.exists():
        return None
    D = pd.read_parquet(fp, columns=["row", "col", "cls"])
    D = D[D.cls.isin(classes)]
    if D.empty:
        return shapely.Polygon()
    r0, r1, c0 = int(D.row.min()), int(D.row.max()), int(D.col.min())
    polys = []
    W = 4096
    for cs in range(c0, int(D.col.max()) + 1, W):
        S = D[(D.col >= cs) & (D.col < cs + W)]
        if S.empty:
            continue
        a = np.zeros((r1 - r0 + 1, W), np.uint8)
        a[S.row.to_numpy() - r0, S.col.to_numpy() - cs] = 1
        tr = Affine(SIN_CELL, 0, SIN_X0 + cs * SIN_CELL, 0, -SIN_CELL, SIN_Y0 - r0 * SIN_CELL)
        for geom, v in rasterio.features.shapes(a, mask=a > 0, transform=tr, connectivity=4):
            polys.append(shapely.geometry.shape(geom))
    G = gpd.GeoSeries(polys, crs=SINU).to_crs(CRS)
    return _union(shapely.make_valid(G.values))


def pieces(geom, kind, extra=None):
    """Explode to polygons, drop slivers, return a GeoDataFrame (area_km2 from the unsimplified piece)."""
    parts = shapely.get_parts(shapely.get_parts(geom))
    parts = parts[shapely.get_type_id(parts) == 3]            # polygons only (drop lines/points from overlays)
    area = shapely.area(parts) / 1e6
    keep = area >= MIN_KM2
    parts, area = parts[keep], area[keep]
    G = gpd.GeoDataFrame({"kind": kind, "area_km2": np.round(area, 2)}, geometry=parts, crs=CRS)
    for k, v in (extra or {}).items():
        G[k] = v
    return G


def _union(geoms):
    # no grid_size: snapping to a grid makes GEOS overlay ~20x slower on these pixel-edged perimeters
    try:
        return shapely.union_all(geoms)
    except shapely.errors.GEOSException:
        return shapely.union_all(geoms, grid_size=1.0)


def diff(a, b):
    try:
        return shapely.difference(a, b)
    except shapely.errors.GEOSException:
        return shapely.difference(a, b, grid_size=1.0)


def emit(layer, year, frames, props):
    G = pd.concat([f for f in frames if len(f)], ignore_index=True) if any(len(f) for f in frames) else None
    if G is None:
        G = gpd.GeoDataFrame({p: [] for p in props}, geometry=[], crs=CRS)
    G = gpd.GeoDataFrame(G, geometry="geometry", crs=CRS)
    G["geometry"] = shapely.simplify(G.geometry.values, SIMPLIFY_M, preserve_topology=True)
    n, mb = _write(G, layer, year, props)
    return n, mb


def main():
    try:
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)    # other heavy jobs run on this laptop
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=int, nargs="+", default=list(range(2012, 2026)))
    ap.add_argument("--layers", nargs="+", default=["new_vs_old", "modis_vs_old", "multisat_gain"])
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    avail, summary = {}, []
    for y in a.years:
        t0 = time.time()
        old = union_perims("old", y)
        joint = union_perims("joint", y) if {"new_vs_old", "multisat_gain"} & set(a.layers) else None
        if "new_vs_old" in a.layers and old is not None and joint is not None:
            fr = [pieces(diff(joint, old), "gained"), pieces(diff(old, joint), "lost")]
            n, mb = emit("new_vs_old", y, fr, ["kind", "area_km2"])
            avail.setdefault("new_vs_old", []).append(y)
            summary += [(y, "gained", fr[0].area_km2.sum()), (y, "lost", fr[1].area_km2.sum())]
            print(f"{y} new_vs_old: {n:,} pieces, {mb:.2f} MB, gained {fr[0].area_km2.sum():,.0f} km2, "
                  f"lost {fr[1].area_km2.sum():,.0f} km2", flush=True)
        if "modis_vs_old" in a.layers and old is not None:
            mall = modis_union(y, ["in", "add", "near", "omit"])
            if mall is not None:
                madd = modis_union(y, ["add"])
                fr = [pieces(diff(mall, old), "modis_not_in_old"), pieces(diff(madd, old), "modis_add_not_in_old")]
                n, mb = emit("modis_vs_old", y, fr, ["kind", "area_km2"])
                avail.setdefault("modis_vs_old", []).append(y)
                summary += [(y, "modis_not_in_old", fr[0].area_km2.sum()), (y, "modis_add_not_in_old", fr[1].area_km2.sum())]
                print(f"{y} modis_vs_old: {n:,} pieces, {mb:.2f} MB, not-in-old {fr[0].area_km2.sum():,.0f} km2, "
                      f"add-not-in-old {fr[1].area_km2.sum():,.0f} km2", flush=True)
        if "multisat_gain" in a.layers and joint is not None:
            singles = {s: union_perims(s, y) for s in SINGLES}
            singles = {s: g for s, g in singles.items() if g is not None}
            if singles:
                su = _union(list(singles.values()))
                sats = ",".join(singles)
                fr = [pieces(diff(joint, su), "joint_only", {"sats": sats}),
                      pieces(diff(su, joint), "singles_only", {"sats": sats})]
                n, mb = emit("multisat_gain", y, fr, ["kind", "area_km2", "sats"])
                avail.setdefault("multisat_gain", []).append(y)
                summary += [(y, "joint_only", fr[0].area_km2.sum()), (y, "singles_only", fr[1].area_km2.sum())]
                print(f"{y} multisat_gain ({sats}): {n:,} pieces, {mb:.2f} MB, joint-only {fr[0].area_km2.sum():,.0f} km2, "
                      f"singles-only {fr[1].area_km2.sum():,.0f} km2", flush=True)
        print(f"  {y} done in {time.time() - t0:.0f} s", flush=True)
    idx_fp = OUT / "index.json"
    idx = json.loads(idx_fp.read_text()) if idx_fp.exists() else {}
    for L, ys in avail.items():
        idx[L] = sorted(set(idx.get(L, [])) | set(ys))
    idx_fp.write_text(json.dumps(idx))
    if summary:
        S = pd.DataFrame(summary, columns=["year", "kind", "km2"]).pivot(index="year", columns="kind", values="km2")
        print(S.round(0).to_string())
        S.round(1).to_csv(OUT.parent / f"compare_diff_summary_{a.years[0]}_{a.years[-1]}.csv")
    print("index:", idx)


if __name__ == "__main__":
    main()
