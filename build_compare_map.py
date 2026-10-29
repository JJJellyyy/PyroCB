"""Data for the atlas comparison page (compare.html next to the main map on Vercel): per year and layer one gzipped GeoJSON
in compare/<layer>_<year>.json.gz, fetched lazily and decompressed by the page (DecompressionStream).

Layers
  old      published Arctic-boreal atlas, final perimeters (2013-2023)
  joint    new multi-satellite atlas v1 (S-NPP + NOAA-20 + NOAA-21 tracked together), final perimeters (2012-2025)
  snpp / n20 / n21   new atlas tracked with one satellite only (years as far as available)
  modis    MCD64A1 v061 burned area (500 m) by agreement class with the joint atlas (pcbcube/mcd64.py):
           in = inside an atlas perimeter, add = within 3 km and burn date in the fire's active window (+-3 d),
           near = within 3 km but outside every window, omit = > 3 km from any atlas fire
Perimeters are simplified (SIMPLIFY_M in EPSG:3571) and rounded to 4 decimals; MODIS cells are polygonised per class
on the MODIS sinusoidal grid.

    python build_compare_map.py [--years 2012 ...] [--layers old joint snpp n20 n21 modis]
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from pcbcube import config as C  # noqa: E402

OUT = HERE / "compare"
AX = ROOT / "INPUT_DATA" / "atlas_extension"
REL = {"joint": AX / "release_multi", "snpp": AX / "release_snpp", "n20": AX / "release_n20", "n21": AX / "release_n21"}
SIMPLIFY_M = 100.0
# MODIS global sinusoidal grid (MCD64A1 500 m)
SIN_X0, SIN_Y0, SIN_CELL = -20015109.354, 10007554.677, 463.312716527778
SINU = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
CLS = {"in": 1, "add": 2, "near": 3, "omit": 4}


def _write(gdf, layer, year, props):
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]
    gdf = gdf.to_crs(4326)
    gdf["geometry"] = shapely.set_precision(gdf.geometry.values, 1e-4, mode="pointwise")   # display only
    gdf = gdf[~gdf.geometry.is_empty]
    gj = json.loads(gdf[props + ["geometry"]].to_json(drop_id=True))
    raw = json.dumps(gj, separators=(",", ":")).encode()
    import gzip
    (OUT / f"{layer}_{year}.json.gz").write_bytes(gzip.compress(raw, 9))       # decompressed in the browser
    old = OUT / f"{layer}_{year}.json"
    if old.exists():
        old.unlink()
    return len(gdf), (OUT / f"{layer}_{year}.json.gz").stat().st_size / 1e6


def perimeters(layer, year):
    if layer == "old":
        fp = glob.glob(str(C.INPUT / "VIIRS" / "**" / f"final_viirs{year}.gpkg"), recursive=True)
        if not fp or year > 2023:
            return None
        G = gpd.read_file(fp[0])
    else:
        fp = REL[layer] / "final_perims" / f"final_viirs{year}.gpkg"
        if not fp.exists():
            return None
        G = gpd.read_file(fp)
    G = G.to_crs(3571)
    G["geometry"] = shapely.simplify(G.geometry.values, SIMPLIFY_M, preserve_topology=True)
    G["start"] = [f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in zip(G.tst_year, G.tst_month, G.tst_day)]
    G["end"] = [f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in zip(G.ted_year, G.ted_month, G.ted_day)]
    G["farea"] = G.farea.round(2)
    props = ["mergid", "farea", "start", "end"]
    for c in ("n_snpp", "n_n20", "n_n21", "sat_first", "lcc_final"):
        if c in G:
            props.append(c)
    if layer == "joint":
        fa = ROOT / "OUTPUT_DATA" / "mcd64" / f"fire_agreement_{year}.parquet"
        if fa.exists():
            A = pd.read_parquet(fa, columns=["mergid", "frac_modis", "add_km2"])
            G = G.merge(A.groupby("mergid").agg(frac_modis=("frac_modis", "max"), add_km2=("add_km2", "sum")).round(3),
                        on="mergid", how="left")
            props += ["frac_modis", "add_km2"]
    return _write(G, layer, year, props)


def modis(year):
    import rasterio.features
    from affine import Affine
    fp = ROOT / "OUTPUT_DATA" / "mcd64" / f"cells_{year}.parquet"
    if not fp.exists():
        return None
    D = pd.read_parquet(fp, columns=["row", "col", "cls"])
    D["c"] = D.cls.map(CLS).astype(np.uint8)
    r0, r1, c0 = int(D.row.min()), int(D.row.max()), int(D.col.min())
    polys, cls = [], []
    W = 4096
    for cs in range(c0, int(D.col.max()) + 1, W):
        S = D[(D.col >= cs) & (D.col < cs + W)]
        if S.empty:
            continue
        a = np.zeros((r1 - r0 + 1, W), np.uint8)
        a[S.row.to_numpy() - r0, S.col.to_numpy() - cs] = S.c.to_numpy()
        tr = Affine(SIN_CELL, 0, SIN_X0 + cs * SIN_CELL, 0, -SIN_CELL, SIN_Y0 - r0 * SIN_CELL)
        for geom, v in rasterio.features.shapes(a, mask=a > 0, transform=tr, connectivity=4):
            polys.append(shapely.geometry.shape(geom))
            cls.append(int(v))
    inv = {v: k for k, v in CLS.items()}
    G = gpd.GeoDataFrame({"cls": [inv[c] for c in cls]}, geometry=polys, crs=SINU)
    G["geometry"] = shapely.simplify(G.geometry.values, SIN_CELL / 2, preserve_topology=True)
    return _write(G, "modis", year, ["cls"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=int, nargs="+", default=list(range(2012, 2026)))
    ap.add_argument("--layers", nargs="+", default=["old", "joint", "snpp", "n20", "n21", "modis"])
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    avail = {}
    for y in a.years:
        for L in a.layers:
            r = modis(y) if L == "modis" else perimeters(L, y)
            if r:
                avail.setdefault(L, []).append(y)
                print(f"{y} {L}: {r[0]:,} features, {r[1]:.1f} MB", flush=True)
    # merge with an existing index (layers built earlier stay listed)
    idx_fp = OUT / "index.json"
    idx = json.loads(idx_fp.read_text()) if idx_fp.exists() else {}
    for L, ys in avail.items():
        idx[L] = sorted(set(idx.get(L, [])) | set(ys))
    idx_fp.write_text(json.dumps(idx))
    print("index:", idx)


if __name__ == "__main__":
    main()
