"""Simplified land cover over each fire perimeter and a 5 km collar.

Source: the 30 m *Land Cover of Canada* rasters (EPSG:3979), 2015 collection for
2013-2017 and 2020 collection for 2018-2023 (`LC_YEAR_SPLIT`). Per fire we read the
raster windowed to (perimeter + `LC_BUFFER_M`), decimated to `LC_CELL_M` with a
majority (`mode`) resample so a 1 km cell takes the commonest 30 m class under it,
mask to the buffered perimeter, polygonise by class and dissolve. Everything is done
in the raster's native CRS - the raster is never reprojected - and only the final
dissolved layer is cast to WGS84 for the map.

Scope follows the map: only fires credited with a pyroCb.
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import shapely
from rasterio.enums import Resampling
from rasterio.windows import from_bounds

from . import config as C

# Land Cover of Canada classes present in the boreal zone: code -> (label, hex).
LC_CLASSES = {
    1:  ("Needleleaf forest",         "#003d00"),
    2:  ("Taiga needleleaf forest",   "#949c70"),
    5:  ("Broadleaf deciduous forest", "#148c3d"),
    6:  ("Mixed forest",              "#5c752b"),
    8:  ("Shrubland",                 "#b38a2d"),
    10: ("Grassland",                 "#e1cf8a"),
    11: ("Polar shrubland-lichen-moss", "#9c7554"),
    12: ("Polar grassland-lichen-moss", "#bad48f"),
    13: ("Polar barren-lichen-moss",  "#408a70"),
    14: ("Wetland",                   "#6ba38a"),
    15: ("Cropland",                  "#e6ae66"),
    16: ("Barren land",               "#a8abae"),
    17: ("Urban",                     "#dc2126"),
    18: ("Water",                     "#4c70a3"),
    19: ("Snow and ice",              "#fffaff"),
}


def _tif(year: int) -> str:
    return C.LC_TIF_2015 if year < C.LC_YEAR_SPLIT else C.LC_TIF_2020


def _poly_only(g):
    """Polygonal part of a geometry (make_valid can hand back a GeometryCollection)."""
    if g.geom_type in ("Polygon", "MultiPolygon"):
        return g
    polys = [p for p in getattr(g, "geoms", []) if p.geom_type in ("Polygon", "MultiPolygon")]
    return shapely.union_all(polys) if polys else shapely.Polygon()


def _polygonise(r, clip_geom, cell, simp):
    """Windowed, decimated (majority) read of `r` over `clip_geom`'s bounds, masked to
    `clip_geom`, polygonised by class and dissolved. Returns list of (cls, geom) in 3979."""
    minx, miny, maxx, maxy = clip_geom.bounds
    ow = max(1, round((maxx - minx) / cell))
    oh = max(1, round((maxy - miny) / cell))
    a = r.read(1, window=from_bounds(minx, miny, maxx, maxy, r.transform),
               out_shape=(oh, ow), resampling=Resampling.mode,
               boundless=True, fill_value=0).astype("uint8")
    tr = rasterio.transform.from_bounds(minx, miny, maxx, maxy, ow, oh)
    m = rasterio.features.geometry_mask([clip_geom], (oh, ow), tr, invert=True)
    a = np.where(m, a, 0).astype("uint8")
    out = []
    for geom, val in rasterio.features.shapes(a, mask=a > 0, transform=tr):
        v = int(val)
        if v not in LC_CLASSES:
            continue
        p = shapely.geometry.shape(geom).simplify(simp)
        if not p.is_valid:
            p = _poly_only(shapely.make_valid(p))
        if not p.is_empty:
            out.append((v, p))
    return out


def landcover_layer(year: int, geoms_3979) -> gpd.GeoDataFrame:
    """Dissolved [cls, res, geometry] land-cover polygons (WGS84) - the COARSE backdrop
    only, at `LC_BG_CELL_M` (10 km) over a wide collar around each fire.

    Everything closer in is a paletted PNG overlay now (`pyrocb/lcraster.py`): finer,
    smaller, GPU-composited, and it skips the slow `features.shapes` + dissolve that
    used to dominate the build. The collar stays vector because at 10 km it is a
    handful of cheap polygons and it has to tile a very large area.
    """
    geoms_3979 = [_poly_only(shapely.make_valid(g)) for g in geoms_3979
                  if g is not None and not g.is_empty]
    geoms_3979 = [g for g in geoms_3979 if not g.is_empty]
    if not geoms_3979:
        return gpd.GeoDataFrame({"cls": [], "res": []},
                                geometry=gpd.GeoSeries([], crs=C.WGS84))

    # The PNG overlays cover perimeter + LC_BUFFER_M, so punch that out of the collar.
    png_union = shapely.union_all([g.buffer(C.LC_BUFFER_M) for g in geoms_3979])
    coarse = []
    with rasterio.open(_tif(year)) as r:
        rbounds = shapely.box(*r.bounds)
        # A wide collar around EACH fire - never a hull over all of them, since fires
        # span the dateline and EPSG:3979 is Canada-only, so a global hull would be a
        # garbage megapolygon and a huge read.
        for g in geoms_3979:
            collar = g.buffer(C.LC_BG_BUFFER_M).difference(png_union).intersection(rbounds)
            if collar.is_empty:
                continue
            for v, poly in _polygonise(r, collar, C.LC_BG_CELL_M, C.LC_BG_SIMPLIFY_M):
                try:
                    p = poly.intersection(collar)
                except Exception:
                    p = poly
                if not p.is_empty:
                    coarse.append((v, p, C.LC_BG_CELL_M))

    # Group by (res, cls) and union with a grid snap - more robust than
    # GeoDataFrame.dissolve, which threw "unable to assign free hole to a shell".
    groups: dict = {}
    for v, g, res in coarse:
        if g is None or g.is_empty:
            continue
        g = _poly_only(shapely.make_valid(g))
        if not g.is_empty:
            groups.setdefault((int(res), int(v)), []).append(g)
    if not groups:
        return gpd.GeoDataFrame({"cls": [], "res": []},
                                geometry=gpd.GeoSeries([], crs=C.WGS84))
    recs = []
    for (res, cls), gs in groups.items():
        try:
            u = shapely.union_all(gs, grid_size=1.0)
        except Exception:
            u = shapely.union_all([shapely.set_precision(x, 5.0) for x in gs])
        u = _poly_only(shapely.make_valid(u))
        if not u.is_empty:
            recs.append((res, cls, u))
    diss = gpd.GeoDataFrame({"res": [r[0] for r in recs], "cls": [r[1] for r in recs]},
                            geometry=[r[2] for r in recs], crs=3979)
    diss["geometry"] = shapely.set_precision(diss.geometry.values, 30.0)
    diss = diss[~diss.geometry.is_empty & diss.geometry.notna()]
    return diss.to_crs(C.WGS84)


def active_fractions(year: int, ids, geoms_3979) -> pd.DataFrame:
    """Per pyroCb event, the land-cover composition of its active-burn polygon,
    sampled at the raster's native 30 m resolution.

    Returns one row per id: top-3 class codes and their fractions of the classified
    burn area, plus `lc_n_px` (30 m cells sampled).
    """
    recs = []
    with rasterio.open(_tif(year)) as r:
        for pid, g in zip(ids, geoms_3979):
            rec = {"pcb_id": pid, "lc_n_px": 0}
            if g is not None and not g.is_empty:
                minx, miny, maxx, maxy = g.bounds
                win = from_bounds(minx, miny, maxx, maxy, r.transform)
                a = r.read(1, window=win, boundless=True, fill_value=0).astype("uint8")
                tr = r.window_transform(win)
                m = rasterio.features.geometry_mask([g], a.shape, tr, invert=True)
                vals = a[m & (a > 0)]
                cnt = np.bincount(vals, minlength=max(LC_CLASSES) + 1).astype(float)
                keep = np.array([c for c in LC_CLASSES if cnt[c] > 0])
                tot = float(cnt[keep].sum()) if len(keep) else 0.0
                if tot > 0:
                    order = keep[np.argsort(cnt[keep])[::-1]][:3]
                    rec["lc_n_px"] = int(tot)
                    for i, c in enumerate(order, 1):
                        rec[f"lc_top{i}_cls"] = int(c)
                        rec[f"lc_top{i}_frac"] = round(float(cnt[c] / tot), 3)
            recs.append(rec)
    return pd.DataFrame(recs)
