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


def landcover_layer(year: int, geoms_3979) -> gpd.GeoDataFrame:
    """Dissolved [cls, geometry] land-cover polygons (WGS84) for the given fire
    geometries (already in EPSG:3979), each grown by `LC_BUFFER_M`."""
    cell, buf, simp = C.LC_CELL_M, C.LC_BUFFER_M, C.LC_SIMPLIFY_M
    parts = []
    with rasterio.open(_tif(year)) as r:
        for g in geoms_3979:
            bg = g.buffer(buf)
            minx, miny, maxx, maxy = bg.bounds
            ow = max(1, round((maxx - minx) / cell))
            oh = max(1, round((maxy - miny) / cell))
            a = r.read(1, window=from_bounds(minx, miny, maxx, maxy, r.transform),
                       out_shape=(oh, ow), resampling=Resampling.mode,
                       boundless=True, fill_value=0).astype("uint8")
            tr = rasterio.transform.from_bounds(minx, miny, maxx, maxy, ow, oh)
            m = rasterio.features.geometry_mask([bg], (oh, ow), tr, invert=True)
            a = np.where(m, a, 0).astype("uint8")
            for geom, val in rasterio.features.shapes(a, mask=a > 0, transform=tr):
                v = int(val)
                if v in LC_CLASSES:
                    parts.append((v, shapely.geometry.shape(geom)))

    if not parts:
        return gpd.GeoDataFrame({"cls": []}, geometry=gpd.GeoSeries([], crs=C.WGS84))
    gdf = gpd.GeoDataFrame({"cls": [p[0] for p in parts]},
                           geometry=[p[1] for p in parts], crs=3979)
    diss = gdf.dissolve("cls").reset_index()
    diss["geometry"] = diss.geometry.simplify(simp)
    diss = diss[~diss.geometry.is_empty & diss.geometry.notna()]
    diss["geometry"] = shapely.set_precision(diss.geometry.values, 30.0)
    return diss.to_crs(C.WGS84)
