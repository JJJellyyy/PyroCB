"""Paths, constants and the caching helper for the pyroCb x FireAtlas pipeline.

Everything tunable lives here so the notebook stays narrative.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- scope
# YEARS = None  ->  every year that has BOTH fire data on disk AND pyroCb records.
#                   That is 2013-2023; 2012 drops out because the pyroCb inventory
#                   does not cover it.
YEARS: list[int] | None = None

# Years for the heavy snapshot/NFP analysis (matching, bands, per-event variables).
# None -> same as the resolved YEARS.
SNAPSHOT_YEARS: list[int] | None = None

# ---------------------------------------------------------------- paths
PERIMS_TPL = str(ROOT / "final_perims" / "final_perims" / "final_viirs{y}.gpkg")
IGN_TPL    = str(ROOT / "ignitions" / "ignitions" / "ign{y}.gpkg")
SNAP_DIR   = str(ROOT / "{y}" / "{y}" / "Snapshot")
NFP_DIR    = str(ROOT / "{y}" / "{y}" / "NFP")
XLSX_FP    = ROOT / "PCBAtlas.xlsx"

OUT_HTML   = ROOT / "fire_pcb_map.html"
OUT_TSV    = ROOT / "pyrocb_fire_attributes.tsv"
OUT_XLSX   = ROOT / "pyrocb_fire_attributes.xlsx"
OUT_DICT   = ROOT / "pyrocb_columns.tsv"
OUT_ACTIVE = ROOT / "pyrocb_active_areas.gpkg"
CACHE_DIR  = ROOT / "cache"

# ---------------------------------------------------------------- CRS
# EPSG:3571 is "WGS 84 / North Pole LAEA Bering Sea": equal-area, centred on lon 180,
# so the dateline - where Chukotka fires sit - is the LEAST distorted part of the
# projection.  All metric work (areas, distances, geometry booleans) happens here.
# (It is NOT EASE-Grid 2.0 North; that is EPSG:6931.)
METRIC_CRS = 3571
WGS84 = 4326

# ---------------------------------------------------------------- geometry
PRECISION_M      = 1.0    # snap grid before every boolean; 1/375 of a VIIRS pixel
SLIVER_KM2       = 0.02   # drop polygon parts below this
MIN_BAND_KM2     = 0.5    # ~3.5 VIIRS pixels; below this a band keeps its row, loses geometry
BAND_SIMPLIFY_M  = 50.0   # band simplification for the map payload (in METRIC_CRS)
PIXEL_AREA_KM2   = 0.141  # exact area FEDS assigns a single VIIRS detection

# Background layers (the ~200k perimeters and ignitions) - much coarser, see notebook.
SIMPLIFY_TOL    = 0.003   # degrees (~330 m, ~1 VIIRS pixel)
COORD_PRECISION = 1e-5    # degrees (~1 m)

# ---------------------------------------------------------------- matching / bands
MATCH_CAP_KM        = 20.0
MATCH_CAP_SENSITIVITY = (10.0, 20.0, 30.0)
MAX_BACKWALK_STEPS  = 30      # 15 days
MATCH_FORWARD_STEPS = 2       # look this many overpasses AHEAD for a fire that blew up
                              # before it was first observed (the nightly ~10 h gap)
NFP_WINDOWS_H       = (6, 12, 24)
# Active-burn polygon: NFP pixel footprints within +/- W hours of onset. Symmetric,
# because the ~10 h nightly blind window means the run that fed the plume is usually
# first seen by the overpass AFTER onset (FINDINGS F3).
ACTIVE_WINDOWS_H    = (6, 12, 24)
ACTIVE_HEADLINE_H   = 12

# Map scope. True -> the map carries only fires credited with a pyroCb (and their
# ignitions), which is the analysis population; False -> every tracked fire, which is
# ~170k polygons and a 27 MB file.
MAP_ONLY_PCB_FIRES  = True
STALE_OBS_H         = 14.0    # no overpass within this of onset -> qc_stale_obs

# Local-solar boundaries of the AM/PM overpass groups, measured from NFP datetime_utc.
LOCAL_AM_END_H = 7.0
LOCAL_PM_END_H = 14.4

# ---------------------------------------------------------------- colour
# Cool -> warm ramp, one colour per year. Also reused (interpolated) for the
# day-of-year gradient.
YEAR_COLORS = {
    2012: "#3b4cc0", 2013: "#5a78d0", 2014: "#7ba3dc", 2015: "#2a9d8f",
    2016: "#3f9b45", 2017: "#7bab2e", 2018: "#c9a227", 2019: "#e08a1e",
    2020: "#e2660f", 2021: "#d1401f", 2022: "#b02020", 2023: "#7d1128",
}
DOY_RAMP = [YEAR_COLORS[y] for y in sorted(YEAR_COLORS)]
# DOY_RANGE is fitted to the whole dataset at build time (see variables.doy_range),
# so it extends automatically if 2026+ are added.  This is only the fallback.
DOY_RANGE_FALLBACK = (90, 300)


# ---------------------------------------------------------------- year discovery
def discover_fire_years() -> list[int]:
    """Years with a final-perimeter file, an ignition file and a Snapshot directory."""
    out = []
    for p in sorted(Path(PERIMS_TPL).parent.glob("final_viirs*.gpkg")):
        m = re.search(r"final_viirs(\d{4})\.gpkg$", p.name)
        if not m:
            continue
        y = int(m.group(1))
        if Path(IGN_TPL.format(y=y)).exists() and Path(SNAP_DIR.format(y=y)).is_dir():
            out.append(y)
    return out


def resolve_years(pcb_years: set[int] | None = None) -> list[int]:
    """YEARS if set, else every fire year that also has pyroCb records."""
    if YEARS is not None:
        return list(YEARS)
    fire = discover_fire_years()
    return [y for y in fire if pcb_years is None or y in pcb_years]


def resolve_snapshot_years(years: list[int]) -> list[int]:
    return list(years) if SNAPSHOT_YEARS is None else [y for y in SNAPSHOT_YEARS if y in years]


# ---------------------------------------------------------------- caching
CODE_VERSION = "1"   # bump to invalidate every cached stage


def params_hash(params: dict) -> str:
    blob = json.dumps(params, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:12]


def _stamp(paths) -> list:
    out = []
    for p in paths:
        p = Path(p)
        try:
            st = p.stat()
            out.append([str(p), int(st.st_mtime), st.st_size])
        except OSError:
            out.append([str(p), None, None])
    return out


def cached(stage: str, year, params: dict, inputs, builder, *, geo: bool = False):
    """Build-or-load one pipeline stage.

    Keyed on CODE_VERSION + a hash of `params` + the (mtime, size) of `inputs`, so
    changing an unrelated constant does not invalidate the stage.
    `builder()` must return a DataFrame (geo=False) or GeoDataFrame (geo=True).
    """
    import pandas as pd

    CACHE_DIR.mkdir(exist_ok=True)
    name = f"{stage}_{year}" if year is not None else stage
    data_fp = CACHE_DIR / f"{name}.parquet"
    meta_fp = CACHE_DIR / f"{name}_meta.json"

    want = {"code": CODE_VERSION, "params": params_hash(params), "inputs": _stamp(inputs)}
    if data_fp.exists() and meta_fp.exists():
        try:
            if json.loads(meta_fp.read_text()) == want:
                if geo:
                    import geopandas as gpd
                    return gpd.read_parquet(data_fp)
                return pd.read_parquet(data_fp)
        except Exception:
            pass  # unreadable cache -> rebuild

    obj = builder()
    if geo:
        obj.to_parquet(data_fp)
    else:
        obj.to_parquet(data_fp, index=False)
    meta_fp.write_text(json.dumps(want))
    return obj


def clear_cache(stage: str | None = None) -> int:
    """Delete cached stages (all, or one stage prefix). Returns files removed."""
    if not CACHE_DIR.exists():
        return 0
    pat = f"{stage}_*" if stage else "*"
    n = 0
    for p in list(CACHE_DIR.glob(pat)):
        if p.is_file():
            p.unlink()
            n += 1
    return n
