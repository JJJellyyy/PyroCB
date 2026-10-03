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
# The fire atlas, inventory and land-cover rasters live ONCE, in the parent project's INPUT_DATA
# (PCB_Analysis_1/INPUT_DATA); this map project no longer keeps its own copies.
DATA = ROOT.parent / "INPUT_DATA"

# ---------------------------------------------------------------- scope
# YEARS = None  ->  every year that has BOTH fire data on disk AND pyroCb records.
#                   That is 2013-2023; 2012 drops out because the pyroCb inventory
#                   does not cover it.
YEARS: list[int] | None = list(range(2013, 2024))   # pinned 2026-10-03: project is 2013-2023 (published atlas); 2024-25 kept on disk, 13-year map in map_backup_2013_2025/

# Years for the heavy snapshot/NFP analysis (matching, bands, per-event variables).
# None -> same as the resolved YEARS.
SNAPSHOT_YEARS: list[int] | None = None

# ---------------------------------------------------------------- paths
PERIMS_TPL = str(DATA / "VIIRS" / "final_perims" / "final_perims" / "final_viirs{y}.gpkg")
IGN_TPL    = str(DATA / "VIIRS" / "ignitions" / "ignitions" / "ign{y}.gpkg")
SNAP_DIR   = str(DATA / "VIIRS" / "{y}" / "{y}" / "Snapshot")
NFP_DIR    = str(DATA / "VIIRS" / "{y}" / "{y}" / "NFP")
# The pyroCb inventory is the parent project's extended one (2013-2025); the map reads it through
# pcbcube.pcb so ids, the boreal filter and the Confirmed-only rule are the same as the dataset's.
import sys as _sys
_sys.path.insert(0, str(ROOT.parent))
from pcbcube import config as _PC          # noqa: E402
XLSX_FP    = Path(_PC.PCB_XLSX)
XLSX_LEGACY = ROOT / "PCBAtlas.xlsx"       # the original 2013-2023 inventory (cache stamp for those years)


def inventory_stamp(year: int) -> Path:
    """File whose mtime/size invalidates a year's caches. The extension years follow the extended
    inventory; 2013-2023 stay stamped with the original so their caches survive its extension."""
    return XLSX_LEGACY if (year < _PC.EXT_FIRST_YEAR and XLSX_LEGACY.exists()) else XLSX_FP

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
BAND_SIMPLIFY_M  = 50.0   # band simplification for the map payload (in METRIC_CRS).
                          # SIMPLIFY ONLY - never a buffer open/close. Douglas-Peucker
                          # deviation is BOUNDED by the tolerance (measured: max 92 m at
                          # 60 m); a +/-60 m close displaced band boundaries by up to
                          # 4.3 km, because bands are thin rings and an open/close
                          # annihilates anything narrower than 2x the radius (F14).
PIXEL_AREA_KM2   = 0.141  # exact area FEDS assigns a single VIIRS detection

# Map payload geometry fidelity. FEDS perimeters are alpha shapes fitted over 375 m
# VIIRS pixels, so the SOURCE has ~375 m granularity; simplifying an order of magnitude
# below that is visually lossless against the real data, and anything finer is just
# storing alpha-shape noise. 0.003 deg (~330 m) used to be the tolerance and was
# plainly visible when zoomed into a fire.
SIMPLIFY_TOL    = 0.0003  # degrees (~33 m, ~1/11 of a VIIRS pixel)
COORD_PRECISION = 1e-5    # degrees (~1.1 m) - ~340x finer than a VIIRS pixel and
                          # well under the simplify tolerance, so the snap never
                          # dominates the error budget. 1e-6 only bought decimetres
                          # nothing can see, at a real cost in bytes.

# ---------------------------------------------------------------- matching / bands
MATCH_CAP_KM        = 20.0
MATCH_CAP_SENSITIVITY = (10.0, 20.0, 30.0)
MAX_BACKWALK_STEPS  = 30      # 15 days
MATCH_FORWARD_STEPS = 2       # look this many overpasses AHEAD for a fire that blew up
                              # before it was first observed (the nightly ~10 h gap)

# "Matched fire v2" (pyrocb/nearby.py): a second, deliberately simple pyroCb-position
# estimate that ignores fire identity (mergid) entirely - just "what was actively
# burning within this radius of the report, right around this time". See that module
# for the full method. Kept separate from MATCH_CAP_KM: that cap governs the
# containment/nearest-fire mergid match, this radius governs a totally different,
# identity-free localisation.
# Computed at every radius in this tuple. Was (3.0, 5.0, 10.0) with 3/5 km dropped for
# being dominated by the atlas's own +/-5.6 km position error rather than anything the
# method could fix (D-55). 15/20 km added back on request (D-66) for side-by-side
# comparison against 10 km, now that NEARBY_METHODS also varies which snapshot(s)
# each radius is measured against - a wider radius may behave differently under
# "before2_after2" than it did alone.
NEARBY_RADII_KM      = (10.0, 15.0, 20.0)
# How many overpasses to try, nearest-in-time first, before giving up on an event.
# A step only has to have detected fire SOMEWHERE to be "available", and VIIRS swaths
# do not cover every location on every pass, so the closest overpass in time can be a
# regional blank. Measured 2013-2023 (under "closest"): trying the 2nd-closest
# recovers 17 of 111 otherwise-missing events and the 3rd recovers 1 more; the 4th
# through 6th recover nothing, so 3 is where the returns stop. Also the step budget
# "before2_after2" (D-66) draws its 2+2 overpasses from, so it cannot be lowered
# below 2 without breaking that method.
NEARBY_MAX_STEPS     = 3
# D-66: four different answers to "which overpass(es) count as 'nearby in time'?",
# shipped side by side (pyrocb/nearby.py has the full reasoning for each) so a
# reviewer can compare them directly rather than the pipeline committing to one.
#   "closest"          nearest overpass, before OR after onset - no causal claim.
#   "before"           nearest AT-OR-BEFORE onset; falls through to AFTER-onset
#                      overpasses only once every before-onset option is exhausted
#                      (i.e. for a blow-up fire, structurally - not a separate flag).
#   "before_after"     nearest before AND nearest after, pooled into one search.
#   "before2_after2"   nearest TWO before AND nearest TWO after, pooled.
NEARBY_METHODS       = ("closest", "before", "before_after", "before2_after2")
# Map-only smoothing (D-67) for the nb2 "active fire used for centering" polygons:
# close_m bridges gaps between separate VIIRS pixel rectangles into solid blobs
# before simp_m simplifies them - plain simplify() alone barely helps a shape that
# IS hundreds of disjoint small rectangles. Needed once before2_after2 (up to 4
# pooled overpasses) and 15/20 km radii made these median 1000+ vertices each, x12
# (radius, method) combos. Map-display geometry only; area_km2/n_pixels in the
# attribute table come from the un-smoothed polygon in pyrocb/nearby.py.
NB2_CLOSE_M          = 200.0
NB2_SIMPLIFY_M       = 30.0
NFP_WINDOWS_H       = (6, 12, 24)
# Active-burn polygon: NFP pixel footprints within +/- W hours of onset. Symmetric,
# because the ~10 h nightly blind window means the run that fed the plume is usually
# first seen by the overpass AFTER onset (FINDINGS F3).
# Active-burn window grid: (hours BEFORE onset, hours AFTER onset). The map slider
# steps over this grid; the GeoPackage carries a polygon for every combination.
# Asymmetric matters - a blow-up fire (F4) is only seen AFTER onset, so a pure-post
# window is the honest choice there.
ACTIVE_STEPS_H     = (0, 6, 12, 18, 24)
ACTIVE_MIN_SPAN_H  = 6
ACTIVE_WINDOWS     = tuple((a, b) for a in ACTIVE_STEPS_H for b in ACTIVE_STEPS_H
                           if a + b >= ACTIVE_MIN_SPAN_H)      # 24 combinations
ACTIVE_HEADLINE    = (12, 12)
# THE ANALYSIS POLYGON. The active-burn window the covariates for the pyroCb model
# (ERA5 temperature/humidity/wind, terrain slope, fire intensity, land cover ...) are
# sampled over, and the window the map's "Analysis polygon" layer opens on. It is a
# window of ACTIVE_WINDOWS, so the map can step to any other one client-side without
# new data - this constant only sets the default and the documented choice.
ANALYSIS_WINDOW    = ACTIVE_HEADLINE
# Windows also written as columns in the attribute table: the headline gets the full
# set of metrics, the rest just area + pixel count.
# (6, 0) is the "burn area in the 6 h before onset" polygon - the geometry the
# pyroCb-info land-cover fractions are sampled over (LC_ACTIVE_WINDOW).
ACTIVE_TABLE_WINDOWS = ((12, 12), (6, 6), (24, 24), (0, 24), (24, 0), (6, 0))

# Localisation. NFP pixels in the window are unioned, near-neighbours bridged across
# swath gaps, then only the cluster containing (or nearest, within the cap) the
# pyroCb's own location is kept. Without this a fire that MERGES with a distant
# neighbour later in the season drags that neighbour's pixels into the "active area",
# because NFP `mergid` is the season-final id (FINDINGS F5). Measured: this drops
# ~48% of in-window pixels across the matched-event set.
ACTIVE_BRIDGE_M      = 1500.0
ACTIVE_ATTACH_CAP_KM = 20.0

# ---------------------------------------------------------------- land cover
# 30 m Land Cover of Canada (EPSG:3979). 2015 collection for the earlier fire
# years, 2020 for LC_YEAR_SPLIT onward. Decimated to LC_CELL_M for the map layer;
# the perimeter is grown by LC_BUFFER_M so the fire's surroundings are included.
# Land cover is Copernicus C3S 300 m ONLY (user decision 2026-10-03), same product as the dataset's
# lc_* columns. The map's C3S overlays and per-pyroCb class mix are written by build_c3s_landcover.py;
# with LC_SOURCE = "c3s" the pipeline's old Land Cover of Canada (NALCMS 30 m) stages return empty, so
# a regenerated year can never bring that product back. "nalcms" restores the old behaviour.
LC_SOURCE     = "c3s"
LC_TIF_2015   = str(DATA / "LANDCOVER" / "landcover-2015-classification.tif")
LC_TIF_2020   = str(DATA / "LANDCOVER" / "landcover-2020-classification.tif")
LC_YEAR_SPLIT = 2018
LC_CELL_M     = 1000.0
LC_BUFFER_M   = 5000.0
LC_SIMPLIFY_M = 2000.0
# Near the fires land cover is a paletted PNG overlay per fire (pyrocb/lcraster.py),
# not polygons: one GPU-composited <img> per fire instead of thousands of canvas
# paths. LC_PNG_CELL_M is deliberately == BURN_CELL_M (375 m, below) and the grid is
# snapped to the same Mercator-origin convention burngrid.py uses (D-64), so a
# land-cover cell and a fire-activity cell at the same fire share exact edges instead
# of merely overlapping - keep the two constants equal if either one changes.
# LC_PNG_CAP_PX is the per-side pixel cap that keeps a continent-sized fire from
# producing a huge image (widens the cell size for that fire only, never the grid's
# alignment - the widened cell still snaps to the same Mercator-origin grid).
LC_PNG_CELL_M = 375.0
LC_PNG_CAP_PX = 640
# Coarse background land cover: away from the fires the map still shows land cover
# rather than OpenStreetMap, but decimated hard so the browser copes. Drawn under
# the 1 km fire layer, over the convex hull of the year's pyroCb fires + this collar.
LC_BG_CELL_M     = 10000.0
LC_BG_BUFFER_M   = 120000.0
LC_BG_SIMPLIFY_M = 6000.0
# The active-burn window whose land-cover composition is reported per pyroCb event
# (top 3 classes + fractions) in the map's info panel: the 6 h before onset.
LC_ACTIVE_WINDOW = (6, 0)
# Blow-up fires (F4) have no pre-onset detections, so the (6, 0) polygon is empty.
# Fall back to the first this-many hours of detected burning that IS available.
LC_BLOWUP_SPAN_H = 6.0

# ---------------------------------------------------------------- satellite imagery
# Post-fire Sentinel-2 true-colour image per fire (pyrocb/satellite.py). Each image is
# a LIVE network fetch (STAC search + windowed COG reads via Microsoft's Planetary
# Computer), unlike everything else in the pipeline. None -> every resolved year
# (matches the YEARS/SNAPSHOT_YEARS None-means-all convention); set to a tuple (e.g.
# (2023,)) to scope it down again for a faster rebuild while iterating.
# Pre-fire images were dropped (D-59): post alone, at every year, was worth more than
# both phases at one year, and the pair cost roughly the same to fetch either way.
# () -> off for every year (D-65): a deliberate, one-line, fully reversible switch -
# already-built images and cache stay on disk untouched, satellite_cached's own year
# gate (`year not in C.SAT_YEARS`) is unconditionally true for an empty tuple, and
# nothing else in the pipeline or the map needs to change to bring it back. Set back
# to None (every year) or a specific tuple to re-enable.
SAT_YEARS: tuple[int, ...] = ()
# The circle drawn "around the whole perimeter": centroid to farthest boundary point,
# times this margin, floored at SAT_RADIUS_MIN_KM so small fires keep visible context.
# Radius varies 2-160 km across 2023's matched fires (D-51), so SAT_OUT_PX - not a
# fixed ground resolution - is what keeps the output image a sane, uniform size.
SAT_RADIUS_MARGIN = 1.25
SAT_RADIUS_MIN_KM = 5.0
SAT_OUT_PX         = 1024
# How many scenes to mosaic per image, closest-date-first, before giving up on
# filling the circle. A small fire's circle fits in one tile and stops immediately
# (the mosaic loop exits as soon as `filled` is 100%); this cap only matters for the
# handful of huge complexes whose circle spans many Sentinel-2 tiles - a 320 km-wide
# circle needs on the order of 15-20 100 km tiles, not a handful.
SAT_MAX_TILES      = 20
SAT_MAX_CLOUD      = 40.0     # scene-level eo:cloud_cover ceiling for a candidate
# Search window in days after the fire's own ted_* (official end): offset forward so
# the scene shows dead ground, not still-active fire/smoke. Candidates are ranked by
# closeness to that anchor date, not by cloud cover (D-52) - so a narrower window
# matters less here than it would under pure cloud-minimisation.
SAT_POST_WINDOW_D  = (10, 75)
# Delivery (D-53): images are JPEG files on disk, not embedded in fire_pcb_map.html -
# at 1024 px true-colour even a circle-cropped PNG is 1-2 MB, and ~70 fires x 2 phases
# would multiply the map's size 10x. SAT_IMG_DIR is relative to ROOT so the map's
# <img src> stays a same-directory relative path wherever the output is copied.
SAT_IMG_DIR        = str(ROOT / "sat_images" / "{y}")
SAT_JPEG_QUALITY   = 82

# -------------------------------------------------------- pyroCb cloud imagery (GOES)
# A snapshot of the actual pyroCb - not the fire - at the moment it happened, from
# GOES-16/17/18/19 visible imagery (pyrocb/goes_cloud.py). Two hard coverage limits,
# both physical, not configurable: GOES only sees the Western Hemisphere, and was not
# operational before Dec 2017. ~252 of 484 events are within reach.
GOES_MIN_YEAR          = 2018
GOES_EXCLUDE_REGIONS   = ("Northern Asia",)   # no geostationary coverage on Planetary Computer
# The circle radius. D-69: centred on the AS-REPORTED point (not v2 - that's what's
# being evaluated, so it can't also be the reference frame), widened from 30 to 60 km
# so every (radius, method) v2 position - up to 20 km from the report - sits well
# inside with room to compare against the actual cloud, rather than at the edge.
GOES_RADIUS_KM         = 60.0
GOES_OUT_PX            = 512
# Search window: earliest DAYLIGHT scene at least GOES_POST_WINDOW_H[0] after onset
# (the pyroCb needs time to develop into a mature, visible cloud) and at most [1]
# after (still the same event, not a later day's unrelated weather).
GOES_POST_WINDOW_H     = (0.0, 12.0)
# Hard cap on candidate-scene READ attempts per event (D-85/D-86 lineage): starting the
# search window at onset (0 h) means a dusk/night onset can walk through many darkening
# candidates before GOES_DAYLIGHT_LOCAL_H/GOES_MIN_DAYLIGHT_REFL reject them, each a real
# network fetch. This bounds the worst-case cost per event without changing which scene
# wins when one is found well inside the cap.
GOES_MAX_READ_ATTEMPTS = 10
# Local-solar-hour band treated as "daylight enough" for a visible-band read, before
# ever fetching a scene - generous, because boreal summer days are long at these
# latitudes. GOES_MIN_DAYLIGHT_REFL is the post-read sanity floor (median reflectance
# of valid pixels) that catches the dawn/dusk edge this coarse hour band lets through.
GOES_DAYLIGHT_LOCAL_H  = (5.0, 21.0)
GOES_MIN_DAYLIGHT_REFL = 0.05
GOES_JPEG_QUALITY      = 82
GOES_IMG_DIR           = str(ROOT / "cloud_images" / "{y}")
# Time series around onset: one frame per offset, so the map can loop through the
# pyroCb's own development instead of showing a single still. Offsets are hours from
# onset; GOES_SERIES_TOL_H is how far from the requested offset a scene may sit and
# still be accepted (full-disk cadence is ~10-15 min, so this is generous). Frames that
# find nothing - usually the late ones, past local dusk - are simply absent, and the
# map's player skips them rather than pretending the series is regular.
GOES_SERIES_OFFSETS_H  = (-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0)
GOES_SERIES_TOL_H      = 1.0
GOES_SERIES_PX         = 384          # smaller than the single still: 7x the files

# ---------------------------------------------------------------- fire-activity grid
# Both animation layers (the +/-24 h window and the whole-fire event) ship as a
# delta-coded, time-sorted grid of burned cells at ONE ground resolution - see
# pyrocb/burngrid.py. 375 m is a single VIIRS pixel, i.e. native.
# LC_PNG_CELL_M (above) is kept equal to this on purpose (D-64) - the land-cover PNG
# overlay grid and this fire-activity grid share the same Mercator-origin snapping,
# so their cells tile against each other exactly rather than merely overlapping.
#
# This replaced a tiered decimation that inflated each kept detection to its dedup
# cell, so the quiet spans between pyroCb onsets were drawn as 3 km blocks - blocks
# that swallowed whole perimeters on small fires. Measured on 2023: the grid is
# 0.64 MB base64 at 375 m against 1.7 MB for the old 1-3 km version.
BURN_CELL_M = 375.0

# Map scope. True -> the map carries only fires credited with a pyroCb (and their
# ignitions), which is the analysis population; False -> every tracked fire, shipped as
# a togglable "All fire activity" layer (off by default in the UI). Measured cost of
# False: ~4 MB gzip/year for perimeters+ignitions (most of it is a handful of huge
# complexes - simplified/precision-snapped small-fire polygons cost almost nothing), so
# the fear of "~170k polygons, 27 MB" from when this was last measured did not hold up.
MAP_ONLY_PCB_FIRES  = False
# A non-pyroCb fire needs at least this many VIIRS pixels to get a baked-in whole-fire
# growth animation when someone clicks its ignition. Below it, an ignition still opens
# the fire panel with basic stats, just no player - most of the ~170k tracked fires are
# 1-5 pixel noise detections with no growth story to show. Cost is dominated by the
# fires that already clear this bar regardless of threshold value (20 vs 50 pixels
# changed 2023's payload by under 5%), so this is chosen for "worth animating", not size.
ANIM_MIN_PIXELS     = 20
# If fewer than this fraction of fires reproduce their n_pixels from the NFP id join,
# the year's NFP labels are stale (2021) and pixels are re-attributed by geometry.
NFP_ID_MIN_MATCH    = 0.95
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
    2024: "#5b1a4a", 2025: "#2d1b4e",     # extension years (the map re-spreads its own RdYlBu ramp over all years)
}
# Magma. Sequential and perceptually uniform, so "later" reads as "brighter" without the
# false category boundaries the old year-colour ramp implied (it was just the per-year
# palette reused, which is qualitative). Drives day-of-burn colouring and the 12 h
# snapshot growth bands; the animation's burn-age ramp (AGE_STOPS in the map) matches.
DOY_RAMP = ["#000004", "#0c0927", "#231151", "#410f75", "#5f187f", "#7b2382",
            "#982d80", "#b73779", "#d3436e", "#eb5760", "#f8765c", "#fd9668",
            "#feb77e", "#fddc9e", "#fcfdbf"]
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
    # by file NAME (not full path), so moving the project folder does not invalidate every cache
    out = []
    for p in paths:
        p = Path(p)
        try:
            st = p.stat()
            out.append([p.name, int(st.st_mtime), st.st_size])
        except OSError:
            out.append([p.name, None, None])
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
