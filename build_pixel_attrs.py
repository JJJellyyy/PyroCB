"""Every pixels.parquet feature on the map: per-cell attributes for the burn grids, fire-level
means for the polygons, and the variable dictionary for the "how was this calculated" lookup.

The map draws every fire as a time-sorted, delta-coded grid of cells (burngrid.py): one grid
per whole fire (`sim` = pyroCb-linked, `simAll` = the "All fire activity" set) and one per pyroCb
detection window (`pix`). This script lines the parquet up with those grids and writes:

    pixattr/meta.json                    every variable: label, unit, group, ramp, scale, and the
                                         dictionary text (definition / calculation / caveats / source)
    pixattr/{year}_{pcb|all}.json        {"sim"|"pix": {key: [offset, n]}} - where each grid's cells sit
    pixattr/v/{col}/{year}_{pcb|all}.bin.gz
                                         ONE variable, every grid of that year, 4-bit codes packed two
                                         per byte (high nibble first): 0..14 = the variable's scale,
                                         15 = no value. Same cell order as the grids' di/dj/dt.
    pixattr/v/spread_bearing/...       the one exception: 8-bit true bearing (spread_ux/uy merged), drawn
                                         as arrows; no fire-level mean
    pixattr/f/{year}_ids.json            the fire (mergid) order of that year's fire-level files
    pixattr/f/{col}/{year}.bin.gz        ONE variable's fire-level mean (mode for categorical) per
                                         fire, uint8 (0..254, 255 = none), on its own fire scale

One file per variable is the whole trick: the browser fetches only the variable that is picked
(~0.05-0.6 MB a year), so shipping all ~155 costs nothing until one is used. 4-bit codes are a
deliberate downscale for display (16 colour steps; patterns stay visible) - they cut size ~40 %
against 8-bit. Re-ordering cells spatially was measured and gains nothing (values are close to
incompressible at 375 m), so cells stay in grid order and need no permutation.

The join is spatial, not by array index: the parquet holds a slightly different pixel set than
the NFP files the grids were built from. A grid cell takes the values of the earliest-detected
parquet pixel inside it; a detection-window entry (cell, time) prefers the pixel with that exact
time. Cells with no parquet pixel get "no value".

Variables and their documentation come from pixel_variable_dictionary.csv (every row not flagged
"not a feature"). Ids, geometry, labels and metadata are left out.

Run:  python build_pixel_attrs.py [year ...]      (default: every year in the map)
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import re
import shutil
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyproj import Geod, Transformer

from pyrocb import config as C

# the dataset lives once, in the parent project (merged layout)
PARQUET = str(C.ROOT.parent / "OUTPUT_DATA" / "datacube" / "pixels.parquet")
DICTIONARY = str(C.ROOT.parent / "OUTPUT_DATA" / "datacube" / "pixel_variable_dictionary.csv")
HTML = "fire_pcb_map.html"
OUT = "pixattr"
NODATA = 255
CHUNK = 40                     # columns per pass (memory: rows x CHUNK float64)

# ---- variable catalogue -------------------------------------------------------------------
CLASS_NAMES = {"needle_evergreen": "Evergreen needleleaf", "needle_deciduous": "Larch (deciduous needleleaf)",
               "broad_deciduous": "Broadleaf", "mixed_forest": "Mixed forest", "shrub": "Shrub",
               "grass_herb": "Grass / herb", "sparse_moss_lichen": "Sparse / lichen / moss",
               "wetland": "Wetland", "cropland": "Cropland", "nonfuel": "Non-fuel"}
BASE_LABELS = {
    "frp": "FRP", "DS": "Pixel size, along-scan", "DT": "Pixel size, along-track", "pix_area_km2": "Pixel area",
    "nb_n_1p5km": "Active neighbours (1.5 km)", "nb_n_5km": "Active neighbours (5 km)",
    "nb_frp_1p5km": "Neighbour FRP (1.5 km)", "nb_frp_5km": "Neighbour FRP (5 km)",
    "nb_elong_5km": "Fire-front elongation (5 km)", "burned_frac_1p5km": "Burned this season (1.5 km)",
    "spread_dist_km": "Spread distance", "spread_dt_h": "Time since predecessor", "spread_rate_kmh": "Spread rate",
    "spread_ux": "Spread direction, grid x", "spread_uy": "Spread direction, grid y",
    "ahead_burned_1p5km": "Already burned 1.5 km ahead", "ahead_burned_3km": "Already burned 3 km ahead",
    "elev": "Elevation", "slope_deg": "Slope", "aspect_sin": "Aspect, eastness", "aspect_cos": "Aspect, northness",
    "tpi_1p5km": "TPI (1.5 km)", "tpi_5km": "TPI (5 km)", "rough_1p5km": "Roughness (1.5 km)",
    "rough_5km": "Roughness (5 km)", "slope_along_spread": "Slope along spread",
    "elev_gain_ahead_1p5km": "Elevation gain 1.5 km ahead", "elev_gain_ahead_3km": "Elevation gain 3 km ahead",
    "slope_ahead_1p5km": "Slope 1.5 km ahead", "slope_ahead_3km": "Slope 3 km ahead",
    "lc_class": "Land-cover class (incl. water)", "lc_group": "Land-cover group",
    "shannon_1p5km": "Shannon diversity (1.5 km)", "shannon_5km": "Shannon diversity (5 km)",
    "simpson_1p5km": "Simpson diversity (1.5 km)", "simpson_5km": "Simpson diversity (5 km)",
    "n_groups_1p5km": "Fuel groups present (1.5 km)", "n_groups_5km": "Fuel groups present (5 km)",
}
BASE_LABELS.update({
    "dist_road_km": "Distance to nearest road or track", "dist_road_public_km": "Distance to nearest public road",
    "aspect_deg": "Aspect (downslope direction)", "tpi_375m": "TPI (375 m, inside the pixel)", "tpi_2km": "TPI (~2 km)",
    "rough_375m": "Roughness (375 m, inside the pixel)", "rough_2km": "Roughness (~2 km)",
    "shannon_2km": "Shannon diversity (~2 km)", "n_groups_2km": "Fuel groups present (~2 km)",
})
WEATHER_LABELS = {
    "t2m_C": "Air temperature (2 m)", "d2m_C": "Dew point (2 m)", "rh2m": "Relative humidity (2 m)", "vpd_kPa": "Vapour pressure deficit",
    "wind10": "Wind speed (10 m)", "wind100": "Wind speed (100 m)", "shear_10_100": "Wind shear 10-100 m", "u10": "Wind, eastward (10 m)",
    "v10": "Wind, northward (10 m)", "blh": "Boundary layer height", "cape": "CAPE", "cin": "CIN", "kx": "K-index",
    "totalx": "Total totals index", "tcwv": "Column water vapour", "tp_mm": "Precipitation (1 h)", "swvl1": "Soil moisture (0-7 cm)",
    "sp_hPa": "Surface pressure", "t850_C": "Temperature 850 hPa", "t700_C": "Temperature 700 hPa", "t500_C": "Temperature 500 hPa",
    "td850_C": "Dew point 850 hPa", "td700_C": "Dew point 700 hPa", "q850": "Specific humidity 850 hPa",
    "q700": "Specific humidity 700 hPa", "q500": "Specific humidity 500 hPa", "wind850": "Wind speed 850 hPa",
    "wind700": "Wind speed 700 hPa", "wind500": "Wind speed 500 hPa", "lapse_850_700": "Temperature drop 850-700 hPa",
    "lapse_700_500": "Temperature drop 700-500 hPa", "haines_mid": "Haines index (mid)", "c_haines": "Continuous Haines index",
    "wind_alignment": "Wind along spread direction",
    "fwi": "Fire Weather Index", "ffmc": "Fine fuel moisture code", "dmc": "Duff moisture code", "dc": "Drought code",
    "isi": "Initial spread index", "bui": "Build-up index", "dsr": "Daily severity rating",
}
GROUP_ORDER = ["Fire behaviour", "Spread", "Terrain", "Access", "Land cover", "Fuel share around", "Fuel change behind",
               "Fuel diversity", "Weather", "Weather change", "Fire weather indices",
               # v1 groups, kept for the original table
               "Burned ground", "Terrain ahead", "Fuel ahead / behind", "Fuel change ahead", "Change along fire path"]
# land-cover group labels, same order as lc_group (and LC_CLASSES in the map)
LC_GROUPS = ["Needleleaf evergreen", "Needleleaf deciduous (larch)", "Broadleaf", "Mixed forest", "Shrub",
             "Grass / herb", "Sparse / lichen / moss", "Wetland", "Cropland", "Non-fuel (urban, bare, ice)"]
# LCCS code -> map land-cover index 0..10 (lc_group order, 10 = water) for lc_class
LCCS_TO_IX = {}
for _ix, _codes in enumerate([(70, 71, 72), (80, 81, 82), (50, 60, 61, 62), (90, 100), (110, 120, 121, 122), (130,),
                              (140, 150, 151, 152, 153), (160, 170, 180), (10, 11, 12, 20, 30, 40),
                              (190, 200, 201, 202, 220), (210,)]):
    for _c in _codes:
        LCCS_TO_IX[_c] = _ix


def _label(col):
    if col.endswith("_dpath"):
        return _label(col[:-6]) + " · Δ along path"
    if col in BASE_LABELS:
        return BASE_LABELS[col]
    m = re.match(r"^(.*)_dtime(12|24)h$", col)
    if m and m.group(1) in WEATHER_LABELS:
        return f"{WEATHER_LABELS[m.group(1)]} · change over {m.group(2)} h"
    m = re.match(r"^(fwi|ffmc|dmc|dc|isi|bui|dsr)_d3d$", col)
    if m:
        return f"{WEATHER_LABELS[m.group(1)]} · change over 3 days"
    if col in WEATHER_LABELS:
        return WEATHER_LABELS[col]
    m = re.match(r"^frac_(.+?)_(3x3|5x5)(_dbehind)?$", col)
    if m and m.group(1) in CLASS_NAMES:
        cls, win = CLASS_NAMES[m.group(1)], {"3x3": "3x3 pixels, ~1.1 km", "5x5": "5x5 pixels, ~1.9 km"}[m.group(2)]
        return f"{cls} share change behind ({win})" if m.group(3) else f"{cls} share ({win})"
    m = re.match(r"frac_(.+?)_(ahead_1p5km_change|ahead_3km_change|ahead_1p5km|ahead_3km|behind_1p5km|1p5km|5km)$", col)
    if m:
        cls, where = CLASS_NAMES.get(m.group(1), m.group(1)), m.group(2)
        return {"1p5km": f"{cls} share (1.5 km)", "5km": f"{cls} share (5 km)",
                "ahead_1p5km": f"{cls} share 1.5 km ahead", "ahead_3km": f"{cls} share 3 km ahead",
                "behind_1p5km": f"{cls} share 1.5 km behind",
                "ahead_1p5km_change": f"{cls} share change, 1.5 km ahead",
                "ahead_3km_change": f"{cls} share change, 3 km ahead"}[where]
    return col


def _group(r):
    col, fam, blk = r["variable"], r["family"], r["block"]
    if fam == "path delta":
        return "Change along fire path"
    if fam == "access":
        return "Access"
    if fam == "fuel change behind":
        return "Fuel change behind"
    if fam == "weather":
        return "Weather"
    if fam == "weather change (same location)":
        return "Weather change"
    if fam == "fire weather index":
        return "Fire weather indices"
    if fam == "fuel share":
        return "Fuel share around"
    if fam in ("fuel ahead", "fuel behind"):
        return "Fuel ahead / behind"
    if fam == "fuel change ahead":
        return "Fuel change ahead"
    if fam == "diversity":
        return "Fuel diversity"
    if col in ("burned_frac_1p5km", "ahead_burned_1p5km", "ahead_burned_3km"):
        return "Burned ground"
    if col in ("lc_class", "lc_group"):
        return "Land cover"
    if blk == "terrain":
        return "Terrain ahead" if fam == "ahead" else "Terrain"
    if blk == "spread" or col == "spread_dt_h":
        return "Spread"
    return "Fire behaviour"


SIGNED_TERRAIN = ("tpi_1p5km", "tpi_5km", "tpi_375m", "tpi_2km", "slope_along_spread", "elev_gain_ahead_1p5km", "elev_gain_ahead_3km")


def _kind_ramp(r):
    """(kind, ramp) - kind: lin | log | div (symmetric about 0) | cat; ramp names match RAMPS in the map."""
    col, fam = r["variable"], r["family"]
    if col in ("lc_class", "lc_group"):
        return "cat", "cat"
    if col.startswith("dist_road"):
        return "log", "index"
    if fam in ("path delta", "fuel change ahead", "fuel change behind", "weather change (same location)") \
            or col in ("spread_ux", "spread_uy", "aspect_sin", "aspect_cos") or col.endswith("_d3d"):
        return "div", "delta"
    if col == "aspect_deg":
        return "lin", "index"
    if col in SIGNED_TERRAIN:
        return "div", "signed"
    if col in ("frp", "nb_frp_1p5km", "nb_frp_5km"):
        return "log", "fire"
    if col == "elev":
        return "lin", "terrain"
    if col.startswith(("slope", "rough")):
        return "lin", "slope"
    if col.startswith("frac_nonfuel"):
        return "lin", "grey"
    if col.startswith("frac_"):
        return "lin", "veg"
    if fam == "diversity" or col in ("DS", "DT", "pix_area_km2", "nb_elong_5km", "spread_dt_h"):
        return "lin", "index"
    return "lin", "fire"


def catalogue():
    d = pd.read_csv(DICTIONARY)
    d = d[~d["auto_flags"].fillna("").str.contains("not a feature")]
    clean = lambda v: "" if pd.isna(v) else str(v)
    out = []
    for _, r in d.iterrows():
        kind, ramp = _kind_ramp(r)
        unit = clean(r["unit"])
        out.append({"col": r["variable"], "label": _label(r["variable"]), "unit": "" if unit == "-" else unit,
                    "kind": kind, "ramp": ramp, "group": _group(r),
                    "def": clean(r["definition"]), "calc": clean(r["calculation"]),
                    "caveats": clean(r["caveats"]), "source": clean(r["source"]),
                    "flags": clean(r["auto_flags"])})
    # spread_ux / spread_uy are one vector: merged into a true bearing, drawn as arrows (D-94)
    out = [m for m in out if m["col"] not in ("spread_ux", "spread_uy")]
    out.append({"col": "spread_bearing", "label": "Spread direction", "unit": "deg", "kind": "dir", "ramp": "fire",
                "group": "Spread", "bits": 8, "rate": "spread_rate_kmh", "lo": 0.0, "hi": 360.0,
                "def": "Direction the fire is moving at the pixel: true bearing, clockwise from north.",
                "calc": "The unit spread vector (spread_ux, spread_uy: this pixel minus its predecessor, in the "
                        "EPSG:3571 grid frame) rotated into true north - both ends projected to lon/lat and the "
                        "geodesic forward azimuth taken (pyproj). Stored at 8 bits (1.4 deg steps).",
                "caveats": "NaN without a predecessor (new lobes). Drawn as arrows on a screen grid: each arrow is "
                           "the mean direction of the cells under it, its length how consistent they are, its "
                           "colour their mean spread rate.",
                "source": "derived from Fire Atlas NFP (spread_ux, spread_uy)", "flags": ""})
    out.sort(key=lambda m: (GROUP_ORDER.index(m["group"]), m["label"]))
    return out


_GEOD = Geod(ellps="WGS84")
_TF_LL = Transformer.from_crs(C.METRIC_CRS, 4326, always_xy=True)


def spread_bearing(x, y, ux, uy):
    """True bearing (deg clockwise from north) of a grid-frame unit vector at grid point (x, y)."""
    out = np.full(len(x), np.nan)
    ok = np.isfinite(ux) & np.isfinite(uy) & ((ux != 0) | (uy != 0))
    if ok.any():
        lon0, lat0 = _TF_LL.transform(x[ok], y[ok])
        lon1, lat1 = _TF_LL.transform(x[ok] + ux[ok] * 100.0, y[ok] + uy[ok] * 100.0)
        az, _, _ = _GEOD.inv(lon0, lat0, lon1, lat1)
        out[ok] = np.mod(az, 360.0)
    return out


# ---- helpers ----------------------------------------------------------------------------------
_TF = Transformer.from_crs(C.METRIC_CRS, 3857, always_xy=True)


def unpack(b64):
    return json.loads(gzip.decompress(base64.b64decode(b64)).decode())


def load_map():
    """The map's DATA object: map_src/data.json (see map_tool.py) if present, else the built HTML."""
    if os.path.exists("map_src/data.json"):
        return json.load(open("map_src/data.json", encoding="utf-8"))
    text = open(HTML, "r", encoding="utf-8").read()
    i = text.index("const DATA=")
    data, _ = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])
    return data


def _prep(col, v):
    """Raw float values ready for scaling (lc_class -> map land-cover index)."""
    # nullable columns (pd.NA, e.g. a flag that some years lack) -> NaN
    v = pd.to_numeric(pd.Series(v), errors="coerce").astype("float64").to_numpy()
    if col == "lc_class":
        out = np.full(v.shape, np.nan)
        ok = np.isfinite(v)
        out[ok] = [LCCS_TO_IX.get(int(c), np.nan) for c in v[ok]]
        return out
    return v


def _fit(m, s, q=(.01, .99)):
    """lo/hi for a variable from a sample of (prepared) values."""
    s = s[np.isfinite(s)]
    if m["kind"] == "cat":
        return 0.0, float(len(m["labels"]) - 1)
    if not len(s):
        return 0.0, 1.0
    if m["kind"] == "log":
        s = np.log10(np.clip(s, 0.1, None))
    lo, hi = float(np.quantile(s, q[0])), float(np.quantile(s, q[1]))
    if m["kind"] == "div":
        a = max(abs(lo), abs(hi)) or 1.0
        return -a, a
    return (lo, hi) if hi > lo else (lo, lo + 1.0)


def _quantise(v, kind, lo, hi):
    v = np.asarray(v, "float64")
    bad = ~np.isfinite(v)
    if kind == "cat":
        q = np.clip(np.rint(np.nan_to_num(v)), lo, hi)
    else:
        if kind == "log":
            v = np.log10(np.clip(v, 0.1, None))
        q = np.rint(np.clip((np.nan_to_num(v) - lo) / (hi - lo), 0, 1) * 254)
    q[bad] = NODATA
    return q.astype(np.uint8)


def pack4(q8):
    """8-bit codes (0..254, 255 none) -> 4-bit codes (0..14, 15 none), two per byte.
    Categorical codes are class indices (<= 10) and pass through unchanged."""
    q = np.where(q8 == NODATA, 15, q8).astype(np.uint8)
    return _nibbles(q)


def pack4_scaled(q8):
    q = np.where(q8 == NODATA, 15, np.rint(q8.astype(np.float32) * 14 / 254)).astype(np.uint8)
    return _nibbles(q)


def _nibbles(q):
    if len(q) % 2:
        q = np.append(q, np.uint8(15))
    return (q[0::2] << 4 | q[1::2]).astype(np.uint8).tobytes()


def _grid_cells(g):
    c = g["c"]
    i = np.cumsum(np.asarray(g["di"], np.int64)) + int(round(g["ox"] / c))
    j = np.cumsum(np.asarray(g["dj"], np.int64)) + int(round(g["oy"] / c))
    t = np.cumsum(np.asarray(g["dt"], np.int64))
    return i, j, t


def grid_rows(grid, px, ix, onset_min=None):
    """Parquet row (index into px) feeding each cell of `grid`, -1 where none. `ix` = candidate rows."""
    n = len(grid["di"])
    out = np.full(n, -1, np.int64)
    if ix is None or not len(ix):
        return out
    i, j, t = _grid_cells(grid)
    mx, my = _TF.transform(px["x"].to_numpy(float)[ix], px["y"].to_numpy(float)[ix])
    c = grid["c"]
    pkey = np.floor(np.asarray(mx) / c).astype(np.int64) * (1 << 32) + np.floor(np.asarray(my) / c).astype(np.int64)
    ptime = px["_min"].to_numpy(np.int64)[ix]
    prel = ptime - onset_min if onset_min is not None else ptime
    o = np.lexsort((ptime, pkey))
    skey, stime = pkey[o], prel[o]
    qkey = i * (1 << 32) + j
    first = np.searchsorted(skey, qkey, side="left")
    ok = (first < len(skey)) & (skey[np.minimum(first, len(skey) - 1)] == qkey)
    pick = np.where(ok, first, -1)
    if onset_min is not None and ok.any():           # a window entry prefers its exact detection time
        last = np.searchsorted(skey, qkey, side="right")
        for k in range(8):
            cand = first + k
            hit = ok & (cand < last) & (stime[np.minimum(cand, len(skey) - 1)] == t)
            pick = np.where(hit & (pick == first), cand, pick)
    good = pick >= 0
    out[good] = ix[o[pick[good]]]
    return out


# ---- build ------------------------------------------------------------------------------------
def fit_scales(meta):
    """Per-cell scale from a sample of the parquet (1st/99th percentile; symmetric for div)."""
    pf = pq.ParquetFile(PARQUET)
    take = list(range(0, pf.num_row_groups, 4))
    cols = [m for m in meta if m["kind"] != "dir"]
    for a in range(0, len(cols), CHUNK):
        chunk = cols[a:a + CHUNK]
        df = pf.read_row_groups(take, columns=[m["col"] for m in chunk]).to_pandas()
        for m in chunk:
            m["lo"], m["hi"] = (round(x, 5) for x in _fit(m, _prep(m["col"], df[m["col"]].to_numpy())))
        del df


def build_year(y, data, meta, fire_means):
    t0 = time.time()
    base = pq.read_table(PARQUET, columns=["mergid", "x", "y", "datetime_utc"], filters=[("year", "=", y)]).to_pandas()
    base["_min"] = base["datetime_utc"].astype("datetime64[ns]").astype("int64") // 60_000_000_000
    by_fire = base.groupby("mergid").indices
    d = data[str(y)]
    pcbp = {f["properties"]["pcb_id"]: f["properties"] for f in unpack(d["pcb"])["features"]}

    # 1. which parquet row feeds each cell, per grid - computed once, reused for every variable
    layout = {"pcb": ({"sim": {}, "pix": {}}, []), "all": ({"sim": {}}, [])}
    for k, v in unpack(d["sim"]).items():
        layout["pcb"][0]["sim"][k] = len(v["g"]["di"])
        layout["pcb"][1].append(grid_rows(v["g"], base, by_fire.get(int(k))))
    for k, g in unpack(d["pix"]).items():
        pp = pcbp.get(k)
        ix = by_fire.get(int(pp["mergid"])) if pp and pp.get("mergid") is not None else None
        onset = pd.Timestamp(pp["datetime"]).value // 60_000_000_000 if pp else None
        layout["pcb"][0]["pix"][k] = len(g["di"])
        layout["pcb"][1].append(grid_rows(g, base, ix, onset))
    for k, v in unpack(d["simAll"]).items():
        layout["all"][0]["sim"][k] = len(v["g"]["di"])
        layout["all"][1].append(grid_rows(v["g"], base, by_fire.get(int(k))))
    rows = {}
    for kind, (idx, parts) in layout.items():
        rows[kind] = np.concatenate(parts) if parts else np.zeros(0, np.int64)
        off, out = 0, {}
        for sec in idx:
            out[sec] = {}
            for key, n in idx[sec].items():
                out[sec][key] = [off, n]
                off += n
        json.dump(out, open(f"{OUT}/{y}_{kind}.json", "w"), separators=(",", ":"))

    # 2. variables in column chunks: per-cell 4-bit files + fire-level means
    mergids = np.array(sorted(by_fire))
    fpos = {m: i for i, m in enumerate(mergids)}
    frow = np.array([fpos[m] for m in base["mergid"].to_numpy()])
    json.dump([int(m) for m in mergids], open(f"{OUT}/f/{y}_ids.json", "w"), separators=(",", ":"))
    nbytes = {"pcb": 0, "all": 0}
    cols = [m for m in meta if m["kind"] != "dir"]
    for a in range(0, len(cols), CHUNK):
        chunk = cols[a:a + CHUNK]
        df = pq.read_table(PARQUET, columns=[m["col"] for m in chunk], filters=[("year", "=", y)]).to_pandas()
        for m in chunk:
            v = _prep(m["col"], df[m["col"]].to_numpy())
            q = _quantise(v, m["kind"], m["lo"], m["hi"])
            os.makedirs(f"{OUT}/v/{m['col']}", exist_ok=True)
            for kind, r in rows.items():
                cell = np.full(len(r), NODATA, np.uint8)
                ok = r >= 0
                cell[ok] = q[r[ok]]
                blob = gzip.compress(pack4(cell) if m["kind"] == "cat" else pack4_scaled(cell), 9, mtime=0)
                open(f"{OUT}/v/{m['col']}/{y}_{kind}.bin.gz", "wb").write(blob)
                nbytes[kind] += len(blob)
            # fire level: mean of the raw values (mode for categorical); scaled once all years are in
            s = pd.Series(v)
            if m["kind"] == "cat":
                fm = s.dropna().groupby(frow[s.notna().to_numpy()]).agg(lambda x: x.mode().iat[0])
            else:
                fm = s.groupby(frow).mean()
            fire_means.setdefault(m["col"], {})[y] = fm.reindex(range(len(mergids))).to_numpy()
        del df
    # spread direction: 8-bit bearing codes (0..253 over 0..360 deg, 255 = none), no fire-level mean
    uv = pq.read_table(PARQUET, columns=["spread_ux", "spread_uy"], filters=[("year", "=", y)]).to_pandas()
    b = spread_bearing(base["x"].to_numpy(float), base["y"].to_numpy(float),
                       uv["spread_ux"].to_numpy(float), uv["spread_uy"].to_numpy(float))
    q = np.where(np.isfinite(b), np.mod(np.rint(np.nan_to_num(b) / 360.0 * 254), 254), NODATA).astype(np.uint8)
    os.makedirs(f"{OUT}/v/spread_bearing", exist_ok=True)
    for kind, r in rows.items():
        cell = np.full(len(r), NODATA, np.uint8)
        ok = r >= 0
        cell[ok] = q[r[ok]]
        blob = gzip.compress(cell.tobytes(), 9, mtime=0)
        open(f"{OUT}/v/spread_bearing/{y}_{kind}.bin.gz", "wb").write(blob)
        nbytes[kind] += len(blob)
    print(f"  {y}: {len(rows['pcb']):>7} pyroCb-grid cells {nbytes['pcb']/1e6:5.1f} MB | "
          f"{len(rows['all']):>7} background cells {nbytes['all']/1e6:5.1f} MB | "
          f"{len(mergids)} fires  [{time.time()-t0:.0f}s]", flush=True)


def write_fire_means(meta, fire_means):
    """Fire means are narrower than pixels, so their scale is the 2nd/98th percentile of the means."""
    for m in meta:
        if m["kind"] == "dir":
            continue
        per_year = fire_means[m["col"]]
        allv = np.concatenate(list(per_year.values()))
        m["flo"], m["fhi"] = (round(x, 5) for x in _fit(m, allv, (.02, .98)))
        os.makedirs(f"{OUT}/f/{m['col']}", exist_ok=True)
        for y, v in per_year.items():
            blob = gzip.compress(_quantise(v, m["kind"], m["flo"], m["fhi"]).tobytes(), 9, mtime=0)
            open(f"{OUT}/f/{m['col']}/{y}.bin.gz", "wb").write(blob)


def main():
    data = load_map()
    years = [int(a) for a in sys.argv[1:] if a.isdigit()] or sorted(int(y) for y in data)
    meta = catalogue()
    for m in meta:
        if m["kind"] == "cat":
            m["labels"] = LC_GROUPS + (["Water"] if m["col"] == "lc_class" else [])
    print(f"{len(meta)} variables from {DICTIONARY}", flush=True)
    for sub in ("v", "f"):
        shutil.rmtree(f"{OUT}/{sub}", ignore_errors=True)
        os.makedirs(f"{OUT}/{sub}", exist_ok=True)
    for old in os.listdir(OUT):                       # the previous all-variables-in-one-file layout
        if old.endswith(".bin.gz") or (old.startswith("fires_") and old.endswith(".json")):
            os.remove(f"{OUT}/{old}")
    fit_scales(meta)
    fire_means = {}
    for y in years:
        build_year(y, data, meta, fire_means)
    write_fire_means(meta, fire_means)
    json.dump({"vars": meta, "nodata": NODATA, "bits": 4,
               "note": "per-cell: 4-bit codes 0..14 on [lo, hi] (categorical: the class index), 15 = none; "
                       "fire means: uint8 0..254 on [flo, fhi], 255 = none"},
              open(f"{OUT}/meta.json", "w"), separators=(",", ":"))
    total = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(OUT) for f in fs)
    print(f"wrote {OUT}/: {total/1e6:.0f} MB, {sum(len(fs) for _, _, fs in os.walk(OUT))} files")


if __name__ == "__main__":
    main()
