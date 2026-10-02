"""Per-cell pixel attributes for the map, built from pixels.parquet.

The map draws every fire as a time-sorted, delta-coded grid of cells (burngrid.py): one
grid per whole fire (`sim` = pyroCb-linked, `simAll` = the "All fire activity" set) and
one per pyroCb +/-24 h window (`pix`). pixels.parquet carries ~170 extra columns per
VIIRS pixel (terrain, land cover, spread, neighbourhood ...). This script lines the two
up and writes, per year, a small quantised attribute array for every grid cell:

    pixattr/meta.json            variable list + how each was quantised (lo/hi/kind)
    pixattr/{year}_pcb.bin.gz    sim + pix grids   (what a pyroCb panel needs)
    pixattr/{year}_pcb.json      {"sim": {mergid: [off, n]}, "pix": {pcb_id: [off, n]}}
    pixattr/{year}_all.bin.gz    simAll grids      (only fetched for background fires)
    pixattr/{year}_all.json      {"sim": {mergid: [off, n]}}
    pixattr/fires_{year}.json    {mergid: [V codes]} one value per variable per FIRE (mean over
                                 all its pixels; mode for categorical), quantised on the
                                 FIRE-level scale meta["vars"][k]["flo"/"fhi"] - what colours
                                 whole perimeters/ignitions on the main map

A block of n cells is V consecutive uint8 arrays of length n (variable-major, same
order as meta["vars"]), in the SAME order as the grid's di/dj/dt, so the browser indexes
attribute k of cell c as  block[k*n + c].  255 = no value.  Offsets are in cells.

The join is spatial, not by array index: the parquet holds a slightly different pixel
set than the NFP files the map grids were built from (a few % fewer pixels, simAll
fires in particular), so index alignment is not possible. A grid cell takes the values
of the earliest-detected parquet pixel inside it; a +/-24 h window entry (cell, time)
prefers the pixel with that exact time. Cells with no parquet pixel get 255.

Run:  python build_pixel_attrs.py [--fires-only] [year ...]   (default: every year in the map)
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyproj import Transformer

from pyrocb import config as C

PARQUET = "pixels.parquet"
HTML = "fire_pcb_map.html"
OUT = "pixattr"
NODATA = 255

# (column, label, unit, kind, group). kind: lin | log | div (diverging about 0) | cat.
# lo/hi come from the 1st/99th percentile of the data itself (see _fit_scales) so a few
# extreme pixels do not flatten the ramp.
VARS = [
    ("frp",                         "FRP",                       "MW",   "log", "Fire"),
    ("spread_rate_kmh",             "Spread rate",               "km/h", "lin", "Fire"),
    ("nb_n_1p5km",                  "Active neighbours (1.5 km)", "px",  "lin", "Fire"),
    ("nb_frp_5km",                  "Neighbour FRP (5 km)",      "MW",   "log", "Fire"),
    ("burned_frac_1p5km",           "Burned fraction (1.5 km)",  "",     "lin", "Fire"),
    ("ahead_burned_3km",            "Burned ahead (3 km)",       "",     "lin", "Fire"),
    ("elev",                        "Elevation",                 "m",    "lin", "Terrain"),
    ("slope_deg",                   "Slope",                     "deg",  "lin", "Terrain"),
    ("tpi_5km",                     "TPI (5 km)",                "m",    "div", "Terrain"),
    ("rough_5km",                   "Roughness (5 km)",          "m",    "lin", "Terrain"),
    ("slope_along_spread",          "Slope along spread",        "",     "div", "Terrain"),
    ("elev_gain_ahead_3km",         "Elevation gain ahead (3 km)", "m",  "div", "Terrain"),
    ("lc_group",                    "Land-cover group",          "",     "cat", "Land cover"),
    ("frac_needle_evergreen_1p5km", "Evergreen needleleaf (1.5 km)", "", "lin", "Land cover"),
    ("frac_needle_deciduous_1p5km", "Larch / deciduous needleleaf (1.5 km)", "", "lin", "Land cover"),
    ("frac_broad_deciduous_1p5km",  "Broadleaf deciduous (1.5 km)", "",  "lin", "Land cover"),
    ("frac_shrub_1p5km",            "Shrub (1.5 km)",            "",     "lin", "Land cover"),
    ("shannon_5km",                 "Land-cover diversity (5 km)", "",   "lin", "Land cover"),
]
# lc_group index -> label, same order as the frac_* columns in the parquet.
LC_GROUPS = ["Needleleaf evergreen", "Needleleaf deciduous (larch)", "Broadleaf deciduous",
             "Mixed forest", "Shrub", "Grass / herb", "Sparse / moss / lichen", "Wetland",
             "Cropland", "Non-fuel (urban, bare, water, ice)"]

_TF = Transformer.from_crs(C.METRIC_CRS, 3857, always_xy=True)


def unpack(b64):
    return json.loads(gzip.decompress(base64.b64decode(b64)).decode())


def load_map():
    text = open(HTML, "r", encoding="utf-8").read()
    i = text.index("const DATA=")
    data, _ = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])
    return data


def _fit_scales(years):
    """lo/hi per variable from a sample of the parquet (1st/99th percentile)."""
    cols = [v[0] for v in VARS]
    pf = pq.ParquetFile(PARQUET)
    take = list(range(0, pf.num_row_groups, 2))
    df = pf.read_row_groups(take, columns=cols).to_pandas()
    meta = []
    for col, label, unit, kind, group in VARS:
        s = df[col].astype("float64").replace([np.inf, -np.inf], np.nan).dropna()
        m = {"col": col, "label": label, "unit": unit, "kind": kind, "group": group}
        if kind == "cat":
            m.update(lo=0, hi=len(LC_GROUPS) - 1, labels=LC_GROUPS)
        elif kind == "log":
            s = np.log10(s.clip(lower=0.1))
            m.update(lo=round(float(s.quantile(.01)), 3), hi=round(float(s.quantile(.99)), 3))
        elif kind == "div":
            a = float(max(abs(s.quantile(.01)), abs(s.quantile(.99))))
            m.update(lo=round(-a, 3), hi=round(a, 3))
        else:
            lo, hi = float(s.quantile(.01)), float(s.quantile(.99))
            if hi <= lo:
                hi = lo + 1.0
            m.update(lo=round(lo, 4), hi=round(hi, 4))
        meta.append(m)
    return meta


def _quantise(vals, m):
    v = vals.astype("float64")
    bad = ~np.isfinite(v)
    if m["kind"] == "cat":
        q = np.clip(np.rint(v), m["lo"], m["hi"])
    else:
        if m["kind"] == "log":
            v = np.log10(np.clip(v, 0.1, None))
        q = np.rint(np.clip((v - m["lo"]) / (m["hi"] - m["lo"]), 0, 1) * 254)
    q[bad] = NODATA
    return q.astype(np.uint8)


def _grid_cells(g):
    c = g["c"]
    i = np.cumsum(np.asarray(g["di"], np.int64)) + int(round(g["ox"] / c))
    j = np.cumsum(np.asarray(g["dj"], np.int64)) + int(round(g["oy"] / c))
    t = np.cumsum(np.asarray(g["dt"], np.int64))
    return i, j, t


def _px_keys(px, c):
    mx, my = _TF.transform(px["x"].to_numpy(float), px["y"].to_numpy(float))
    i = np.floor(np.asarray(mx) / c).astype(np.int64)
    j = np.floor(np.asarray(my) / c).astype(np.int64)
    return i, j


def _pick(order_key, order_t, q_key, q_t):
    """Index into the sorted pixel arrays for each query cell: exact (key, t) if q_t is
    given and present, else the first (earliest) pixel of the key, else -1."""
    first = np.searchsorted(order_key, q_key, side="left")
    ok = (first < len(order_key)) & (order_key[np.minimum(first, len(order_key) - 1)] == q_key)
    idx = np.where(ok, first, -1)
    if q_t is not None and ok.any():
        last = np.searchsorted(order_key, q_key, side="right")
        # exact-time refinement inside the key's run (runs are short: <= a few pixels)
        for off in range(0, 8):
            cand = first + off
            m = ok & (cand < last) & (order_t[np.minimum(cand, len(order_key) - 1)] == q_t)
            idx = np.where(m & (idx == first), cand, idx)
    return idx


def block_for(grid, px, vals_q, onset_ms=None):
    """uint8 [V, n] attribute block for one grid, from this fire's parquet pixels."""
    n = len(grid["di"])
    out = np.full((len(VARS), n), NODATA, np.uint8)
    if px is None or not len(px):
        return out, 0
    i, j, t = _grid_cells(grid)
    pi, pj = _px_keys(px, grid["c"])
    pkey = pi * (1 << 32) + pj
    ptime = px["_min"].to_numpy(np.int64)
    # window grids are relative to the pyroCb onset, whole-fire grids to the fire's own
    # first detection; only the window case can match a time exactly.
    ptime_rel = (ptime - onset_ms) if onset_ms is not None else ptime
    o = np.lexsort((ptime, pkey))
    skey, stime = pkey[o], ptime_rel[o]
    qkey = i * (1 << 32) + j
    idx = _pick(skey, stime, qkey, t if onset_ms is not None else None)
    hit = idx >= 0
    rows = o[idx[hit]]
    out[:, hit] = vals_q[:, rows]
    return out, int(hit.sum())


def build_year(y, data, meta):
    t0 = time.time()
    cols = ["mergid", "x", "y", "datetime_utc", "pcb_id"] + [v[0] for v in VARS]
    px = pq.read_table(PARQUET, columns=cols, filters=[("year", "=", y)]).to_pandas()
    px["_min"] = (px["datetime_utc"].astype("datetime64[ns]").astype("int64") // 60_000_000_000)
    vals_q = np.stack([_quantise(px[m["col"]].to_numpy(), m) for m in meta])
    px = px.reset_index(drop=True)
    by_fire = px.groupby("mergid").indices

    d = data[str(y)]
    pcbp = {f["properties"]["pcb_id"]: f["properties"]
            for f in unpack(d["pcb"])["features"]}

    def sub(ix):
        s = px.iloc[ix]
        return s, vals_q[:, ix]

    files = {"pcb": ({"sim": {}, "pix": {}}, []), "all": ({"sim": {}}, [])}
    off = {"pcb": 0, "all": 0}
    hit = {"pcb": [0, 0], "all": [0, 0]}

    def add(kind, section, key, grid, ix, onset_ms=None):
        if ix is None:
            blk = np.full((len(VARS), len(grid["di"])), NODATA, np.uint8)
            h = 0
        else:
            s, vq = sub(ix)
            s = s.reset_index(drop=True)
            blk, h = block_for(grid, s, vq, onset_ms)
        idx_json, chunks = files[kind]
        idx_json[section][key] = [off[kind], blk.shape[1]]
        chunks.append(np.ascontiguousarray(blk).tobytes())
        off[kind] += blk.shape[1]
        hit[kind][0] += h
        hit[kind][1] += blk.shape[1]

    for k, v in unpack(d["sim"]).items():
        add("pcb", "sim", k, v["g"], by_fire.get(int(k)))
    for k, g in unpack(d["pix"]).items():
        # pcb_id-tagged rows are only a subset of the +/-24 h window's pixels, so draw
        # from every pixel of the fire that made the pyroCb instead.
        pp = pcbp.get(k)
        ms = pd.Timestamp(pp["datetime"]).value // 60_000_000_000 if pp else None
        add("pcb", "pix", k, g, by_fire.get(int(pp["mergid"])) if pp and pp.get("mergid") is not None else None, ms)
    for k, v in unpack(d["simAll"]).items():
        add("all", "sim", k, v["g"], by_fire.get(int(k)))

    os.makedirs(OUT, exist_ok=True)
    msg = []
    for kind, (idx_json, chunks) in files.items():
        raw = b"".join(chunks)
        with open(f"{OUT}/{y}_{kind}.bin.gz", "wb") as f:
            f.write(gzip.compress(raw, 9))
        with open(f"{OUT}/{y}_{kind}.json", "w") as f:
            json.dump(idx_json, f, separators=(",", ":"))
        h, n = hit[kind]
        msg.append(f"{kind}: {n:>8} cells, {100*h/max(n,1):5.1f}% matched, "
                   f"{os.path.getsize(f'{OUT}/{y}_{kind}.bin.gz')/1e6:5.1f} MB gz")
    print(f"  {y}  " + " | ".join(msg) + f"  [{time.time()-t0:.0f}s]", flush=True)


def _fire_means(y):
    """DataFrame: one row per mergid, one column per variable (mean, or mode for cat)."""
    cols = ["mergid"] + [v[0] for v in VARS]
    df = pq.read_table(PARQUET, columns=cols, filters=[("year", "=", y)]).to_pandas()
    for c in cols[1:]:
        df[c] = df[c].astype("float64").replace([np.inf, -np.inf], np.nan)
    g = df.groupby("mergid")
    out = g[[c for c in cols[1:]]].mean()
    for col, _l, _u, kind, _g in VARS:
        if kind == "cat":
            out[col] = df.dropna(subset=[col]).groupby("mergid")[col].agg(lambda s: s.mode().iat[0])
    return out


def build_fires(years, meta):
    """Fire-level aggregates -> pixattr/fires_{year}.json, and flo/fhi into meta."""
    means = {y: _fire_means(y) for y in years}
    allm = pd.concat(means.values())
    for m in meta:
        s = allm[m["col"]].dropna()
        if m["kind"] == "cat":
            m["flo"], m["fhi"] = m["lo"], m["hi"]
            continue
        if m["kind"] == "log":
            s = np.log10(s.clip(lower=0.1))
        lo, hi = float(s.quantile(.02)), float(s.quantile(.98))
        if m["kind"] == "div":
            lo, hi = -max(abs(lo), abs(hi)), max(abs(lo), abs(hi))
        if hi <= lo:
            hi = lo + 1.0
        m["flo"], m["fhi"] = round(lo, 4), round(hi, 4)
    for y, df in means.items():
        cols = []
        for m in meta:
            fm = dict(m, lo=m["flo"], hi=m["fhi"])
            cols.append(_quantise(df[m["col"]].to_numpy(), fm))
        q = np.stack(cols, axis=1)
        with open(f"{OUT}/fires_{y}.json", "w") as f:
            json.dump({str(int(k)): r.tolist() for k, r in zip(df.index, q)}, f, separators=(",", ":"))
        print(f"  fires_{y}.json: {len(df)} fires", flush=True)


def main():
    args = sys.argv[1:]
    fires_only = "--fires-only" in args
    data = load_map()
    years = [int(a) for a in args if a.isdigit()] or sorted(int(y) for y in data)
    os.makedirs(OUT, exist_ok=True)
    if fires_only:
        meta = json.load(open(f"{OUT}/meta.json"))["vars"]
    else:
        meta = _fit_scales(years)
    if not fires_only:
        for y in years:
            build_year(y, data, meta)
    build_fires(years, meta)
    with open(f"{OUT}/meta.json", "w") as f:
        json.dump({"vars": meta, "nodata": NODATA}, f, separators=(",", ":"))


if __name__ == "__main__":
    main()
