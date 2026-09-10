"""Generate 'PCB_MATCHING_CLEAN.ipynb'.

A separate, minimal notebook with one job: match every pyroCb event to the FEDS
Snapshot fire (`mergid`) that was actively burning near it, using the radius-ladder
algorithm below, and export the matched events with their active-burn geometry and
every Snapshot attribute (fireid..mergedate, no lcc_final). Unlike 'Data exploration
.ipynb' (D-20: analysis lives in the package, notebook narrates), this file is
deliberately self-contained -- the matching logic is written out in the cells
themselves, not hidden behind a function call, so it can be read and checked
line-by-line. It only reuses long-established, already-tested primitives from
`pyrocb` (CRS constants, step parsing, geometry repair) -- never rounds a
coordinate or a timestamp, and documents every column's aggregation rule inline.

Run:  python make_pcb_matching_clean_notebook.py
Then open and run PCB_MATCHING_CLEAN.ipynb in Jupyter.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
def md(t):   cells.append(nbf.v4.new_markdown_cell(t.strip("\n")))
def code(t): cells.append(nbf.v4.new_code_cell(t.strip("\n")))


md(r"""
# PyroCb <-> FEDS Snapshot matching (clean, minimal, full precision)

One job: for every pyroCb event, find the nearest **actively burning** fire in the
FEDS Snapshot atlas around the time of onset, using the radius-ladder algorithm
below, and export that event with the fire's active-burn geometry and every
Snapshot attribute available for it (`fireid` .. `mergedate`, everything except
`lcc_final` -- which is not a Snapshot column at all, see the note in step 6).

No coordinate or timestamp is rounded anywhere in this file. The PCB source
coordinates themselves are stored to 0.1 deg (~11 km) in `PCBAtlas.xlsx` -- that is
a property of the input data, not something this notebook does; everything computed
here keeps full `float64`/`datetime64[ns]` precision.

**Standalone.** This notebook imports no first-party module -- not even `pyrocb`.
Everything it needs from that package (CRS constants, the half-day step schedule,
`read_pcb`, geometry repair) is copied inline in step 1b, verbatim, so the whole
matching pipeline can be read, run and audited from this one file with nothing to
jump to. Only external packages (`numpy`/`pandas`/`geopandas`/`shapely`) are imported.

## The matching algorithm (as specified)

1. Start from the pyroCb's own **coordinates and onset timestamp**.
2. Within **10 km**, look for the **closest snapshot in time** (before OR after
   onset, whichever is nearer). If any fire (`mergid`) is within 10 km in that
   snapshot -> **matched**.
3. If not: look at the **second-closest** snapshot in time, but **only if** it is
   within **+/-12 h** of onset. If it has a fire within 10 km -> **matched**.
4. Still nothing? Step the radius up: **12.5, 15, 17.5, 20 km**, repeating 2-3 at
   each radius (same two candidate snapshots -- they depend only on time, not on
   the search radius, so they don't need to be re-picked).
5. Nothing within 20 km at either candidate snapshot -> the event is **unmatched**
   (`N/A`) -- no active burning area could be associated with it.
6. For every match, store the **delta** between onset and the matched snapshot's
   nominal time, and the **radius** at which the match was found, alongside every
   Snapshot attribute for the matched fire.

### What "the active burning area for that PCB" means here

The buffer is not just a hit test -- it defines the area. Once the ladder settles on
a radius R, the geometry stored for the event is **every fire polygon in the matched
snapshot that reaches into the R km buffer, unioned and clipped to that buffer**, not
the whole fire.

That distinction is the difference between a usable answer and a nonsense one on big
fires. A FEDS snapshot polygon is **cumulative season-to-date burned area**, not the
actively-burning front (this project measures that directly: max-snapshot `n_pixels`
over final is exactly 1.0 for every fire; the burning front is carried separately by
`flinelen`/`isactive`/`t_inactive`). So a large complex's footprint here spans up to
~300 km, and its representative point can land **150 km** from the pyroCb that
matched it -- a different lobe of the fire, weeks of spread away. Clipping to the
buffer the ladder actually searched keeps the geometry, and any position derived from
it, where the event is.

The whole-fire footprint is not discarded: it ships as a second GeoPackage layer
(`matched_fire_full`), so the two can be compared directly.

**Consequence, stated plainly:** the `snap_*` attribute columns describe the **whole
matched fire** (that is the only scale FEDS publishes them at), while `geometry` and
`clip_area_km2` describe **only the burn inside the buffer**. Do not divide one by
the other and call it a density.

Two more things worth being explicit about before any code runs:

* **Only the second candidate carries the +/-12 h cap.** The first (closest-in-time)
  candidate is tried unconditionally, whatever its lag -- it's "the closest we've
  got", and its actual lag is stored (`cand1_dt_h`) so a large one is visible rather
  than silently accepted. The cap is checked with `abs(dt_h) <= 12.0`, applied to
  the candidate's **nominal** overpass time (see step 7 below) -- verified by
  assertion in the checks cell, since the prompt asked for this to be double-checked.
* **"Closest snapshot" is picked globally by time**, not restricted to the fixed
  AM/PM pair around onset -- it's whichever existing Snapshot file (this year, or
  the adjacent year near a Dec/Jan boundary) has the nearest nominal overpass time,
  so a real data gap is handled correctly instead of assuming the schedule is dense.
""")

md("## 1. Setup")
code(r"""
# External packages only -- no first-party import. See step 1b for the small set
# of pyrocb primitives (CRS constants, step parsing, geometry repair) copied
# inline so this notebook has nothing else to depend on.
import warnings; warnings.filterwarnings("ignore")

import datetime as dt
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely

pd.set_option("display.width", 160)
pd.set_option("display.max_columns", 60)
""")

md(r"""
## 1b. Inlined `pyrocb` primitives

Copied verbatim from `pyrocb/config.py`, `pyrocb/timebins.py` and `pyrocb/io.py` --
not reimplemented or approximated, so this notebook's numbers stay identical to the
rest of the project. Kept to exactly the handful this algorithm needs: the two CRS
codes, the Snapshot file-path template, the half-day step schedule (parsing a step
name, its nominal UTC instant, and the local-solar AM/PM convention), geometry
repair, and the PCB loader.
""")
code(r"""
# ---- paths + CRS (pyrocb/config.py) ----------------------------------------
# Assumes this notebook runs from the repo root (same directory as PCBAtlas.xlsx
# and the per-year data folders) -- the same convention "Data exploration.ipynb"
# already relies on for its own relative output paths.
ROOT = Path.cwd()
SNAP_DIR = str(ROOT / "{y}" / "{y}" / "Snapshot")
NFP_DIR = str(ROOT / "{y}" / "{y}" / "NFP")
XLSX_FP = ROOT / "NRL_pyroCb_inventory_2013-2024_updated_20260319_with_auto_count.xlsx"
METRIC_CRS = 3571   # WGS 84 / North Pole LAEA Bering Sea - equal-area, true km distances
WGS84 = 4326
PRECISION_M = 1.0   # snap grid before every boolean; 1/375 of a VIIRS pixel

# ---- half-day step schedule (pyrocb/timebins.py) ---------------------------
# FEDS bins observations into two per-day VIIRS overpass groups, AM/PM, named by
# the calendar date + which half. NOMINAL_LOCAL_H is only the nominal local-solar
# hour of each group's overpass, used to convert a step name to an approximate
# UTC instant anywhere on Earth (via longitude).
STEP_RE = re.compile(r"^(\d{4})(\d{2})(\d{2})(AM|PM)$")
NOMINAL_LOCAL_H = {"AM": 2.0, "PM": 13.0}

def parse_step(step):
    m = STEP_RE.match(step)
    if not m:
        raise ValueError(f"bad step {step!r}")
    y, mo, d, ap = m.groups()
    return dt.date(int(y), int(mo), int(d)), ap

def step_sort_key(step):
    d, ap = parse_step(step)
    return (d, 0 if ap == "AM" else 1)

def make_step(d, ampm):
    return f"{d:%Y%m%d}{ampm}"

def norm_lon(lon):
    '''Normalise longitude to (-180, 180].'''
    lon = np.asarray(lon, dtype=float)
    out = ((lon + 180.0) % 360.0) - 180.0
    return np.where(out == -180.0, 180.0, out)

def nominal_utc(step, lon):
    '''Approximate UTC instant of a step's overpass at a given longitude.'''
    d, ap = parse_step(step)
    return (pd.Timestamp(d) + pd.Timedelta(hours=NOMINAL_LOCAL_H[ap])
            - pd.Timedelta(hours=float(norm_lon(lon)) / 15.0))

# ---- geometry repair (pyrocb/io.py) ----------------------------------------
# Never `buffer(0)`: on a bowtie it silently deletes the smaller loop and on some
# MultiPolygons drops whole lobes, i.e. it changes area without saying so - fatal
# when the output *is* an area. `make_valid` preserves area.
_POLY = {"Polygon", "MultiPolygon"}

def _polygonal(g):
    '''Keep only polygonal parts; make_valid can emit GeometryCollections.'''
    if g is None or g.is_empty:
        return shapely.Polygon()
    if g.geom_type in _POLY:
        return g
    parts = [p for p in shapely.get_parts(g) if p.geom_type in _POLY]
    if not parts:
        return shapely.Polygon()
    return parts[0] if len(parts) == 1 else shapely.union_all(parts)

def ensure_valid(geoms):
    '''Repair invalid / non-polygonal geometry. Cheap: no-ops when input is clean.'''
    g = np.asarray(geoms, dtype=object)
    bad = ~shapely.is_valid(g)
    if bad.any():
        g = g.copy()
        g[bad] = shapely.make_valid(g[bad])
    types = shapely.get_type_id(g)          # 3 = Polygon, 6 = MultiPolygon
    odd = ~np.isin(types, (3, 6))
    if odd.any():
        g = g.copy()
        g[odd] = [_polygonal(x) for x in g[odd]]
    return g

# ---- Snapshot step listing (pyrocb/io.py) ----------------------------------
def list_steps(year):
    '''Sorted snapshot step names for a year, excluding the *_FL fire-line files.'''
    d = SNAP_DIR.format(y=year)
    names = [os.path.basename(p)[:-5] for p in glob.glob(os.path.join(d, "*.gpkg"))]
    names = [n for n in names if not n.endswith("_FL")]
    return tuple(sorted(names, key=step_sort_key))

# ---- PCB loader (pyrocb/io.py) ---------------------------------------------
def _make_pcb_ids(d):
    '''PCB_YYYY_MMDD_HHMM, with a deterministic a/b/c suffix where a timestamp repeats.'''
    base = "PCB_" + d["datetime"].dt.strftime("%Y_%m%d_%H%M")
    out = base.copy()
    for key, idx in base.groupby(base).groups.items():
        if len(idx) == 1:
            continue
        order = d.loc[idx].sort_values(["lat", "lon"]).index
        for i, ix in enumerate(order):
            out.loc[ix] = f"{key}{chr(ord('a') + i)}"
    return out

def read_pcb():
    '''The NRL pyroCb inventory, cleaned, with a stable `pcb_id`. No caching layer
    here (unlike pyrocb.io.read_pcb) - 912 rows from one Excel file is cheap enough
    to just recompute every run, keeping this notebook dependency-free.

    The 2013-2024 refresh is global (`country`, not the old boreal-only `region`)
    and its sheet ships an auto-generated count/QC block appended below the real
    rows ("TOTAL"/"COUNTRY"/"PYROCB COUNT" etc., every field but fire_name blank or
    non-numeric). Renaming keeps `region` the name used everywhere below. lat/lon
    need an explicit numeric coercion first - the summary rows leave the column
    object-dtyped (mixed floats and strings like "PROBABLE"), and a bare dropna
    only catches an actual NaN, not a string; coercing first turns every such row
    into one, and the existing dropna on datetime/lat/lon removes it (checked: it
    already does, since every summary row's datetime is also unparseable).'''
    d = pd.read_excel(XLSX_FP, sheet_name=0)
    if "country" in d.columns and "region" not in d.columns:
        d = d.rename(columns={"country": "region"})
    d["datetime"] = pd.to_datetime(d["datetime"], errors="coerce")
    d["lat"] = pd.to_numeric(d["lat"], errors="coerce")
    d["lon"] = pd.to_numeric(d["lon"], errors="coerce")
    d = d.dropna(subset=["datetime", "lat", "lon"]).copy()
    d["year"] = d["datetime"].dt.year
    d["doy"] = d["datetime"].dt.dayofyear
    d = d.sort_values(["datetime", "lat", "lon"]).reset_index(drop=True)
    for c in d.columns:
        if d[c].dtype == object:
            d[c] = d[c].astype("string")
    d["pcb_id"] = _make_pcb_ids(d)
    return d
""")

md(r"""
## 2. Algorithm constants

`RADII_KM` is the exact ladder from the spec. `TIME_CAP_H` is the +/-12 h cap on
the *second* candidate only (see the note above).
""")
code(r"""
RADII_KM = [10.0, 12.5, 15.0, 17.5, 20.0]   # search radii, in the order tried
TIME_CAP_H = 12.0                           # cap on the 2nd (not 1st) candidate
ANALYSIS_RADIUS_KM = 5.0                    # stage 3: window around the corrected centre
ANALYSIS_TIME_CAP_H = 12.0                  # stage 3: pixels must ALSO be within +/-12h of onset
V4_N_STEPS = 4                              # v4: how many time-closest snapshots to pool
V4_TAU_H = 12.0                             # v4: time-decay constant, one overpass interval
""")

md(r"""
## 3. Load the pyroCb events

`read_pcb()` (step 1b) is a verbatim copy of the loader `pyrocb.matching` uses.
Nothing is rounded or reprojected here: `datetime` stays `datetime64[ns]` (naive
UTC, as supplied), `lat`/`lon` stay the raw `float64` degrees from `PCBAtlas.xlsx`.
All 18 source columns are kept so the exported table still reads as
pyroCb-events-first.
""")
code(r"""
pcb = read_pcb()
print(f"{len(pcb)} pyroCb events, {pcb['year'].min()}-{pcb['year'].max()}")
pcb[["pcb_id", "datetime", "lat", "lon"]].head()
""")

md(r"""
## 4. Snapshot schema: what "fireid .. mergedate" means, and how each field is handled

A raw Snapshot `.gpkg` (`{step}.gpkg`, *not* the reduced `{step}_FL.gpkg` variant)
has exactly 32 attribute columns, `fireid` first and `mergedate` last -- confirmed
directly against the file schema, not assumed. `lcc_final` is not among them; it
lives only in the season-final-perimeter product, a different file this notebook
never reads, so "no lcc_final" is automatically satisfied rather than something to
filter out.

A single `mergid` (the stable fire id) can be spread over many `fireid` rows in one
snapshot -- up to dozens, for a fire built from several merged components that
haven't been geometrically dissolved together yet. Reporting "the" Snapshot
attributes for a matched fire therefore means **combining its component rows**, and
how each column combines depends on what kind of quantity it is. This is the
complete, column-by-column policy -- applied explicitly in step 9, nothing left
implicit:

| raw column(s) | kind | rule | output column(s) |
|---|---|---|---|
| `fireid` | per-component id, not a fire-level value | list every contributing id | `snap_fireids`, `snap_n_components` |
| `mergid` | stable fire id | unchanged (join key) | `mergid` |
| `ftype` | categorical, **varies per component** (checked: e.g. `forest-peat` + `Crop` + `grassland:peat` in one real merge group) | area-weighted mode (component with the largest `farea`) **and** the full unique set, so nothing is silently dropped | `snap_ftype_dominant`, `snap_ftype_all` |
| `isactive` | 0/1 flag | max (any component active). Snapshots are documented as "active fires only" -- verified empirically (see step 8) that it is 1 for every row seen | `snap_isactive` |
| `t_inactive` | hours since last active step, per component | min (freshest component) | `snap_t_inactive` |
| `isignition` | 0/1 flag | max. **Deviation from `pyrocb.io.read_snapshot`**, which sums this by default -- summing a boolean flag across components conflates "1 ignition-stage component" with "3", so max is used here instead | `snap_isignition` |
| `invalid` | 0/1 flag | max, defensively. Empirically always 0 (checked over 10,014 rows / 60 files) -- kept rather than dropped in case a future file has one | `snap_invalid` |
| `n_pixels`, `n_newpixels` | extensive counts | sum | `snap_n_pixels`, `snap_n_newpixels` |
| `farea`, `fperim`, `flinelen` | extensive (km2 / km) | sum -- components are geometrically disjoint (not yet merged), so summed perimeter/line-length is exact, not an approximation | `snap_farea_km2`, `snap_fperim_km`, `snap_flinelen_km` |
| `duration` | days burning, per component | max (oldest/longest-running component) | `snap_duration_days` |
| `pixden` | ratio (pixels / area) | **recomputed** as `sum(n_pixels) / sum(farea)` rather than averaged. **Deviation from `read_snapshot`**, which means the per-component ratios -- biased when components differ in size; recomputing from the summed numerator/denominator is the statistically correct way to combine a density | `snap_pixden` |
| `meanFRP` (W/m2) | intensive (area-normalised) | area-weighted mean, `sum(meanFRP_i * farea_i) / sum(farea_i)`. **Deviation from both `read_snapshot` and `pyrocb.spread`**, which plain-average it -- a small hot component would otherwise pull the mean as hard as a large one | `snap_meanFRP_Wm2` |
| `FRP95` (MW) | 95th-percentile statistic, per component | max of the per-component P95s. The true complex-level P95 needs pixel-level data this file doesn't carry, so max is a conservative proxy (matches existing convention) | `snap_FRP95_MW` |
| `avgspread` (m) | intensive (per-step advance), **empty in ~25% of raw rows** | area-weighted mean **over the non-null components only** -- a plain `np.average` propagates one missing component's NaN across the whole fire, which emptied this column for 79-100% of multi-component fires. Component count that contributed is reported alongside | `snap_avgspread_m`, `snap_n_comp_with_spread` |
| `spread95` (m) | 95th-percentile statistic | max, same reasoning as `FRP95` | `snap_spread95_m` |
| `tst_year/month/day/ampm` | per-component start half-day | converted to a real UTC instant via `timebins.nominal_utc` at *that component's own* representative-point longitude, then **min** across components (earliest-starting component = the complex's ignition) | `snap_tst_utc` |
| `ted_year/month/day/ampm` | per-component last-active half-day | same conversion, **max** across components (most recently active) | `snap_ted_utc` |
| `ted_doy` | day-of-year of `ted` | **not** aggregated from the raw column -- recomputed fresh from `snap_ted_utc` so the two stay consistent by construction | `snap_ted_doy` |
| `lake_area`, `lake_border`, `lake_border_tot` | extensive (km2 / km) | sum | `snap_lake_area_km2`, `snap_lake_border_km`, `snap_lake_border_tot_km` |
| `lake_no` | count | sum (same disjoint-components convention as `farea`; a lake straddling two components would be double-counted, a known, documented approximation shared with `read_snapshot`) | `snap_lake_no` |
| `mergedate` | date this *component* was merged **into** another fireid, **varies per component**, and is **null** for a component that never merged into anything (checked directly: every size-1 group is 100% null; a merged group typically has exactly one null member -- the surviving/absorbing fireid, which nothing merges *it* into) | max of the non-null values (latest merge = when the currently-visible complex became complete; `NaT` if the fire never merged at all), plus the full non-null unique set | `snap_mergedate_latest`, `snap_mergedate_all` |

`avgspread`/`spread95` ship as strings in the raw file (checked) and are coerced to
numeric before any of this runs.
""")

md(r"""
## 5. Which years have Snapshot data

Every step file that actually exists on disk, per year -- used both to build the
candidate-snapshot search pool and to read a specific matched step later. This is
`list_steps` (step 1b), unmodified from `io.list_steps`: it globs `{step}.gpkg` and
explicitly excludes the reduced `{step}_FL.gpkg` variant (checked against its
source -- an earlier pass in this investigation accidentally read an `_FL` file
and got a misleadingly tiny 2-column schema, which is exactly the mistake this
exclusion prevents).
""")
code(r"""
# Cache per year so repeated candidate searches don't re-glob the directory.
_steps_by_year = {}
def steps_for_year(y):
    if y not in _steps_by_year:
        _steps_by_year[y] = list_steps(y)   # () if the year has no Snapshot dir
    return _steps_by_year[y]

for y in range(pcb["year"].min(), pcb["year"].max() + 1):
    print(y, len(steps_for_year(y)), "steps")
""")

md(r"""
## 6. Step 1-2: the two time-closest candidate snapshots

For a given onset time and longitude, rank every step from the onset's own year plus
its immediate neighbours (covers a pyroCb landing right at a Dec/Jan boundary) by
how close its **nominal** overpass time is to onset, and keep the closest two.
`nominal_utc` (step 1b) already encodes the local-solar AM/PM schedule -- reused
as-is rather than re-derived. Ties (equal `abs(dt_h)`) prefer the earlier step, a
one-line, deterministic tie-break that is not expected to bind in practice.

`dt_h` is signed: negative = snapshot before onset, positive = after -- this *is*
the delta the spec asks to store, before any match is even attempted.
""")
code(r"""
def closest_two_steps(onset_utc, lon, year, n=2):
    '''[(step, dt_h), ...], length 0-n, nearest-in-time first.

    n=2 is the ladder (steps 2-3). v4 calls the same function with n=V4_N_STEPS, so
    both use one definition of "closest in time" rather than two that could drift.'''
    pool = [s for y in (year - 1, year, year + 1) for s in steps_for_year(y)]
    if not pool:
        return []
    dt_h = [(nominal_utc(s, lon) - onset_utc).total_seconds() / 3600.0 for s in pool]
    order = sorted(range(len(pool)), key=lambda i: (abs(dt_h[i]), dt_h[i]))
    return [(pool[i], dt_h[i]) for i in order[:n]]

# Smoke-test on the first event.
r0 = pcb.iloc[0]
closest_two_steps(r0["datetime"], r0["lon"], int(r0["year"]))
""")

md(r"""
## 7. Step 2-4: "is there fire within R km" at a candidate snapshot

Two small, cached readers:

* `load_raw(year, step)` -- the untouched per-component rows of one Snapshot file,
  reprojected to `METRIC_CRS` (EPSG:3571, equal-area -- distances below are true
  km) and geometry-repaired with `ensure_valid` (step 1b; never `buffer(0)` -
  that silently changes area, which matters when this notebook reports an
  active-burn *geometry*). Cached because several events can share a step.
* `dissolved(year, step)` -- the same file collapsed to one **union** geometry per
  `mergid` (components of an unfinished merge are physically the same fire, so the
  nearest-fire search should see them as one shape), with a `shapely.STRtree` for
  fast nearest-neighbour queries. Geometry-only at this stage; the full attribute
  aggregation (step 9's table) happens later, and only for the one fire that
  actually wins a match -- no need to aggregate every fire in every candidate file.
""")
code(r"""
_raw_cache, _dis_cache = {}, {}

def load_raw(year, step):
    key = (year, step)
    if key not in _raw_cache:
        fp = os.path.join(SNAP_DIR.format(y=year), f"{step}.gpkg")
        g = gpd.read_file(fp)
        if g.crs is None or g.crs.to_epsg() != METRIC_CRS:
            g = g.to_crs(METRIC_CRS)
        for c in ("avgspread", "spread95"):     # ship as strings in the raw file
            g[c] = pd.to_numeric(g[c], errors="coerce")
        g["geometry"] = ensure_valid(g.geometry.values)
        _raw_cache[key] = g
    return _raw_cache[key]

def dissolved(year, step):
    key = (year, step)
    if key not in _dis_cache:
        g = load_raw(year, step)
        d = g[["mergid", "geometry"]].dissolve(by="mergid").reset_index()
        tree = shapely.STRtree(d.geometry.values)
        _dis_cache[key] = (d, tree)
    return _dis_cache[key]

def nearest_mergid(pt_metric, year, step):
    '''(mergid, distance_km) of the fire nearest `pt_metric` in this snapshot.'''
    d, tree = dissolved(year, step)
    if not len(d):
        return None, np.inf
    k = tree.nearest(pt_metric)
    dist_km = float(shapely.distance(pt_metric, d.geometry.values[k])) / 1000.0
    return int(d["mergid"].values[k]), dist_km
""")

md(r"""
## 8. Run the radius ladder over every event

PCB points are projected to `METRIC_CRS` once, up front, at full precision (no
rounding of the projected coordinates either). For each event: get its two
candidates (step 6), find the nearest fire at each (step 7) **once** -- the
distances don't change as the radius steps up, only which radius first admits
them does -- then walk the ladder exactly as specified: at each radius, candidate 1
first, candidate 2 only if it passes the +/-12 h cap.
""")
code(r"""
pts_metric = (gpd.GeoSeries(shapely.points(pcb["lon"].values, pcb["lat"].values), crs=WGS84)
              .to_crs(METRIC_CRS))

rows = []
for i, (_, r) in enumerate(pcb.iterrows()):
    onset, lon, year, pt = r["datetime"], r["lon"], int(r["year"]), pts_metric.values[i]
    cands = closest_two_steps(onset, lon, year)
    c1 = cands[0] if len(cands) > 0 else None
    c2 = cands[1] if len(cands) > 1 else None
    m1 = nearest_mergid(pt, year, c1[0]) if c1 else (None, np.inf)
    m2 = nearest_mergid(pt, year, c2[0]) if c2 else (None, np.inf)
    cand2_eligible = bool(c2 and abs(c2[1]) <= TIME_CAP_H)

    match = None
    for radius in RADII_KM:
        if c1 and m1[1] <= radius:                       # step 2: closest snapshot
            match = dict(step=c1[0], dt_h=c1[1], radius=radius, dist_km=m1[1],
                         mergid=m1[0], rank=1)
            break
        if cand2_eligible and m2[1] <= radius:            # step 3: 2nd, capped at 12h
            match = dict(step=c2[0], dt_h=c2[1], radius=radius, dist_km=m2[1],
                         mergid=m2[0], rank=2)
            break

    rows.append(dict(
        pcb_id=r["pcb_id"],
        match_status="matched" if match else "unmatched",
        cand1_step=c1[0] if c1 else None, cand1_dt_h=c1[1] if c1 else np.nan,
        cand1_dist_km=m1[1],
        cand2_step=c2[0] if c2 else None, cand2_dt_h=c2[1] if c2 else np.nan,
        cand2_dist_km=m2[1], cand2_eligible_12h=cand2_eligible,
        match_step=match["step"] if match else None,
        match_dt_h=match["dt_h"] if match else np.nan,          # the requested delta
        search_radius_km=match["radius"] if match else np.nan,  # the requested radius
        match_rank=match["rank"] if match else pd.NA,
        # The ladder's own answer: the NEAREST fire, which is what decided whether
        # this event matched at all. Which fire gets CREDITED (and supplies every
        # snap_* column) is decided separately, on burn inside the buffer - see
        # local_active_burn(). Both are kept so the two rules stay auditable.
        ladder_mergid=match["mergid"] if match else pd.NA,
        ladder_dist_km=match["dist_km"] if match else np.nan,
    ))

match_tbl = pd.DataFrame(rows)
print(match_tbl["match_status"].value_counts())
match_tbl.head()
""")

md(r"""
## 9. Aggregate the matched fire's Snapshot attributes (the table from step 4)

For every **matched** event, pull that fire's raw component rows out of the already
-cached snapshot file (`load_raw`, no re-read) and combine them into the single
`snap_*` row using exactly the rule assigned to each column above -- one clearly
labelled block per rule, nothing bundled into a generic aggfunc dict so each choice
stays visible and checkable.
""")
code(r"""
def aggregate_matched_fire(year, step, mergid):
    '''One dict of every snap_* column (step 4's table) for this fire at this step.'''
    g = load_raw(year, step)
    c = g[g["mergid"] == mergid]
    farea = c["farea"].to_numpy(dtype=float)
    out = {}

    # -- identity: list every contributing component, don't collapse to one --
    out["snap_fireids"] = ",".join(str(v) for v in sorted(c["fireid"]))
    out["snap_n_components"] = len(c)

    # -- categorical, varies per component: area-weighted dominant + full set --
    out["snap_ftype_dominant"] = c.loc[c["farea"].idxmax(), "ftype"]
    out["snap_ftype_all"] = ";".join(sorted(c["ftype"].unique()))

    # -- 0/1 flags: max --
    out["snap_isactive"] = int(c["isactive"].max())
    out["snap_isignition"] = int(c["isignition"].max())
    out["snap_invalid"] = int(c["invalid"].max())

    # -- extensive counts / areas / lengths: sum (components are disjoint) --
    for src, dst in [("n_pixels", "snap_n_pixels"), ("n_newpixels", "snap_n_newpixels"),
                      ("farea", "snap_farea_km2"), ("fperim", "snap_fperim_km"),
                      ("flinelen", "snap_flinelen_km"), ("lake_area", "snap_lake_area_km2"),
                      ("lake_border", "snap_lake_border_km"),
                      ("lake_border_tot", "snap_lake_border_tot_km"),
                      ("lake_no", "snap_lake_no")]:
        out[dst] = float(c[src].sum())

    # -- per-component extremes: min / max --
    out["snap_t_inactive"] = float(c["t_inactive"].min())
    out["snap_duration_days"] = float(c["duration"].max())

    # -- ratio, recomputed from summed parts rather than averaged --
    out["snap_pixden"] = out["snap_n_pixels"] / out["snap_farea_km2"] if out["snap_farea_km2"] else np.nan

    # -- intensive quantities: area-weighted mean over the components that
    # -- actually HAVE a value. 25% of raw `avgspread` cells are empty in the
    # -- source, and np.average propagates a single NaN over the whole fire - which
    # -- wiped out snap_avgspread_m for 79-100% of MULTI-component fires (i.e.
    # -- precisely the big merge complexes) while sparing 92% of single-component
    # -- ones. Missing components are dropped from numerator and denominator
    # -- instead; all-missing still yields NaN, and the count of components that
    # -- did contribute is reported so a thin mean is visible rather than implied.
    w = farea if farea.sum() > 0 else np.ones_like(farea)

    def wmean(col):
        v = c[col].to_numpy(dtype=float)
        ok = np.isfinite(v) & np.isfinite(w)
        if not ok.any():
            return np.nan
        return float(np.average(v[ok], weights=w[ok]) if w[ok].sum() > 0 else v[ok].mean())

    out["snap_meanFRP_Wm2"] = wmean("meanFRP")
    out["snap_avgspread_m"] = wmean("avgspread")
    out["snap_n_comp_with_spread"] = int(np.isfinite(c["avgspread"].to_numpy(dtype=float)).sum())

    # -- 95th-percentile statistics: max of the per-component values --
    out["snap_FRP95_MW"] = float(c["FRP95"].max())
    out["snap_spread95_m"] = float(c["spread95"].max())

    # -- per-component half-day stamps -> real UTC instants, at each component's --
    # -- own representative-point longitude, then min (tst) / max (ted) --
    def half_day_utc(sub, prefix):
        lon = sub.geometry.representative_point().to_crs(WGS84).x.to_numpy()
        steps = [make_step(pd.Timestamp(int(y), int(m), int(d)).date(), ap)
                 for y, m, d, ap in zip(sub[f"{prefix}_year"], sub[f"{prefix}_month"],
                                        sub[f"{prefix}_day"], sub[f"{prefix}_ampm"])]
        return [nominal_utc(s, lo) for s, lo in zip(steps, lon)]

    out["snap_tst_utc"] = min(half_day_utc(c, "tst"))
    out["snap_ted_utc"] = max(half_day_utc(c, "ted"))
    out["snap_ted_doy"] = int(out["snap_ted_utc"].dayofyear)   # recomputed, not aggregated

    # -- per-component merge date: max (latest) + full set. NaN for any component --
    # -- that never merged INTO another (checked: a lone fireid == mergid is 100% --
    # -- NaN here; in a merged group it's the one surviving/absorbing component -- --
    # -- so `dropna` before sorting, and an all-NaN group (no merge yet) reports --
    # -- None rather than crashing on a mixed str/NaN sort --
    md_dates = pd.to_datetime(c["mergedate"], format="%Y-%m-%d-%p", errors="coerce")
    out["snap_mergedate_latest"] = md_dates.max()  # NaT if the fire never merged
    seen = sorted(c["mergedate"].dropna().unique())
    out["snap_mergedate_all"] = ";".join(seen) if seen else None

    # -- the WHOLE fire's footprint at this step: union of its components, full --
    # -- precision, no simplification. This is the fire, NOT "the burn near this --
    # -- pyroCb" - local_active_burn() below produces that. --
    out["fire_geometry"] = shapely.union_all(c.geometry.values)
    return out


def local_active_burn(year, step, pt_metric, radius_km):
    '''The burn actually around this pyroCb: every fire polygon in this snapshot
    that reaches into the search buffer, unioned and clipped to that buffer.

    This, not the whole fire, is the geometry the event gets. Clipping matters
    enormously on big fires: a FEDS snapshot polygon is CUMULATIVE season-to-date
    burned area, so a large complex spans up to ~300 km here and its
    representative point can sit 150 km from the pyroCb - a different lobe of the
    fire, weeks of spread away. The buffer is what the ladder actually searched,
    so the burn inside it is what the ladder actually found.

    Identity-free, like the search itself: if two fires reach into the buffer both
    are kept (mergids recorded), rather than dropping one because the other
    happened to be a few metres nearer.'''
    g = load_raw(year, step)
    buf = shapely.buffer(pt_metric, radius_km * 1000.0)
    hit = g[shapely.dwithin(pt_metric, g.geometry.values, radius_km * 1000.0)].copy()
    if not len(hit):
        return {"snap_clip_geometry": None, "clip_area_km2": 0.0, "n_fires_in_buffer": 0,
                "mergids_in_buffer": None, "mergid": None, "match_dist_km": np.nan,
                "dominant_clip_km2": np.nan, "fb_x": np.nan, "fb_y": np.nan}

    # Per-FIRE clipped area, so "which fire is this pyroCb's" can be decided on how
    # much of each is actually burning in the buffer.
    per_fire = (hit[["mergid", "geometry"]].dissolve(by="mergid").reset_index())
    per_fire["clip_km2"] = shapely.area(
        shapely.intersection(per_fire.geometry.values, buf)) / 1e6
    per_fire["dist_km"] = shapely.distance(pt_metric, per_fire.geometry.values) / 1000.0

    # The fire the event is credited to is the DOMINANT one in the buffer, not merely
    # the nearest. Nearest-wins let a 0.11 km2 speck 4.0 km away beat a 15.1 km2 fire
    # 6.2 km away (PCB_2023_0630_0620) and a 0.42 km2 sliver beat a 41.7 km2 fire
    # (PCB_2022_0708_2210) - attaching that speck's FRP, area and duration to the
    # pyroCb. Ties break on distance. The LADDER decision is untouched: whether this
    # event matched at all, and at which radius, still comes from the nearest fire
    # (cand*_dist_km), so the search rule and the credit rule stay separable.
    win = per_fire.sort_values(["clip_km2", "dist_km"],
                               ascending=[False, True]).iloc[0]

    clipped = shapely.intersection(shapely.union_all(hit.geometry.values), buf)
    mids = sorted(hit["mergid"].unique())
    pos = burn_position(clipped)
    return {"snap_clip_geometry": clipped,
            "clip_area_km2": float(shapely.area(clipped)) / 1e6,
            "n_fires_in_buffer": len(mids),
            "mergids_in_buffer": ",".join(str(int(v)) for v in mids),
            "mergid": int(win["mergid"]),
            "match_dist_km": float(win["dist_km"]),
            "dominant_clip_km2": float(win["clip_km2"]),
            # Fallback only - the shipped position comes from the pixels (stage 2).
            "fb_x": float(shapely.get_x(pos)), "fb_y": float(shapely.get_y(pos))}


def burn_position(clipped):
    '''Where to plot this pyroCb: the AREA-WEIGHTED centroid of the burn in the
    buffer, snapped back onto the burn if it lands off it.

    `representative_point()` (GEOS PointOnSurface) was used first and is the wrong
    tool: it only guarantees a point INSIDE the polygon and is weighted by nothing
    at all, so a bigger or denser lobe does not pull it. 140 of 395 clipped burns
    here are multi-part (up to 7 lobes) - exactly the case where that matters - and
    it sat a median 0.86 km, worst case 11.9 km, from the weighted centre.

    The plain centroid fixes the weighting but reintroduces D-49's problem: a
    concave or multi-lobed burn can put its centre of mass on unburned ground (7.6%
    of events here). Snapping the centroid to the nearest point on the burn keeps
    the weighting AND guarantees the marker sits on ground that actually burned.'''
    if clipped is None or shapely.is_empty(clipped):
        return None
    cen = shapely.centroid(clipped)
    if shapely.within(cen, clipped):
        return cen
    return shapely.get_point(shapely.shortest_line(clipped, cen), 0)


# ---------------------------------------------------------------- NFP: the pixels
_NFP_COLS = ["mergid", "y", "x", "DS", "DT", "frp", "sat", "datetime_utc"]
_nfp_cache = {}

def load_nfp(year, step):
    '''The VIIRS pixels newly detected burning at ONE overpass, in METRIC_CRS.

    NFP is the incremental product: consecutive steps share no coordinates, so one
    step's file is exactly "what was seen burning then" - unlike the Snapshot
    polygon, which is cumulative and cannot tell today's front from ground that
    burned a fortnight ago. DS/DT are the along-scan / along-track pixel dimensions
    in km, i.e. each detection's real footprint.

    Only the ONE matched step is read, never the whole year, so this stays cheap.
    Pixels are selected spatially here (distance to the report), never by mergid, so
    the stale-id repair the full pipeline needs does not apply.'''
    key = (year, step)
    if key not in _nfp_cache:
        fp = os.path.join(NFP_DIR.format(y=year), step + "_NFP.txt")
        if not os.path.exists(fp):
            _nfp_cache[key] = None
            return None
        try:
            d = pd.read_csv(fp, usecols=_NFP_COLS)
        except (ValueError, pd.errors.EmptyDataError):
            _nfp_cache[key] = None
            return None
        if d.empty:
            _nfp_cache[key] = None
            return None
        d["datetime_utc"] = pd.to_datetime(d["datetime_utc"])
        # x/y are EPSG:3571 metres everywhere EXCEPT 2021, which ships lon/lat.
        # Detected from the magnitudes rather than trusted from the spec (the same
        # check pyrocb/io.py makes), so every year lands in one frame.
        if d["x"].abs().max() <= 180 and d["y"].abs().max() <= 90:
            g = gpd.GeoSeries(shapely.points(d["x"].to_numpy(float), d["y"].to_numpy(float)),
                              crs=WGS84).to_crs(METRIC_CRS)
            d["x"], d["y"] = g.x.to_numpy(), g.y.to_numpy()
        _nfp_cache[key] = d
    return _nfp_cache[key]


def _pts_of(d):
    return shapely.points(d["x"].to_numpy(float), d["y"].to_numpy(float))


def _footprints(d):
    '''Union of each detection's real DS x DT footprint (km -> half-width in m).'''
    x, y = d["x"].to_numpy(float), d["y"].to_numpy(float)
    hx, hy = d["DS"].to_numpy(float) * 500.0, d["DT"].to_numpy(float) * 500.0
    return ensure_valid([shapely.union_all(
        shapely.box(x - hx, y - hy, x + hx, y + hy), grid_size=PRECISION_M)])[0]


def active_burn(year, step, pt_metric, radius_km, onset_utc, fallback_xy):
    '''Stages 2 and 3: find where the fire really was, then cut the analysis window.

    STAGE 2 - the position. The plain (unweighted) mean of the detection CENTRES
    inside the ladder's own buffer. Unweighted over centres is exactly "more pixels
    pull it closer": ten detections one side and one the other lands the point ten
    times nearer the ten. This replaces taking a centroid of the Snapshot polygon,
    which cannot work at all - that polygon is cumulative burned area carrying no
    pixel density whatsoever, and measured here it is a median 3.8x (p90 21x) larger
    than what was actually burning in the same circle at the same overpass, so its
    centre of mass tracks old cold ground instead of the front.
    Not FRP-weighted, deliberately: this is the centre of the burning AREA and a cold
    flank counts as much as the head - the same call D-49 made for v2.

    STAGE 3 - the analysis polygon. Re-centre on that corrected position and take
    every active pixel within ANALYSIS_RADIUS_KM of it, as the union of their real
    footprints. That is the burn credited with producing this pyroCb and the geometry
    downstream analysis should sample. Re-centring is the point: the ladder buffer is
    centred on the atlas coordinate, which is rounded to 0.1 deg (~11 km), so a
    window around the report is not the same window as one around the fire.

    No detections in the buffer (cloud, or a fire burning but not newly detected that
    overpass) -> no active-burn polygon, and the position falls back to the Snapshot
    polygon's weighted centroid. Flagged in v3_pos_source, never silently mixed.'''
    out = {"geometry": None, "active_area_km2": 0.0, "active_n_pixels": 0,
           "active_frp_sum": np.nan, "active_frp_p95": np.nan, "active_det_dt_h": np.nan,
           "active_max_abs_dt_h": np.nan, "v3_n_det_buffer": 0, "v3_pos_source": "polygon_fallback",
           "v3_x": fallback_xy[0], "v3_y": fallback_xy[1]}
    d = load_nfp(year, step)
    if d is None:
        return out
    pts = _pts_of(d)
    inb = d[shapely.dwithin(pt_metric, pts, radius_km * 1000.0)]
    if not len(inb):
        return out

    # Stage 2: the corrected position.
    cx, cy = float(inb["x"].mean()), float(inb["y"].mean())
    out.update(v3_x=cx, v3_y=cy, v3_pos_source="detections", v3_n_det_buffer=int(len(inb)))

    # Stage 3: the re-centred analysis window, capped at +/-12h of onset MAXIMUM -
    # a pixel can sit inside the 5 km circle and still be a much later (or earlier)
    # detection than what actually produced this pyroCb; the radius alone does not
    # bound time. For v3 this rarely binds (one overpass is already close to onset
    # by construction of the ladder), but it is enforced explicitly rather than
    # assumed, and checked in the verification cell below.
    off_h_all = (d["datetime_utc"] - onset_utc).dt.total_seconds().to_numpy() / 3600.0
    win = d[shapely.dwithin(shapely.points(cx, cy), pts, ANALYSIS_RADIUS_KM * 1000.0)
            & (np.abs(off_h_all) <= ANALYSIS_TIME_CAP_H)]
    if not len(win):
        return out
    geom = _footprints(win)
    off_h = (win["datetime_utc"] - onset_utc).dt.total_seconds().to_numpy() / 3600.0
    out.update(geometry=geom,
               active_area_km2=float(shapely.area(geom)) / 1e6,
               active_n_pixels=int(len(win)),
               active_frp_sum=float(win["frp"].sum()),
               active_frp_p95=float(np.percentile(win["frp"].to_numpy(float), 95)),
               active_det_dt_h=float(np.median(off_h)),
               active_max_abs_dt_h=float(np.max(np.abs(off_h))))   # for the +/-12h check below
    return out


def v4_burn(year, onset_utc, lon, pt_metric, radius_km):
    '''v4: pool the V4_N_STEPS time-closest overpasses, weighting pixels by time and
    normalising each step's contribution.

    v3 uses ONE overpass, so it sees one 12-hourly slice of a fire and nothing else.
    v4 widens that to the four closest in time and asks a sharper question: where was
    the burning that this pyroCb came out of, given that some of these looks are much
    closer to onset than others?

    Two weights, doing two different jobs:

    * TIME DECAY, per pixel: `exp(-|dt| / V4_TAU_H)` on the pixel's REAL
      `datetime_utc` (NFP carries it), not the step's nominal hour - so the decay is
      exact rather than quantised to half-days. tau = 12 h (one overpass interval)
      gives a pixel 1 h from onset ~6.6x the pull of one 24 h away.
    * PER-STEP NORMALISATION: each pixel's weight is divided by the number of pixels
      its own step contributed. Without it a step's influence would scale with how
      much it happened to detect, so a massive burn 24 h later would simply outvote a
      small one an hour before onset by sheer count. With it, each step contributes
      its own centroid with a weight set by time alone, and "a few pixels just before
      onset" can outweigh "a huge area a day later" - which is the intent.

    Equivalent formulation, and how to read the result: v4 is the time-weighted mean
    of the four per-step centroids. Steps with no pixels in the buffer simply drop
    out; if none of the four has any, v4 has no answer (NaN) rather than a guess.'''
    out = {"v4_geometry": None, "v4_x": np.nan, "v4_y": np.nan, "v4_area_km2": 0.0,
           "v4_n_pixels": 0, "v4_n_det": 0, "v4_n_steps_used": 0, "v4_steps": None,
           "v4_det_dt_h": np.nan, "v4_frp_sum": np.nan, "v4_max_abs_dt_h": np.nan}
    cands = closest_two_steps(onset_utc, lon, year, n=V4_N_STEPS)
    if not cands:
        return out

    xs, ys, ws, used, pooled = [], [], [], [], []
    for step, _dt_nominal in cands:
        d = load_nfp(year, step)
        if d is None:
            continue
        inb = d[shapely.dwithin(pt_metric, _pts_of(d), radius_km * 1000.0)]
        if not len(inb):
            continue
        # Real per-pixel offset from onset, then decay; divide by this step's own
        # count so the step contributes its centroid, not its population.
        dt_h = (inb["datetime_utc"] - onset_utc).dt.total_seconds().to_numpy() / 3600.0
        w = np.exp(-np.abs(dt_h) / V4_TAU_H) / len(inb)
        xs.append(inb["x"].to_numpy(float))
        ys.append(inb["y"].to_numpy(float))
        ws.append(w)
        used.append(step)
        pooled.append(inb)
    if not used:
        return out

    x = np.concatenate(xs); y = np.concatenate(ys); w = np.concatenate(ws)
    if w.sum() <= 0:
        return out
    cx, cy = float(np.average(x, weights=w)), float(np.average(y, weights=w))
    allp = pd.concat(pooled, ignore_index=True)
    out.update(v4_x=cx, v4_y=cy, v4_n_det=int(len(allp)),
               v4_n_steps_used=len(used), v4_steps=",".join(used))

    # Same stage-3 idea as v3: re-centre, then take the analysis window - but pooled
    # over the steps v4 actually used, so the polygon matches the position's evidence.
    # Capped at +/-12h of onset MAXIMUM: v4 pools up to V4_N_STEPS=4 overpasses, and
    # in a data-sparse stretch the 4th-closest can be well over a day from onset -
    # the radius alone would happily include a much later (or earlier) detection
    # just because it fell within 5 km. This is where the cap actually binds, unlike
    # v3's single-overpass case.
    off_h_allp = (allp["datetime_utc"] - onset_utc).dt.total_seconds().to_numpy() / 3600.0
    win = allp[shapely.dwithin(shapely.points(cx, cy), _pts_of(allp), ANALYSIS_RADIUS_KM * 1000.0)
              & (np.abs(off_h_allp) <= ANALYSIS_TIME_CAP_H)]
    if not len(win):
        return out
    geom = _footprints(win)
    off_h = (win["datetime_utc"] - onset_utc).dt.total_seconds().to_numpy() / 3600.0
    out.update(v4_geometry=geom,
               v4_area_km2=float(shapely.area(geom)) / 1e6,
               v4_n_pixels=int(len(win)),
               v4_max_abs_dt_h=float(np.max(np.abs(off_h))),
               v4_det_dt_h=float(np.median(off_h)),
               v4_frp_sum=float(win["frp"].sum()))
    return out
""")

md(r"""
`mergedate` ships as e.g. `2023-5-16-PM` -- year, non-padded month/day, then AM/PM
where `%p` needs upper-case `AM`/`PM` to parse; checked against the real values
above before writing the `format=` string.
""")
code(r"""
# Confirm the mergedate format assumption directly against real data before it's
# relied on above (pandas silently returns NaT on a mismatched format string).
# mergedate itself is legitimately null for a component that never merged into
# another (see step 4's table), so the check is "every *non-null* value parses",
# not "every value is non-null".
_first_match = match_tbl.index[match_tbl["match_status"] == "matched"][0]
_probe = load_raw(int(pcb.loc[_first_match, "year"]),
                   match_tbl.loc[_first_match, "match_step"])
_parsed = pd.to_datetime(_probe["mergedate"], format="%Y-%m-%d-%p", errors="coerce")
_nonnull = _probe["mergedate"].notna()
assert (_parsed.notna() == _nonnull).all(), "mergedate format assumption failed"
print("mergedate format OK,", _parsed.notna().sum(), "/", _nonnull.sum(), "non-null rows parsed")
""")

md(r"""
## 10. Assemble the final table

Run step 9 for every matched event, then join: full PCB record + matching
diagnostics (step 8) + every `snap_*` attribute (step 9, whole-fire) + the local
active-burn geometry clipped to the buffer that found it, and the whole matched
fire's footprint alongside it as `fire_geometry`. Unmatched events keep their PCB and
diagnostic columns with `snap_*` fields empty -- they are genuinely `N/A`, not
dropped from the table.
""")
code(r"""
matched_idx = match_tbl.index[match_tbl["match_status"] == "matched"]
agg_rows = []
for i in matched_idx:
    row = match_tbl.loc[i]
    year = int(pcb.loc[i, "year"])
    # Stage 1 result -> the dominant fire in the buffer + its cumulative footprint.
    loc = local_active_burn(year, row["match_step"], pts_metric.values[i],
                            float(row["search_radius_km"]))
    # Stages 2 and 3 -> the pixel-weighted position, then the 5 km analysis window
    # re-centred on it. This is what `geometry` becomes.
    act = active_burn(year, row["match_step"], pts_metric.values[i],
                      float(row["search_radius_km"]), pcb.loc[i, "datetime"],
                      (loc["fb_x"], loc["fb_y"]))
    # Aggregate the fire the buffer says this pyroCb belongs to, not the nearest one.
    # v4: the same question asked over the 4 time-closest overpasses, time-weighted.
    v4 = v4_burn(year, pcb.loc[i, "datetime"], pcb.loc[i, "lon"], pts_metric.values[i],
                 float(row["search_radius_km"]))
    agg_rows.append({"__idx": i,
                     **aggregate_matched_fire(year, row["match_step"], loc["mergid"]),
                     **loc, **act, **v4})
agg_tbl = pd.DataFrame(agg_rows).set_index("__idx")

final = pcb.join(match_tbl.drop(columns="pcb_id").set_index(pcb.index)).join(agg_tbl)
final = gpd.GeoDataFrame(final, geometry="geometry", crs=METRIC_CRS)

# The shipped v3 position, in WGS84, computed ONCE here. Everything downstream (the
# map marker, the GOES star) reads these columns instead of deriving a position of
# its own - that duplication is what let the map plot a point the notebook never
# endorsed. NaN stays NaN: an unmatched event has no position.
for _v in ("v3", "v4"):
    _ok = final[f"{_v}_x"].notna()
    final[f"{_v}_lat"] = np.nan
    final[f"{_v}_lon"] = np.nan
    if _ok.any():
        _ll = gpd.GeoSeries(shapely.points(final.loc[_ok, f"{_v}_x"].to_numpy(float),
                                           final.loc[_ok, f"{_v}_y"].to_numpy(float)),
                            crs=METRIC_CRS).to_crs(WGS84)
        final.loc[_ok, f"{_v}_lat"] = _ll.y.to_numpy()
        final.loc[_ok, f"{_v}_lon"] = _ll.x.to_numpy()
print(final.shape)
final.columns.tolist()
""")

md("## 11. Verification -- double- and triple-checks")
code(r"""
m = final[final["match_status"] == "matched"]

# 1) every match respects its own reported radius
assert (m["match_dist_km"] <= m["search_radius_km"] + 1e-6).all(), "distance exceeds its search radius"

# 2) the +/-12h cap actually held for every rank-2 match (the point the prompt
#    explicitly asked to double-check)
r2 = m[m["match_rank"] == 2]
assert (r2["match_dt_h"].abs() <= TIME_CAP_H + 1e-9).all(), "a rank-2 match exceeds the 12h cap"
print(f"rank-2 matches: {len(r2)}, max |dt_h| = {r2['match_dt_h'].abs().max() if len(r2) else float('nan'):.3f} h (cap {TIME_CAP_H})")

# 3) radius is always one of the ladder rungs, never something else
assert m["search_radius_km"].isin(RADII_KM).all()

# 4) independently re-derive snap_n_pixels / snap_pixden from a fresh read, for a
#    random sample, and compare -- catches any bug in aggregate_matched_fire()
chk = m.sample(min(15, len(m)), random_state=0)
for i, row in chk.iterrows():
    g = load_raw(int(pcb.loc[i, "year"]), row["match_step"])
    c = g[g["mergid"] == row["mergid"]]
    assert row["snap_n_pixels"] == c["n_pixels"].sum()
    exp_pixden = c["n_pixels"].sum() / c["farea"].sum()
    assert abs(row["snap_pixden"] - exp_pixden) < 1e-9
print(f"re-derivation check passed on {len(chk)} random matched events")

# 5) no duplicate events, one row per pcb_id
assert final["pcb_id"].is_unique

# 6) full precision preserved: datetime keeps sub-day components where present,
#    lat/lon/dt_h/dist_km are float64, not rounded ints
assert final["datetime"].dtype.kind == "M"
for c in ["lat", "lon", "match_dt_h", "match_dist_km", "snap_farea_km2"]:
    assert final[c].dtype == np.float64, c

# 7) THE POSITION. It must sit inside the ladder buffer that found the event: the
#    stage-2 position is the mean of detection centres drawn from that buffer, and a
#    mean of points inside a disc is inside the disc. This is the invariant that was
#    violated when the position came from the whole fire's footprint - it sat up to
#    150 km from its own pyroCb (worst case 2023-G90288).
_pos = shapely.points(m["v3_x"].to_numpy(float), m["v3_y"].to_numpy(float))
_pts = pts_metric.values[m.index.to_numpy()]
_off_km = shapely.distance(_pos, _pts) / 1000.0
assert (_off_km <= m["search_radius_km"].to_numpy() + 1e-6).all(), "v3 position escaped its buffer"

# 8) THE CUMULATIVE CLIP stays inside its buffer (area can not exceed the disc).
#    Note it pools every fire reaching in, while snap_farea_km2 is only the credited
#    one, so a clip legitimately exceeds that single fire's area.
assert (m["clip_area_km2"] <= np.pi * m["search_radius_km"] ** 2 + 1e-6).all(), "clip exceeds its buffer"

# 9) THE ANALYSIS POLYGON is the stage-3 window: active pixel footprints within
#    ANALYSIS_RADIUS_KM of the CORRECTED centre - so it is bounded by that radius
#    plus one pixel's half-width, NOT by the ladder buffer around the report.
act = m[m["active_n_pixels"] > 0]
_ac = shapely.points(act["v3_x"].to_numpy(float), act["v3_y"].to_numpy(float))
_slack = 2000.0                      # generous: a VIIRS pixel is ~0.4-1.6 km across
assert shapely.covered_by(
    act.geometry.values,
    shapely.buffer(_ac, ANALYSIS_RADIUS_KM * 1000.0 + _slack)).all(), "analysis window escaped its radius"
assert (act["active_area_km2"] <= np.pi * (ANALYSIS_RADIUS_KM + 2) ** 2).all(), "analysis area too large"

# 10) Every event with detections must be positioned FROM those detections, and every
#     fallback must genuinely have had none - the two must never be silently mixed.
assert ((m["v3_pos_source"] == "detections") == (m["v3_n_det_buffer"] > 0)).all(), "position source mislabelled"

_det = m[m["v3_pos_source"] == "detections"]
print(f"v3 position offset from report: median {np.median(_off_km):.2f} km, "
      f"max {_off_km.max():.2f} km (ladder max {max(RADII_KM)} km)")
print(f"position from VIIRS pixels: {len(_det)}/{len(m)}; "
      f"{int((m['v3_pos_source'] == 'polygon_fallback').sum())} fell back to the polygon centroid")
print(f"active burn in the {ANALYSIS_RADIUS_KM:g} km window: median "
      f"{act['active_area_km2'].median():.1f} km2 / {int(act['active_n_pixels'].median())} px "
      f"| cumulative clip for comparison: {m['clip_area_km2'].median():.1f} km2")
print(f"credited fire differs from the merely nearest one for "
      f"{int((m['mergid'] != m['ladder_mergid']).sum())} events; "
      f"{int((m['n_fires_in_buffer'] > 1).sum())} events have >1 fire in the buffer")

# 11) v4 must obey the same geometric invariants, and its position must be a
#     genuine weighted mean of the pixels it pooled - so it lies inside the ladder
#     buffer too (a weighted mean of points in a disc is in the disc).
v4 = m[m["v4_x"].notna()]
if len(v4):
    _v4p = shapely.points(v4["v4_x"].to_numpy(float), v4["v4_y"].to_numpy(float))
    _v4off = shapely.distance(_v4p, pts_metric.values[v4.index.to_numpy()]) / 1000.0
    assert (_v4off <= v4["search_radius_km"].to_numpy() + 1e-6).all(), "v4 position escaped its buffer"
    assert (v4["v4_n_steps_used"] <= V4_N_STEPS).all(), "v4 pooled more steps than allowed"
    _v4g = v4[v4["v4_n_pixels"] > 0]
    assert shapely.covered_by(
        _v4g["v4_geometry"].values,
        shapely.buffer(shapely.points(_v4g["v4_x"].to_numpy(float),
                                      _v4g["v4_y"].to_numpy(float)),
                       ANALYSIS_RADIUS_KM * 1000.0 + _slack)).all(), "v4 window escaped its radius"
    _v3v4 = shapely.distance(
        shapely.points(v4["v3_x"].to_numpy(float), v4["v3_y"].to_numpy(float)), _v4p) / 1000.0
    print(f"v4: {len(v4)}/{len(m)} events, median {v4['v4_n_steps_used'].median():.0f} steps pooled, "
          f"median {v4['v4_area_km2'].median():.1f} km2 / {int(v4['v4_n_pixels'].median())} px")
    print(f"v4 vs v3 position: median {np.median(_v3v4):.2f} km apart, max {_v3v4.max():.2f} km")

print("all checks passed")
final["match_status"].value_counts()
""")

md(r"""
## 12. Spot-check a few matched events by eye

Every intermediate value the algorithm used (both candidates' steps, deltas and
distances, which one won and at what radius) is kept in the table specifically so a
match can be audited by eye without re-running anything.
""")
code(r"""
cols = ["pcb_id", "datetime", "lat", "lon", "mergid", "match_status",
        "cand1_step", "cand1_dt_h", "cand1_dist_km",
        "cand2_step", "cand2_dt_h", "cand2_dist_km", "cand2_eligible_12h",
        "match_step", "match_dt_h", "search_radius_km", "match_dist_km", "match_rank"]
final[cols].sample(8, random_state=2)
""")

md(r"""
## 13. Unmatched events

Real cases where no snapshot within 20 km, at either candidate time, has an active
fire -- worth a look before treating them as "no data".
""")
code(r"""
u = final[final["match_status"] == "unmatched"]
print(f"{len(u)} / {len(final)} events unmatched ({100*len(u)/len(final):.1f}%)")
u[["pcb_id", "datetime", "lat", "lon", "region",
   "cand1_step", "cand1_dt_h", "cand1_dist_km",
   "cand2_step", "cand2_dt_h", "cand2_dist_km", "cand2_eligible_12h"]]
""")

md(r"""
## 14. Export

Geometry in `METRIC_CRS` (EPSG:3571, equal-area), matching every other geometry
output in this project (`pyrocb_active_areas.gpkg`, final perimeters, ignitions).
No column is rounded on the way out -- `to_csv`/`to_file` write the native
`float64`/`datetime64` precision as-is.
""")
code(r"""
# Delete first: to_file APPENDS layers to an existing GeoPackage, so a rerun that
# renames or drops a layer would silently leave the previous run's version behind,
# still readable (and first in the layer list) for anyone who opens the file without
# naming a layer. Starting from a clean file makes the output match this code exactly.
if os.path.exists("pcb_matching_clean.gpkg"):
    os.remove("pcb_matching_clean.gpkg")

# Layer 1 (default): one row per pyroCb, geometry = the LOCAL active burn (the
# matched snapshot clipped to the search buffer). Layer 2: the same rows carrying
# the WHOLE matched fire footprint instead, so nothing is discarded and the two can
# be compared directly. One GeoPackage holds both.
_geoms = ["geometry", "v4_geometry", "snap_clip_geometry", "fire_geometry"]
# Layer 1 (default) - THE analysis polygon: active VIIRS pixel footprints within
# ANALYSIS_RADIUS_KM of the corrected, pixel-weighted centre (stages 2-3).
final.drop(columns=[c for c in _geoms if c != "geometry"]).to_file(
    "pcb_matching_clean.gpkg", layer="pcb_active_burn", driver="GPKG")
# Layer 2 - the cumulative Snapshot polygon clipped to the ladder buffer. Kept for
# comparison; it is ~3.8x the active area and must never be read as "actively burning".
(final[final["snap_clip_geometry"].notna()]
    .drop(columns=[c for c in _geoms if c != "snap_clip_geometry"])
    .set_geometry("snap_clip_geometry")
    .to_file("pcb_matching_clean.gpkg", layer="snapshot_clip", driver="GPKG"))
# Layer 3 - v4's analysis window: the same construction pooled over the 4
# time-closest overpasses, time-weighted and normalised per step.
(final[final["v4_geometry"].notna()]
    .drop(columns=[c for c in _geoms if c != "v4_geometry"])
    .set_geometry("v4_geometry")
    .to_file("pcb_matching_clean.gpkg", layer="pcb_active_burn_v4", driver="GPKG"))
# Layer 4 - the whole matched fire, unclipped.
(final[final["fire_geometry"].notna()]
    .drop(columns=[c for c in _geoms if c != "fire_geometry"])
    .set_geometry("fire_geometry")
    .to_file("pcb_matching_clean.gpkg", layer="matched_fire_full", driver="GPKG"))
final.drop(columns=_geoms).to_csv("pcb_matching_clean.tsv", sep="\t", index=False)
print("wrote pcb_matching_clean.gpkg (layers pcb_active_burn, snapshot_clip,"
      " pcb_active_burn_v4, matched_fire_full) and pcb_matching_clean.tsv:", final.shape)
""")

nb["cells"] = cells
nbf.write(nb, "PCB_MATCHING_CLEAN.ipynb")
print("wrote PCB_MATCHING_CLEAN.ipynb:", len(cells), "cells")
