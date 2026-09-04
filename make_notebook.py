"""Generate 'Data exploration.ipynb'.

The notebook narrates; the analysis lives in the `pyrocb` package so it can be
imported, tested and profiled.  Running the notebook writes:

    fire_pcb_map.html              interactive map, per-year toggles + growth bands
    pyrocb_fire_attributes.tsv     one row per pyroCb event, ~100 columns
    pyrocb_fire_attributes.xlsx
    pyrocb_columns.tsv             data dictionary

Run:  python make_notebook.py
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
def md(t):   cells.append(nbf.v4.new_markdown_cell(t.strip("\n")))
def code(t): cells.append(nbf.v4.new_code_cell(t.strip("\n")))


md(r"""
# Boreal FireAtlas × pyroCb — growth bands and the event attribute table

Each pyroCb event is associated with the fire that produced it, that fire's footprint is
split at **every** pyroCb it generated, and the result is exported as one row per event.

| Dataset | Files | What it is |
|---|---|---|
| **pyroCb events** | `PCBAtlas.xlsx` | NRL inventory (Peterson et al. 2025), boreal subset. One row per **event**; `lat`/`lon` = approximate **fire** location, stored to 0.1°. 2013–2023. |
| **Final perimeters** | `final_perims/.../final_viirs{Y}.gpkg` | Season-final burned area per fire. EPSG:4326. |
| **Ignitions** | `ignitions/.../ign{Y}.gpkg` | First-detection footprint per fire. EPSG:3571. |
| **Snapshots** | `{Y}/{Y}/Snapshot/YYYYMMDD{AM,PM}.gpkg` | 12-hourly cumulative state, **active fires only**. EPSG:3571. |
| **NFP** | `{Y}/{Y}/NFP/YYYYMMDD{AM,PM}_NFP.txt` | 12-hourly **incremental** new fire pixels with a real `datetime_utc`. EPSG:3571. |

## The one thing to know: `mergid` is the join key

`ignitions.mergid` == `final.fireid` == `final.mergid` == `NFP.mergid` == `snapshot.mergid`,
a perfect 1:1 (10,707 fires in 2023).

* **Snapshot `fireid` is unsafe** — ~5% of snapshot rows carry a `fireid` absent from the final
  file, and matching on it can attach a small wrong polygon.
* **NFP `fireid` is a per-file row counter** (`== fid == 0..N-1`), not a fire id at all.
* `mergid` groups the separately-ignited fires that later grew together: 1,444 of 12,151 2023
  ignition rows are non-primary members of a merge complex, and one `mergid` absorbs **115**.

## `farea` is cumulative burned area, not actively-burning area

A perimeter does **not** shrink when parts of the fire stop burning. Checked over all 608
snapshots × 10,707 fires (2023): max-snapshot `n_pixels` ÷ final is **exactly 1.0 for every
fire**, and `farea` agrees with the final within 0.01% for 99.92% of them. The attribute table
says the same — `n_pixels` is *"Number of **total** fire pixels"*, `n_newpixels` is *"at current
time step"*, and the burning front is carried separately by **`flinelen`**, `isactive` and
`t_inactive`. The 9 exceptions are alpha-shape artefacts on 2–14 pixel fires.

The analysis lives in the `pyrocb` package; this notebook drives it.
""")

md("## 1. Setup")
code(r"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import pyrocb
from pyrocb import (config as C, io, timebins as T, matching, bands, active,
                    variables, pipeline)

pd.set_option("display.max_columns", 60)
pd.set_option("display.width", 200)
print("pyrocb package loaded")
""")

md(r"""
## 2. Configuration

`YEARS = None` resolves to *years with fire data on disk **and** pyroCb records*. That is
2013–2023 — **2012 drops out on its own**, because the pyroCb inventory does not cover it.

`SNAPSHOT_YEARS` gates the heavy snapshot/NFP analysis (matching, bands, per-event variables).
Set it to `None` to run every year. Everything is cached per (stage, year), so going from
`[2023]` to all years reuses 2023 and only computes the rest.

All metric work happens in **EPSG:3571 — "WGS 84 / North Pole LAEA Bering Sea"**: equal-area,
and centred on lon 180, so the dateline (where the Chukotka fires sit) is the *least* distorted
part of the projection. It is **not** EASE-Grid 2.0 North; that is EPSG:6931.
""")
code(r"""
pcb_all = io.read_pcb()
YEARS    = C.resolve_years(set(pcb_all["year"]))      # map years
SNAP_YRS = C.resolve_snapshot_years(YEARS)            # heavy-analysis years

print(f"pyroCb inventory : {len(pcb_all)} events, {pcb_all.year.min()}-{pcb_all.year.max()}")
print(f"map years        : {YEARS}")
print(f"analysis years   : {SNAP_YRS}   (set C.SNAPSHOT_YEARS = None for all)")
print(f"match cap        : {C.MATCH_CAP_KM} km | sliver {C.SLIVER_KM2} km2 | snap grid {C.PRECISION_M} m")
pd.crosstab(pcb_all.year, pcb_all.region, margins=True)
""")

md(r"""
### Event IDs

`pcb_id = PCB_YYYY_MMDD_HHMM` (UTC). This **collides**: 484 events give only 449 distinct
timestamps, because a pyroCb outbreak blows up several separate fires in one afternoon —
`PCB_2023_0729_2220` is 2023WB7, 2023WB18 *and* MWF-043-2023, hundreds of km apart. Colliding
timestamps get a deterministic `a`/`b`/`c` suffix ordered by (lat, lon), so ids are unique and
stable across rebuilds.
""")
code(r"""
print("unique pcb_id:", pcb_all.pcb_id.nunique(), "of", len(pcb_all), "events")
pcb_all[pcb_all.pcb_id.str.contains("2023_0729_2220")][["pcb_id", "fire_name", "region", "lat", "lon"]]
""")

md(r"""
## 3. Matching each pyroCb to a fire

**Containment first, then nearest active fire within `MATCH_CAP_KM`.**

Matching against the *snapshot* rather than the season-final perimeter matters: the final
footprint of a big fire can be 100 km across, so a point inside it says nothing about whether
the fire had reached that spot on the day, and a fire that ignited weeks later would still
swallow the point. The snapshot at onset is the fire as it actually was.

Events unmatched against the snapshot get a second chance against the final perimeters,
restricted to fires whose start–end interval brackets onset — a fire obscured by cloud vanishes
from a snapshot even though it is burning.

The floor on match quality is the atlas: `lat`/`lon` are stored to 0.1°, a worst-case
displacement of ~7.8 km, so `match_dist_km` above ~5 km is provisional.
""")
code(r"""
matches = {y: pipeline.match_cached(y) for y in SNAP_YRS}
m = pd.concat(matches.values(), ignore_index=True)
print(m.match_method.value_counts().to_string())
print(f"\ndistinct fires: {m.mergid.dropna().nunique()}")
m[["pcb_id", "fire_name", "region", "match_method", "match_dist_km", "mergid"]].head(10)
""")

md(r"""
## 4. Concentric growth bands

For a fire with N distinct pyroCb half-day steps we build N+1 bands:

```
S_0     = V(snap(t_0))               # burned before the first pyroCb
S_i     = S_{i-1} ∪ V(snap(t_i))     # cumulative, monotone BY CONSTRUCTION
S_final = S_{N-1} ∪ V(final_perim)
B_i     = S_i \ S_{i-1}
```

Building a running **union** rather than differencing raw snapshots is what makes the bands
disjoint, non-negative, and exactly tiling `S_final`. Consecutive snapshots are not perfectly
nested (re-polygonisation loses ~0.17% of the earlier polygon on average), so plain differencing
would give negative-area bands and rings that punch holes in each other on the map.

Geometry hygiene: `make_valid` (**never** `buffer(0)` — it silently deletes loops and changes
area, fatal when the output *is* an area), then both operands snapped to a common 1 m grid before
every boolean, which is what makes the difference topologically exact. 1 m is 1/375 of a VIIRS pixel.

### Which step counts as "before onset"

FEDS bins observations into two **local-solar** overpass groups, not calendar halves: measured
against real NFP timestamps, AM falls at local 0.5–7.0 h and PM at local 7.1–14.4 h, with **no
observation between local 14.4 h and 00.5 h**. pyroCb onset peaks at local ~15.0 h — inside that
blind window.

So we use the observations themselves (**Tier A**): the latest step whose *median* detection time
for that fire precedes onset. The median, not the max, so one stray late pixel cannot disqualify
an overpass whose bulk clearly precedes onset. `obs_lag_h` records how stale the resulting
pre-band is.
""")
code(r"""
band_geom, band_ev, active_geom = {}, {}, {}
for y in SNAP_YRS:
    band_geom[y], band_ev[y] = pipeline.bands_cached(y, matches[y])
bd = pd.concat(band_geom.values(), ignore_index=True) if band_geom else None
print(f"\n{len(bd)} bands over {bd.mergid.nunique()} fires")
bd.groupby("mergid").n_bands.first().value_counts().sort_index().rename_axis("bands per fire").to_frame("fires")
""")

md(r"""
### Band integrity

The bands must be pairwise disjoint and sum to `S_final`. This is an assertion, not a hope —
if it fails, every area column downstream is wrong.
""")
code(r"""
import shapely
overlap = 0
for mid, g in bd.groupby("mergid"):
    gs = [x for x in g.geometry if x is not None and not x.is_empty]
    if gs and abs(sum(x.area for x in gs) - shapely.union_all(gs).area) > max(5e4, 1e-3 * sum(x.area for x in gs)):
        overlap += 1
chk = bd.groupby("mergid").agg(s=("area_km2", "sum"), c=("cum_area_km2", "first"))
chk["rel_err"] = (chk.s - chk.c).abs() / chk.c.clip(lower=1e-9)

print(f"fires with overlapping bands : {overlap} of {bd.mergid.nunique()}")
print(f"band sum vs S_final          : max rel err {chk.rel_err.max():.2e}, fires >0.1% off: {(chk.rel_err > 1e-3).sum()}")
print(f"negative-area bands          : {int((bd.area_km2 < 0).sum())}")
print(f"sliver area dropped          : {bd.groupby('mergid').qc_sliver_dropped_km2.first().sum():.3f} km2")
""")

md(r"""
## 5. Per-event variables

One row per pyroCb event — **every** event, including unmatched ones (null fire columns,
`match_method="unmatched"`). The unmatched set is non-random, so dropping it silently would bias
the sample toward well-detected boreal forest fires.

**Deduplication caveat, worth reading before you plot anything.** Two pyroCb on the same fire in
the same half-day share a band, so every *step-derived* column (`area_pre_km2`, `area_post_km2`,
`frac_burned_before`, `band_id`) is **identical** for them — genuinely unresolvable at 12-hourly
sampling. `dedup_collapsed` flags those rows. The NFP-derived columns use exact `datetime_utc`
and *do* separate them. Without this in mind, a scatter of `area_pre_km2` shows vertical stacks
of identical values and looks like a bug.
""")
code(r"""
_ev = []
for y in SNAP_YRS:
    e = pipeline.events_cached(y, matches[y], band_ev[y])
    active_geom[y] = pipeline.active_cached(y, e)
    s = active.summarise(active_geom[y])
    if len(s):
        e = e.merge(s, on="pcb_id", how="left")
    _ev.append(e)
events = pd.concat(_ev, ignore_index=True)
active_all = pd.concat([active_geom[y] for y in SNAP_YRS if len(active_geom[y])],
                       ignore_index=True)
print(f"{len(events)} events x {events.shape[1]} columns")
print(f"dedup_collapsed rows: {int(events.dedup_collapsed.fillna(False).sum())}")
events[["pcb_id", "fire_name", "mergid", "pcb_index", "n_pcb_total", "n_pcb_before", "n_pcb_after",
        "n_ign_total", "n_ign_before", "n_ign_after", "area_pre_km2", "area_post_km2",
        "frac_burned_before", "obs_lag_h"]].head(12)
""")

md(r"""
### Does one pyroCb beget more?

`n_pcb_before` / `n_pcb_after` count the pyroCb on the *same fire* either side of each event, so
the table can answer that question directly.
""")
code(r"""
mt = events[events.mergid.notna()]
per_fire = mt.groupby("mergid").agg(n_pcb=("n_pcb_total", "first"),
                                    area=("fire_farea_km2", "first"),
                                    n_ign=("n_ign_total", "first"),
                                    merge_span=("ign_merge_span_days", "first"))
print(per_fire.n_pcb.value_counts().sort_index().rename_axis("pyroCb on the fire").to_frame("fires").to_string())
print(f"\nSpearman(final fire area, number of pyroCb) = {per_fire.area.corr(per_fire.n_pcb, method='spearman'):.3f}")
print(f"Spearman(ignitions in complex, n pyroCb)    = {per_fire.n_ign.corr(per_fire.n_pcb, method='spearman'):.3f}")
per_fire.sort_values("n_pcb", ascending=False).head(8)
""")

md(r"""
### Export

`pyrocb_fire_attributes.tsv` / `.xlsx`, plus `pyrocb_columns.tsv` — a data dictionary naming every
column's units, definition and source file, for handing to collaborators.
""")
code(r"""
_ = pipeline.export(events, active_all)
""")

md(r"""
### Caveats to carry into any analysis

1. **Merged fires break the concentric premise.** Bands assume one spreading front, but the
   biggest fires — exactly the pyroCb ones — are the most likely to be merge complexes.
   `n_ign_total` and `ign_merge_span_days` tell you when `band_rate_km2_h` is not a spread rate.
2. **`frac_burned_before` is biased low, correlated with time of day** — worst for early-afternoon
   pyroCb, so not random. Use `obs_lag_h` as a covariate or restrict to `obs_lag_h < 6 h`.
3. **0.1° coordinate rounding** produces plausible-looking but wrong band stacks with no visible
   symptom. Treat `match_dist_km > 5 km` as provisional.
4. **FRP sums are not energy.** Headline `frp_p95_prev24h`; always carry `n_obs_prev{W}h` so
   "zero because no overpass" stays distinguishable from "zero because no fire".
5. **Unmatched events are non-random** — report them with a regional breakdown, never as a total.
""")
code(r"""
print("unmatched by region:")
print(events[events.match_method == "unmatched"].region.value_counts().to_string())
lag = events.obs_lag_h.dropna()
print(f"\nobs_lag_h: median {lag.median():.1f} h, {(lag < 6).sum()}/{len(lag)} events under 6 h")
print(f"match_dist_km > 5 km (provisional): {(events.match_dist_km > 5).sum()}")
print(f"merge complexes (n_ign_total > 1): {(events.n_ign_total > 1).sum()} events")
""")

md(r"""
## 6. Background map layers

Perimeters and ignitions for **every** map year, thinned for the browser. `doy` (day of year) is
precomputed on every feature so the map can colour by season without parsing dates in JS.

At ~330 m (one VIIRS pixel) the 200 k background perimeters compress to a workable size. The
growth bands get a far finer 50 m tolerance — there are only ~180 of them and the shape of a run
is the whole point of that layer.
""")
code(r"""
import base64, gzip, json, geopandas as gpd, shapely

def doy_of(df, pre="tst"):
    d = pd.to_datetime(dict(year=df[f"{pre}_year"], month=df[f"{pre}_month"],
                            day=df[f"{pre}_day"]), errors="coerce")
    return d.dt.dayofyear

def pack(gj):
    raw = json.dumps(gj, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode(), len(raw), len(gzip.compress(raw, 9))

def smooth_for_map(gser, close_m, simp_m):
    # Map-only geometry reduction. The active-burn polygon is a union of thousands of
    # individual VIIRS pixel rectangles (~125 kB of coordinates each). The map only
    # needs to show roughly where; the exact sampling geometry lives in the GeoPackage.
    # Close the pixel lattice into solid blobs (buffer +/- close_m), then simplify.
    g = gser.buffer(close_m, join_style="mitre").buffer(-close_m, join_style="mitre")
    return g.simplify(simp_m, preserve_topology=True)

payloads, doys = {}, []
for y in YEARS:
    per = io.read_final(y, metric=False)
    per["geometry"] = per.geometry.simplify(C.SIMPLIFY_TOL, preserve_topology=True)
    per["geometry"] = shapely.set_precision(per.geometry.values, C.COORD_PRECISION)
    per = per[~per.geometry.is_empty & per.geometry.notna()].copy()
    per["doy"] = doy_of(per)
    per_out = per[["mergid", "farea", "duration", "n_pixels", "lcc_final", "doy", "geometry"]] \
        .rename(columns={"farea": "farea_km2", "duration": "duration_days"}).round({"farea_km2": 3})
    doys.append(per_out[["doy", "farea_km2"]].dropna())   # feeds the DOY-window fit below

    ign = io.read_ignitions(y).to_crs(C.WGS84)
    ign["doy"] = doy_of(ign)
    ign["n_ign_total"] = ign.groupby("mergid")["mergid"].transform("size")
    ign["geometry"] = ign.geometry.representative_point()
    ign_out = ign[["mergid", "ign_index", "n_ign_total", "ftype", "n_pixels",
                   "farea", "doy", "geometry"]].rename(columns={"farea": "farea_km2"})

    p = pcb_all[pcb_all.year == y]
    ev_y = events[events.year == y] if len(events) else pd.DataFrame()
    if len(ev_y):
        p = ev_y
    pg = gpd.GeoDataFrame(p.copy(), geometry=gpd.points_from_xy(p.lon, p.lat), crs=C.WGS84)
    for c in pg.columns:
        if pd.api.types.is_datetime64_any_dtype(pg[c]):
            pg[c] = pg[c].dt.strftime("%Y-%m-%d %H:%M")
    pg = pg[[c for c in ("pcb_id", "datetime", "doy", "fire_name", "region", "status",
                         "max_inject_alt", "trop_relative_alt", "pcb_min_bt", "mergid",
                         "match_method", "match_dist_km", "pcb_index", "n_pcb_total",
                         "n_pcb_before", "n_pcb_after", "n_ign_total", "area_pre_km2",
                         "area_post_km2", "frac_burned_before", "obs_lag_h", "geometry")
                if c in pg.columns]]

    # Fires credited with a pyroCb, for the "only fires with a pyroCb" filter.
    linked = set(ev_y.loc[ev_y.match_method != "unmatched", "mergid"].dropna().astype(int)) \
        if len(ev_y) else set()
    per_out["pcb_linked"] = per_out["mergid"].isin(linked)
    ign_out["pcb_linked"] = ign_out["mergid"].isin(linked)   # mergid, NOT fireid

    if C.MAP_ONLY_PCB_FIRES:
        # Ship only the analysis population. Essential now ignitions are DOM star
        # markers - 180k SVG icons would lock the browser; ~1.5k is nothing.
        per_out = per_out[per_out.pcb_linked].copy()
        ign_out = ign_out[ign_out.pcb_linked].copy()

    b64 = {}
    b64["per"], _, gp = pack(json.loads(per_out.to_json(drop_id=True)))
    b64["ign"], _, gi = pack(json.loads(ign_out.to_json(drop_id=True)))
    b64["pcb"], _, _  = pack(json.loads(pg.to_json(drop_id=True)))

    bg = band_geom.get(y)
    if bg is not None and len(bg):
        bg = bg.copy()
        bg["geometry"] = smooth_for_map(bg.geometry, 60.0, C.BAND_SIMPLIFY_M * 3)
        bg = gpd.GeoDataFrame(bg, geometry="geometry", crs=C.METRIC_CRS).to_crs(C.WGS84)
        bg = bg[~bg.geometry.is_empty & bg.geometry.notna()]
        bg["doy"] = pd.to_numeric(bg["step_hi"].str[4:8].str.slice(0, 2), errors="coerce")
        bg["doy"] = pd.to_datetime(bg["step_hi"].str[:8], format="%Y%m%d",
                                   errors="coerce").dt.dayofyear
        b64["bnd"], _, _ = pack(json.loads(bg.drop(columns=[c for c in ("step_lo",) if c in bg]).to_json(drop_id=True)))
    else:
        b64["bnd"], _, _ = pack({"type": "FeatureCollection", "features": []})

    ag = active_geom.get(y)
    if ag is not None and len(ag):
        ag = ag[ag.window_h == C.ACTIVE_HEADLINE_H].copy()
        # +/- 200 m closes the ~375 m pixel lattice into solid blobs before simplifying.
        ag["geometry"] = smooth_for_map(ag.geometry, 200.0, 150.0)
        ag = gpd.GeoDataFrame(ag, geometry="geometry", crs=C.METRIC_CRS).to_crs(C.WGS84)
        ag = ag[~ag.geometry.is_empty & ag.geometry.notna()]
        ag["doy"] = pd.to_datetime(ag["datetime"]).dt.dayofyear
        ag["datetime"] = pd.to_datetime(ag["datetime"]).dt.strftime("%Y-%m-%d %H:%M")
        for c in ("t_first", "t_last"):
            ag[c] = pd.to_datetime(ag[c]).dt.strftime("%Y-%m-%d %H:%M")
        b64["act"], _, _ = pack(json.loads(ag.to_json(drop_id=True)))
    else:
        b64["act"], _, _ = pack({"type": "FeatureCollection", "features": []})

    payloads[y] = b64
    print(f"  {y}: {len(per_out):>6} perims  {len(ign_out):>6} ignitions  {len(pg):>3} pyroCb"
          f"  {0 if bg is None else len(bg):>3} bands   packed {(gp+gi)/1e6:5.2f} MB", flush=True)

# Fit the DOY window to the whole dataset, weighted by BURNED AREA. Weighting matters:
# an unweighted percentile over fire counts is dragged to (16, 349) by tiny January crop
# fires, which would waste the entire gradient on months with no boreal burning. The
# area-weighted 1-99% span is (69, 281) - 98% of all burned area - and we widen it to
# cover every pyroCb so no event clamps. Refitting on new years extends it automatically.
_d = pd.concat(doys).sort_values("doy")
_cw = _d.farea_km2.cumsum() / _d.farea_km2.sum()
lo = int(_d.doy[_cw >= 0.01].iloc[0])
hi = int(_d.doy[_cw <= 0.99].iloc[-1])
DOY_RANGE = (min(lo, int(pcb_all.doy.min())), max(hi, int(pcb_all.doy.max())))
print(f"\nDOY colour window (area-weighted, whole dataset): {DOY_RANGE}"
      f"  [{100 * _d.farea_km2[(_d.doy >= DOY_RANGE[0]) & (_d.doy <= DOY_RANGE[1])].sum() / _d.farea_km2.sum():.1f}% of burned area]")
""")

md(r"""
## 7. The map

`fire_pcb_map.html` — one self-contained file.

* **Year chips** toggle years; each year's payload is decompressed and turned into Leaflet layers
  only the first time you switch it on, which is why a 12-year file opens fast.
* **Colour by: Year / Day-of-year.** The DOY ramp reuses the same cool→warm palette, interpolated
  over a window fitted to the *whole* dataset — so it extends automatically as years are added,
  and colours stay comparable between years. Select a single year to watch the season progress.
* **Growth bands** draw the pyroCb split. Band order is encoded in **lightness** (pale core =
  burned before the first pyroCb, dark = late runs) so hue stays free for the colour mode.
  The outermost band — "everything after the last pyroCb, at unknown times" — is drawn fainter,
  because it is not a growth episode of known duration.
* Hovering a band highlights every band of the same fire and lists its events.
* pyroCb markers are labelled with `pcb_index`, which is what makes a 12-event fire readable.
""")
code(r'''
DATA_JS   = json.dumps({str(y): payloads[y] for y in YEARS})
YEARS_JS  = json.dumps(YEARS)
COLORS_JS = json.dumps({str(y): C.YEAR_COLORS[y] for y in YEARS})
RAMP_JS   = json.dumps(C.DOY_RAMP)

HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Boreal fires &amp; pyroCb __Y0__-__Y1__</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css">
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<style>
  html,body{margin:0;height:100%}
  #map{position:absolute;inset:0;background:#aad3df}
  .panel{background:#fff;padding:9px 11px;border-radius:6px;font:12px/1.4 system-ui,sans-serif;
         box-shadow:0 1px 5px rgba(0,0,0,.4);max-width:252px}
  .panel h4{margin:0 0 6px;font-size:12px}
  .panel .sub{color:#666;margin:2px 0 6px}
  .chips{display:flex;flex-wrap:wrap;gap:4px}
  .chip{border:1px solid #bbb;background:#fff;color:#333;border-radius:11px;
        padding:2px 8px;font:11px system-ui;cursor:pointer;user-select:none}
  .chip.on{color:#fff;border-color:transparent;font-weight:600}
  .mini{display:flex;gap:6px;margin:7px 0 4px}
  .mini button{flex:1;border:1px solid #bbb;background:#f7f7f7;border-radius:4px;
               padding:3px 0;font:11px system-ui;cursor:pointer}
  .layers{border-top:1px solid #e3e3e3;margin-top:8px;padding-top:7px}
  .layers label{display:flex;align-items:center;gap:6px;margin:3px 0;cursor:pointer}
  .seg{display:flex;gap:0;margin:4px 0 2px}
  .seg button{flex:1;border:1px solid #bbb;background:#f7f7f7;padding:3px 0;
              font:11px system-ui;cursor:pointer}
  .seg button:first-child{border-radius:4px 0 0 4px}
  .seg button:last-child{border-radius:0 4px 4px 0;border-left:none}
  .seg button.on{background:#333;color:#fff;border-color:#333}
  .legend .row{display:flex;align-items:center;gap:6px;margin:2px 0}
  .legend i{width:13px;height:13px;display:inline-block;border-radius:50%;flex:none;border:1px solid #666}
  .legend i.sq{border-radius:2px}
  #bar{height:11px;border-radius:2px;margin:3px 0 1px}
  #barlab{display:flex;justify-content:space-between;color:#666;font-size:10px}
  .pcblbl{font:700 10px system-ui;color:#fff;text-align:center;line-height:15px;
          text-shadow:0 0 2px #000,0 0 2px #000}
  .ignstar{position:relative}
  .ignstar svg{display:block}
  .ignstar span{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
     font:700 8px system-ui;color:#fff;text-shadow:0 0 2px #000,0 0 2px #000;pointer-events:none}
  .leaflet-popup-content{font:12px/1.45 system-ui,sans-serif;max-height:340px;overflow:auto}
  #busy{position:absolute;z-index:9999;top:12px;left:50%;transform:translateX(-50%);
        background:#222;color:#fff;padding:6px 14px;border-radius:14px;font:12px system-ui;display:none}
</style>
</head>
<body>
<div id="map"></div>
<div id="busy">Loading&hellip;</div>
<script>
const DATA=__DATA__, YEARS=__YEARS__, YCOL=__COLORS__, RAMP=__RAMP__;
const DOY0=__DOY0__, DOY1=__DOY1__, ORPHAN_KM=__ORPHAN__;

async function unpack(b64){
  const bin = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
  const s = new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"));
  return JSON.parse(await new Response(s).text());
}
const esc = s => s==null ? "" : String(s).replace(/[&<>]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
const rows = l => l.filter(r => r[1]!=null && r[1]!=="" && r[1]!==false)
  .map(r=>`<tr><td style="color:#666;padding-right:8px;vertical-align:top">${esc(r[0])}</td><td>${esc(r[1])}</td></tr>`).join("");
const busy = document.getElementById("busy");

// ---- colour ---------------------------------------------------------------
function hex2rgb(h){return [1,3,5].map(i=>parseInt(h.slice(i,i+2),16));}
function rgb2hex(r){return "#"+r.map(v=>Math.max(0,Math.min(255,Math.round(v))).toString(16).padStart(2,"0")).join("");}
function ramp(t){                       // t in [0,1] across the RAMP stops
  t = Math.max(0, Math.min(1, t));
  const x = t*(RAMP.length-1), i = Math.min(RAMP.length-2, Math.floor(x)), f = x-i;
  const a = hex2rgb(RAMP[i]), b = hex2rgb(RAMP[i+1]);
  return rgb2hex(a.map((v,k)=>v+(b[k]-v)*f));
}
const doyColor = d => ramp(((d==null?DOY0:d)-DOY0)/Math.max(1,DOY1-DOY0));
function lighten(hex, L){               // push a colour toward white / black
  const c = hex2rgb(hex);
  return L>=0 ? rgb2hex(c.map(v=>v+(255-v)*L)) : rgb2hex(c.map(v=>v*(1+L)));
}
let MODE = "year";
const hueOf = (y,f) => MODE==="year" ? YCOL[y] : doyColor(f.properties.doy);

// 5-point star, black outline, hue fill, ignition number centred on it.
const STAR = "M10,1 L12.12,7.09 L18.56,7.22 L13.42,11.11 L15.29,17.28 "
           + "L10,13.6 L4.71,17.28 L6.58,11.11 L1.44,7.22 L7.88,7.09 Z";
function starIcon(hue, n){
  return L.divIcon({className:"", iconSize:[20,20], iconAnchor:[10,10], html:
    `<div class="ignstar"><svg width="20" height="20" viewBox="0 0 20 20">`
    + `<path d="${STAR}" fill="${hue}" stroke="#111" stroke-width="1.4"`
    + ` stroke-linejoin="round"/></svg><span>${n==null?"":n}</span></div>`});
}

// ---- map ------------------------------------------------------------------
const map = L.map("map").setView([62,-105], 3);
L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
  {maxZoom:19, attribution:"&copy; OpenStreetMap"}).addTo(map);

const PANES = ["per","act","bnd","ign","pcb","orp"];
const R = {};
PANES.forEach((k,i)=>{ map.createPane(k); map.getPane(k).style.zIndex = 410+i*10;
                       R[k] = L.canvas({pane:k}); });
map.createPane("ignlbl"); map.getPane("ignlbl").style.zIndex = 445;   // stars: above bands, below pyroCb
map.createPane("lbl"); map.getPane("lbl").style.zIndex = 470;

const ALL = ["per","perL","act","bnd","ign","ignL","pcb","pcbC","orp"];
const cache = {}, byFire = {};
const activeYears = new Set([YEARS[YEARS.length-1]]);
const showKind = {per:true, act:true, bnd:true, ign:true, pcb:true, orp:true};
let onlyPcbFires=false, hideOrphans=false;

function buildLayers(y, gj){
  const perOpts = f => ({
    pane:"per", renderer:R.per, filter:f,
    style: ft => ({color:hueOf(y,ft), weight:1, fillColor:hueOf(y,ft),
                   fillOpacity:(ft.properties.farea_km2||0)>100?0.55:0.3}),
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>Fire #${esc(p.mergid)}</b> ${y}<table>${rows([
        ["Area (km²)",p.farea_km2],["Duration (d)",p.duration_days],
        ["VIIRS pixels",p.n_pixels],["Land cover",p.lcc_final],["Day of year",p.doy],
        ["Linked to a pyroCb",p.pcb_linked?"yes":"no"]])}</table>`);}
  });
  const ignOpts = f => ({
    pane:"ignlbl", filter:f,
    pointToLayer:(ft,ll)=>L.marker(ll,{pane:"ignlbl", icon:starIcon(hueOf(y,ft),
      ft.properties.ign_index)}),
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>Ignition ${esc(p.ign_index)} of ${esc(p.n_ign_total)}</b>`
        +` &middot; fire #${esc(p.mergid)} &middot; ${y}<table>${rows([
        ["Order within fire",p.ign_index+" of "+p.n_ign_total+" (by ignition time)"],
        ["Day of year",p.doy],["Fuel",p.ftype],["Area (km²)",p.farea_km2],
        ["VIIRS pixels",p.n_pixels]])}</table>`);}
  });
  const pcbOpts = f => ({
    pane:"pcb", renderer:R.pcb, filter:f,
    pointToLayer:(ft,ll)=>{
      const p=ft.properties, n=p.pcb_index;
      return n ? L.marker(ll,{pane:"lbl", icon:L.divIcon({className:"", iconSize:[15,15],
               html:`<div class="pcblbl" style="background:${hueOf(y,ft)};border:1.5px solid #111;
                     border-radius:50%;width:15px;height:15px">${n}</div>`})})
               : L.circleMarker(ll,{pane:"pcb",renderer:R.pcb,radius:6,color:"#111",
                   weight:1.5,fillColor:hueOf(y,ft),fillOpacity:1});},
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>${esc(p.pcb_id)}</b><br>${esc(p.fire_name||"unnamed")}<table>${rows([
        ["Date/time (UTC)",p.datetime],["Region",p.region],["Day of year",p.doy],
        ["Max injection alt (m)",p.max_inject_alt],["Rel. tropopause (m)",p.trop_relative_alt],
        ["Min BT (°C)",p.pcb_min_bt],["Match",p.match_method],["Dist. to fire (km)",p.match_dist_km],
        ["Fire (mergid)",p.mergid],["pyroCb # on this fire",p.pcb_index+" of "+p.n_pcb_total],
        ["pyroCb before / after",p.n_pcb_before+" / "+p.n_pcb_after],
        ["Ignitions in complex",p.n_ign_total],
        ["Area before (km²)",p.area_pre_km2],["Area after (km²)",p.area_post_km2],
        ["Fraction burned before",p.frac_burned_before],
        ["Obs lag (h)",p.obs_lag_h]])}</table>`);}
  });
  const bndOpts = {
    pane:"bnd", renderer:R.bnd,
    style: ft => {
      const p=ft.properties, n=Math.max(1,p.n_bands-1);
      const t = n>0 ? p.band_i/n : 0.5;                 // 0 = core, 1 = outermost
      const hue = MODE==="year" ? YCOL[y] : doyColor(p.doy);
      // band order in LIGHTNESS so hue stays free for the colour mode
      const fill = MODE==="year" ? lighten(hue, 0.72-1.05*t) : hue;
      return {color: MODE==="year" ? "#fff" : "#fff",
              weight: MODE==="year" ? 0.8 : 0.8+1.4*t,
              fillColor: fill,
              fillOpacity: p.is_post_final_band ? 0.32 : (MODE==="year"?0.8:0.55+0.3*t)};
    },
    onEachFeature:(ft,l)=>{const p=ft.properties;
      (byFire[p.mergid] = byFire[p.mergid] || []).push(l);
      l.bindPopup(`<b>Fire #${esc(p.mergid)} — band ${p.band_i+1} of ${p.n_bands}</b><table>${rows([
        ["Burned between",(p.step_hi==="final"?"last pyroCb → season end":"start → "+p.step_hi)],
        ["Area (km²)",p.area_km2],["Fire total (km²)",p.cum_area_km2],
        ["After last pyroCb",p.is_post_final_band?"yes (unknown duration)":""],
        ["Snapshot regression (km²)",p.qc_regress_km2]])}</table>`);
      l.on("mouseover",()=>{(byFire[p.mergid]||[]).forEach(x=>x.setStyle({weight:2.4,color:"#111"}));});
      l.on("mouseout", ()=>{(byFire[p.mergid]||[]).forEach(x=>bndLayerStyle(x));});}
  };
  const bndLayerStyle = l => l.setStyle(bndOpts.style(l.feature));

  // Active burn area: union of VIIRS pixel footprints within +/-12 h of onset. This is
  // the polygon the land-cover and ERA5 covariates get sampled over, so it is drawn as a
  // bright outline on top rather than a fill - you need to see what it covers.
  const actOpts = {
    pane:"act", renderer:R.act,
    style: ft => ({color:"#00e5ff", weight:1.6, fill:true, fillColor:"#00e5ff",
                   fillOpacity:0.18, dashArray:"4 3"}),
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>Active burn at onset &plusmn;${p.window_h} h</b><br>${esc(p.pcb_id)}<table>${rows([
        ["Fire (mergid)",p.mergid],["Onset (UTC)",p.datetime],
        ["Active area (km²)",p.area_km2],["Raw pixel footprint (km²)",p.footprint_km2],
        ["VIIRS detections",p.n_pixels],["Distinct overpasses",p.n_obs],
        ["Detections before / after onset",p.n_before_onset+" / "+p.n_after_onset],
        ["FRP p95 (MW)",p.frp_p95],["First / last obs",p.t_first+" → "+p.t_last]])}</table>`);}
  };

  return {
    per:  L.geoJSON(gj.per, perOpts(undefined)),
    perL: L.geoJSON(gj.per, perOpts(f=>f.properties.pcb_linked===true)),
    bnd:  L.geoJSON(gj.bnd, bndOpts),
    act:  L.geoJSON(gj.act, actOpts),
    ign:  L.geoJSON(gj.ign, ignOpts(undefined)),
    ignL: L.geoJSON(gj.ign, ignOpts(f=>f.properties.pcb_linked===true)),
    pcb:  L.geoJSON(gj.pcb, pcbOpts(undefined)),
    pcbC: L.geoJSON(gj.pcb, pcbOpts(f=>f.properties.match_method!=="unmatched")),
    orp:  L.geoJSON(gj.pcb, {pane:"orp", renderer:R.orp,
            filter:f=>f.properties.match_method==="unmatched",
            pointToLayer:(f,ll)=>L.circleMarker(ll,{pane:"orp",renderer:R.orp,radius:11,
              color:"#e11",weight:2,fill:false}),
            onEachFeature:(f,l)=>{const p=f.properties;
              l.bindPopup(`<b>No ${y} fire within ${ORPHAN_KM} km</b><table>${rows([
                ["Event",p.pcb_id],["Fire name",p.fire_name],["Region",p.region],
                ["Nearest (km)",p.match_dist_km]])}</table>`);}}),
  };
}

function wanted(){
  const w = new Set();
  if (showKind.per) w.add(onlyPcbFires?"perL":"per");
  if (showKind.bnd) w.add("bnd");
  if (showKind.act) w.add("act");
  if (showKind.ign) w.add(onlyPcbFires?"ignL":"ign");
  if (showKind.pcb) w.add(hideOrphans?"pcbC":"pcb");
  if (showKind.orp && !hideOrphans) w.add("orp");
  return w;
}

async function ensureYear(y){
  if (cache[y]) return cache[y];
  const d = DATA[y];
  const [per,ign,pcb,bnd,act] = await Promise.all([d.per,d.ign,d.pcb,d.bnd,d.act].map(unpack));
  cache[y] = buildLayers(y,{per,ign,pcb,bnd,act});
  return cache[y];
}

let refreshing=false, pending=false;
async function refresh(){
  if (refreshing){ pending=true; return; }
  refreshing=true;
  do {
    pending=false;
    const need = YEARS.filter(y=>activeYears.has(y) && !cache[y]);
    if (need.length){ busy.textContent=`Loading ${need.join(", ")}…`; busy.style.display="block"; }
    const w = wanted();
    for (const y of YEARS){
      if (activeYears.has(y)){
        const c = await ensureYear(y);
        for (const k of ALL){
          if (w.has(k)){ if(!map.hasLayer(c[k])) c[k].addTo(map); }
          else if (map.hasLayer(c[k])) map.removeLayer(c[k]);
        }
      } else if (cache[y]){
        for (const k of ALL) if (map.hasLayer(cache[y][k])) map.removeLayer(cache[y][k]);
      }
    }
    busy.style.display="none";
    document.querySelectorAll(".chip[data-year]").forEach(c=>{
      const on = activeYears.has(+c.dataset.year);
      c.classList.toggle("on", on);
      c.style.background = on ? YCOL[c.dataset.year] : "#fff";
    });
  } while (pending);
  refreshing=false;
}

function recolour(){
  for (const y of Object.keys(cache)){
    const c = cache[y];
    ["per","perL"].forEach(k=>c[k].setStyle(ft=>({color:hueOf(y,ft),fillColor:hueOf(y,ft),
      weight:1, fillOpacity:(ft.properties.farea_km2||0)>100?0.55:0.3})));
    ["ign","ignL"].forEach(k=>c[k].eachLayer(l=>{
      const el=l.getElement&&l.getElement(), pth=el&&el.querySelector("path");
      if (pth) pth.setAttribute("fill", hueOf(y,l.feature));
    }));
    ["pcb","pcbC"].forEach(k=>c[k].eachLayer(l=>{
      if (l.setStyle) l.setStyle({fillColor:hueOf(y,l.feature)});
      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = hueOf(y,l.feature);
    }));
    c.bnd.eachLayer(l=>{
      const p=l.feature.properties, n=Math.max(1,p.n_bands-1), t=n>0?p.band_i/n:0.5;
      const hue = MODE==="year" ? YCOL[y] : doyColor(p.doy);
      l.setStyle({color:"#fff", weight: MODE==="year"?0.8:0.8+1.4*t,
        fillColor: MODE==="year" ? lighten(hue,0.72-1.05*t) : hue,
        fillOpacity: p.is_post_final_band?0.32:(MODE==="year"?0.8:0.55+0.3*t)});
    });
  }
  drawLegend();
}

// ---- controls -------------------------------------------------------------
const panel = L.control({position:"topright"});
panel.onAdd = () => {
  const d = L.DomUtil.create("div","panel");
  d.innerHTML = `<h4>Years</h4><div class="sub">Click a year to show / hide it.</div>
    <div class="chips">${YEARS.map(y=>`<button class="chip" data-year="${y}">${y}</button>`).join("")}</div>
    <div class="mini"><button id="allY">All</button><button id="noneY">None</button></div>
    <div class="layers"><b style="font-size:11px">Colour by</b>
      <div class="seg"><button id="mYear" class="on">Year</button><button id="mDoy">Day of year</button></div>
      <div class="sub">Pick one year to watch a season progress.</div></div>
    <div class="layers">
      <label><input type="checkbox" data-kind="per" checked> Fire perimeters</label>
      <label><input type="checkbox" data-kind="bnd" checked> Growth bands (split at pyroCb)</label>
      <label><input type="checkbox" data-kind="act" checked> Active burn area at onset (&plusmn;12 h)</label>
      <label><input type="checkbox" data-kind="ign" checked> Fire ignitions</label>
      <label><input type="checkbox" data-kind="pcb" checked> pyroCb events</label>
      <label><input type="checkbox" data-kind="orp" checked> Flag unmatched pyroCb</label>
    </div>
    <div class="layers">
      <label><input type="checkbox" id="onlyPcb"> Only fires with a pyroCb</label>
      <label><input type="checkbox" id="hideOrp"> Hide unmatched pyroCb</label>
    </div>`;
  L.DomEvent.disableClickPropagation(d); L.DomEvent.disableScrollPropagation(d);
  d.querySelectorAll(".chip").forEach(c=>c.onclick=()=>{
    const y=+c.dataset.year; activeYears.has(y)?activeYears.delete(y):activeYears.add(y); refresh();});
  d.querySelector("#allY").onclick  = ()=>{YEARS.forEach(y=>activeYears.add(y)); refresh();};
  d.querySelector("#noneY").onclick = ()=>{activeYears.clear(); refresh();};
  d.querySelectorAll("input[data-kind]").forEach(i=>i.onchange=()=>{
    showKind[i.dataset.kind]=i.checked; refresh();});
  d.querySelector("#onlyPcb").onchange = e=>{onlyPcbFires=e.target.checked; refresh();};
  d.querySelector("#hideOrp").onchange = e=>{hideOrphans=e.target.checked; refresh();};
  const my=d.querySelector("#mYear"), mdy=d.querySelector("#mDoy");
  my.onclick  = ()=>{MODE="year"; my.classList.add("on"); mdy.classList.remove("on"); recolour();};
  mdy.onclick = ()=>{MODE="doy";  mdy.classList.add("on"); my.classList.remove("on"); recolour();};
  return d;
};
panel.addTo(map);

const legend = L.control({position:"bottomright"});
let legendDiv;
legend.onAdd = () => { legendDiv = L.DomUtil.create("div","panel legend"); drawLegend(); return legendDiv; };
legend.addTo(map);

function drawLegend(){
  if (!legendDiv) return;
  const common = `<div class="row"><i class="sq" style="background:#888"></i>Fire perimeter</div>
    <div class="row"><svg width="15" height="15" viewBox="0 0 20 20" style="flex:none">
      <path d="${STAR}" fill="#888" stroke="#111" stroke-width="1.4" stroke-linejoin="round"/></svg>
      Ignition (labelled with its order within the fire)</div>
    <div class="row"><i style="background:#888;border:1.5px solid #111"></i>pyroCb (labelled with its # on the fire)</div>
    <div class="row"><i style="border:2px solid #e11;background:none"></i>unmatched pyroCb</div>
    <div class="row"><i class="sq" style="border:2px solid #00e5ff;background:rgba(0,229,255,.18)"></i>active burn at onset (&plusmn;12 h)</div>
    <div class="sub" style="margin-top:5px">Growth bands: pale core = burned before the first
    pyroCb, darker = later runs; the faintest outer band is after the last pyroCb.</div>`;
  if (MODE==="year"){
    legendDiv.innerHTML = `<h4>Boreal fires &amp; pyroCb</h4>${common}
      <div class="sub" style="margin:6px 0 3px">Colour = year</div>
      <div class="chips">${YEARS.map(y=>
        `<span class="chip on" style="background:${YCOL[y]};cursor:default">${y}</span>`).join("")}</div>`;
  } else {
    const stops = Array.from({length:24},(_,i)=>ramp(i/23)).join(",");
    const mid = Math.round((DOY0+DOY1)/2);
    const md = d => { const t=new Date(2001,0,1); t.setDate(d);
      return t.toLocaleString("en",{month:"short"}); };
    legendDiv.innerHTML = `<h4>Boreal fires &amp; pyroCb</h4>${common}
      <div class="sub" style="margin:6px 0 3px">Colour = day of year</div>
      <div id="bar" style="background:linear-gradient(90deg,${stops})"></div>
      <div id="barlab"><span>${md(DOY0)} (${DOY0})</span><span>${md(mid)}</span><span>${md(DOY1)} (${DOY1})</span></div>`;
  }
}

refresh();
</script>
</body>
</html>
"""

HTML = (HTML.replace("__DATA__", DATA_JS).replace("__YEARS__", YEARS_JS)
            .replace("__COLORS__", COLORS_JS).replace("__RAMP__", RAMP_JS)
            .replace("__DOY0__", str(DOY_RANGE[0])).replace("__DOY1__", str(DOY_RANGE[1]))
            .replace("__ORPHAN__", str(int(C.MATCH_CAP_KM)))
            .replace("__Y0__", str(YEARS[0])).replace("__Y1__", str(YEARS[-1])))
C.OUT_HTML.write_text(HTML, encoding="utf-8")
print(f"Wrote {C.OUT_HTML.name}  ({C.OUT_HTML.stat().st_size/1e6:.1f} MB)")
''')

md(r"""
## 8. Notes

* **2012 is absent by construction** — `YEARS = None` keeps only years present in the pyroCb
  inventory, and the inventory starts in 2013.
* **Ignitions are filtered on `mergid`, not `fireid`.** Ignition `fireid` is unique per row, so
  filtering on it keeps only the *primary* ignition of each merged fire and drops the other 1,444
  (2023). This was a real bug in the earlier version of this map.
* **Half-day AM/PM stamps are local-solar overpass groups**, not calendar halves, and there is a
  ~10 h nightly observation gap that contains the pyroCb peak. Every step-derived quantity
  inherits that resolution; `obs_lag_h` says how much.
* To run every year, set `C.SNAPSHOT_YEARS = None` and re-run — cached years are reused.
* `C.clear_cache()` wipes the cache; `C.clear_cache("bands")` just one stage.
""")

nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}

# Guard: every code cell must parse. A `\"\"\"` docstring written inside a
# code(r\"\"\"...\"\"\") block silently terminates the cell early and truncates the
# notebook - it has bitten us twice (ymd, smooth_for_map). Catch it here, at
# generation time, rather than discovering it as a stale notebook three runs later.
import ast as _ast
for _i, _c in enumerate(cells):
    if _c["cell_type"] != "code":
        continue
    try:
        _ast.parse(_c["source"])
    except SyntaxError as _e:
        raise SystemExit(
            f"code cell {_i} does not parse: {_e}\n"
            f"  (a triple-quoted string inside a code(r\"\"\"...\"\"\") block ends the cell -\n"
            f"   use # comments instead)\n---\n{_c['source'][:400]}") from None

nbf.write(nb, "Data exploration.ipynb")
print("wrote Data exploration.ipynb  (%d cells, all parse)" % len(cells))
