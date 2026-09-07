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
                    variables, pipeline, landcover, firesim)

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
    lf = pipeline.landcover_fractions_cached(y, active_geom[y], e)
    if len(lf) and lf.shape[1] > 1:
        e = e.merge(lf, on="pcb_id", how="left")
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

payloads, doys, spread_by_year = {}, [], {}
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

    # "Blew up before it was ever observed": the +/-12h active-burn window has
    # detections, but every one of them comes from the overpass AFTER onset - the fire
    # existed and was already burning, it just wasn't seen until the next look (F4).
    if "active_n_pixels_pre12_post12" in pg.columns:
        pg["active_is_blowup"] = ((pg["active_n_pixels_pre12_post12"] > 0)
                                   & (pg["active_n_before_onset_pre12_post12"] == 0))

    # Slider / animation data: one localisation pass gives the 24-window area lookup
    # and the per-detection table the map filters by time.
    agrid, px = active.slider_data(y, ev_y) if len(ev_y) else ({}, {})
    if agrid:
        pg["agrid"] = pg["pcb_id"].map(agrid)

    # POSITION ACCURACY. The atlas stores lat/lon to 0.1 deg for 93% of events - about
    # 11 km of latitude, 5 km of longitude at 60N - so the plotted pyroCb dot is only
    # ever a neighbourhood. The matched fire, on the other hand, is known to one VIIRS
    # pixel. Ship BOTH: `loc_*` is the reported point, `fire_*` the centroid of the
    # 6 h-pre-onset active burn, and `loc_rounded` says whether the report was snapped.
    ag_y = active_geom.get(y)
    if ag_y is not None and len(ag_y):
        lo_, hi_ = C.LC_ACTIVE_WINDOW
        cen = ag_y[(ag_y.pre_h == lo_) & (ag_y.post_h == hi_)].copy()
        cen = cen[~cen.geometry.is_empty & cen.geometry.notna()]
        if len(cen):
            c4326 = gpd.GeoSeries(cen.geometry.values, crs=C.METRIC_CRS).centroid.to_crs(C.WGS84)
            cmap = dict(zip(cen["pcb_id"], zip(c4326.y.round(5), c4326.x.round(5))))
            pg["fire_lat"] = pg["pcb_id"].map(lambda k: (cmap.get(k) or (None, None))[0])
            pg["fire_lon"] = pg["pcb_id"].map(lambda k: (cmap.get(k) or (None, None))[1])
    _r = lambda v: bool(abs(v * 10 - round(v * 10)) < 1e-9) if pd.notna(v) else False
    pg["loc_rounded"] = [(_r(a) and _r(b)) for a, b in zip(p["lat"], p["lon"])]

    pg = pg[[c for c in ("pcb_id", "datetime", "doy", "fire_name", "region", "status",
                         "max_inject_alt", "trop_relative_alt", "pcb_min_bt", "mergid",
                         "match_method", "match_dist_km", "pcb_index", "n_pcb_total",
                         "n_pcb_before", "n_pcb_after", "n_ign_total", "area_pre_km2",
                         "area_post_km2", "frac_burned_before", "obs_lag_h",
                         "active_area_km2_pre12_post12", "active_n_pixels_pre12_post12",
                         "active_n_obs_pre12_post12", "active_n_before_onset_pre12_post12",
                         "active_n_after_onset_pre12_post12", "active_frp_p95_pre12_post12",
                         "active_detach_km_pre12_post12", "active_is_blowup", "agrid",
                         "lc_top1_cls", "lc_top1_frac", "lc_top2_cls", "lc_top2_frac",
                         "lc_top3_cls", "lc_top3_frac", "lc_n_px",
                         "fire_lat", "fire_lon", "loc_rounded",
                         "geometry")
                if c in pg.columns]]

    # Fires credited with a pyroCb, for the "only fires with a pyroCb" filter.
    linked = set(ev_y.loc[ev_y.match_method != "unmatched", "mergid"].dropna().astype(int)) \
        if len(ev_y) else set()
    per_out["pcb_linked"] = per_out["mergid"].isin(linked)
    ign_out["pcb_linked"] = ign_out["mergid"].isin(linked)   # mergid, NOT fireid

    # Country/region per FIRE, inherited from the pyroCb events matched to it, so the
    # region toggles filter perimeters and ignitions as well as the events themselves.
    if len(ev_y):
        rmap = (ev_y.dropna(subset=["mergid"])
                    .assign(mergid=lambda d: d["mergid"].astype(int))
                    .groupby("mergid")["region"].first().to_dict())
    else:
        rmap = {}
    per_out["region"] = per_out["mergid"].map(rmap)
    ign_out["region"] = ign_out["mergid"].map(rmap)

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
        # Simplify only - NO buffer-close. The close was for the old active-burn pixel
        # lattice; bands are snapshot alpha shapes already, and +/-60 m of buffering
        # was rounding real concave detail off their boundaries.
        bg["geometry"] = bg.geometry.simplify(C.BAND_SIMPLIFY_M, preserve_topology=True)
        bg = gpd.GeoDataFrame(bg, geometry="geometry", crs=C.METRIC_CRS).to_crs(C.WGS84)
        bg = bg[~bg.geometry.is_empty & bg.geometry.notna()]
        bg["region"] = bg["mergid"].map(rmap)
        bg["doy"] = pd.to_numeric(bg["step_hi"].str[4:8].str.slice(0, 2), errors="coerce")
        bg["doy"] = pd.to_datetime(bg["step_hi"].str[:8], format="%Y%m%d",
                                   errors="coerce").dt.dayofyear
        b64["bnd"], _, _ = pack(json.loads(bg.drop(columns=[c for c in ("step_lo",) if c in bg]).to_json(drop_id=True)))
    else:
        b64["bnd"], _, _ = pack({"type": "FeatureCollection", "features": []})

    # +/-24 h cell grids (localised, 375 m, delta-coded) -> the slider / animation.
    b64["pix"], _, gpx = pack(px)

    # Land cover over each pyroCb fire + 5 km, decimated to 1 km, dissolved by class.
    lc = pipeline.landcover_cached(y, linked) if linked else None
    if lc is not None and len(lc):
        lc = lc.copy()
        # Tier-aware reduction for the browser: the sub-1 km tiers carry far more
        # vertices than the map can show, so simplify them ~75 m and snap ~22 m.
        fine = lc["res"] < 1000
        lc.loc[fine, "geometry"] = lc.loc[fine, "geometry"].simplify(7e-4, preserve_topology=True)
        lc["geometry"] = shapely.set_precision(lc.geometry.values, 2e-4)
        lc = lc[~lc.geometry.is_empty & lc.geometry.notna()]
        b64["lc"], _, glc = pack(json.loads(lc.to_json(drop_id=True)))
    else:
        b64["lc"], _, glc = pack({"type": "FeatureCollection", "features": []})

    # Land cover close to each fire: paletted PNG overlays (see pyrocb/lcraster.py).
    # Left un-gzipped: the payload is already base64 PNG, so packing it would just add
    # a second base64 pass over incompressible bytes.
    lcpng = pipeline.lcpng_cached(y, linked) if linked else {}
    gpng = sum(len(v["d"]) for v in lcpng.values())
    b64["lcpng"] = lcpng

    # Whole-fire evolution: 375 m cell grid, first-burn time per cell.
    sim = pipeline.firesim_cached(y, events, linked) if (len(ev_y) and linked) else {}
    b64["sim"], _, gsim = pack(sim)

    # 12-hourly spread-rate series per fire, for the animation's colour modes.
    spr = pipeline.spread_steps_cached(y, linked) if linked else {}
    spread_by_year[str(y)] = spr
    b64["spr"], _, gspr = pack(spr)

    payloads[y] = b64
    _nsim = sum(v["n_cell"] for v in sim.values())
    _npix = sum(len(v["di"]) for v in px.values())
    print(f"  {y}: {len(per_out):>6} perims  {len(ign_out):>6} ign  {len(pg):>3} pyroCb"
          f"  {0 if bg is None else len(bg):>3} bands"
          f" | pix {_npix:>7} cells {gpx/1e6:4.2f}MB"
          f" | sim {len(sim):>3} fires {_nsim:>7} cells {gsim/1e6:4.2f}MB"
          f" | spr {gspr/1e6:4.2f}MB"
          f" | lc {glc/1e6:4.2f}MB poly + {gpng/1e6:4.2f}MB png({len(lcpng)})",
          flush=True)

from pyrocb import spread as _SP
SPREAD_CLAMP = _SP.clamp(spread_by_year)
print(f"\nspread ramp clamp (p2..p98 non-zero): "
      f"avgspread {SPREAD_CLAMP['av'][0]:.0f}..{SPREAD_CLAMP['av'][1]:.0f} m/12h, "
      f"spread95 {SPREAD_CLAMP['s95'][0]:.0f}..{SPREAD_CLAMP['s95'][1]:.0f} m/12h")

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
SPREAD_CLAMP_JS = json.dumps(SPREAD_CLAMP)

# Month-by-month pyroCb counts per (year, region), so the animate panel's month ticks
# and its running total follow whichever countries are toggled on.
_mc = (pcb_all.assign(m=pcb_all.datetime.dt.month)
       .groupby(["year", "region", "m"]).size())
MONTH_COUNTS_JS = json.dumps({
    str(y): {r: {int(m): int(_mc.get((y, r, m), 0)) for m in range(3, 11)}
             for r in sorted(pcb_all.region.dropna().unique())}
    for y in YEARS})
REGIONS_JS = json.dumps(sorted(pcb_all.region.dropna().unique().tolist()))
print("regions:", pcb_all.region.value_counts().to_dict())

# Tiny global index for the "find a fire" search box: one row per event, every year,
# so a search hits regardless of which years are currently loaded.
_si = (events if len(events) else pcb_all).copy()
_si_cols = {"pcb_id": "id", "fire_name": "nm", "mergid": "mg", "year": "y",
            "lat": "la", "lon": "lo", "datetime": "dt"}
_si = _si[[c for c in _si_cols if c in _si.columns]].rename(columns=_si_cols)
if "dt" in _si:
    _si["dt"] = pd.to_datetime(_si["dt"]).dt.strftime("%Y-%m-%d")
_si["nm"] = _si["nm"].astype("string").fillna("")
_si["mg"] = pd.to_numeric(_si.get("mg"), errors="coerce")
SEARCH_JS = json.dumps([{k: (None if (isinstance(v, float) and pd.isna(v)) else v)
                         for k, v in r.items()}
                        for r in _si.round({"la": 4, "lo": 4}).to_dict("records")])

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
         box-shadow:0 1px 5px rgba(0,0,0,.4);max-width:252px;
         max-height:calc(100vh - 28px);overflow-y:auto;overflow-x:hidden}
  .leaflet-right .panel{clear:both}
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
  details.lcdd{margin-top:6px;border:1px solid #ddd;border-radius:4px}
  details.lcdd>summary{cursor:pointer;padding:4px 7px;color:#444;font-weight:600;
    list-style:none;user-select:none}
  details.lcdd>summary::-webkit-details-marker{display:none}
  details.lcdd>summary::before{content:"\25B8";color:#888;margin-right:5px}
  details.lcdd[open]>summary::before{content:"\25BE"}
  details.lcdd .body{padding:3px 7px 5px;max-height:190px;overflow:auto;border-top:1px solid #eee}
  #bar{height:11px;border-radius:2px;margin:3px 0 1px}
  #barlab{display:flex;justify-content:space-between;color:#666;font-size:10px}
  .pcblbl{font:700 10px system-ui;color:#fff;text-align:center;line-height:15px;
          text-shadow:0 0 2px #000,0 0 2px #000}
  /* categorical raster: never interpolate between class colours */
  .lcimg{image-rendering:pixelated;image-rendering:crisp-edges}
  .ignstar{position:relative}
  .ignstar svg{display:block}
  .ignstar span{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
     font:700 8px system-ui;color:#fff;text-shadow:0 0 2px #000,0 0 2px #000;pointer-events:none}
  .leaflet-popup-content{font:12px/1.45 system-ui,sans-serif;max-height:340px;overflow:auto}
  #busy{position:absolute;z-index:9999;top:48px;left:50%;transform:translateX(-50%);
        background:#222;color:#fff;padding:6px 14px;border-radius:14px;font:12px system-ui;display:none}
  /* find-a-fire search */
  #search{position:absolute;z-index:9997;top:10px;left:50%;transform:translateX(-50%);
        width:min(360px,70vw);font:12px system-ui}
  #search input{width:100%;box-sizing:border-box;padding:6px 10px;border:1px solid #999;
        border-radius:16px;box-shadow:0 1px 5px rgba(0,0,0,.35);outline:none}
  #search ul{list-style:none;margin:3px 0 0;padding:0;background:#fff;border:1px solid #ccc;
        border-radius:8px;box-shadow:0 2px 12px rgba(0,0,0,.3);max-height:46vh;overflow:auto;display:none}
  #search li{padding:5px 10px;cursor:pointer;border-bottom:1px solid #f0f0f0}
  #search li:last-child{border-bottom:none}
  #search li:hover,#search li.sel{background:#eef3ff}
  #search li b{font-weight:700}
  #search li span{color:#666}
  /* global playback bar */
  /* animate section, inside the control panel (not floating over the map) */
  .anim .grow{display:flex;align-items:center;gap:7px;margin:5px 0 2px}
  .anim #gplay{border:1px solid #bbb;background:#f7f7f7;border-radius:5px;min-width:26px;
        padding:2px 6px;cursor:pointer;font-size:11px}
  .anim #glabel{font-weight:700;white-space:nowrap;font-size:11px}
  .anim #gsub{color:#666;white-space:nowrap;margin-left:auto;font-size:10px}
  /* labels ABOVE the track so the slider thumb never covers them */
  .anim #gticks{display:flex;justify-content:space-between;margin:2px 0 1px;
        font-size:9px;color:#666;line-height:1.15}
  .anim #gticks span{flex:1;text-align:center;white-space:nowrap}
  .anim #gticks span em{font-style:normal;display:block;color:#999;font-size:8.5px}
  .anim #gticks span.hot{color:#7d1128;font-weight:700}
  .anim #gticks span.hot em{color:#a5405a}
  .anim #gscrub{width:100%;margin:0;display:block}
  #notice{position:absolute;z-index:9999;top:48px;left:50%;transform:translateX(-50%);
        max-width:420px;background:#fff3cd;color:#664d03;border:1px solid #ffda6a;
        padding:8px 14px;border-radius:8px;font:12px/1.4 system-ui;display:none;
        box-shadow:0 2px 10px rgba(0,0,0,.25)}
  #notice.info{background:#e7f1ff;color:#084298;border-color:#9ec5fe}
  #notice b{display:block;margin-bottom:2px}
  #notice .x{float:right;cursor:pointer;font-weight:700;margin-left:8px}
  #info{position:absolute;z-index:9998;left:12px;bottom:12px;width:310px;max-height:70vh;
        overflow:auto;background:#fff;border:1px solid #ccc;border-radius:8px;
        box-shadow:0 2px 12px rgba(0,0,0,.3);font:12px/1.45 system-ui;display:none}
  #info .hd{background:#222;color:#fff;padding:7px 10px;border-radius:7px 7px 0 0;
        font-weight:700;position:relative}
  #info .hd .x{position:absolute;right:8px;top:6px;cursor:pointer;font-weight:700}
  #info .bd{padding:8px 10px}
  #info table{border-collapse:collapse;width:100%}
  #info td{padding:1px 6px 1px 0;vertical-align:top}
  #info .sl{margin:8px 0 4px;padding-top:6px;border-top:1px solid #eee}
  #info .sl label{display:flex;align-items:center;gap:6px;margin:3px 0}
  #info .sl input[type=range]{flex:1}
  #info .rd{font-weight:700;margin-top:4px}
  #info .warn{color:#8a5a00}
  #info .seg{display:inline-flex;gap:0;margin:2px 0 2px 6px;vertical-align:middle;flex-wrap:wrap}
  #info .seg button,#info #play{border:1px solid #bbb;background:#f7f7f7;
     padding:2px 7px;font:11px system-ui;cursor:pointer}
  #info .seg button:first-child{border-radius:4px 0 0 4px}
  #info .seg button:last-child{border-radius:0 4px 4px 0}
  #info .seg button+button{border-left:none}
  #info .seg button.on{background:#333;color:#fff;border-color:#333}
  #info #play{border-radius:4px;min-width:32px;font-size:12px}
  #info .player{display:flex;align-items:center;gap:8px;margin:6px 0 1px}
  #info #cur{font-weight:700;white-space:nowrap}
  #info #curabs{color:#666;font-size:11px;white-space:nowrap}
  #info #scrubwrap{position:relative;margin:7px 0 4px;padding-top:16px}
  #info #scrub{width:100%;margin:0;display:block;position:relative;z-index:1}
  #info #ticks{position:absolute;left:0;right:0;top:0;pointer-events:none}
  #info #ticks b{position:absolute;top:0;width:14px;height:14px;margin-left:-7px;
     background:#7d1128;color:#fff;font:700 8px/14px system-ui;text-align:center;
     border:1px solid #fff;border-radius:3px;cursor:pointer;pointer-events:auto;
     box-shadow:0 0 2px rgba(0,0,0,.5)}
  #info #ticks b.staggered::after{content:"";position:absolute;left:6px;top:13px;width:1px;
     background:rgba(125,17,40,.45);height:var(--stem,13px)}
  #info #ticks b:hover{background:#e11;z-index:3;transform:scale(1.25)}
  #info #loopCustom{width:46px;margin-left:6px;font:11px system-ui;padding:2px 4px;
     border:1px solid #bbb;border-radius:4px}
  @keyframes pcbpulse{0%{transform:scale(.3);opacity:.9}100%{transform:scale(3.4);opacity:0}}
  .pcbpulse{box-sizing:border-box;width:22px;height:22px;border:3px solid #7d1128;
     border-radius:50%;animation:pcbpulse 1.1s ease-out forwards}
</style>
</head>
<body>
<div id="map"></div>
<div id="search"><input type="text" placeholder="find a fire &mdash; pyroCb id, fire name, or merge id" autocomplete="off"><ul></ul></div>
<div id="busy">Loading&hellip;</div>
<div id="notice"></div>
<div id="info"></div>

<script>
const DATA=__DATA__, YEARS=__YEARS__, YCOL=__COLORS__, RAMP=__RAMP__;
const DOY0=__DOY0__, DOY1=__DOY1__, ORPHAN_KM=__ORPHAN__;
const SPREAD_CLAMP=__SPREADCLAMP__, MONTH_COUNTS=__MONTHCOUNTS__, SEARCH=__SEARCH__;
const REGIONS=__REGIONS__;
const MONTH_ABBR=["","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
// day-of-year for the 1st of a given month (non-leap; fires are Mar-Oct so leap day is moot)
const MONTH_DOY=[0,1,32,60,91,121,152,182,213,244,274,305,335,366];

async function unpack(b64){
  const bin = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
  const s = new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"));
  return JSON.parse(await new Response(s).text());
}
const esc = s => s==null ? "" : String(s).replace(/[&<>]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
const rows = l => l.filter(r => r[1]!=null && r[1]!=="" && r[1]!==false)
  .map(r=>`<tr><td style="color:#666;padding-right:8px;vertical-align:top">${esc(r[0])}</td><td>${esc(r[1])}</td></tr>`).join("");
const busy = document.getElementById("busy");
const noticeEl = document.getElementById("notice");
function showNotice(kind, html){
  noticeEl.className = kind;
  noticeEl.innerHTML = `<span class="x" title="dismiss">&times;</span>${html}`;
  noticeEl.style.display = "block";
  noticeEl.querySelector(".x").onclick = () => noticeEl.style.display="none";
}
function hideNotice(){ noticeEl.style.display = "none"; }

// pyroCb -> its per-detection pixels (widest window, localised to the fire that made
// the pyroCb), keyed by the globally unique pcb_id so a click works no matter which
// years are toggled on. Each entry: {o:[lat,lon], dlat:[], dlon:[], hw:[], hh:[], t:[]}.
const byPcbPix = {};
// mergid -> whole-fire series {t0, dur_h, pcb_h, pcb_ix, pcb_ll, pcb_id, g, n_raw, n_cell}.
const byFireSim = {};
// mergid -> 12-hourly spread series {t:[iso], av:[m/12h], s95:[m/12h]}, decoded lazily.
const bySpread = {};
const _spDec = new WeakMap();
function decodeSpread(s){
  if (!s || !s.t || !s.t.length) return null;
  let d = _spDec.get(s);
  if (d) return d;
  const n = s.t.length, ms = new Float64Array(n);
  for (let k = 0; k < n; k++) ms[k] = Date.parse(s.t[k] + "Z");
  d = {n, ms, av:s.av, s95:s.s95};
  _spDec.set(s, d);
  return d;
}

// ---- fire-activity grid ---------------------------------------------------
// Payloads are delta-coded, TIME-SORTED cell grids in Web Mercator (burngrid.py).
// Time-sorted is the whole trick: every view the map offers - cumulative up to t, a
// sliding window, a static -pre..+post window - is a CONTIGUOUS index range. So the
// canvas can draw the newly-revealed slice on top of what is already there instead of
// rebuilding thousands of Leaflet paths sixty times a second.
const R_MERC = 6378137.0, MERC_K = 1/(2*Math.PI*R_MERC);
const _decoded = new WeakMap();
function decodeGrid(p){
  if (!p) return null;
  let d = _decoded.get(p);
  if (d) return d;
  const n = p.di.length;
  const ix = new Int32Array(n), iy = new Int32Array(n), t = new Float32Array(n);
  let a = 0, b = 0, c = 0;
  for (let k = 0; k < n; k++){
    a += p.di[k]; b += p.dj[k]; c += p.dt[k];
    ix[k] = a; iy[k] = b; t[k] = c/60;              // minutes -> hours
  }
  d = {n, ix, iy, t, ox:p.ox, oy:p.oy, c:p.c};
  _decoded.set(p, d);
  return d;
}
// t is ascending: first index with t > v  /  first index with t >= v
function ubound(t, n, v){ let lo=0, hi=n; while(lo<hi){const m=(lo+hi)>>1; if(t[m]<=v) lo=m+1; else hi=m;} return lo; }
function lbound(t, n, v){ let lo=0, hi=n; while(lo<hi){const m=(lo+hi)>>1; if(t[m]< v) lo=m+1; else hi=m;} return lo; }

// Animation colour modes. "age" is a PURPLE ramp (old -> new) on purpose: red/yellow
// reads as fire *intensity*, which this is not. "av"/"s95" DO use a red->yellow ramp,
// there deliberately meaning speed - deep red = fast spread, yellow = slow - from the
// 12-hourly avgspread / spread95 snapshot attributes.
const AGE_STOPS    = ["#2c0b52","#43217a","#5c34a6","#7a51c7","#9c78dd","#c1a6ec","#ecdcfb"];
const SPREAD_STOPS = ["#ffe98a","#ffcf5c","#ffab33","#ff7a1f","#f24a12","#c81e0c","#7d0000"];
function _bandArray(stops, n){
  const out = [];
  for (let k = 0; k < n; k++){
    const x = k/(n-1)*(stops.length-1), i = Math.min(stops.length-2, x|0), f = x-i;
    const a = hex2rgb(stops[i]), b = hex2rgb(stops[i+1]);
    out.push(rgb2hex(a.map((v,q)=>v+(b[q]-v)*f)));
  }
  return out;
}
const AGE_BANDS = _bandArray(AGE_STOPS, 20), SPREAD_BANDS = _bandArray(SPREAD_STOPS, 20);
let BURN_COLOR_MODE = "age";           // "age" | "av" | "s95"
const BURN_MODE_LABEL = {age:"burn age", av:"avg spread", s95:"spread p95"};

// Subclasses L.Renderer rather than L.Layer so Leaflet's own container positioning
// and zoom-animation transform bookkeeping come for free - the same machinery L.Canvas
// uses. Like L.Canvas we translate the context by -bounds.min and then draw in LAYER
// coordinates, which is exactly where the mercator->pixel transform lands.
const BurnLayer = L.Renderer.extend({
  initialize(opts){
    L.Renderer.prototype.initialize.call(this, opts);
    this._g=null; this._range=[0,0]; this._drawn=null; this._span=1; this._t0=0;
    this._absT0=0; this._sp=null;         // grid t=0 in ms UTC; decoded spread series
  },
  _initContainer(){
    const c = this._container = document.createElement("canvas");
    c.style.pointerEvents = "none";
    this._ctx = c.getContext("2d");
  },
  _destroyContainer(){
    L.Util.cancelAnimFrame(this._redrawRequest);
    delete this._ctx;
    L.DomUtil.remove(this._container);
    L.DomEvent.off(this._container);
    delete this._container;
  },
  setGrid(grid, span, t0, absT0, spread){
    this._g = grid; this._span = Math.max(1e-6, span||1); this._t0 = t0||0;
    this._absT0 = absT0||0; this._sp = decodeSpread(spread);
    this._range = [0,0]; this._drawn = null; this._blank();
  },
  setContext(span, t0, absT0, spread){
    this._span = Math.max(1e-6, span||1); this._t0 = t0||0;
    this._absT0 = absT0||0; this._sp = decodeSpread(spread);
  },
  repaint(){ if (!this._g) return;
    const r = this._range; this._drawn = null; this._blank();
    this._paint(r[0], r[1]); this._drawn = [r[0], r[1]]; },
  clear(){ this._g = null; this._sp = null; this._range=[0,0]; this._drawn=null; this._blank(); },
  _onZoomEnd(){},                       // Renderer binds this; only L.Canvas defines it
  count(){ return this._range[1] - this._range[0]; },
  setRange(i0, i1){
    if (!this._g) return 0;
    i1 = Math.max(i0, i1);
    const d = this._drawn;
    if (d && d[0] === i0 && i1 >= d[1]) this._paint(d[1], i1);   // extend: draw only the new slice
    else { this._blank(); this._paint(i0, i1); }
    this._range = [i0, i1]; this._drawn = [i0, i1];
    return i1 - i0;
  },
  _blank(){
    if (!this._ctx || !this._bounds) return;
    const b = this._bounds, s = b.getSize();
    this._ctx.clearRect(b.min.x, b.min.y, s.x, s.y);
  },
  _update(){
    if (this._map._animatingZoom && this._bounds) return;
    L.Renderer.prototype._update.call(this);
    const b = this._bounds, size = b.getSize(), c = this._container;
    L.DomUtil.setPosition(c, b.min);
    c.width = size.x; c.height = size.y;
    c.style.width = size.x + "px"; c.style.height = size.y + "px";
    this._ctx.setTransform(1, 0, 0, 1, 0, 0);
    this._ctx.translate(-b.min.x, -b.min.y);      // draw in layer coordinates
    // Web Mercator metres -> layer point is a pure scale+offset, so a cell index
    // becomes a pixel with one multiply-add. No per-feature projection, ever.
    const scale = 256 * Math.pow(2, this._map.getZoom()), o = this._map.getPixelOrigin();
    this._A  = scale * MERC_K;
    this._Bx = 0.5*scale - o.x;
    this._By = 0.5*scale - o.y;
    const r = this._range;
    this._drawn = null; this._blank(); this._paint(r[0], r[1]); this._drawn = [r[0], r[1]];
  },
  _paint(i0, i1){
    const g = this._g;
    if (!g || !this._ctx || !this._bounds || i1 <= i0) return;
    const ctx = this._ctx, b = this._bounds;
    const S = g.c * this._A, PX = g.ox*this._A + this._Bx, PY = -g.oy*this._A + this._By;
    const w = Math.max(1, S + 0.6);
    const x0 = b.min.x - w, x1 = b.max.x, y0 = b.min.y - w, y1 = b.max.y;
    const ix = g.ix, iy = g.iy, t = g.t, span = this._span, t0 = this._t0;
    const mode = BURN_COLOR_MODE, sp = this._sp;
    const bands = (mode === "age" || !sp) ? AGE_BANDS : SPREAD_BANDS;
    const NB = bands.length;
    // spread modes: a forward pointer into the 12-hourly series (cells are time-sorted)
    const spOn = mode !== "age" && sp;
    const cl = spOn ? SPREAD_CLAMP[mode] : null;
    const spArr = spOn ? (mode === "av" ? sp.av : sp.s95) : null;
    let spi = 0;
    if (spOn){ const ms0 = this._absT0 + t[i0]*3600000;
      while (spi + 1 < sp.n && ms0 >= sp.ms[spi + 1]) spi++; }
    ctx.globalAlpha = 0.66;
    let band = -1;
    for (let k = i0; k < i1; k++){
      let tt;
      if (spOn){
        const ms = this._absT0 + t[k]*3600000;
        while (spi + 1 < sp.n && ms >= sp.ms[spi + 1]) spi++;
        tt = (spArr[spi] - cl[0]) / Math.max(1e-6, cl[1] - cl[0]);
      } else {
        tt = (t[k]-t0)/span;
      }
      const bnd = Math.min(NB-1, Math.max(0, (tt*NB)|0));
      if (bnd !== band){ band = bnd; ctx.fillStyle = bands[bnd]; }
      const x = ix[k]*S + PX, y = PY - (iy[k]+1)*S;
      if (x < x0 || y < y0 || x > x1 || y > y1) continue;
      ctx.fillRect(x, y, w, w);
    }
    ctx.globalAlpha = 1;
  },
  gridBounds(){
    const g = this._g;
    if (!g || !g.n) return null;
    let a=g.ix[0],A=a,b=g.iy[0],B=b;
    for (let k=1;k<g.n;k++){
      if(g.ix[k]<a)a=g.ix[k]; if(g.ix[k]>A)A=g.ix[k];
      if(g.iy[k]<b)b=g.iy[k]; if(g.iy[k]>B)B=g.iy[k];
    }
    const ll = (i,j) => { const mx=g.ox+i*g.c, my=g.oy+j*g.c;
      return [ (2*Math.atan(Math.exp(my/R_MERC)) - Math.PI/2) * 180/Math.PI,
               mx/R_MERC * 180/Math.PI ]; };
    return L.latLngBounds(ll(a,b), ll(A+1,B+1));
  },
});
// Growth bands straight off the 375 m burn grid: every cell coloured by the day it
// first burned. Needs no new payload (it reuses byFireSim) and is finer than the
// polygon bands, which only split at pyroCb onsets.
const SnapBandLayer = BurnLayer.extend({
  initialize(opts){ BurnLayer.prototype.initialize.call(this, opts); this._gs = []; },
  setGrids(list){ this._gs = list || []; this._blank(); this._paintAll(); },
  setRange(){ return 0; },
  _update(){
    if (this._map._animatingZoom && this._bounds) return;
    L.Renderer.prototype._update.call(this);
    const b = this._bounds, size = b.getSize(), c = this._container;
    L.DomUtil.setPosition(c, b.min);
    c.width = size.x; c.height = size.y;
    c.style.width = size.x + "px"; c.style.height = size.y + "px";
    this._ctx.setTransform(1, 0, 0, 1, 0, 0);
    this._ctx.translate(-b.min.x, -b.min.y);
    const scale = 256 * Math.pow(2, this._map.getZoom()), o = this._map.getPixelOrigin();
    this._A = scale * MERC_K; this._Bx = 0.5*scale - o.x; this._By = 0.5*scale - o.y;
    this._blank(); this._paintAll();
  },
  _paintAll(){
    if (!this._ctx || !this._bounds || !this._gs.length) return;
    const ctx = this._ctx, b = this._bounds;
    ctx.globalAlpha = 0.8;
    for (const it of this._gs){
      const g = it.g, S = g.c*this._A, PX = g.ox*this._A + this._Bx, PY = -g.oy*this._A + this._By;
      const w = Math.max(1, S + 0.6);
      const x0 = b.min.x - w, x1 = b.max.x, y0 = b.min.y - w, y1 = b.max.y;
      let band = -1;
      for (let k = 0; k < g.n; k++){
        const doy = it.doy0 + g.t[k]/24;
        const bnd = Math.max(0, Math.min(DOY_BANDS.length-1,
          ((doy - DOY0) / Math.max(1, DOY1 - DOY0) * DOY_BANDS.length)|0));
        if (bnd !== band){ band = bnd; ctx.fillStyle = DOY_BANDS[bnd]; }
        const x = g.ix[k]*S + PX, y = PY - (g.iy[k]+1)*S;
        if (x < x0 || y < y0 || x > x1 || y > y1) continue;
        ctx.fillRect(x, y, w, w);
      }
    }
    ctx.globalAlpha = 1;
  },
});
const DOY_BANDS = Array.from({length:48}, (_,i) => ramp(i/47));
const fireRegion = {};                 // mergid -> region, for filtering the snap bands
let snapBand = null;
function rebuildSnapBands(){
  if (!snapBand) return;
  if (!(showKind.bnd && BAND_SRC === "snap")){ snapBand.setGrids([]); return; }
  const out = [];
  for (const y of YEARS){
    if (!activeYears.has(y) || !cache[y]) continue;
    for (const mid in byFireSim){
      const f = byFireSim[mid];
      if (!f || !f.g) continue;
      if (String(f.t0).slice(0,4) !== String(y)) continue;
      const r = fireRegion[mid];
      if (r != null && showRegion[r] === false) continue;
      const ms = Date.parse(String(f.t0).replace(" ","T") + "Z");
      const doy0 = (ms - Date.UTC(new Date(ms).getUTCFullYear(),0,1)) / 86400000 + 1;
      const gg = decodeGrid(f.g);
      if (gg) out.push({g:gg, doy0});
    }
  }
  snapBand.setGrids(out);
}

let burn = null;                       // the single BurnLayer instance (created with the map)
function clearHighlight(){ if (burn) burn.clear(); }
let popLayer = null;                   // "a pyroCb popped" markers during whole-fire animation
function clearPops(){
  if (popLayer){ map.removeLayer(popLayer); popLayer = null; }
}
const fmtDur = h => { const d = Math.floor(h/24);
  return d ? d+" d "+Math.round(h-24*d)+" h" : Math.round(h)+" h"; };

// Show a grid over an inclusive time range; returns {n, nb} (cells, cells before onset).
// `o` carries: span, rampT0 (hours), absT0Ms (grid t=0 in ms UTC), mergid, fit.
function showGrid(payload, lo, hi, o){
  o = o || {};
  const span = o.span||1, rampT0 = o.rampT0||0;
  const g = decodeGrid(payload);
  if (!burn) return {n:0, nb:0};
  if (!g){ burn.clear(); return {n:0, nb:0}; }
  const spr = o.mergid!=null ? bySpread[o.mergid] : null;
  if (burn._g !== g) burn.setGrid(g, span, rampT0, o.absT0Ms||0, spr);
  else burn.setContext(span, rampT0, o.absT0Ms||0, spr);
  const i0 = lbound(g.t, g.n, lo), i1 = ubound(g.t, g.n, hi);
  const n = burn.setRange(i0, i1);
  const z = Math.min(Math.max(lbound(g.t, g.n, 0), i0), i1);
  if (o.fit && FIT_ON_SELECT && n){
    const bb = burn.gridBounds();
    if (bb) try { map.fitBounds(bb.pad(0.25), {maxZoom:11}); } catch(e){}
  }
  return {n, nb: z - i0};
}

// t=0 of a +/-24 h window grid in ms UTC = the event's onset.
const winAbsT0 = p => Date.parse(String(p.datetime).replace(" ", "T") + "Z");

// Static wrapper used by popup mode and the panel's Static tab.
function renderActiveWindow(p, pre, post){
  const {n, nb} = showGrid(byPcbPix[p.pcb_id], -pre, post,
    {span:48, rampT0:-24, absT0Ms:winAbsT0(p), mergid:p.mergid, fit:true});
  if (p.active_is_blowup){
    showNotice("warn", `<b>&#9888; Blow-up fire</b> No VIIRS detections before onset &mdash; `
      +`first seen by the overpass <i>after</i> the pyroCb (FINDINGS F4). This is post-onset burning.`);
  } else if (!n){
    showNotice("info", `<b>No detections</b> in the &minus;${pre} h &hellip; +${post} h window.`);
  } else { hideNotice(); }
  return {n, nb};
}

// ---- bottom-left detail panel: static window slider + animation --------
let DETAIL_MODE = "panel";            // "panel" | "popup"
let FIT_ON_SELECT = false;            // pan/zoom the map to the event on select (opt-in)
let PANEL_MODE  = "static";           // "static" | "animate"
let winPre = 12, winPost = 12;        // sticky; grid {0,6,12,18,24}, sum >= 6
let loopSec = 5, animStyle = "cumulative", slideW = 6;   // slideW = sliding-window width (h)
const infoEl = document.getElementById("info");
L.DomEvent.disableClickPropagation(infoEl); L.DomEvent.disableScrollPropagation(infoEl);
let animReq = null;
function stopAnim(){ if (animReq){ cancelAnimationFrame(animReq); animReq = null; } }
function closeInfo(){ stopAnim(); infoEl.style.display = "none"; clearHighlight(); clearPops();
  hideNotice(); if (typeof updateGBar === "function") updateGBar(); }

function openInfoPanel(p){
  stopAnim();
  const dt = [
    ["Date/time (UTC)",p.datetime],["Region",p.region],["Fire (mergid)",p.mergid],
    ["pyroCb # on fire",p.pcb_index+" of "+p.n_pcb_total],
    ["Match",p.match_method+(p.match_dist_km!=null?` (${p.match_dist_km} km)`:"")],
    ["Fraction burned before",p.frac_burned_before],["Obs lag (h)",p.obs_lag_h],
  ];
  const seg = (id,a,b) => `<div class="seg"><button id="${id}A" class="on">${a}</button>`
    +`<button id="${id}B">${b}</button></div>`;
  infoEl.innerHTML =
    `<div class="hd">${esc(p.pcb_id)}<span class="x" title="close">&times;</span></div>`
    +`<div class="bd"><div style="color:#666;margin-bottom:3px">${esc(p.fire_name||"unnamed fire")}`
    +(p.active_is_blowup?` &nbsp;<span class="warn">&#9888;&#65039; blow-up</span>`:"")
    +`</div><table>${rows(dt)}</table>`
    +(lcTop3Html(p) ? `<div class="sl"><b>Land cover &mdash; `
        +(p.active_is_blowup?`first 6 h of burning <span class="warn">(blow-up)</span>`:`6 h pre-onset burn`)+`</b>`
        +`<div style="color:#666;font-size:11px">top 3 classes &middot; native 30 m`
        +(p.lc_n_px!=null?` &middot; ${p.lc_n_px} cells`:"")+`</div>${lcTop3Html(p)}</div>` : "")
    +`<div class="sl"><b>Active burn</b> ${seg("pm","Static","Animate")}`
    +`<div id="paneStatic">`
      +`<label>before <input type="range" id="slPre" min="0" max="24" step="6" value="${winPre}">`
      +`<span id="lbPre">${winPre}</span> h</label>`
      +`<label>after &nbsp;<input type="range" id="slPost" min="0" max="24" step="6" value="${winPost}">`
      +`<span id="lbPost">${winPost}</span> h</label></div>`
    +`<div id="paneAnim" hidden>`
      +`<label>scope <div class="seg"><button id="scA" class="on">&pm;24 h window</button>`
        +`<button id="scB"${byFireSim[p.mergid]?"":" disabled title='no fire track for this event'"}>Whole fire</button></div></label>`
      +`<div class="player"><button id="play" title="play / pause">&#9654;</button>`
        +`<span id="cur"></span><span id="curabs"></span></div>`
      +`<div id="scrubwrap"><div id="ticks"></div>`
        +`<input type="range" id="scrub" min="0" max="1000" step="1" value="0"`
        +` title="drag to scrub through the event"></div>`
      +`<label>loop <span class="seg" id="loopSeg">`
        + [3,5,10,15,20].map(s=>`<button data-s="${s}">${s}s</button>`).join("") + `</span>`
        +`<input type="number" id="loopCustom" min="1" max="120" step="1" placeholder="s"`
        +` title="custom loop length, 1-120 s"></label>`
      +`<label>${seg("as","Cumulative","Sliding")}</label>`
      +`<label id="slwRow" hidden>window <span class="seg" id="slwSeg">`
        + [6,24,48,168].map(s=>`<button data-h="${s}">${s>=168?"1 wk":s+" h"}</button>`).join("") + `</span></label>`
      +`<label>colour <span class="seg" id="cmSeg">`
        +`<button data-m="age" class="on">burn age</button>`
        +`<button data-m="av"${bySpread[p.mergid]?"":" disabled title='no spread series for this fire'"}>avg spread</button>`
        +`<button data-m="s95"${bySpread[p.mergid]?"":" disabled"}>spread p95</button></span></label>`
      +`<div id="tierkey" style="color:#666;font-size:11px"></div></div>`
    +`<div class="rd" id="rd"></div>`
    +`<div style="color:#666;font-size:11px">VIIRS active fire, 375 m grid, this fire only.`
    +` Steps snap to the 12 h VIIRS overpass cadence.</div>`
    +`</div></div>`;
  infoEl.style.display = "block";
  infoEl.querySelector(".x").onclick = closeInfo;

  const $ = id => infoEl.querySelector(id);
  const rd = $("#rd");
  const setReadout = (loH, hiH, n, nb) => {
    let km2 = null;
    if (loH === -winPre && hiH === winPost) km2 = (p.agrid||{})[winPre+"_"+winPost];
    rd.innerHTML = n
      ? `${km2!=null?km2.toFixed(1)+" km&sup2; &middot; ":""}${n} cells `
        +`(${nb} before / ${n-nb} after onset)`
      : `no detections in this window`;
  };

  // static
  const slPre = $("#slPre"), slPost = $("#slPost");
  const applyStatic = src => {
    let a = +slPre.value, b = +slPost.value;
    if (a + b < 6){ if (src==="pre") b = 6 - a; else a = 6 - b; slPre.value=a; slPost.value=b; }
    winPre = a; winPost = b;
    $("#lbPre").textContent = a; $("#lbPost").textContent = b;
    const {n, nb} = renderActiveWindow(p, a, b);
    setReadout(-a, b, n, nb);
  };
  slPre.oninput = () => applyStatic("pre");
  slPost.oninput = () => applyStatic("post");

  // ---- animate: a scrubbable cursor, not a captive loop ---------------------
  // `animU` in [0,1] is the single source of truth. Playing advances it by elapsed
  // time; the scrubber and the pyroCb ticks set it directly. So pausing, scrubbing
  // and jumping to an event all work without special cases.
  const fsim = byFireSim[p.mergid];
  let animScope = "win";                 // "win" | "fire"
  let animU = 0, lastTs = 0, prevCur = null, snapPts = [];
  const scrub = $("#scrub"), playBtn = $("#play");
  const spanOf = () => (animScope === "fire" && fsim) ? Math.max(1e-6, fsim.dur_h) : 48;
  const t0Of   = () => (animScope === "fire" && fsim) ? 0 : -24;
  const absT0Of = () => (animScope === "fire" && fsim)
    ? Date.parse(String(fsim.t0).replace(" ", "T") + "Z") : winAbsT0(p);
  const gridOf  = () => (animScope === "fire" && fsim) ? fsim.g : byPcbPix[p.pcb_id];
  const curOf  = u => t0Of() + u * spanOf();
  // Everything the map shows is a 12-hourly VIIRS overpass state, so snap the cursor
  // to the 12 h grid (plus the exact pyroCb onset times, which are worth landing on).
  const buildSnap = () => {
    const t0 = t0Of(), sp = spanOf(), pts = [];
    for (let h = t0; h <= t0 + sp + 1e-6; h += 12) pts.push(Math.min(t0 + sp, h));
    if (animScope === "fire" && fsim) fsim.pcb_h.forEach(h => { if (h >= t0 && h <= t0+sp) pts.push(h); });
    snapPts = [...new Set(pts.map(v => +v.toFixed(3)))].sort((a,b)=>a-b);
  };
  const snapCur = v => {
    if (!snapPts.length) return v;
    let best = snapPts[0], bd = Math.abs(v - best);
    for (const q of snapPts){ const d = Math.abs(v - q); if (d < bd){ bd = d; best = q; } }
    return best;
  };
  const gridCtx = () => ({span:spanOf(), rampT0:t0Of(), absT0Ms:absT0Of(), mergid:p.mergid});

  function pause(){ stopAnim(); lastTs = 0; playBtn.innerHTML = "&#9654;"; }
  function play(){
    if (animReq) return;
    if (animU >= 0.999){ animU = 0; prevCur = null; }     // restart from a finished run
    lastTs = 0; playBtn.innerHTML = "&#10073;&#10073;";
    animReq = requestAnimationFrame(frame);
  }

  // pyroCb markers over the scrubber: numbered, and vertically staggered when several
  // events fall close together in time so the labels never sit on top of each other.
  const drawTicks = () => {
    const tr = $("#ticks");
    tr.innerHTML = "";
    tr.style.height = "13px";
    if (animScope !== "fire" || !fsim) return;
    const sp = spanOf();
    const items = fsim.pcb_h.map((h, i) => ({h, i, x: 100 * h / sp})).sort((a,b) => a.x - b.x);
    // Greedy row packing: drop each marker into the FIRST row whose last marker is far
    // enough left. Round-robin used to collide again as soon as a 12-event fire bunched
    // several onsets into one afternoon.
    const GAP = 5.0;                     // % of track width one marker needs
    const rowLastX = [];
    let maxRow = 0;
    items.forEach(it => {
      let r = rowLastX.findIndex(x => it.x - x >= GAP);
      if (r < 0){ r = rowLastX.length; }
      rowLastX[r] = it.x; it.row = r;
      if (r > maxRow) maxRow = r;
    });
    items.forEach(it => {
      const b = document.createElement("b");
      b.textContent = fsim.pcb_ix ? fsim.pcb_ix[it.i] : (it.i + 1);
      b.title = fsim.pcb_id[it.i] + "  (" + fmtDur(it.h) + ")";
      b.style.left = it.x.toFixed(3) + "%";
      b.style.top = (it.row * 13) + "px";
      if (it.row){ b.classList.add("staggered");
        b.style.setProperty("--stem", ((maxRow - it.row) * 13 + 13) + "px"); }
      b.onclick = ev => { ev.preventDefault(); ev.stopPropagation();
        pause(); animU = Math.min(1, Math.max(0, it.h / sp)); prevCur = null; render(); };
      tr.appendChild(b);
    });
    tr.style.height = (14 + maxRow * 13) + "px";
    const wrap = infoEl.querySelector("#scrubwrap");
    if (wrap) wrap.style.paddingTop = (16 + maxRow * 13) + "px";
  };

  // dark-red dot at every pyroCb whose onset has been reached; a pulsing ring on the
  // one that just fired this frame (flashIdx).
  const renderPops = (cur, flashIdx) => {
    if (!popLayer) popLayer = L.layerGroup().addTo(map);
    popLayer.clearLayers();
    if (animScope !== "fire" || !fsim) return;
    for (let i = 0; i < fsim.pcb_h.length; i++){
      if (fsim.pcb_h[i] > cur) continue;
      const ll = fsim.pcb_ll[i];
      L.circleMarker(ll, {pane:"pop", renderer:R.pop, radius:4, color:"#fff",
        weight:1.5, fillColor:"#7d1128", fillOpacity:1}).addTo(popLayer);
      if (i === flashIdx)
        L.marker(ll, {pane:"pop", icon:L.divIcon({className:"", iconSize:[22,22],
          iconAnchor:[11,11], html:'<div class="pcbpulse"></div>'})}).addTo(popLayer);
    }
  };

  function render(){
    const cur = snapCur(curOf(animU));
    const lo = animStyle === "cumulative" ? t0Of() : cur - slideW;
    const absMs = absT0Of() + cur * 3600e3;
    const abs = isNaN(absMs) ? "" : new Date(absMs).toISOString().slice(0,16).replace("T"," ") + " UTC";
    if (animScope === "fire" && fsim){
      const n = showGrid(fsim.g, lo, cur, gridCtx()).n;
      $("#cur").textContent = fmtDur(cur) + " / " + fmtDur(fsim.dur_h);
      $("#curabs").textContent = abs;
      rd.innerHTML = `${n} of ${fsim.n_cell} cells &middot; coloured by ${BURN_MODE_LABEL[BURN_COLOR_MODE]}`;
      if (prevCur === null || cur < prevCur) prevCur = cur;   // init / loop wrap: no flash
      let flash = -1;
      for (let i = 0; i < fsim.pcb_h.length; i++)
        if (fsim.pcb_h[i] > prevCur && fsim.pcb_h[i] <= cur) flash = i;
      prevCur = cur;
      renderPops(cur, flash);
      if (flash >= 0)
        rd.innerHTML = `<b style="color:#7d1128">&#9733; pyroCb ${fsim.pcb_ix?("#"+fsim.pcb_ix[flash]+" "):""}`
          + `&mdash; ${esc(fsim.pcb_id[flash])}</b> &middot; ${fmtDur(cur)}`;
    } else {
      const {n, nb} = showGrid(byPcbPix[p.pcb_id], lo, cur, gridCtx());
      $("#cur").textContent = (cur >= 0 ? "+" : "") + cur.toFixed(0) + " h";
      $("#curabs").textContent = abs;
      setReadout(null, null, n, nb);
    }
    scrub.value = Math.round(animU * 1000);
  }

  const frame = ts => {
    if (lastTs) animU = (animU + (ts - lastTs) / (loopSec * 1000)) % 1;
    lastTs = ts;
    render();
    animReq = requestAnimationFrame(frame);
  };
  playBtn.onclick = () => animReq ? pause() : play();
  scrub.oninput = () => { pause(); animU = (+scrub.value)/1000; prevCur = null; render(); };

  const pair = (aId,bId,set,va,vb) => {
    $(aId).onclick = () => { set(va); $(aId).classList.add("on"); $(bId).classList.remove("on"); };
    $(bId).onclick = () => { set(vb); $(bId).classList.add("on"); $(aId).classList.remove("on"); };
  };
  // single-select button group: <span class="seg"> of <button data-{attr}="v">
  const radio = (seg, attr, cur, onPick) => {
    const bs = [...seg.querySelectorAll("button")];
    const mark = v => bs.forEach(b => b.classList.toggle("on", +b.dataset[attr] === v));
    bs.forEach(b => b.onclick = () => { const v = +b.dataset[attr]; onPick(v); mark(v); });
    mark(cur);
    return mark;
  };
  const markLoop = radio($("#loopSeg"), "s", loopSec, v => { loopSec = v; cust.value = ""; });
  const cust = $("#loopCustom");
  cust.oninput = () => {
    const v = Math.round(+cust.value);
    if (v >= 1 && v <= 120){ loopSec = v; markLoop(-1); }   // -1 clears every preset
  };
  radio($("#slwSeg"), "h", slideW, v => { slideW = v; render(); });
  pair("#asA","#asB", v => { animStyle = v; $("#slwRow").hidden = v !== "sliding"; render(); },
       "cumulative", "sliding");

  const AGE_KEY = `Colour = burn age: <span style="color:${AGE_STOPS[0]}">&#9632;</span> first`
    + ` &rarr; <span style="color:${AGE_STOPS[6]};text-shadow:0 0 1px #999">&#9632;</span> last.`;
  const SPD_KEY = m => `Colour = ${BURN_MODE_LABEL[m]}: `
    + `<span style="color:${SPREAD_STOPS[0]};text-shadow:0 0 1px #999">&#9632;</span> slow`
    + ` &rarr; <span style="color:${SPREAD_STOPS[6]}">&#9632;</span> fast`
    + ` (${SPREAD_CLAMP[m][0].toFixed(0)}&ndash;${SPREAD_CLAMP[m][1].toFixed(0)} m / 12 h).`;
  const setColourKey = () => {
    $("#tierkey").innerHTML = (BURN_COLOR_MODE === "age" ? AGE_KEY : SPD_KEY(BURN_COLOR_MODE))
      + ` Click a numbered <span style="color:#7d1128">&#9670;</span> to jump to that pyroCb.`;
  };
  if (BURN_COLOR_MODE !== "age" && !bySpread[p.mergid]) BURN_COLOR_MODE = "age";
  const cmBtns = [...$("#cmSeg").querySelectorAll("button")];
  cmBtns.forEach(b => {
    b.classList.toggle("on", b.dataset.m === BURN_COLOR_MODE);
    b.onclick = () => {
      if (b.disabled) return;
      BURN_COLOR_MODE = b.dataset.m;
      cmBtns.forEach(x => x.classList.toggle("on", x === b));
      setColourKey();
      burn.setContext(spanOf(), t0Of(), absT0Of(), bySpread[p.mergid]); burn.repaint();
    };
  });

  const setScope = s => {
    if (s === "fire" && !fsim) s = "win";
    animScope = s; prevCur = null; pause();
    $("#scA").classList.toggle("on", s === "win"); $("#scB").classList.toggle("on", s === "fire");
    if (s !== "fire") clearPops();
    buildSnap(); drawTicks(); setColourKey();
    // Show the whole extent first, so you can see what you are about to animate;
    // pressing play rewinds to the start.
    showGrid(gridOf(), t0Of(), t0Of() + spanOf(), {...gridCtx(), fit:true});
    animU = 1; render();
  };
  $("#scA").onclick = () => setScope("win");
  $("#scB").onclick = () => setScope("fire");

  $("#pmA").onclick = () => { PANEL_MODE="static"; pause(); clearPops();
    $("#pmA").classList.add("on"); $("#pmB").classList.remove("on");
    $("#paneStatic").hidden=false; $("#paneAnim").hidden=true; applyStatic();
    if (typeof updateGBar === "function") updateGBar(); };
  $("#pmB").onclick = () => { PANEL_MODE="animate";
    $("#pmB").classList.add("on"); $("#pmA").classList.remove("on");
    $("#paneStatic").hidden=true; $("#paneAnim").hidden=false;
    $("#slwRow").hidden = animStyle!=="sliding";
    setScope(animScope);
    if (typeof updateGBar === "function") updateGBar(); };

  (PANEL_MODE==="animate" ? $("#pmB") : $("#pmA")).onclick();
}

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

const PANES = ["lc","per","bnd","ign","pcb","orp"];
const R = {};
PANES.forEach((k,i)=>{ map.createPane(k); map.getPane(k).style.zIndex = 405+i*10;
                       R[k] = L.canvas({pane:k}); });
map.getPane("lc").style.pointerEvents = "none";                      // land cover: backdrop only, never eats a click
map.createPane("ignlbl"); map.getPane("ignlbl").style.zIndex = 452;   // stars: above bands, below pyroCb
map.createPane("lbl"); map.getPane("lbl").style.zIndex = 480;
map.createPane("hl");  map.getPane("hl").style.zIndex = 490;          // fire-activity canvas: on top of all
map.getPane("hl").style.pointerEvents = "none";                      // ...but transparent to clicks, so a new pyroCb stays selectable
map.createPane("pop"); map.getPane("pop").style.zIndex = 495;        // "a pyroCb popped" markers, above the burn canvas
map.getPane("pop").style.pointerEvents = "none";
R.pop = L.canvas({pane:"pop"});
burn = new BurnLayer({pane:"hl", padding:0.15}).addTo(map);
snapBand = new SnapBandLayer({pane:"bnd", padding:0.15}).addTo(map);

const ALL = ["lc","per","bnd","ign","pcb","pcbC","orp"];
const cache = {}, byFire = {};
const activeYears = new Set([YEARS[YEARS.length-1]]);
const showKind = {lc:false, per:true, bnd:true, ign:true, pcb:true, orp:true};
const showRegion = {};                 // country toggles, all on to start
REGIONS.forEach(r => showRegion[r] = true);
const OUTLINE = {per:false, bnd:false}; // draw as border only, no fill
let BAND_SRC = "pcb";                  // "pcb" = split at pyroCb | "snap" = 12 h snapshot
let PCB_POS  = "report";               // "report" = atlas lat/lon | "fire" = matched-burn centroid
let SHOW_ERR = false;                  // draw the +/-0.1 deg report uncertainty
const regionOK = ft => { const r = ft.properties.region;
  return r == null ? true : showRegion[r] !== false; };

// Land Cover of Canada palette (classes present in the boreal zone).
const LC_COL = {1:"#003d00",2:"#949c70",5:"#148c3d",6:"#5c752b",8:"#b38a2d",10:"#e1cf8a",
  11:"#9c7554",12:"#bad48f",13:"#408a70",14:"#6ba38a",15:"#e6ae66",16:"#a8abae",
  17:"#dc2126",18:"#4c70a3",19:"#fffaff"};
const LC_LAB = {1:"Needleleaf forest",2:"Taiga needleleaf",5:"Broadleaf deciduous",
  6:"Mixed forest",8:"Shrubland",10:"Grassland",11:"Polar shrub-lichen-moss",
  12:"Polar grass-lichen-moss",13:"Polar barren-lichen-moss",14:"Wetland",15:"Cropland",
  16:"Barren land",17:"Urban",18:"Water",19:"Snow & ice"};

// Land-cover composition of a pyroCb's 6 h-before-onset active-burn polygon (30 m).
function lcSwatch(c){return `<span style="display:inline-block;width:12px;height:12px;`
  +`border-radius:2px;margin-right:5px;vertical-align:-2px;background:${LC_COL[c]||'#999'}"></span>`;}
function lcTop3Html(p){
  const a=[];
  for(let i=1;i<=3;i++){const c=p["lc_top"+i+"_cls"], f=p["lc_top"+i+"_frac"];
    if(c!=null&&f!=null) a.push(`<div style="margin:2px 0">${lcSwatch(c)}`
      +`${esc(LC_LAB[c]||("class "+c))} <b>${Math.round(f*100)}%</b></div>`);}
  return a.join("");
}
function lcTop3Text(p){
  const a=[];
  for(let i=1;i<=3;i++){const c=p["lc_top"+i+"_cls"], f=p["lc_top"+i+"_frac"];
    if(c!=null&&f!=null) a.push((LC_LAB[c]||("class "+c))+" "+Math.round(f*100)+"%");}
  return a.join(", ");
}
let hideOrphans=false;

// One definition of the perimeter / band look, shared by first draw and recolour().
// "outline only" keeps the boundary and drops the fill, which is what you want when
// bands or perimeters are stacked over land cover.
function perStyle(y, ft){
  const hue = hueOf(y, ft);
  if (OUTLINE.per) return {color:hue, weight:1.6, opacity:1, fill:false};
  return {color:hue, weight:1, opacity:1, fill:true, fillColor:hue,
          fillOpacity:(ft.properties.farea_km2||0)>100?0.55:0.3};
}
function bndStyle(y, ft){
  const p = ft.properties, n = Math.max(1, p.n_bands-1);
  const t = n>0 ? p.band_i/n : 0.5;                  // 0 = core, 1 = outermost
  const hue = MODE==="year" ? YCOL[y] : doyColor(p.doy);
  if (OUTLINE.bnd)
    return {color:hue, weight:1.1+1.1*t, opacity:1, fill:false};
  // band order in LIGHTNESS so hue stays free for the colour mode
  const fill = MODE==="year" ? lighten(hue, 0.72-1.05*t) : hue;
  return {color:"#fff", weight: MODE==="year" ? 0.8 : 0.8+1.4*t, opacity:1,
          fill:true, fillColor:fill,
          fillOpacity: p.is_post_final_band ? 0.32 : (MODE==="year"?0.8:0.55+0.3*t)};
}

function buildLayers(y, gj){
  const perOpts = {
    pane:"per", renderer:R.per,
    style: ft => perStyle(y, ft),
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>Fire #${esc(p.mergid)}</b> ${y}<table>${rows([
        ["Area (km²)",p.farea_km2],["Duration (d)",p.duration_days],
        ["VIIRS pixels",p.n_pixels],["Land cover",p.lcc_final],
        ["Day of year",p.doy]])}</table>`);}
  };
  const ignOpts = {
    pane:"ignlbl",
    pointToLayer:(ft,ll)=>L.marker(ll,{pane:"ignlbl", icon:starIcon(hueOf(y,ft),
      ft.properties.ign_index)}),
    onEachFeature:(ft,l)=>{const p=ft.properties;
      l.bindPopup(`<b>Ignition ${esc(p.ign_index)} of ${esc(p.n_ign_total)}</b>`
        +` &middot; fire #${esc(p.mergid)} &middot; ${y}<table>${rows([
        ["Order within fire",p.ign_index+" of "+p.n_ign_total+" (by ignition time)"],
        ["Day of year",p.doy],["Fuel",p.ftype],["Area (km²)",p.farea_km2],
        ["VIIRS pixels",p.n_pixels]])}</table>`);}
  };
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
      if (DETAIL_MODE==="popup"){
        l.bindPopup(`<b>${esc(p.pcb_id)}</b><br>${esc(p.fire_name||"unnamed")}`
          +(p.active_is_blowup?` <span title="Blew up before it was ever observed">&#9888;&#65039;</span>`:"")
          +`<table>${rows([
          ["Date/time (UTC)",p.datetime],["Region",p.region],["Day of year",p.doy],
          ["Max injection alt (m)",p.max_inject_alt],["Rel. tropopause (m)",p.trop_relative_alt],
          ["Min BT (°C)",p.pcb_min_bt],["Match",p.match_method],["Dist. to fire (km)",p.match_dist_km],
          ["Fire (mergid)",p.mergid],["pyroCb # on this fire",p.pcb_index+" of "+p.n_pcb_total],
          ["pyroCb before / after",p.n_pcb_before+" / "+p.n_pcb_after],
          ["Ignitions in complex",p.n_ign_total],
          ["Area before (km²)",p.area_pre_km2],["Area after (km²)",p.area_post_km2],
          ["Fraction burned before",p.frac_burned_before],["Obs lag (h)",p.obs_lag_h],
          ["Active burn ±12h (km²)",p.active_area_km2_pre12_post12],
          ["Active-burn detections",p.active_n_pixels_pre12_post12],
          ["  before / after onset",p.active_n_before_onset_pre12_post12+" / "+p.active_n_after_onset_pre12_post12],
          ["Land cover (6 h pre-onset)",lcTop3Text(p)],
          ["Blew up before observed",p.active_is_blowup?"yes":""]])}</table>`);
      }
      l.on("click", () => DETAIL_MODE==="panel" ? openInfoPanel(p)
                                               : renderActiveWindow(p, 12, 12));}
  });
  const bndOpts = {
    pane:"bnd", renderer:R.bnd,
    style: ft => bndStyle(y, ft),
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

  return {
    // Land cover: a paletted PNG per fire (native-ish, GPU-composited) over a coarse
    // 10 km polygon backdrop for everything beyond them.
    lc:   (() => {
            const grp = L.layerGroup();
            L.geoJSON(gj.lc, {pane:"lc", renderer:R.lc, interactive:false,
              style: ft => ({stroke:false, fill:true,
                             fillColor:LC_COL[ft.properties.cls]||"#999", fillOpacity:1})}).addTo(grp);
            for (const k in (gj.lcpng||{})){
              const o = gj.lcpng[k];
              L.imageOverlay(o.d, o.b, {pane:"lc", opacity:1, interactive:false,
                className:"lcimg"}).addTo(grp);
            }
            return grp;
          })(),
    per:  L.geoJSON(gj.per, perOpts),
    bnd:  L.geoJSON(gj.bnd, bndOpts),
    ign:  L.geoJSON(gj.ign, ignOpts),
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
  if (showKind.lc) w.add("lc");
  if (showKind.per) w.add("per");
  if (showKind.bnd && BAND_SRC === "pcb") w.add("bnd");
  if (showKind.ign) w.add("ign");
  if (showKind.pcb) w.add(hideOrphans?"pcbC":"pcb");
  if (showKind.orp && !hideOrphans) w.add("orp");
  return w;
}

async function ensureYear(y){
  if (cache[y]) return cache[y];
  const d = DATA[y];
  const [per,ign,pcb,bnd,pix,lc,sim,spr] = await Promise.all(
    [d.per,d.ign,d.pcb,d.bnd,d.pix,d.lc,d.sim,d.spr].map(unpack));
  Object.assign(byPcbPix, pix);          // pcb_id -> +/-24 h cell grid
  Object.assign(byFireSim, sim);         // mergid  -> whole-fire series
  Object.assign(bySpread, spr||{});      // mergid  -> 12-hourly spread series
  (per.features||[]).forEach(f => { if (f.properties.region != null)
    fireRegion[String(f.properties.mergid)] = f.properties.region; });
  // lcpng is already-base64 PNG, so it is embedded raw - gzipping it would only add
  // a second base64 pass over incompressible data.
  cache[y] = buildLayers(y,{per,ign,pcb,bnd,lc,lcpng:d.lcpng||{}});
  return cache[y];
}

let refreshing=false, pending=false;
async function refresh(){
  if (refreshing){ pending=true; return; }
  refreshing=true;
  if (infoEl.style.display !== "block"){ clearHighlight(); clearPops(); hideNotice(); }
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
  if (typeof updateGBar === "function") updateGBar();
  forceFilter();                       // re-mask newly added layers (region + date)
  rebuildSnapBands();
  repositionPcb();
}

function recolour(){
  for (const y of Object.keys(cache)){
    const c = cache[y];
    c.per.setStyle(ft => perStyle(y, ft));
    c.bnd.setStyle(ft => bndStyle(y, ft));
    c.ign.eachLayer(l=>{
      const el=l.getElement&&l.getElement(), pth=el&&el.querySelector("path");
      if (pth) pth.setAttribute("fill", hueOf(y,l.feature));
    });
    ["pcb","pcbC"].forEach(k=>c[k].eachLayer(l=>{
      if (l.setStyle) l.setStyle({fillColor:hueOf(y,l.feature)});
      else if (l.getElement && l.getElement()) l.getElement().firstChild.style.background = hueOf(y,l.feature);
    }));
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
    <details class="layers anim" id="animBox"><summary style="cursor:pointer;font-weight:600;font-size:11px">Animate</summary>
      <div class="seg" id="gmodeSeg" style="margin-top:4px">
        <button data-g="off" class="on">Off</button><button data-g="year">Years</button>
        <button data-g="month">Months</button></div>
      <div id="gbody" hidden>
        <div id="gticks"></div>
        <input type="range" id="gscrub" min="0" max="1000" step="1" value="1000">
        <div class="grow"><button id="gplay" title="play / pause">&#9654;</button>
          <span id="glabel"></span><span id="gsub"></span></div>
        <label style="font-size:10px;color:#666">step <span class="seg" id="gspdSeg">
          <button data-s="0.4">fast</button><button data-s="0.9" class="on">med</button>
          <button data-s="2">slow</button></span></label>
      </div></details>
    <div class="layers"><b style="font-size:11px">Region</b>
      ${REGIONS.map(r=>`<label><input type="checkbox" data-region="${r}" checked> ${r}</label>`).join("")}
    </div>
    <div class="layers"><b style="font-size:11px">Colour by</b>
      <div class="seg"><button id="mYear" class="on">Year</button><button id="mDoy">Day of year</button></div>
      <div class="sub">Pick one year to watch a season progress.</div></div>
    <div class="layers">
      <label><input type="checkbox" data-kind="lc"> Land cover (&asymp;120 m over fires)</label>
      <label><input type="checkbox" data-kind="per" checked> Fire perimeters</label>
      <label style="margin-left:18px"><input type="checkbox" id="perOutline"> outline only</label>
      <label><input type="checkbox" data-kind="bnd" checked> Growth bands</label>
      <label style="margin-left:18px"><input type="checkbox" id="bndOutline"> outline only</label>
      <div class="seg" id="bndSrcSeg" style="margin-left:18px">
        <button data-b="pcb" class="on">split at pyroCb</button><button data-b="snap">12 h snapshot</button></div>
      <label><input type="checkbox" data-kind="ign" checked> Fire ignitions</label>
      <label><input type="checkbox" data-kind="pcb" checked> pyroCb events</label>
      <label><input type="checkbox" data-kind="orp" checked> Flag unmatched pyroCb</label>
    </div>
    <div class="layers">
      <label><input type="checkbox" id="hideOrp"> Hide unmatched pyroCb</label>
    </div>
    <div class="layers"><b style="font-size:11px">pyroCb position</b>
      <div class="seg"><button id="posRep" class="on">As reported</button><button id="posFire">Matched fire</button></div>
      <label><input type="checkbox" id="showErr"> show &plusmn;0.1&deg; report uncertainty</label>
      <div class="sub">93 % of atlas coordinates are rounded to 0.1&deg; (~11 km N-S).
      "Matched fire" re-plots each event on the centroid of its own 6 h pre-onset burn,
      which is known to one 375 m VIIRS pixel.</div></div>
    <div class="layers"><b style="font-size:11px">pyroCb click</b>
      <div class="seg"><button id="dPanel" class="on">Info panel</button><button id="dPopup">Popup</button></div>
      <label><input type="checkbox" id="fitSel"> Zoom map to event on select</label></div>`;
  L.DomEvent.disableClickPropagation(d); L.DomEvent.disableScrollPropagation(d);
  d.querySelectorAll(".chip").forEach(c=>c.onclick=()=>{
    const y=+c.dataset.year; activeYears.has(y)?activeYears.delete(y):activeYears.add(y); refresh();});
  d.querySelector("#allY").onclick  = ()=>{YEARS.forEach(y=>activeYears.add(y)); refresh();};
  d.querySelector("#noneY").onclick = ()=>{activeYears.clear(); refresh();};
  d.querySelectorAll("input[data-kind]").forEach(i=>i.onchange=()=>{
    showKind[i.dataset.kind]=i.checked; refresh(); if (i.dataset.kind==="lc") drawLegend();});
  d.querySelectorAll("input[data-region]").forEach(i=>i.onchange=()=>{
    showRegion[i.dataset.region]=i.checked; forceFilter(); updateGBar();});
  d.querySelector("#hideOrp").onchange = e=>{hideOrphans=e.target.checked; refresh();};
  d.querySelector("#fitSel").onchange = e=>{FIT_ON_SELECT=e.target.checked;};
  d.querySelector("#perOutline").onchange = e=>{OUTLINE.per=e.target.checked; forceFilter();};
  d.querySelector("#bndOutline").onchange = e=>{OUTLINE.bnd=e.target.checked; forceFilter();};
  d.querySelectorAll("#bndSrcSeg button").forEach(b=>b.onclick=()=>{
    BAND_SRC=b.dataset.b;
    d.querySelectorAll("#bndSrcSeg button").forEach(x=>x.classList.toggle("on",x===b));
    refresh();});
  const my=d.querySelector("#mYear"), mdy=d.querySelector("#mDoy");
  my.onclick  = ()=>{MODE="year"; my.classList.add("on"); mdy.classList.remove("on"); recolour();};
  mdy.onclick = ()=>{MODE="doy";  mdy.classList.add("on"); my.classList.remove("on"); recolour();};
  const pr=d.querySelector("#posRep"), pf=d.querySelector("#posFire");
  pr.onclick = ()=>{PCB_POS="report"; pr.classList.add("on"); pf.classList.remove("on"); repositionPcb();};
  pf.onclick = ()=>{PCB_POS="fire";   pf.classList.add("on"); pr.classList.remove("on"); repositionPcb();};
  d.querySelector("#showErr").onchange = e=>{SHOW_ERR=e.target.checked; repositionPcb();};
  const dp=d.querySelector("#dPanel"), du=d.querySelector("#dPopup");
  dp.onclick = ()=>{DETAIL_MODE="panel"; dp.classList.add("on"); du.classList.remove("on");
    Object.keys(cache).forEach(y=>cache[y].pcb.eachLayer(l=>l.unbindPopup&&l.unbindPopup()));};
  du.onclick = ()=>{DETAIL_MODE="popup"; du.classList.add("on"); dp.classList.remove("on");
    closeInfo();
    Object.keys(cache).forEach(y=>rebindPopups(y));};
  // animate controls
  d.querySelectorAll("#gmodeSeg button").forEach(b=>b.onclick=()=>{
    d.querySelectorAll("#gmodeSeg button").forEach(x=>x.classList.toggle("on",x===b));
    setAnimMode(b.dataset.g);});
  d.querySelectorAll("#gspdSeg button").forEach(b=>b.onclick=()=>{
    gStepSec=+b.dataset.s;
    d.querySelectorAll("#gspdSeg button").forEach(x=>x.classList.toggle("on",x===b));});
  return d;
};

// Re-attach pyroCb popups after a switch to "popup" mode (they are only bound at
// build time when the mode is already "popup").
function rebindPopups(y){
  const c = cache[y]; if (!c) return;
  c.pcb.eachLayer(l=>{
    const p = l.feature.properties;
    l.bindPopup(`<b>${esc(p.pcb_id)}</b><br>${esc(p.fire_name||"unnamed")}`
      +(p.active_is_blowup?` <span title="blew up before observed">&#9888;&#65039;</span>`:"")
      +`<table>${rows([
      ["Date/time (UTC)",p.datetime],["Region",p.region],["Fire (mergid)",p.mergid],
      ["pyroCb # on fire",p.pcb_index+" of "+p.n_pcb_total],
      ["Match",p.match_method],["Dist. to fire (km)",p.match_dist_km],
      ["Fraction burned before",p.frac_burned_before],["Obs lag (h)",p.obs_lag_h],
      ["Active burn ±12h (km²)",p.active_area_km2_pre12_post12],
      ["Active-burn detections",p.active_n_pixels_pre12_post12],
      ["  before / after onset",p.active_n_before_onset_pre12_post12+" / "+p.active_n_after_onset_pre12_post12],
      ["Land cover (6 h pre-onset)",lcTop3Text(p)]
      ])}</table>`);
  });
}
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
    <div class="sub" style="margin-top:5px">${BAND_SRC==="snap"
      ? `Growth bands = <b>12 h snapshot</b>, 375 m cells coloured by the day each cell
         first burned (same ramp as day-of-year below).`
      : `Growth bands: pale core = burned before the first pyroCb, darker = later runs;
         the faintest outer band is after the last pyroCb.`}
    Click a pyroCb for its active-burn window &amp; animation.</div>`
    + (PCB_POS==="fire" ? `<div class="sub" style="color:#7d1128">pyroCb plotted on the
        matched fire's 6 h pre-onset burn centroid (375 m), not the atlas point.</div>` : "")
    + (showKind.lc ? `<details class="lcdd"><summary>Land cover</summary><div class="body">`
        + [1,2,5,6,8,10,11,12,13,14,15,16,18].map(c =>
            `<div class="row"><i class="sq" style="background:${LC_COL[c]}"></i>${LC_LAB[c]}</div>`
          ).join("")
        + `<div class="sub" style="margin-top:4px">&asymp;120 m raster over each fire + 5 km,`
        + ` 10 km polygons beyond. Exact 30 m composition is in the pyroCb info panel.</div>`
        + `</div></details>` : "");
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

// ---- find-a-fire search -------------------------------------------------
(function(){
  const box = document.getElementById("search");
  const inp = box.querySelector("input"), ul = box.querySelector("ul");
  L.DomEvent.disableClickPropagation(box); L.DomEvent.disableScrollPropagation(box);
  let hits = [], sel = -1;
  const hide = () => { ul.style.display = "none"; sel = -1; };
  const query = q => {
    q = q.trim().toLowerCase();
    if (q.length < 2) return [];
    const num = q.replace(/[^0-9]/g, "");
    return SEARCH.filter(r =>
      (r.id && r.id.toLowerCase().includes(q)) ||
      (r.nm && r.nm.toLowerCase().includes(q)) ||
      (num && r.mg != null && String(r.mg).includes(num))
    ).slice(0, 12);
  };
  const draw = () => {
    ul.innerHTML = hits.map((r,i) =>
      `<li data-i="${i}" class="${i===sel?"sel":""}"><b>${esc(r.id||"")}</b> `
      + `<span>&middot; ${esc(r.nm||"unnamed")} &middot; fire ${r.mg??"?"} &middot; ${r.y}</span></li>`
    ).join("");
    ul.style.display = hits.length ? "block" : "none";
    ul.querySelectorAll("li").forEach(li =>
      li.onclick = () => pick(hits[+li.dataset.i]));
  };
  async function pick(r){
    if (!r) return;
    inp.value = r.id || r.nm || String(r.mg); hide();
    if (!activeYears.has(r.y)){ activeYears.add(r.y); await refresh(); }
    else await refresh();
    if (r.la != null && r.lo != null) map.setView([r.la, r.lo], Math.max(map.getZoom(), 7));
    const c = cache[r.y]; if (!c) return;
    let feat = null;
    c.pcb.eachLayer(l => { if (l.feature && l.feature.properties.pcb_id === r.id) feat = l.feature; });
    if (!feat && r.mg != null)
      c.pcb.eachLayer(l => { if (l.feature && +l.feature.properties.mergid === +r.mg) feat = feat || l.feature; });
    if (feat) (DETAIL_MODE === "panel" ? openInfoPanel(feat.properties)
                                       : renderActiveWindow(feat.properties, 12, 12));
  }
  inp.oninput = () => { hits = query(inp.value); sel = -1; draw(); };
  inp.onfocus = () => { if (inp.value.trim().length >= 2){ hits = query(inp.value); draw(); } };
  inp.onkeydown = e => {
    if (e.key === "Escape"){ inp.value = ""; hide(); inp.blur(); }
    else if (e.key === "ArrowDown"){ sel = Math.min(hits.length-1, sel+1); draw(); e.preventDefault(); }
    else if (e.key === "ArrowUp"){ sel = Math.max(0, sel-1); draw(); e.preventDefault(); }
    else if (e.key === "Enter"){ pick(hits[sel >= 0 ? sel : 0]); }
  };
  document.addEventListener("click", e => { if (!box.contains(e.target)) hide(); });
})();

// ---- animation engine (lives in the control panel, not over the map) ----
// The mode is an EXPLICIT user choice, never inferred from how many years happen to be
// selected. Inferring it was the year-sweep bug: sweeping back to one year flipped the
// mode to "month", which reset the animation and cleared the year set mid-play.
const gplay = () => document.getElementById("gplay");
const gscrub = () => document.getElementById("gscrub");
const gticks = () => document.getElementById("gticks");
let gMode = "off", gU = 1, gReq = null, gLastTs = 0, gStepSec = 0.9;
let dateFilterDoy = null, _lastFilterDoy = -1;
let yearsBeforeAnim = null;            // restored when the year animation is switched off
const MONTH_LO = 3, MONTH_HI = 10;
const gSteps = () => gMode === "year" ? YEARS.length : (MONTH_HI + 1 - MONTH_LO) + 1;
// step index 0..n-1 from the 0..1 scrub position
const gIdx = () => Math.min(gSteps() - 1, Math.floor(gU * gSteps() * 0.999999));
// pyroCb counts for one year, summed over the countries currently toggled on
function monthCounts(y){
  const per = MONTH_COUNTS[String(y)] || {};
  const out = {};
  for (let m = MONTH_LO; m <= MONTH_HI; m++){
    let n = 0;
    for (const r of REGIONS) if (showRegion[r]) n += ((per[r] || {})[m] || 0);
    out[m] = n;
  }
  return out;
}
const animYear = () => [...activeYears].sort((a,b)=>a-b)[0] ?? YEARS[YEARS.length-1];

// pyroCb position. The atlas rounds 93 % of coordinates to 0.1 deg (~11 km N-S), so
// "as reported" is a neighbourhood, not a point. "Matched fire" moves each event onto
// the centroid of its own 6 h pre-onset active burn, which is good to one 375 m pixel.
let errLayer = null;
function repositionPcb(){
  if (errLayer){ map.removeLayer(errLayer); errLayer = null; }
  if (SHOW_ERR) errLayer = L.layerGroup().addTo(map);
  for (const y of Object.keys(cache)){
    const c = cache[y];
    ["pcb","pcbC","orp"].forEach(k => c[k] && c[k].eachLayer(l => {
      const p = l.feature.properties, g = l.feature.geometry.coordinates;
      const useFire = PCB_POS === "fire" && p.fire_lat != null && p.fire_lon != null;
      if (l.setLatLng) l.setLatLng(useFire ? [p.fire_lat, p.fire_lon] : [g[1], g[0]]);
    }));
    if (errLayer && PCB_POS === "report" && activeYears.has(+y) && c.pcb)
      c.pcb.eachLayer(l => {
        const p = l.feature.properties, g = l.feature.geometry.coordinates;
        if (!p.loc_rounded || !regionOK(l.feature)) return;
        L.rectangle([[g[1]-0.05, g[0]-0.05],[g[1]+0.05, g[0]+0.05]],
          {pane:"orp", renderer:R.orp, color:"#e11", weight:1, opacity:0.7,
           dashArray:"3 3", fill:false, interactive:false}).addTo(errLayer);
      });
  }
}

// Visibility filter: country toggles, plus (in month animation) a day-of-year cursor.
// recolour() rebuilds every style from scratch first, then hidden features are knocked
// to zero opacity. Keyed on (doy, region set) so it is a no-op when nothing changed.
function filterKey(){ return dateFilterDoy + "|" + REGIONS.filter(r=>showRegion[r]).join(","); }
function forceFilter(){ _lastFilterDoy = "\\u0000"; applyFilters(); }
function applyFilters(){
  const key = filterKey();
  if (key === _lastFilterDoy) return;
  _lastFilterDoy = key;
  const doy = dateFilterDoy;
  recolour();
  const hide = ft => !regionOK(ft)
    || (doy != null && !(ft.properties.doy != null && ft.properties.doy <= doy));
  for (const y of Object.keys(cache)){
    const c = cache[y];
    ["per","bnd","orp"].forEach(k => c[k] && c[k].eachLayer(l => {
      if (l.setStyle) l.setStyle(hide(l.feature) ? {opacity:0, fillOpacity:0} : {opacity:1});
    }));
    ["pcb","pcbC","ign"].forEach(k => c[k] && c[k].eachLayer(l => {
      const off = hide(l.feature);
      if (l.setOpacity) l.setOpacity(off ? 0 : 1);
      else if (l.setStyle) l.setStyle({opacity: off?0:1, fillOpacity: off?0:1});
    }));
  }
  if (errLayer) repositionPcb();
}

function gStop(){ if (gReq){ cancelAnimationFrame(gReq); gReq = null; } gLastTs = 0;
  const b = gplay(); if (b) b.innerHTML = "&#9654;"; }
function gFrame(ts){
  if (gLastTs){
    gU += (ts - gLastTs) / (gStepSec * gSteps() * 1000);
    if (gU >= 1){ gU = 1; gLastTs = 0; gRender(); gStop(); return; }
  }
  gLastTs = ts; gRender(); gReq = requestAnimationFrame(gFrame);
}
function gPlayPause(){
  if (gReq){ gStop(); return; }
  if (gU >= 0.999) gU = 0;
  gLastTs = 0; gplay().innerHTML = "&#10073;&#10073;"; gReq = requestAnimationFrame(gFrame);
}

// Switch animation mode. "off" restores whatever the year chips said beforehand.
function setAnimMode(mode){
  gStop();
  if (gMode === "year" && mode !== "year" && yearsBeforeAnim){
    activeYears.clear(); yearsBeforeAnim.forEach(y => activeYears.add(y));
    yearsBeforeAnim = null; refresh();
  }
  gMode = mode; gU = 1;
  dateFilterDoy = null;
  if (mode === "year" && !yearsBeforeAnim) yearsBeforeAnim = [...activeYears];
  const body = document.getElementById("gbody");
  if (body) body.hidden = (mode === "off");
  forceFilter();
  buildTicks();
  if (mode !== "off") gRender();
}

function buildTicks(){
  const t = gticks(); if (!t) return;
  if (gMode === "year"){
    t.innerHTML = YEARS.map(y => `<span data-i="${YEARS.indexOf(y)}">${String(y).slice(2)}</span>`).join("");
  } else if (gMode === "month"){
    const mc = monthCounts(animYear());
    let h = "";
    for (let m = MONTH_LO; m <= MONTH_HI; m++)
      h += `<span data-i="${m - MONTH_LO}">${MONTH_ABBR[m]}<em>${mc[m]}</em></span>`;
    h += `<span data-i="${MONTH_HI + 1 - MONTH_LO}">all<em>&nbsp;</em></span>`;
    t.innerHTML = h;
  } else t.innerHTML = "";
}

let _gRendering = false;
function gRender(){
  if (_gRendering) return;               // refresh() -> updateGBar() can re-enter mid-frame
  _gRendering = true;
  try { _gRender(); } finally { _gRendering = false; }
}
function _gRender(){
  const sc = gscrub(), glabel = document.getElementById("glabel"),
        gsub = document.getElementById("gsub");
  if (!sc || gMode === "off") return;
  sc.value = Math.round(gU * 1000);
  const i = gIdx();
  [...gticks().children].forEach(sp => sp.classList.toggle("hot", +sp.dataset.i === i));
  if (gMode === "year"){
    const cy = YEARS[i];
    let changed = false;
    YEARS.forEach(y => { const want = y <= cy;
      if (want !== activeYears.has(y)){ want ? activeYears.add(y) : activeYears.delete(y); changed = true; } });
    if (changed) refresh();
    glabel.textContent = YEARS[0] + " – " + cy;
    gsub.textContent = (i + 1) + " of " + YEARS.length + " years";
  } else {
    const y = animYear(), mc = monthCounts(y), m = MONTH_LO + i;
    const done = m > MONTH_HI;
    dateFilterDoy = done ? null : (MONTH_DOY[m + 1] - 1);
    applyFilters();
    let cum = 0; for (let k = MONTH_LO; k <= Math.min(m, MONTH_HI); k++) cum += mc[k];
    glabel.textContent = done ? "whole of " + y : MONTH_ABBR[m] + " " + y;
    gsub.textContent = cum + " pyroCb by then";
  }
}

// Kept as the single hook the rest of the code calls when state moves under the bar.
function updateGBar(){
  if (gMode === "month") buildTicks();
  if (gMode !== "off") gRender();
}

// Wire the controls once the panel exists.
function wireAnim(){
  const b = gplay(), sc = gscrub();
  if (!b || !sc) return;
  b.onclick = gPlayPause;
  sc.oninput = () => { gStop(); gU = (+sc.value)/1000; gRender(); };
}

refresh().then(() => { wireAnim(); updateGBar(); });
</script>
</body>
</html>
"""

HTML = (HTML.replace("__DATA__", DATA_JS).replace("__YEARS__", YEARS_JS)
            .replace("__COLORS__", COLORS_JS).replace("__RAMP__", RAMP_JS)
            .replace("__SPREADCLAMP__", SPREAD_CLAMP_JS).replace("__MONTHCOUNTS__", MONTH_COUNTS_JS)
            .replace("__SEARCH__", SEARCH_JS).replace("__REGIONS__", REGIONS_JS)
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
