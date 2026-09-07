# CHANGELOG

Patch notes for the pyroCb × FireAtlas project. Newest first.

Rationale for each choice lives in [DECISIONS.md](DECISIONS.md); data discoveries in
[FINDINGS.md](FINDINGS.md).

---

## v0.16 — 2026-09-07 — positional accuracy, country filter, band/perimeter styling

### Added
* **Country filter.** Toggle Canada / Northern Asia / Alaska. Region is carried on every
  layer — events, perimeters, ignitions, growth bands and the snapshot bands — by
  inheriting each fire's region from the pyroCb matched to it.
* **Outline-only** toggles for growth bands *and* fire perimeters (border, no fill).
* **Growth bands from the 12 h snapshot** as an alternative to split-at-pyroCb. Rendered
  straight off the existing 375 m burn grid, coloured by day-of-year — no new payload,
  and finer than the polygon bands.
* **pyroCb position: "as reported" vs "matched fire"**, plus an optional ±0.1° report
  uncertainty box. See Accuracy below.

### Accuracy (the headline)
* **Reported pyroCb coordinates are 0.1°-rounded for ~93 % of events** — ±5.6 km N–S,
  ±2.8 km E–W at 60°N (FINDINGS **F13**). That is unrecoverable from the source. So the
  map now also ships `fire_lat`/`fire_lon`: the centroid of each event's own 6 h-pre-onset
  active burn, which is located to **one 375 m VIIRS pixel**, and can plot events there
  instead. `loc_rounded` flags which events were snapped.
* **Perimeter simplification 0.003° → 0.0003°** (~330 m → ~33 m) and the coordinate snap
  1e-5° → 1e-6° (~1.1 m → ~0.11 m). FEDS perimeters are alpha shapes over 375 m pixels,
  so ~33 m is an order of magnitude below the source's own granularity — visually lossless,
  where 330 m was plainly visible zoomed in.
* **Growth bands: the ±60 m buffer-close removed** (simplify-only at 50 m). This turned
  out to be the *largest* positional error in the map, far bigger than the perimeter
  issue we set out to fix: bands are thin rings, and a `buffer(+r).buffer(-r)` close
  deletes anything narrower than `2r`. Measured worst-case displacement was **4.3 km**
  (v0.15), versus **92 m** for simplify-only at the same payload scale. Shrinking the
  close radius does not help — any close jumps the max error to kilometres. See F14.

### Fixed
* **Year sweep deselected / failed to select years.** The playback mode was *inferred*
  from how many years were active, so sweeping back to a single year flipped the mode to
  "months" mid-play, resetting the animation and clearing the year set. Mode is now an
  explicit choice (Off / Years / Months) and never changes itself.
* **The playback bar is out of the middle of the map** — it is an "Animate" section in the
  control panel under the year chips, as requested.
* **Month tick labels no longer read as dates.** The slider used to sit over them; ticks
  are now above the track, with the count on its own muted sub-line (`Jun` / `40`) so it
  can't be misread as "Jun 40".
* **Month "pyroCb by then" follows the country toggles.**
* **Whole-fire animation tick labels crowd less**: greedy row packing (first row with
  space) instead of round-robin, with the stem drawn to the right length.

---

## v0.15 — 2026-09-07 — animation colour modes, media scrubber, global playback bar, search

### Added
* **Animation colour modes** (`colour` toggle in the Animate panel):
  * **burn age** — now a **purple** ramp (old → new). Red/yellow was intuitively misleading
    because those colours read as fire *intensity*.
  * **avg spread** / **spread p95** — a red→yellow ramp meaning *speed* (deep red = fast,
    yellow = slow), from the 12-hourly `avgspread` / `spread95` Snapshot attributes
    (`pyrocb/spread.py`, cache stage `spread`). Global ramp clamp = p2–p98 of the non-zero
    values. Disabled for fires with no spread series.
* **Media-player scrubber.** The fire-animation slider is now a real transport: drag to any
  moment, pause, study detail. The pyroCb ticks over it are **numbered** (`pcb_index`) and
  **stagger vertically** when several events fall close together in time.
* **Custom loop length** already existed (1–120 s); the presets stay 3/5/10/15/20 s.
* **Global playback bar** (bottom-centre), contextual:
  * ≥ 2 years shown → steps by **year**, cumulative reveal, current year highlighted.
  * exactly 1 year → steps by **month** (Mar–Oct); each tick is labelled with that month's
    pyroCb count, and playing reveals perimeters / ignitions / events month by month
    (`MONTH_COUNTS`, `applyDateFilter`).
  * a fire's Animate panel open → the bar hides; that panel's scrubber is the control.
* **Find-a-fire search** (top-centre): match on our `pcb_id`, the fire name, or the merge id;
  selecting a result loads that year, flies to the fire and opens its panel. Ships a tiny
  global index so it works before a year is loaded.

### Changed
* **All animation / timeline scrubbing snaps to the 12 h VIIRS overpass cadence** (plus the
  exact pyroCb onset times). NFP is actually overpass-resolution, not 12-hourly, but the
  honest cadence to expose is the overpass — sub-overpass scrubbing showed motion that isn't
  in the data.

### Notes
`spread.py` reads ~1200 Snapshot files per year on a cold cache (~90 s/year, then cached).
Spread modes apply per cell by binning its first-burn time into the covering 12 h step.

A performance/resolution pass. Everything got **finer and smaller at the same time**, because
all three hot paths were shipping the wrong representation.

### Fixed
* **Giant squares covering whole perimeters.** The whole-fire animation used to decimate
  detections by proximity to a pyroCb onset and then *inflate each survivor to its dedup
  cell* — so the quiet spans between events were drawn as 1 km and 3 km blocks, which on a
  small fire swallowed the entire perimeter. Gone: fire activity is now a **375 m grid
  everywhere** (one VIIRS pixel).

### Added / changed
* **Gridded burn-time payload** (`pyrocb/burngrid.py`). Detections snap to a 375 m grid built
  in EPSG:3857, reduced to the time each cell first burned (whole-fire) or to unique
  (cell, observation) pairs (±24 h window), sorted by time and delta-encoded.
  2023: **393 k cells at 375 m in 0.64 MB** vs 276 k decimated points at 1–3 km in 1.7 MB —
  **2.7× smaller and ~8× finer**. Because entries are time-sorted, every view the map offers
  is a *contiguous index range*.
* **Direct-to-canvas renderer** (`BurnLayer`, an `L.Renderer` subclass). Each frame used to
  rebuild a GeoJSON FeatureCollection and hand it to Leaflet, which created a `Path` object
  per box. Now a cell index becomes a canvas pixel with one multiply-add and a `fillRect`,
  and because the visible slice is contiguous, **advancing the animation draws only the newly
  revealed cells** instead of redrawing everything.
* **Colour = burn age.** Cells are shaded on an ember→flame ramp keyed to absolute time, so
  even a paused cumulative view shows the order the ground burned in. Free: cells are already
  time-sorted, so the fill colour changes a handful of times per repaint.
* **Media-player scrubber.** The animation is now a scrubbable cursor, not a captive loop:
  drag the slider to any moment, pause, and step through detail. pyroCb onsets appear as
  clickable diamonds above the scrubber — click one to jump straight to that event. Absolute
  UTC timestamp shown alongside elapsed time.
* **Custom loop duration**, 1–120 s, next to the 3/5/10/15/20 s presets.
* **Land cover is now a paletted PNG overlay per fire** (`pyrocb/lcraster.py`), warped to
  EPSG:3857 so Leaflet's placement is exact, with a 3×3 majority filter that is invisible at
  map scale and takes ~35 % off the PNG. **~120 m over every fire + 5 km for ~1.9 MB** (2023),
  against 1 km polygons for 1.1 MB — ~8× finer, GPU-composited, and it deletes the slow
  `features.shapes` + dissolve + simplify path that dominated build time and threw topology
  exceptions. The 10 km polygon backdrop beyond the fires stays.

### Notes
* `FIRESIM_*` and `LC_FINE_*` config are gone, replaced by `BURN_CELL_M`, `LC_PNG_MIN_M`,
  `LC_PNG_CAP_PX`.
* The exact 30 m land-cover composition per event (`active_fractions`) is untouched — that is
  the science number and it was always read at native resolution.

---

## v0.13 — 2026-09-06 — animation window/loop controls, native LC in the perimeter, pop flash

### Added
* **Sliding-window width** for the whole-fire animation: 6 h / 24 h / 48 h / 1 week, shown
  only when the *Sliding* style is picked (`slideW`, replaces the fixed `SLIDE_W = 6`).
* **Loop lengths 15 s and 20 s** (now 3 / 5 / 10 / 15 / 20 s).
* **"A pyroCb just popped" flash.** During whole-fire animation, every pyroCb whose onset the
  cursor has passed gets a dark-red dot; the one firing this frame gets a pulsing ring and the
  readout names the event. Onsets are also ticked on the progress bar (already in v0.12).
* **Finer land cover inside the fire perimeter.** A third LC tier clipped to the perimeter
  itself (no collar): **≈90 m**, but **only for fires whose bbox is under `LC_FINE_MAX_KM2`
  (150 km²)** — large perimeters at any sub-km resolution blow up the payload. Larger fires
  keep just the 1 km tier. Every LC feature carries `res` (90 | 1000 | 10000).
  Measured on the way here: true 30 m polygons over every pyroCb perimeter = ~14 MB
  compressed for 2023 alone (57 MB whole file); 90/300 m over all of them still ~7 MB
  (31 MB file). The **exact 30 m composition still reaches the user** through
  `active_fractions` (info panel + `pyrocb_fire_attributes.tsv`); only the cosmetic map fill
  is limited. A real 30 m *visual* needs a raster/tile path — IDEAS I-3. `LC_FINE_*` config.

### Fixed
* **Loop buttons were glitchy** — the old 2-button `.seg` plus a bolted-on "10 s" button never
  kept a single consistent `.on` state. Replaced with one `radio()` helper over a proper
  button group; same helper drives the new window-width control.

### Note
30 m vector land cover over every perimeter is ~14 MB compressed for 2023 (57 MB whole file);
even 90/300 m over all fires is ~7 MB. So the inside-perimeter tier is 90 m for small fires
only and the whole file lands at 16.8 MB (≈ v0.12). The 30 m numbers are still exact in the
info panel and attribute table; a closer *visual* needs a raster/tile path — IDEAS I-3.

---

## v0.12 — 2026-09-06 — whole-fire animation, LC-fraction fixes, legend dropdown

### Added
* **Animate the entire fire event.** The info-panel *Animate* tab gains a scope toggle:
  *±24 h window* (as before) or **Whole fire** — first VIIRS detection to last, cumulative or
  a sliding 6 h front, over the chosen 3/5/10 s loop, with pyroCb onsets marked on the
  progress bar. Data is adaptively decimated (`pyrocb/firesim.py`, cache stage `firesim`):
  native resolution within ±24 h of any pyroCb onset, 1 km / 6 h between 24–72 h, 3 km / 24 h
  in the quiet spans, capped at `FIRESIM_MAX_PTS` per fire. Tier is colour-coded in the
  animation. `FIRESIM_*` constants in config.
* **Zoom-to-event is now opt-in.** New "Zoom map to event on select" checkbox (default off);
  selecting a pyroCb no longer moves the map unless it's ticked.

### Changed
* **Land-cover legend is a proper collapsible dropdown** (caret, bordered, scrolls) instead
  of the inline `<details>`.
* **Blow-up land-cover fractions now have a defined source.** A blow-up has no burning in the
  6 h before onset, so its `(6, 0)` polygon is empty; the composition now falls back to the
  first 6 h of detected burning **at or after onset** (`active.first_burn_polys`), same
  construction — localised VIIRS pixel footprints only. The panel labels which case applies.
* Confirmed and locked in: the per-event land-cover fractions are sampled over **only the
  actual detected-burning pixel footprints** in the 6 h window (localised to the pyroCb's own
  cluster), at native 30 m — not a hull, not the perimeter, not the 1 km map layer.

---

## v0.11 — 2026-09-06 — land cover as backdrop, per-event LC fractions, panel fixes

### Added
* **Land-cover top 3 per pyroCb.** The info panel and the marker popup now list the three
  commonest *Land Cover of Canada* classes and their fractions within the **6 h-before-onset**
  active-burn polygon (`ACTIVE_TABLE_WINDOWS` gains `(6, 0)`), sampled at the raster's native
  **30 m** — not the 1 km map layer. Columns `lc_top{1,2,3}_cls` / `_frac` / `lc_n_px` are
  also written to `pyrocb_fire_attributes.tsv`/`.xlsx`.
  (`landcover.active_fractions`, `pipeline.landcover_fractions_cached`, cache stage `lcfrac`.)
* **Coarse land-cover backdrop.** Away from the fires the map now shows 10 km land cover
  instead of OpenStreetMap: a `LC_BG_BUFFER_M` (120 km) collar around **each** fire (never a
  hull over all of them — fires cross the dateline and EPSG:3979 is Canada-only), clipped to
  the raster, with the 1 km footprints punched out. `LC_BG_CELL_M` / `LC_BG_BUFFER_M` /
  `LC_BG_SIMPLIFY_M` in config; `res` (1000 | 10000) on every LC feature.

### Changed
* **Land cover is now an opaque backdrop**, drawn under everything, `pointer-events:none`
  (never eats a click, no popup). OpenStreetMap only shows through where there is no land
  cover. All other layers sit above it unchanged.
* **Removed the "Active burn area at onset (±12 h)" layer** — redundant now that clicking a
  pyroCb gives the window slider + animation. The ±12 h polygons are still in
  `pyrocb_active_areas.gpkg` and the attribute table.
* Land-cover legend swatches moved into a collapsed `<details>`; control panel and legend
  cap at viewport height and scroll internally, so a long legend no longer hides the panel
  behind it.

### Fixed
* **Could not select a second pyroCb without reloading.** The active-burn highlight canvas
  (pane `hl`, z-index 490) sat above the markers and its filled polygons swallowed the click.
  The `hl` pane is now `pointer-events:none` and the highlight layer `interactive:false`; the
  info panel also stops click/scroll propagation to the map.

---

## v0.10 — 2026-09-06 — land-cover context layer

### Added
* **Land cover** (`pyrocb/landcover.py`) — the 30 m *Land Cover of Canada* rasters
  (EPSG:3979; 2015 collection for 2013–2017, 2020 for 2018–2023) read windowed to each
  pyroCb fire's perimeter **+ 5 km**, majority-resampled to **1 km**, masked, polygonised
  by class and dissolved. New toggleable map layer under everything else, coloured by the
  15-class boreal palette, with a legend. Payload cost after 11 m coordinate snapping:
  ~0.1 MB/year, 0.6 MB for 2023 — the whole map went 9.9 → 11.3 MB.
* `LC_TIF_2015` / `LC_TIF_2020` / `LC_CELL_M` / `LC_BUFFER_M` / `LC_SIMPLIFY_M` in config.

### Note
This is the *context* layer only. Per-active-burn-polygon land-cover fractions for the
attribute table are the next step.

---

## v0.9 — 2026-09-06 — active burn localised, window slider + animation

### Fixed
* **Active-burn area was several fires, not one.** `NFP.mergid` is the season-final id, so a
  big fire's detections already include every neighbour that merged into it later; a ±12 h
  *time* filter did nothing about that. For the six largest 2023 pyroCb fires the "active
  area" spanned 130–210 km with the median detection 22–121 km from the pyroCb. Now the
  in-window detections are clustered and only the one at the pyroCb is kept
  (`active._localise`, KD-tree connected components). Headline active area: median 92 → **31
  km²**, max 15,000+ → **1,280 km²**; ~48 % of in-window pixels were dropped as detached.
  New columns `active_detach_km_*`, `n_clusters`. See FINDINGS F12.

### Added
* **Active-burn window slider.** Click a pyroCb → the bottom-left panel carries independent
  *before* / *after* sliders over {0, 6, 12, 18, 24} h (min 6 h total). The map redraws the
  VIIRS detections in that window live; the readout shows area (24-window lookup, rasterised
  server-side) and the before/after detection split.
* **Animation mode.** A Static / Animate toggle in the same panel. Animate sweeps a time
  cursor −24 h → +24 h on a 3 / 5 / 10 s loop, either **cumulative** (the fire growing from
  −24 h) or a **sliding 6 h window** (the front moving). Pure client-side from the shipped
  per-detection pixels — no extra payload.
* **Bottom-left detail panel** replaces the marker popup by default (toggle back to popups
  in the control panel).
* `active_area_km2_*` now for asymmetric windows `pre{A}_post{B}`; the GeoPackage carries the
  five in `ACTIVE_TABLE_WINDOWS` (12/12, 6/6, 24/24, 0/24 for blow-ups, 24/0).

### Changed
* `active.py` reworked: one localisation pass per event, shared by the polygons, the
  24-window area grid and the pixel export (`active.slider_data`).
* `make_notebook.py` now aborts if a generated code cell fails to `ast.parse` — a `"""` in a
  `code(r"""…""")` block has truncated the notebook twice.

---

## v0.8 — 2026-09-04 — click-to-highlight active burn area, blow-up flag

### Added
* **Click a pyroCb marker to highlight its active-burn area.** Builds a one-off overlay from
  that event's own `pyrocb_active_areas` polygon (±12 h window) in a dedicated top pane, and
  fits the map to it, instead of relying on the always-on cyan outline alone.
* **Blow-up flag.** `active_is_blowup` — every VIIRS detection in the ±12 h window postdates
  onset, i.e. the fire was first observed *after* its own pyroCb (F4: the ~10 h nightly
  overpass gap). Computed once in Python so the map and attribute table agree; a click on a
  flagged event pops a banner explaining the highlighted area is post-onset, not pre-onset.
  2023: 3 events, including `PCB_2023_0625_0000` — the `2023-R90877` type case from F4.
  A second banner covers the plain "no detections in this window" case.

### Removed
* **"Only fires with a pyroCb" toggle.** Redundant since v0.7 wired `MAP_ONLY_PCB_FIRES` into
  the payload itself — the layers it switched between were already identical. Also dropped
  the now-always-true "Linked to a pyroCb" popup row on perimeters/ignitions.

---

## v0.7 — 2026-09-04 — all-year map lands, ignition stars, export & build fixes

The v0.6 all-years run computed every stage but **never produced a working map or the
GeoPackage** — three defects in the last mile, now fixed.

### Fixed
* **Map generation crashed before writing HTML.** The DOY-colour window is fitted over a
  `doys` list that the payload loop never appended to → `pd.concat([])` →
  `ValueError: No objects to concatenate`. Every prior "all years" map on disk was the
  **stale 2023-only build** left behind by this crash — which is why bands and numbered
  pyroCb only showed for 2023. Now populated in the loop; the map carries all of 2013–2023.
* **`pyrocb_active_areas.gpkg` was never written.** `pipeline.export()` took an `active`
  argument and ignored it. Now writes the `active_burn` layer (1,188 rows = 396 events × 3
  windows; empty-window rows kept with null geometry). This is the sampling geometry for the
  30 m land-cover / ERA5 work, so it not existing was the real blocker.
* **`make_notebook.py` shipped a truncated notebook twice.** A `"""` docstring written
  inside a `code(r"""…""")` block ends the block early; the generator then wrote a short
  notebook and `_nb_run.py` was extracted from it, so edits silently didn't take. Added an
  `ast.parse` guard over every code cell at generation time — it now aborts loudly instead.

### Added
* **Ignitions are numbered stars.** Black-outlined 5-point SVG stars (hue = active colour
  mode) in a new `ignlbl` DOM pane, replacing the canvas dots, each labelled with
  **`ign_index`** — the ignition's chronological rank within its merge complex (`mergid`),
  so a 40-ignition complex like `mergid 7783` reads 1…40 in assembly order. `ign_t` /
  `ign_index` are defined once in `io.read_ignitions` and shared with `variables`.

### Changed
* **`MAP_ONLY_PCB_FIRES` is now actually applied** to the payload (it only flagged rows
  before). Essential with DOM star markers — 180 k SVG icons would lock the browser; the
  analysis population is ~1.5 k. Non-pyroCb fires and their ignitions no longer shipped.
* **Map geometry smoothed for size.** `smooth_for_map` closes the VIIRS pixel lattice
  (`buffer ±`) before simplifying, so the active-burn polygons — a union of thousands of
  375 m boxes, ~125 kB each — collapse to smooth blobs. **26 MB → 10 MB.** Full-resolution
  geometry is untouched in `pyrocb_active_areas.gpkg`.

---

## v0.6 — 2026-09-04 — active burn area, all years, 2021 data repair

### Added
* **Active burn area at pyroCb onset** (`pyrocb/active.py`) — the polygon of ground actually
  detected burning around each event, built by unioning every VIIRS detection's DS×DT footprint
  over a **symmetric** ±6/12/24 h window. Symmetric because the ~10 h nightly blind window means
  the run that fed the plume is usually first observed *after* onset (F3). Exported as
  **`pyrocb_active_areas.gpkg`** (WGS84) so it drops straight into QGIS or `rasterstats` for the
  30 m land-cover fractions and ERA5 sampling. Also folded into the event table as
  `active_area_km2_{6,12,24}h` and friends.
  2023 sanity: ±12 h median 92.8 km², i.e. 7.1% of the fire's final area (p90 39%).
* `MAP_ONLY_PCB_FIRES` — the map now ships only fires credited with a pyroCb and their ignitions,
  which is the analysis population, instead of ~170 k background polygons.
* Active-burn layer on the map, drawn as a cyan outline over the growth bands.
* NFP reader now carries `DS`/`DT` (pixel footprint dimensions).
* Per-year NFP id validation (`NFP_ID_MIN_MATCH`) with geometric repair and the
  `qc_nfp_repaired` / `qc_nfp_pixel_exact` flags.
* Stable event schema across years (`_SCHEMA_DEFAULTS`), so per-year frames concatenate.

### Fixed
* **Ignition-day blow-ups were unmatched.** A fire that ignites and blows up inside the nightly
  observation gap is first seen by the *next* overpass, so it exists in no snapshot at or before
  onset. Matching now searches up to `MATCH_FORWARD_STEPS` (2) overpasses ahead. Type case
  **`2023-R90877`** (`PCB_2023_0625_0000`): was unmatched with the nearest active fire 86 km away;
  now matches `mergid 11918` at **5.24 km** — the same fire its 06-28 sibling matched. That fire's
  first detection came **10.3 h after** the pyroCb already carrying 117 pixels, so it was plainly
  burning at onset, just unseen. 2023 unmatched 40 → 38. See F4.
* **The 2021 NFP files are from a different processing run** — caught by the standing `n_pixels`
  assertion, which failed for 2021 and no other year. Two defects: coordinates ship as **lon/lat**
  instead of EPSG:3571 metres (every other year is metres), and `mergid` comes from a **stale id
  space** (18,931 ids vs the final product's 11,688). The pixel *total* agrees to 0.04%, so the
  detections are all there and correctly located — only the labels are wrong. Now auto-detected
  and reprojected, then re-attributed by point-in-polygon: 7.0% → **68.3%** of fires exact,
  83.1% within ±2 pixels. 2021 is **flagged, not faked** — treat its NFP-derived columns
  (`n_pixels_*`, `frp_*`, `active_*`) as approximate. See F11.

### Changed
* `SNAPSHOT_YEARS = None` — the heavy analysis now covers all of 2013–2023.

---

## v0.5 — 2026-09-04 — version control baseline

### Added
* `CHANGELOG.md`, `DECISIONS.md`, `FINDINGS.md`.
* Git repository with LFS for `*.html` / `*.xlsx` / `*.gpkg`, so every historical map stays
  retrievable. ~4.3 GB of source archives excluded via `.gitignore`.

---

## v0.4 — 2026-09-04 — the `pyrocb` package: matching, growth bands, attribute table

The big one. Analysis moved out of notebook string literals into an importable package
(`config`, `io`, `timebins`, `matching`, `bands`, `variables`, `pipeline` — ~1,300 lines).

### Added
* **pyroCb → fire matching** against the 12-hourly snapshot (containment first, then nearest
  active fire within 20 km), with a dormant-fire rescue tier and a 10/20/30 km cap sensitivity.
  2023: 48 inside, 71 nearest, 40 unmatched, 68 distinct fires.
* **Concentric growth bands** — each fire's footprint split at *every* pyroCb it produced, built
  as a running cumulative union so bands are provably disjoint and exactly tiling. 180 bands over
  68 fires in 2023. Verified: 0 overlapping fires, max relative error 2e-4, 0 negative-area
  bands, 1.68 km² of slivers dropped.
* **Per-event attribute table**, 159 events × 101 columns → `pyrocb_fire_attributes.tsv` / `.xlsx`
  plus `pyrocb_columns.tsv` (data dictionary). Includes the pyroCb-sequence variables
  (`n_pcb_before` / `n_pcb_after` / `pcb_index`), ignition counts by merge complex
  (`n_ign_before/after/total`), exact NFP pixel and FRP splits with 6/12/24 h windows, snapshot
  state at onset, and QC columns (`obs_lag_h`, `qc_stale_obs`, `qc_final_smaller`, …).
* **Day-of-year colour mode** on the map, with the ramp window fitted area-weighted over the whole
  dataset → (69, 283), covering 98.2% of burned area and extending automatically with new years.
* **Growth-band map layer**, band order encoded in lightness so hue stays free for the colour mode.
* Per-(stage, year) parquet caching keyed on config hash + input mtimes. Warm rerun: 0.3 s.
* `pyarrow` dependency (parquet + GeoParquet).

### Fixed
* **Ignitions were filtered on `fireid` instead of `mergid`** — since ignition `fireid` is unique
  per row, the "only fires with a pyroCb" toggle kept only each merged fire's *primary* ignition
  and dropped the other 1,444 (2023). Count went 68 → 518. See F1.
* Boolean columns became the strings `"True"`/`"False"` on parquet round-trip.
* `EPSG:3571` was mislabelled "EASE-Grid 2.0 North" in a comment; it is *North Pole LAEA Bering
  Sea* (EASE-Grid 2.0 North is 6931).

### Changed
* **2012 dropped from the map** — now falls out of `YEARS = None` resolving to "years with fire
  data *and* pyroCb records", rather than being hard-coded.
* Snapshot reads split into cheap validation at read time and precision-snapping only where a
  boolean needs it: 8.0 s → 0.7 s per snapshot.

---

## v0.3 — 2026-09-03 — all years, OSM only, filter toggles

### Added
* **All 12 years (2012–2023)** with per-year chips, `All`/`None`, and lazy per-year
  decompression — 169,439 fires, 184,154 ignitions, 484 pyroCb events in one 28 MB file.
* Per-year colour ramp (cool → warm, 2012 blue → 2023 crimson).
* **"Only fires with a pyroCb"** and **"Hide pyroCb >20 km from a fire"** toggles.
* Runtime smoke test (`smoke.js`) driving the map's JS against stubbed Leaflet/DOM.

### Changed
* **OpenStreetMap only** — CARTO and Esri basemaps and the basemap switcher removed.
* Perimeter simplification relaxed to ~330 m (one VIIRS pixel) to fit 12 years in one file.

### Fixed
* Legend colour swatches were being blanked by the year-chip restyler
  (`.chip` → `.chip[data-year]`).
* Year clicks made during a load were dropped; `refresh()` now coalesces them.

---

## v0.2 — 2026-09-03 — documented notebook, higher resolution

### Added
* `Data exploration.ipynb` — a documented notebook generated by `make_notebook.py`.
* pyroCb-to-perimeter proximity analysis with green/amber/red overlap classes.
* Dataset boundary outlines (convex + concave hulls). *Removed again in v0.3 as redundant.*

### Changed
* Perimeter simplification 1 km → ~3 m, coordinates snapped at ~1 cm; gzip+base64 embedding with
  in-browser `DecompressionStream` keeps the file at 15 MB.

### Fixed
* **pyroCb coordinates were described as the plume position.** `PCB readme.pdf` Table 1 defines
  `Lat`/`Lon` as the *approximate location of the **fire*** linked to the pyroCb. The "plume
  advection" explanation for non-overlapping points was wrong and was removed; the real causes are
  the 0.1° storage precision (~11 km) and genuine gaps in the fire atlas.
* Terminology: "pulse" → "event" throughout, matching the readme (`pulse_code` separately flags
  single- vs multi-pulse).

---

## v0.1 — 2026-09-03 — first map

* Standalone `fire_pcb_map_2023.html`: 2023 final perimeters, ignitions and pyroCb events with
  layer toggles, canvas rendering and click popups.
