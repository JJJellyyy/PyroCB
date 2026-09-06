# CHANGELOG

Patch notes for the pyroCb × FireAtlas project. Newest first.

Rationale for each choice lives in [DECISIONS.md](DECISIONS.md); data discoveries in
[FINDINGS.md](FINDINGS.md).

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
