# CHANGELOG

Patch notes for the pyroCb × FireAtlas project. Newest first.

Rationale for each choice lives in [DECISIONS.md](DECISIONS.md); data discoveries in
[FINDINGS.md](FINDINGS.md).

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
