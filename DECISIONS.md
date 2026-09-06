# DECISIONS

Every methodological choice, why it was made, and what it rules out. Append here whenever a
choice is made — the point is that a reviewer (or you in six months) can reconstruct *why* the
numbers are what they are.

Format: **D-n — decision.** Context → choice → rationale → consequence/alternative rejected.

Data-driven discoveries live in [FINDINGS.md](FINDINGS.md); release history in
[CHANGELOG.md](CHANGELOG.md).

---

## Scope and identity

### D-1 — Join everything on `mergid`
`fireid` is the pre-merge id and is unsafe (see F1: incomplete in snapshots, a meaningless row
counter in NFP). `mergid` is 1:1 across all four products.
**Rejected:** joining on `fireid`, which silently drops merged components and can attach a
polygon 92 km from the real fire.

### D-2 — `YEARS = None` means "years with fire data *and* pyroCb records"
Resolves to 2013–2023, so **2012 drops out on its own** for the stated reason (no pyroCb
coverage) rather than being hard-coded. Adding 2024+ data needs no code change.
**Rejected:** `YEARS = [2013, …, 2023]` literal, which silently rots.

### D-3 — `SNAPSHOT_YEARS` gates the heavy analysis separately from the map
The map is cheap per year; snapshot/NFP matching is not. Two knobs let the map cover all years
while the expensive pipeline runs on a subset. `None` → same as `YEARS`.

### D-4 — `pcb_id = PCB_YYYY_MMDD_HHMM` + deterministic `a/b/c` suffix
The requested format collides on outbreak days (F8). Suffix ordered by (lat, lon) so ids are
unique *and* stable across rebuilds.
**Rejected:** appending `mergid` — informative, but the id would change whenever matching
parameters change, which defeats the purpose of a stable key.

---

## Coordinate systems and geometry

### D-5 — All metric work in EPSG:3571
"WGS 84 / North Pole LAEA Bering Sea": equal-area (so `area/1e6` is a defensible km²) and centred
on lon 180, making the dateline — where Chukotka fires sit — the *least* distorted part of the
projection. Three of four inputs are already native 3571, so booleans operate on bit-identical
source coordinates. Convert to 4326 only for the map payload.
**Rejected:** working in degrees. A metre-scale sliver threshold is impossible at 50–72°N
(0.00001° is 1.1 m in latitude but 0.35 m in longitude at 70°N), and antimeridian differences
produce world-spanning garbage.

### D-6 — `make_valid`, never `buffer(0)`
`buffer(0)` silently deletes the smaller loop of a bowtie and can drop whole MultiPolygon lobes —
it changes area without saying so, which is fatal when the output *is* an area. `make_valid`
preserves area. Cost is negligible: 0 of 1413 snapshot geometries are invalid.

### D-7 — Snap both boolean operands to a common 1 m grid
`set_precision(1.0, mode="valid_output")` before every difference makes the operation
topologically exact and removes hairline slivers. 1 m is 1/375 of a VIIRS pixel, so it discards
nothing real. Applied only to the few fires that get banded — it costs ~2 s per 1400 geometries,
so it is deliberately *not* run over whole snapshots (that alone took read time from 8 s to 0.7 s).

### D-8 — Two sliver thresholds at different levels
`SLIVER_KM2 = 0.02` per exploded polygon part (observed slivers after snapping are < 1e-4 km²),
so a genuine 40 km² band is not discarded for carrying a 2 m² hair. `MIN_BAND_KM2 = 0.5`
(≈ 3.5 VIIRS pixels) on the band total, for whether to emit geometry.
**A sub-threshold band keeps its row and loses only its geometry**, so band index stays aligned
with event index.

---

## Matching pyroCb to fires

### D-9 — Match against the *snapshot*, not the season-final perimeter
The final footprint of a big fire can be 100 km across; a point inside it says nothing about
whether the fire had reached that spot on the day, and a fire that ignited weeks later would
still swallow the point. The snapshot at onset is the fire as it actually was.

### D-10 — Containment first, then nearest active fire within a cap
`MATCH_CAP_KM = 20`, ~2.5× the worst-case 0.1° rounding displacement (~7.8 km). A 10/20/30 km
sensitivity is reported every run; the nearest fire never changes with the cap, only acceptance.
**Consequence:** `match_dist_km > 5 km` should be treated as provisional.

### D-11 — Three rescue tiers, in order
1. **Forward search** (`_rescue_new_fire`, up to 2 overpasses) for fires that blew up inside the
   nightly blind window before being observed — see F4. Recovered 2 of 40 in 2023.
2. **Dormant search** against final perimeters whose start–end interval brackets onset, for fires
   hidden by cloud that overpass. Rescued **0** in 2023 — which is itself the evidence that the
   remaining unmatched events are genuinely absent from the atlas (F10).
3. Otherwise `unmatched`, with fire columns null.

### D-12 — Keep unmatched events in the export
They are non-random (F10), so silently dropping them biases the sample toward well-detected
boreal forest fires. Flagged via `match_method`, never deleted.

---

## Growth bands

### D-13 — Split at *every* pyroCb, not just the first
A fire with N distinct pyroCb half-day steps yields N+1 concentric bands. Measured for 2023:
51 of 70 fires have a single step (a simple 2-band split), only 3 exceed 3 steps (max 9 → 10
bands), so the rings stay legible.
**Rejected:** splitting only at the first or the largest pyroCb — loses the pulse-by-pulse growth
that motivates the whole exercise.

### D-14 — Build a running cumulative *union*, never difference raw snapshots
`S_i = S_{i-1} ∪ V(snap(t_i))`, then `B_i = S_i \ S_{i-1}`. Consecutive snapshots are not
perfectly nested (F5), so plain differencing yields negative-area bands and rings that punch
holes in each other.
**Guarantees:** bands disjoint, non-negative, exactly tiling `S_final`; `area_pre` monotone in
event index; `frac_burned_before ∈ [0,1]` by construction; high fill opacity is safe on the map.
**Verified for 2023:** 0 overlapping fires, max relative error 2e-4, 0 negative bands.
**But:** `qc_regress_km2` / `qc_final_smaller` record what the union papered over, so the
correction is auditable rather than invisible.

### D-15 — Pick the pre-band step from NFP observation times (Tier A), not a clock
`step_before = max{ s : t_med(mergid, s) < onset }`, using the **median** detection time so one
stray late pixel cannot disqualify an overpass whose bulk precedes onset. Needs no local-time
model at all. Tier B (local-solar clock, using the *fire's* centroid longitude, not the atlas's
rounded lon) only if the fire has no NFP; Tier C walks back over missing files.
**2023 result:** Tier A resolved 119/119 events, 0 back-walks.

### D-16 — Deduplicate bands, keep events
Two events on the same fire in the same half-day share a band, so every step-derived column is
**identical** for them — genuinely unresolvable at 12-hourly sampling. Both events keep their own
row; `dedup_collapsed` flags them. NFP-derived columns still separate them.
**Documented loudly** because otherwise a scatter of `area_pre_km2` shows vertical stacks of
identical values and looks like a bug.

---

## Presentation

### D-17 — Colour by year *or* day-of-year, with the DOY window fitted to the whole dataset
Fitted **area-weighted**: an unweighted percentile over fire counts gives (16, 348), dragged out
by tiny January crop fires, wasting the entire gradient on months with no boreal burning. The
area-weighted 1–99% span is (69, 281), widened to cover every pyroCb → **(69, 283)**, 98.2% of
all burned area. Refitting on new years extends it automatically.
**Rejected:** per-year fitting — maximises within-year contrast but makes colours incomparable
between years.

### D-18 — Band order in lightness, hue reserved for the colour mode
Year mode: hue = year colour, lightness stepped 0.78 (pale pre-pyroCb core) → 0.34 (dark late
runs). DOY mode: hue = day-of-year, so band order moves to stroke weight. The outermost band
("everything after the last pyroCb, at unknown times") is drawn fainter so it cannot be misread
as a growth episode of known duration.

### D-19 — Lazy per-year decompression in the map
Each year's payload is gzip+base64 embedded and unpacked only when the year is switched on. That
is what lets a 10 MB file open quickly.

### D-24 — The map ships only the pyroCb analysis population
`MAP_ONLY_PCB_FIRES` pre-filters the payload to fires credited with a pyroCb and their ignitions
(~1.5 k ignitions across 11 years), not every tracked fire (~170 k). Forced by the move to DOM
star markers for ignitions — 170 k `<svg>` icons lock the browser; the background fire population
was never the point of this map anyway.
**Rejected:** keeping the full background and only styling the pyroCb fires — fine for canvas
dots, impossible for labelled markers.

### D-25 — Ignitions are numbered stars; the number is order within the merge complex
`ign_index` is the ignition's 1-based chronological rank within its `mergid` (ties on `fireid`),
deliberately parallel to `pcb_index`. For a merge complex the stars then read 1…N in the order
the complex assembled, making visible the structure that `n_ign_total` / `ign_merge_span_days`
only summarise. Defined once in `io.read_ignitions` (with `ign_t`) so the map and the attribute
table cannot disagree.
**Rejected:** numbering by `fireid` (an opaque global id, not meaningful per fire); labelling
only multi-ignition complexes (inconsistent with how every pyroCb is labelled).

### D-26 — Smooth map geometry, keep full resolution in the GeoPackage
The active-burn polygon is a union of thousands of 375 m VIIRS pixel rectangles — ~125 kB of
coordinates each, 14 MB for one year. `smooth_for_map` closes the lattice (`buffer ±200 m`) into
solid blobs before simplifying to 150 m; bands get `buffer ±60 m` + 150 m. The map only needs to
show *roughly where*; the exact sampling geometry — what land-cover fractions and ERA5 are
computed against — is untouched in `pyrocb_active_areas.gpkg`. Cut the map 26 MB → 10 MB.

### D-27 — Drop the "only fires with a pyroCb" toggle
Made redundant by D-24: once the payload itself is pre-filtered to the analysis population,
the "linked" and "unlinked" layers it toggled between are identical, so the checkbox did
nothing. Removed the toggle, the duplicate `perL`/`ignL` layers, and the now-always-true
"Linked to a pyroCb" popup row.

### D-28 — Highlight the active-burn area by click, not as a permanent overlay
Clicking a pyroCb marker builds a one-off highlight layer from that event's own
`pyrocb_active_areas` feature (indexed by the globally-unique `pcb_id`, so it works
regardless of which years are toggled on) and fits the map to it, rather than relying on
the always-on cyan `act` layer alone. Two things fall out of doing this per-click instead of
per-layer-style: (1) a genuinely empty window (no detections) and (2) a **blow-up** event —
every detection in the ±12 h window postdates onset (F4) — can each get their own notice
banner, which a shared layer style cannot express. `active_is_blowup` is computed once in
Python (`active_n_pixels_12h > 0 & active_n_before_onset_12h == 0`) rather than re-derived in
JS, so the map and the attribute table can never disagree on which events are flagged.
**Rejected:** styling blow-up fires differently in the always-on `act` layer — cannot show
the "no data at all" case (empty geometry has no style to apply), and clutters the base map
with a distinction only relevant when you're already looking at one specific event.

### D-29 — Localise the active-burn area to the pyroCb's own fire cluster
The window filter is on `NFP.mergid`, which is the season-final id, so it silently includes
every fire that merged in later (F12). After the time filter the detections are clustered and
only the one holding the pyroCb (or the nearest within `ACTIVE_ATTACH_CAP_KM`, for the
0.1°-rounded location landing in a gap) is kept. Clustering is a KD-tree connected-component
join, not a shapely buffer/unbuffer — the latter over thousands of parts was the runtime.
**Rejected:** using the onset snapshot polygon as the spatial mask — it is *also* dissolved by
the final `mergid`, so for the merge cases it is the same over-large blob.

### D-30 — The map slider ships pixels, not 24 polygons
Independent before/after sliders imply 24 (pre,post) windows. Precomputing a smoothed polygon
for each, per event, is ~90 MB and half an hour of `union_all`. Instead the map ships the
localised per-detection footprints once (~1 MB gz for all years) and the client draws the
boxes for whatever window — or animation frame — is showing. Exact union polygons are kept
only for the five `ACTIVE_TABLE_WINDOWS` (GeoPackage + `act` layer); the slider's area readout
uses a 100 m rasterised area (sub-percent error).
**Rejected:** turf.js `union` in the browser — too slow to run on every slider tick.

### D-31 — Animation is a client-side sweep, not pre-rendered frames
"Active-area GIF" = a `requestAnimationFrame` loop moving a time cursor over the same shipped
pixels, cumulative or sliding-window, 3/5/10 s loop. No frames, tiles or video in the payload;
it reuses one Leaflet layer via `clearLayers()`/`addData()` so it stays smooth.

### D-32 — Land cover: read in the raster's own CRS, decimate on read, polygonise
The *Land Cover of Canada* rasters are EPSG:3979, 190001 × 160001, uint8. We never
reproject the raster - the fire perimeters go to 3979 instead - and we never read it at
full 30 m: `rasterio` reads windowed to (perimeter + 5 km) with `out_shape` set for 1 km
and `Resampling.mode`, so a 1 km cell is the commonest 30 m class under it. Then mask,
`rasterio.features.shapes`, dissolve by class, simplify 2 km. 2015 collection for fire
years < `LC_YEAR_SPLIT` (2018), 2020 collection after.
**Rejected:** shipping a PNG image overlay per fire - crisper, but 255 `L.imageOverlay`
elements, and the polygon layer is only ~1 MB for all years after coordinate snapping.
**Rejected:** one union bbox per year - the 2023 pyroCb fires span BC to Quebec to Yukon,
a 5300 × 7900 km box that is 99 % nodata; per-fire windows keep the reads tiny.

### D-20 — Analysis lives in a package, not in notebook string literals
`make_notebook.py` emits the notebook; the logic sits in `pyrocb/` so it can be imported, tested,
linted and profiled. The notebook narrates and calls it.

### D-21 — Cache one parquet per (stage, year)
Keyed on `CODE_VERSION` + a hash of the *relevant* config subset + input mtimes, so changing an
unrelated constant does not invalidate a stage. Every stage is per-year and independent, so
widening `SNAPSHOT_YEARS` reuses what is already computed. Warm reruns: 0.3 s.

### D-23 — Validate the NFP id join per year; repair by geometry, flag, never fake
The `n_pixels` identity is checked for every year rather than assumed. Below
`NFP_ID_MIN_MATCH` (0.95) the year's labels are treated as stale and pixels are re-attributed by
point-in-polygon against the final perimeters (2021 only — see F11). The result is flagged with
`qc_nfp_repaired` / `qc_nfp_pixel_exact` rather than being presented as exact, and the hard
assertion still runs for the sound years.
**Rejected:** relaxing the assertion globally, which would have hidden the defect; and dropping
2021, which discards 63 pyroCb events over a real fire season.

---

## Version control

### D-22 — Git with LFS for the large binary outputs
Maps are ~27 MB of base64, which git deltas poorly. `*.html`, `*.xlsx` and `*.gpkg` go to LFS so
every historical map stays retrievable without bloating the pack. Source data (~4.3 GB of yearly
archives) is **not** tracked — it is a fixed external input, listed in `.gitignore`.
