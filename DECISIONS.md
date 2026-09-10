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

### D-48 — Matched fire v2: an identity-free nearby-snapshot centroid, alongside D-9–D-12
D-9–D-12 commit to one fire per event (`mergid`) and a strict before-onset snapshot. That
commitment is also a risk: if the atlas coordinate is off, or two fires burn shoulder-to-shoulder,
attaching to "the" fire can point at the wrong lobe of ground. `pyrocb/nearby.py` adds a second,
independent position that makes no such commitment: clip every fire polygon in the snapshot
(any mergid) to a `NEARBY_RADIUS_KM = 3` circle around the report, in whichever snapshot is
closest in time to onset — before OR after, since this is not a causal match — union what
survives, and take the centroid. No containment tier, no distance cap, no rescue passes:
nothing burning nearby means no v2 position (`fire2_lat`/`fire2_lon` = NaN), which is itself
informative rather than a gap to paper over.
**Kept separate from D-9–D-12, not a replacement.** v1 answers "which fire produced this
pyroCb" and is exported as `mergid`/`match_*`; v2 answers "what was on fire right here" and is
exported as `fire2_*`. The map's "pyroCb position" toggle (D-43) carries all three — as
reported, matched fire (v1), matched fire v2 — so a reviewer can compare them per event.

### D-49 — v2 averages NFP **detections**, not the snapshot hull, and does so unweighted
The first cut of D-48 took the area centroid of the Snapshot polygon clipped to the radius.
Measured over all years, that centroid fell on **unburned ground for 5.4 / 8.0 / 8.2 %** of events
at 3 / 5 / 10 km, by up to **3.1 km** — because a snapshot polygon is an alpha-shape hull over
*cumulative* burned area, and the centre of mass of a concave or two-lobed hull need not lie
inside it (36 % of 10 km hulls are multi-part). So it was the middle of a hull, not of a fire.
Replaced with the **unweighted mean of the centres of the NFP detections** at the nearest
overpass within the radius. NFP is incremental, so one step's detections are exactly the pixels
newly seen burning then. Re-measured: the miss is now **0.00 px (median) at 3 and 5 km**, 0.29 px
at 10 km (1 px = 375 m), and what remains is mostly the mean landing in a sub-pixel gap in the
lattice rather than off the fire; only 14 / 16 / 36 % are more than one pixel out. Coverage barely
moved (295→259, 338→319, 380→373 of 484) — the losses are events with cumulative burn nearby but
nothing actually detected at the closest overpass, which *should* be unmatched under this rule.
**Unweighted, not FRP-weighted:** this is the centre of the detected burning *area*; a cold flank
counts as much as the head. FRP-weighting would answer a different question and FRP is already
carried per event by `active.py`.
**Still visible at 10 km:** a third of cases are >1 px out because a 10 km circle routinely spans
two separate fronts and any mean lands between them. That is the radius telling you it is too
wide, not a defect of the estimator.

### D-70 — v3: a third, standalone matching pass with an explicit radius ladder, in its own notebook
D-9–D-12 (v1) commit to one fire via containment-then-cap; D-48/D-49 (v2) refuse to commit to any
identity at all. v3 is a third position, specified end-to-end by the user rather than derived:
within 10 km, take the time-closest snapshot (before or after onset, whichever is nearer); if it
has no active fire, try the 2nd-closest **only if** that candidate is within +/-12 h of onset;
still nothing, widen to 12.5 / 15 / 17.5 / 20 km, retrying the *same* two time-candidates at each
step (they depend only on onset time, not on the search radius); nothing within 20 km at either
candidate is `N/A`. `search_radius_km` (the ladder rung) and `match_dt_h` (signed, onset to the
matched snapshot's nominal time) are stored for every match, as specified.
Shipped as a **separate, deliberately non-package notebook**, `PCB_MATCHING_CLEAN.ipynb` +
`make_pcb_matching_clean_notebook.py` — unlike every other analysis in this repo (D-20: logic
lives in `pyrocb/`, the notebook narrates), here the matching and attribute-aggregation code is
written directly in the notebook's own cells, on request, so it can be read and checked line by
line without jumping to a package file. It reuses only long-established primitives (CRS
constants, `timebins.nominal_utc`/step parsing, `io.ensure_valid`) — nothing that rounds a
coordinate or a timestamp.
**Full Snapshot attribute set, not the `_SNAP_KEEP` subset v1 exports.** Every one of the 32 raw
columns (`fireid`..`mergedate`) is carried through for the matched fire, aggregated across its
merge-in-progress components with an explicit, documented rule per column (sum for extensive
quantities, max/min for per-component extremes, area-weighted mean rather than plain mean for
`meanFRP`/`avgspread`, `pixden` recomputed from summed n_pixels/farea rather than averaged) —
see the notebook's own step-4 table for the full rationale per column, including two
deviations from `pyrocb.io.read_snapshot`/`pyrocb.spread`'s existing conventions (`isignition`
maxed rather than summed; `pixden`/`meanFRP`/`avgspread` recomputed rather than plain-averaged).
**Found and fixed while building it:** `mergedate` is null for any component that never merged
into another fireid (a lone unmerged fire is 100% null; a merged group typically has exactly one
null member — the surviving/absorbing component, since nothing merges *into* it) — the first
aggregation pass assumed every value was present and crashed on `sorted()` over a mixed str/None
column; fixed to `dropna()` first, `NaT`/`None` when a fire has no merge history yet.
**Validated against real data before shipping**, not just code-reviewed: run standalone against
all 484 events (not executed as a notebook — this environment has no `nbconvert`/`nbclient`, so
every cell was run as one script and then, separately, executed cell-by-cell straight out of the
written `.ipynb` via `nbformat` to prove the shipped file itself runs end to end). Result: 395/484
matched (391 at the base 10 km rung, 12 needed the capped 2nd candidate, 4 needed to widen past
10 km), 89 genuinely `N/A` (real gaps: nothing burning within 20 km at either candidate time, up
to hundreds of km away — not a bug). Assertions confirmed by construction and re-checked against
independently re-derived values on a random sample: every match's distance is within its own
reported radius, every 2nd-candidate match is within the 12 h cap (max observed 11.12 h),
`snap_isactive`/`snap_invalid` are 1/0 for every matched row (consistent with Snapshot's "active
fires only" documentation), every matched geometry is valid, `snap_tst_utc <= snap_ted_utc`
always. Outputs (`pcb_matching_clean.gpkg`/`.tsv`) are checked-in run artefacts of that
validated pass, not hand-written.

### D-71 — v3 gets its own marker on the GOES cloud imagery: a green star, drawn last
`goes_cloud.py`'s comparison markers (D-69) were built for v2's (radius, method) grid — a dot per
combination, colour = radius, number = method. v3 (D-70) is a single position per event, not a
grid, so it does not fit that scheme and gets a distinct shape instead: a filled green 5-point
star, collected from optional `v3_lat`/`v3_lon` columns (missing/unmatched just skips it, same as
a missing v2 dot) and drawn **last of everything, including the reported crosshair** — the one
marker explicitly meant to be the most visible thing on the image, not competing for a colour.
It still goes through the same declutter/rosette pass as every other marker (D-69), so it is
never fully hidden behind an agreeing dot. Regenerating the images did **not** need a full
`Data exploration.ipynb` rebuild: `pipeline.run()`'s per-(stage, year) parquet caches (D-21) for
matching/bands/events/nearby are untouched by this change, so only `goes_cloud_cached`'s own
cache key needed a version bump to force new imagery — a ~6-year, 245-image GOES rebuild that
would otherwise cost a ~10-12 min cold run instead ran in well under two minutes.

### D-72 — v3 joins the map as a third "pyroCb position" option, mirrored on v2's UI but simplified for a single result
`fire_pcb_map.html`'s "pyroCb position" panel gets a third button, "Matched fire v3", alongside
"As reported" and "Matched fire v2". Structurally it copies v2's pattern end to end: a
`v3_lat`/`v3_lon` position that `repositionPcb()` moves markers onto, an `nb3` pane/canvas
layered just above `nb2`, and the same 3 extra layers as v2 (search circle, matched-fire polygon,
as-reported&rarr;matched arrow) with matching toggle checkboxes. The one deliberate difference:
v2 exposes a radius/method **picker** because it computes 12 different positions per event; v3 is
a single deterministic result per event (D-70's ladder already picked the radius), so there is
nothing to pick — its "search circle" toggle draws each event's own `v3_search_radius_km`
(whichever ladder rung actually matched it) rather than a global slider value.
**Data path:** `pcb_matching_clean.gpkg` (D-70's finished run artefact) is read once near the top
of `Data exploration.ipynb`, merged onto `events` by `pcb_id` (representative point of the
matched geometry -> `v3_lat`/`v3_lon`, same reasoning as D-49's centroid-vs-representative-point
fix), and its matched-fire polygon is simplified/reprojected into its own `nb3` GeoJSON layer,
exactly paralleling `nb2`'s "active fire used for centering" treatment. **Not recomputed here** —
`PCB_MATCHING_CLEAN.ipynb` stays the single source of truth for its own algorithm (D-73); the map
only visualises a finished result, so a v3 rerun and a map rebuild are two independent steps.
**Verified before shipping:** the generated JS was extracted from the built notebook and checked
with `node --check` (catches nothing at the Python-string level, since embedded JS is opaque text
to `ast.parse`), and every new `id`/`querySelector` pair was cross-checked 1:1 by grep — both
before spending the ten-plus minutes on a real rebuild to find out.

### D-73 — PCB_MATCHING_CLEAN.ipynb is standalone: no `pyrocb` import, everything inlined
Every other analysis notebook in this repo imports `pyrocb` (D-20: logic lives in the package).
`PCB_MATCHING_CLEAN.ipynb` deliberately breaks that pattern a second way (on top of D-70's
"logic lives IN the notebook, not behind a function call"): it imports no first-party module at
all, only `numpy`/`pandas`/`geopandas`/`shapely`. The handful of `pyrocb` primitives it still
needs — the two CRS codes, the Snapshot path template, the half-day step schedule
(`parse_step`/`make_step`/`norm_lon`/`nominal_utc`), geometry repair (`ensure_valid`), and the PCB
loader (`read_pcb`) — are copied into the notebook's own "1b" cell **verbatim**, not
reimplemented, so its numbers stay identical to the rest of the project while the file itself
has nothing else to depend on. `read_pcb`'s parquet-cache wrapper (`config.cached`) is dropped —
484 rows from one Excel file is cheap enough to just recompute every run, and pulling in the
cache layer would have reintroduced a `pyrocb` dependency for no benefit.

### D-74 — v3 stores the burn INSIDE the search buffer, not the whole matched fire
The first v3 build attached the matched fire's entire snapshot polygon to the event, and put the
map marker on that polygon's representative point. On small fires this is invisible (the
footprint is barely bigger than the buffer). On big fires it is nonsense, because a FEDS snapshot
polygon is **cumulative season-to-date burned area**, not the burning front: matched footprints
here span up to **298 km**, so 100 of 395 events had their v3 marker plotted **>20 km** from their
own pyroCb and 37 of them **>50 km**. Worst case `PCB_2023_0602_2020` / 2023-G90288 — the pyroCb
sits *inside* the fire (`match_dist_km` = 0) and the marker landed **150 km away**, on a different
lobe, weeks of spread from the event.

The ladder was never wrong: for both G90288 events the matched fire is the *only* fire within
10 km. What was wrong is what got attached afterwards. So the buffer is now treated as defining
the area, not just testing for a hit: the stored geometry is every fire polygon reaching into the
resolved R km buffer, unioned and **clipped to that buffer** (G90288: 3,610 km² whole footprint →
104 km² of burn actually near the pyroCb). Any position derived from it is then inside R km of
the report by construction — asserted in the notebook's checks cell, so this regression cannot
come back silently.
**Identity-free inside the buffer**, matching D-49's reasoning for v2: if two fires reach in, both
are kept and their mergids recorded (`n_fires_in_buffer`, `mergids_in_buffer`) rather than
dropping one because the other was a few metres nearer. The nearest fire still supplies the
`snap_*` identity.
**Nothing is discarded:** the whole-fire footprint ships as a second GeoPackage layer,
`matched_fire_full`, beside the default `pcb_active_burn`.
**The caveat this creates, stated in the notebook:** `snap_*` attributes describe the whole
matched fire (the only scale FEDS publishes them at) while `geometry`/`clip_area_km2` describe
only the buffer interior. They must not be divided into each other.

### D-75 — Area-weighted snapshot means skip missing components instead of returning NaN
`avgspread` is empty in ~25 % of raw Snapshot rows. `np.average` propagates a single NaN across
the whole weighted mean, so one blank component emptied `snap_avgspread_m` for the entire fire —
which hit **79-100 % of multi-component fires** against **8 % of single-component ones**. The
missingness is therefore strongly correlated with fire size, i.e. exactly the fires the analysis
is about, and it looked like a data gap rather than an arithmetic artefact. The weighted mean now
runs over the components that have a value, dropping the others from both numerator and
denominator; all-missing still yields NaN, and `snap_n_comp_with_spread` reports how many
components actually backed the number so a thin mean is visible rather than implied.
**Rejected:** filling blanks with zero — that is a measured "no spread" claim the source never
made, and it would drag every big fire's mean toward zero.

### D-76 — v3's position comes from the VIIRS pixels, not from any centroid of the Snapshot polygon
Three separate failures, found by pulling the thread on "more pixels do not pull the centroid".

1. **The point was weighted by nothing.** `representative_point()` is GEOS *PointOnSurface*: it
   guarantees only that the point is inside the polygon. 140 of 395 clipped burns are multi-part
   (up to 7 lobes) — precisely where a bigger lobe should win — and it sat a median 0.86 km, worst
   case **11.9 km**, from the area-weighted centre.
2. **The credited fire was the nearest, not the dominant one.** A 0.11 km² speck 4.0 km away beat
   a 15.1 km² fire 6.2 km away (`PCB_2023_0630_0620`), and a 0.42 km² sliver beat a 41.7 km² fire
   (`PCB_2022_0708_2210`) — attaching the speck's FRP, area and duration to the pyroCb. Now the
   fire with the most burn **inside the buffer** is credited, ties broken on distance. The ladder
   decision (did this match, at what radius) still uses the nearest fire, so search and credit stay
   separable and both are recorded.
3. **The root cause: no centroid of that polygon can work.** A Snapshot polygon is *cumulative*
   burned area. Measured against v2's detections in the *same* circle at the *same* overpass, the
   clipped polygon is a median **3.8×** larger (p75 7.4×, p90 21×) — median 60.8 km² against
   14.2 km². It carries no pixel density at all, and most of its mass is ground that went cold days
   ago, so its centre of mass tracks old burn rather than the front.

So the position is now the plain mean of the **VIIRS detection centres** in the ladder buffer at the
matched overpass — unweighted over centres, which *is* "more pixels pull it closer", and the same
call D-49 made for v2 (centre of burning area; a cold flank counts as much as the head). Validation:
on `VQ006-23` the marker went from 6.59 km off to **0.40 km**, landing on v2's independently-computed
answer to two decimals; `VQ009-23` 5.44 → 3.44 km against v2's 3.5 km. Two separately-implemented
estimators agreeing is worth more than either one's internal checks.
**Then the analysis window is re-cut around that corrected centre** (`ANALYSIS_RADIUS_KM` = 5 km),
not around the report — the atlas coordinate is rounded to 0.1° (~11 km), so a window centred on it
is not a window centred on the fire. Median active burn 11.5 km² / 68 px, against 59.5 km² for the
cumulative clip it replaces.
**Fallback, flagged not hidden:** 9 of 395 events have no detections in range (cloud, or burning but
not newly detected that overpass). Those keep the polygon's area-weighted centroid, snapped onto the
burn, and say so in `v3_pos_source`.
**Layers, so nothing is lost:** `pcb_active_burn` (the analysis polygon), `snapshot_clip` (the
cumulative clip), `pcb_active_burn_v4`, `matched_fire_full`.
**One source of truth:** the notebook now ships `v3_lat`/`v3_lon` and the map reads them. The map
used to derive its own position, which is exactly how it displayed a point the notebook never
computed.

### D-77 — v4: four overpasses, per-pixel time decay, normalised per step
v3 sees a single 12-hourly look. v4 pools the **4 time-closest** overpasses and weights each pixel
`exp(-|Δt| / 12 h)` on its **real `datetime_utc`** — NFP carries it, so the decay is exact rather
than quantised to half-days — then **divides by that step's own pixel count**.
The normalisation is the load-bearing part. Without it a step's influence scales with how much it
happened to detect, so a massive burn 24 h later simply outvotes a small one an hour before onset by
count. With it each step contributes its own centroid at a weight set by timing alone, which is the
stated intent: a few pixels just before onset can outweigh a huge area a day later. Equivalently,
v4 is the time-weighted mean of the four per-step centroids.
τ = 12 h (one overpass interval) gives a pixel 1 h from onset ~6.6× the pull of one 24 h away.
Coverage 394/395 events, median 4 steps pooled, median 24.9 km² / 193 px — larger than v3's window
because it pools four looks. v4 sits a median 1.01 km from v3 (max 9.59 km), so it is a genuine
second opinion rather than a restatement.
**Shipped alongside v3, not replacing it**, as a fourth "pyroCb position" option on the map with its
own overlays, in amber against v3's green.

### D-78 — v4 gets its own star on the GOES cloud imagery too, purple, alongside v3's green
Same reasoning as D-70/D-71 for v3: a single-position match doesn't fit the v2 (radius, method) dot
grid, so it gets its own shape. `goes_cloud.py`'s marker code was generalised from one hardcoded
`v3_point` to a `stars` dict keyed `"v3"/"v4"`, both batched through the same projection/declutter
pipeline and drawn last (v4 on top of v3, both on top of the reported crosshair), so neither can be
silently hidden by an agreeing v2 dot or by each other. Colour reuses `RADIUS_COLOR`'s existing
purple (168, 85, 247) — v4's own legend line disambiguates it from the v2 "10 km" dot, and the
shapes differ (star vs. filled circle) regardless.
`pipeline.goes_cloud_cached`'s position fingerprint (D-74's fix for stale imagery surviving a moved
star) now hashes v3 **and** v4's coordinates together, so a future change to either invalidates every
affected image — the exact bug D-74 fixed for v3 alone would otherwise still be open for v4.

### D-79 — Supervisor demo pass: hide, don't delete, behind a single "Show advanced options" toggle
For a clean first look, the map now opens with only the features a supervisor needs visible:
Years, Layers, Colour by, pyroCb position (v4 only), Land cover in the info panel. Everything
else — Animate, Regions, pyroCb click (panel/popup + zoom-to-select), the "As reported"/v2/v3
position options, and most of the animation panel (Static mode, the ±24 h window scope,
Sliding style + window-size) — still exists, still works, just starts hidden.
Mechanism: every element that should start hidden gets `class="extOnly"`; a single
`SHOW_EXTENDED` flag + `applyExtVisibility(root)` (`querySelectorAll(".extOnly").forEach(el =>
el.hidden = !SHOW_EXTENDED)`) toggles them all at once from one checkbox at the top of the
panel. Re-run after every `openInfoPanel()` build (that DOM is rebuilt per click) and after the
control panel is first attached. **Rejected:** removing the code — this is a working-session
map that gets iterated on again; deleting v2/v3/static/sliding would mean re-writing them
later. **Rejected:** per-feature toggles — one switch is what "clean for a demo" means; a wall
of checkboxes defeats the purpose.

### D-80 — The analysis-polygon map layer is fully removed, not hidden
Unlike everything else in D-79, the analysis-polygon layer (`AnalysisLayer`, `rebuildAnalysis`,
its pane, its panel section, its legend line, the `ANALYSIS_WIN`/`__ANLWIN__` template
plumbing) is deleted outright rather than hidden behind the toggle. It was a reviewer-facing
diagnostic for the covariate-sampling window, unrelated to the v3/v4 matched-fire re-centred
windows the map is actually about, and not something worth explaining in a demo. `pyrocb/
pipeline.py`'s `ACTIVE_TABLE_WINDOWS` output, `pyrocb_active_areas.gpkg` and the attribute
table are untouched — this only removes the map's own extra canvas layer.

### D-81 — pyroCb position defaults to v4; the info panel opens on Animate, scoped to the fire
`PCB_POS` now starts `"fire4"` (was `"report"`) and `PANEL_MODE` starts `"animate"` (was
`"static"`) — v4 is the map's best current match and the animation is the more informative
view of an event, so those are what a first-time viewer sees without touching anything.

### D-82 — Info-panel popup shows only id / fire name / datetime up front; everything else is one click away
`openInfoPanel()`'s header used to be a full `<table>` of provenance fields (region, mergid,
match method, v2/v3/v4 summaries, satellite + GOES imagery) shown unconditionally. That table,
plus the satellite/GOES imagery, now sits inside a closed-by-default `<details class=
"moreDetails">` under the header — a supervisor sees the pyroCb id, fire name and UTC datetime
immediately, and can expand for the rest. Land cover is pulled *out* of that dropdown and shown
directly below it (top 5 classes now, was top 3 — `landcover.active_fractions` takes an
`n_top` parameter, default 3, called with `n_top=5` from `pipeline.landcover_cached`), with the
"6 h pre-onset" window label dropped from the panel text (the window itself is unchanged;
`LC_ACTIVE_WINDOW` is still `(6, 0)` with the blow-up fallback from D-37) — the fixed 6 h
disclaimer.

### D-83 — Animation defaults to Whole fire / Cumulative 1 week; ±24 h window is an extended feature; loop range sliders added
The per-event animation scope picker used to default to a fixed ±24 h window. It now offers
"Whole fire" (unchanged) and a new "Cumulative 1 week" scope by default, with the old ±24 h
window pushed behind the advanced toggle. "Cumulative 1 week" reuses the fire's own whole-fire
grid and time axis (same `fsim.g`, same t=0-at-fire-start absolute anchor as "Whole fire") but
bounds `t0Of()`/`spanOf()` to a 168 h window centred on the *clicked event's own onset*
(`myFireH()`), clipped to `[0, fsim.dur_h]` so it never runs off either end of the fire's
observed record. `usesFire()` generalises the old `animScope === "fire" && fsim` guard to cover
both fire-grid scopes; every scattered check in `render()`/`drawTicks()`/`drawIgnTicks()`/
`renderPops()`/`buildSnap()` now calls it instead of comparing scope strings directly, and the
per-tick x-position math was corrected to subtract the scope's own `t0Of()` (it happened to be
0 for "fire", so this only mattered once "week" existed). New loop start/stop range sliders
(`loopLoU`/`loopHiU`, fractions of the current scope's span) let a reviewer scrub the playback
loop down to the interesting stretch instead of always looping the whole scope; `setScope()`
resets the range to `[0,1]` on every scope change, since the span itself just moved.

### D-84 — GOES cloud imagery: earliest scene from onset (0 h), not 2 h after, with a hard read-attempt cap
`GOES_POST_WINDOW_H` moves from `(2.0, 12.0)` to `(0.0, 12.0)` — the earliest usable scene at
or after onset, rather than giving the storm 2 h to mature first, matches what the map now
promises ("nearest scene to onset"). Starting the search at 0 h means a dusk/night onset can
walk through many darkening candidates before `GOES_DAYLIGHT_LOCAL_H`/
`GOES_MIN_DAYLIGHT_REFL` reject them, each a real network fetch; `GOES_MAX_READ_ATTEMPTS = 10`
caps the candidate-scene loop in `_build_image()` so one bad event can't blow up the whole
rebuild's runtime. The scene that wins when found well inside the cap is unaffected.

---

## The analysis polygon

**Removed from the map in the presentation build (D-83) — this section is kept as the historical
record of why it was built the way it was, not as a description of anything still on screen.**

### D-50 — Covariates are sampled over the NFP active-burn window, and the map draws it from the cell grid
The geometry every per-event variable (ERA5 temperature/humidity/wind, terrain slope, land cover,
fire intensity) is sampled over is `active.py`'s localised NFP pixel-footprint union for a
(pre, post) window around onset — the same object that already backs the land-cover fractions
(D-34) — with `C.ANALYSIS_WINDOW` naming the default (±12 h).
**Rejected: shipping it as a polygon layer.** The exact polygons are 10.3 MB gzipped for 2023
alone across the six table windows; a single window over all years is still 1.8–4.7 MB. Reducing
them with a ±400 m buffer-close gets to 0.9 MB but **inflates area by +5.0 km² median and up to
53 km² (p95)** by filling the lattice gaps — unacceptable when the drawn shape *is* the sampling
geometry. Simplify-only at 150 m is honest (median −0.09 km²) but does not pay for itself.
**Chosen:** draw it client-side from `byPcbPix`, the per-event 375 m detection grid the map
already ships for the window slider (0.10 MB for all of 2023). Cell times are hours relative to
that event's own onset, so any window is a contiguous index range — **zero new payload, no
simplification at all, and every window in `ACTIVE_WINDOWS` is reachable from a slider.**
This reverses D-35 (which dropped a static active-burn layer) now that there is a reason to want
one permanently on: it is the model's sampling footprint, not decoration.

---

## Pre/post satellite imagery

### D-51 — A circle sized as a PIXEL BUDGET, not a fixed ground resolution
"Around the whole perimeter" taken literally means the enclosing radius (centroid to farthest
boundary point) sets the circle, and that varies 2–160 km across 2023's 69 matched fires (median
22 km, `sub.encl_km.describe()`). A fixed ground resolution (e.g. native 10 m) would make the
biggest fire's image 30,000+ px on a side. `SAT_OUT_PX = 1024` is a pixel cap instead: small
fires render near-native, huge ones downsample — still enough to see the scar's shape, which is
the point, not survey-grade detail.

### D-52 — Candidates ranked by closeness to the anchor date, not by cloud cover
First cut sorted STAC results by ascending `eo:cloud_cover` across the whole search window.
Measured failure: a fire ending 2023-07-22 at 57.7°N, searched `ted+10d … ted+120d`, picked its
globally-clearest scene on **2023-11-10** — by which point the ground was already snow-covered,
which erases the burn scar completely (visually confirmed). The window has to be wide because
weather is unpredictable, but the *clearest* day in a wide window can easily be the *wrong
season*. Fixed by sorting candidates by `|scene_date − anchor_date|` instead (anchor = ignition
for pre, official end for post) and taking the closest one under the cloud ceiling. Re-measured on
the same fire: pre landed the day *before* ignition, post landed 14 days after end, both 100 %
coverage. **Rejected:** narrowing the window instead — boreal fires that end in October have no
snow-free window far enough out to reliably beat cloud cover, so proximity-first is the only
rule that degrades gracefully (picks the least-bad clear-ish scene) rather than catastrophically
(picks a technically-clear but seasonally-useless one).

### D-53 — Images are JPEG files on disk, keyed by mergid, not embedded in the map HTML
A circle-cropped 1024 px true-colour PNG is 1.5–2 MB. At ~70 fires × 2 phases that is 200+ MB —
embedding it the way every other layer works (gzip + base64 inside `fire_pcb_map.html`) would
multiply the map's size 10×, for a feature most sessions never open. Flattening the transparent
circle edge onto a solid background and re-encoding as JPEG (no georeferencing needed — this is
a standalone photo, not a map overlay, so the earlier EPSG:3857 reprojection step was dropped
too) brings a typical image to 100–200 KB. Files live under `sat_images/{year}/{mergid}_{pre,post}.jpg`,
next to `fire_pcb_map.html`; the map ships only a small per-fire metadata dict (path, date, cloud
%, coverage, scene count) and loads the actual image lazily via a plain `<img src>` the moment a
viewer opens that fire's dropdown in the info panel — never before.
**Consequence:** the map is no longer a single self-contained file when this layer has data
(2023 only, `C.SAT_YEARS`). That is an explicit trade, not an oversight — the alternative was not
"one bigger file", it was "an unusable one".

### D-54 — Coverage is reported, not hidden, when a fire is bigger than the mosaic budget
A handful of 2023 complexes exceed 100 km enclosing radius (max 160 km, mergid 7783 = 14,966 km²)
— a circle wider than a dozen Sentinel-2 tiles. `SAT_MAX_TILES = 20` (raised from an initial 6
after measuring 18–32 % coverage on the worst case) still only reaches ~45–50 % coverage for the
single largest fire; the rest of the circle stays transparent/background rather than being
silently cropped or stretched to look complete. Every image ships its `coverage` fraction and
`n_scenes` count, shown in the map's caption, so a reviewer can see "this one is incomplete"
rather than trusting a partial mosaic as if it were the whole picture. Small/medium fires (the
large majority) are unaffected: the mosaic loop exits as soon as one or two tiles fill 100 %, so
raising the tile cap costs those fires nothing.

### D-55 — v2 map simplified to one radius, and "matched fire v1" dropped from the map
Three changes made together after reviewing the shipped map:
1. **v2's mean is already density-weighted, by construction, with no extra weighting needed.**
   The position is `x.mean(), y.mean()` over the raw NFP detection points inside the radius -
   each detection is one term in the sum, so a cluster of 10 detections pulls the mean 10x as
   hard as 1 isolated detection on the opposite side. Re-verified two ways: (a) a synthetic
   9-vs-1 point check confirming the mean lands at the analytically correct 9:1 split, not the
   midpoint; (b) the shipped position, decoded from `fire_pcb_map.html`, reproduced by
   recomputing `x.mean()/y.mean()` straight from the raw NFP files - agreement to 0.44 m. There
   is no separate "weight by density" step to add on top; doing so would double-count what plain
   averaging already does.
2. **Radius cut to 10 km only** (D-49's 3/5/10 km selector removed). Measured coverage at 3 km
   was 53 % of events, 5 km 70 % - both dominated by the atlas's own position error (±5.6 km N-S)
   rather than by anything the method could fix, so the smaller radii could not be trusted and the
   selector was pure clutter. 10 km reaches 77 %.
3. **"Matched fire v1" (the mergid active-burn centroid, D-43) removed from the map's position
   toggle.** It answered the same question as v2 - "where did this fire actually burn near onset"
   - using a NARROWER, mergid-scoped input (D-29's own-cluster localisation) where v2 uses a
   plain radius. Once v2 moved to real NFP detections (D-49) the two were redundant, and shipping
   both cost the reviewer a three-way comparison for very little contrast. v1's `fire_lat`/
   `fire_lon` computation is deleted from the map notebook entirely (D-43's decision record is
   left in place as history, not edited, per this file's own convention). Nothing changes for
   `active.py`/the attribute table - this is a map-payload change only.

**What was NOT the cause of the "some points don't move / show no circle" symptom that prompted
this review:** JSON validity. `NaN` fields (the majority of events, since most have no v2 match at
a given radius) were checked end-to-end - `geopandas.to_json()` already emits `null` for them, not
a bare `NaN` token that would break `JSON.parse` - so that was never the mechanism. The far more
likely explanation is coverage itself: at the old default radius (3 km), only 53 % of events had
ANY v2 position, so roughly half of all pyroCb markers legitimately sat still and drew no circle -
expected behaviour for "nothing detected nearby," easily misread as broken. 10-km-only raises
that to 77 %; the remaining 23 % are still expected to stay put by design, not a bug.

### D-57 — v2 tries the next-closest overpass when the closest one imaged nothing here
Investigating "why did 48 of 159 events not move at all?" turned up two separate things.
**A bug.** The bracket search that was supposed to find the nearest available overpass
"before OR after" onset had a forward walk that started AT the anchor step instead of one step
past it. `timebins._walk_back` deliberately returns its own start step when that step is
available (its contract is "at-or-before"), so with both walks starting at the anchor, `before`
and `after` collapsed to the same step whenever the anchor had data anywhere in the dataset —
the common case. The forward direction was never actually explored, and the whole feature
silently degraded to "just use `step_containing`". Fixed by starting the forward walk at
`next_step(anchor)`. Measured: **17 of 484 events move to a closer overpass, 0 to a farther one**
(median time-to-overpass 7.9 h → 5.1 h among those; best case 8.9 h closer).
**A design gap the bug was hiding.** A step is "available" if it detected fire *somewhere* — but
VIIRS swaths do not cover every location on every pass, so the nearest-in-time overpass can be a
regional blank over this particular fire. Type case: PCB_2023_0527_0720, whose chosen overpass had
its nearest detection **262 km away**, while the next overpass (3 h *after* onset rather than 10 h
before) carried 1141 detections whose closest sat **190 m** from the pyroCb. So the step choice now
falls through: try the closest overpass, and if it imaged nothing inside the widest radius, try the
next-closest, up to `NEARBY_MAX_STEPS = 3`. Same tiered-rescue shape as D-11.
**Combined effect: 373 → 391 of 484 events (77.1 % → 80.8 %).** The remaining 93 have nothing
within 10 km across ~3 days of overpasses — genuinely no detected fire near the reported point,
which is information, not a failure: they keep their as-reported position by design.
**Cap justified by measurement, not taste:** the 2nd-closest overpass recovers 17 events, the 3rd
recovers 1, the 4th through 6th recover 0.

### D-58 — The v2 search circle is drawn around the REPORTED point, not around the result
It was first built around the v2 centroid. That is backwards: `nearby.py` balls the KD-tree
around the **atlas lat/lon**, so the reported point is the centre of the search and the centroid
is its output. Drawing the circle around the output implies the detections were gathered around
the answer — and it hides the thing the circle is actually useful for, which is seeing whether
the detections that moved the point sat near the edge of the search area (i.e. whether a slightly
larger radius would have changed the answer). Together with the D-56 arrow, the honest picture is:
circle fixed on the report, arrow running from the report to wherever the detections pulled it.
Implementation note: read the origin from `feature.geometry.coordinates`, not from the marker —
`PCB_POS` moves markers with `setLatLng`, which does not write back to the feature, so the
geometry always still holds the as-reported coordinates.

### D-56 — Displacement arrow: as-reported to matched-fire-v2, in the report's own frame
A dashed circle and a teal detection footprint show *where* v2 landed; neither shows *how far*
or *which direction* from the atlas point. The arrow toggle draws both in one glance: a
polyline plus a CSS-rotated triangle at the v2 end, bearing computed with the standard
initial-bearing formula (accurate; the arrows are short enough that great-circle vs rhumb
line is not a visible difference). No new geometry payload - both endpoints are already-shipped
per-event fields, so this is pure client-side drawing, gated to `PCB_POS === "fire2"` like the
circle/footprint toggles it sits next to.

### D-59 — Satellite imagery: drop the pre-fire image, ship post-fire for every year
D-53 shipped pre+post pairs scoped to 2023 only (`SAT_YEARS = (2023,)`), one merged fire at a
time. Revisited once the method (D-52 date-proximity search) and delivery (D-53 external JPEG
files) were both proven out: **pre-fire dropped, post-fire extended to every resolved year**
(`SAT_YEARS = None`). Rationale: the map already carries the "before" picture as vector layers —
the final perimeter, the growth bands, the fire's own ignitions — laid over a live basemap a
reviewer can already compare against; a static pre-fire photo mostly repeated that. The postfire
photo is the one thing nothing else on the map shows (the actual burn scar), and it is the one
that answers "did this really burn where the atlas/FEDS says it did" - worth having for every
matched fire, not just one year's. Cost-neutral to fetch: one search+mosaic per fire either way,
so dropping pre roughly offsets extending post from 69 fires to ~275.

### D-61 — PyroCb cloud imagery needs a different satellite entirely, and only covers ~half the atlas
`satellite.py` answers "did this burn here"; it cannot answer "what did the pyroCb actually look
like", because Sentinel-2 revisits any point every ~5 days and a pyroCb exists for hours - the
probability of a pass landing near one specific event's timestamp is essentially zero. That needs
a geostationary satellite (same hemisphere imaged every 5-15 min), which is a structurally
different, smaller-coverage source: `pyrocb/goes_cloud.py` uses GOES-16/17/18/19 (`goes-cmi` on
Planetary Computer), and two hard limits come with it, neither fixable in code -
**no geostationary coverage of Siberia** on Planetary Computer (Northern Asia's ~24% of events
excluded outright) and **GOES-16 not operational before Dec 2017** (2013-2017 excluded outright).
Ceiling: ~252 of 484 events (2018-2023, Canada/Alaska).
**Band: visible (C02), not infrared, reversing the first instinct.** IR (10.3 µm "clean window")
is the textbook way to spot a pyroCb's cold overshooting top and works at night, so it was the
obvious pick - and was tried first. Measured against it: at the 20 km radius asked for, IR's 2 km
pixels give an ~11×22 px native window, which upsamples to a blur with no discernible structure.
The single visible band at 0.5 km gives ~4x the linear pixel density and reads as actual cloud
texture. Switched despite losing night coverage, because that loss turned out to be nearly free:
98% of the 252 GOES-eligible events fall in local daytime, because pyroCb convection is itself
driven by daytime heating - the same physics behind the method is why the "downside" barely
applies to this specific dataset. No colour synthesis needed either way - one band, grayscale,
percentile-stretched.
**Search: closest scene in time, unconditionally** - there is no cloud-cover property to rank by
here (cloud is the subject, not noise to minimise). Measured on a real event: the chosen scene
landed 0.7 minutes from the pyroCb's own onset.
**Centred on the v2 matched-fire position** (`fire2_r10_lat/lon`) when one exists, the as-reported
point otherwise - the same fallback the map's own position toggle already uses, so this needed no
new judgement call about which coordinate counts as "the fire".
**Delivery** follows D-53/D-59 exactly: JPEG file per event under `config.GOES_IMG_DIR`, only
small metadata embedded in the map payload, loaded lazily per pcb_id (not per mergid - this is
about the moment, not the fire as a whole).

### D-62 — GOES imagery moved to "earliest daylight scene 2-12 h after onset", widened to 30 km
D-61 shipped "closest scene to onset, +/-30 min" at a 20 km radius. Revised on request, for two
independent reasons:
1. **Timing.** The exact minute the atlas timestamps as onset is often too early for the pyroCb to
   have built into a distinct, recognisable cloud yet, and a meaningful share of onsets fall at
   local dusk regardless of season. Searching `GOES_POST_WINDOW_H = (2, 12)` gives the storm time
   to mature, and restricting to `GOES_DAYLIGHT_LOCAL_H` (checked via `local_solar_hour`, before
   any fetch) plus a post-read reflectance floor (`GOES_MIN_DAYLIGHT_REFL`, for the dawn/dusk edge
   the coarse hour band lets through) keeps every scene actually readable. Candidates are taken
   **earliest first**, not clearest or closest-to-a-preference: "2-12 h after onset" should mean
   what it says, not silently drift to whichever hour in the window happened to look best - the
   same reasoning as D-52's date-proximity search for the post-fire imagery.
2. **Radius**: 20 km -> 30 km, for more context around the storm now that the image is no longer
   pinned to the literal instant of onset - by 2-12 h later the cloud (and, in the metadata, `dt_h`
   itself) has had room to move.
`dt_h` in the metadata replaces the old `dt_min`: now routinely 2-12, not a fraction of a minute.

### D-63 — The analysis-polygon map layer is cut to 5 km of the event's own v2 centroid
D-50 shipped the analysis polygon as the full NFP active-burn window, localised to the pyroCb's own
fire cluster (D-29) but otherwise unbounded - for a large complex that cluster can still span tens
of km, most of it not actually relevant to any one pyroCb sitting near one edge of it. Cut on
request to `ANL_RADIUS_KM = 5` of that event's own `fire2_r10_lat/lon` (D-49's identity-free
detection mean), computed client-side in EPSG:3857 (the burn grid's native storage CRS - D-40) with
the standard Mercator-inflation correction (`/cos(lat)`) so 5 km stays true on the ground rather
than on the map. An event with no v2 position has no centroid to cut around and shows nothing.
**Scope: this map layer only.** `pyrocb_active_areas.gpkg` and the attribute table's land-cover
fractions still use the full, uncut fire-cluster geometry (D-34) - narrowing THAT would mean
re-deriving every downstream covariate for every event, a much heavier change than a map filter,
and is deliberately left as a later decision rather than done as a side effect of a map request.
The map's own copy says so explicitly, so the two are never mistaken for describing the same
geometry.

### D-64 — Land-cover PNG overlay resampled to 375 m, snapped to the fire-activity grid
Was 120 m, resampled in the source CRS (EPSG:3979) at each fire's own bounding-box-relative
origin, then reprojected to EPSG:3857 with an auto-computed transform (`calculate_default_transform`)
sized to fit that one image. Two consequences: no two fires' land-cover grids lined up with
anything, and the resolution didn't match `burngrid.py`'s 375 m fire-activity grid the whole-fire
animation draws, so overlaying the two never quite agreed pixel-for-pixel.
Fixed on both counts. `LC_PNG_CELL_M` is now literally `BURN_CELL_M` (375 m, one VIIRS pixel) - not
independently chosen and coincidentally similar, the same number for the same reason. The output
grid is snapped the way `burngrid.encode` grids the fire-activity layer: cell edges fall on exact
multiples of the cell size measured from the EPSG:3857 Mercator origin (`floor(x/cell)*cell`), not
from each fire's own bbox. Because that's an infinite tiling anchored at a fixed point, not a
window-relative grid, two rasters built with the same cell size share cell boundaries regardless of
which piece of the world either one covers - a land-cover cell and a fire-activity cell at the same
fire now tile against each other exactly, not merely overlap. Implementation: one direct
`rasterio.warp.reproject` from a decimated (not full 30 m native) read straight to the
Mercator-snapped destination transform, replacing the old two-step resample-then-reproject - fewer
moving parts as well as more correct. `LC_PNG_CAP_PX` still widens the cell for a continent-sized
fire, but the widened cell still snaps to the same origin.
**Scope: this map layer only**, same caveat as D-63 - the land-cover FRACTION columns in the
attribute table (`lc_top1_cls` etc.) still sample at native 30 m and are untouched by this.

### D-65 — Satellite (Sentinel-2 post-fire) imagery switched off, kept as a one-line restore
`SAT_YEARS` set to `()` rather than deleting or gutting `satellite.py`/its pipeline wiring. The
existing year-gate in `satellite_cached` (`year not in C.SAT_YEARS`) already treats an empty tuple
as "skip every year" with no code change needed elsewhere, and every image and cache entry already
on disk is untouched - re-enabling later is `SAT_YEARS = None` (every year) or a specific tuple,
not a rebuild-from-scratch. Left off purely to keep a same-day build fast: with it off, this build
touches no Sentinel-2 network calls at all, only the (already-cached) GOES imagery and the
land-cover regrid, which is local computation.

### D-66 — v2 ships 4 snapshot-selection strategies x 3 radii, side by side, for comparison
The whole point of v2 (D-49) was always "a second opinion", but until now it committed to ONE
opinion: closest overpass, one radius. With no ground truth to validate against, the honest move
is to stop picking a winner and ship several genuinely different answers so a reviewer can compare
them directly - not four cosmetic variations, four different *causal claims* about which
overpass(es) should count as "near" a pyroCb in time:
- `closest` - nearest in time, no direction preference (the original D-49/D-52 behaviour, kept
  as one of four rather than replaced).
- `before` - nearest AT-OR-BEFORE onset; only reaches for an after-onset overpass once every
  before-onset option is exhausted. This is closer to `matching.py`'s causal stance (D-9-12) than
  `closest` is, applied without the mergid commitment.
- `before_after` - nearest before AND nearest after, POOLED - captures the transition across
  onset rather than picking a side.
- `before2_after2` - nearest two on each side, pooled - the widest net, useful for judging whether
  a wider time window changes the answer much.
Radii went from 10 km alone to (10, 15, 20) for the same reason - a comparison needs more than
one point per axis. 3/5 km stay dropped (D-55: dominated by the atlas's own position error).
**Implementation.** The old single-step-per-event search became a per-step KD-tree cache (one tree
per overpass, built once and shared across every event/method/radius that touches it, since a
step's tree does not depend on who's asking) plus a per-method "tiers of step-groups" list -
`before_after`/`before2_after2` pool several steps' detections into one search before computing
the centroid, rather than picking a single step. Regression-tested against the exact 9-vs-1
density-weighting check from D-55 to confirm the centroid math survived the rewrite unchanged, and
against real 2023 data to confirm the four methods diverge exactly where they should (9 events
where `closest` and `before` pick different steps, all because `closest` took a nearer AFTER-onset
step that `before` correctly held out for) and pool detections correctly (median n_pixels roughly
doubles `closest` -> `before_after` -> `before2_after2`, tracking the overpass count).
**Scope.** All 12 (radius, method) combinations ship to both the attribute table (96
`fire2_r{R}_{M}_*` columns) and the map (a radius selector next to a new method selector, both
driving the same `v2Key()` the position toggle, circle, footprint and arrow already used).

### D-67 — nb2 map geometry needed real simplification once D-66 landed
D-66's before2_after2 pools up to 4 overpasses' worth of separate VIIRS pixel rectangles, across
12 (radius, method) combos instead of 1 - shipped un-smoothed (the earlier "small enough,
radius-capped" assumption held for one 10 km combo, not for this), the first full build came to
35.1 MB, up from 16.3 MB, almost entirely from the nb2 "active fire used for centering" polygon
layer (measured: 87.6 MB raw across all years, the largest single payload key by far). Plain
`simplify()` barely helped (median vertices 1264->699 even at a 150 m tolerance) because the shape
IS hundreds of disjoint small rectangles, not a curve with redundant vertices to remove. Fixed with
the SAME close-then-simplify treatment `smooth_for_map` already applies to the other pixel-lattice
layers: buffer +/-`NB2_CLOSE_M` (200 m) to merge nearby detections into solid blobs, then
`simplify(NB2_SIMPLIFY_M)` (30 m) those. Cut 2023's nb2 payload from its share of 87.6 MB raw down
to 9.14 MB (0.83 MB gzipped), zero polygons lost to the process; final map 19.2 MB. Map-display
geometry only - `area_km2`/`n_pixels` in the data are computed from the un-smoothed polygon in
`pyrocb/nearby.py`, so only the drawn shape changed, never the reported numbers.

### D-68 — The static active-burn window snaps to real overpasses, not fixed hour marks
Was a pair of 0-24 h sliders in 6 h steps, driving a precomputed area lookup at the FIXED grid
`config.ACTIVE_STEPS_H` (0/6/12/18/24). Real overpasses do not land on round hour marks - they
follow the AM/PM local-solar schedule (`timebins.py`), which varies by longitude - so "12 h before"
was frequently not actually the nearest look. Replaced with two 3-way controls (none / 1st / 2nd)
per side, each button showing the EVENT'S OWN real hour offset for that snapshot
(`timebins.snapshot_offsets_h`, new). That helper reuses `nearby.py`'s exact anchor convention (b1
IS `step_containing`'s anchor, not one step earlier - the same "anchor counts as the nearest
before" rule the v2 `before` method uses) so the slider's "1st before" and v2's `before` tier mean
the identical overpass for the same event - verified directly (`snapshot_offsets_h`'s b1 matched
`nearby._dt_hours` on the same anchor to 4 decimal places on a real event).
**Trade-off, accepted deliberately:** the precise km² readout (`setReadout`'s `agrid` lookup) only
matches the FIXED 0/6/12/18/24 grid, so it will now rarely show - cell COUNT still does (the
existing graceful fallback, unchanged). The grid DISPLAY itself (`showGrid`/`lbound`/`ubound`) was
always window-agnostic and needed no change; only the precise-area number is affected, and it was
explicitly what the user asked to trade for schedule-accurate snapping.

### D-69 — GOES cloud imagery: centred on the report, every v2 position marked and numbered
Was one crosshair at a single v2 position (or the report as fallback), 30 km radius. Reworked into
a comparison tool: centred on the AS-REPORTED point unconditionally (it has to be the fixed
reference something else is measured against, so it can't also be a v2 result), radius widened to
60 km (D-66's widest match radius is 20 km; 60 km leaves 3x margin so nothing sits at the edge),
and every valid `fire2_r{R}_{M}_lat/lon` gets its own marker - colour = radius (10/15/20 km =
purple/blue/yellow), number = method (1 closest, 2 before, 3 before+after, 4 before2+after2) - so
a reviewer sees at a glance where each (radius, method) combination would have placed the fire
relative to the actual convective cloud. A fixed legend (same corner, every image) spells out the
scheme once rather than requiring it be memorised.
**No mosaicking needed to zoom out**, unlike Sentinel-2's 110 km UTM tiles (D-51) - GOES full disk
is one contiguous raster covering the whole hemisphere, so widening the radius is just a bigger
read window on the same scene, tested well past 60 km with no coverage gap.
**Bugs found and fixed during this build**, both from extending a single-marker design to many:
1. The legend was invisible - drawn on the square canvas BEFORE the circular vignette mask was
   applied, and a square's corners fall outside its inscribed circle, so the mask silently erased
   it. Fixed by drawing the legend after the vignette paste, directly on the final image.
2. Closest/before frequently resolve to the identical position (D-66: 91% of 2023 events) -
   INCLUDING matching the reported point itself for a tight match - so markers routinely stacked
   exactly on top of each other with only the last one drawn ever visible, silently hiding
   agreement instead of showing it. Fixed with a small declutter pass: points within 18 px of each
   other fan into a small ring around their shared position, the reported-point crosshair included
   and always drawn last (on top), so agreement is now visible as a tight rosette rather than one
   marker standing in for twelve.
**Hardening, not a live bug:** `_candidate_scenes`'s daylight check now `.tz_convert("UTC")`
before stripping the tz label, rather than a bare `.tz_localize(None)` that would have silently
kept a non-UTC wall-clock if a catalogue ever returned one - STAC guarantees UTC so this never
actually fired, but the fix was zero-cost so there was no reason to leave the assumption unstated
in the code. Checked `satellite.py` for the same pattern: it subtracts two tz-AWARE timestamps
directly rather than stripping tz first, which is correct regardless of source tz - no change needed there.

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

### D-33 — Land cover is an opaque backdrop; OSM only shows where LC does not
Set the LC pane to `fillOpacity:1`, draw it under every other layer and give it
`pointer-events:none` (no popup, never intercepts a click). Two resolutions in one layer:
1 km over each fire + 5 km, and 10 km over a 120 km collar around **each** fire with the
1 km footprints punched out (`res` = 1000 | 10000). Everything else draws over it unchanged.
**Rejected:** a hull over all of a year's fires for the coarse pass - they cross the
dateline and EPSG:3979 is Canada-only, so the hull becomes a garbage megapolygon and the
10 km read blows up memory. Per-fire collars keep every read a ~240 km box.
**Rejected:** keeping LC semi-transparent over OSM - the user wants LC to read as the base
map, and stacking two basemaps just muddies both.

### D-34 — Per-event land-cover fractions come from the 6 h-pre-onset polygon at native 30 m
The info-panel "top 3 land covers" is sampled over the `(6, 0)` active-burn polygon (the
6 h before onset - the ground that fed the plume), read at the raster's full 30 m, not from
the decimated 1 km map layer. Added `(6, 0)` to `ACTIVE_TABLE_WINDOWS` so that polygon
exists everywhere; `active_fractions` does a windowed native read + `geometry_mask` +
`bincount`. Fractions are over classified cells in the boreal palette only (nodata excluded).
**Rejected:** the ±12 h headline polygon - it includes post-onset spread, which is not the
fuel that drove the pyroCb. **Rejected:** reusing the 1 km layer - a ~10 km² polygon is
only ~10 cells at 1 km, far too coarse for a composition.

### D-40 — Fire activity ships as a gridded burn-time array, not as detections
Snap NFP detections to a 375 m grid (one VIIRS pixel) built in **EPSG:3857**, keep the time
each cell first burned (whole-fire) or unique (cell, observation) pairs (±24 h window), sort
by time, delta-encode. Measured on 2023: 393 k cells at 375 m = **0.64 MB** base64, against
1.7 MB for 276 k points at 1–3 km. Finer *and* smaller, because a cell index delta is a tiny
integer and gzip eats runs of those, while a lat/lon box is six floats.
Two consequences that the rest of the design leans on:
1. **Mercator, not lat/lon.** Web Mercator is the display projection, so a cell index becomes
   a canvas pixel with one multiply-add. Storing lat/lon would force a projection per cell
   per frame.
2. **Time-sorted means every view is a contiguous index range** — cumulative to *t*, a
   sliding window, a static −pre..+post window. That is what makes incremental drawing
   possible (D-41).
**Rejected:** binary typed arrays instead of JSON — measured *larger* after gzip (0.75 MB vs
0.64 MB), because gzip compresses runs of small decimal integers better than raw int32.
**Rejected:** the old proximity-tiered decimation — it inflated each surviving detection to
its dedup cell, so quiet spans drew as 3 km blocks that covered whole small perimeters.

### D-41 — The animation draws to a canvas directly, and only draws what is new
`BurnLayer` subclasses `L.Renderer` (not `L.Layer`) so Leaflet's container positioning and
zoom-animation transforms come for free — the same machinery `L.Canvas` uses — and, like
`L.Canvas`, the context is translated by `-bounds.min` so we draw in layer coordinates.
Advancing the cursor extends the visible index range on the right, so the new slice is
painted *on top of what is already there*; only a rewind, a style change or a map move
forces a full repaint.
**Rejected:** rebuilding a GeoJSON FeatureCollection per frame (what v0.9–v0.13 did) — it
allocated a Leaflet `Path` per box sixty times a second.
**Rejected:** plain `L.Layer` with hand-rolled positioning — got the pan right but needed
Leaflet's `_updateTransform` bookkeeping to survive zoom animation.

### D-42 — Land cover near the fires is a paletted PNG overlay, not polygons
Supersedes D-39. One 8-bit paletted PNG per fire, warped to EPSG:3857 (so `L.imageOverlay`'s
linear placement is exact rather than stretched), ~120 m ground resolution, with a 3×3
majority filter before encoding. 2023: **~1.9 MB for all 69 fires at ~120 m**, against 1.1 MB
for 1 km polygons — ~8× finer for under 2× the bytes, one GPU-composited `<img>` per fire
instead of thousands of canvas paths, and it deletes the `features.shapes` + dissolve +
simplify path that dominated build time and threw topology exceptions.
**Measured on the way:** RGBA PNG is 36 MB (palette compression destroyed); paletted at a
512 px cap is 3.4 MB; the majority filter alone takes ~35 % off. The 10 km polygon backdrop
beyond the fires stays vector — it is a handful of cheap shapes over a very large area.
**Still true:** `active_fractions` reads the raster at native **30 m** for the attribute
table and the info panel. The raster tier is presentation; the science number never changed.

### D-39 — Land cover map layer: finer inside SMALL perimeters, but not true 30 m vector
*(superseded by D-42 — kept for the measurements)*
Three tiers, painter-ordered coarse→fine: 10 km collar, 1 km over perimeter + 5 km, and an
inside-perimeter tier at **≈90 m for fires whose bbox is under `LC_FINE_MAX_KM2` (150 km²)**
only. Larger fires keep the 1 km tier. Sub-km features are simplified again in the payload
(~75 m) and snapped ~22 m. Each feature keeps `res` (m).
**Measured and rejected — true 30 m vector:** polygonising the raster at 30 m over every
pyroCb perimeter, dissolved by class, gave a **14 MB** compressed 2023 land-cover payload
(**57 MB** whole file). 90/300 m over *all* fires was still ~7 MB (31 MB file). Only 90 m
over small fires keeps the whole file at ~17 MB (≈ the v0.12 baseline).
**Not a loss for the science:** `active_fractions` still reads the raster at native 30 m for
the attribute table and the info-panel top-3 — only the cosmetic map fill is limited.
**Deferred to I-3:** a raster/tile render path if the *visual* layer really needs 30 m
everywhere.

### D-43 — Plot the *fire*, not the reported pyroCb point, when position matters
93 % of atlas pyroCb coordinates are exact 0.1° multiples (F13) — a ±5.6 km N–S / ±2.8 km
E–W neighbourhood at 60°N. Rather than pretend otherwise, the map ships both positions:
the reported point (default, with an optional ±0.1° box) and `fire_lat`/`fire_lon`, the
centroid of that event's own 6 h-pre-onset active burn, which VIIRS locates to one 375 m
pixel. `loc_rounded` says which events were snapped.
**Rejected:** silently re-plotting everything on the fire centroid — the reported location
is the published record and some events have no pre-onset burn to centre on.
**Rejected:** interpolating a "better" pyroCb point from the plume geometry — we do not
have the plume, only the atlas row.

### D-44 — Reduce map geometry with `simplify` only; never a morphological close
FEDS perimeters are alpha shapes fitted over 375 m pixels, so the *source* is granular at
375 m. `SIMPLIFY_TOL` 0.003° (~330 m, ~1 pixel — visible when zoomed in) → **0.0003°
(~33 m)**; `COORD_PRECISION` stays 1e-5° (~1.1 m, ~340× finer than a pixel — 1e-6 only
bought invisible decimetres at a real cost in bytes).
Growth bands are **simplify-only at 50 m**, with the ±60 m buffer-close deleted. Douglas–
Peucker error is *bounded by its tolerance* (measured max 92 m at 60 m); a close is not —
bands are thin rings, and `buffer(+r).buffer(-r)` annihilates anything narrower than `2r`,
which measured **4.3 km** of worst-case displacement (F14). Reducing the radius does not
help; any close jumps the max error to kilometres.
**Rejected:** shipping unsimplified perimeters — 12.3 MB for 2023 alone, and the extra
vertices are alpha-shape noise, not real boundary detail.
**Rejected:** a small (±15–25 m) close to suppress boolean crenellation — measured, and it
still cost 2.5–3.1 km of max error for no payload saving.
**Rule going forward:** morphological open/close is valid only on the pixel-lattice
layers, where fusing disjoint 375 m squares is the point.

### D-46 — Report-uncertainty box is per-event, sized to that event's coordinate grid
The atlas nominal precision is 0.1°, but 2013–2016 mix in some 0.01° coordinates. Gating
the box on "exactly 0.1°-rounded" hid it for the finer ones (looked like a bug for 2014).
Now every event ships `loc_prec` (0.1° or 0.01°) and the box is drawn at half that cell,
so it always appears. Whole-number coordinates are treated as 0.1°, not 1° — they are
almost certainly 0.1-grid values that landed on an integer, and a ~110 km box would lie.
The box also obeys the same visibility as the marker (`pcbShown`): region filter, month
cursor, hide-unmatched, layer toggle.

### D-47 — The control panel is dropdowns; one summary tag per section
The panel had grown past a laptop viewport (items clipped off the bottom). Every section
is a `<details>` now, with a one-word state tag in the summary (`Regions … all`,
`Colour by … day of year`) so the collapsed panel still tells you what is set. Years and
Layers stay open; the rest collapse. The legend collapses too.

### D-45 — Animation mode is explicit, never inferred from the year selection
Inferring "years vs months" from `activeYears.size` was the year-sweep bug: sweeping back
to one year flipped the mode, which reset the animation and cleared the year set. The mode
is now a three-way user choice (Off / Years / Months); turning it off restores the year
selection that was in place before the sweep started.
**Rejected:** keeping the inference with guards — every guard added another state to get
wrong, and the user has to be able to say "animate months" while several years are loaded.

### D-40 — Animation burn-age ramp is purple; spread modes are the red→yellow one
Red/yellow on active-fire pixels reads as *intensity*. Burn age (old vs new cell) is not
intensity, so it gets a dark→light **purple** ramp. The red→yellow ramp is reassigned to
the two *spread-rate* colour modes (`avgspread`, `spread95` from the 12-hourly Snapshot
attributes), where "hot = fast" is the intended reading. Deep red = fast, yellow = slow.
**Rejected:** a computed local front-crossing rate per cell — the user pointed out the atlas
already carries `avgspread` / `spread95` per step; using the published number is both less
work and defensible.

### D-41 — Every animation scrubber snaps to the 12 h overpass cadence
NFP `datetime_utc` is real per-overpass time (F6), not 12-hourly, so the raw scrubber let
you slide between overpasses and *see* motion that is pure interpolation. The cursor now
snaps to a 12 h grid (∪ the exact pyroCb onset times). This is the "12 h snapshot" the user
asked for — the honest observation cadence, even though the underlying timestamps are finer.

### D-42 — The global playback bar is contextual, and does not fight the year chips
One bottom-centre transport that changes meaning with the selection: ≥ 2 years → cumulative
year sweep; 1 year → month-by-month reveal (Mar–Oct, ticks labelled with pyroCb counts);
a fire's Animate panel open → the bar hides and that panel's scrubber takes over. It only
*drives* `activeYears` / the date filter once the user presses play or drags it (`gDriving`);
until then it just mirrors state, so manual chip selection is never overridden silently.
**Rejected:** folding the per-fire animation into the global bar — the panel scrubber has
pyroCb ticks, colour modes and loop options bound to one event; duplicating that in a global
widget with three unrelated modes would be worse, not better.

### D-60 — Ignitions get their own row, below the scrubber, not merged with the pyroCb ticks
The whole-fire animation scrubber already carried numbered pyroCb badges above the track
(`drawTicks`). Adding ignitions to that same row risked exactly what it looks like it would:
a merge complex can carry dozens of them (19 on one 2023 fire) interleaved in time with the
pyroCb onsets, and mixing two different kinds of event into one row of marks defeats the purpose
of marking them at all. Ignitions instead get a second row, symmetric but below the track
(`drawIgnTicks`, `#igticks`, amber dots vs the pyroCb row's dark-red numbered squares) - same
greedy row-packing algorithm so a cluster of ignitions still fans out legibly, but visually and
spatially impossible to confuse with "a pyroCb happened here". Left unlabelled on purpose (unlike
the pyroCb badges, which show `pcb_index`): which pyroCb is which matters for cross-referencing
the map and the attribute table; which ignition is which mostly doesn't, and 19 numbered badges
in one row would be unreadable regardless of row-packing. Order and exact date are one hover or
click away. Timing comes from `io.read_ignitions`' `ign_t`/`ign_index`, threaded through
`firesim.fire_series` as `ign_h`/`ign_ix`/`ign_ll` on the same t0-relative clock as `pcb_h`, so
no new per-year payload was needed - `pyrocb/firesim.py` already ships one payload per fire.
**Shipped invisible the first time, for a by-now-familiar reason:** `pipeline.firesim_cached`'s
cache key (`{"cell": C.BURN_CELL_M, "onsets": onsets}`) didn't change when `fire_series`'s output
schema did, so every year's firesim payload kept serving from the pre-`ign_h` cache and the row
silently had nothing to draw - same class of bug as the FRP fields in D-(spread.py v3), same fix:
bump the cache's version tag (`"v": 2`) so a schema change forces a rebuild instead of quietly
being served stale. Verified post-fix: `ign_h` present in a fresh build, real ignitions rendering.

### D-36 — Whole-fire animation is adaptively decimated, tiered on time-to-pyroCb
NOTE: superseded in v0.14 by the uniform 375 m grid (see CHANGELOG v0.14 / D-39 area).
Kept for history. Original text:
Shipping every NFP detection of a big merge complex (10 k+, months long) would bloat the
map. Instead, per detection, take the minimum time gap to any of the fire's pyroCb onsets:
≤ 24 h keeps native footprint + time; ≤ 72 h dedups to a 1 km grid / 6 h bin; beyond that,
3 km / 24 h. If the fire still exceeds `FIRESIM_MAX_PTS` the B/C grids widen until it fits.
So the resolution is highest exactly where the science is (around each event) and coarse in
the quiet spans, which is what the user asked for. Cached as JSON (it's a dict, not a frame).
**Rejected:** 12-hourly Snapshot perimeters — lighter, but they're cumulative alpha-shape
hulls, so you lose the front structure and the merge moments the animation is *for*.
**Rejected:** one fixed decimation — either too heavy near events or too coarse to watch a
run develop.

### D-37 — Blow-up land-cover fractions: first 6 h of burning at/after onset
The `(6, 0)` "6 h before onset" polygon is empty for a blow-up (F4 — no pre-onset
detections). Rather than leave those events without a composition, sample the first
`LC_BLOWUP_SPAN_H` (6 h) of detected burning at or after onset — same construction (union of
localised VIIRS pixel footprints), just shifted in time. The map labels which window was
used. **Rejected:** the ±12 h headline polygon as the fallback — it includes far more
post-onset spread than the 6 h the user specified. **Rejected:** no fallback — silently
drops ~30 % of matched events (the blow-ups) from any land-cover analysis, and they are
not a random subset.

### D-38 — Per-event land-cover fractions come from detected pixels only, at 30 m
Locked in after an explicit check: `active_fractions` masks the raster to the *union of the
VIIRS pixel footprints* in the window (localised to the pyroCb's cluster by
`active._localise`), read at the raster's native 30 m. Not the fire perimeter, not a convex
hull, not the decimated 1 km map layer. `lc_n_px` reports the 30 m cell count so a thin
sample is visible.

### D-35 — Drop the static "active burn at onset" map layer
The ±12 h active-burn outline layer is redundant once a pyroCb click opens the window
slider + animation (D-30, D-31), which shows the same geometry interactively for any
window. The polygons remain in `pyrocb_active_areas.gpkg` and as attribute-table columns.
**Rejected:** keeping it as a default-off layer - one more toggle for a strictly worse view.

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
