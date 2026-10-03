# FINDINGS

Things we learned from the data itself — schema gotchas, failure modes, and the checks
that caught them. Every entry is backed by a measurement against the real files, not by
reading documentation. **Add to this file whenever the data surprises you.**

Format: what we believed → what is actually true → how we know → what we did about it.

---

## F1. `mergid` is the only valid join key

**Believed:** `fireid` identifies a fire across products.

**Actually:** `fireid` is the ID *before* merging and is unsafe in two of the four products.
The stable key is `mergid`:

```
ignitions.mergid == final.fireid == final.mergid == NFP.mergid == snapshot.mergid
```

Perfect 1:1 — 10,707 fires in 2023, intersection 10,707, both differences 0.

* **Snapshot `fireid` is incomplete.** Only ~95% of snapshot rows carry a `fireid` present in
  the final file; the rest are components absorbed by a merge. Matching final fire 12406 by
  `fireid` at `20230715PM` returns a 13.3 km² fragment **92 km away**, while `mergid` returns
  5 rows including the real 450.9 km² main body.
* **NFP `fireid` is not a fire id at all** — it equals `fid` and the row index, `0..N-1`, in all
  606 files of 2023. Using it would attribute pixels essentially at random.

**Consequence:** a real bug shipped in the first map — ignitions were filtered on `fireid`, which
is unique per row, so the "only fires with a pyroCb" toggle kept only each merged fire's *primary*
ignition and dropped the other 1,444 (2023). After the fix the count went 68 → 518.

---

## F2. `farea` is cumulative burned area, not actively-burning area

**Believed (natural but wrong):** a perimeter shrinks as parts of the fire burn out.

**Actually:** it never shrinks. Checked across all 608 snapshots × 10,707 fires (2023):

| check | result |
|---|---|
| max-snapshot `n_pixels` ÷ final `n_pixels` | **exactly 1.0 for every fire** (std 0.0) |
| `farea` within 0.01% of final | 10,698 / 10,707 (99.92%) |
| snapshot area exceeding final | 9 fires (6 by >1%) |

The attribute table agrees: `n_pixels` is *"Number of **total** fire pixels"* while `n_newpixels`
is *"at current time step"*, and the burning front is carried separately by **`flinelen`**
(*"Length of **active** fire line"*), `isactive` and `t_inactive`. Burnt-out area stays inside the
perimeter — which is why the pixel count can never decrease.

The 9 exceptions are alpha-shape fitting artefacts on *tiny* fires: identical pixel counts, only
the polygon differs. `mergid 16049` has 2 pixels, snapshot 0.141 km² (exactly one VIIRS pixel)
vs final 0.0027 km² — a collapsed sliver. pyroCb fires are the largest in the dataset, so this
cannot reach our analysis population.

---

## F3. AM/PM are local-solar overpass groups, and there is a 10-hour nightly blind window

**Believed:** AM/PM split the calendar day at 12:00 UTC.

**Actually:** they are the two VIIRS overpass groups in *local solar* time. Measured against real
NFP `datetime_utc` (`local = utc + lon/15`):

* AM pixels fall at local **0.5–7.0 h** (median ~2.0)
* PM pixels fall at local **7.1–14.4 h** (median ~12.5)
* **No observation at all between local 14.4 h and 00.5 h**

A single AM/PM file therefore spans ~24 h of UTC, not 12. And **pyroCb onset peaks at local
15.0 h** (n=484: p5 11.8, p25 13.6, p50 15.0, p75 16.6) — squarely inside the blind window.

**Consequences, both load-bearing:**

* **Straddle.** 19.3% of 2023 events have ≥1 same-day-PM pixel detected *after* onset; 10.9% have
  the *majority* after. Fixed by choosing the pre-band step from the NFP **median** observation
  time per (mergid, step) rather than from a clock — one stray late pixel then cannot disqualify
  an overpass whose bulk precedes onset.
* **`frac_burned_before` is biased low and the bias is correlated with time of day**, worst for
  early-afternoon pyroCb. Not random noise. `obs_lag_h` records the staleness; use it as a
  covariate or restrict to `obs_lag_h < 6 h`. Measured for 2023: median 3.0 h, 94/119 under 6 h.

---

## F4. A fire can blow up *before* it is ever observed — matching must look forward

**Believed:** the fire that produced a pyroCb is present in the snapshot at or before onset.

**Actually:** not if it ignited and blew up inside the nightly blind window (F3). Its first look
is the **next** overpass.

**Type case — `2023-R90877` (`PCB_2023_0625_0000`):**

| | |
|---|---|
| pyroCb onset | 2023-06-25 00:00 UTC = local **15.5 h** (inside the blind window) |
| step containing onset | `20230624PM` — fire absent; nearest active fire **86 km** away → unmatched |
| fire `mergid 11918` first detected | 2023-06-25 **10:17 UTC**, i.e. **10.3 h after onset** |
| pixels on that first look | **117** — so it was plainly burning at onset, just unseen |
| final size | 135.7 km², 2023-06-25 → 07-15 |

The same fire's *second* pyroCb (06-28) matched normally at 5.5 km, which is what made the miss
obvious: two events with the same `fire_name` disagreeing.

**Fix:** `matching._rescue_new_fire` searches up to `MATCH_FORWARD_STEPS` (2) overpasses ahead for
events still unmatched, tagging them `inside_newfire` / `nearest_newfire` with `match_fwd_steps`.
Recovered 2 of 40 unmatched events in 2023 at 5.24 km.

**Deliberately scoped to matching only.** The pre-band still uses the last observation strictly
*before* onset, so an ignition-day blow-up correctly gets `band0_undefined` and `area_pre = 0` —
it really had burned nothing observable beforehand.

---

## F5. Snapshots contain only *active* fires, and are not perfectly nested

* `isactive == 1` for every row of every snapshot checked; a fire that goes quiet **vanishes from
  the file** rather than appearing with a smaller area. So walking backwards to a fire's last
  appearance is *correct*, not an approximation — its cumulative footprint is unreported, not
  shrunk. (`bands._fetch` does this; 0 back-walks were needed in 2023.)
* Total snapshot area is far below the season total — 97,459 km² at `20230715PM` vs 233,171 km²
  for the 2023 season, because burnt-out fires have dropped out.
* Between consecutive steps the earlier polygon loses on average **0.17%** of its area (max 50%
  in a sampled outlier) to re-polygonisation. Small, but enough that plain differencing produces
  negative-area bands. Hence the cumulative-union construction.
* A single snapshot file holds **70 duplicate-`mergid` rows** (the not-yet-executed components of
  a merge), so `dissolve(by="mergid")` is mandatory. Validation that it is right: dissolved
  `farea` 97,459 vs summed polygon area 97,457 km² — 0.002% agreement.
* `mergedate` values point to the *future* relative to the snapshot: `mergid` is a retroactively
  back-filled forward pointer, because the season was reprocessed after the fact.

---

## F6. NFP is incremental, exact, and reproduces the final product

* Consecutive periods share **zero** rounded coordinates (checked on 10 sampled pairs) — the
  files are per-step increments, not cumulative.
* Summing NFP rows per `mergid` over a whole year reproduces `final.n_pixels` **exactly for
  100% of the 10,707 fires**. This is the strongest end-to-end check available on the entire
  join chain, and it is asserted in `variables.build_events`.
* The remaining 4,938 NFP rows belong to `mergid`s absent from the final product, each with
  exactly 1 pixel — single-detection fires filtered out downstream.
* One VIIRS pixel = **0.141 km²** exactly (confirmed as the `farea` of a 1-pixel ignition row).
* `datetime_utc` is a real timestamp, so NFP-derived variables are the only ones that can
  separate two pyroCb falling in the same half-day.

---

## F7. Type and unit traps

| Field | Trap |
|---|---|
| `avgspread`, `spread95` (snapshot) | stored as **strings** (`'0'`, `'166.83…'`) — silently poison any aggregation |
| `mergedate` (snapshot) | string like `'2023-9-23-PM'`, not a datetime |
| `fperim` (snapshot **and** ignitions) | **mixed km/m regimes** — values ~1–10 alongside ~1000+. Do not use without checking |
| `fire_name` (PCBAtlas) | mixes ints (Quebec fires are numbered: `172`, `260`) with strings (`EWF-031-2023`) — breaks parquet round-trips unless cast |
| NFP column order | `y` comes **before** `x` |
| CRS | snapshots / ignitions / NFP are **EPSG:3571**; final perimeters are **EPSG:4326** |
| EPSG:3571 | is *"WGS 84 / North Pole LAEA Bering Sea"*, **not** EASE-Grid 2.0 North (that is 6931) |
| Snapshot glob | must exclude `*_FL.gpkg` (fire lines), or you read each step twice |
| NFP step gaps | are real — 2012 is missing `20120127PM` and `20120218PM` |

---

## F8. `pcb_id` collides on outbreak days

`PCB_YYYY_MMDD_HHMM` gives only **449 distinct values for 484 events**. 26 timestamps carry
several events, because a pyroCb outbreak blows up several *separate* fires in one afternoon:
`PCB_2023_0729_2220` is 2023WB7, 2023WB18 **and** MWF-043-2023, hundreds of km apart.
Resolved with a deterministic `a`/`b`/`c` suffix ordered by (lat, lon).

---

## F9. Merged fires strain the concentric-band model

Bands assume one spreading front. But `mergid` groups separately-ignited fires that later grew
together, and the biggest fires — exactly the pyroCb producers — are the most likely to be
complexes:

* 1,444 of 12,151 (11.9%) 2023 ignition rows are non-primary members of a merge complex.
* `mergid 3283` (West Siberia) absorbs **115** separate ignitions between 2023-04-09 and 05-23.
* `mergid 7783`, the largest 2023 fire (15,054 km², 12 pyroCb), merges **three separately named
  fires** — 2023-G90288, SS008-23 and FS001-23 — with **40 ignitions**.

For those, rings are several fronts joining, not one front advancing, and `band_rate_km2_h` is
not a spread rate. Always report `n_ign_total` and `ign_merge_span_days` beside band statistics.

---

## F10. Unmatched pyroCb are non-random

38 of 159 (2023) have no fire within 20 km. They are **not** a dormancy artefact — retrying
against the final perimeters with a timing constraint rescued **zero** of them, so the fires are
genuinely absent from the atlas. They cluster in two places:

* the **western Canadian Cordillera** (BC interior — `2023-R/K/C…` ids, Crater Creek at 49.2°N),
  where steep terrain and persistent cloud defeat VIIRS detection;
* **southern Siberia / the Mongolian steppe** (~47–52°N), where the fuel is grassland the boreal
  product does not map.

Dropping them biases the sample toward well-detected boreal forest fires. Report unmatched counts
**with a regional breakdown**, never as a bare total.


---

## F11. The 2021 NFP files are from a different processing run

Caught by the standing `n_pixels` assertion (F6), which failed for 2021 and no other year.
Two independent defects, one root cause — the 2021 NFP was not re-issued after the season was
reprocessed:

**a) Coordinates are lon/lat, not EPSG:3571 metres.** 2021 `x` spans −179…180 and `y` 50.6…71.4;
every other year spans ±4.5 million metres. The attribute table documents EPSG:3571 for all
years, so this is only detectable from the magnitudes. `io._read_nfp_raw` now sniffs it
(`|x| ≤ 180 and |y| ≤ 90`) and reprojects, setting `nfp_reprojected`.

**b) `mergid` is from a stale id space.** 2021 NFP carries **18,931** distinct ids against the
final product's 11,688, and `final.mergid` is *not* a subset of them — yet the pixel total agrees
to 0.04% (1,652,778 vs 1,653,384). The detections are all present and correctly located; only
the labels are wrong.

**Repair:** since the perimeters are built from these pixels, re-attributing by point-in-polygon
against the final perimeters recovers most of the mapping — 7.0% → **68.3%** of fires exact,
**83.1% within ±2 pixels**, median difference 0. Not perfect: pixels on perimeter edges and in
overlapping polygons cannot be resolved this way.

**Treatment:** 2021 is flagged, not faked. `qc_nfp_repaired` marks every 2021 row and
`qc_nfp_pixel_exact` marks the individual rows that still reconcile. The hard assertion still
runs for the ten sound years. **Treat 2021 NFP-derived columns — the `n_pixels_*`, `frp_*` and
`active_*` families — as approximate.**

*Lesson worth keeping: the assertion earned its place. Without it, 2021 would have produced
plausible-looking pixel counts and active-burn polygons that were silently attributed to the
wrong fires.*

---

## F12. "Active burn area" was pulling in fires that only merged later

**Believed:** filtering `NFP[mergid == fire] ` to ±12 h of onset gives the ground burning
around the pyroCb.

**Actually:** for a big fire it gave *several* fires. `NFP.mergid` is the season-**final**
id (F5), so a fire's row set already contains every neighbour that merged into it weeks
later. A narrow *time* filter does nothing about that — the neighbours were also burning
that afternoon. Measured on the six largest 2023 pyroCb fires, ±12 h of onset:

| event | detections | median dist to pyroCb | cloud diagonal |
|---|---|---|---|
| PCB_2023_0705_1940 | 470 | **121 km** | 211 km |
| PCB_2023_0922_2020 | 4708 | 39 km | 198 km |
| PCB_2023_0716_2220 | 3107 | 37 km | 128 km |

`PCB_2023_0922_2020` is the clean case: 4 snapshot rows with 4 distinct `fireid`, `mergedate`
**2023-09-23** — the day *after* the pyroCb. Four fires 48–123 km apart, dissolved into one
because they eventually touched.

There is also a second mechanism: `PCB_2023_0705_1940` is a *single* `fireid` the atlas tracks
as one 9,448 km², 35-part blob spanning 200 km. Not a merge — just a sprawling complex.

**Fix (`active._localise`).** After the time filter, cluster the detections (KD-tree
connected components at `ACTIVE_BRIDGE_M` + a pixel diagonal) and keep only the cluster the
pyroCb sits in — or the nearest within `ACTIVE_ATTACH_CAP_KM` if the 0.1°-rounded location
lands in a gap. `n_clusters` and `detach_km` record what happened.

Result across the matched set: **~48 % of in-window pixels dropped**; headline active area
median 92 → **31 km²**, max 15,000+ → **1,280 km²**. `detach_km > 0` for 65 % of events
(rounding puts the point just off a pixel; snap is usually < 1 km), `> 5 km` for only 27 —
those are the genuinely weak localisations and are flagged.

---

## F13. The atlas pyroCb coordinate is a neighbourhood, not a point — and the fire is not

Measured over all 484 events in `PCBAtlas.xlsx`:

| | share |
|---|---|
| `lat` an exact multiple of 0.1° | **92.6 %** |
| `lon` an exact multiple of 0.1° | **93.4 %** |
| both rounded (`loc_rounded`) | ~93 % |

Only ~7 % of events carry a second decimal. So for most events the plotted pyroCb dot
carries **±0.05° of quantisation — about 5.6 km N–S, and 2.8 km E–W at 60°N** (half-cell),
i.e. a worst-case displacement near 11 km / 5 km. No processing can recover that; it is
how the source was recorded.

**The matched fire is three orders of magnitude better.** VIIRS active-fire detections
are located to one 375 m pixel, so the centroid of an event's own 6 h-pre-onset active
burn (`active_area_km2_pre6_post0`) is a genuinely accurate position for *the fire that
made the plume*. The map now ships both (`fire_lat` / `fire_lon` beside the reported
point) and can plot either, with the ±0.1° report box drawable on demand.

**Consequence for analysis:** never compute a distance, a land-cover sample or a
terrain covariate from the reported pyroCb lat/lon. Use the matched fire geometry.
`match_dist_km` is measured against the rounded point and so inherits the same floor —
which is exactly why anything above ~5 km was already flagged provisional.

### Map-side geometry fidelity

Separately, the *map* was simplifying perimeters at 0.003° (~330 m). Measured deviation
from the unsimplified geometry over 2023's 69 pyroCb fires (Hausdorff, metres):

| perimeter simplify | max dev | median dev | 2023 payload |
|---|---|---|---|
| 0.003° (~334 m) — shipped through v0.15 | **585 m** | 335 m | 0.55 MB |
| 0.001° (~111 m) | 199 m | 124 m | 0.84 MB |
| **0.0003° (~33 m) — v0.16** | **56 m** | 33 m | 1.55 MB |
| 0.0001° (~11 m) | 14 m | 11 m | 2.50 MB |

The old tolerance was displacing boundaries by up to **585 m — more than one full VIIRS
pixel**, which is why it was visible when zoomed in. FEDS perimeters are alpha shapes
fitted over 375 m pixels, so the source itself is granular at 375 m: ~33 m is an order of
magnitude below that (bounded worst case 56 m, ~1/7 of a pixel) and anything finer just
stores alpha-shape noise at 2–8× the bytes. The coordinate snap stays at 1e-5° (~1.1 m),
already ~340× finer than a pixel; 1e-6° bought only invisible decimetres.

Growth bands were fixed separately and far more dramatically — see **F14**.

---

## F14. Buffer open/close destroys thin geometry — never use it to "tidy" growth bands

Growth bands are **thin rings** (the area burned between two successive cumulative
unions). A `buffer(+r).buffer(-r)` close annihilates any feature narrower than `2r` and
merges anything separated by less than `2r`, so on rings it does not smooth — it deletes.

Measured on 2023's 183 bands, deviation from the unsmoothed geometry (Hausdorff, metres):

| band map treatment | payload | max dev | median dev |
|---|---|---|---|
| close ±60 m + simplify 150 m *(shipped through v0.15)* | 1.53 MB | **4304 m** | 150 m |
| close ±25 m + simplify 25 m | 3.39 MB | **2487 m** | 25 m |
| close ±15 m + simplify 20 m | 3.60 MB | **3051 m** | 20 m |
| **simplify 25 m, no close** | 3.36 MB | **42 m** | 25 m |
| **simplify 60 m, no close** | 2.24 MB | **92 m** | 60 m |

Note the pattern: adding *any* close jumps the maximum error to kilometres while barely
changing the payload, and shrinking the close radius does not help. Douglas–Peucker
simplification, by contrast, is **bounded by its tolerance** — 60 m of tolerance gives
92 m of worst-case error.

**So the band layer is simplify-only at `BAND_SIMPLIFY_M`, and the close is gone.** The
±60 m close had been inherited from the old active-burn *pixel lattice* layer, where
closing a grid of disjoint 375 m squares into a blob was the entire point. It was never
appropriate for bands, and it had been silently displacing their boundaries by up to
4.3 km — a far larger error than the 330 m perimeter simplification we set out to fix.

**General rule for this project:** morphological open/close is only ever valid on the
pixel-lattice layers. Everywhere else, reduce vertices with `simplify`, whose error is
bounded and reportable.


---

## F15. Land Cover of Canada maps the burn scar, not the fuel — Copernicus C3S does not

`lc_compare/compare_landcover.py`: 20,000 VIIRS fire pixels in Canada (pixels.parquet),
Canada 30 m majority class over a ~330 m window vs the Copernicus C3S 300 m class of the fire
year, both mapped onto 9 shared fuel classes (`lc_compare/results.json`,
`lc_compare/confusion_canada_vs_c3s.csv`).

| | agreement |
|---|---|
| overall | **47%** |
| fire **after** the Canada map's epoch (2017 vs 2015 map; 2021–23 vs 2020 map) | **65%** |
| fire **in or before** the epoch (2013–16 vs 2015; 2018–20 vs 2020) | **35%** |
| C3S fire year vs C3S year before the fire | **99.5%** |

Agreement collapses exactly when the Canada map was made *after* the fire: it then shows the
burn scar. Canada calls **31%** of fire pixels grassland (C3S:
0.8%) and 43% needleleaf (C3S:
75%); where Canada does say needleleaf, C3S agrees
91% of the time. C3S is annual and barely
moves across the burn (99.5 % same-year vs year-before), so it describes the fuel that burned.
It is also global (Alaska, Russia), where the Canada product has nothing. Caveats: 300 m vs
30 m, and C3S stops at 2022 (2023 fires use the 2022 map). Decision: D-87.
