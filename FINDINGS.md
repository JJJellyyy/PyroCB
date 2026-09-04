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
