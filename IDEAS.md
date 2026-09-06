# IDEAS

Parked ideas and bigger directions that are **not** on the current work list. Nothing here
is decided — each entry is a sketch plus a rough feasibility read. Do not start any of these
without an explicit go-ahead.

Release history: [CHANGELOG.md](CHANGELOG.md) · rationale: [DECISIONS.md](DECISIONS.md) ·
data discoveries: [FINDINGS.md](FINDINGS.md).

---

## I-1 — Animate a whole fire event — DONE (v0.12)

Shipped as the *Whole fire* scope in the info-panel Animate tab. Went with **adaptively
decimated NFP detections** (`pyrocb/firesim.py`), not the snapshot-perimeter idea sketched
here — the pixel footprints keep the front structure and the merge moments, which is the
point. Native res within ±24 h of a pyroCb, 1 km/6 h to ±72 h, 3 km/24 h between. See D-36.

---

## I-3 — Country-wide land cover with a resolution gradient

**Ask (2026-09-06):** display land cover for the *entire country* at a very rough scale and as
fine as possible around pyroCb fires — a smooth resolution gradient rather than the current
two fixed tiers (1 km over fires + 5 km, 10 km per-fire collar). "Probably many underlying
optimisations."

**Sketch / options:**
- **Add a third, national tier.** One decimated read of the whole raster at ~25–50 km
  (`3800×3200` cells at 50 km — a few seconds, tens of MB RAM), dissolved + hard-simplified,
  as an always-on base. Then the existing 10 km collar, then 1 km over fires. Cheap, and
  the payload stays small because coarse polygons have few vertices.
- **True gradient:** cell size as a smooth function of distance to the nearest pyroCb
  (e.g. `cell = clip(0.3·d_km, 1 km, 40 km)`), one variable-resolution polygonisation.
  Cleaner visually, more code; needs a distance-banded read loop.
- **Raster tiles instead of polygons:** pre-render the coarse national layer to a handful of
  PNG `L.imageOverlay`s (or a tiny tile pyramid embedded as data URIs). Crisper at speed,
  but a different rendering path and base64 PNGs are heavy.
- **Optimisations to look at:** read the raster's internal overviews (`rasterio` `overview_level`)
  instead of decimating band 1; `rio warp`/`gdal_translate` a one-off coarse GeoTIFF checked
  into the repo so the build doesn't touch the 2 GB file; topology-preserving simplification
  (`shapely.simplify` vs `topojson`) to cut vertices harder; quantise coordinates more
  aggressively for the coarse tier (100 m snap is invisible at 25 km).

**Recommendation:** start with the third fixed national tier (small, low-risk), measure the
payload, then decide whether the smooth gradient is worth it.

**Progress (v0.13):** the *inward* direction was attempted and measured. True 30 m vector
land cover over every pyroCb perimeter = **14 MB** compressed for 2023 alone (57 MB whole
file); 90/300 m over all fires still ~7 MB. Shipped a compromise: **90 m inside small
perimeters only** (bbox < 150 km²), whole file back to ~17 MB (D-39). The exact 30 m
composition is still in the panel/table via `active_fractions`.

**Done (v0.14) — the raster path.** Pre-rendering each fire's land cover to a paletted PNG
`L.imageOverlay` was the answer, exactly as sketched here: **~120 m over every fire + 5 km
for ~1.9 MB** (2023), ~8× finer than the 1 km polygons it replaced, GPU-composited, and much
faster to build. See D-42. Measured limits: RGBA kills palette compression (36 MB); a 512 px
cap is 3.4 MB; a 3×3 majority filter is worth ~35 %.

**Still open:**
- **The national rough tier** — one ~25–50 km pass over the whole raster as an always-on
  base, so land cover is visible between fires rather than only around them. Now cheap to
  add as one more PNG per region.
- **A smooth `cell = f(distance-to-pyroCb)` gradient** instead of the current fixed tiers
  (~120 m PNG near fires, 10 km polygons beyond).
- Pushing the per-fire PNG below 120 m: the cap is `LC_PNG_MIN_M` / `LC_PNG_CAP_PX`; 60 m
  measured at ~3.8 MB for 2023 if that trade is ever worth it.

**Status:** raster path shipped; national tier and gradient still open.

---

## I-2 — Standalone VIIRS FireAtlas visualiser (not an HTML)

**Ask (2026-09-06):** since we have the full VIIRS FireAtlas data, we could write our own
program to visualise fire evolution rather than fighting the browser.

**Sketch.** A small local app — Python + `pydeck`/`datashader`/`matplotlib`-animation, or a
`deck.gl` desktop shell — that reads the Snapshot/NFP products directly and renders:
- time-scrubbable fire growth for any `mergid` or region, at full pixel resolution;
- FRP as a colour/height channel;
- pyroCb onsets and injection heights on the timeline;
- export to MP4/GIF.

**Why park it.** The HTML is the shared deliverable and is close to "great". A native tool is
a parallel track with its own packaging/its own audience question. Revisit once the HTML is
locked.

**Status:** explicitly do not start without asking.
