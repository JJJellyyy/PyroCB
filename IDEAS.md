# IDEAS

Parked ideas and bigger directions that are **not** on the current work list. Nothing here
is decided — each entry is a sketch plus a rough feasibility read. Do not start any of these
without an explicit go-ahead.

Release history: [CHANGELOG.md](CHANGELOG.md) · rationale: [DECISIONS.md](DECISIONS.md) ·
data discoveries: [FINDINGS.md](FINDINGS.md).

---

## I-1 — Animate a whole fire event (first pyroCb → last, or first ignition → last burning pixel)

**Ask (2026-09-06):** selecting a pyroCb, play the entire event — not just the ±24 h active
window — from the first ignition of the complex to the last detected burning pixel, looping in
a chosen number of seconds.

**Feasibility: doable in the HTML, but not from the per-detection pixel arrays we ship now.**
- The current animation reads `byPcbPix` — NFP detections in the widest **±24 h** window,
  localised to the pyroCb's own cluster. A whole-event animation needs the fire's *entire*
  detection history (weeks to months), which for the biggest merge complexes is 10 k+ pixels
  and would bloat the payload if shipped raw per event.
- **Better source: the 12-hourly Snapshot perimeters.** One simplified cumulative polygon per
  half-day step, first detection → season end: ~60–120 frames for even a long-lived fire,
  a few KB each after simplification. Cumulative growth animation, pyroCb onsets marked on the
  timeline. Far lighter than pixels and it already exists on disk.
- New work: a `snapshot_series(mergid)` pipeline stage (cache per year), a per-fire payload
  keyed by `mergid`, and a second animation mode in the info panel ("This window" vs "Whole
  fire"). Loop-length control reused. Estimate: similar size to the v0.9 slider/animation work.

**Status:** waiting on a go-ahead. If yes, do the snapshot-series version, not raw pixels.

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
