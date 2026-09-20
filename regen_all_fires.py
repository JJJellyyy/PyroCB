"""One-off: rebuild per/ign/sim for every year with MAP_ONLY_PCB_FIRES=False and
ANIM_MIN_PIXELS-gated animation eligibility (config.py), then splice the three
payloads into the existing fire_pcb_map.html DATA object in place.

Deliberately does NOT go through make_notebook.py - that generator writes a Jupyter
notebook, not the HTML file itself, and fire_pcb_map.html has since diverged from it
with a round of hand-written UI work this script must not disturb. pipeline.run_year()
gives the same `ev` (events) DataFrame make_notebook.py's own per-year loop builds, up
to the land-cover/v3 columns that per/ign/sim never touch - so this reproduces exactly
what a full regen would produce for these three keys, without running the rest of it.
"""
from __future__ import annotations

import base64
import gzip
import json
import time

import pandas as pd
import shapely

from pyrocb import config as C
from pyrocb import io
from pyrocb import pipeline

YEARS = list(range(2013, 2024))


def doy_of(df, pre="tst"):
    d = pd.to_datetime(dict(year=df[f"{pre}_year"], month=df[f"{pre}_month"],
                             day=df[f"{pre}_day"]), errors="coerce")
    return d.dt.dayofyear


def pack(gj):
    raw = json.dumps(gj, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode(), len(raw)


def build_year(y):
    t0 = time.time()
    _, _, ev, _ = pipeline.run_year(y, verbose=False)

    per = io.read_final(y, metric=False)
    per["geometry"] = per.geometry.simplify(C.SIMPLIFY_TOL, preserve_topology=True)
    per["geometry"] = shapely.set_precision(per.geometry.values, C.COORD_PRECISION)
    per = per[~per.geometry.is_empty & per.geometry.notna()].copy()
    per["doy"] = doy_of(per)
    per_out = per[["mergid", "farea", "duration", "n_pixels", "lcc_final", "doy", "geometry"]] \
        .rename(columns={"farea": "farea_km2", "duration": "duration_days"}).round({"farea_km2": 3})

    ign = io.read_ignitions(y).to_crs(C.WGS84)
    ign["doy"] = doy_of(ign)
    ign["n_ign_total"] = ign.groupby("mergid")["mergid"].transform("size")
    ign["geometry"] = ign.geometry.representative_point()
    ign_out = ign[["mergid", "ign_index", "n_ign_total", "ftype", "n_pixels",
                   "farea", "doy", "geometry"]].rename(columns={"farea": "farea_km2"})

    linked = set(ev.loc[ev["match_method"] != "unmatched", "mergid"].dropna().astype(int)) \
        if len(ev) else set()
    per_out["pcb_linked"] = per_out["mergid"].isin(linked)
    ign_out["pcb_linked"] = ign_out["mergid"].isin(linked)

    rmap = (ev.dropna(subset=["mergid"]).assign(mergid=lambda d: d["mergid"].astype(int))
              .groupby("mergid")["region"].first().to_dict()) if len(ev) else {}
    per_out["region"] = per_out["mergid"].map(rmap)
    ign_out["region"] = ign_out["mergid"].map(rmap)

    sim = pipeline.firesim_cached(y, ev, linked) if len(ev) else {}

    per_b64, per_raw = pack(json.loads(per_out.to_json(drop_id=True)))
    ign_b64, ign_raw = pack(json.loads(ign_out.to_json(drop_id=True)))
    sim_b64, sim_raw = pack(sim)

    n_anim_only = sum(1 for v in sim.values() if not v.get("pcb_h"))
    print(f"  {y}: {len(per_out):>6} perims ({per_out.pcb_linked.sum():>3} pcb-linked)  "
          f"{len(ign_out):>6} ign  |  sim {len(sim):>4} fires "
          f"({n_anim_only:>4} anim-only, no pyroCb)  "
          f"per {per_raw/1e6:5.2f}MB ign {ign_raw/1e6:5.2f}MB sim {sim_raw/1e6:5.2f}MB  "
          f"[{time.time()-t0:5.1f}s]", flush=True)
    return per_b64, ign_b64, sim_b64


def main():
    new_payloads = {}
    for y in YEARS:
        new_payloads[str(y)] = build_year(y)

    path = "fire_pcb_map.html"
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    i = text.index("const DATA=")
    data, end = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])

    for y, (per_b64, ign_b64, sim_b64) in new_payloads.items():
        data[y]["per"] = per_b64
        data[y]["ign"] = ign_b64
        data[y]["sim"] = sim_b64

    new_data_json = json.dumps(data, separators=(",", ":"))
    new_text = text[:i] + "const DATA=" + new_data_json + text[i + len("const DATA=") + end:]

    import shutil
    shutil.copy2(path, path + ".bak2")
    with open(path, "w", encoding="utf-8") as f:
        f.write(new_text)
    print(f"\nwrote {path} (backup at {path}.bak2)")


if __name__ == "__main__":
    main()
