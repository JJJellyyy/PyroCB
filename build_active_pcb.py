"""'Active pyroCb' pixels for the map: DATA[year]["apc"], straight from the dataset's labels.

The tabular dataset (pcbcube, LABEL_MODE = "active_buffer") labels as positive every VIIRS pixel - of
ANY fire - within 20 km of the pyroCb's fixed (matched) coordinate that was detected in the overpass
just before onset or just after it. This writes exactly those pixels per event, so the map shows the
dataset's positives, not a re-derivation:

    {pcb_id: {"c": [lon, lat],            fixed coordinate (pcb_events lon_c/lat_c)
              "r": 20,                    buffer radius, km
              "mg": credited mergid, "wk": weak_match, "src": pos_source,
              "s": [step_before, step_after],
              "fires": [mergid, ...],     fires with active pixels (credited one first)
              "p": [[dlon, dlat, hx, hy, side, fire_ix], ...]}}
                                          dlon/dlat: offset from "c" in 1e-5 deg; hx/hy: footprint half
                                          size in 1e-5 deg (from DS/DT); side 0 = before onset, 1 = after

    python build_active_pcb.py [year ...]      then   python map_tool.py join
"""
from __future__ import annotations

import base64
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DC = HERE.parent / "OUTPUT_DATA" / "datacube"
DATA = HERE / "map_src" / "data.json"
COLS = ["lon", "lat", "DS", "DT", "fire_key", "mergid", "si", "step", "label", "pcb_id", "pcb_step_offset"]


def pack(o):
    return base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def year_payload(y: int, ev: pd.DataFrame) -> dict:
    pix = pd.read_parquet(DC / "_cache" / f"pixels_{y}.parquet", columns=COLS, filters=[("label", "=", 1)])
    if pix.empty:
        return {}
    pix = pix.assign(pcb_id=pix["pcb_id"].astype(str).str.split(";")).explode("pcb_id")
    out = {}
    for pid, g in pix.groupby("pcb_id"):
        if pid not in ev.index:
            continue
        e = ev.loc[pid]
        lon0, lat0 = float(e["lon_c"]), float(e["lat_c"])
        cred = int(e["mergid"]) if pd.notna(e["mergid"]) else None
        fires = sorted(set(g["mergid"].astype(int)), key=lambda m: (m != cred, m))
        fix = {m: i for i, m in enumerate(fires)}
        coslat = np.cos(np.radians(g["lat"].to_numpy(float)))
        hy = g["DT"].to_numpy(float) / 2 / 111.32 * 1e5
        hx = g["DS"].to_numpy(float) / 2 / (111.32 * coslat) * 1e5
        p = np.column_stack([np.rint((g["lon"].to_numpy(float) - lon0) * 1e5), np.rint((g["lat"].to_numpy(float) - lat0) * 1e5),
                             np.rint(hx), np.rint(hy), g["pcb_step_offset"].to_numpy(int),
                             [fix[int(m)] for m in g["mergid"]]]).astype(int)
        steps = g.sort_values("si").drop_duplicates("si")
        before = steps.loc[steps["pcb_step_offset"] == 0, "step"]
        after = steps.loc[steps["pcb_step_offset"] == 1, "step"]
        out[pid] = {"c": [round(lon0, 5), round(lat0, 5)], "r": 20, "mg": cred,
                    "wk": bool(e["weak_match"]), "src": str(e["pos_source"]),
                    "s": [before.iloc[0] if len(before) else None, after.iloc[0] if len(after) else None],
                    "fires": fires, "p": p.tolist()}
    return out


def main(years=None):
    ev = pd.read_parquet(DC / "pcb_events.parquet").set_index("pcb_id")
    data = json.loads(DATA.read_text(encoding="utf-8"))
    years = years or sorted(int(y) for y in data)
    for y in years:
        apc = year_payload(y, ev[ev["onset_utc"].dt.year == y])
        data[str(y)]["apc"] = pack(apc)
        n = sum(len(v["p"]) for v in apc.values())
        print(f"  {y}: {len(apc)} events, {n:,} active-pyroCb pixels, {len(data[str(y)]['apc']) / 1e6:.2f} MB", flush=True)
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
    print(f"data.json now {DATA.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
