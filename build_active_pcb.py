"""'Active pyroCb' pixels for the map: DATA[year]["apc"], straight from the dataset's labels.

The tabular dataset (pcbcube, LABEL_MODE = "active_buffer") labels as positive every VIIRS pixel - of
ANY fire - within 20 km of the pyroCb's fixed (matched) coordinate that was detected in the overpass
just before onset or just after it. This writes exactly those pixels per event, so the map shows the
dataset's positives, not a re-derivation:

    {pcb_id: {"c": [lon, lat],            fixed coordinate (pcb_events lon_c/lat_c)
              "r": 10,                    buffer radius, km (the label set's)
              "mg": credited mergid, "wk": weak_match, "src": pos_source,
              "s": [step_before, step_after],
              "fires": [mergid, ...],     fires with active pixels (credited one first)
              "p": [[dlon, dlat, hx, hy, side, fire_ix], ...]}}
                                          dlon/dlat: offset from "c" in 1e-5 deg; hx/hy: footprint half
                                          size in 1e-5 deg (from DS/DT); side 0 = overpass(es) at/before onset, 1 = after

    python build_active_pcb.py [year ...] [--labels NAME]      then   python map_tool.py join

--labels NAME draws a label set made by `python -m pcbcube.relabel` (OUTPUT_DATA/datacube/labels/NAME) instead
of the labels baked into the pixel tables - try another radius / overpass rule on the map without a rebuild.
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
GEO = ["pix_id", "lon", "lat", "DS", "DT", "fire_key", "mergid", "si", "step"]


def pack(o):
    return base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def year_payload(y: int, ev: pd.DataFrame, labels: pd.DataFrame | None = None, radius_km: float = 20) -> dict:
    if labels is None:
        pix = pd.read_parquet(DC / "_cache" / f"pixels_{y}.parquet", columns=COLS, filters=[("label", "=", 1)])
    else:
        geo = pd.read_parquet(DC / "_cache" / f"pixels_{y}.parquet", columns=GEO)
        pix = geo.merge(labels[labels["label"] == 1], on="pix_id")
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
                             np.rint(hx), np.rint(hy), (g["pcb_step_offset"].to_numpy(int) >= 1).astype(int),
                             [fix[int(m)] for m in g["mergid"]]]).astype(int)
        steps = g.sort_values("si").drop_duplicates("si")
        before = steps.loc[steps["pcb_step_offset"] <= 0, "step"]          # 0, -1, ... = at/before onset
        after = steps.loc[steps["pcb_step_offset"] >= 1, "step"]
        out[pid] = {"c": [round(lon0, 5), round(lat0, 5)], "r": float(radius_km), "mg": cred,
                    "wk": bool(e["weak_match"]), "src": str(e["pos_source"]),
                    "s": [" + ".join(before) or None, " + ".join(after) or None],
                    "fires": fires, "p": p.tolist()}
    return out


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def add_new_positions(d: dict, ev: pd.DataFrame) -> int:
    """Put the dataset's matched ('new') pyroCb coordinate on each event of the map payload, replacing the
    old v2/v3/v4 position columns: np_lat/np_lon, np_src, np_r (ladder radius), np_move_km, np_wk, np_n (px in range),
    np_mg (credited mergid). Unmatched events keep the reported point (np_* = null)."""
    pcb = unpack(d["pcb"])
    n = 0
    for f in pcb["features"]:
        p = f["properties"]
        for k in [k for k in p if k.startswith(("v3_", "v4_", "fire2_"))]:
            del p[k]
        e = ev.loc[p["pcb_id"]] if p["pcb_id"] in ev.index else None
        ok = e is not None and bool(e["matched"])
        p["np_lat"] = round(float(e["lat_c"]), 5) if ok else None
        p["np_lon"] = round(float(e["lon_c"]), 5) if ok else None
        p["np_src"] = str(e["pos_source"]) if ok else None
        p["np_r"] = float(e["radius_km"]) if ok else None
        p["np_move_km"] = round(float(e["move_km"]), 2) if ok else None
        p["np_wk"] = bool(e["weak_match"]) if ok else None
        p["np_n"] = int(e["n_pix"]) if ok else None
        p["np_mg"] = int(e["mergid"]) if ok else None
        n += ok
    d["pcb"] = pack(pcb)
    empty = pack({"type": "FeatureCollection", "features": []})
    d["nb2"] = d["nb3"] = d["nb4"] = empty                 # the old matched-fire v2/v3/v4 overlays are gone
    return n


def dem_at(lon: float, lat: float) -> float | None:
    """Copernicus GLO-90 elevation at a point (nearest cell), None if the tile is missing."""
    import math
    import rasterio
    sys.path.insert(0, str(HERE.parent))
    from pcbcube.dem import tile_path
    fp = tile_path(math.floor(lat), math.floor(lon))
    if not fp.exists():
        return None
    with rasterio.open(fp) as r:
        v = next(r.sample([(lon, lat)]))[0]
    return None if (v is None or not np.isfinite(v) or v < -500) else float(v)


def add_heights(d: dict, inv: pd.DataFrame, fire_elev: dict) -> int:
    """Injection-height variables for colouring pyroCbs on the map:
       trop_alt (m ASL, MERRA-2 tropopause), fire_elev_m (median DEM elevation of the event's active pyroCb
       pixels, else the DEM at the plotted point), inj_above_fire = max_inject_alt - fire_elev_m,
       inj_frac_trop = (max_inject_alt - fire_elev_m) / (trop_alt - fire_elev_m) (1 = reached the tropopause),
       inj_lb (injection height is a lower bound)."""
    pcb = unpack(d["pcb"])
    n = 0
    for f in pcb["features"]:
        p = f["properties"]
        pid = p["pcb_id"]
        r = inv.loc[pid] if pid in inv.index else None
        trop = float(r["trop_alt"]) if r is not None and pd.notna(r["trop_alt"]) else None
        inj = float(r["max_inject_alt"]) if r is not None and pd.notna(r["max_inject_alt"]) else None
        el, src = fire_elev.get(pid, (None, None))
        if el is None:
            lon, lat = (p["np_lon"], p["np_lat"]) if p.get("np_lat") is not None else f["geometry"]["coordinates"]
            el, src = dem_at(float(lon), float(lat)), "dem_at_point"
        p["trop_alt"] = round(trop) if trop is not None else None
        p["fire_elev_m"] = round(el) if el is not None else None
        p["fire_elev_src"] = src
        p["inj_above_fire"] = round(inj - el) if inj is not None and el is not None else None
        p["inj_frac_trop"] = (round((inj - el) / (trop - el), 3)
                              if inj is not None and el is not None and trop is not None and trop > el else None)
        p["inj_lb"] = bool(r["inject_is_lower_bound"]) if r is not None and "inject_is_lower_bound" in r else False
        n += p["inj_above_fire"] is not None
    d["pcb"] = pack(pcb)
    return n


def main(years=None, label_set=None):
    ev = pd.read_parquet(DC / "pcb_events.parquet").set_index("pcb_id")
    ev = ev[~ev.index.duplicated()]
    data = json.loads(DATA.read_text(encoding="utf-8"))
    years = years or sorted(int(y) for y in data)
    sys.path.insert(0, str(HERE.parent))
    from pcbcube import config as PC
    from pcbcube.relabel import set_name
    if not label_set:                                   # default: the config's own label set (made by the build)
        label_set = set_name(PC.LABEL_MODE, PC.ACTIVE_PCB_RADIUS_KM, tuple(PC.ACTIVE_PCB_STEP_OFFSETS),
                             PC.LABEL_MAX_DT_H, PC.LABEL_MAX_DIST_KM)
    labels, radius = None, float(PC.ACTIVE_PCB_RADIUS_KM)
    if label_set:
        d = DC / "labels" / label_set
        labels = pd.read_parquet(d / "pixel_labels.parquet")
        st = json.loads((d / "settings.json").read_text())
        if st["mode"] != "active_buffer":
            sys.exit(f"label set {label_set} is mode {st['mode']}; the map layer draws active_buffer sets")
        radius = st["radius_km"]
        print(f"label set {label_set}: radius {radius} km, overpasses {st['offsets']}")
    from pcbcube.pcb import read_pcb
    inv = read_pcb().drop_duplicates("pcb_id").set_index("pcb_id")
    for y in years:
        apc = year_payload(y, ev[ev["onset_utc"].dt.year == y], labels, radius)
        data[str(y)]["apc"] = pack(apc)
        nn = add_new_positions(data[str(y)], ev)
        # fire altitude = median elevation of the event's active pyroCb pixels (the dataset's terrain values)
        el = pd.read_parquet(DC / "_cache" / f"pixels_{y}.parquet", columns=["pix_id", "elev"])
        lab = (labels[labels["label"] == 1][["pix_id", "pcb_id"]] if labels is not None else
               pd.read_parquet(DC / "_cache" / f"pixels_{y}.parquet", columns=["pix_id", "pcb_id", "label"])
               .query("label == 1")[["pix_id", "pcb_id"]])
        lab = lab.merge(el, on="pix_id").assign(pcb_id=lambda t: t["pcb_id"].astype(str).str.split(";")).explode("pcb_id")
        fe = {k: (float(v), "active_pixels_median") for k, v in lab.groupby("pcb_id")["elev"].median().items() if pd.notna(v)}
        nh = add_heights(data[str(y)], inv, fe)
        n = sum(len(v["p"]) for v in apc.values())
        print(f"  {y}: {nh} with height-above-fire, {nn} events with a new coordinate, {len(apc)} events, {n:,} active-pyroCb pixels, {len(data[str(y)]['apc']) / 1e6:.2f} MB", flush=True)
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
    print(f"data.json now {DATA.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("years", type=int, nargs="*")
    ap.add_argument("--labels")
    a = ap.parse_args()
    main(a.years or None, a.labels)
