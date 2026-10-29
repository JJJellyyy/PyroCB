"""PyroCb-region polygons for the HTML map: outline layer + a per-fire tag for the "only inside the polygons" filter.

Polygons: OUTPUT_DATA/pcb_regions/pcb_regions.geojson (python -m pcbcube.pcbregion) -> copied to pcb_regions.json next to
the map (loaded on demand, like permafrost.json).
Tag rg on every fire feature of map_src/data.json (per, perAll, ign, ignAll, bnd, sb): the pcb_region (1 North America,
2 Eurasia) holding the MAJORITY of the fire's VIIRS pixels (dataset tag OUTPUT_DATA/pcb_regions/pcb_region_tag.parquet), else 0.
pyroCb points (pcb) get the polygon at their plotted position (new coordinate, else the report).

    python build_region_tags.py        then   python map_tool.py join
"""
from __future__ import annotations

import base64
import gzip
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from pcbcube import config as C                                    # noqa: E402
from pcbcube import pcbregion                                      # noqa: E402
from pcbcube.atlas import to_metric                                # noqa: E402

DATA = HERE / "map_src" / "data.json"
FIRE_KEYS = ("per", "perAll", "ign", "ignAll", "bnd", "sb")


def unpack(b):
    return json.loads(gzip.decompress(base64.b64decode(b)).decode())


def pack(o):
    return base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def fire_regions() -> dict:
    """{(year, mergid): rg} - region with the majority of the fire's pixels (0 if most pixels are outside)."""
    P = pq.read_table(C.OUT / "pixels.parquet", columns=["pix_id", "year", "mergid"]).to_pandas()
    T = pd.read_parquet(pcbregion.TAG_FP)
    P = P.merge(T[["pix_id", "pcb_region"]], on="pix_id", how="left")
    cnt = P.groupby(["year", "mergid", "pcb_region"]).size().unstack(fill_value=0)
    share_in = 1 - cnt.get(0, 0) / cnt.sum(axis=1)
    best = cnt.drop(columns=[0], errors="ignore").idxmax(axis=1)
    rg = np.where(share_in >= 0.5, best, 0).astype(int)
    return {(int(y), int(m)): int(r) for (y, m), r in zip(cnt.index, rg)}


def main():
    src = pcbregion.OUT / "pcb_regions.geojson"
    shutil.copyfile(src, HERE / "pcb_regions.json")
    rg = fire_regions()
    data = json.loads(DATA.read_text(encoding="utf-8"))
    tot = {}
    for y in sorted(data, key=int):
        for key in FIRE_KEYS:
            if key not in data[y]:
                continue
            fc = unpack(data[y][key])
            for f in fc["features"]:
                p = f["properties"]
                p["rg"] = rg.get((int(y), int(p["mergid"])), 0) if p.get("mergid") is not None else 0
                tot.setdefault(key, [0, 0])[0] += p["rg"] > 0; tot[key][1] += 1
            data[y][key] = pack(fc)
        if "pcb" in data[y]:
            fc = unpack(data[y]["pcb"])
            if fc["features"]:
                lon = np.array([f["properties"].get("np_lon") or f["geometry"]["coordinates"][0] for f in fc["features"]], float)
                lat = np.array([f["properties"].get("np_lat") or f["geometry"]["coordinates"][1] for f in fc["features"]], float)
                x, yy = to_metric(lon, lat)
                for f, r in zip(fc["features"], pcbregion.tag(x, yy)):
                    f["properties"]["rg"] = int(r)
                    tot.setdefault("pcb", [0, 0])[0] += r > 0; tot["pcb"][1] += 1
            data[y]["pcb"] = pack(fc)
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    for k, (a, n) in tot.items():
        print(f"  {k:<7} {a:>7,} of {n:>7,} features inside a polygon ({a / max(n, 1):.0%})")
    print("wrote", DATA, "and pcb_regions.json")


if __name__ == "__main__":
    main()
