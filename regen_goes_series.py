"""Fetch the GOES time series (one frame per C.GOES_SERIES_OFFSETS_H) for every
eligible pyroCb event, then splice it into fire_pcb_map.html's DATA as a new
per-year "goesq" payload.

Kept separate from regen_all_fires.py because this one is network-bound (a STAC
search plus a windowed read per frame) rather than CPU-bound, and because the JPEGs
it writes live next to the existing stills in cloud_images/ rather than inside the
HTML - only the small metadata index is embedded.

Frames that find no usable daylight scene are absent from the list, so the map's
player must treat the series as irregular and use each frame's own offset_h.
"""
from __future__ import annotations

import base64
import gzip
import json
import os
import shutil
import time

from pyrocb import config as C
from pyrocb import goes_cloud as G
from pyrocb import pipeline

YEARS = [y for y in range(2013, 2024) if y >= C.GOES_MIN_YEAR]
CACHE = "goes_series_cache.json"


def main():
    done = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    for y in YEARS:
        if str(y) in done:
            print(f"{y}: cached ({len(done[str(y)])} events)", flush=True)
            continue
        t0 = time.time()
        _, _, ev, _ = pipeline.run_year(y, verbose=False)
        series = G.event_series(y, ev, C.GOES_IMG_DIR.format(y=y))
        done[str(y)] = series
        json.dump(done, open(CACHE, "w"))
        print(f"{y}: {len(series)} events in {time.time()-t0:.0f}s", flush=True)

    path = "fire_pcb_map.html"
    text = open(path, "r", encoding="utf-8").read()
    i = text.index("const DATA=")
    data, end = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])

    total = 0
    for y in data:
        s = done.get(y, {})
        total += sum(len(v["frames"]) for v in s.values())
        raw = json.dumps(s, separators=(",", ":")).encode()
        data[y]["goesq"] = base64.b64encode(gzip.compress(raw, 9)).decode()

    new = text[:i] + "const DATA=" + json.dumps(data, separators=(",", ":")) \
        + text[i + len("const DATA=") + end:]
    shutil.copy2(path, path + ".goesbak")
    open(path, "w", encoding="utf-8").write(new)
    print(f"\nspliced {total} frames into {path} (backup at {path}.goesbak)")


if __name__ == "__main__":
    main()
