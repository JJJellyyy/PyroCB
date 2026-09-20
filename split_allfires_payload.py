"""Split each year's per/ign payload into the pyroCb-linked analysis population and
the "All fire activity" background set, so the map can load the first eagerly and the
second only if someone actually turns the toggle on.

Shipping every tracked fire in the SAME payload cost ~15 s of Leaflet layer
construction at startup (10.7k polygons + 12.2k DOM star markers for 2023 alone),
every time, whether or not the toggle was ever touched. After this split:

    per  / ign     pyroCb-linked only, full precision  - unpacked at load, as before
    perAll/ignAll  everything else, coarsely simplified - untouched until toggled on

Non-linked perimeters are background context, never the analysis geometry, so they are
simplified 8x coarser and snapped to a ~220 m grid: 2023 goes from 3.80 MB to 0.87 MB
gzipped with no visible difference at the zoom levels they are ever seen at.

Works off the payloads already embedded in fire_pcb_map.html (every feature already
carries pcb_linked), so it needs no pipeline run.
"""
from __future__ import annotations

import base64
import gzip
import json
import shutil

import geopandas as gpd
import shapely

from pyrocb import config as C

COARSE_TOL = C.SIMPLIFY_TOL * 8          # ~265 m
COARSE_PREC = 2e-3                       # ~220 m


def unpack(b64):
    return json.loads(gzip.decompress(base64.b64decode(b64)).decode())


def pack(obj):
    raw = json.dumps(obj, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode(), len(raw)


def split(fc, coarsen):
    linked = [f for f in fc["features"] if f["properties"].get("pcb_linked")]
    other = [f for f in fc["features"] if not f["properties"].get("pcb_linked")]
    if coarsen and other:
        g = gpd.GeoDataFrame.from_features(other, crs=C.WGS84)
        g["geometry"] = g.geometry.simplify(COARSE_TOL, preserve_topology=True)
        g["geometry"] = shapely.set_precision(g.geometry.values, COARSE_PREC)
        g = g[~g.geometry.is_empty & g.geometry.notna()]
        other = json.loads(g.to_json(drop_id=True))["features"]
    mk = lambda fs: {"type": "FeatureCollection", "features": fs}
    return mk(linked), mk(other)


def main():
    path = "fire_pcb_map.html"
    text = open(path, "r", encoding="utf-8").read()
    i = text.index("const DATA=")
    data, end = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])

    tot_before = tot_after = 0
    for y in sorted(data):
        for key, coarsen in (("per", True), ("ign", False)):
            fc = unpack(data[y][key])
            before = len(data[y][key])
            keep, rest = split(fc, coarsen)
            data[y][key], _ = pack(keep)
            data[y][key + "All"], _ = pack(rest)
            after = len(data[y][key]) + len(data[y][key + "All"])
            tot_before += before; tot_after += after
            if key == "per":
                print(f"  {y} per: {len(keep['features']):>3} linked + "
                      f"{len(rest['features']):>6} other | "
                      f"{before/1e6:5.2f} MB -> {len(data[y][key])/1e6:4.2f} (eager) + "
                      f"{len(data[y][key+'All'])/1e6:4.2f} (lazy) MB b64", flush=True)

    new = text[:i] + "const DATA=" + json.dumps(data, separators=(",", ":")) \
        + text[i + len("const DATA=") + end:]
    shutil.copy2(path, path + ".splitbak")
    open(path, "w", encoding="utf-8").write(new)
    print(f"\nper+ign base64 total: {tot_before/1e6:.1f} MB -> {tot_after/1e6:.1f} MB")
    print(f"wrote {path} (backup at {path}.splitbak)")


if __name__ == "__main__":
    main()
