"""Second payload split: whole-fire animation grids, and dropping the dead v2 overlays.

sim -> sim (pyroCb-linked fires) + simAll (everything else). A non-linked fire's grid
is only reachable by clicking its ignition, which requires the "All fire activity"
toggle, so simAll rides along with perAll/ignAll instead of being parsed at startup.
That is ~97% of the sim payload for a typical year.

nb2 (matched-fire-v2 circles/polygons) is emptied: the v2 position option was removed
from the UI, so PCB_POS can never be "fire2" and rebuildNb2() returns immediately -
the payload was pure dead weight.
"""
import base64, gzip, json, shutil

def unpack(b): return json.loads(gzip.decompress(base64.b64decode(b)).decode())
def pack(o):
    raw = json.dumps(o, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw, 9)).decode()

path = "fire_pcb_map.html"
text = open(path, encoding="utf-8").read()
i = text.index("const DATA=")
data, end = json.JSONDecoder().raw_decode(text[i + len("const DATA="):])

EMPTY = pack({"type": "FeatureCollection", "features": []})
for y in sorted(data):
    per = unpack(data[y]["per"])                 # linked-only after the first split
    linked = {str(f["properties"]["mergid"]) for f in per["features"]}
    sim = unpack(data[y]["sim"])
    keep = {k: v for k, v in sim.items() if k in linked}
    rest = {k: v for k, v in sim.items() if k not in linked}
    before = len(data[y]["sim"])
    data[y]["sim"] = pack(keep)
    data[y]["simAll"] = pack(rest)
    nb2_before = len(data[y].get("nb2", ""))
    data[y]["nb2"] = EMPTY
    print(f"  {y} sim: {len(keep):>3} linked + {len(rest):>5} other | "
          f"{before/1e6:4.2f} -> {len(data[y]['sim'])/1e6:4.2f} eager + "
          f"{len(data[y]['simAll'])/1e6:4.2f} lazy MB  | nb2 {nb2_before/1e6:4.2f} -> 0", flush=True)

new = text[:i] + "const DATA=" + json.dumps(data, separators=(",", ":")) + text[i + len("const DATA=") + end:]
shutil.copy2(path, path + ".simbak")
open(path, "w", encoding="utf-8").write(new)
print("\nwrote", path)
