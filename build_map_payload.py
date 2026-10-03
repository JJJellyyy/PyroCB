"""Regenerate the map's DATA payload (all years incl. the 2024-25 extension) WITHOUT touching the hand-edited
map_src/head.html and tail.html.

The notebook `Data exploration.ipynb` writes the whole HTML from a template that is older than the hand edits,
so it is never run to produce fire_pcb_map.html any more. This script executes its code cells in one namespace
with two patches and then leaves the notebook's own HTML in a TEMP file:

  * the matched-fire v3/v4 layers are read from BOTH pcb_matching_clean.gpkg (2013-2023, original algorithm run)
    and pcb_matching_new_years.gpkg (2024-25, same algorithm - see run_pcb_matching_new_years.py)
  * OUT_HTML points at `cache/_payload_probe.html`

`merge_payload.py` then splits that file and folds the new DATA / global constants into map_src.

    python build_map_payload.py 2014 2015 2022 2023 2024 2025     # the years to regenerate

Cells 18-22 (narrative QC and pipeline.export, which would overwrite the all-years attribute tables with
a subset) are skipped.
"""
import json
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
NB = HERE / "Data exploration.ipynb"
PROBE = HERE / "cache" / "_payload_probe.html"


def patch(i: int, src: str) -> str:
    if i == 6:   # matched-fire v3/v4 layers: original file + the extension-years file
        src = src.replace(
            'v3_gdf = gpd.read_file(C.ROOT / "pcb_matching_clean.gpkg", layer="pcb_active_burn",',
            'def _both(layer, **kw):\n'
            '    parts = [gpd.read_file(C.ROOT / f, layer=layer, **kw)\n'
            '             for f in ("pcb_matching_clean.gpkg", "pcb_matching_new_years.gpkg") if (C.ROOT / f).exists()]\n'
            '    return pd.concat(parts, ignore_index=True)\n'
            'v3_gdf = _both("pcb_active_burn",')
        src = src.replace('v4_gdf = gpd.read_file(C.ROOT / "pcb_matching_clean.gpkg", layer="pcb_active_burn_v4",',
                          'v4_gdf = _both("pcb_active_burn_v4",')
        assert "_both(" in src and src.count("_both(") == 3
    return src


SKIP = {18, 20, 22}


def main(years):
    nb = json.load(open(NB, encoding="utf-8"))
    ns = {"__name__": "__map_payload__"}
    t0 = time.time()
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code" or i in SKIP:
            continue
        src = patch(i, "".join(c["source"]))
        if i == 2:
            src += f"\nC.OUT_HTML = C.CACHE_DIR / '_payload_probe.html'\nC.YEARS = {years!r}\n"
        t1 = time.time()
        try:
            exec(compile(src, f"<nb cell {i}>", "exec"), ns)
        except Exception:
            print(f"--- FAILED in notebook cell {i} ---", flush=True)
            raise
        print(f"cell {i} ran in {time.time() - t1:.0f}s", flush=True)
    print(f"BUILD-PAYLOAD-DONE {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]])
