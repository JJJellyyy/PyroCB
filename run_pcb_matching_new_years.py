"""Run the map's own v3/v4/v5 matched-fire notebook (make_pcb_matching_clean_notebook.py) for the
EXTENSION years only (2024-2025) and write `pcb_matching_new_years.gpkg/.tsv` next to the original
`pcb_matching_clean.gpkg`, which stays untouched.

Same technique as run_pcb_matching_clean.py (execute the notebook's code cells in one namespace), with
four source patches so the cells read the merged project's layout:
  * SNAP_DIR / NFP_DIR  -> ../INPUT_DATA/VIIRS/{y}/{y}/...
  * read_pcb            -> pcbcube's extended inventory (Confirmed boreal, stable ids), years >= 2024 only,
                           reduced to the 18 columns the original layers carry
  * load_nfp            -> swaps DS/DT for the extension release (it ships them swapped)
  * output names        -> pcb_matching_new_years.*
"""
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
SRC = HERE / "make_pcb_matching_clean_notebook.py"
VIIRS = (HERE.parent / "INPUT_DATA" / "VIIRS").as_posix()
KEEP = ["datetime", "fire_name", "region", "status", "lat", "lon", "pcb_min_bt", "trop_alt", "max_inject_alt",
        "trop_relative_alt", "min_smoke_mass", "max_smoke_mass", "plume_code", "pulse_code", "notes",
        "year", "doy", "pcb_id"]


def patch(cell: str) -> str:
    cell = cell.replace('SNAP_DIR = "{y}/{y}/Snapshot"', f'SNAP_DIR = "{VIIRS}/{{y}}/{{y}}/Snapshot"')
    cell = cell.replace('NFP_DIR = "{y}/{y}/NFP"', f'NFP_DIR = "{VIIRS}/{{y}}/{{y}}/NFP"')
    cell = cell.replace(
        '    d = pd.read_excel(XLSX_FP, sheet_name=0)\n',
        '    from pcbcube.pcb import read_pcb as _rp\n'
        f'    d = _rp()[lambda t: t["year"] >= 2024][{KEEP!r}].reset_index(drop=True)\n')
    cell = cell.replace('    d["pcb_id"] = _make_pcb_ids(d)\n', '')           # keep pcbcube's stable ids
    cell = cell.replace(
        '        d["datetime_utc"] = pd.to_datetime(d["datetime_utc"])\n        # x/y are EPSG:3571 metres everywhere EXCEPT 2021',
        '        d["datetime_utc"] = pd.to_datetime(d["datetime_utc"])\n'
        '        if year >= 2024:\n            d[["DS", "DT"]] = d[["DT", "DS"]].to_numpy()\n'
        '        # x/y are EPSG:3571 metres everywhere EXCEPT 2021')
    cell = cell.replace("pcb_matching_clean", "pcb_matching_new_years")
    return cell


def main():
    text = SRC.read_text(encoding="utf-8")
    cells = [patch(c) for c in re.findall(r'code\(r"""(.*?)"""\)', text, re.DOTALL)]
    print(f"extracted {len(cells)} code cells", flush=True)
    ns = {"__name__": "__pcb_matching_new_years__"}
    t0 = time.time()
    for i, cell in enumerate(cells, 1):
        t1 = time.time()
        try:
            exec(compile(cell, f"<cell {i}>", "exec"), ns)
        except Exception:
            print(f"\n--- FAILED in cell {i}/{len(cells)} ---")
            raise
        print(f"  cell {i}/{len(cells)} ran in {time.time() - t1:5.1f}s", flush=True)
    print(f"\nall {len(cells)} cells ran in {time.time() - t0:.0f}s total")


if __name__ == "__main__":
    import os
    os.chdir(HERE)
    main()
