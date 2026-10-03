"""Drop years from the map (2026-10-03: back to 2013-2023), keeping the hand-edited JS.

Removes DATA[year] from map_src/data.json and rewrites the global constants in tail.html for the remaining
years: YEARS, YCOL (RdYlBu reversed over the years left), MONTH_COUNTS, SEARCH, SPREAD_CLAMP. The full
13-year map is saved in map_backup_2013_2025/.

    python rollback_map_years.py 2024 2025      then   python map_tool.py join
"""
import base64
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from merge_payload import _decode_after, year_colours, unpack   # noqa: E402
from pyrocb import spread as SP                                    # noqa: E402

DATA = HERE / "map_src" / "data.json"
TAIL = HERE / "map_src" / "tail.html"


def main(drop):
    drop = {str(y) for y in drop}
    data = json.loads(DATA.read_text(encoding="utf-8"))
    for y in drop:
        data.pop(y, None)
    years = sorted(int(y) for y in data)
    t = TAIL.read_text(encoding="utf-8")
    _, a, b = _decode_after(t, ", YEARS=")
    t = t[:a] + json.dumps(years) + t[b:]
    _, a, b = _decode_after(t, "YCOL=")
    t = t[:a] + json.dumps(year_colours(years)) + t[b:]
    mc, a, b = _decode_after(t, "MONTH_COUNTS=")
    t = t[:a] + json.dumps({k: v for k, v in mc.items() if k not in drop}) + t[b:]
    se, a, b = _decode_after(t, "SEARCH=")
    t = t[:a] + json.dumps([r for r in se if str(r["y"]) not in drop]) + t[b:]
    _, a, b = _decode_after(t, "const SPREAD_CLAMP=")
    t = t[:a] + json.dumps(SP.clamp({y: unpack(data[y]["spr"]) for y in data})) + t[b:]
    DATA.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8", newline="\n")
    TAIL.write_text(t, encoding="utf-8", newline="\n")
    print(f"map years now {years}; data.json {DATA.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]])
