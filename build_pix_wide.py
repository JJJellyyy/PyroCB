"""Widen each pyroCb's detection grid (map_src/data.json `pix`) from +/-24 h to -72 h .. +24 h and add
snap_b3_h / snap_b4_h to the pyroCb properties, so the static "before" view can offer the 3rd and 4th
overpass before onset (nominal half-day steps put those 24-48 h before onset, outside the old grid).

Same construction as pyrocb.active.slider_data / _cores - NFP detections of the event's fire,
clustered, keeping only the cluster that holds the pyroCb - just with an asymmetric, wider time
window. Entries are unique (cell, time), so any window is still a contiguous index range.

Run:  python build_pix_wide.py [year ...]     then   python map_tool.py join
"""
import base64, gzip, json
import geopandas as gpd, numpy as np, pandas as pd
from pyrocb import config as C, io, pipeline, burngrid
from pyrocb import timebins as T
from pyrocb.active import _localise

PRE_H, POST_H = 72.0, 24.0
DATA = "map_src/data.json"
unpack = lambda b: json.loads(gzip.decompress(base64.b64decode(b)).decode())
pack = lambda o: base64.b64encode(gzip.compress(json.dumps(o, separators=(",", ":")).encode(), 9)).decode()


def wide_cells(year, ev):
    ev = ev[(ev["year"] == year) & ev["mergid"].notna()].copy()
    ev["mergid"] = ev["mergid"].astype("int64")
    pts = gpd.GeoSeries(gpd.points_from_xy(ev["lon"], ev["lat"]), crs=C.WGS84).to_crs(C.METRIC_CRS).values
    by = {k: v for k, v in io.read_nfp_year(year).groupby("mergid")}
    out = {}
    for (_, e), pt in zip(ev.iterrows(), pts):
        g = by.get(e["mergid"])
        if g is None or not len(g):
            continue
        t = pd.Timestamp(e["datetime"])
        tt = g["datetime_utc"].values
        sel = (tt >= np.datetime64(t - pd.Timedelta(hours=PRE_H))) & (tt <= np.datetime64(t + pd.Timedelta(hours=POST_H)))
        if not sel.any():
            continue
        s = g[sel]
        x, y = s["x"].to_numpy(float), s["y"].to_numpy(float)
        dx, dy = s["DS"].to_numpy(float) * 500.0, s["DT"].to_numpy(float) * 500.0
        mask, _, _ = _localise(x, y, dx, dy, pt)
        if not mask.any():
            continue
        rel = (pd.to_datetime(s["datetime_utc"]) - t).dt.total_seconds().to_numpy() / 3600.0
        enc = burngrid.encode(x[mask], y[mask], rel[mask] * 60.0, collapse_first=False)
        if enc is not None:
            out[e["pcb_id"]] = enc
    return out


def main(years=None):
    data = json.load(open(DATA, encoding="utf-8"))
    for y, d in sorted(data.items()):
        if years and int(y) not in years:
            continue
        _, _, ev, _ = pipeline.run_year(int(y), verbose=False)
        new = wide_cells(int(y), ev)
        old = unpack(d["pix"])
        n_old, n_new = sum(len(v["di"]) for v in old.values()), sum(len(v["di"]) for v in new.values())
        # events the old grid had but the rebuild lost would silently lose their static view
        missing = sorted(set(old) - set(new))
        d["pix"] = pack({**old, **new})          # keep any old-only entry rather than drop it
        pcb = unpack(d["pcb"])
        look = ev.set_index("pcb_id")
        for f in pcb["features"]:
            p = f["properties"]
            if p["pcb_id"] in look.index:
                r = look.loc[p["pcb_id"]]
                if isinstance(r, pd.DataFrame):
                    r = r.iloc[0]
                offs = T.snapshot_pre_offsets_h(pd.Timestamp(r["datetime"]), float(r["lon"]), 4)
                p["snap_b3_h"], p["snap_b4_h"] = round(offs[2], 2), round(offs[3], 2)
                # b1/b2 must agree with what the map already carries
                if p.get("snap_b1_h") is not None:
                    assert abs(offs[0] - p["snap_b1_h"]) < 0.02 and abs(offs[1] - p["snap_b2_h"]) < 0.02, p["pcb_id"]
        d["pcb"] = pack(pcb)
        print(f"{y}: pix cells {n_old} -> {n_new}, events {len(old)} -> {len(new)}, lost {missing}", flush=True)
    json.dump(data, open(DATA, "w", encoding="utf-8", newline="\n"), separators=(",", ":"))


if __name__ == "__main__":
    import sys
    main([int(a) for a in sys.argv[1:]] or None)
