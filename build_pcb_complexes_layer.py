"""Burning complexes of every kept pyroCb for the HTML map (tick box "PyroCb burning complexes").

From pcbcube.complexes (python -m pcbcube.complexes): per pyroCb and overpass offset 0 ... -4 the complex of still-burning pixels
with the strongest pull on the pyroCb (pixels linked with at most one undetected pixel between them; weights 1 / 1 / 0.5 / 0.25 for
ages 0 / 12 / 24 / 36 h). One feature per pyroCb x offset x age class: the union of those pixels' footprints, closed by CLOSE_M and
simplified for display. Properties: pcb_id, off, age, weight, step, year, after_prev (overpass within GAP_DAYS after an earlier
pyroCb of the same atlas fire) and the complex metrics (area_w, frp_new, frp_w, n_pix_active, share_credited_power,
frp_new_credited_fire, pull_other_ratio, n_complexes_pull10).
Output pcb_complexes.json next to the map (GeoJSON, WGS84, loaded on demand).

    python build_pcb_complexes_layer.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pcbcube import prepulse                                       # noqa: E402
from pcbcube.atlas import StepIndex                                # noqa: E402
from pcbcube.complexes import GEOM_FP, OUT_FP                      # noqa: E402

CLOSE_M, SIMPLIFY_M = 0.0, 60.0          # geometry is already closed in pcbcube.complexes
MET = ["area_unified", "unify_radius_km", "cen_dist_km", "area_closed_w", "area_closed", "area_w", "frp_new", "frp_w", "n_pix_active", "share_credited_power", "frp_new_credited_fire", "pull_other_ratio",
       "n_complexes_pull10"]


def main():
    R = pd.read_parquet(OUT_FP)
    R = R[R["kept"]]
    G = pd.read_parquet(GEOM_FP).merge(R[["pcb_id", "offset", "step", "year", "si", "fire_key"] + MET], on=["pcb_id", "offset"])
    st = StepIndex()
    ev = prepulse.events(prepulse.GAP_DAYS)
    on = ev.groupby("fire_key")["onset_utc"].apply(lambda s: np.sort(s.to_numpy("datetime64[ns]")))
    onset = ev.set_index("pcb_id")["onset_utc"]

    def after_prev(r):
        o = on.get(r.fire_key)
        if o is None:
            return False
        prev = o[o < onset[r.pcb_id].to_datetime64()]
        t = st.t0[int(r.si)].to_datetime64()
        return bool(len(prev) and t >= prev.max() and (t - prev.max()) / np.timedelta64(1, "h") / 24 <= prepulse.GAP_DAYS)

    G["after_prev"] = [after_prev(r) for r in G.itertuples()]
    geoms = []
    for w in G["wkb"]:
        g = shapely.from_wkb(w)
        geoms.append(shapely.simplify(g, SIMPLIFY_M, preserve_topology=True))
    gdf = gpd.GeoDataFrame(G.drop(columns=["wkb", "si", "fire_key"]).rename(columns={"offset": "off"}), geometry=geoms, crs=3571).to_crs(4326)
    gdf["geometry"] = shapely.transform(gdf.geometry.values, lambda a: np.round(a, 4))
    gdf = gdf[~gdf.geometry.is_empty].sort_values(["off", "age"], ascending=[True, False])      # oldest first, onset / newest on top
    for c in MET:
        gdf[c] = gdf[c].astype(float).round(3)
    # position option "centre of the onset polygon": small lookup {pcb_id: [lat, lon]} (all matched pyroCbs with an onset complex)
    from pyproj import Transformer
    c = R[(R["offset"] == 0)].dropna(subset=["x_cen"]).drop_duplicates("pcb_id")
    lon, lat = Transformer.from_crs(3571, 4326, always_xy=True).transform(c["x_cen"].to_numpy(), c["y_cen"].to_numpy())
    import json
    (HERE / "pcb_cpxpos.json").write_text(json.dumps({p: [round(float(la), 5), round(float(lo), 5)] for p, la, lo in zip(c["pcb_id"], lat, lon)}),
                                          encoding="utf-8")
    out = HERE / "pcb_complexes.json"
    out.write_text(gdf.to_json(drop_id=True), encoding="utf-8")
    print(f"{len(gdf):,} features ({gdf.pcb_id.nunique()} pyroCbs); {out.stat().st_size / 1e6:.1f} MB -> {out}")


if __name__ == "__main__":
    main()
