"""PyroCb bands v2 for the HTML map (tick box "PyroCb bands v2").

From pcbcube.bands_v2 (S-NPP-only buffered atlas + MODIS reference, 2012-2025):
  kind "c"  per kept event x overpass offset 0 ... -6 the closed outline of the event's complex (observed detections of that
            pass only, footprints closed over 2-pixel gaps), with the band metrics (size, power, local 20 km totals, growth vs the
            previous pass, MODIS reference, conditions at the detections, dominant fuel)
  kind "e"  one point per kept event at the corrected pyroCb position, with the 7-pass trajectory (h before onset, area, power,
            MODIS burned area) and why it has no bands (S-NPP offline) when dropped
Output pcb_bands_v2.json next to the map (GeoJSON, WGS84, loaded on demand).

    python build_pcb_bands_v2_layer.py
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
from pcbcube.bands_v2 import OUT                                    # noqa: E402

SIMPLIFY_M = 60.0
FUEL = {"frac_needle_evergreen_5x5": "evergreen needleleaf", "frac_needle_deciduous_5x5": "deciduous needleleaf (larch)",
        "frac_broad_deciduous_5x5": "broadleaf deciduous", "frac_mixed_forest_5x5": "mixed forest", "frac_shrub_5x5": "shrub",
        "frac_grass_herb_5x5": "grass / herbaceous", "frac_sparse_moss_lichen_5x5": "sparse / moss / lichen",
        "frac_wetland_5x5": "wetland", "frac_cropland_5x5": "cropland", "frac_nonfuel_5x5": "non-fuel"}
ROUND = {"area_fp": 2, "area_closed": 2, "area_unified": 2, "frp_mw": 0, "frp_dens_max": 1, "rg_km": 2, "d_pcb_km": 2,
         "pull_other_ratio": 3, "area_fp_r20": 2, "frp_mw_r20": 0, "area_fp_r50": 2, "frp_mw_r50": 0, "d_other_km": 1,
         "mcd64_km2_cpx": 2, "mcd64_add_km2_cpx": 2, "mcd64_km2_r20": 2, "mcd14_frp_cpx": 0, "fwi": 1, "ffmc": 1, "isi": 1,
         "bui": 1, "dc": 0, "t2m_C": 1, "rh2m": 0, "wind10": 1, "vpd_kPa": 2, "cape": 0, "c_haines": 1, "blh": 0, "elev": 0,
         "slope_deg": 1, "h_before_onset": 1}
INTS = ["n_pix", "n_complexes_pull10", "mcd14_n_cpx", "mcd14_n_r20", "n_fire_ids"]


def main():
    E = pd.read_parquet(OUT / "events.parquet")
    E = E[E["kept"]].copy()
    B = pd.read_parquet(OUT / "bands.parquet")
    G = pd.read_parquet(OUT / "bands_geom.parquet")
    yr = E.set_index("pcb_id").onset_utc.dt.year

    # complexes
    fr = B[list(FUEL)].to_numpy(float)
    ok = np.isfinite(fr).any(axis=1)
    j = np.where(ok, np.nanargmax(np.where(np.isfinite(fr), fr, -1), axis=1), 0)
    B["fuel_main"] = np.where(ok, [list(FUEL.values())[k] for k in j], None)
    B["fuel_share"] = np.where(ok, np.round(fr[np.arange(len(fr)), j] * 100), np.nan)
    for v in ("area_fp", "frp_mw"):                                   # growth factor vs the previous pass
        B[f"{v}_x_prev"] = np.round(10 ** B[f"{v}_lr_prev"], 2)
    B["t_pass"] = pd.to_datetime(B.t_pass).dt.strftime("%Y-%m-%d %H:%M")
    keep = ["pcb_id", "offset", "step", "t_pass", "picked_nearest", "fuel_main", "fuel_share", "area_fp_x_prev", "frp_mw_x_prev"] \
        + list(ROUND) + INTS
    Bc = B[keep].copy()
    for c, n in ROUND.items():
        Bc[c] = Bc[c].astype(float).round(n)
    C = G.merge(Bc, on=["pcb_id", "offset"])
    geoms = [shapely.simplify(shapely.from_wkb(w), SIMPLIFY_M, preserve_topology=True) for w in C.pop("geometry")]
    gc = gpd.GeoDataFrame(C.rename(columns={"offset": "off"}), geometry=geoms, crs=3571).to_crs(4326)
    gc["kind"] = "c"
    gc["year"] = gc.pcb_id.map(yr).astype(int)

    # events with their trajectory
    T = B.sort_values(["pcb_id", "offset"])
    traj = {p: [[int(r.offset), *(None if pd.isna(x) else round(float(x), n) for x, n in
                                 ((r.h_before_onset, 1), (r.area_fp, 2), (r.area_closed, 2), (r.frp_mw, 0), (r.mcd64_km2_cpx, 2)))]
                for r in g.itertuples()] for p, g in T.groupby("pcb_id")}
    E["onset"] = E.onset_utc.dt.strftime("%Y-%m-%d %H:%M")
    E["year"] = E.onset_utc.dt.year
    E["reason"] = np.select([E.window_gap, E.window_month_edge], ["S-NPP offline (gap day) in the window",
                                                                  "lost month-edge night pass in the window"], "")
    E["traj"] = E.pcb_id.map(traj)
    ge = gpd.GeoDataFrame(E[["pcb_id", "year", "onset", "region", "fire_name", "status", "move_km", "mergid", "dropped_snpp",
                             "reason", "window_partial_unverified", "no_preonset_look", "max_inject_alt", "traj"]],
                          geometry=gpd.points_from_xy(E.lon_c, E.lat_c), crs=4326)
    ge["move_km"] = ge.move_km.round(1)
    ge["kind"] = "e"
    ge["mergid"] = ge.mergid.astype("Int64")

    gdf = pd.concat([gc.sort_values("off"), ge], ignore_index=True)       # oldest pass first, onset on top, points last
    gdf = gpd.GeoDataFrame(gdf, geometry="geometry", crs=4326)
    gdf["geometry"] = shapely.transform(gdf.geometry.values, lambda a: np.round(a, 4))
    gdf = gdf[~gdf.geometry.is_empty]
    out = HERE / "pcb_bands_v2.json"
    out.write_text(gdf.to_json(drop_id=True, na="drop"), encoding="utf-8")
    print(f"{(gdf.kind == 'c').sum():,} complexes, {(gdf.kind == 'e').sum()} events ({int(ge.dropped_snpp.sum())} dropped for S-NPP); "
          f"{out.stat().st_size / 1e6:.1f} MB -> {out}")


if __name__ == "__main__":
    main()
