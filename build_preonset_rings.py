"""Pre-onset 12 h ACTIVE AREAS of every kept pyroCb for the HTML map (tick box "Pre-onset 12 h rings").

Per kept pyroCb (pcbcube.prepulse.events) and per overpass at offsets 0 (onset overpass) ... -4 (48 h before), colours as in the
power-before-onset figures (fig5b):
  geometry  = the credited fire's 12 h growth ring UNITED with the footprints (DS x DT) of all its new VIIRS pixels at that overpass
              (about half of a fire's new pixels lie inside the old perimeter and add no ring area - they still count as burning)
  kept parts: only the separate parts (polygons) that come within NEAR_KM (20 km) of the pyroCb coordinate (corrected position);
              a part may extend far beyond 20 km, but parts entirely farther away (other branches of a merged fire) are dropped
  after_prev: overpasses within GAP_DAYS after an earlier pyroCb of the same fire are KEPT (they were dropped before) and flagged
              (dashed outline on the map); the power-before-onset figures still leave them out.
Properties: pcb_id, off, step, year, after_prev, area_local (km2, kept parts), area_ring_local (growth-ring part of it), n_pix_local
and frp_local (pixels in the kept parts), n_pix_fire and frp_fire (all the fire's pixels at that overpass), n_parts_dropped.
Display outline: footprints closed by CLOSE_M and simplified (areas computed before smoothing).
Output pcb_preonset.json next to the map (GeoJSON, WGS84, loaded on demand).

    python build_preonset_rings.py        (no map_tool join needed: the file is fetched at run time)
"""
from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pcbcube import config as C                                    # noqa: E402
from pcbcube import prepulse                                       # noqa: E402
from pcbcube.atlas import StepIndex                                # noqa: E402
from pcbcube.geometry import pixel_footprints                      # noqa: E402
from pcbcube.rings import rings_path                               # noqa: E402

OFFSETS = [0, -1, -2, -3, -4]
NEAR_KM = 20.0
CLOSE_M = 150.0
SIMPLIFY_M = 60.0


def main():
    ev = prepulse.events(prepulse.GAP_DAYS)
    on = ev.groupby("fire_key")["onset_utc"].apply(lambda s: np.sort(s.to_numpy("datetime64[ns]")))
    kept = ev[ev["kept"] & ev["fire_key"].notna()]
    st = StepIndex()
    feats, tot = [], dict(bands=0, empty_near=0, after_prev=0, parts_dropped=0)
    for y in sorted(kept["onset_utc"].dt.year.unique()):
        K = kept[kept["onset_utc"].dt.year == y]
        G = gpd.read_parquet(rings_path(int(y)))
        G = G[G["fire_key"].isin(K["fire_key"].astype("int64")) & (G["ring_area_km2"] > 0)].set_index(["fire_key", "si"])
        sis = sorted({int(s) + o for s in K["si_grid"] for o in OFFSETS})
        P = pq.read_table(C.OUT / "pixels.parquet", columns=["fire_key", "si", "x", "y", "DS", "DT", "frp", "pix_area_km2"],
                          filters=[("year", "=", int(y)), ("si", "in", sis)]).to_pandas()
        P = P[P["fire_key"].isin(K["fire_key"].astype("int64"))]
        P["frp_mw"] = P["frp"] * P["pix_area_km2"]
        Pg = {k: g for k, g in P.groupby(["fire_key", "si"])}
        for _, e in K.iterrows():
            fk, sg = int(e["fire_key"]), int(e["si_grid"])
            prev = on[fk][on[fk] < e["onset_utc"].to_datetime64()]
            zone = shapely.Point(e["x"], e["y"]).buffer(NEAR_KM * 1000.0)
            for k in OFFSETS:
                si = sg + k
                t = st.t0[si].to_datetime64()
                after = bool(len(prev) and t >= prev.max() and (t - prev.max()) / np.timedelta64(1, "h") / 24 <= prepulse.GAP_DAYS)
                ring = G.loc[(fk, si), "geometry"] if (fk, si) in G.index else None
                q = Pg.get((fk, si))
                fp = pixel_footprints(q.x, q.y, q.DS, q.DT) if q is not None and len(q) else None
                parts_in = [g for g in (ring, fp) if g is not None and not g.is_empty]
                if not parts_in:
                    continue
                act = shapely.union_all(parts_in)
                parts = list(shapely.get_parts(act))
                keep = [p for p in parts if p.intersects(zone)]
                tot["parts_dropped"] += len(parts) - len(keep)
                if not keep:
                    tot["empty_near"] += 1
                    continue
                local = shapely.union_all(keep)
                ring_local = ring.intersection(local).area / 1e6 if ring is not None else 0.0
                if q is not None and len(q):
                    inl = shapely.within(shapely.points(q.x.to_numpy(), q.y.to_numpy()), local.buffer(1.0))
                    n_loc, f_loc, n_all, f_all = int(inl.sum()), float(q.frp_mw[inl].sum()), len(q), float(q.frp_mw.sum())
                else:
                    n_loc = n_all = 0; f_loc = f_all = 0.0
                disp = shapely.union_all([p.buffer(CLOSE_M, quad_segs=2).buffer(-CLOSE_M, quad_segs=2) for p in keep])
                disp = shapely.simplify(disp if not disp.is_empty else local, SIMPLIFY_M, preserve_topology=True)
                feats.append(dict(geometry=disp, pcb_id=e["pcb_id"], off=k, step=str(st.steps[si]), year=int(y), after_prev=after,
                                  area_local=round(local.area / 1e6, 2), area_ring_local=round(ring_local, 2),
                                  n_pix_local=n_loc, frp_local=round(f_loc, 1), n_pix_fire=n_all, frp_fire=round(f_all, 1),
                                  n_parts_dropped=len(parts) - len(keep)))
                tot["bands"] += 1; tot["after_prev"] += after
    gdf = gpd.GeoDataFrame(feats, geometry="geometry", crs=3571).to_crs(4326)
    gdf["geometry"] = shapely.transform(gdf.geometry.values, lambda a: np.round(a, 4))
    gdf = gdf[~gdf.geometry.is_empty]
    out = HERE / "pcb_preonset.json"
    out.write_text(gdf.to_json(drop_id=True), encoding="utf-8")
    loc_share = (gdf["frp_local"] / gdf["frp_fire"].where(gdf["frp_fire"] > 0)).median()
    print(f"{len(gdf):,} active-area bands of {kept['pcb_id'].nunique()} kept pyroCbs; {tot}; median share of the fire's power in the "
          f"kept (<= {NEAR_KM:g} km) parts {loc_share:.0%}; {out.stat().st_size / 1e6:.1f} MB -> {out}")
    print(gdf.groupby("off").size().to_dict())


if __name__ == "__main__":
    main()
