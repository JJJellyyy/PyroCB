"""Comparison-page layers for the v2 product: S-NPP-only atlas on the 100 km buffered domain (release_snpp_buf) + MODIS
fusion (pcbcube.fusion) + per-pixel weather / terrain (pcbcube.retag). Written next to the other layers as
compare/<layer>_<year>.json.gz and merged into compare/index.json (see build_compare_map.py).

Layers
  snppbuf    S-NPP buffered final perimeters; popup: per-fire medians over the fire's VIIRS detections of FWI, FFMC,
             ISI, BUI, DC (pcbcube.fwi_grid), 2 m temperature, RH, 10 m wind, VPD, CAPE (ERA5 cube), elevation, slope
             (Copernicus GLO-90) and the dominant land-cover group (C3S) - values at the detection times and places
  fused      VIIRS perimeter U MODIS-added cells per fire (fused_<year>.gpkg), same medians + MODIS added area,
             fused_id (fires bridged by MODIS share it) and the number of bridges
  modisadd   only the MODIS-added part of each fire (fused - VIIRS perimeter)
  modisonly  MODIS-only fires (MCD64A1 cells no VIIRS fire reached)
  bridges    points where MODIS connects two VIIRS fires (first contact cell), with the two mergids
Medians skip pixels outside the ERA5 (30 Mar - 31 Oct) / FWI (Mar - Oct) seasons.

    python build_compare_v2.py [--years 2012 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from build_compare_map import OUT, SIMPLIFY_M, _write  # noqa: E402
from pcbcube import config as C  # noqa: E402
from pcbcube.landcover_data import GROUP_NAMES  # noqa: E402

REL = C.INPUT / "atlas_extension" / "release_snpp_buf"
FUS = ROOT / "OUTPUT_DATA" / "fusion" / "release_snpp_buf"
TAG = ROOT / "OUTPUT_DATA" / "pixtags" / "release_snpp_buf"
MED = {"fwi": "fwi", "ffmc": "ffmc", "isi": "isi", "bui": "bui", "dc": "dc", "t2m_C": "t2m", "rh2m": "rh",
       "wind10": "wind10", "vpd_kPa": "vpd", "cape": "cape", "elev": "elev", "slope_deg": "slope"}


def fire_summary(year) -> pd.DataFrame:
    """Per mergid: medians over its VIIRS detections + dominant land-cover group + MODIS cell counts."""
    P = pd.read_parquet(TAG / f"pixels_{year}.parquet", columns=["pix_id", "mergid", "cls"])
    parts = [pd.read_parquet(TAG / f"fwi_{year}.parquet", columns=["pix_id", "fwi", "ffmc", "isi", "bui", "dc"]),
             pd.read_parquet(TAG / f"era5_{year}.parquet", columns=["pix_id", "t2m_C", "rh2m", "wind10", "vpd_kPa", "cape"]),
             pd.read_parquet(TAG / f"terrain_{year}.parquet", columns=["pix_id", "elev", "slope_deg"]),
             pd.read_parquet(TAG / f"landcover_{year}.parquet", columns=["pix_id", "lc_group"])]
    D = P
    for X in parts:
        D = D.merge(X, on="pix_id", how="left")
    V = D[D.cls == "viirs"]
    S = V.groupby("mergid")[list(MED)].median().rename(columns=MED)
    lc = V.dropna(subset=["lc_group"]).groupby("mergid").lc_group.agg(lambda s: int(s.value_counts().idxmax()))
    S["lc_main"] = lc.map(lambda g: GROUP_NAMES[g].replace("_", " "))
    M = D[D.cls != "viirs"].groupby(["mergid", "cls"]).size().unstack(fill_value=0)
    for c in ("modis_in", "modis_add"):
        S[f"n_{c}"] = M[c] if c in M else 0
    S["n_viirs"] = V.groupby("mergid").size()
    for c in ("fwi", "ffmc", "isi", "bui", "dc", "t2m", "rh", "wind10", "vpd", "elev", "slope"):
        S[c] = S[c].round(1)
    S["cape"] = S.cape.round(0)
    return S.reset_index()


def build_year(year):
    out = {}
    S = fire_summary(year)
    summ = [c for c in S.columns if c != "mergid"]

    # S-NPP buffered perimeters
    G = gpd.read_file(REL / "final_perims" / f"final_viirs{year}.gpkg").to_crs(3571)
    G["start"] = [f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in zip(G.tst_year, G.tst_month, G.tst_day)]
    G["end"] = [f"{y:04d}-{m:02d}-{d:02d}" for y, m, d in zip(G.ted_year, G.ted_month, G.ted_day)]
    G["farea"] = G.farea.round(2)
    G = G.merge(S, on="mergid", how="left")
    Pv = G[["mergid", "geometry"]].dissolve("mergid")          # exact VIIRS perimeter per fire, for the MODIS-added part
    Gs = G.copy()
    Gs["geometry"] = shapely.simplify(Gs.geometry.values, SIMPLIFY_M, preserve_topology=True)
    out["snppbuf"] = _write(Gs, "snppbuf", year, ["mergid", "farea", "start", "end", "lcc_final", "domain_edge"] + summ)

    # fused perimeters, MODIS-added parts, MODIS-only fires
    F = gpd.read_file(FUS / f"fused_{year}.gpkg")
    L = pd.read_parquet(FUS / f"links_{year}.parquet")
    nb = pd.concat([L.mergid_a, L.mergid_b]).value_counts()
    Fv = F[F.source == "VIIRS"].copy()
    Fv["n_bridges"] = Fv.mergid.map(nb).fillna(0).astype(int)
    for c in ("area_viirs", "area_modis_add", "area_fused"):
        Fv[c] = Fv[c].round(2)
    Fv["start"] = [f"{int(y):04d}-{int(m):02d}-{int(d):02d}" for y, m, d in zip(Fv.tst_year, Fv.tst_month, Fv.tst_day)]
    Fv["end"] = [f"{int(y):04d}-{int(m):02d}-{int(d):02d}" for y, m, d in zip(Fv.ted_year, Fv.ted_month, Fv.ted_day)]
    Fv = Fv.merge(S, on="mergid", how="left")
    A = Fv[Fv.area_modis_add > 0.01][["mergid", "fused_id", "area_modis_add", "geometry"]].copy()
    vg = Pv.reindex(A.mergid).geometry.to_numpy()

    def diff(a, b):
        if b is None:
            return a
        try:
            return shapely.difference(a, b)
        except shapely.errors.GEOSException:
            return shapely.difference(shapely.make_valid(a), shapely.make_valid(b))
    A["geometry"] = [diff(a, b) for a, b in zip(A.geometry.to_numpy(), vg)]
    Fs = Fv.copy()
    Fs["geometry"] = shapely.simplify(Fs.geometry.values, SIMPLIFY_M, preserve_topology=True)
    out["fused"] = _write(Fs, "fused", year, ["mergid", "fused_id", "area_viirs", "area_modis_add", "area_fused",
                                              "n_bridges", "start", "end"] + summ)
    A["geometry"] = shapely.simplify(A.geometry.values, SIMPLIFY_M / 2, preserve_topology=True)
    out["modisadd"] = _write(A, "modisadd", year, ["mergid", "fused_id", "area_modis_add"])
    Mo = F[F.source == "MCD64A1"][["mergid", "area_fused", "geometry"]].copy()
    Mo["area_fused"] = Mo.area_fused.round(2)
    Mo["geometry"] = shapely.simplify(Mo.geometry.values, SIMPLIFY_M / 2, preserve_topology=True)
    out["modisonly"] = _write(Mo, "modisonly", year, ["area_fused"])

    # bridges as points
    if len(L):
        B = gpd.GeoDataFrame(L.assign(date=(pd.Timestamp(year, 1, 1) + pd.to_timedelta(L.t_first - 1, unit="D")).dt.strftime("%Y-%m-%d")),
                             geometry=gpd.points_from_xy(L.lon, L.lat), crs=4326).to_crs(3571)
        out["bridges"] = _write(B, "bridges", year, ["mergid_a", "mergid_b", "n_contacts", "date"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=int, nargs="+", default=list(range(2012, 2026)))
    a = ap.parse_args()
    avail = {}
    for y in a.years:
        r = build_year(y)
        for L_, v in r.items():
            if v:
                avail.setdefault(L_, []).append(y)
        print(y, ", ".join(f"{k} {v[0]:,} ({v[1]:.1f} MB)" for k, v in r.items() if v), flush=True)
    idx_fp = OUT / "index.json"
    idx = json.loads(idx_fp.read_text()) if idx_fp.exists() else {}
    for L_, ys in avail.items():
        idx[L_] = sorted(set(idx.get(L_, [])) | set(ys))
    idx_fp.write_text(json.dumps(idx))


if __name__ == "__main__":
    main()
