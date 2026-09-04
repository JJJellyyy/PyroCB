"""Readers for the four FireAtlas products, plus the geometry validity pipeline.

The single join key across all of them is `mergid`:
    ignitions.mergid == final.fireid == final.mergid == NFP.mergid == snapshot.mergid
Snapshot `fireid` is incomplete (~5% of rows are absent from the final file) and can
attach a wrong minor polygon; NFP `fireid` is just a per-file row counter.  Neither is
ever used here.
"""
from __future__ import annotations

import glob
import os
from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from . import config as C
from .timebins import step_sort_key

_POLY = {"Polygon", "MultiPolygon"}


# ---------------------------------------------------------------- geometry validity
def ensure_valid(geoms):
    """Repair invalid / non-polygonal geometry. Cheap: no-ops when input is clean.

    Never `buffer(0)`: on a bowtie it silently deletes the smaller loop and on some
    MultiPolygons drops whole lobes, i.e. it changes area without saying so.  That is
    fatal when the output *is* an area.  `make_valid` preserves area.

    Measured: 0 of 1413 snapshot geometries and 4 of 10707 final perimeters are
    invalid, so both branches below almost always skip.
    """
    g = np.asarray(geoms, dtype=object)
    bad = ~shapely.is_valid(g)
    if bad.any():
        g = g.copy()
        g[bad] = shapely.make_valid(g[bad])
    types = shapely.get_type_id(g)          # 3 = Polygon, 6 = MultiPolygon
    odd = ~np.isin(types, (3, 6))
    if odd.any():
        g = g.copy()
        g[odd] = [_polygonal(x) for x in g[odd]]
    return g


def snap(geoms):
    """Snap to the 1 m grid, for geometries about to enter a boolean.

    Snapping both operands of a difference to a common grid is what makes it
    topologically exact and kills hairline slivers; 1 m is 1/375 of a VIIRS pixel.
    Costs ~2 s per 1400 geometries, so it is applied only to the few fires that
    actually get banded, never to a whole snapshot.
    """
    g = shapely.set_precision(np.asarray(geoms, dtype=object), C.PRECISION_M,
                              mode="valid_output")
    return ensure_valid(g)


def _polygonal(g):
    """Keep only polygonal parts; make_valid can emit GeometryCollections."""
    if g is None or g.is_empty:
        return shapely.Polygon()
    if g.geom_type in _POLY:
        return g
    parts = [p for p in shapely.get_parts(g) if p.geom_type in _POLY]
    if not parts:
        return shapely.Polygon()
    return parts[0] if len(parts) == 1 else shapely.union_all(parts)


def drop_slivers(g, min_km2: float = C.SLIVER_KM2):
    """Explode, drop parts under `min_km2`, re-collect. Returns (geom, dropped_km2)."""
    if g is None or g.is_empty:
        return g, 0.0
    keep, dropped = [], 0.0
    for p in shapely.get_parts(g):
        a = p.area / 1e6
        if a >= min_km2:
            keep.append(p)
        else:
            dropped += a
    if not keep:
        return shapely.Polygon(), dropped
    return (keep[0] if len(keep) == 1 else shapely.union_all(keep)), dropped


def km2(g) -> float:
    return 0.0 if g is None or g.is_empty else g.area / 1e6


# ---------------------------------------------------------------- steps
@lru_cache(maxsize=32)
def list_steps(year: int) -> tuple[str, ...]:
    """Sorted snapshot step names for a year, excluding the *_FL fire-line files."""
    d = C.SNAP_DIR.format(y=year)
    names = [os.path.basename(p)[:-5] for p in glob.glob(os.path.join(d, "*.gpkg"))]
    names = [n for n in names if not n.endswith("_FL")]
    return tuple(sorted(names, key=step_sort_key))


@lru_cache(maxsize=32)
def steps_set(year: int) -> frozenset:
    return frozenset(list_steps(year))


# ---------------------------------------------------------------- snapshots
_SNAP_NUM = ["farea", "fperim", "flinelen", "duration", "pixden", "meanFRP", "FRP95",
             "lake_area", "lake_border", "lake_border_tot"]
_SNAP_INT = ["n_pixels", "n_newpixels", "isactive", "t_inactive", "isignition",
             "invalid", "lake_no"]


@lru_cache(maxsize=8)
def read_snapshot(year: int, step: str) -> gpd.GeoDataFrame:
    """One snapshot, validated and dissolved by `mergid`, in METRIC_CRS.

    Dissolving is mandatory: a single file can hold 70 rows sharing a mergid (the
    not-yet-executed components of a merge).  Aggregating by mergid is what makes
    snapshot n_pixels reproduce the cumulative NFP count.
    """
    fp = os.path.join(C.SNAP_DIR.format(y=year), f"{step}.gpkg")
    g = gpd.read_file(fp)
    if g.crs is None or g.crs.to_epsg() != C.METRIC_CRS:
        g = g.to_crs(C.METRIC_CRS)

    # avgspread / spread95 ship as strings and would poison any aggregation.
    for c in ("avgspread", "spread95"):
        if c in g.columns:
            g[c] = pd.to_numeric(g[c], errors="coerce")
    g = g.drop(columns=[c for c in ("mergedate", "ftype") if c in g.columns])

    g["geometry"] = ensure_valid(g.geometry.values)

    agg = {c: "sum" for c in _SNAP_NUM + _SNAP_INT if c in g.columns}
    agg.update({c: "mean" for c in ("meanFRP", "pixden", "avgspread") if c in g.columns})
    agg.update({c: "max" for c in ("FRP95", "spread95", "duration", "isactive") if c in g.columns})
    agg.update({c: "min" for c in ("t_inactive",) if c in g.columns})
    out = g.dissolve(by="mergid", aggfunc=agg)
    out["n_components"] = g.groupby("mergid").size()
    return out.reset_index()


# ---------------------------------------------------------------- final / ignitions
@lru_cache(maxsize=16)
def read_final(year: int, metric: bool = True) -> gpd.GeoDataFrame:
    """Season-final perimeters keyed by `mergid`.

    In this product `fireid == mergid` for every row (asserted), so we keep `mergid`
    as the single key and drop `fireid` to remove any temptation to join on it.
    """
    g = gpd.read_file(C.PERIMS_TPL.format(y=year))
    assert (g["fireid"] == g["mergid"]).all(), f"final {year}: fireid != mergid"
    g = g.drop(columns=["fireid"]).to_crs(C.METRIC_CRS if metric else C.WGS84)
    if metric:
        g["geometry"] = ensure_valid(g.geometry.values)
    return g


@lru_cache(maxsize=16)
def read_ignitions(year: int) -> gpd.GeoDataFrame:
    """Ignition footprints in METRIC_CRS. One row per `fireid`; `mergid` groups the
    separately-ignited fires that later grew together (1,444 of 12,151 rows in 2023).

    Adds two shared columns so the map and `variables.ignition_vars` agree:
      ign_t      nominal ignition instant (tst_* half-day at the polygon's own
                 centroid longitude -> local solar -> UTC)
      ign_index  1-based chronological rank of the ignition within its `mergid`,
                 ties (same half-day) broken on fireid. Its max == n_ign_total.
    """
    from . import timebins as T

    g = gpd.read_file(C.IGN_TPL.format(y=year))
    if g.crs is None or g.crs.to_epsg() != C.METRIC_CRS:
        g = g.to_crs(C.METRIC_CRS)

    lon = g.geometry.representative_point().to_crs(C.WGS84).x.values
    d = pd.to_datetime(dict(year=g["tst_year"], month=g["tst_month"], day=g["tst_day"]),
                       errors="coerce")
    h = g["tst_ampm"].map(T.NOMINAL_LOCAL_H).fillna(7.0).values
    g["ign_t"] = d + pd.to_timedelta(h - T.norm_lon(lon) / 15.0, unit="h")
    g["ign_index"] = (g.sort_values(["ign_t", "fireid"])
                        .groupby("mergid").cumcount().add(1)
                        .reindex(g.index))
    return g


# ---------------------------------------------------------------- NFP
# DS/DT are the along-scan / along-track pixel dimensions in km; needed to give each
# detection its real footprint when building the active-burn polygon.
_NFP_COLS = ["mergid", "y", "x", "DS", "DT", "frp", "sat", "datetime_utc"]


def _read_nfp_raw(year: int) -> pd.DataFrame:
    """All new-fire-pixel detections for a year, thinned and cached, as-shipped.

    NFP is incremental (consecutive steps share no coordinates) and carries a real
    `datetime_utc`, which is what makes exact before/after splits possible.
    """
    d = C.NFP_DIR.format(y=year)
    files = sorted(glob.glob(os.path.join(d, "*_NFP.txt")))

    def build():
        frames = []
        for p in files:
            step = os.path.basename(p).replace("_NFP.txt", "")
            try:
                df = pd.read_csv(p, usecols=_NFP_COLS)
            except (ValueError, pd.errors.EmptyDataError):
                continue
            if df.empty:
                continue
            df["step"] = step
            frames.append(df)
        out = pd.concat(frames, ignore_index=True)
        out["datetime_utc"] = pd.to_datetime(out["datetime_utc"])
        out["mergid"] = out["mergid"].astype("int64")
        for c in ("frp", "DS", "DT"):
            out[c] = out[c].astype("float32")
        # The attribute table says x/y are EPSG:3571 metres, and they are - except in
        # 2021, which ships lon/lat. Detect it from the magnitudes rather than trusting
        # the spec, and reproject so every year lands in the same frame.
        if out["x"].abs().max() <= 180 and out["y"].abs().max() <= 90:
            print(f"    NFP {year}: coordinates are lon/lat, not EPSG:{C.METRIC_CRS} "
                  f"- reprojecting", flush=True)
            g = gpd.GeoSeries(shapely.points(out["x"].to_numpy(float),
                                             out["y"].to_numpy(float)),
                              crs=C.WGS84).to_crs(C.METRIC_CRS)
            out["x"], out["y"] = g.x.to_numpy(), g.y.to_numpy()
            out["nfp_reprojected"] = True
        else:
            out["nfp_reprojected"] = False
        return out

    return C.cached("nfp", year, {"cols": _NFP_COLS}, files[:1] + files[-1:], build)


def read_nfp_year(year: int) -> pd.DataFrame:
    """NFP detections with `mergid` validated against the final perimeters.

    The id join is checked, not assumed: if fewer than NFP_ID_MIN_MATCH of fires
    reproduce their `n_pixels`, the year's labels are stale and pixels are
    re-attributed by geometry instead. `nfp_repaired(year)` reports which happened.
    """
    files = sorted(glob.glob(os.path.join(C.NFP_DIR.format(y=year), "*_NFP.txt")))

    def build():
        raw = _read_nfp_raw(year)
        if nfp_id_health(year) >= C.NFP_ID_MIN_MATCH:
            raw = raw.copy()
            raw["nfp_repaired"] = False
            return raw
        fixed = _repair_nfp_by_geometry(year, raw)
        fixed["nfp_repaired"] = True
        print(f"    NFP {year}: stale ids in the shipped files - re-attributed "
              f"{len(fixed):,}/{len(raw):,} pixels by geometry", flush=True)
        return fixed

    return C.cached("nfpv", year, {"cols": _NFP_COLS, "min": C.NFP_ID_MIN_MATCH},
                    files[:1] + files[-1:], build)


def nfp_repaired(year: int) -> bool:
    return bool(read_nfp_year(year)["nfp_repaired"].iloc[0])


def nfp_id_health(year: int) -> float:
    """Fraction of final fires whose NFP pixel count matches `n_pixels` exactly."""
    nfp = _read_nfp_raw(year)
    fin = read_final(year)
    s = nfp.groupby("mergid").size()
    j = fin[["mergid", "n_pixels"]].set_index("mergid").join(s.rename("nfp")).fillna(0)
    return float((j["n_pixels"] == j["nfp"]).mean())


def _repair_nfp_by_geometry(year: int, nfp: pd.DataFrame) -> pd.DataFrame:
    """Re-attribute NFP pixels to fires spatially, ignoring the file's own `mergid`.

    Needed for 2021, where the NFP files were not re-issued after the season was
    reprocessed: they carry 18,931 distinct ids against the final product's 11,688, and
    `final.mergid` is not even a subset of them - yet the pixel TOTAL agrees to 0.04%.
    So the detections are all present and correctly located; only the labels are stale.

    Since the perimeters are themselves built from these pixels, a point-in-polygon
    assignment against the final perimeters recovers the intended mapping.
    """
    fin = read_final(year)
    pts = shapely.points(nfp["x"].to_numpy(float), nfp["y"].to_numpy(float))
    tree = shapely.STRtree(fin.geometry.values)
    pi, gi = tree.query(pts, predicate="within")
    out = nfp.copy()
    out["mergid_file"] = out["mergid"]
    out["mergid"] = pd.NA
    out.iloc[pi, out.columns.get_loc("mergid")] = fin["mergid"].values[gi]
    # a pixel landing in two overlapping perimeters keeps the first; drop unassigned
    out = out[out["mergid"].notna()].copy()
    out["mergid"] = out["mergid"].astype("int64")
    out = out[~out.index.duplicated(keep="first")]
    return out


def nfp_obs_index(year: int) -> pd.DataFrame:
    """Per (mergid, step) observation times: t_min, t_med, t_max, n.

    This is what Tier A of `timebins.step_before` consumes - it lets us pick the last
    overpass that genuinely precedes onset without any local-time assumption.
    """
    def build():
        nfp = read_nfp_year(year)
        g = nfp.groupby(["mergid", "step"])["datetime_utc"]
        out = g.agg(t_min="min", t_max="max", n="size").reset_index()
        out["t_med"] = g.median().values
        return out

    files = sorted(glob.glob(os.path.join(C.NFP_DIR.format(y=year), "*_NFP.txt")))
    return C.cached("obs", year, {"v": 1}, files[:1] + files[-1:], build)


# ---------------------------------------------------------------- pyroCb atlas
def read_pcb() -> pd.DataFrame:
    """The NRL pyroCb inventory, cleaned, with a stable `pcb_id`.

    `lat`/`lon` are the APPROXIMATE FIRE location (readme Table 1), not the plume -
    but stored to 0.1 deg (~11 km), which sets a floor on match quality.
    """
    def build():
        d = pd.read_excel(C.XLSX_FP, sheet_name=0)
        d["datetime"] = pd.to_datetime(d["datetime"], errors="coerce")
        d = d.dropna(subset=["datetime", "lat", "lon"]).copy()
        d["year"] = d["datetime"].dt.year
        d["doy"] = d["datetime"].dt.dayofyear
        d = d.sort_values(["datetime", "lat", "lon"]).reset_index(drop=True)
        # fire_name mixes ints (Quebec fires are numbered: 172, 260) with strings
        # (EWF-031-2023); normalise every object column so the cache round-trips.
        for c in d.columns:
            if d[c].dtype == object:
                d[c] = d[c].astype("string")
        d["pcb_id"] = _make_pcb_ids(d)
        return d

    return C.cached("pcb_clean", None, {"v": 2}, [C.XLSX_FP], build)


def _make_pcb_ids(d: pd.DataFrame) -> pd.Series:
    """PCB_YYYY_MMDD_HHMM, with a deterministic a/b/c suffix where a timestamp repeats.

    484 events collapse to only 449 distinct timestamps: pyroCb outbreaks blow up
    several separate fires in the same afternoon (PCB_2023_0729_2220 is 2023WB7,
    2023WB18 and MWF-043-2023, hundreds of km apart).  Ordering by (lat, lon) makes
    the suffix stable across rebuilds.
    """
    base = "PCB_" + d["datetime"].dt.strftime("%Y_%m%d_%H%M")
    out = base.copy()
    for key, idx in base.groupby(base).groups.items():
        if len(idx) == 1:
            continue
        order = d.loc[idx].sort_values(["lat", "lon"]).index
        for i, ix in enumerate(order):
            out.loc[ix] = f"{key}{chr(ord('a') + i)}"
    return out
