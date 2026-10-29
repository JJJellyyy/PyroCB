"""Data for the 3D fire viewer (fire3d.html): the N largest fires of an atlas release, each as a terrain block with its
burn history.

Per fire, fire3d/<year>_<mergid>/ holds (binary files gzipped, little-endian, decompressed in the browser):
  meta.json   ids, area, first detection (t0, UTC), duration, grid (nx, ny, res_m, origin in EPSG:3571), pyroCbs, counts
  dem.bin.gz   int16 elevation (m) on the grid, Copernicus GLO-90 DSM averaged into the grid cells; cells without a local
              tile (outside our download) take the nearest value (share in meta dem_filled_frac)
  burn.bin.gz  uint16 smoothed first-burn time per grid cell: 0 = not burned, else 1 + round(4 * hours since t0);
              inverse-distance mean (power 2) of the 8 nearest first-burn points, only inside the fused perimeter
              (VIIRS final perimeter U MODIS-added cells, pcbcube.fusion). DISPLAY interpolation, not data.
  src.bin.gz   uint8 per grid cell: source of the nearest first-burn point (1 VIIRS, 2 MODIS-added) - 0 outside
  pix.bin.gz   per point float32 x, y (m from the grid's top-left corner, x east, y south), float32 t (h since t0),
              float32 value (VIIRS: FRP density MW/km2 from the atlas NFP; MODIS cells: MCD14ML FRP in MW where a
              Terra/Aqua detection matched, else NaN - reference only), uint8 kind (1 S-NPP, 2 NOAA-20, 3 NOAA-21,
              4 MODIS-added cell, 5 MODIS cell inside the VIIRS perimeter); arrays stored one after the other, with
              float32 t_next (h, next detection within SAME_M of the point, inf = none) before kind.
First-burn points: per 375 m bin the earliest VIIRS detection, plus every MODIS cell joined to the fire with its fused
burn time (pcbcube.fusion t_burn). Grid: bbox of the fused perimeter + 8 %, cell = max(90 m, longest side / MAXDIM).

    python build_fire3d.py [--release release_multi] [--top 10]
"""
from __future__ import annotations

import argparse
import glob
import gzip
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
from pcbcube import config as C  # noqa: E402

OUT = HERE / "fire3d"
DEM_DIR = C.INPUT / "TOPOLOGY" / "cop90"
MAXDIM = 1280
MARGIN = 0.08
K = 8
BIN_M = 375.0
SAME_M = 250.0                                         # 'same spot' for re-detections (~2/3 of a VIIRS I-band pixel)
SAT = {"SNPP": 1, "NOAA20": 2, "NOAA21": 3}


def top_fires(rel, n):
    rows = []
    for fp in sorted(glob.glob(str(rel / "final_perims" / "final_viirs*.gpkg"))):
        G = gpd.read_file(fp, ignore_geometry=True, columns=["mergid", "farea"])
        G["year"] = int(Path(fp).stem[-4:])
        rows.append(G)
    return pd.concat(rows).sort_values("farea", ascending=False).head(n)


def dem_grid(x0, y0, nx, ny, res):
    """Copernicus GLO-90 tiles averaged onto the EPSG:3571 grid (top-left x0, y0)."""
    import rasterio
    from rasterio.merge import merge
    from rasterio.warp import Resampling, reproject, transform_bounds
    lo0, la0, lo1, la1 = transform_bounds(3571, 4326, x0, y0 - ny * res, x0 + nx * res, y0, densify_pts=50)
    tiles = []
    for la in range(int(np.floor(la0)), int(np.floor(la1)) + 1):
        for lo in range(int(np.floor(lo0)), int(np.floor(lo1)) + 1):
            ns, ew = ("N" if la >= 0 else "S"), ("E" if lo >= 0 else "W")
            f = DEM_DIR / f"Copernicus_DSM_COG_30_{ns}{abs(la):02d}_00_{ew}{abs(lo):03d}_00_DEM.tif"
            if f.exists():
                tiles.append(rasterio.open(f))
    out = np.zeros((ny, nx), np.float32)
    if not tiles:
        return out.astype(np.int16), 1.0
    mos, tr = merge(tiles, nodata=0)
    for t in tiles:
        t.close()
    reproject(mos[0], out, src_transform=tr, src_crs="EPSG:4326", src_nodata=0,
              dst_transform=rasterio.Affine(res, 0, x0, 0, -res, y0), dst_crs="EPSG:3571", dst_nodata=0,
              resampling=Resampling.average)
    hole = out == 0                                    # tiles outside our download (e.g. across a border): nearest value, display only
    if hole.any() and (~hole).any():
        from scipy.ndimage import distance_transform_edt
        _, (ri, ci) = distance_transform_edt(hole, return_indices=True)
        out = out[ri, ci]
    return np.round(out).astype(np.int16), float(hole.mean())


def build(year, mergid, rel, release):
    import rasterio.features
    from pyproj import Transformer
    from scipy.spatial import cKDTree
    key = f"{year}_{mergid}"
    odir = OUT / key
    odir.mkdir(parents=True, exist_ok=True)
    fus = ROOT / "OUTPUT_DATA" / "fusion" / release

    N = pd.read_parquet(rel / str(year) / "compact" / "nfp.parquet", columns=["x", "y", "frp", "sat", "datetime_utc", "mergid"])
    N = N[N.mergid == mergid].copy()
    N["dt"] = pd.to_datetime(N.datetime_utc)
    t0 = N.dt.min()
    N["t"] = (N.dt - t0) / pd.Timedelta(hours=1)

    Fz = gpd.read_file(fus / f"fused_{year}.gpkg") if (fus / f"fused_{year}.gpkg").exists() else None
    if Fz is not None and (Fz.mergid == mergid).any():
        per = Fz[Fz.mergid == mergid].geometry.union_all()
        D = pd.read_parquet(fus / f"cells_{year}.parquet", columns=["x", "y", "mergid", "cls", "t_burn", "t_src", "af_frp"])
        Mc = D[(D.mergid == mergid) & D.cls.isin(["modis_add", "modis_in"])].copy()
        Mc["t"] = (pd.Timestamp(year, 1, 1) + pd.to_timedelta(Mc.t_burn - 1, unit="D") - t0) / pd.Timedelta(hours=1)
        fused_id = int(Fz.loc[Fz.mergid == mergid, "fused_id"].iloc[0])
        area_v = float(Fz.loc[Fz.mergid == mergid, "area_viirs"].sum())
        area_m = float(Fz.loc[Fz.mergid == mergid, "area_modis_add"].sum())
    else:                                              # no fusion output: VIIRS final perimeter only
        P = gpd.read_file(rel / "final_perims" / f"final_viirs{year}.gpkg").to_crs(3571)
        per = P[P.mergid == mergid].geometry.union_all()
        Mc = pd.DataFrame(columns=["x", "y", "cls", "t", "af_frp"])
        fused_id = mergid
        area_v, area_m = per.area / 1e6, 0.0

    bx0, by0, bx1, by1 = per.bounds
    w, h = bx1 - bx0, by1 - by0
    bx0, bx1, by0, by1 = bx0 - MARGIN * w, bx1 + MARGIN * w, by0 - MARGIN * h, by1 + MARGIN * h
    res = max(90.0, max(bx1 - bx0, by1 - by0) / MAXDIM)
    nx, ny = int(np.ceil((bx1 - bx0) / res)), int(np.ceil((by1 - by0) / res))
    x0, y0 = bx0, by1                                  # top-left

    dem, dem_filled = dem_grid(x0, y0, nx, ny, res)

    # first-burn points: earliest VIIRS detection per 375 m bin + MODIS cells
    b = pd.DataFrame({"bx": np.floor(N.x / BIN_M), "by": np.floor(N.y / BIN_M), "x": N.x, "y": N.y, "t": N.t})
    fb = b.sort_values("t").drop_duplicates(["bx", "by"])
    px = np.r_[fb.x.to_numpy(), Mc.x.to_numpy(float)]
    py = np.r_[fb.y.to_numpy(), Mc.y.to_numpy(float)]
    pt = np.r_[fb.t.to_numpy(), Mc.t.to_numpy(float)]
    ps = np.r_[np.ones(len(fb), np.uint8), np.where(Mc.cls.to_numpy() == "modis_add", 2, 1).astype(np.uint8)]
    ok = np.isfinite(pt)
    px, py, pt, ps = px[ok], py[ok], pt[ok], ps[ok]

    mask = rasterio.features.rasterize([(per, 1)], out_shape=(ny, nx), transform=rasterio.Affine(res, 0, x0, 0, -res, y0),
                                       fill=0, dtype=np.uint8, all_touched=True).astype(bool)
    jj, ii = np.nonzero(mask)
    cx, cy = x0 + (ii + 0.5) * res, y0 - (jj + 0.5) * res
    tree = cKDTree(np.c_[px, py])
    d, j = tree.query(np.c_[cx, cy], k=min(K, len(px)))
    d, j = d.reshape(len(cx), -1), j.reshape(len(cx), -1)
    wgt = 1.0 / np.maximum(d, 1.0) ** 2
    tt = (wgt * pt[j]).sum(1) / wgt.sum(1)
    burn = np.zeros((ny, nx), np.uint16)
    burn[jj, ii] = (1 + np.round(4 * np.clip(tt, 0, None))).clip(1, 65535).astype(np.uint16)
    src = np.zeros((ny, nx), np.uint8)
    src[jj, ii] = ps[j[:, 0]]

    # points for blocks / point cloud
    kind_v = N.sat.map(SAT).fillna(1).to_numpy(np.uint8)
    kind_m = np.where(Mc.cls.to_numpy() == "modis_add", 4, 5).astype(np.uint8)
    X = np.r_[N.x.to_numpy(), Mc.x.to_numpy(float)] - x0
    Y = y0 - np.r_[N.y.to_numpy(), Mc.y.to_numpy(float)]
    T = np.r_[N.t.to_numpy(), Mc.t.to_numpy(float)]
    V = np.r_[N.frp.to_numpy(float), Mc.af_frp.to_numpy(float)]
    Kd = np.r_[kind_v, kind_m]
    o = np.argsort(T, kind="stable")                   # ascending time: the viewer finds windows by binary search
    X, Y, T, V, Kd = X[o], Y[o], T[o], V[o], Kd[o]
    # time of the next detection at the same spot (any point within SAME_M, strictly later): the viewer's "active front"
    # shows each spot by its LAST detection, hiding a detection once a newer one there is visible. inf = none.
    TN = np.full(len(T), np.inf, np.float32)
    pr = cKDTree(np.c_[X, Y]).query_pairs(SAME_M, output_type="ndarray")
    if len(pr):
        a, b2 = pr[:, 0], pr[:, 1]
        sw = T[a] > T[b2]
        a, b2 = np.where(sw, b2, a), np.where(sw, a, b2)  # a = earlier, b2 = later
        later = T[b2] > T[a] + 0.01                      # same scan: neither replaces the other
        np.minimum.at(TN, a[later], T[b2[later]].astype(np.float32))

    def wr(name, *arrs):
        (odir / name).write_bytes(gzip.compress(b"".join(np.ascontiguousarray(a).tobytes() for a in arrs), 6))
    wr("dem.bin.gz", dem.astype("<i2"))
    wr("burn.bin.gz", burn.astype("<u2"))
    wr("src.bin.gz", src)
    wr("pix.bin.gz", X.astype("<f4"), Y.astype("<f4"), T.astype("<f4"), V.astype("<f4"), TN.astype("<f4"), Kd)

    # pyroCbs inside the fire: within 10 km of the fused perimeter, between first detection - 1 d and last + 1 d
    from pcbcube.pcb import read_pcb
    Pc = read_pcb()
    Pc = Pc[Pc.year == year]
    tf = Transformer.from_crs(4326, 3571, always_xy=True)
    qx, qy = tf.transform(Pc.lon.to_numpy(), Pc.lat.to_numpy())
    near = shapely.dwithin(per, shapely.points(qx, qy), 10000)
    th = ((pd.to_datetime(Pc.datetime) - t0) / pd.Timedelta(hours=1)).to_numpy()
    sel = near & (th >= -24) & (th <= T.max() + 24)
    pcbs = [dict(id=i, h=round(float(hh), 2), x=round(float(a - x0), 1), y=round(float(y0 - c), 1), status=s,
                 inj=None if pd.isna(z) else float(z))
            for i, hh, a, c, s, z in zip(Pc.pcb_id[sel], th[sel], qx[sel], qy[sel], Pc.status[sel], Pc.max_inject_alt[sel])]

    lon, lat = Transformer.from_crs(3571, 4326, always_xy=True).transform(x0 + nx * res / 2, y0 - ny * res / 2)
    farea = float(gpd.read_file(rel / "final_perims" / f"final_viirs{year}.gpkg", ignore_geometry=True,
                                columns=["mergid", "farea"]).set_index("mergid").farea[mergid])
    meta = dict(key=key, year=year, mergid=int(mergid), fused_id=fused_id, release=release, farea_km2=round(farea, 1),
                area_fused_km2=round(per.area / 1e6, 1), area_viirs_km2=round(area_v, 1), area_modis_add_km2=round(area_m, 1), t0=t0.strftime("%Y-%m-%d %H:%M"), dur_h=round(float(T.max()), 2),
                lon=round(lon, 3), lat=round(lat, 3), grid=dict(nx=nx, ny=ny, res_m=round(res, 3), x0=round(x0, 1), y0=round(y0, 1)),
                dem_min=int(dem[dem != 0].min()) if (dem != 0).any() else 0, dem_max=int(dem.max()), dem_filled_frac=round(dem_filled, 4),
                n_viirs=int((Kd <= 3).sum()), n_modis_add=int((Kd == 4).sum()), n_modis_in=int((Kd == 5).sum()),
                n_points=int(len(T)), pcbs=pcbs)
    (odir / "meta.json").write_text(json.dumps(meta))
    size = sum(f.stat().st_size for f in odir.iterdir()) / 1e6
    print(f"{key}: {nx}x{ny} @ {res:.0f} m, {len(T):,} points, {len(pcbs)} pyroCbs, {size:.1f} MB", flush=True)
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", default="release_multi")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--only", nargs="*", default=None, help="year_mergid keys to (re)build")
    a = ap.parse_args()
    rel = C.INPUT / "atlas_extension" / a.release
    T = top_fires(rel, a.top)
    OUT.mkdir(exist_ok=True)
    idx = []
    for y, m in zip(T.year, T.mergid):
        k = f"{y}_{m}"
        if a.only and k not in a.only:
            mp = OUT / k / "meta.json"
            if mp.exists():
                idx.append(json.loads(mp.read_text()))
            continue
        idx.append(build(int(y), int(m), rel, a.release))
    keep = ["key", "year", "mergid", "farea_km2", "area_fused_km2", "t0", "dur_h", "lon", "lat", "n_points"]
    (OUT / "index.json").write_text(json.dumps([{k: d[k] for k in keep} | {"n_pcb": len(d["pcbs"])} for d in idx]))


if __name__ == "__main__":
    main()
