"""Compare Land Cover of Canada (30 m, NRCan) with Copernicus C3S land cover (300 m) at VIIRS fire pixels.

Samples fire pixels from pixels.parquet that fall in Canada, takes the Canada majority class over a
~330 m window (so both products describe the same ground area) and the C3S class of the fire year
(2023 uses 2022, the last C3S year available), maps both onto shared broad fuel classes and writes
agreement statistics to lc_compare/results.json + confusion matrix CSV.
"""
import json, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq, rasterio
from pyproj import Transformer
from rasterio.windows import Window

ROOT = "C:/Users/jelle/Desktop/PCB_Analysis_1/PCB_exploratory"
DATA = "C:/Users/jelle/Desktop/PCB_Analysis_1/INPUT_DATA"
C3S = "C:/Users/jelle/Desktop/PCB_Analysis_1/INPUT_DATA/LANDCOVER_C3S/lccs_class_{y}.tif"
CAN = {2015: DATA + "/LANDCOVER/landcover-2015-classification.tif", 2020: DATA + "/LANDCOVER/landcover-2020-classification.tif"}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
CLS = ["Needleleaf", "Broadleaf deciduous", "Mixed", "Shrub", "Grass/herb", "Lichen/sparse", "Wetland", "Cropland", "Non-fuel"]
CAN_MAP = {1:0, 2:0, 5:1, 6:2, 8:3, 10:4, 11:5, 12:5, 13:5, 14:6, 15:7, 16:8, 17:8, 18:8, 19:8}
C3S_MAP = {10:7,11:7,12:7,20:7,30:7,40:7, 50:1, 60:1,61:1,62:1, 70:0,71:0,72:0, 80:0,81:0,82:0, 90:2, 100:2,
           110:3,120:3,121:3,122:3, 130:4, 140:5,150:5,151:5,152:5,153:5, 160:6,170:6,180:6,
           190:8,200:8,201:8,202:8,210:8,220:8}

rng = np.random.default_rng(0)
px = pq.read_table("C:/Users/jelle/Desktop/PCB_Analysis_1/redundant/old_map_project_pixel_table/pixels.parquet", columns=["lon", "lat", "year"]).to_pandas()
px = px[(px.lon > -141) & (px.lon < -52) & (px.lat > 41)]
print("fire pixels in Canada bbox:", len(px))
px = px.sample(min(N * 2, len(px)), random_state=0)
tf = Transformer.from_crs(4326, 3979, always_xy=True)
rows = []
for yr_group, tif in ((2015, CAN[2015]), (2020, CAN[2020])):
    sub = px[(px.year < 2018) if yr_group == 2015 else (px.year >= 2018)]
    with rasterio.open(tif) as r:
        xs, ys = tf.transform(sub.lon.values, sub.lat.values)
        for (lon, lat, y), x, yv in zip(sub[["lon", "lat", "year"]].itertuples(index=False), xs, ys):
            row, col = r.index(x, yv)
            if row < 5 or col < 5 or row > r.height - 6 or col > r.width - 6:
                continue
            w = r.read(1, window=Window(col - 5, row - 5, 11, 11))
            v = w[w > 0]
            if not len(v):
                continue
            cls = np.bincount(v).argmax()
            rows.append((lon, lat, y, int(cls), int(w[5, 5])))
d = pd.DataFrame(rows, columns=["lon", "lat", "year", "can_maj", "can_pt"]).head(N)
print("sampled with Canada coverage:", len(d))
c3 = np.zeros(len(d), np.uint8)
for y, idx in d.groupby("year").groups.items():
    with rasterio.open(C3S.format(y=min(int(y), 2022))) as r:
        pos = d.loc[idx]
        c3[d.index.get_indexer(idx)] = [v[0] for v in r.sample(list(zip(pos.lon, pos.lat)))]
d["c3s"] = c3
# C3S one year BEFORE the fire (pre-fire fuel, before the burn can show in the annual map)
c3p = np.zeros(len(d), np.uint8)
for y, idx in d.groupby("year").groups.items():
    with rasterio.open(C3S.format(y=max(2012, min(int(y) - 1, 2022)))) as r:
        pos = d.loc[idx]
        c3p[d.index.get_indexer(idx)] = [v[0] for v in r.sample(list(zip(pos.lon, pos.lat)))]
d["c3s_prev"] = c3p
d["can_g"] = d.can_maj.map(CAN_MAP); d["c3s_g"] = d.c3s.map(C3S_MAP); d["c3p_g"] = d.c3s_prev.map(C3S_MAP)
d = d.dropna(subset=["can_g", "c3s_g", "c3p_g"]).astype({"can_g": int, "c3s_g": int, "c3p_g": int})
cm = pd.crosstab(d.can_g, d.c3s_g).reindex(index=range(9), columns=range(9), fill_value=0)
cm.index = [f"CAN {c}" for c in CLS]; cm.columns = [f"C3S {c}" for c in CLS]
cm.to_csv(ROOT + "/lc_compare/confusion_canada_vs_c3s.csv")
agree = float((d.can_g == d.c3s_g).mean())
# forest vs non-forest and needleleaf/broad/mixed (forest type) agreement
forest = lambda s: s.isin([0, 1, 2])
res = {"n": int(len(d)), "overall_agreement": agree,
       "forest_vs_nonforest_agreement": float((forest(d.can_g) == forest(d.c3s_g)).mean()),
       "per_class": {}}
for i, c in enumerate(CLS):
    nc, n3 = int((d.can_g == i).sum()), int((d.c3s_g == i).sum())
    both = int(((d.can_g == i) & (d.c3s_g == i)).sum())
    res["per_class"][c] = {"canada_n": nc, "c3s_n": n3, "both": both,
        "share_canada": nc / len(d), "share_c3s": n3 / len(d),
        "c3s_agrees_when_canada_is_class": both / nc if nc else None}
# Canada map epoch: 2015 map for fires 2013-2017, 2020 map for 2018-2023. "fire before epoch" = the
# map was made AFTER the fire, so it shows the burn scar, not the fuel that burned.
ep = np.where(d.year < 2018, 2015, 2020)
before = d.year <= ep - 0  # fire year <= epoch year
res["agreement_canada_vs_c3s_same_year"] = {"fire_not_after_map_epoch": float((d.can_g == d.c3s_g)[d.year <= ep].mean()), "fire_after_map_epoch": float((d.can_g == d.c3s_g)[d.year > ep].mean())}
res["agreement_canada_vs_c3s_prefire_year"] = float((d.can_g == d.c3p_g).mean())
res["agreement_c3s_same_year_vs_prefire_year"] = float((d.c3s_g == d.c3p_g).mean())
res["by_year_agreement"] = {int(y): float((g.can_g == g.c3s_g).mean()) for y, g in d.groupby("year")}
json.dump(res, open(ROOT + "/lc_compare/results.json", "w"), indent=1)
d.to_csv(ROOT + "/lc_compare/sample.csv", index=False)
print(json.dumps({k: v for k, v in res.items() if k != "per_class"}, indent=1))
print((d.groupby("year").apply(lambda g: pd.Series({"n":len(g),"can_vs_c3s":(g.can_g==g.c3s_g).mean(),"can_vs_c3s_prefire":(g.can_g==g.c3p_g).mean(),"c3s_vs_prefire":(g.c3s_g==g.c3p_g).mean()}))).round(3).to_string())
print(cm.to_string())
for c, v in res["per_class"].items(): print(f"{c:20s} CAN {v['share_canada']:.3f}  C3S {v['share_c3s']:.3f}  agree|CAN {v['c3s_agrees_when_canada_is_class']}")
