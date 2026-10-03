"""Replace the map's Land Cover of Canada layers with Copernicus C3S land cover.

Why (see lc_compare/): the Canada product is a 2015/2020 snapshot, so for fires that burned BEFORE its
epoch it maps the burn scar (31 % "grassland" under fire pixels vs 0.8 % in C3S), it stops at the
border (no Alaska / Russia), and C3S is already what pixels.parquet's lc_class / lc_group use. C3S
is annual, so each fire gets the map of its own year (2023 uses 2022, the last published).

What it rewrites in fire_pcb_map.html's DATA (via map_src/data.json):

  lcpng   per fire (every pyroCb-linked fire, any country): paletted PNG overlay on the 375 m grid
  lc      the old 10 km polygon backdrop is emptied (it only existed for Canada)
  pcb     lc_top{1..5}_cls / _frac / lc_n_px per pyroCb: the class mix under the active-burn cells
          in the 6 h before onset (first 6 h after for blow-ups) - same footprint definition as before

Class indices are the pixels.parquet lc_group order + 1 (0 = transparent / no data); 11 = water:
  1 needleleaf evergreen  2 needleleaf deciduous (larch)  3 broadleaf  4 mixed  5 shrub  6 grass/herb
  7 sparse/lichen  8 wetland  9 cropland  10 non-fuel (urban, bare, ice)  11 water

Run:  python build_c3s_landcover.py [year ...]   then   python map_tool.py join
"""
import base64, gzip, io, json, math
import numpy as np
import rasterio
from PIL import Image

C3S = "C:/Users/jelle/Desktop/PCB_Analysis_1/INPUT_DATA/LANDCOVER_C3S/lccs_class_{y}.tif"
DATA = "map_src/data.json"
R = 6378137.0
CELL_M, CAP_PX, COLLAR_M = 375.0, 640, 5000.0

CLS = {1: (70, 71, 72), 2: (80, 81, 82), 3: (50, 60, 61, 62), 4: (90, 100), 5: (110, 120, 121, 122),
       6: (130,), 7: (140, 150, 151, 152, 153), 8: (160, 170, 180), 9: (10, 11, 12, 20, 30, 40),
       10: (190, 200, 201, 202, 220), 11: (210,)}
LUT = np.zeros(256, np.uint8)
for k, codes in CLS.items():
    for c in codes:
        LUT[c] = k
# keep in step with LC_CLASSES in the map's JS
PALETTE = {1: "#1b5e20", 2: "#9ccc3c", 3: "#f2b705", 4: "#6aa84f", 5: "#c27c3a", 6: "#d9c07a",
           7: "#bdbdbd", 8: "#3a7ca5", 9: "#d95f9e", 10: "#4b4b4b", 11: "#9ecae1"}


def unpack(b64):
    return json.loads(gzip.decompress(base64.b64decode(b64)).decode())


def pack(obj):
    return base64.b64encode(gzip.compress(json.dumps(obj, separators=(",", ":")).encode(), 9)).decode()


def merc_to_ll(x, y):
    return np.degrees(2 * np.arctan(np.exp(np.asarray(y) / R)) - np.pi / 2), np.degrees(np.asarray(x) / R)


class C3SYear:
    def __init__(self, year):
        self.r = rasterio.open(C3S.format(y=min(year, 2022)))
        self.res = self.r.res[0]
        self.west, self.north = self.r.bounds.left, self.r.bounds.top

    def classes(self, lat, lon):
        """LUT class index for arrays of lat/lon (0 where outside the raster / no data)."""
        lat, lon = np.asarray(lat, float), np.asarray(lon, float)
        col = np.floor((lon - self.west) / self.res).astype(int)
        row = np.floor((self.north - lat) / self.res).astype(int)
        out = np.zeros(lat.shape, np.uint8)
        ok = (row >= 0) & (row < self.r.height) & (col >= 0) & (col < self.r.width)
        if not ok.any():
            return out
        r0, r1, c0, c1 = row[ok].min(), row[ok].max() + 1, col[ok].min(), col[ok].max() + 1
        win = self.r.read(1, window=((r0, r1), (c0, c1)))
        out[ok] = LUT[win[row[ok] - r0, col[ok] - c0]]
        return out


def grid_cells(g):
    """(x, y) Web Mercator centres + minutes-relative-t (hours) of a payload grid."""
    i, j, t = (np.cumsum(np.asarray(g[k], np.int64)) for k in ("di", "dj", "dt"))
    return g["ox"] + (i + .5) * g["c"], g["oy"] + (j + .5) * g["c"], t / 60.0


def overlay(c3, g):
    x, y, _ = grid_cells(g)
    half = g["c"] / 2 + COLLAR_M / max(0.15, math.cos(math.radians(merc_to_ll(0, y.mean())[0])))
    x0, x1, y0, y1 = x.min() - half, x.max() + half, y.min() - half, y.max() + half
    cell = g["c"]
    w, h = int(math.ceil((x1 - x0) / cell)), int(math.ceil((y1 - y0) / cell))
    s = min(1.0, CAP_PX / max(w, h))
    w, h = max(1, int(w * s)), max(1, int(h * s))
    xs = x0 + (np.arange(w) + .5) * (x1 - x0) / w
    ys = y1 - (np.arange(h) + .5) * (y1 - y0) / h
    lat, lon = merc_to_ll(*np.meshgrid(xs, ys))
    idx = c3.classes(lat, lon)
    if not idx.any():
        return None
    im = Image.fromarray(idx, "P")
    pal = [0, 0, 0] * 256
    for k, hexcol in PALETTE.items():
        pal[k * 3:k * 3 + 3] = [int(hexcol[i:i + 2], 16) for i in (1, 3, 5)]
    im.putpalette(pal)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True, transparency=0)
    (s_, w_), (n_, e_) = merc_to_ll(x0, y0), merc_to_ll(x1, y1)
    return {"d": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(),
            "b": [[float(s_[()]), float(w_[()])], [float(n_[()]), float(e_[()])]]}


def top_classes(c3, g, blowup):
    x, y, t = grid_cells(g)
    sel = (t >= -6) & (t <= 0)
    if not sel.any():
        sel = (t >= 0) & (t <= 6)
    if not sel.any():
        sel = t <= t.min() + 6
    lat, lon = merc_to_ll(x[sel], y[sel])
    idx = c3.classes(lat, lon)
    idx = idx[idx > 0]
    if not len(idx):
        return {}
    cnt = np.bincount(idx, minlength=12)
    order = np.argsort(-cnt)[:5]
    out = {"lc_n_px": float(len(idx))}
    for n, k in enumerate(order, 1):
        if cnt[k]:
            out[f"lc_top{n}_cls"] = float(k)
            out[f"lc_top{n}_frac"] = round(float(cnt[k]) / len(idx), 3)
    return out


def main(years=None):
    data = json.load(open(DATA, encoding="utf-8"))
    for y, d in sorted(data.items()):
        if years and int(y) not in years:
            continue
        c3 = C3SYear(int(y))
        sim = unpack(d["sim"])
        lcpng = {}
        for mid, v in sim.items():
            o = overlay(c3, v["g"])
            if o:
                lcpng[mid] = o
        pix = unpack(d["pix"])
        pcb = unpack(d["pcb"])
        n_top = 0
        for f in pcb["features"]:
            p = f["properties"]
            for i in range(1, 6):
                p[f"lc_top{i}_cls"] = p[f"lc_top{i}_frac"] = None
            p["lc_n_px"] = None
            g = pix.get(p.get("pcb_id"))
            if g:
                tc = top_classes(c3, g, p.get("active_is_blowup"))
                p.update(tc)
                n_top += bool(tc)
        d["lcpng"] = lcpng
        d["lc"] = pack({"type": "FeatureCollection", "features": []})
        d["pcb"] = pack(pcb)
        print(f"{y}: {len(lcpng)} overlays ({sum(len(o['d']) for o in lcpng.values())/1e3:.0f} kB), "
              f"{n_top}/{len(pcb['features'])} events with a class mix", flush=True)
    json.dump(data, open(DATA, "w", encoding="utf-8", newline="\n"), separators=(",", ":"))


if __name__ == "__main__":
    import sys
    main([int(a) for a in sys.argv[1:]] or None)
