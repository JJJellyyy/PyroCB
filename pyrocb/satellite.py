"""Post-fire Sentinel-2 imagery: a true-colour circle around each merged perimeter
some weeks after its official end, with the perimeter traced on top.

Source is Sentinel-2 L2A's "visual" (true-colour) asset via Microsoft's Planetary
Computer STAC API - free, no account, no API key. This makes live network calls per
fire, unlike the rest of the pipeline; `config.SAT_YEARS` (default: every resolved
year) can still scope it down for a faster rebuild while iterating.

A pre-fire image used to ship alongside this one. Dropped (D-59): comparing "before"
against "after" is the reviewer's job with the perimeter/growth-band layers already on
the map, and every fire getting a post image at every year was worth more than both
phases at one year alone.

**The circle.** Centred on the perimeter's own centroid, radius = the farthest the
perimeter reaches from that centroid (`hausdorff_distance` to a point IS that
distance) times `SAT_RADIUS_MARGIN`, floored at `SAT_RADIUS_MIN_KM` so small fires
still get context. This is "around the whole perimeter" taken literally, which means
radius varies enormously with fire size (2-160 km across 2023's matched fires alone),
which is why output resolution scales with it (`SAT_OUT_PX`, a pixel BUDGET) rather
than staying at a fixed ground resolution - see D-51.

**Search window, and why by DATE not by cloud.** Searches forward from the official
end (`ted_*`), offset by `SAT_POST_WINDOW_D[0]` days so the scene is dead ground, not
still-smoking ground. Candidates are sorted by closeness to that anchor date, not by
cloud cover - a boreal fire ending in October can have its globally-clearest scene in
a 120-day window land in November, by which point the ground is already under snow.
Proximity keeps the scene seasonally honest even at the cost of some cloud (D-52).

**Mosaicking.** A single Sentinel-2 tile is 110 km on a side, so a large fire's circle
can straddle two or three tiles (or two UTM zones). Up to `SAT_MAX_TILES` candidates
are read, closest-date-first, each filling in only the pixels still empty - i.e. the
closest-in-time scene wins every pixel it has, later scenes only patch gaps. Coverage
is reported per image; a real gap (no scene available) stays transparent.

**Delivery: files on disk, not embedded.** Even after circle-cropping, a true-colour
1024 px image is ~1.5-2 MB as PNG. Across every matched fire in the dataset that is
easily 100+ MB - embedding that in `fire_pcb_map.html` the way every other layer works
would multiply the map's size several times over for a feature most sessions never
open. Each image is instead flattened onto a solid background (no georeferencing is
needed - this is a standalone photo, not a map overlay) and written as a JPEG file
under `config.SAT_IMG_DIR`; the map loads it lazily via a plain `<img src>` only when
a viewer opens that fire's image panel. See D-53.
"""
from __future__ import annotations

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import planetary_computer as pc
import pystac_client
import rasterio
import shapely
from PIL import Image, ImageDraw
from rasterio.enums import Resampling
from rasterio.warp import transform, transform_bounds
from rasterio.windows import from_bounds

from . import config as C
from . import io as pio

_STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
_COLLECTION = "sentinel-2-l2a"
_BG = (17, 17, 17)   # flattened background outside the circle / where no scene reached


def fire_images(year: int, mergids, out_dir: str, verbose: bool = True) -> dict:
    """{str(mergid): meta} for every mergid that got a post-fire image.

    meta = {"file": path relative to config.ROOT, "date": iso, "cloud": pct,
            "coverage": 0..1, "n_scenes": int}
    Writes one JPEG per fire into `out_dir` as a side effect.
    """
    want = {int(m) for m in mergids}
    fin = pio.read_final(year, metric=True)
    fin = fin[fin["mergid"].isin(want)]
    if not len(fin):
        return {}
    os.makedirs(out_dir, exist_ok=True)

    cat = pystac_client.Client.open(_STAC_URL, modifier=pc.sign_inplace)
    out: dict[str, dict] = {}
    for i, r in enumerate(fin.itertuples(), 1):
        geom = r.geometry
        cen = shapely.centroid(geom)
        radius_m = max(shapely.hausdorff_distance(shapely.points(cen.x, cen.y), geom)
                       * C.SAT_RADIUS_MARGIN, C.SAT_RADIUS_MIN_KM * 1000.0)
        bbox = gpd.GeoSeries([shapely.buffer(cen, radius_m)], crs=C.METRIC_CRS) \
                  .to_crs(C.WGS84).total_bounds

        t1 = pd.Timestamp(int(r.ted_year), int(r.ted_month), int(r.ted_day))
        post_dt = (f"{(t1 + pd.Timedelta(days=C.SAT_POST_WINDOW_D[0])).date()}/"
                   f"{(t1 + pd.Timedelta(days=C.SAT_POST_WINDOW_D[1])).date()}")

        mid = int(r.mergid)
        post = _build_image(cat, geom, cen, radius_m, bbox, post_dt, t1,
                            os.path.join(out_dir, f"{mid}_post.jpg"))
        if post:
            out[str(mid)] = post
        if verbose:
            tag = "-" if post is None else f"{post['coverage']*100:.0f}% cov, {post['cloud']:.1f}% cloud"
            print(f"    sat {i}/{len(fin)} mergid {mid}: {tag}", flush=True)
    if verbose:
        print(f"  {year}: satellite imagery for {len(out)}/{len(fin)} fires", flush=True)
    return out


def _search(cat, bbox, datetime_range: str, anchor: pd.Timestamp):
    """Candidates under the cloud ceiling, sorted by closeness to `anchor` - NOT by
    cloud cover. A scene chosen purely for minimum cloud can drift to the far end of
    a wide search window and land in a different season (see module docstring).
    Deduped to one item per (date, tile), keeping the lower-cloud copy: the catalogue
    carries both an original and a reprocessed copy of many scenes.
    """
    items = cat.search(collections=[_COLLECTION], bbox=list(bbox),
                       datetime=datetime_range,
                       query={"eo:cloud_cover": {"lt": C.SAT_MAX_CLOUD}},
                       limit=100).items()
    best: dict[tuple, object] = {}
    for it in items:
        k = (it.datetime.date(), it.properties.get("s2:mgrs_tile"))
        cur = best.get(k)
        if cur is None or it.properties.get("eo:cloud_cover", 100) < cur.properties.get("eo:cloud_cover", 100):
            best[k] = it
    anchor_utc = anchor.tz_localize("UTC") if anchor.tzinfo is None else anchor
    return sorted(best.values(), key=lambda it: abs((it.datetime - anchor_utc).total_seconds()))


def _build_image(cat, geom, cen, radius_m: float, bbox_wgs84, datetime_range: str,
                 anchor: pd.Timestamp, out_path: str):
    """One mosaic + perimeter overlay, flattened and written as a JPEG file."""
    candidates = _search(cat, bbox_wgs84, datetime_range, anchor)
    if not candidates:
        return None

    out_px = C.SAT_OUT_PX
    out_crs = None
    win_transform = None
    bx = None
    canvas = np.zeros((3, out_px, out_px), dtype=np.uint8)
    filled = np.zeros((out_px, out_px), dtype=bool)
    used = []

    for it in candidates[:C.SAT_MAX_TILES]:
        if filled.all():
            break
        try:
            href = it.assets["visual"].href
            with rasterio.open(href) as src:
                if out_crs is None:
                    out_crs = src.crs
                    bx = transform_bounds("EPSG:4326", out_crs, *bbox_wgs84)
                    win_transform = rasterio.transform.from_bounds(*bx, out_px, out_px)
                b2 = transform_bounds(out_crs, src.crs, *bx)
                win = from_bounds(*b2, transform=src.transform)
                data = src.read([1, 2, 3], window=win, out_shape=(3, out_px, out_px),
                                resampling=Resampling.average, boundless=True, fill_value=0)
        except Exception:
            continue
        valid = data.any(axis=0) & ~filled
        if not valid.any():
            continue
        for b in range(3):
            canvas[b][valid] = data[b][valid]
        filled |= valid
        used.append(it)

    if not used:
        return None

    coverage = float(filled.mean())
    px_size = (bx[2] - bx[0]) / out_px
    cx, cy = transform("EPSG:3571", out_crs, [cen.x], [cen.y])
    cen_px = ~win_transform * (cx[0], cy[0])
    yy, xx = np.mgrid[0:out_px, 0:out_px]
    on_ground = (np.hypot(xx - cen_px[0], yy - cen_px[1]) * px_size <= radius_m) & filled

    img = Image.fromarray(np.moveaxis(canvas, 0, -1), "RGB")
    _draw_perimeter(img, geom, out_crs, win_transform)
    # Flatten to a solid background where the mosaic never reached (circle edge,
    # unfilled tile gaps) - a standalone photo has no map to blend transparently with.
    bg = Image.new("RGB", img.size, _BG)
    bg.paste(img, mask=Image.fromarray((on_ground * 255).astype(np.uint8)))
    bg.save(out_path, "JPEG", quality=C.SAT_JPEG_QUALITY, optimize=True)

    best = used[0]
    return {
        "file": os.path.relpath(out_path, C.ROOT).replace(os.sep, "/"),
        "date": best.datetime.date().isoformat(),
        "cloud": round(float(best.properties.get("eo:cloud_cover", 0.0)), 1),
        "coverage": round(coverage, 3), "n_scenes": len(used),
    }


def _draw_perimeter(img: Image.Image, geom, out_crs, win_transform,
                    dash: int = 10, gap: int = 6, width: int = 3,
                    color=(255, 40, 40)) -> None:
    """Burn the fire perimeter into `img` as a red dashed line, in place."""
    per = gpd.GeoSeries([geom], crs=C.METRIC_CRS).to_crs(out_crs).iloc[0]
    parts = shapely.get_parts(per) if per.geom_type == "MultiPolygon" else [per]
    draw = ImageDraw.Draw(img)
    for part in parts:
        pts = [(~win_transform * (x, y)) for x, y in part.exterior.coords]
        for i in range(len(pts) - 1):
            x0, y0 = pts[i]
            x1, y1 = pts[i + 1]
            seg = float(np.hypot(x1 - x0, y1 - y0))
            if seg < 1e-6:
                continue
            n = max(1, int(seg / (dash + gap)))
            for k in range(n + 1):
                t0 = k * (dash + gap) / seg
                if t0 >= 1.0:
                    break
                t1 = min(1.0, t0 + dash / seg)
                draw.line([(x0 + (x1 - x0) * t0, y0 + (y1 - y0) * t0),
                          (x0 + (x1 - x0) * t1, y0 + (y1 - y0) * t1)],
                         fill=color, width=width)
