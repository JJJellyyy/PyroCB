"""PyroCb cloud imagery: a GOES visible-band snapshot of each event's storm, with the
reported point and every v2 matched-fire position marked on it for comparison.

The post-fire imagery (`satellite.py`) answers "did this really burn here"; this
answers a different question - "what did the actual pyroCb look like" - and needs a
different kind of satellite entirely. Sentinel-2 revisits any point roughly every 5
days, so the odds of a pass landing near one specific pyroCb are close to zero. A
geostationary satellite images the same hemisphere every 5-15 minutes, so it can
actually catch the event, not just the neighbourhood it happened in.

**Source: GOES-16/17/18/19 (`goes-cmi` on Planetary Computer), visible red band C02.**
Two structural coverage gaps, both permanent, not implementation gaps:

* **Region.** GOES only covers the Western Hemisphere. Planetary Computer carries no
  geostationary satellite over Siberia, so `region == "Northern Asia"` events - about
  a quarter of the atlas - get no image, ever, from this source.
* **Time.** GOES-16 was not operational until December 2017, so 2013-2017 (`GOES_MIN_YEAR`)
  is out regardless of location.

**Band choice: visible, not infrared, and not a big deal that it's day-only.** The
obvious "works day and night" choice is the 10.3 um clean-IR band (brightness
temperature - literally how meteorologists spot a pyroCb's cold overshooting top). It
was tried first and measured against the alternative: IR's 2 km pixels give a blurry
native window at the radii this ships at, while the single visible band (C02, 0.5 km)
reads as actual cloud texture - no colour synthesis needed, one band is enough.

**Search: earliest daylight scene 2-12 h after onset, not the instant of onset.** A
pyroCb often has not yet built into a distinct, recognisable cloud at the exact minute
the atlas timestamps as onset, and a meaningful fraction of onsets land at local
dusk/night regardless. `GOES_POST_WINDOW_H = (2, 12)` gives the storm time to mature;
`GOES_DAYLIGHT_LOCAL_H` skips candidates before they are ever downloaded (cheap - no
fetch); `GOES_MIN_DAYLIGHT_REFL` is a post-read floor for the dawn/dusk edge that
coarse hour band lets through. Earliest candidate that clears both wins.

**Markers (D-69).** Centred on the AS-REPORTED point (not a v2 position - the whole
point is comparing v2 candidates against a fixed reference), at `GOES_RADIUS_KM` = 60,
wide enough that every (radius, method) combination's matched position - up to 20 km
from the report - sits comfortably inside with room around it. Every valid
`fire2_r{R}_{M}_lat/lon` gets a small dot: colour by radius (`RADIUS_COLOR` - 10 km
purple, 15 km blue, 20 km yellow), number by method (`METHOD_NUMBER` - 1 closest,
2 before, 3 before+after, 4 before2+after2), so a reviewer can see at a glance which
combinations land where relative to the actual convective cloud. The reported point
itself keeps the plain red crosshair used elsewhere on the map. A small fixed legend
is burned into a corner of every image (same position every time, so it reads as part
of the frame rather than clutter).

**v3/v4 markers (D-70, D-78).** `PCB_MATCHING_CLEAN.ipynb`'s two single-position
matches - the radius-ladder match (v3, green star) and the 4-overpass time-weighted
match (v4, purple star) - each get their own star, drawn last and in that order (v4
on top) so they are always the most visible marks on the image, never hidden by an
agreeing v2 dot or by each other (the declutter rosette still separates them if they
coincide). Missing/unmatched events simply carry no star for that version, same as a
missing v2 dot.

**No mosaicking needed to zoom out.** Unlike Sentinel-2 (110 km UTM tiles, D-51),
GOES full disk is ONE contiguous raster covering the whole hemisphere - "zoom out"
here is just a bigger read window on the same scene, not a multi-tile stitch. Tested
to well past `GOES_RADIUS_KM`; the only failure mode already handled is a genuine
edge-of-visible-disk point (D-69's projection-domain guard), which no amount of
compositing would fix - that point is outside what the satellite can see, from any
scene of it.

**Delivery**: same as `satellite.py` - a JPEG file per event under `config.GOES_IMG_DIR`,
metadata only in the map payload, loaded lazily via `<img src>`.
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
from PIL import Image, ImageDraw, ImageFont
from rasterio.enums import Resampling
from rasterio.warp import transform, transform_bounds
from rasterio.windows import from_bounds

from . import config as C
from . import nearby as N
from . import timebins as T

_STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
_COLLECTION = "goes-cmi"
_BAND = "C02_0.5km"
_BG = (8, 8, 16)

# D-69: colour = radius, number = method - fixed across every image so the scheme
# only has to be learned once.
RADIUS_COLOR = {10.0: (168, 85, 247), 15.0: (59, 130, 246), 20.0: (250, 204, 21)}   # purple, blue, yellow
METHOD_NUMBER = {"closest": 1, "before": 2, "before_after": 3, "before2_after2": 4}
_REPORTED_COLOR = (255, 40, 40)
# v3/v4 (D-70, D-78): PCB_MATCHING_CLEAN.ipynb's single-position matches - not a
# (radius, method) grid like v2, so each gets its own shape (star) rather than
# competing for a dot colour. Drawn last, on top of everything including the
# reported crosshair, so neither is ever hidden by an agreeing v2 marker.
_V3_COLOR = (34, 197, 94)     # green  - v3, the radius-ladder match
_V4_COLOR = (168, 85, 247)    # purple - v4, the 4-overpass time-weighted match
_STAR_COLORS = {"v3": _V3_COLOR, "v4": _V4_COLOR}


def event_images(year: int, events: pd.DataFrame, out_dir: str, verbose: bool = True) -> dict:
    """{pcb_id: meta} for every event GOES can plausibly have seen.

    `events` needs pcb_id, datetime, lat, lon, region, and every `fire2_r{R}_{M}_lat
    /lon` column `nearby.summarise` produces (missing ones are just skipped as
    markers), plus optional `v3_lat`/`v3_lon` and `v4_lat`/`v4_lon` (PCB_MATCHING_CLEAN's
    two single-position matches - missing/NaN just skips that star for that event,
    same as a missing fire2_* skips its dot).
    meta = {"file": path relative to config.ROOT, "date": iso scene timestamp,
            "dt_h": hours after onset the scene was taken, "coverage": 0..1,
            "n_markers": how many (radius, method) positions got drawn}.
    """
    if year < C.GOES_MIN_YEAR:
        return {}
    sub = events[~events["region"].isin(C.GOES_EXCLUDE_REGIONS)].copy()
    if not len(sub):
        return {}
    os.makedirs(out_dir, exist_ok=True)

    cat = pystac_client.Client.open(_STAC_URL, modifier=pc.sign_inplace)
    out: dict[str, dict] = {}
    for i, r in enumerate(sub.itertuples(), 1):
        bbox = _bbox_wgs84(r.lat, r.lon, C.GOES_RADIUS_KM)
        markers, stars = _collect_markers(r)
        meta = _build_image(cat, r.lat, r.lon, markers, stars, bbox, float(r.lon),
                            pd.Timestamp(r.datetime), os.path.join(out_dir, f"{r.pcb_id}_cloud.jpg"))
        if meta:
            out[r.pcb_id] = meta
        if verbose:
            tag = "-" if meta is None else f"+{meta['dt_h']:.1f} h, {meta['n_markers']} markers, {meta['coverage']*100:.0f}% cov"
            print(f"    goes {i}/{len(sub)} {r.pcb_id}: {tag}", flush=True)
    if verbose:
        print(f"  {year}: pyroCb cloud imagery for {len(out)}/{len(sub)} eligible events", flush=True)
    return out


def _bbox_wgs84(lat: float, lon: float, radius_km: float):
    """WGS84 bounds of a circle, computed via the equal-area METRIC_CRS so `radius_km`
    is a true ground distance (a plain lat/lon box would distort badly this far north)."""
    cen = gpd.GeoSeries([shapely.points(lon, lat)], crs=C.WGS84).to_crs(C.METRIC_CRS).iloc[0]
    return gpd.GeoSeries([shapely.buffer(cen, radius_km * 1000.0)], crs=C.METRIC_CRS) \
              .to_crs(C.WGS84).total_bounds


def _collect_markers(r) -> tuple[list[dict], dict]:
    """(v2 markers, star points) for this event row.

    v2 markers: every valid (radius, method) position, as {"lat", "lon",
    "radius_km", "method"}. star points: {"v3": {"lat","lon"} | None, "v4": ...},
    from the row's `v3_lat`/`v3_lon` and `v4_lat`/`v4_lon` columns - None (never
    drawn) if either is missing or NaN."""
    out = []
    for radius in C.NEARBY_RADII_KM:
        for method in C.NEARBY_METHODS:
            pre = f"fire2_{N._key(radius, method)}"
            lat = getattr(r, f"{pre}_lat", None)
            lon = getattr(r, f"{pre}_lon", None)
            if lat is not None and lon is not None and pd.notna(lat) and pd.notna(lon):
                out.append({"lat": float(lat), "lon": float(lon),
                           "radius_km": radius, "method": method})

    stars = {}
    for v in ("v3", "v4"):
        vlat = getattr(r, f"{v}_lat", None)
        vlon = getattr(r, f"{v}_lon", None)
        stars[v] = ({"lat": float(vlat), "lon": float(vlon)}
                    if vlat is not None and vlon is not None and pd.notna(vlat) and pd.notna(vlon)
                    else None)
    return out, stars


def _candidate_scenes(cat, bbox, t_onset: pd.Timestamp, lon: float):
    """FULL DISK scenes in [onset+GOES_POST_WINDOW_H], EARLIEST first, restricted to
    scenes whose local solar hour at `lon` falls inside GOES_DAYLIGHT_LOCAL_H - a
    cheap check with no fetch, so obviously-dark candidates never cost a read."""
    win_lo = t_onset + pd.Timedelta(hours=C.GOES_POST_WINDOW_H[0])
    win_hi = t_onset + pd.Timedelta(hours=C.GOES_POST_WINDOW_H[1])
    items = list(cat.search(collections=[_COLLECTION], bbox=list(bbox),
                            datetime=f"{win_lo.isoformat()}Z/{win_hi.isoformat()}Z",
                            query={"goes:image-type": {"eq": "FULL DISK"}}).items())
    items.sort(key=lambda it: it.datetime)
    lo_h, hi_h = C.GOES_DAYLIGHT_LOCAL_H
    out = []
    for it in items:
        # STAC datetimes are always UTC per spec, but .tz_convert first (a no-op if
        # already UTC) rather than a bare .tz_localize(None), which would silently
        # keep whatever OTHER timezone's wall-clock a non-conformant catalogue
        # returned instead of actually converting to UTC.
        scene_utc = pd.Timestamp(it.datetime).tz_convert("UTC").tz_localize(None)
        lh = float(T.local_solar_hour(scene_utc, lon))
        if lo_h <= lh <= hi_h:
            out.append(it)
    return out


def _build_image(cat, cen_lat: float, cen_lon: float, markers: list[dict], stars: dict,
                 bbox_wgs84, lon: float, t_onset: pd.Timestamp, out_path: str):
    """Earliest daylight scene 2-12 h after onset, circle-cropped, marked up and
    written as a JPEG. `cen_lat/lon` is the AS-REPORTED point (D-69) - the fixed
    reference every marker is compared against, not a v2 result. `stars`
    (D-70, D-78) is {"v3": {"lat","lon"}|None, "v4": ...} - each present entry is
    one single-position match drawn as its own star."""
    out_px = C.GOES_OUT_PX
    star_keys = [k for k in ("v3", "v4") if stars.get(k) is not None]
    it = bx = win_transform = cen_px = star_px = marker_px = src_meta = data = None

    # A candidate's bbox, or any one of the marker points, can fail transform for a
    # location beyond that satellite's actual visible geostationary disk, even when
    # the STAC search's coarse bounding-box match let the item through. Everything
    # that touches this candidate's projection - bbox, read, AND every marker
    # (batched into one transform call) - has to happen inside this one try block:
    # a candidate that passes the bbox/read can still fail on a marker, and falling
    # through to the next candidate only works if nothing has been committed yet.
    # cen (index 0), then each present star in star_keys order, then every v2 marker.
    n_lead = 1 + len(star_keys)
    all_lon = [cen_lon] + [stars[k]["lon"] for k in star_keys] + [m["lon"] for m in markers]
    all_lat = [cen_lat] + [stars[k]["lat"] for k in star_keys] + [m["lat"] for m in markers]
    attempts = 0
    for cand in _candidate_scenes(cat, bbox_wgs84, t_onset, lon):
        attempts += 1
        if attempts > C.GOES_MAX_READ_ATTEMPTS:
            break
        cand = pc.sign(cand)
        try:
            with rasterio.open(cand.assets[_BAND].href) as src:
                cbx = transform_bounds("EPSG:4326", src.crs, *bbox_wgs84)
                if not all(np.isfinite(cbx)):
                    continue
                pxs, pys = transform("EPSG:4326", src.crs, all_lon, all_lat)
                if not (np.all(np.isfinite(pxs)) and np.all(np.isfinite(pys))):
                    continue
                cwin = from_bounds(*cbx, transform=src.transform)
                cdata = src.read(1, window=cwin, out_shape=(out_px, out_px),
                                 resampling=Resampling.bilinear, boundless=True,
                                 fill_value=src.nodata if src.nodata is not None else -1)
                tags = src.tags(bidx=1)
                cscale = float(tags.get("scale_factor", 1.0))
                coffset = float(tags.get("add_offset", 0.0))
                cnodata = src.nodata
        except Exception:
            continue
        cvalid = cdata != cnodata if cnodata is not None else np.ones_like(cdata, dtype=bool)
        if not cvalid.any():
            continue
        # Post-read sanity floor: GOES_DAYLIGHT_LOCAL_H is a coarse hour band, and a
        # scene right at its dawn/dusk edge can still be too dark to be useful even
        # though the check above let it through.
        crefl = cdata.astype("float32")[cvalid] * cscale + coffset
        if float(np.median(crefl)) < C.GOES_MIN_DAYLIGHT_REFL:
            continue
        it, data, bx = cand, cdata, cbx
        win_transform = rasterio.transform.from_bounds(*cbx, out_px, out_px)
        px_all = [(~win_transform * (x, y)) for x, y in zip(pxs, pys)]
        cen_px, lead, marker_px = px_all[0], px_all[1:n_lead], px_all[n_lead:]
        star_px = dict(zip(star_keys, lead))
        src_meta = (cscale, coffset, cnodata)
        break
    if it is None:
        return None
    scale, offset, nodata = src_meta
    # Methods frequently agree closely (D-66: closest/before landed identically in
    # 91% of 2023 events) with the report itself - without this, agreeing markers
    # (the reported point included) stack exactly on top of each other and only the
    # last-drawn one is ever visible, silently hiding the agreement instead of
    # showing it. Fan overlapping markers into a small rosette around their shared
    # position so every one stays visible. `cen_px` (the TRUE reported position) is
    # kept separately for the vignette-circle maths below; only the VISUAL crosshair
    # (and v3 star) position moves. v3 rides in the same lead slot used above.
    declut = _declutter([cen_px] + [star_px[k] for k in star_keys] + list(marker_px))
    cen_px_visual, lead_visual, marker_px = declut[0], declut[1:n_lead], declut[n_lead:]
    star_px_visual = dict(zip(star_keys, lead_visual))

    valid = data != nodata if nodata is not None else np.ones_like(data, dtype=bool)
    coverage = float(valid.mean())

    # Raw digital number -> reflectance factor (CMI product convention), then a mild
    # percentile stretch + gamma: linear reflectance leaves convective tops looking
    # flat grey: this is a legibility choice, not a scientific one - the underlying
    # value is a straight linear radiometric quantity, no synthesis or blending.
    refl = data.astype("float32") * scale + offset
    lo_r = float(np.percentile(refl[valid], 2)) if valid.any() else 0.0
    hi_r = max(lo_r + 0.05, float(np.percentile(refl[valid], 99.5)) if valid.any() else 0.7)
    norm = np.clip((refl - lo_r) / (hi_r - lo_r), 0.0, 1.0) ** 0.85
    gray = (norm * 255).astype(np.uint8)
    gray[~valid] = 0

    img = Image.fromarray(gray, "L").convert("RGB")
    for (px, py), m in zip(marker_px, markers):
        _mark_matched(img, (px, py), m["radius_km"], m["method"])
    _mark_reported(img, cen_px_visual)
    for k in ("v3", "v4"):                 # drawn last of all (D-70/D-78), v4 topmost
        if k in star_px_visual:
            _mark_star(img, star_px_visual[k], _STAR_COLORS[k])

    # Circular vignette, matching satellite.py's presentation of the fire imagery.
    yy, xx = np.mgrid[0:out_px, 0:out_px]
    px_size = (bx[2] - bx[0]) / out_px
    on_circle = (np.hypot(xx - cen_px[0], yy - cen_px[1]) * px_size <= C.GOES_RADIUS_KM * 1000.0) & valid

    bg = Image.new("RGB", img.size, _BG)
    bg.paste(img, mask=Image.fromarray((on_circle * 255).astype(np.uint8)))
    # Legend AFTER the vignette paste, in the square canvas's corner - the corner
    # sits outside the inscribed circle, so drawing it before the mask (like the
    # markers, which are always well inside the circle) would just erase it.
    _draw_legend(bg)
    bg.save(out_path, "JPEG", quality=C.GOES_JPEG_QUALITY, optimize=True)

    t_utc = t_onset.tz_localize("UTC") if t_onset.tzinfo is None else t_onset
    return {
        "file": os.path.relpath(out_path, C.ROOT).replace(os.sep, "/"),
        "date": it.datetime.strftime("%Y-%m-%d %H:%M"),
        "dt_h": round((it.datetime - t_utc).total_seconds() / 3600.0, 2),
        "coverage": round(coverage, 3),
        "n_markers": len(markers),
        "has_v3": stars.get("v3") is not None,
        "has_v4": stars.get("v4") is not None,
    }


def _declutter(px_list: list[tuple[float, float]], min_sep: float = 18.0) -> list[tuple[float, float]]:
    """Greedily cluster points within `min_sep` px of each other, then fan any
    cluster of 2+ into a small evenly-spaced ring around its centroid so every
    marker stays visible even when the true positions coincide almost exactly."""
    clusters: list[list[int]] = []
    centroids: list[tuple[float, float]] = []
    for i, (x, y) in enumerate(px_list):
        for c, (cx, cy) in enumerate(centroids):
            if np.hypot(x - cx, y - cy) <= min_sep:
                clusters[c].append(i)
                n = len(clusters[c])
                centroids[c] = ((cx * (n - 1) + x) / n, (cy * (n - 1) + y) / n)
                break
        else:
            clusters.append([i])
            centroids.append((x, y))

    out = list(px_list)
    ring_r = min_sep * 0.55
    for idxs, (cx, cy) in zip(clusters, centroids):
        if len(idxs) == 1:
            continue
        for k, i in enumerate(idxs):
            ang = 2 * np.pi * k / len(idxs)
            out[i] = (cx + ring_r * np.cos(ang), cy + ring_r * np.sin(ang))
    return out


def _mark_reported(img: Image.Image, px_py, size: int = 7) -> None:
    """Red crosshair at the AS-REPORTED point - same style D-58's map circle uses for
    the reported point, so the two stay visually consistent."""
    px, py = px_py
    draw = ImageDraw.Draw(img)
    draw.line([(px - size, py), (px + size, py)], fill=_REPORTED_COLOR, width=2)
    draw.line([(px, py - size), (px, py + size)], fill=_REPORTED_COLOR, width=2)


def _mark_matched(img: Image.Image, px_py, radius_km: float, method: str, r: int = 8) -> None:
    """Small filled dot for one (radius, method) v2 position: colour = radius,
    number = method (D-69)."""
    px, py = px_py
    color = RADIUS_COLOR.get(radius_km, (200, 200, 200))
    n = METHOD_NUMBER.get(method, 0)
    draw = ImageDraw.Draw(img)
    draw.ellipse([px - r, py - r, px + r, py + r], fill=color, outline=(0, 0, 0), width=1)
    font = _font(11)
    tb = draw.textbbox((0, 0), str(n), font=font)          # (left, top, right, bottom)
    draw.text((px - (tb[0] + tb[2]) / 2, py - (tb[1] + tb[3]) / 2),
             str(n), fill=(0, 0, 0), font=font)


def _mark_star(img: Image.Image, px_py, color: tuple, r: int = 10) -> None:
    """Filled 5-point star at a single-position match (D-70 v3 / D-78 v4) - a
    single position, not a (radius, method) grid, so a distinct shape rather than
    another dot in the v2 colour/number scheme. Colour picks which version."""
    px, py = px_py
    pts = []
    for k in range(10):
        ang = -np.pi / 2 + k * np.pi / 5      # start pointing straight up
        rad = r if k % 2 == 0 else r * 0.42   # outer / inner points alternate
        pts.append((px + rad * np.cos(ang), py + rad * np.sin(ang)))
    draw = ImageDraw.Draw(img)
    draw.polygon(pts, fill=color, outline=(0, 0, 0), width=1)


_FONT_CACHE: dict[int, "ImageFont.ImageFont"] = {}


def _font(size: int):
    if size not in _FONT_CACHE:
        try:
            _FONT_CACHE[size] = ImageFont.load_default(size=size)
        except TypeError:                    # Pillow < 10.1: no size arg
            _FONT_CACHE[size] = ImageFont.load_default()
    return _FONT_CACHE[size]


def _draw_legend(img: Image.Image, margin: int = 8) -> None:
    """Fixed legend in the bottom-left corner, same on every image: radius colours
    and method numbers (D-69), plus the reported-point symbol."""
    draw = ImageDraw.Draw(img)
    font = _font(12)
    lines = [("Reported:", None, "+ red"),
             ("Match v3 (ladder):", ("star", _V3_COLOR), None),
             ("Match v4 (4-overpass):", ("star", _V4_COLOR), None)] + \
        [(f"{r:g} km:", RADIUS_COLOR[r], None) for r in sorted(RADIUS_COLOR)] + \
        [("Method:", None, "1 closest  2 before  3 before+after  4 before2+after2")]
    h = img.height
    lh = 15
    y = h - margin - lh * len(lines)
    for label, swatch, extra in lines:
        x = margin
        if isinstance(swatch, tuple) and swatch[0] == "star":
            star_color = swatch[1]
            draw.polygon([(x + 5, y + 3), (x + 6, y + 7), (x + 10, y + 7), (x + 7, y + 10),
                          (x + 8, y + 14), (x + 5, y + 11), (x + 2, y + 14), (x + 3, y + 10),
                          (x, y + 7), (x + 4, y + 7)], fill=star_color, outline=(0, 0, 0))
            x += 15
        elif swatch is not None:
            draw.ellipse([x, y + 3, x + 10, y + 13], fill=swatch, outline=(0, 0, 0))
            x += 15
        draw.text((x, y), label, fill=(255, 255, 255), font=font,
                  stroke_width=1, stroke_fill=(0, 0, 0))
        if extra:
            tb = draw.textbbox((0, 0), label, font=font)
            draw.text((x + tb[2] - tb[0] + 5, y), extra, fill=(255, 255, 255), font=font,
                     stroke_width=1, stroke_fill=(0, 0, 0))
        y += lh
