# pyroCb × FireAtlas

Associating every **pyrocumulonimbus (pyroCb)** event in the boreal zone with the fire that
produced it, and with the fire's state — footprint, active-burn area, ignition history — at the
moment of onset. The endpoint is one clean row per pyroCb event carrying ~50 environmental
predictors (ERA5, topography, hydrology, 30 m land cover, …) for downstream analysis.

## Layout

```
pyrocb/                    importable pipeline (see below)
make_notebook.py           emits "Data exploration.ipynb" — the narrative + build driver
Data exploration.ipynb     run top-to-bottom to reproduce everything
CHANGELOG.md               what changed, newest first
DECISIONS.md               every methodological choice and why  (D-1 … D-26)
FINDINGS.md                what the data itself told us          (F1 … F11)

pyrocb_fire_attributes.tsv / .xlsx   one row per pyroCb event  (484 × 122)
pyrocb_columns.tsv                   data dictionary for the above
pyrocb_active_areas.gpkg             active-burn polygon per (event, ±6/12/24 h window)
fire_pcb_map.html                    interactive map (latest)
maps/                                dated snapshots of the map, one per release
```

### The `pyrocb` package

| module | role |
|---|---|
| `config`   | paths, constants, year resolution, the parquet cache helper |
| `io`       | readers for the four FireAtlas products + the geometry validity pipeline |
| `timebins` | half-day (AM/PM overpass) step handling; the "snapshot before onset" rule |
| `matching` | pyroCb event → fire (`mergid`); containment, then nearest active fire |
| `bands`    | concentric growth bands, split at every pyroCb, cumulative-union construction |
| `active`   | active-burn polygon at onset, from NFP pixel footprints |
| `variables`| the per-event attribute table |
| `pipeline` | per-year orchestration with caching, then export |

## Reproducing

```bash
python make_notebook.py           # regenerate the notebook from source
python -m pip install -r requirements.txt   # geopandas, pyogrio, shapely>=2.1, pyarrow, openpyxl, nbformat
# then run "Data exploration.ipynb", or:
jupyter nbconvert --to notebook --execute "Data exploration.ipynb"
```

Cold build ≈ 10–12 min (all years); warm reruns seconds (per-(stage, year) parquet cache in
`cache/`, keyed on input mtime + the relevant config hash). `C.clear_cache()` to wipe.

Scope knobs in `pyrocb/config.py`: `YEARS` (map; `None` → years with fire data *and* pyroCb
records = 2013–2023) and `SNAPSHOT_YEARS` (heavy analysis; `None` → same as `YEARS`).

## Source data — not in the repo

The raw FireAtlas / FEDS archives (~30 GB) and the 30 m land-cover rasters (~2 GB each) are
external inputs and `.gitignore`d. Expected layout:

```
{year}/{year}/Snapshot/YYYYMMDD{AM,PM}.gpkg   12-hourly cumulative perimeters (EPSG:3571)
{year}/{year}/NFP/YYYYMMDD{AM,PM}_NFP.txt      new fire pixels, real datetime_utc
final_perims/final_perims/final_viirs{year}.gpkg   season-final perimeters (EPSG:4326)
ignitions/ignitions/ign{year}.gpkg                 first-detection footprints (EPSG:3571)
PCBAtlas.xlsx                                       NRL pyroCb inventory (Peterson et al. 2025)
landcover-YYYY-classification.tif                   30 m Canada land cover
```

Attribute definitions: `fireatlas_attribute_table_v1.pdf`, `PCB readme.pdf`.

## Editing the map

`fire_pcb_map.html` is one 47 MB file. Edit its parts instead:

```bash
python map_tool.py split   # -> map_src/head.html (CSS + markup), tail.html (JS), data.json (payload)
python map_tool.py join    # -> fire_pcb_map.html
node map_check.js          # syntax-check the JS
```

Data builders write into `map_src/data.json` or `pixattr/`: `build_c3s_landcover.py` (land cover),
`build_pix_wide.py` (pyroCb detection grids), `build_pixel_attrs.py` (pixel variables). Deploy
the map with `pixattr/` and `sat_images/` next to it (it fetches them at runtime).

## Version control

Git with **LFS** for `*.html`, `*.gpkg`, `*.xlsx`, `*.pdf` so the map (rebuilt often, ~10–26 MB)
and the GeoPackage don't bloat the pack. Every release drops a dated copy in `maps/`.
