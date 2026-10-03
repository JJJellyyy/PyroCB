The interactive map as deployed on 2026-10-03 with 2013-2025 (2024-25 from the extended fire atlas), saved
when the project was rolled back to 2013-2023. Includes the pyroCb inventory extension, Copernicus land
cover for every year, flowing arrows, the oklab GIF fix. tail.html here also already contains the
"Active pyroCb pixels" layer code (its data was never built for this version).

To view: open fire_pcb_map.html over http (python -m http.server) - needs ./pixattr (copied here, 13-year scales), ../sat_images, ../cloud_images.
To restore: copy map_src/* back to ../map_src/ and run python map_tool.py join.
