Tampa Bay Watershed — Rainfall Extremes Trend Analysis
Replicates key analyses from:
Wright, D.B., Bosma, C.D., & Lopez-Cantu, T. (2019). "U.S. Hydrologic Design Standards Insufficient Due to Large Increases in Frequency of Rainfall Extremes." Geophysical Research Letters, 46, 8144-8153. https://doi.org/10.1029/2019GL083235

Scope: All GHCN-Daily stations within the Tampa Bay estuary watershed boundary
Study period: 1900–2025 (extending the paper's primary analysis window)

Analyses replicated:
   1. GHCN station retrieval and spatial filter to Tampa Bay watershed
   2. Exceedance counting relative to NOAA Atlas 14 IDF estimates (Section 3.1)
   3. Negative binomial regression trend analysis (Section 3.1)
   4. Rainstorm cluster identification (Section 3.2 / Section 2.4)
   5. Design-vs-observed ARI comparison (Section 3.3 / Section 2.5)
   6. Publication-quality figures for all analyses

Dependencies (install once before sourcing):
   install.packages(c("tidyverse", "sf", "terra", "rnoaa", "MASS", "pscl", "zoo", "lubridate", "ggplot2", "patchwork", "scales", "httr", "jsonlite", "tigris", "units"))

Data sources fetched automatically at runtime:
   1. GHCN-Daily via rnoaa::ghcnd_*()
   2. NOAA Atlas 14 IDF values via NOAA Precipitation Frequency Data Server (PFDS) REST API  <https://hdsc.nws.noaa.gov/pfds/>
   3. Tampa Bay watershed boundary via internal shapefile, fallback is the USGS StreamStats or NHD boundary (a local GeoJSON fallback is provided if API is unavailable)
