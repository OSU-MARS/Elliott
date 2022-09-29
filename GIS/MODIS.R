## MCD18A2 - MODIS/Terra+Aqua Photosynthetically Active Radiation Daily/3-Hour L3 Global 1km SIN Grid
# https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MCD18A2/ for data
#   wget with OAuth token
#     https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD18A2/2001/001/MCD18A2.A2001001.h09v04.006.20171771945.hdf
#     https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD18A2/2001/002/MCD18A2.A2001002.h09v04.006.20171771946.hdf
#     ...
#   or use data download tool when it's working
#     Tiles in the tool are drawn inaccurately with the Elliott crossing the h08v04-h09v04 boundary.
# https://wiki.earthdata.nasa.gov/display/DAS/HEG%3A++HDF-EOS+to+GeoTIFF+Conversion+Tool provides HDF to GeoTiff conversion but 
#   requires Java and fragile setup. It's likely more robust, at least for HDF v5 files, to manipulate them in R with terra or
#   in QGIS. Reprojection from the MODIS to EPSG:6557 in terra dramatically modifies cell shapes, so maintain MODIS CRS and use
#   zonal statistics.
library(dplyr)
library(terra)

dayTile = merge(rast(file.path(getwd(), "GIS/MODIS/MCD18A2.A2001001.h09v04.006.20171771945.hdf")), # GMT_nn00_PAR in W/m²
                rast(file.path(getwd(), "GIS/MODIS/MCD18A2.A2001001.h08v04.006.20171771945.hdf")))
#writeRaster(dayTile, file.path(getwd(), "GIS/MODIS/2001.01.01.tif"), overwrite = TRUE)
dayTile = merge(rast(file.path(getwd(), "GIS/MODIS/MCD18A2.A2001002.h09v04.006.20171771946.hdf")),
                rast(file.path(getwd(), "GIS/MODIS/MCD18A2.A2001002.h08v04.006.20171771945.hdf")))
#writeRaster(dayTile, file.path(getwd(), "GIS/MODIS/2001.01.02.tif"), overwrite = TRUE)

dayElliott = crop(project(dayTile, "EPSG:6556"), ext(c(102000, 138000, 190000, 226000) + 10000 * c(-1, 1, -1, 1)))

dayElliott$PAR = dayElliott$GMT_0000_PAR + dayElliott$GMT_0300_PAR + dayElliott$GMT_0600_PAR + dayElliott$GMT_0900_PAR + dayElliott$GMT_1200_PAR + dayElliott$GMT_1500_PAR  + dayElliott$GMT_1800_PAR + dayElliott$GMT_2100_PAR

writeRaster(dayElliott$PAR, file.path(getwd(), "GIS/MODIS/MCD18A2.A2001001.h09v04.006.20171771945.tif"), overwrite = TRUE)
#writeRaster(dayElliott$PAR, file.path(getwd(), "GIS/MODIS/MCD18A2.A2001002.h09v04.006.20171771946.tif"), overwrite = TRUE)
