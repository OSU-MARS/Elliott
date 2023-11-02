library(dplyr)
library(lidR)
library(future)
library(readxl)
library(sf)
library(terra)

chunkIndex = 18
setwd(file.path(getwd(), "../../.."))

## 2015-16 cruise plot metrics
# Plot metrics are memory compact, so memory constraint on worker count (12 workers -> ~20 GB DDR). However, compute times
# are easily long.
# workers  5950X runtime, hours
#          1000 plots  18312 plots, estimated
#  1                   200
#  2       
#  4                   51
#  8       1.7         31                      two concurrent eight worker jobs
# 12                   22
plotStartTime = Sys.time()
plan(multisession, workers = 8)

chunkSize = 1000

# lidR doesn't translate across CRSes so plot centers and point clouds must be in the same CRS.
elliottNormalized = readLAScatalog("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz normalized")
trees2016 = read_xlsx("trees/Elliott final cruise records 2015-16.xlsx", sheet = "CRUISERECS")
cruisePlots2016 = st_transform(st_read("GIS/Trees/2015-16 cruise/CruisePlots_All_20151211.gpkg"), crs = st_crs(6557)) %>%
  filter(PltInteger %in% unique(trees2016$PlotID))

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, nrow(cruisePlots2016))
cat(paste0("\nProcessing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", endIndex - startIndex, " plots...\n")) # leading newline because readLAScatalog() doesn't write trailing newline

plotStdmetrics = plot_metrics(elliottNormalized, .stdmetrics, cruisePlots2016[startIndex:endIndex, ], radius = 11.284)
write_sf(plotStdmetrics, paste0("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics/Elliott stdmetrics plot chunk ", chunkIndex, ".gpkg"), delete_layer = TRUE)

cat(paste0("Standard point cloud metrics for ", nrow(cruisePlots2016), " cruise plots in ", format(Sys.time() - plotStartTime), ".\n"))

# deprecated due to low lidR performance, reliability, and lack of non-normalized support
# Use Get-GridMetrics in Clouds instead.
# if (gridMetrics)
# {
#   ## forest level pixel metrics: 20 m ABA grid (iLand aligned within 1 mm)
#   gridStartTime = Sys.time()
#   plan(multisession, workers = 8)
#   
#   elliottNormalized = readLAScatalog("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz normalized") # EPSG:6557 + EPSG:8228
#   elliottAbaGrid = rast(file.path(getwd(), "GIS/Trees/Elliott ABA grid 20 m.tif")) # EPSG:6556 + no vertical CRS
#   chunkBoundary = read_sf(file.path(getwd(), "GIS/Trees/Elliott ABA grid 20 m.gpkg"), layer = "chunk grid 1 km")[chunkIndex, ] # EPSG:6556
#   chunkCells = trim(crop(elliottAbaGrid, chunkBoundary)) # trim() excludes rectangular no data regions in boundary chunks but won't follow (pixelated) curves
#   chunkCells = project(chunkCells, st_crs(elliottNormalized)$wkt) # convert to EPSG:6557 + EPSG:8228 to match LiDAR tiles
#   
#   cat(paste0("Processing chunk ", chunkIndex, " with ", prod(dim(chunkCells)[1:2]), " cells...\n"))
#   chunkStdmetrics = pixel_metrics(elliottNormalized, .stdmetrics, chunkCells)
#   #chunkGrid = st_transform(read_sf(file.path(getwd(), "GIS/Trees/Elliott ABA grid 20 m.gpkg"), layer = "grid 20 m"), crs = st_crs(6557)) %>% filter(chunk == chunkIndex)
#   #chunkStdmetrics = polygon_metrics(elliottNormalized, .stdmetrics, chunkGrid[1, ]) # not implemented in lidR
#   writeRaster(chunkStdmetrics, sprintf("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics/Elliott stdmetrics 20 m chunk %03i.tif", chunkIndex), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
#   
#   cat(paste0("Standard point cloud metrics at 20 m resolution in ", format(Sys.time() - gridStartTime), ".\n"))
# }
