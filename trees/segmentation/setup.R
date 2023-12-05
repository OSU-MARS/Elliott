#.libPaths(.libPaths()[2])
#install.packages(c("arrow", "caret", "data.table", "FNN", "lidR", "ranger", "sf", "terra", "tidyterra"))
library(data.table)
library(dplyr)
library(future)
library(ggplot2)
library(lidR)
library(magrittr)
library(patchwork)
library(readxl)
library(sf)
library(stringr)
library(terra)
library(tibble)
library(tidyterra)
library(tidyr)

metrics_rgbn = function(returnNumber, r, g, b, nir)
{
  firstReturns = which(returnNumber == 1)
  return(data.table(r = sqrt(mean(r[firstReturns]^2)), g = sqrt(mean(g[firstReturns]^2)), b = sqrt(mean(b[firstReturns]^2)), nir = sqrt(mean(nir[firstReturns]^2))))
}

metrics_rgbnc = function(classification, returnNumber, r, g, b, nir)
{
  firstReturns = which((classification == 1) & (returnNumber == 1))
  return(data.table(r = sqrt(mean(r[firstReturns]^2)), g = sqrt(mean(g[firstReturns]^2)), b = sqrt(mean(b[firstReturns]^2)), nir = sqrt(mean(nir[firstReturns]^2))))
}

#metrics_zero = function(returnNumber, z) # lidR testing function
#{
#  return(tibble(n1 = 0, n2 = 0))
#}

metrics_zn = function(returnNumber, z)
{
  return(data.table(n1 = sum(returnNumber == 1), n2 = sum(returnNumber == 2), n3 = sum(returnNumber == 3), n4 = sum(returnNumber == 4), n5 = sum(returnNumber == 5),
                    z1 = max((returnNumber == 1) * z), z2 = max((returnNumber == 2) * z), z3 = max((returnNumber == 3) * z), z4 = max((returnNumber == 4) * z), z5 = max((returnNumber == 5) * z)))
}

remove_redundant_tiles = function(directoryToRetain, directoryToPrune)
{
  tileNamesEligibleToObsolete = list.files(directoryToPrune, pattern = "\\.laz$")
  tileNamesToRetain = list.files(directoryToRetain, pattern = "\\.laz$")
  
  tileNamesToObsolete = intersect(tileNamesToRetain, tileNamesEligibleToObsolete)
  cat(paste0("Removing ", length(tileNamesToObsolete), " tiles from ", directoryToPrune, "...\n"))
  for (tileName in tileNamesToObsolete)
  {
    #cat(paste0("Remove", file.path(directoryToPrune, tileName)))
    file.remove(file.path(directoryToPrune, tileName))
  }
}

segment_cloud = function(pointCloud, isNormalized = FALSE)
{
  chm = rasterize_canopy(pointCloud, res = 1, algorithm = p2r(subcircle = 0)) # compute the canopy height model using a pixel of 0.5 units, in this case feet
  #plot(st_geometry(standBuffer30m))
  #plot(chm,add=TRUE)
  
  if (isNormalized)
  {
    chmHeightQuantiles = quantile(values(chm), probs = c(0.1, 0.85, 0.95), na.rm = TRUE, names = FALSE) # point density dependent: 0.1, 0.65, 0.95? for 2009 flight?
    names(chmHeightQuantiles) = c("minimumTreeHeight", "q85", "q95")
  
    get_window_size = function(z)
    {
      return(z * log((chmHeightQuantiles[3] - chmHeightQuantiles[2])) / 25 + (chmHeightQuantiles[3] - chmHeightQuantiles[2]) / 10)
    }
    localMaxima = locate_trees(chm, lmf(get_window_size, hmin = chmHeightQuantiles[1]), uniqueness = "incremental") # identifies the tree tops from the CHM using the local maximum filter with a window of 6 feet and height >8 ft
    #localMaxima = locate_trees(chm,lmf(chmHeightQuantiles[3]/25,hmin=chmHeightQuantiles[1]),uniqueness = "incremental")
  } else {
    localMaxima = locate_trees(chm, lmf(get_window_size, hmin = chmHeightQuantiles[1]), uniqueness = "incremental")
  }
  #dim(localMaxima)
  #plot(localMaxima, pch=20, add=TRUE)
  
  minimumTreeHeight = chmHeightQuantiles["minimumTreeHeight"] # ft
  maxCrownDiameter = chmHeightQuantiles["q95"] / 8 # ft
  segmentationAlgorithm = dalponte2016(chm, localMaxima, th_tree = minimumTreeHeight, max_cr = maxCrownDiameter)
  #segmentationAlgorithm = watershed(chm, th_tree=15)
  #segmentationAlgorithm = li2012(dt1=1, dt2=2, Zu=30, hmin=10, speed_up = 8)
  #segmentationAlgorithm = silva2016(chm,localMaxima, max_cr_factor = 0.5, exclusion = 0.3, ID="treeID")
  
  # TODO: how to avoid  GDAL message: Setting GeoTIFF 1.1 requires libgeotiff >= 1.6?
  pointCloud = segment_trees(pointCloud, algorithm = segmentationAlgorithm) # segment point cloud: identifies individual tree crowns
  return(list(standPoints = pointCloud, chm = chm, minimumTreeHeight = minimumTreeHeight))
}

stdmetrics_z_simple = function(z) # reduced from https://github.com/r-lidar/lidR/blob/master/R/metrics_stdmetrics.R
{
  n = length(z)

  probs = seq(0, 1.0, by = 0.1)
  zq = as.list(stats::quantile(z, probs))
  names(zq) = paste0("zq", 100 * probs)

  return(c(n = length(z), zq, zmean = mean(z), zsd = sd(z)))
}

# setup
# 2021 tile naming convention: s<easting>w<northing> where easting ≥ 03450 and northing ≥ 06360 increment in steps of 30 (units are 100 ft)
dataPath = "E:/Elliott/GIS/DOGAMI/2009 OLC South Coast"
#dataPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"

dtm2021nv5 = rast("GIS/DOGAMI/2021 OLC Coos County/DTM_NV5_2021.tif")
elliottStands = st_read("GIS/Planning/ESRF_Stands062022Fixed.shp") # June 2022 stand boundaries
elliottStands = st_transform(elliottStands, 6557) # override EPSG:2992 to 6557 - same projection and units, just different EPSG
elliottLidarTiles = readLAScatalog(file.path(dataPath, "points"), progress = FALSE) # progress = FALSE turns off plotting of tiles during clip_roi() calls, many warnings: only 2 bytes until point block

# check assignment of stand IDs
# STD_ID: Oregon Department of Forestry-Oregon Department of State Lands stand ID
# OSU_SID: Oregon State University revised stand ID, typically a 10 prefix added to ODF-DSL STD_ID
if (nrow(elliottStands) != length(unique(elliottStands$OSU_SID)))
{
  stop("Stand IDs are not unique to each inventory polygon.")
}

#plot(dtm2021nv5)
#plot(st_geometry(elliottStands), add=TRUE)
#head(elliottStands)
#st_crs(elliottStands)

#stands.mis = st_read("./People/MARS/Strimbu/TreeSeg2/Results/StandsMissed.shp") #Stands that were inventoried and not inventoried
#crashed = stands.mis[which(is.na(stands.mis$Freq)),]
#plot(st_geometry(crashed), add=TRUE, col="red")
#crashed = st_transform(crashed, 6557)
#elliottStands.id = elliottStands[,c("StandID", "geometry", "Area_ft2","OSU_SID","T4I")]
#elliottStands.id = elliottStands.id[,-(1)]
#head(elliottStands.id)
#st_crs(elliottStands.id)


## tile level pixel metrics, 5950X: get_lidr_threads() = 16
dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileName = "s04020w06690.laz" # 40-45 year old plantations + 130+ reserve stands
tileName = "s04260w06720.laz" # 4, 5, 13, 46, 52 year plantations + 95 & 120 year reserves
tile = readLAS(file.path(dataPath, "pointz", tileName))
tile = filter_poi(tile, (Classification != 7L & (Classification != 18L)))

#metricsZero = pixel_metrics(tile, ~metrics_zero(ReturnNumber, Z), res = 3.28084 * 1) # 1.4 min with data.table(), 5.3 min with tibble() -> 3.8x penalty
#metrics = pixel_metrics(tile, ~stdmetrics(X, Y, Z, Intensity, ReturnNumber, Classification), dz = 3.28084 * 1, res = 3.28084 * 1) # 11.6 min / tile @ 1 m resolution with dz = 1 foot -> 126 MB GeoTIFF
#metricsCtrl = pixel_metrics(tile, .stdmetrics_ctrl, res = 3.28084 * 1) # only n and area but just a few seconds, error with ~stdmetrics_ctrl(X, Y, Z)
#metricsReturns = pixel_metrics(tile, .stdmetrics_rn, res = 3.28084 * 1) # 1.0 min -> 15 MB GeoTIFF
#metricsRgbn = pixel_metrics(tile, ~metrics_rgbn(ReturnNumber, R, G, B, NIR), res = 1) # 20.3 min -> 64 MB GeoTIFF
#metricsZ = pixel_metrics(tile, .stdmetrics_z, dz = 3.28084 * 1, res = 3.28084 * 1) # 5.53 min -> 71 MB GeoTIFF, z mean, σ, quantiles, skew, kurtosis, entropy, probabilities, error with ~stdmetrics_z(X, Y, Z)
#metricsZ = pixel_metrics(tile, ~stdmetrics_z_simple(Z), res = 3.28084 * 1) # 2.5 min -> 34 MB GeoTIFF
#metricsZn = pixel_metrics(tile, ~metrics_zn(ReturnNumber, Z), res = 3.28084 * 1) # 2.4 min -> 16.9 MB GeoTIFF

metricsStartTime = Sys.time()
#metricsRgbn = pixel_metrics(tile, ~metrics_rgbn(ReturnNumber, R, G, B, NIR), res = 1) # 20.3 min -> 64 MB GeoTIFF
metricsRgbnc = pixel_metrics(tile, ~metrics_rgbnc(Classification, ReturnNumber, R, G, B, NIR), res = 1) # ~21 min -> 64 MB GeoTIFF
format(Sys.time() - metricsStartTime)


tileBaseName = str_remove(tileName, "\\.laz$")
writeRaster(metrics, file.path(dataPath, "metrics 1m", paste0(tileBaseName, ".tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsCtrl, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " ctrl.tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsReturns, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " rn.tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsRgbn, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " rgbn.tif")), datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsRgbnc, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " rgbnc.tif")), datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsZ, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " Z.tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
writeRaster(metricsZn, file.path(dataPath, "metrics 1m", paste0(tileBaseName, " Zn.tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)


## remove old versions of .laz files which have been updated
remove_redundant_tiles("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz RGBIR", "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz")


## testing
# stand id   2021 tiles                   notes
# 443        s03960w07050, s03990w07050   7.5 ha, 150 y natural regen, N central
# 605        s03990w06930, s04020w06930   40 ha, young plantation, steep, S central
# 1874       s03510w06480, s03540w06480   1.1 ha, natural regen, flat, SW corner
tileBaseName = "s03960w07050"
tile = readLAS(file.path(dataPath, "Points", paste0(tileBaseName, "_las.las")))
#table(tile$Classification)
#nonNormalizedSegmentation = segment_cloud(tile)
#nonNormalizedTrees = get_trees(nonNormalizedSegmentation)
#st_write(tileTrees, file.path(dataPath, "segmentation", paste0(tileName, " trees non-normalized.gpkg")), append = FALSE) # only 17 with default Dalponte

tileDsm = rasterize_canopy(tile, res = 1) # get tile's digital surface model at 1 foot resolution using default p2r() algorith,m
tileDtm = project(resample(dtm2021nv5, tileDsm, threads = TRUE), crs(tileDsm)) # bilinear interpolation by default, use projection to close difference in vertical CRS TODO: consider cubic
tileChm = tileDsm - tileDtm

ggplot() +
  geom_spatraster(aes(fill = Z), tileDsm) +
  coord_sf(datum = st_crs(tileChm)) +
  labs(fill = "canopy\nelevation, ft") +
  scale_fill_distiller(palette = "Greys") +
ggplot() +
  geom_spatraster(aes(fill = Z), tileChm) +
  coord_sf(datum = st_crs(tileChm)) +
  labs(fill = "tree\nheight, ft") +
  scale_fill_viridis_c() +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1)

normalizedTile = classify_noise(tile, ivf(res = 3, n = 1)) # identifies noise using IVF method
dtm = crop(dtm2021nv5, st_buffer(st_as_sfc(st_bbox(tile)), 30 + 3))
normalizedTile = normalize_height(normalizedTile, dtm)
normalizedTile = filter_poi(normalizedTile, Classification == 1L, Z >= 0, Z <= 400) # exclude ground points and high outliers, Z range in feet

normalizedSegmentation = segment_cloud(normalizedTile, isNormalized = TRUE)
normalizedTrees = get_trees(normalizedSegmentation)
st_write(normalizedTrees, file.path(dataPath, "segmentation", paste0(tileName, " trees normalized.gpkg")), append = FALSE) # 31466 trees

get_trees = function(segmentation)
{
  allTreeMetrics = crown_metrics(segmentation$standPoints, func = .stdmetrics, geom = "convex") %>% # crown polygons and various statistics including the percentiles
    rename(crownArea = area) # ft², projected
  allTreeCentroids = st_centroid(st_geometry(allTreeMetrics))
  st_geometry(allTreeMetrics) = allTreeCentroids
  allTreeCentroids = st_coordinates(allTreeCentroids)
  
  treesNearestNeighbor3 = get.knn(allTreeCentroids, 
                                  k = 3, 
                                  algorithm = "brute")
  colnames(treesNearestNeighbor3$nn.index) = c("nn1", "nn2", "nn3")
  colnames(treesNearestNeighbor3$nn.dist) = c("nn1distance", "nn2distance", "nn3distance")
  
  neighborMergeDistance = 7 # ft
  tileTrees = bind_cols(allTreeMetrics,
                        as_tibble(treesNearestNeighbor3$nn.index), 
                        as_tibble(treesNearestNeighbor3$nn.dist)) %>% 
    mutate(STD_ID = 0L, # TODO: join by location to get stand ID
           treeID = as.integer(treeID),
           species = "DF", # default all trees to Douglas-fir
           xCentroid = allTreeCentroids[,"X"], yCentroid = allTreeCentroids[,"Y"],
           x1 = xCentroid[nn1], x2 = xCentroid[nn2], x3 = xCentroid[nn3],
           y1 = yCentroid[nn1], y2 = yCentroid[nn2], y3 = yCentroid[nn3],
           z1 = zmax[nn1], z2 = zmax[nn2], z3 = zmax[nn3],
           isShorterThanNearNeighbor = ((nn1distance < neighborMergeDistance) & (zmax < z1)) |
             ((nn2distance < neighborMergeDistance) & (zmax < z2)) |
             ((nn3distance < neighborMergeDistance) & (zmax < z3))) %>%
    relocate(STD_ID, treeID, species, xCentroid, yCentroid, zmax, crownArea)
  
  return(tileTrees)
}


## problems with lack of parallel support in terra
#library(furrr)
#options(future.debug = TRUE)
#plan(multisession, workers = 8)
boundingBoxes = future_map(elliottLidarTiles$filename, function(tilePath)
{
  boundingBox = merge(st_as_sfc(st_bbox(readLASheader(tilePath))),
                      tibble(name = stringr::str_remove(fs::path_file(tilePath), "_las.las")))
  return(boundingBox)
}, .options = furrr_options(seed = TRUE))

st_write(st_sf(bind_rows(boundingBoxes)), file.path(dataPath, "Points", "tile index.gpkg"), append = TRUE)

#standIDs = c(1874, 1875)
#future_map(standIDs, function(standID)
#{
#  process_stand(standID)
#})
#ggsave(file.path(dataPath, "segmentation", "604 trees.tif"), device = tiff(compression = "lzw"), width = 30, height = 30, units = "cm", dpi = 600)


## comparative stdmetrics for Clouds unit test
psmeCells = raster::raster(file.path(getwd(), "../tools/Clouds/UnitTests/PSME ABA grid cells.tif"))
psme146 = readLAS(file.path(getwd(), "../tools/Clouds/UnitTests/PSME LAS 1.4 point type 6.las"))
psmeMetrics = pixel_metrics(psme146, ~stdmetrics(X, Y, Z, Intensity, ReturnNumber, Classification, dz = 1, th = 920.52 + 2), zmin = min(psme146$Z), psmeCells)
psmeMetrics[15, 9]
psmeMetrics[15, 10]

## no longer used
#library(FNN) # fast nearest neighbor::get.knn() for cleaning segmentation
# normalize_stand = function(standID)
# {
#   standBoundary = (elliottStands %>% filter(STD_ID == standID))[1]
#   #standBuffer30m = st_buffer(standBoundary, 30) # m
#   #units::set_units(st_area(standBoundary), acres)
#   #plot(st_geometry(standBuffer30m), col="green")
#   #plot(st_geometry(standBoundary), add=TRUE, col="red")
#   
#   ## point cloud normalization and canopy height model creation
#   segmentationBufferWidth = 30 # ft TODO: is a 30 foot radius is too small for correct segmentation of large trees along the stand boundary?
#   standAndBufferPoints = clip_roi(elliottLidarTiles, st_buffer(standBoundary, segmentationBufferWidth))
#   #standBuffer30m = st_transform(standBuffer30m, st_crs(standPoints)) # TODO: why clip a second time?
#   #standPoints = clip_roi(standPoints, st_buffer(standBoundary, segmentationBufferWidth))
#   #standPoints = las_update(standPoints) # update the header
#   
#   standAndBufferPoints = classify_noise(standAndBufferPoints, ivf(res = 3, n = 1)) #identifies noise using ivf method
#   #cat("noise classification...", fill = TRUE)
#   #table(standPoints@data$Classification)
#   #standPoints.sor=classify_noise(standPoints, sor(treeID=10,m=9)) # identifies noise using sor method
#   #table(standPoints@data$Classification,standPoints.sor@data$Classification)
#   
#   dtm = crop(dtm2021nv5, st_buffer(standBoundary, segmentationBufferWidth + 3)) # crop the large DTM to the point cloud extent, adding additional buffer width so all the points normalize_height() needs to process lie over the DTM TODO: why crop from the large DTM?
#   standAndBufferPoints = normalize_height(standAndBufferPoints, dtm)
#   standAndBufferPoints = filter_poi(standAndBufferPoints, Classification == 1L, Z >= 0, Z <= 400) # exclude ground points and high outliers, Z range in feet
#   #plot(decimate_points(standPoints, random(1)))
#   #rm(standPoints)
#   #gc()
#   
#   return(standAndBufferPoints)
# }
# 
# process_stand = function(standID, plotSegmentation = FALSE)
# {
#   ## initial segmentation
#   normalizationStartTime = Sys.time()
#   standAndBufferPoints = normalize_stand(standID)
#   cat(sprintf("%d normalization: %.1f s", standID, difftime(Sys.time(), normalizationStartTime, units = "secs")), fill = TRUE)
#   
#   segmentationStartTime = Sys.time()
#   segmentation = segment_cloud(standAndBufferPoints)
#   minimumTreeHeight = segmentation$minimumTreeHeight
#   standAndBufferPoints = segmentation$standPoints
#   #rm(segmentation)
#   cat(sprintf("%d initial segmentation: %.1f s.", standID, difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
#   #length(standPoints@data$treeID[(standPoints@data$treeID =="NA")])
#   #sum(is.na(standPoints@data$treeID))
#   
#   ## look for trees with outlier points above them, remove outliers if found
#   highPointStartTime = Sys.time()
#   outlierAboveThreshold = 5 # ft
#   outlierAboveTrees = standAndBufferPoints@data %>% filter(is.na(treeID) == FALSE) %>% 
#     group_by(treeID) %>% 
#     reframe(quantiles = c(0.995, 1.0), z = quantile(Z, probs = quantiles, names = FALSE)) %>%
#     pivot_wider(id_cols = "treeID", names_prefix = "q", names_from = "quantiles", values_from = "z") %>%
#     filter((q1 - q0.995) > outlierAboveThreshold)
#   
#   if (nrow(outlierAboveTrees) > 0)
#   {
#     for (treeIndex in 1:nrow(outlierAboveTrees))
#     {
#       treeID = outlierAboveTrees$treeID[treeIndex]
#       zThreshold = outlierAboveTrees$q0.995[treeIndex]
#       standAndBufferPoints@data$Classification[which((standAndBufferPoints@data$Z >= zThreshold) & (standAndBufferPoints@data$treeID == treeID))] = 18L # point class: high
#     }
#     
#     standAndBufferPoints = standAndBufferPoints[(standAndBufferPoints@data$Classification != 7L) & (standAndBufferPoints@data$Classification != 18L),] # remove outlier points
#     standAndBufferPoints = las_update(standAndBufferPoints)
#   }
#   cat(sprintf("%d upper quantile checking: found %d trees with large upper quantile differences in %.1f s.", standID, nrow(outlierAboveTrees), difftime(Sys.time(), highPointStartTime, units = "secs")), fill = TRUE)
#   
#   if (nrow(outlierAboveTrees) > 0)
#   {
#     segmentationStartTime = Sys.time()
#     segmentation = segment_cloud(standAndBufferPoints)
#     standAndBufferPoints = segmentation$standPoints
#     cat(sprintf("%d resegmentation: %.1f s.", standID, difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
#   }
#   
#   #nr.trees.final = length(treeIDs)
#   #TPA = nr.trees.final / st_area(standBuffer30m) * 45360
#   #cat(sprintf("segmentation time: %.1f s.", difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
#   
#   allTreeMetrics = crown_metrics(standAndBufferPoints, func = .stdmetrics, geom = "convex") %>% # crown polygons and various statistics including the percentiles
#     rename(crownArea = area) # ft², projected
#   #allTreeMetrics = st_centroid(allTreeMetrics, of_largest_polygon = FALSE, ensure_within = TRUE) # not needed as crown_metrics() returns a point geometry with centroid coordinates TODO: clear warning st_centroid assumes attributes are constant over geometries
#   #treeCentroidXY = st_coordinates(treeCentroids)
#   #allTreeMetrics %<>% mutate(Xcentroid = treeCentroidXY[,"X"], Ycentroid = treeCentroidXY[,"Y"])
#   
#   
#   ## remove segmented trees whose tops are less than a threshold distance apart (like 7 ft) or are shorter than minimumTreeHeight
#   neighborEliminationStartTime = Sys.time()
#   minimumCrownArea = 15 # ft²
#   acceptedTreeMetrics = allTreeMetrics %>% filter(crownArea > minimumCrownArea, zmax > minimumTreeHeight)
#   acceptedTreeCentroids = st_centroid(st_geometry(acceptedTreeMetrics)) # reduce trees' crown polygons to their centroids (layer's geometry type changes from polygon to point)
#   st_geometry(acceptedTreeMetrics) = acceptedTreeCentroids
#   
#   # find three nearest neighbors based on distance between trees' centroids
#   # $nn.index, which identifies the row index of the tree in the st geometry/data frame, is not the tree number. The matching of index with treeID should be done explicitly.
#   acceptedTreeCentroids = st_coordinates(acceptedTreeCentroids) # flatten point geometry to XY matrix
#   treesNearestNeighbor3 = get.knn(acceptedTreeCentroids, 
#                                   k = 3, 
#                                   algorithm = "brute")
#   colnames(treesNearestNeighbor3$nn.index) = c("nn1", "nn2", "nn3")
#   colnames(treesNearestNeighbor3$nn.dist) = c("nn1distance", "nn2distance", "nn3distance")
#   
#   neighborMergeDistance = 7 # ft
#   standAndBufferTrees = bind_cols(acceptedTreeMetrics,
#                                   as_tibble(treesNearestNeighbor3$nn.index), 
#                                   as_tibble(treesNearestNeighbor3$nn.dist)) %>% 
#     mutate(STD_ID = as.integer(standID),
#            treeID = as.integer(treeID),
#            species = "DF", # default all trees to Douglas-fir
#            xCentroid = acceptedTreeCentroids[,"X"], yCentroid = acceptedTreeCentroids[,"Y"],
#            x1 = xCentroid[nn1], x2 = xCentroid[nn2], x3 = xCentroid[nn3], # while they don't need to be stored explicitly, unpack coordinates of each tree's three nearest neighbors for now
#            y1 = yCentroid[nn1], y2 = yCentroid[nn2], y3 = yCentroid[nn3],
#            z1 = zmax[nn1], z2 = zmax[nn2], z3 = zmax[nn3],
#            isShorterThanNearNeighbor = ((nn1distance < neighborMergeDistance) & (zmax < z1)) |
#              ((nn2distance < neighborMergeDistance) & (zmax < z2)) |
#              ((nn3distance < neighborMergeDistance) & (zmax < z3))) %>%
#     relocate(STD_ID, treeID, species, xCentroid, yCentroid, zmax, crownArea)
#   standAndBufferTrees %<>% filter(crownArea > minimumCrownArea, # ft² TODO: vary minimum crown projection area with tree height as a 15 ft² fixed area seems high for short trees and small for tall ones?
#                                   isShorterThanNearNeighbor == FALSE) %>% 
#     select(-isShorterThanNearNeighbor) # for now, update point cloud is not updated merged tree IDs
#   # TODO: update crown area, centroid, and so on based on merge of dropped trees?
#   #nr.standAndBufferTrees = nrow(standAndBufferTrees)
#   #sum(standAndBufferTrees$nn1Dist < neighborMergeDistance)
#   #sum(standAndBufferTrees$FP)
#   cat(sprintf("%d neighbor elimination time: %.1f s.", standID, difftime(Sys.time(), neighborEliminationStartTime, units = "secs")), fill = TRUE)
#   
#   ## find x and y coordinate of tallest point within the crown
#   # Very time consuming as currently implemented. TODO: is group_by() and summarize(Z == acceptedTreeMetrics$zmax) more performant?
#   #trees.max.xy = as.data.frame(matrix(0, length(standMetrics$treeID), 3))
#   #names(trees.max.xy) = c("Xmax","Ymax", "Zmax")
#   #  
#   #time.xymax = Sys.time()
#   #for (kk in 1:length(standMetrics$treeID)) 
#   #{
#   #  treePoints = filter_poi(standPoints, treeID == standMetrics$treeID[kk])
#   #  trees.max.xy[kk,] = treePoints@data[which.max(treePoints@data$Z), c("X","Y","Z")]
#   #}
#   #cat(sprintf("coordinates of highest point in trees: %.1f s.", difftime(Sys.time(),time.xymax,units = "secs")))
#   #  
#   #allTreeMetrics = cbind(allTreeMetrics, trees.max.xy)
#   #rm(tree.max.xy)
#   
#   #plot(chm)
#   #plot(st_geometry(standMetrics), add=TRUE)
#   #plot(st_geometry(allTreeMetrics), add=TRUE, pch=20)
#   #sum(is.na(standMetrics$treeID))
#   
#   # since segmentation area was buffered beyond stand boundary, clip segmentation to stand
#   # Not strictly necessary (there is duplicate segmentation wherever a buffer projects into another stand) but a
#   # per stand segmentation yielding per stand data files is simple.
#   #coordinates(standAndBufferTrees) = ~Xcentroid + Ycentroid # trees 
#   #standAndBufferTrees = st_as_sf(standAndBufferTrees)
#   #st_crs(standAndBufferTrees) = st_crs(treeMetrics)
#   #standAndBufferTrees = st_transform(standAndBufferTrees, st_crs(standBoundary))
#   #st_crs(standBoundary)
#   standBoundary = st_transform((elliottStands %>% filter(STD_ID == standID))[1], st_crs(standAndBufferTrees)) # essentially a no-op transform: both the stand polygon and trees are EPSG:6557 but st_intersection() sees a mismatch and errors out because a vertical CRS isn't set on the stand (even through st_intersection() is a horizontal operation)
#   standTrees = st_intersection(standAndBufferTrees, standBoundary) # TODO: how to avoid warning: attribute variables are assumed to be spatially constant throughout all geometries?
#   
#   # check plot (not directly useful in larger stands but can be ggsaved() at high resolution)
#   if (plotSegmentation)
#   {
#     chmExtent = ext(segmentation$chm)
#     pointSize = 1.5 / as.numeric(units::set_units(st_area(standBoundary), "acres"))
#     standPlot = ggplot() +
#       geom_spatraster(aes(fill = Z), segmentation$chm) +
#       geom_sf(aes(), allTreeMetrics, color = "cyan2", fill = "transparent") +
#       geom_point(aes(x = X, y = Y), as_tibble(acceptedTreeCentroids), color = "cyan2", shape = 16, size = pointSize) +
#       geom_sf(aes(), standTrees, color = "red", shape = 16, size = pointSize) +
#       guides(fill = "none") +
#       labs(fill = "height, ft") +
#       coord_sf(datum = st_crs(segmentation$chm)) + # geom_spatraster() defaults to showing coordinates in degrees rather than using raster's CRS
#       scale_fill_viridis_c(na.value = "transparent") +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       theme_void()
#     # export to uncompressed .tiff since ggplot() ignores compression settings
#     tiffPath = file.path(dataPath, "segmentation", "604 trees.tiff") # ggsave() doesn't recognize .tiff
#     coordSFaspectRatio = (chmExtent[2] - chmExtent[1]) / (chmExtent[4] - chmExtent[3])
#     ggsave(standPlot, filename = tiffPath, width = 30 * coordSFaspectRatio, height = 30, units = "cm", dpi = 600)
#     # convert .tiff to GeoTIFF
#     standPlotTiff = rast(tiffPath)
#     standPlotPanel1 = ggplot_build(standPlot)$layout$panel_params[[1]]
#     ext(standPlotTiff) = c(standPlotPanel1$x_range, standPlotPanel1$y_range)
#     crs(standPlotTiff) = standPlotPanel1$crs$wkt
#     #standPlotTiff = raster::stack(tiffPath)
#     #standPlotPanel1 = ggplot_build(standPlot)$layout$panel_params[[1]]
#     #raster::extent(standPlotTiff) = c(standPlotPanel1$x_range, standPlotPanel1$y_range)
#     #raster::projection(standPlotTiff) = crs(segmentation$chm)
#     geoTiffPath = file.path(dataPath, "segmentation", "604 trees.tif") # GDAL doesn't recognize .tif
#     writeRaster(standPlotTiff, geoTiffPath, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9", "PHOTOMETRIC=RGB"), datatype = "INT1U", overwrite = TRUE)
#   }
#   
#   #hist(standMetrics$zmax)
#   #hist(standTrees$zmax)
#   #poly.inv.f$StandID = standBoundary$StandID
#   #plot(st_geometry(standAndBufferTrees), add = TRUE, pch = 20)
#   #plot(st_geometry(poly.inv.f), add = TRUE, pch = 20)
#   
#   #allTreeToNNindices = match(standTrees$treeID, treeMetrics$treeID) # unclear if this check remains relevant
#   #if (all(treeMetrics$treeID[allTreeToNNindices] == standTrees$treeID) == FALSE) # check validity of match operation
#   #{
#   #  stop("Index")
#   #}
#   
#   # write segmented trees
#   # TODO; where is best to convert from EPSG:6557 (ft) to 6556 (m)?
#   st_write(standTrees, file.path(dataPath, "segmentation", paste0(standID, " trees.gpkg")), append = FALSE)
#   writeRaster(segmentation$chm, file.path(dataPath, "segmentation", paste0(standID, " chm.tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"))
#   writeLAS(standAndBufferPoints, file.path(dataPath, "segmentation", paste0(standID, " points.laz")))
#   #st_write(allTreeMetrics, file.path(dataPath, "segmentation", paste0(standID, " metrics.gpkg"), append = FALSE) # crowns for the final segmented trees
#   
#   #plot(chm)
#   #plot(st_geometry(standMetrics), add=TRUE)
#   #plot(st_geometry(standTrees), pch=20, add=TRUE)
# }
# 
