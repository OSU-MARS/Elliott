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

segmentationOptions = tibble(includeInvestigatory = FALSE)


# check assignment of stand IDs
# STD_ID: Oregon Department of Forestry-Oregon Department of State Lands stand ID
# OSU_SID: Oregon State University revised stand ID, typically a 10 prefix added to ODF-DSL STD_ID
if (nrow(elliottStands) != length(unique(elliottStands$OSU_SID)))
{
  stop("Stand IDs are not unique to each inventory polygon.")
}


## tile level pixel metrics from lidR, 5950X: get_lidr_threads() = 16
if (segmentationOptions$includeInvestigatory)
{
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
  
  metrics_zero = function(returnNumber, z) # lidR testing function
  {
    return(tibble(n1 = 0, n2 = 0))
  }
  
  metrics_zn = function(returnNumber, z)
  {
    return(data.table(n1 = sum(returnNumber == 1), n2 = sum(returnNumber == 2), n3 = sum(returnNumber == 3), n4 = sum(returnNumber == 4), n5 = sum(returnNumber == 5),
                      z1 = max((returnNumber == 1) * z), z2 = max((returnNumber == 2) * z), z3 = max((returnNumber == 3) * z), z4 = max((returnNumber == 4) * z), z5 = max((returnNumber == 5) * z)))
  }
  
  stdmetrics_z_simple = function(z) # reduced from https://github.com/r-lidar/lidR/blob/master/R/metrics_stdmetrics.R
  {
    n = length(z)
    
    probs = seq(0, 1.0, by = 0.1)
    zq = as.list(stats::quantile(z, probs))
    names(zq) = paste0("zq", 100 * probs)
    
    return(c(n = length(z), zq, zmean = mean(z), zsd = sd(z)))
  }
  
  dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
  tileName = "s04020w06690.laz" # 40-45 year old plantations + 130+ reserve stands
  tileName = "s04260w06720.laz" # 4, 5, 13, 46, 52 year plantations + 95 & 120 year reserves
  tile = readLAS(file.path(dataPath, "pointz", tileName))
  tile = filter_poi(tile, (Classification != 7L & (Classification != 18L)))
  
  metricsZero = pixel_metrics(tile, ~metrics_zero(ReturnNumber, Z), res = 3.28084 * 1) # 1.4 min with data.table(), 5.3 min with tibble() -> 3.8x penalty
  metrics = pixel_metrics(tile, ~stdmetrics(X, Y, Z, Intensity, ReturnNumber, Classification), dz = 3.28084 * 1, res = 3.28084 * 1) # 11.6 min / tile @ 1 m resolution with dz = 1 foot -> 126 MB GeoTIFF
  metricsCtrl = pixel_metrics(tile, .stdmetrics_ctrl, res = 3.28084 * 1) # only n and area but just a few seconds, error with ~stdmetrics_ctrl(X, Y, Z)
  metricsReturns = pixel_metrics(tile, .stdmetrics_rn, res = 3.28084 * 1) # 1.0 min -> 15 MB GeoTIFF
  metricsRgbn = pixel_metrics(tile, ~metrics_rgbn(ReturnNumber, R, G, B, NIR), res = 1) # 20.3 min -> 64 MB GeoTIFF
  metricsZ = pixel_metrics(tile, .stdmetrics_z, dz = 3.28084 * 1, res = 3.28084 * 1) # 5.53 min -> 71 MB GeoTIFF, z mean, σ, quantiles, skew, kurtosis, entropy, probabilities, error with ~stdmetrics_z(X, Y, Z)
  metricsZ = pixel_metrics(tile, ~stdmetrics_z_simple(Z), res = 3.28084 * 1) # 2.5 min -> 34 MB GeoTIFF
  metricsZn = pixel_metrics(tile, ~metrics_zn(ReturnNumber, Z), res = 3.28084 * 1) # 2.4 min -> 16.9 MB GeoTIFF
  
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
}

## remove old versions of .laz files which have been updated
if (segmentationOptions$includeInvestigatory)
{
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
  
  remove_redundant_tiles("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz RGBIR", "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz")
}

## testing
# stand id   2021 tiles                   notes
# 443        s03960w07050, s03990w07050   7.5 ha, 150 y natural regen, N central
# 605        s03990w06930, s04020w06930   40 ha, young plantation, steep, S central
# 1874       s03510w06480, s03540w06480   1.1 ha, natural regen, flat, SW corner
if (segmentationOptions$includeInvestigatory)
{
  # setup
  # 2021 tile naming convention: s<easting>w<northing> where easting ≥ 03450 and northing ≥ 06360 increment in steps of 30 (units are 100 ft)
  dataPath = "E:/Elliott/GIS/DOGAMI/2009 OLC South Coast"
  #dataPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
  
  dtm2021nv5 = rast("GIS/DOGAMI/2021 OLC Coos County/DTM_NV5_2021.tif")
  elliottStands = st_read("GIS/Planning/ESRF_Stands062022Fixed.shp") # June 2022 stand boundaries
  elliottStands = st_transform(elliottStands, 6557) # override EPSG:2992 to 6557 - same projection and units, just different EPSG
  elliottLidarTiles = readLAScatalog(file.path(dataPath, "points"), progress = FALSE) # progress = FALSE turns off plotting of tiles during clip_roi() calls, many warnings: only 2 bytes until point block
  
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
  #segment_cloud = function(pointCloud, isNormalized = FALSE)
  #{
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
  #
  #boundingBoxes = future_map(elliottLidarTiles$filename, function(tilePath)
  #{
  #  boundingBox = merge(st_as_sfc(st_bbox(readLASheader(tilePath))),
  #                      tibble(name = stringr::str_remove(fs::path_file(tilePath), "_las.las")))
  #  return(boundingBox)
  #}, .options = furrr_options(seed = TRUE))
  #
  #st_write(st_sf(bind_rows(boundingBoxes)), file.path(dataPath, "Points", "tile index.gpkg"), append = TRUE)
}

## comparative stdmetrics for Clouds unit test
if (segmentationOptions$includeInvestigatory)
{
  psmeCells = raster::raster(file.path(getwd(), "../tools/Clouds/UnitTests/PSME ABA grid cells.tif"))
  psme146 = readLAS(file.path(getwd(), "../tools/Clouds/UnitTests/PSME LAS 1.4 point type 6.las"))
  psmeMetrics = pixel_metrics(psme146, ~stdmetrics(X, Y, Z, Intensity, ReturnNumber, Classification, dz = 1, th = 920.52 + 2), zmin = min(psme146$Z), psmeCells)
}