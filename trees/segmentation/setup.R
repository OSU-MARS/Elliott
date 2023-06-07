library(dplyr)
library(FNN) # fast nearest neighbor::get.knn() for cleaning segmentation
library(lidR)
library(magrittr)
library(sf)
library(terra)
library(tibble)
library(tidyr)

normalize_stand = function(standID)
{
  standBoundary = (elliottStands %>% filter(STD_ID == standID))[1]
  #standBuffer30m = st_buffer(standBoundary, 30) # m
  #units::set_units(st_area(standBoundary), acres)
  #plot(st_geometry(standBuffer30m), col="green")
  #plot(st_geometry(standBoundary), add=TRUE, col="red")
  
  ## point cloud normalization and canopy height model creation
  segmentationBufferWidth = 30 # ft TODO: is a 30 foot radius is too small for correct segmentation of large trees along the stand boundary?
  standAndBufferPoints = clip_roi(elliottLidarTiles, st_buffer(standBoundary, segmentationBufferWidth))
  #standBuffer30m = st_transform(standBuffer30m, st_crs(standPoints)) # TODO: why clip a second time?
  #standPoints = clip_roi(standPoints, st_buffer(standBoundary, segmentationBufferWidth))
  #standPoints = las_update(standPoints) # update the header
  
  standAndBufferPoints = classify_noise(standAndBufferPoints, ivf(res = 3, n = 1)) #identifies noise using ivf method
  #cat("noise classification...", fill = TRUE)
  #table(standPoints@data$Classification)
  #standPoints.sor=classify_noise(standPoints, sor(treeID=10,m=9)) # identifies noise using sor method
  #table(standPoints@data$Classification,standPoints.sor@data$Classification)
  
  dtm = crop(dtm2021nv5, st_buffer(standBoundary, segmentationBufferWidth + 3)) # crop the large DTM to the point cloud extent, adding additional buffer width so normalize_height() TODO: why?
  standAndBufferPoints = normalize_height(standAndBufferPoints, dtm)
  standAndBufferPoints = filter_poi(standAndBufferPoints, Classification == 1L, Z >= 0, Z <= 400) # exclude ground points and high outliers, Z range in feet
  #plot(decimate_points(standPoints, random(1)))
  #rm(standPoints)
  #gc()
  
  return(standAndBufferPoints)
}

process_stand = function(standID)
{
  ## initial segmentation
  normalizationStartTime = Sys.time()
  standAndBufferPoints = normalize_stand(standID)
  cat(sprintf("%d normalization: %.1f s", standID, difftime(Sys.time(), normalizationStartTime, units = "secs")), fill = TRUE)
  
  segmentationStartTime = Sys.time()
  segmentation = segment_stand(standAndBufferPoints)
  minimumTreeHeight = segmentation$minimumTreeHeight
  standAndBufferPoints = segmentation$standPoints
  rm(segmentation)
  cat(sprintf("%d initial segmentation: %.1f s.", standID, difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
  #length(standPoints@data$treeID[(standPoints@data$treeID =="NA")])
  #sum(is.na(standPoints@data$treeID))
  
  ## look for trees with outlier points above them, remove outliers if found
  highPointStartTime = Sys.time()
  outlierAboveThreshold = 5 # ft
  outlierAboveTrees = standAndBufferPoints@data %>% filter(is.na(treeID) == FALSE) %>% 
    group_by(treeID) %>% 
    reframe(quantiles = c(0.995, 1.0), z = quantile(Z, probs = quantiles, names = FALSE)) %>%
    pivot_wider(id_cols = "treeID", names_prefix = "q", names_from = "quantiles", values_from = "z") %>%
    filter((q1 - q0.995) > outlierAboveThreshold)
  
  if (nrow(outlierAboveTrees) > 0)
  {
    for (treeIndex in 1:nrow(outlierAboveTrees))
    {
      treeID = outlierAboveTrees$treeID[treeIndex]
      zThreshold = outlierAboveTrees$q0.995[treeIndex]
      standAndBufferPoints@data$Classification[which((standAndBufferPoints@data$Z >= zThreshold) & (standAndBufferPoints@data$treeID == treeID))] = 18 # point class: high
    }
    
    standAndBufferPoints = standAndBufferPoints[standAndBufferPoints@data$Classification != 18,] # remove outlier points
    standAndBufferPoints = las_update(standAndBufferPoints)
  }
  cat(sprintf("%d upper quantile checking: found %d trees with large upper quantile differences in %.1f s.", standID, nrow(outlierAboveTrees), difftime(Sys.time(), highPointStartTime, units = "secs")), fill = TRUE)
  
  if (nrow(outlierAboveTrees) > 0)
  {
    segmentationStartTime = Sys.time()
    segmentation = segment_stand(standAndBufferPoints)
    standAndBufferPoints = segmentation$standPoints
    cat(sprintf("%d resegmentation: %.1f s.", standID, difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
  }
  
  #nr.trees.final = length(treeIDs)
  #TPA = nr.trees.final / st_area(standBuffer30m) * 45360
  #cat(sprintf("segmentation time: %.1f s.", difftime(Sys.time(), segmentationStartTime, units = "secs")), fill = TRUE)
  
  allTreeMetrics = crown_metrics(standAndBufferPoints, func = .stdmetrics, geom = "convex") %>% # crown polygons and various statistics including the percentiles
    rename(crownArea = area) # ft², projected
  #allTreeMetrics = st_centroid(allTreeMetrics, of_largest_polygon = FALSE, ensure_within = TRUE) # not needed as crown_metrics() returns a point geometry with centroid coordinates TODO: clear warning st_centroid assumes attributes are constant over geometries
  #treeCentroidXY = st_coordinates(treeCentroids)
  #allTreeMetrics %<>% mutate(Xcentroid = treeCentroidXY[,"X"], Ycentroid = treeCentroidXY[,"Y"])
  
  
  ## remove segmented trees whose tops are less than a threshold distance apart (like 7 ft) or are shorter than minimumTreeHeight
  neighborEliminationStartTime = Sys.time()
  minimumCrownArea = 15 # ft²
  acceptedTreeMetrics = allTreeMetrics %>% filter(crownArea > minimumCrownArea, zmax > minimumTreeHeight)
  acceptedTreeCentroids = st_centroid(st_geometry(acceptedTreeMetrics)) # reduce trees' crown polygons to their centroids (layer's geometry type changes from polygon to point)
  st_geometry(acceptedTreeMetrics) = acceptedTreeCentroids
  #plot(chm)
  #plot(st_geometry(treeMetrics.t), add = TRUE, pch = 20)
  
  # find three nearest neighbors based on distance between trees' centroids
  # $nn.index, which identifies the row index of the tree in the st geometry/data frame, is not the tree number. The matching of index with treeID should be done explicitly.
  acceptedTreeCentroids = st_coordinates(acceptedTreeCentroids) # flatten point geometry to XY matrix
  treesNearestNeighbor3 = get.knn(acceptedTreeCentroids, 
                                  k = 3, 
                                  algorithm = "brute")
  colnames(treesNearestNeighbor3$nn.index) = c("nn1", "nn2", "nn3")
  colnames(treesNearestNeighbor3$nn.dist) = c("nn1distance", "nn2distance", "nn3distance")
  
  neighborMergeDistance = 7 # ft
  standAndBufferTrees = bind_cols(acceptedTreeMetrics,
                                  as_tibble(treesNearestNeighbor3$nn.index), 
                                  as_tibble(treesNearestNeighbor3$nn.dist)) %>% 
    mutate(STD_ID = as.integer(standID),
           treeID = as.integer(treeID),
           species = "DF", # default all trees to Douglas-fir
           xCentroid = acceptedTreeCentroids[,"X"], yCentroid = acceptedTreeCentroids[,"Y"],
           x1 = xCentroid[nn1], x2 = xCentroid[nn2], x3 = xCentroid[nn3], # while they don't need to be stored explicitly, unpack coordinates of each tree's three nearest neighbors for now
           y1 = yCentroid[nn1], y2 = yCentroid[nn2], y3 = yCentroid[nn3],
           z1 = zmax[nn1], z2 = zmax[nn2], z3 = zmax[nn3],
           isShorterThanNearNeighbor = ((nn1distance < neighborMergeDistance) & (zmax < z1)) |
             ((nn2distance < neighborMergeDistance) & (zmax < z2)) |
             ((nn3distance < neighborMergeDistance) & (zmax < z3))) %>%
    relocate(STD_ID, treeID, species, xCentroid, yCentroid, zmax, crownArea)
  standAndBufferTrees %<>% filter(crownArea > minimumCrownArea, # ft² TODO: vary minimum crown projection area with tree height as a 15 ft² fixed area seems high for short trees and small for tall ones?
                                  isShorterThanNearNeighbor == FALSE) %>% 
    select(-isShorterThanNearNeighbor) # for now, update point cloud is not updated merged tree IDs
  # TODO: update crown area, centroid, and so on based on merge of dropped trees?
  #nr.standAndBufferTrees = nrow(standAndBufferTrees)
  #sum(standAndBufferTrees$nn1Dist < neighborMergeDistance)
  #sum(standAndBufferTrees$FP)
  cat(sprintf("%d neighbor elimination time: %.1f s.", standID, difftime(Sys.time(), neighborEliminationStartTime, units = "secs")), fill = TRUE)
  
  ## find x and y coordinate of tallest point within the crown
  # Very time consuming as currently implemented. TODO: is group_by() and summarize(Z == acceptedTreeMetrics$zmax) more performant?
  #trees.max.xy = as.data.frame(matrix(0, length(standMetrics$treeID), 3))
  #names(trees.max.xy) = c("Xmax","Ymax", "Zmax")
  #  
  #time.xymax = Sys.time()
  #for (kk in 1:length(standMetrics$treeID)) 
  #{
  #  treePoints = filter_poi(standPoints, treeID == standMetrics$treeID[kk])
  #  trees.max.xy[kk,] = treePoints@data[which.max(treePoints@data$Z), c("X","Y","Z")]
  #}
  #cat(sprintf("coordinates of highest point in trees: %.1f s.", difftime(Sys.time(),time.xymax,units = "secs")))
  #  
  #allTreeMetrics = cbind(allTreeMetrics, trees.max.xy)
  #rm(tree.max.xy)
  
  #plot(chm)
  #plot(st_geometry(standMetrics), add=TRUE)
  #plot(st_geometry(allTreeMetrics), add=TRUE, pch=20)
  #sum(is.na(standMetrics$treeID))
  
  # since segmentation area was buffered beyond stand boundary, clip segmentation to stand
  # Not strictly necessary (there is duplicate segmentation wherever a buffer projects into another stand) but a
  # per stand segmentation yielding per stand data files is simple.
  #coordinates(standAndBufferTrees) = ~Xcentroid + Ycentroid # trees 
  #standAndBufferTrees = st_as_sf(standAndBufferTrees)
  #st_crs(standAndBufferTrees) = st_crs(treeMetrics)
  #standAndBufferTrees = st_transform(standAndBufferTrees, st_crs(standBoundary))
  #st_crs(standBoundary)
  standBoundary = st_transform((elliottStands %>% filter(STD_ID == standID))[1], st_crs(standAndBufferTrees)) # essentially a no-op transform: both the stand polygon and trees are EPSG:6557 but st_intersection() sees a mismatch and errors out because a vertical CRS isn't set on the stand (even through st_intersection() is a horizontal operation)
  standTrees = st_intersection(standAndBufferTrees, standBoundary) # TODO: how to avoid warning: attribute variables are assumed to be spatially constant throughout all geometries?
  
  #hist(standMetrics$zmax)
  #hist(standTrees$zmax)
  #poly.inv.f$StandID = standBoundary$StandID
  #plot(st_geometry(standAndBufferTrees), add = TRUE, pch = 20)
  #plot(st_geometry(poly.inv.f), add = TRUE, pch = 20)
  
  #allTreeToNNindices = match(standTrees$treeID, treeMetrics$treeID) # unclear if this check remains relevant
  #if (all(treeMetrics$treeID[allTreeToNNindices] == standTrees$treeID) == FALSE) # check validity of match operation
  #{
  #  stop("Index")
  #}
  
  # write segmented trees
  # TODO; where is best to convert from EPSG:6557 (ft) to 6556 (m)?
  st_write(standTrees, file.path(dataPath, "segmentation", paste0(standID, " trees.gpkg")), append = FALSE)
  #st_write(allTreeMetrics, file.path(dataPath, "segmentation", paste0(standID, " metrics.gpkg"), append = FALSE) # crowns for the final segmented trees
  
  #plot(chm)
  #plot(st_geometry(standMetrics), add=TRUE)
  #plot(st_geometry(standTrees), pch=20, add=TRUE)
}

segment_stand = function(standAndBufferPoints)
{
  chm = rasterize_canopy(standAndBufferPoints, res = 1, algorithm = p2r(subcircle = 0)) # compute the canopy height model using a pixel of 0.5 units, in this case meters
  #plot(st_geometry(standBuffer30m))
  #plot(chm,add=TRUE)
  
  chmHeightQuantiles = quantile(values(chm), probs = c(0.1, 0.85, 0.95), na.rm = TRUE, names = FALSE) # point density dependent: 0.1, 0.65, 0.95? for 2009 flight?
  names(chmHeightQuantiles) = c("minimumTreeHeight", "q85", "q95")

  get_window_size <- function(x)
  { 
    return(x * log((chmHeightQuantiles[3] - chmHeightQuantiles[2])) / 25 + (chmHeightQuantiles[3] - chmHeightQuantiles[2]) / 10)
  }
  localMaxima = find_trees(chm, lmf(get_window_size, hmin = chmHeightQuantiles[1]), uniqueness = "incremental") # identifies the tree tops from the CHM using the local maximum filter with a window of 6 feet and height >8 ft
  #localMaxima = find_trees(chm,lmf(chmHeightQuantiles[3]/25,hmin=chmHeightQuantiles[1]),uniqueness = "incremental")
  #dim(localMaxima)
  #plot(localMaxima, pch=20, add=TRUE)
  
  minimumTreeHeight = chmHeightQuantiles["minimumTreeHeight"] # ft
  maxCrownDiameter = chmHeightQuantiles["q95"] / 8 # ft
  segmentationAlgorithm = dalponte2016(chm, localMaxima, th_tree = minimumTreeHeight, max_cr = maxCrownDiameter) #create an R object that stores the algorithm to be used for tree segmentation
  #segmentationAlgorithm = watershed(chm, th_tree=15)
  #segmentationAlgorithm = li2012(dt1=1, dt2=2, Zu=30, hmin=10, speed_up = 8)
  #segmentationAlgorithm = silva2016(chm,localMaxima, max_cr_factor = 0.5, exclusion = 0.3, ID="treeID")
  
  # TODO: how to avoid  GDAL message: Setting GeoTIFF 1.1 requires libgeotiff >= 1.6?
  standAndBufferPoints = segment_trees(standAndBufferPoints, algorithm = segmentationAlgorithm) # segment point cloud: identifies individual tree crowns
  return(list(standPoints = standAndBufferPoints, minimumTreeHeight = minimumTreeHeight))
}

# setup
#dataPath = "E:/Elliott/GIS/DOGAMI/2009 OLC South Coast"
dataPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"

dtm2021nv5 = rast("GIS/NV5/DTM_NV5_2021.tif")
elliottStands = st_read("GIS/Planning/ESRF_Stands062022Fixed.shp") # June 2022 stand boundaries
elliottStands = st_transform(elliottStands, 6557) # override EPSG:2992 to 6557 - same projection and units, just different EPSG
elliottLidarTiles = readLAScatalog(file.path(dataPath, "Points"), progress = FALSE) # progress = FALSE turns off plotting of tiles during clip_roi() calls, many warnings: only 2 bytes until point block

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


## testing
#process_stand(1874) # one among many small stands for testing

#library(furrr)
#options(future.debug = TRUE)
#plan(multisession, workers = 2)
#
#standIDs = c(1874, 1875)
#future_map(standIDs, function(standID)
#{
#  process_stand(standID)
#})
