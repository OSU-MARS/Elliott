# some investigatory code here assumes segmentation/setup.R
library(caret)
library(dplyr)
library(FNN)
library(purrr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(progressr)
library(ranger)
library(rsample)
library(sf)
library(stringr)
library(terra)
library(tidyr)

options(cli.progress_format_iterator = "{cli::pb_bar} {cli::pb_percent} | cross validation fold {cli::pb_current} of {cli::pb_total}, {cli::pb_elapsed_clock} elapsed, {prettyunits::pretty_sec(as.numeric(cli::pb_eta_raw))} remaining")
theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             panel.border = element_blank(), 
                             plot.title = element_text(size = 10)))

# ranger performance maxima: 
# Zen 3 + DDR4-3200: one thread per core (default of two threads per core is slower and bogs the UX)
# Zen 5 + DDR5-5600: TBD
treetopOptions = tibble(fitRandomForest = FALSE,
                        folds = 2,
                        repetitions = 2,
                        neighborhoodBufferWidth = 25, # in CRS units, so feet
                        rangerThreads = 0.5 * future::availableCores(),
                        dsmCellSize = 1.5, # feet
                        tileSize = 3000, # ft
                        verticalExaggeration = 25, # DSM multiplier
                        includeInvestigatory = FALSE,
                        includeSetup = FALSE)

localMaximaPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta/local maxima"
acceptedTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted"
candidateTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops"
tileCrs = NULL

fit_ranger_treetop = function(trainingMaxima, neighborhoodMaxima = trainingMaxima, mtry, minNodeSize, sampleFraction, classWeights = NULL, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
{
  progressBar = progressor(steps = folds * repetitions)
  
  if ((folds == 1) & (repetitions == 1))
  {
    start = Sys.time()
    allFit = ranger(treetop ~ ., data = trainingMaxima %>% select(-all_of(accuracyVariables)), mtry = mtry, min.node.size = minNodeSize, sample.fraction = sampleFraction, class.weights = classWeights, num.threads = treetopOptions$rangerThreads)
    allFitPredicted = predict(allFit, trainingMaxima %>% select(-treetop))$predictions

    tileMergePoints = get_merge_points(trainingMaxima, neighborhoodMaxima)
    allFitPredicted[tileMergePoints$mergePointIndices] = "merge"
    allFitPredicted[tileMergePoints$treetopIndices] = "yes"
    
    allFitMetrics = get_treetop_accuracy(allFitPredicted, trainingMaxima)
    allFitMetrics$fitTimeInS = as.numeric(difftime(Sys.time(), start, units = "secs"))
    progressBar()
    return(allFitMetrics %>% mutate(repetition = 1, fold = 1))
  }
  
  # dataFold = splitsAndFits$splits[[1]]
  fitFunction = function(dataFold)
  {
    start = Sys.time()
    foldTrainingData = analysis(dataFold)
    fit = ranger(treetop ~ ., data = foldTrainingData %>% select(-all_of(accuracyVariables)), mtry = mtry, min.node.size = minNodeSize, sample.fraction = sampleFraction, class.weights = classWeights, num.threads = treetopOptions$rangerThreads)
    validationData = assessment(dataFold)
    fitPredicted = predict(fit, validationData %>% select(-treetop))$predictions

    tileMergePoints = get_merge_points(validationData, neighborhoodMaxima, treetopClassification = fitPredicted)
    fitPredicted[tileMergePoints$mergePointIndices] = "merge"
    fitPredicted[tileMergePoints$treetopIndices] = "yes"
    fitMetrics = get_treetop_accuracy(fitPredicted, validationData)
    fitMetrics$fitTimeInS = as.numeric(difftime(Sys.time(), start, units = "secs"))
    progressBar()
    return(fitMetrics)
  }
  
  # use random cross validation as the training dataset lacks the spatial extent to block by stands
  # ranger is expected to use all cores, so little advantage to future_map() rather than map() is assumed
  splitsAndFits = vfold_cv(trainingMaxima, v = folds, repeats = repetitions) %>% 
    mutate(fit = map(splits, fitFunction)) %>% 
    select(-splits) %>%
    unnest(fit)
  if ("id2" %in% names(splitsAndFits))
  {
    splitsAndFits %<>% rename(repetition = id, fold = id2) %>% mutate(repetition = as.integer(str_remove(repetition, "Repeat")),
                                                                      fold = as.integer(str_remove(fold, "Fold")))
  } else { # 
    splitsAndFits %<>% rename(fold = id) %>% mutate(repetition = 1,
                                                    fold = as.integer(str_remove(fold, "Fold"))) %>%
      relocate(repetition, fold)
  }
  return(splitsAndFits)
}

get_merge_points = function(tileMaxima, neighborhoodMaxima, treetopClassification = tileMaxima$treetop)
{
  ## unsupervised postprocessing of local maxima classified by random forest to generate and resolve merge clusters
  # Merge clusters are groups of local maxima believed to correspond a single treetop. LiDAR point clouds yield clusters of local maxima
  # treetop candidates when
  # - Hits at identical heights occur in adjacent cells, either due to a treetop lying near the cell boundary or the LiDAR spot size being
  #   comparable to the cell size.
  # - Flat tops present multiple maxima whose elevations are not distinguishable within the flight data's vertical accuracy. This may be
  #   due to habit, branch architecture, and leaf shape (multistemmed hardwoods, for example) or suppression (likely made visible by overstory
  #   mortality).
  # - Trees movement between flight lines due to wind sway, resulting in top displacements of two or more cells.
  # - Misalignment of flight lines causes a tree's top to appear in multiple locations within the merged data.
  # - Broken tops exposing a whorl of branches as a circle of local maxima, possibly with a central maxima at the main stem.
  # Ideally, processing errors which produce copies of trees' upper portions are rejected as noise. These are difficult to identify in general,
  # however, and more difficult to locate using the reduced set of information available from a digital surface model compared to a point 
  # cloud.
  #
  # In postprocessing of classified local maxima, merge clusters occur where
  # - A local maxima is classified as a merge point. If other maxima are close enough a cluster is formed, possibly including other merge
  #   points or linking with additional local maxima captured by adjacent merge points. This case captures errors where the classifier puts up
  #   a merge point but mislabels other points which should be merged as treetops or incorrectly discards local maxima.
  # - Two or more nearby maxima classified as treetops are close enough to each other to be deemed a cluster. This case acts to reduce 
  #   oversegmentation by capturing classification errors where the classifier should have generated a merge pair, triplet, or larger group.
  # 
  # The approach used here is
  # - Reclassify treetops which pair with nearby treetops as candidate merge points. (Triples and larger groups are possible but unlikely.)
  # - Find sufficiently close neighbors of all merge points, whether designated by classification or obtained from treetop pairing. 
  # - Add any neighbors not already present in the merge point set and establish their neighbors within the set. Neighbors outside of the
  #   set are excluded to prevent further cluster growth.
  # - Extract merge clusters from the neighbor connectivity graph and find the clusters' xyz centroids. Clusters typically have limited
  #   variation in height (z), either due to flat topping or limited changes in treetop elevation with wind sway and misalingment between
  #   different flight lines' point cloud strips.
  # - Identify and remove outliers from clusters. Clusters with a single merge point are reverted to treetops, converted to treetops if 
  #   initially labeled as merge points, or dropped if not initially labeled as a treetop or merge point.
  # - Yield treetop positions as the xyz averages of clusters with two or more merge points.
  #
  # Merge clustering around the initial classification is controlled primarily by the definition of neighbor distances, both for treetop
  # pairing and cluster inclusion. Because neighbors are found by FNN's implementation of k-nearest neighbors (kNN), Euclidean distances 
  # are used. And, because clusters' vertical extent tends to be limited, distinctions between two- (xy) and and three-dimensional (xyz) 
  # neighborhoods are constrained. Currently, kNN searches are done in two dimensions and decisions on whether to include kNN identified
  # neighbors in a merge cluster are based on separate consideration of points' xyz positions relative to the cluster's centroid. Since
  # merge point clustering is presumed needed to seed crown segmentation, only the candidate cluster points and geometric information
  # passed through from the digital surface model and point clouds is available for decision making.

  # find cluster convertable treetop pairs (triplets, ...)
  treetopTileIndices = which(treetopClassification == "yes")
  treetopPoints = tileMaxima[treetopTileIndices, ]
  neighborhoodTreetops = neighborhoodMaxima %>% filter(treetop == "yes")
  treetopMergeKnn = get.knnx(neighborhoodTreetops[, c("x", "y")], treetopPoints[, c("x", "y")], k = 5)
  treetopMergeKnn = tibble(tile = treetopPoints$tile, id = treetopPoints$id, uniqueID = treetopPoints$uniqueID, sourceID = treetopPoints$sourceID,
                           treetop = treetopPoints$treetop, x = treetopPoints$x, y = treetopPoints$y,
                           radius = treetopPoints$radius, dsmZ = treetopPoints$dsmZ, cmmZ = treetopPoints$cmmZ, height = treetopPoints$height,
                           neighborDistanceThreshold = get_merge_distance(treetopPoints, treetopMerge = TRUE),
                           neighborhoodIndex = treetopMergeKnn$nn.index[, 1], # nn.dist[, 1] is self since get.knnx(neighborhood, tile) is an overlapping query
                           neighbor1distance = treetopMergeKnn$nn.dist[, 2],
                           neighbor2distance = treetopMergeKnn$nn.dist[, 3],
                           neighbor3distance = treetopMergeKnn$nn.dist[, 4],
                           neighbor4distance = treetopMergeKnn$nn.dist[, 5],
                           neighbors = (neighbor1distance <= neighborDistanceThreshold) + (neighbor2distance <= neighborDistanceThreshold) + (neighbor3distance <= neighborDistanceThreshold) + (neighbor4distance <= neighborDistanceThreshold),
                           # gather IDs of neighbors
                           neighbor1uniqueID = if_else(neighbor1distance <= neighborDistanceThreshold, neighborhoodTreetops$uniqueID[treetopMergeKnn$nn.index[, 2]], NA_integer_), # NAs propagate through indexing, so IDs are NA for neighbors beyond mergeable range
                           neighbor2uniqueID = if_else(neighbor2distance <= neighborDistanceThreshold, neighborhoodTreetops$uniqueID[treetopMergeKnn$nn.index[, 3]], NA_integer_),
                           neighbor3uniqueID = if_else(neighbor3distance <= neighborDistanceThreshold, neighborhoodTreetops$uniqueID[treetopMergeKnn$nn.index[, 4]], NA_integer_),
                           neighbor4uniqueID = if_else(neighbor4distance <= neighborDistanceThreshold, neighborhoodTreetops$uniqueID[treetopMergeKnn$nn.index[, 5]], NA_integer_),
                           # find indices of neighbors
                           neighbor1neighborhoodIndex = match(neighbor1uniqueID, neighborhoodTreetops$uniqueID),
                           neighbor2neighborhoodIndex = match(neighbor2uniqueID, neighborhoodTreetops$uniqueID),
                           neighbor3neighborhoodIndex = match(neighbor3uniqueID, neighborhoodTreetops$uniqueID),
                           neighbor4neighborhoodIndex = match(neighbor4uniqueID, neighborhoodTreetops$uniqueID))
  #table(treetopMergeKnn$neighbors) # likely ~99.5% singleton treetops which don't contribute merge points
  #treetopMergeKnn %>% filter(id %in% c(24704, 25078, 24432, 25893, 26179)) %>% mutate(neighbor1uniqueID = neighbor1uniqueID - 423006840000000) %>% select(tile, id, neighbors, neighborDistanceThreshold, neighbor1uniqueID, neighbor1distance)
  
  # find merge point neighborhoods
  mergePointTileIndices = which(treetopClassification == "merge")
  if (length(mergePointTileIndices) > 0)
  {
    mergePoints = tileMaxima[mergePointTileIndices, ]
    mergePointKnn = get.knnx(neighborhoodMaxima[, c("x", "y")], mergePoints[, c("x", "y")], k = 5) # could also use cutree(hclust(dist(xy))) but that's O(N²)
    #range(mergePointKnn$nn.dist[, 1]) # should be [0, 0] as all merge points should be neighborhoodMaxima
    mergePointKnn = tibble(tile = mergePoints$tile, id = mergePoints$id, uniqueID = mergePoints$uniqueID, sourceID = mergePoints$sourceID,
                           treetop = mergePoints$treetop, x = mergePoints$x, y = mergePoints$y,
                           radius = mergePoints$radius, dsmZ = mergePoints$dsmZ, cmmZ = mergePoints$cmmZ, height = mergePoints$height,
                           mergeDistanceThreshold = get_merge_distance(mergePoints), # manual tuning from tile review
                           neighborhoodIndex = mergePointKnn$nn.index[, 1], # nn.dist[, 1] is self since get.knnx(neighborhood, tile) is an overlapping query
                           neighbor1distance = mergePointKnn$nn.dist[, 2],
                           neighbor2distance = mergePointKnn$nn.dist[, 3],
                           neighbor3distance = mergePointKnn$nn.dist[, 4],
                           neighbor4distance = mergePointKnn$nn.dist[, 5],
                           # indices of mergeable neighbors within neighborhoodMaxima
                           # TODO: should distance threshold be adjusted as a function of a neighbor's height relative to a merge point? 
                           #       whether a neighbor has the same or different source ID?
                           #       with neighborhood rugosity?
                           #       should lower members of merge clusters sometimes be removed?
                           neighbor1neighborhoodIndex = if_else(neighbor1distance <= mergeDistanceThreshold, mergePointKnn$nn.index[, 2], NA_integer_),
                           neighbor2neighborhoodIndex = if_else(neighbor2distance <= mergeDistanceThreshold, mergePointKnn$nn.index[, 3], NA_integer_),
                           neighbor3neighborhoodIndex = if_else(neighbor3distance <= mergeDistanceThreshold, mergePointKnn$nn.index[, 4], NA_integer_),
                           neighbor4neighborhoodIndex = if_else(neighbor4distance <= mergeDistanceThreshold, mergePointKnn$nn.index[, 5], NA_integer_),
                           neighbors = (neighbor1distance < mergeDistanceThreshold) + (neighbor2distance < mergeDistanceThreshold) + (neighbor3distance < mergeDistanceThreshold) + (neighbor4distance < mergeDistanceThreshold),
                           # gather IDs of neighbors
                           neighbor1uniqueID = neighborhoodMaxima$uniqueID[neighbor1neighborhoodIndex], # NAs propagate through indexing, so IDs are NA for neighbors beyond mregeable range
                           neighbor2uniqueID = neighborhoodMaxima$uniqueID[neighbor2neighborhoodIndex],
                           neighbor3uniqueID = neighborhoodMaxima$uniqueID[neighbor3neighborhoodIndex],
                           neighbor4uniqueID = neighborhoodMaxima$uniqueID[neighbor4neighborhoodIndex])
    #table(mergePointKnn$neighbors)
    #mergePointKnn %>% filter(id %in% c(59103, 59104, 59494)) %>% mutate(neighbor1uniqueID = neighbor1uniqueID - 423006840000000, neighbor2uniqueID = neighbor2uniqueID - 423006840000000) %>% select(tile, id, neighbors, mergeDistanceThreshold, neighbor1uniqueID, neighbor1distance, neighbor2uniqueID, neighbor2distance)
  } else {
    mergePointKnn = tibble(uniqueID = numeric(), neighbors = integer())
  }
  
  # find merge clusters
  # Most merge clusters are simple in the sense they consist of a pair (~75%) or triplet (~20%) of points with matching sets of neighbors.
  # More complex clusters regularly occur where merge points capture a nearby point which, in turn, connects to one or more additional points
  # which may connect to further points. Assembling # these clusters requires traversal of their connectivity graph, which is done here with 
  # the Matrix and igraph packages (code modified substantially from https://stackoverflow.com/questions/47322126/merging-list-with-common-elements).
  mergeClusterInitiatingPoints = bind_rows(mergePointKnn %>% filter(neighbors > 0),
                                           treetopMergeKnn %>% filter(neighbors > 0))
  mergePointUniqueIDs = unique(c(mergeClusterInitiatingPoints$uniqueID,
                                 na.omit(mergeClusterInitiatingPoints$neighbor1uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor2uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor3uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor4uniqueID)))
  addedMergePointIDs = setdiff(mergePointUniqueIDs, mergeClusterInitiatingPoints$uniqueID)

  mergePoints = bind_rows(mergeClusterInitiatingPoints %>% mutate(isInitiating = TRUE), 
                          neighborhoodMaxima %>% filter(uniqueID %in% addedMergePointIDs) %>% mutate(isInitiating = FALSE))

  adjacencyMatrixFrame = pivot_longer(mergePoints %>% select(uniqueID, neighbor1uniqueID, neighbor2uniqueID, neighbor3uniqueID, neighbor4uniqueID) %>% mutate(neighbor0uniqueID = uniqueID) %>% rename(mergePointUniqueID = uniqueID),
                                      cols = c("neighbor0uniqueID", "neighbor1uniqueID", "neighbor2uniqueID", "neighbor3uniqueID", "neighbor4uniqueID"),
                                      names_pattern = "neighbor(.)uniqueID", names_to = "neighbor",
                                      values_to = "clusterMemberUniqueID") %>%
    filter(is.na(clusterMemberUniqueID) == FALSE) %>%
    mutate(row = match(mergePointUniqueID, mergePoints$uniqueID), column = match(clusterMemberUniqueID, mergePoints$uniqueID))
  adjacencyMatrix = Matrix::sparseMatrix(i = adjacencyMatrixFrame$row, j = adjacencyMatrixFrame$column, x = TRUE, dimnames = list(mergePoints$uniqueID, mergePoints$uniqueID))
  connectivityMatrix = Matrix::tcrossprod(adjacencyMatrix, boolArith = TRUE) > 0
  mergeClusterComponents = igraph::components(igraph::graph_from_adjacency_matrix(connectivityMatrix))
  
  mergePoints = left_join(mergePoints,
                          tibble(uniqueID = as.numeric(names(mergeClusterComponents$membership)), mergeClusterNumber = as.integer(mergeClusterComponents$membership)),
                          by = join_by("uniqueID")) %>%
    group_by(mergeClusterNumber) %>%
    # calculate merge cluster geometry
    # Not currently weighted by DSM z within each source ID, though doing so would likely increase accuracy in merging upper elevation wind sway clusters.
    mutate(clusterCentroidX = mean(x),
           clusterCentroidY = mean(y), 
           clusterCentroidDsmZ = mean(dsmZ),
           clusterDistance = sqrt((x - clusterCentroidX)^2 + (y - clusterCentroidY)^2 + (dsmZ - clusterCentroidDsmZ)^2),
           clusterMaxDsmZ = max(dsmZ),
           clusterMaxOriginatingDsmZ = max(if_else(isInitiating, dsmZ, -Inf)),
           # ±14 cm vegetation accuracy in feet + height based wind motion tolerance, manual height coefficient tune from tile inspection
           clusterThickness = 2 * 0.46 + 0.008 * max(height),
           # enable rejection of non-initiating points significantly above the initiating points
           # Makes clusters less likely to climb branches of adjacent trees.
           clusterMaxDsmZ = if_else((clusterMaxDsmZ - clusterThickness) > clusterMaxOriginatingDsmZ, clusterMaxOriginatingDsmZ, clusterMaxDsmZ),
           clusterMeanDistance = mean(clusterDistance)) %>%
    # TODO: refine outlier removal
    # If a noise point's included in a cluster, typically due to not being classified as such, .
    filter((dsmZ >= (clusterMaxDsmZ - clusterThickness)) & (dsmZ <= clusterMaxDsmZ), # within height range
           clusterDistance < 2 * clusterMeanDistance) %>% # within distance of centroid
    # recalculate merge cluster geometry after outlier removal
    mutate(clusterCentroidX = mean(x),
           clusterCentroidY = mean(y), 
           clusterCentroidDsmZ = mean(dsmZ),
           clusterDistance = sqrt((x - clusterCentroidX)^2 + (y - clusterCentroidY)^2 + (dsmZ - clusterCentroidDsmZ)^2),
           clusterMeanDistance = mean(clusterDistance),
           clusterPoints = n())
    # leave grouped
    # likely contains some clusters with one point due to all neighbors being rejected as cluster suitable
  #table(mergePoints$clusterPoints)
  #mergePoints %>% filter(clusterPoints == 1) %>% select(id, treetop, isInitiating, neighbors, neighbor1distance, neighbor2distance, neighbor3distance)
  #mergePoints %>% ungroup() %>% summarize(mergePoints = n(), clusterlessPoints = sum(is.na(mergeClusterNumber))) # likely fewer points than the initial set due to outlier removal
  #mergePoints %>% filter(id %in% c(24704, 25078)) %>% mutate(neighbor1uniqueID = neighbor1uniqueID - 423006840000000) %>% select(mergeClusterNumber, clusterCentroidX, clusterCentroidY, tile, id, x, y, neighborDistanceThreshold, neighbor1uniqueID, neighbor1distance, neighbor2uniqueID, neighbor2distance)
  #mergePoints %>% filter(mergeClusterNumber == 638) %>% mutate(neighbor1uniqueID = neighbor1uniqueID - 423006840000000) %>% select(mergeClusterNumber, clusterCentroidX, clusterCentroidY, tile, id, x, y, neighborDistanceThreshold, neighbor1uniqueID, neighbor1distance, neighbor2uniqueID, neighbor2distance)
  
  # summarize treetops out of merge clusters
  # Treetop IDs are ID of the lowest, on tile local maxima present in the merge cluster. Including off tile maxima potentially yields non-unique 
  # treetop IDs within tiles as taking an off tile ID could collide with another treetop on the tile.
  # Treetops which lie in other tiles are excluded by checking against tile extents to maintain tile boundaries. Processing of the adjacent tiles 
  # will also find these merge clusters.
  tileNames = unique(tileMaxima$tile)
  tileGridIndices = get_tile_grid_indices(tileNames)
  mergePointsInOnTileClusters = left_join(mergePoints %>% # still grouped by mergeClusterNumber
                                            filter(clusterPoints > 1) %>% # exclude singletons
                                            mutate(tileGridIndexX = floor(clusterCentroidX / treetopOptions$tileSize),
                                                   tileGridIndexY = floor(clusterCentroidY / treetopOptions$tileSize)),
                                            tileGridIndices %>% rename(centroidTile = tile),
                                          by = join_by(tileGridIndexX, tileGridIndexY)) %>%
                                  filter(is.na(tile) == FALSE)
  treetopsFromMergePoints = mergePointsInOnTileClusters %>% # still grouped by mergeClusterNumber
    filter(clusterPoints > 1) %>% # exclude singletons
    summarize(id = min(if_else(tile == centroidTile[1], id, NA_integer_), na.rm = TRUE), # all clusters should have at least one merge point in tileMaxima
              tile = centroidTile[1],
              sourceID = if_else(length(unique(sourceID)) == 1, sourceID[1], as.integer(0)), 
              x = clusterCentroidX[1],
              y = clusterCentroidY[1],
              radius = max(radius), dsmZ = clusterCentroidDsmZ[1], cmmZ = mean(cmmZ), height = mean(height),
              maxima = n(),
              # leave ring statistics unpopulated as they're not well defined for treetops obtained from merge clusters
              sourceIDs = length(unique(sourceID)), 
              mergeClusterNumber = mergeClusterNumber[1],
              .groups = "drop")

  # identify merge points which need reclassification as treetops because they're in single point clusters
  # Treetops with neighbors in treetopMergeKnn but which do not aggregate into multi-point clusters do not need to be listed here
  # as they're already classified as treetops.
  # Singleton merge points which do not appear to merit conversion to treetops could be rejected here.
  singletonMergePointIDs = unique(c((mergePointKnn %>% filter(neighbors == 0))$uniqueID, # initiating merge points without any neighbors within mergable range
                                    (mergePoints %>% filter(tile %in% tileNames, clusterPoints == 1, treetop == "merge"))$uniqueID)) # initiating merge points with in range neighbors rejected by cluster filtering
  singletonMergePointTileIndices = match(singletonMergePointIDs, tileMaxima$uniqueID)
  #mergePointKnn %>% filter(neighbors == 0, id %in% c(24704, 25078, 24432, 25893, 26179))
  #sum(is.na(singletonMergePointTileIndices))
  
  # identify on tile local maxima which are merge points so that caller can reclassify any non-merge points that were clustered in
  # Contains NAs for merge points which are in neighborhoodMaxima but not in tileMaxima.
  mergePointTileIndices = match((mergePoints %>% filter(tile %in% tileNames, clusterPoints > 1))$uniqueID, tileMaxima$uniqueID)
  #sum(is.na(mergePointTileIndices))
  
  return(list(treetops = treetopsFromMergePoints, mergePointIndices = mergePointTileIndices, treetopIndices = singletonMergePointTileIndices))
}

get_merge_distance = function(tileMaxima, treetopMerge = FALSE)
{
  # TODO: include radius penalty for same source ID?
  return(treetopOptions$dsmCellSize * (1 + 1/60 * tileMaxima$height + if_else(tileMaxima$height < 125, 0, 1/40 * (tileMaxima$height - 125))))
  #if (treetopMerge)
  #{
  #  return(treetopOptions$dsmCellSize * (1 + 1/60 * tileMaxima$height + if_else(tileMaxima$height < 125, 0, 1/40 * (tileMaxima$height - 125))))
  #} else {
  #  return(treetopOptions$dsmCellSize + 0.1 * treetopOptions$verticalExaggeration + tileMaxima$height^0.45)
  #}
  
  #inclusionRadius = tibble(height = seq(0, 250), 
  #                         linear = treetopOptions$dsmCellSize * (1 + 1/60 * height + if_else(height < 125, 0, 1/40 * (height - 125))),
  #                         power = height^0.5)
  #ggplot() +
  #  geom_line(aes(x = linear, y = height, color = "linear"), inclusionRadius) +
  #  geom_line(aes(x = power, y = height, color = "power"), inclusionRadius) +
  #  labs(x = "radius, feet", y = "height, feet", color = NULL) +
  #  scale_x_continuous(breaks = seq(0, 20, by = 3.28084), labels = scales::number_format(accuracy = 0.01))
}

get_tile_extent = function(tileName, bufferWidth = 0)
{
  xTileMin = 100 * as.integer(str_sub(tileName, 2, 6))
  yTileMin = 100 * as.integer(str_sub(tileName, 8, 12))
  return(list(xMin = xTileMin - bufferWidth,
              xMax = xTileMin + treetopOptions$tileSize + bufferWidth,
              yMin = yTileMin - bufferWidth,
              yMax = yTileMin + treetopOptions$tileSize + bufferWidth))
}

get_tile_grid_indices = function(tileNames)
{
  return(tibble(tile = tileNames,
                tileGridIndexX = 100 * as.integer(str_sub(tileNames, 2, 6)) / treetopOptions$tileSize,
                tileGridIndexY = 100 * as.integer(str_sub(tileNames, 8, 12)) / treetopOptions$tileSize))
}

get_treetop_accuracy = function(predicted, validationData)
{
  predictedBinary = forcats::fct_collapse(predicted, no = c("no", "noise", "maybe noise"), yes = c("yes", "merge"))
  expectedBinary = forcats::fct_collapse(validationData$treetop, no = c("no", "noise", "maybe noise"), yes = c("yes", "merge"))
  confusionMatrix = confusionMatrix(predictedBinary, expectedBinary, positive = "yes", mode = "everything")
  confusionSubmatrix = confusionMatrix(predicted, validationData$treetop)
  
  classAccuracy = diag(confusionMatrix$table) / rowSums(confusionMatrix$table)
  subclassAccuracy = diag(confusionSubmatrix$table) / rowSums(confusionSubmatrix$table)
  classCounts = table(predicted)
  expectedClassCounts = table(validationData$treetop)
  
  # can't readily compute distance and height error on mismatches here
  # - with random k-fold cross validation the status, on average, of s/k nearby maxima is unknown (where s is the tile sampling fraction)
  # - with spatially blocked cross validation the window to use to exclude edge effects on the validation data must be known
  return(tibble(n = nrow(validationData),
                treetop = classCounts["yes"],
                merge = classCounts["merge"],
                noise = classCounts["noise"] + classCounts["maybe noise"],
                overallAccuracy = confusionMatrix$overall["Accuracy"],
                treetopAccuracy = classAccuracy["yes"], 
                nonTreetopAccuracy = classAccuracy["no"], 
                overallSubaccuracy = confusionSubmatrix$overall["Accuracy"],
                expectedTreetop = expectedClassCounts["yes"],
                expectedMerge = expectedClassCounts["merge"],
                expectedNoise = expectedClassCounts["noise"] + expectedClassCounts["maybe noise"],
                subclassAccuracy = list(subclassAccuracy),
                confusionMatrix = list(confusionMatrix), 
                confusionSubmatrix = list(confusionSubmatrix)))
}

get_treetop_eligible_maxima = function(tileName, minimumHeight = 3.28084, acceptedTileName = NULL) # default to minimum height of 1 m, but in feet
{
  localMaxima = vect(file.path(localMaximaPath, paste0(tileName, ".gpkg")))
  if (is.null(tileCrs))
  {
    tileCrs <<- crs(localMaxima)
  } else {
    if (same.crs(tileCrs, crs(localMaxima)) == FALSE)
    {
      stop(paste0("Tile ", localMaximaPath, " has crs ", crs(localMaxima), " rather than ", tileCrs))
    }
  }
  # use minimumHeight to exclude groundcover maxima and maxima from sensor noise or error
  # for now, also exclude treetop candidates with so few adjacent points ring 1 or ring 2 has no data since it's 1) it's difficult to tell if these are maxima, 2) it's unlikely they're actually maxima, and 3) these likely comprise < 0.01% of maximas
  localMaxima = subset(localMaxima, (localMaxima$height >= minimumHeight) & (is.na(localMaxima$ring1mean) == FALSE) & (is.na(localMaxima$ring2mean) == FALSE))
  
  slopeAspect = rast(file.path(localMaximaPath, "../slopeAspect", paste0(tileName, ".tif")))
  localMaximaSlopeAspect = terra::extract(slopeAspect, localMaxima, method = "simple") # nearest neighbor
  localMaxima$dsmSlope = localMaximaSlopeAspect$dsmSlope
  localMaxima$cmmSlope3 = localMaximaSlopeAspect$cmmSlope3
  
  localMaximaXY = geom(localMaxima)[, c("x", "y")]
  localMaxima = as_tibble(localMaxima) %>% 
    mutate(x = localMaximaXY[, "x"], 
           y = localMaximaXY[, "y"], 
           uniqueID = 1000000 * as.integer(str_c(str_sub(tileName, 2, 6), str_sub(tileName, 8, 12))) + localMaxima$id)
  
  if (is.null(acceptedTileName))
  {
    return(localMaxima)
  }
  
  acceptedTreetopsFilePath = file.path(acceptedTreetopsPath, paste0(acceptedTileName, ".gpkg"))
  if (file.exists(acceptedTreetopsFilePath) == FALSE)
  {
    stop(paste0("Tile ", acceptedTileName, " lacks an accepted treetop GeoPackage at '", acceptedTreetopsFilePath, "'."))
  }
  
  localMaxima$treetop = factor("no", levels = c("no", "yes", "merge", "noise", "maybe noise"))
  
  availableTruthLayers = vector_layers(acceptedTreetopsFilePath)
  if ("treetops" %in% availableTruthLayers)
  {
    acceptedTreetops = vect(acceptedTreetopsFilePath, layer = "treetops")
    isTreetopKnn = get.knnx(localMaximaXY, geom(acceptedTreetops)[, c("x", "y")], k = 1)
    isTreetopKnn = tibble(index = isTreetopKnn$nn.index[, 1], distance = isTreetopKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    localMaxima$treetop[isTreetopKnn$index] = "yes"
  }
  if ("merge treetops" %in% availableTruthLayers)
  {
    mergePoints = vect(acceptedTreetopsFilePath, layer = "merge treetops")
    isMergeKnn = get.knnx(localMaximaXY, geom(mergePoints)[, c("x", "y")], k = 1)
    isMergeKnn = tibble(index = isMergeKnn$nn.index[, 1], distance = isMergeKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    localMaxima$treetop[isMergeKnn$index] = "merge"
  }
  if ("noise points" %in% availableTruthLayers)
  {
    noisePoints = vect(acceptedTreetopsFilePath, layer = "noise points")
    noisePointKnn = get.knnx(localMaximaXY, geom(noisePoints)[, c("x", "y")], k = 1)
    noisePointKnn = tibble(index = noisePointKnn$nn.index[, 1], distance = noisePointKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    localMaxima$treetop[noisePointKnn$index] = "noise"
  }
  if ("maybe noise points" %in% availableTruthLayers)
  {
    maybeNoisePoints = vect(acceptedTreetopsFilePath, layer = "maybe noise points")
    maybeNoisePointKnn = get.knnx(localMaximaXY, geom(maybeNoisePoints)[, c("x", "y")], k = 1)
    maybeNoisePointKnn = tibble(index = maybeNoisePointKnn$nn.index[, 1], distance = maybeNoisePointKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    localMaxima$treetop[maybeNoisePointKnn$index] = "maybe noise"
  }
  
  return(localMaxima)
}

get_treetop_error = function(predicted, tileMaxima)
{
  falsePositives = tileMaxima[which((predicted %in% c("yes", "merge")) & (tileMaxima$treetop %in% c("yes", "merge") == FALSE)), ]
  if (nrow(falsePositives) > 0)
  {
    falseNegatives = tileMaxima[which((predicted %in% c("yes", "merge") == FALSE) & (tileMaxima$treetop %in% c("yes", "merge"))), ]
    falseKnn = get.knnx(falseNegatives[, c("x", "y")], falsePositives[, c("x", "y")], k = 1)
    treetopError = tibble(distanceMae = mean(falseKnn$nn.dist[, 1]),
                          heightMae = mean(abs(falsePositives$height - falseNegatives$height[falseKnn$nn.index[, 1]])))
  } else {
    treetopError = tibble(distanceMae = 0, heightMae = 0)
  }
  
  return(treetopError)
}

get_treetop_eligible_neighborhood = function(tileName, tileMaxima)
{
  xTile = as.integer(str_sub(tileName, 2, 6))
  yTile = as.integer(str_sub(tileName, 8, 12))
  tileExtentBuffered = get_tile_extent(tileName, treetopOptions$neighborhoodBufferWidth)
  
  northwestMaxima = NULL
  northMaxima = NULL
  northeastMaxima = NULL
  eastMaxima = NULL
  westMaxima = NULL
  southwestMaxima = NULL
  southMaxima = NULL
  southeastMaxima = NULL
  
  tileNameNorthwest = sprintf("s%05dw%05d", xTile - 0.01 * treetopOptions$tileSize, yTile + 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameNorthwest, ".gpkg"))))
  {
    northwestMaxima = get_treetop_eligible_maxima(tileNameNorthwest) %>% filter(x >= tileExtentBuffered$xMin, y <= tileExtentBuffered$yMax)
  }
  tileNameNorth = sprintf("s%05dw%05d", xTile, yTile + 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameNorth, ".gpkg"))))
  {
    northMaxima = get_treetop_eligible_maxima(tileNameNorth) %>% filter(y <= tileExtentBuffered$yMax)
  }
  tileNameNortheast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile + 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameNortheast, ".gpkg"))))
  {
    northeastMaxima = get_treetop_eligible_maxima(tileNameNortheast) %>% filter(x <= tileExtentBuffered$xMax, y <= tileExtentBuffered$yMax)
  }

  tileNameWest = sprintf("s%05dw%05d", xTile - 0.01 * treetopOptions$tileSize, yTile)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameWest, ".gpkg"))))
  {
    westMaxima = get_treetop_eligible_maxima(tileNameWest) %>% filter(x >= tileExtentBuffered$xMin)
  }
  tileNameEast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameEast, ".gpkg"))))
  {
    eastMaxima = get_treetop_eligible_maxima(tileNameEast) %>% filter(x <= tileExtentBuffered$xMax)
  }
  
  tileNameSouthwest = sprintf("s%05dw%05d", xTile - 0.01 * treetopOptions$tileSize, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSouthwest, ".gpkg"))))
  {
    southwestMaxima = get_treetop_eligible_maxima(tileNameSouthwest) %>% filter(x >= tileExtentBuffered$xMin, y >= tileExtentBuffered$yMin)
  }
  tileNameSouth = sprintf("s%05dw%05d", xTile, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSouth, ".gpkg"))))
  {
    southMaxima = get_treetop_eligible_maxima(tileNameSouth) %>% filter(y >= tileExtentBuffered$yMin)
  }
  tileNameSoutheast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSoutheast, ".gpkg"))))
  {
    southeastMaxima = get_treetop_eligible_maxima(tileNameSoutheast) %>% filter(x <= tileExtentBuffered$xMax, y >= tileExtentBuffered$yMin)
  }
  
  return(bind_rows(northwestMaxima, northMaxima, northeastMaxima, eastMaxima, tileMaxima, westMaxima, southwestMaxima, southMaxima, southeastMaxima))
}

get_treetop_predictors = function(treetopEligibleMaxima)
{
  localMaximaXY = treetopEligibleMaxima[, c("x", "y")]
  
  neighborKnn = get.knn(localMaximaXY, k = 50)
  neighborKnn = tibble(index1 = neighborKnn$nn.index[, 1], distance1 = neighborKnn$nn.dist[, 1],
                       index2 = neighborKnn$nn.index[, 2], distance2 = neighborKnn$nn.dist[, 2],
                       index3 = neighborKnn$nn.index[, 3], distance3 = neighborKnn$nn.dist[, 3],
                       index4 = neighborKnn$nn.index[, 4], distance4 = neighborKnn$nn.dist[, 4],
                       index5 = neighborKnn$nn.index[, 5], distance5 = neighborKnn$nn.dist[, 5],
                       meanDistance2 = rowMeans(neighborKnn$nn.dist[, 1:2]),
                       meanDistance3 = rowMeans(neighborKnn$nn.dist[, 1:3]),
                       meanDistance4 = rowMeans(neighborKnn$nn.dist[, 1:4]),
                       meanDistance5 = rowMeans(neighborKnn$nn.dist[, 1:5]),
                       meanDistance10 = rowMeans(neighborKnn$nn.dist[, 1:10]),
                       meanDistance20 = rowMeans(neighborKnn$nn.dist[, 1:20]),
                       meanDistance50 = rowMeans(neighborKnn$nn.dist[, 1:50]))
  
  treetopPredictors = as_tibble(treetopEligibleMaxima) %>% 
    mutate(neighbor1sourceID = sourceID[neighborKnn$index1], neighbor1distance = neighborKnn$distance1, neighbor1dsmZ = dsmZ[neighborKnn$index1],
           neighbor2sourceID = sourceID[neighborKnn$index2], neighbor2distance = neighborKnn$distance2, neighbor2dsmZ = dsmZ[neighborKnn$index2],
           neighbor3sourceID = sourceID[neighborKnn$index3], neighbor3distance = neighborKnn$distance3, neighbor3dsmZ = dsmZ[neighborKnn$index3],
           neighbor4sourceID = sourceID[neighborKnn$index4], neighbor4distance = neighborKnn$distance4, neighbor4dsmZ = dsmZ[neighborKnn$index4],
           neighbor5sourceID = sourceID[neighborKnn$index5], neighbor5distance = neighborKnn$distance5, neighbor5dsmZ = dsmZ[neighborKnn$index5],
           neighborDistance2Mean = neighborKnn$meanDistance2, neighborDistance3Mean = neighborKnn$meanDistance3, neighborDistance4Mean = neighborKnn$meanDistance4, 
           neighborDistance5Mean = neighborKnn$meanDistance5, neighborDistance10Mean = neighborKnn$meanDistance10, neighborDistance20Mean = neighborKnn$meanDistance20, neighborDistance50Mean = neighborKnn$meanDistance50,
           # convert all distances from feet to m
           #across(where(is.double), ~0.3048 * .x),
           deltaCmm = dsmZ - cmmZ,
           mean1delta = dsmZ - ring1mean, 
           mean2delta = dsmZ - ring2mean, 
           mean3delta = dsmZ - ring3mean, 
           mean4delta = dsmZ - ring4mean, 
           mean5delta = dsmZ - ring5mean,
           mean1deltaNormalized = mean1delta / height, 
           mean2deltaNormalized = mean2delta / height, 
           mean3deltaNormalized = mean3delta / height, 
           mean4deltaNormalized = mean4delta / height, 
           mean5deltaNormalized = mean5delta / height,
           neighbor1differentSourceID = sourceID != neighbor1sourceID, 
           neighbor2differentSourceID = sourceID != neighbor2sourceID, 
           neighbor3differentSourceID = sourceID != neighbor3sourceID, 
           neighbor4differentSourceID = sourceID != neighbor4sourceID, 
           neighbor5differentSourceID = sourceID != neighbor5sourceID, 
           neighbor1prominence = dsmZ - neighbor1dsmZ, 
           neighbor2prominence = dsmZ - neighbor2dsmZ, 
           neighbor3prominence = dsmZ - neighbor3dsmZ, 
           neighbor4prominence = dsmZ - neighbor4dsmZ, 
           neighbor5prominence = dsmZ - neighbor5dsmZ, 
           neighbor1prominenceNormalized = neighbor1prominence / height,
           neighbor2prominenceNormalized = neighbor2prominence / height,
           neighbor3prominenceNormalized = neighbor3prominence / height,
           neighbor4prominenceNormalized = neighbor4prominence / height,
           neighbor5prominenceNormalized = neighbor5prominence / height,
           neighborDifferentSourceID2 = neighbor1differentSourceID + neighbor2differentSourceID,
           neighborDifferentSourceID3 = neighbor1differentSourceID + neighbor2differentSourceID + neighbor3differentSourceID,
           neighborDifferentSourceID4 = neighbor1differentSourceID + neighbor2differentSourceID + neighbor3differentSourceID + neighbor4differentSourceID,
           neighborDifferentSourceID5 = neighbor1differentSourceID + neighbor2differentSourceID + neighbor3differentSourceID + neighbor4differentSourceID + neighbor5differentSourceID,
           neighborDistance5Variance = 0.2 * ((neighbor1distance - neighborDistance5Mean)^2 + (neighbor2distance - neighborDistance5Mean)^2 + (neighbor3distance - neighborDistance5Mean)^2 + (neighbor4distance - neighborDistance5Mean)^2 + (neighbor5distance - neighborDistance5Mean)^2),
           neighborProminence5Mean = 0.2 * (neighbor1prominence + neighbor2prominence + neighbor3prominence + neighbor4prominence + neighbor5prominence),
           prominence1 = dsmZ - ring1max, 
           prominence2 = dsmZ - ring2max, 
           prominence3 = dsmZ - ring3max, 
           prominence4 = dsmZ - ring4max, 
           prominence5 = dsmZ - ring5max,
           prominence1normalized = prominence1 / height, 
           prominence2normalized = prominence2 / height, 
           prominence3normalized = prominence3 / height, 
           prominence4normalized = prominence4 / height, 
           prominence5normalized = prominence5 / height,
           prominenceMean = 0.2 * (prominence1 + prominence2 + prominence3 + prominence4 + prominence5),
           prominenceMeanNormalized = prominenceMean / height,
           prominenceVariance = 0.2 * ((prominence1 - prominenceMean)^2 + (prominence2 - prominenceMean)^2 + (prominence3 - prominenceMean)^2 + (prominence4 - prominenceMean)^2 + (prominence5 - prominenceMean)^2),
           prominenceVarianceNormalized = prominenceVariance / height,
           prominenceNeighbor5Variance = 0.2 * ((neighbor1prominence - neighborProminence5Mean)^2 + (neighbor2prominence - neighborProminence5Mean)^2 + (neighbor3prominence - neighborProminence5Mean)^2 + (neighbor4prominence - neighborProminence5Mean)^2 + (neighbor5prominence - neighborProminence5Mean)^2),
           prominenceNeighbor5VarianceNormalized = prominenceNeighbor5Variance / height,
           range1 = ring1max - ring1min, 
           range2 = ring2max - ring2min, 
           range3 = ring3max - ring3min, 
           range4 = ring4max - ring4min, 
           range5 = ring5max - ring5min,
           rangeMean = 0.2 * (range1 + range2 + range3 + range4 + range5),
           rangeMeanNormalized = rangeMean / height,
           range1normalized = range1 / height, 
           range2normalized = range2 / height, 
           range3normalized = range3 / height, 
           range4normalized = range4 / height,
           range5normalized = range5 / height,
           rangeVariance = 0.2 * ((range1 - rangeMean)^2 + (range2 - rangeMean)^2 + (range3 - rangeMean)^2 + (range4 - rangeMean)^2 + (range5 - rangeMean)^2),
           rangeVarianceNormalized = rangeVariance / height,
           rangeNeighbor2 = pmax(neighbor1dsmZ, neighbor2dsmZ) - pmin(neighbor1dsmZ, neighbor2dsmZ),
           rangeNeighbor3 = pmax(neighbor1dsmZ, neighbor2dsmZ, neighbor3dsmZ) - pmin(neighbor1dsmZ, neighbor2dsmZ, neighbor3dsmZ),
           rangeNeighbor5 = pmax(neighbor1dsmZ, neighbor2dsmZ, neighbor3dsmZ, neighbor4dsmZ, neighbor5dsmZ) - pmin(neighbor1dsmZ, neighbor2dsmZ, neighbor3dsmZ, neighbor4dsmZ, neighbor5dsmZ),
           rangeNeighbor2Normalized = rangeNeighbor2 / height,
           rangeNeighbor3Normalized = rangeNeighbor3 / height,
           rangeNeighbor5Normalized = rangeNeighbor5 / height,
           netProminence = prominence1 + prominence2 + prominence3 + prominence4 + prominence5,
           netProminenceNeighbor2 = neighbor1prominence + neighbor2prominence,
           netProminenceNeighbor3 = neighbor1prominence + neighbor2prominence + neighbor3prominence,
           netProminenceNeighbor5 = neighbor1prominence + neighbor2prominence + neighbor3prominence + neighbor4prominence + neighbor5prominence,
           netProminenceNormalized = netProminence / height,
           netProminenceNeighbor2Normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance)),
           netProminenceNeighbor3Normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance)),
           netProminenceNeighbor4Normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance + neighbor4prominence / neighbor4distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance + 1 / neighbor4distance)),
           netProminenceNeighbor5Normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance + neighbor4prominence / neighbor4distance + neighbor5prominence / neighbor5distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance + 1 / neighbor4distance + 1 / neighbor5distance)),
           netRange = range1 + range2 + range3 + range4 + range5,
           netRangeNormalized = netRange / height,
           varianceNormalized1 = ring1variance / height, 
           varianceNormalized2 = ring2variance / height, 
           varianceNormalized3 = ring3variance / height, 
           varianceNormalized4 = ring4variance / height, 
           varianceNormalized5 = ring5variance / height)

  return(treetopPredictors)
}

## treetop random forests
# 15 tile load (2.6M maxima, single threaded): ~30 s 9900X + DDR5-5600, ~41 s 5950X + DDR4-3200
loadStart = Sys.time()
s4268maxima = get_treetop_predictors(bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810"), # truthed tiles
                                               get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840"), 
                                               get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810"),
                                               get_treetop_eligible_maxima("s04170w06870"), # adjacent tiles for neighborhood calculations
                                               get_treetop_eligible_maxima("s04200w06870"), # 15 warnings on z coordinate loss on read, one for each tile
                                               get_treetop_eligible_maxima("s04230w06870"), # 2 additional warnings from s04200w06840 and s04230w06810 truth
                                               get_treetop_eligible_maxima("s04170w06840"),
                                               get_treetop_eligible_maxima("s04230w06840"),
                                               get_treetop_eligible_maxima("s04260w06840"),
                                               get_treetop_eligible_maxima("s04170w06810"),
                                               get_treetop_eligible_maxima("s04260w06840"),
                                               get_treetop_eligible_maxima("s04170w06780"),
                                               get_treetop_eligible_maxima("s04200w06780"),
                                               get_treetop_eligible_maxima("s04230w06780"),
                                               get_treetop_eligible_maxima("s04260w06780"))) %>% 
  filter(is.na(cmmSlope3) == FALSE, # is.na(dsmSlope) == FALSE not needed as dsmSlope is not used as a predictor variable
         ((420000 - treetopOptions$neighborhoodBufferWidth) < x) & (x < (426000 + treetopOptions$neighborhoodBufferWidth)), ((680000 - treetopOptions$neighborhoodBufferWidth) < y) & (y < (686000 + treetopOptions$neighborhoodBufferWidth))) %>% # window neighboring tiles, all of s04230w06840 is still included
  select(-neighbor1sourceID, -neighbor2sourceID, -neighbor3sourceID, -neighbor4sourceID, -neighbor5sourceID, # drop most non-portable values (id, sourceID, x, y, dsmZ, and cmmZ are excluded from training by predictor variable selection)
         -ring1max, -ring2max, -ring3max, -ring4max, -ring5max, -ring1mean, -ring2mean, -ring3mean, -ring4mean, -ring5mean, -ring1min, -ring2min, -ring3min, -ring4min, -ring5min, # drop elevations
         -neighbor1dsmZ, -neighbor2dsmZ, -neighbor3dsmZ, -neighbor4dsmZ, -neighbor5dsmZ) 
Sys.time() - loadStart
sum(is.na(s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(-dsmSlope))) # check for incomplete cases in training data
#colSums(is.na(s4268maxima)))
#(sum(s4268maxima$neighbor1distance <= s4268maxima$neighbor2distance) + sum(s4268maxima$neighbor2distance <= s4268maxima$neighbor3distance) + sum(s4268maxima$neighbor3distance <= s4268maxima$neighbor4distance) + sum(s4268maxima$neighbor4distance <= s4268maxima$neighbor5distance)) / nrow(s4268maxima) # check neighbor distance sort; should be exactly 4

# variable selection
if (treetopOptions$includeInvestigatory)
{
  # PCA
  predictorPca = prcomp(~ ., s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(-treetop, -tile, -id, -uniqueID, -sourceID, -x, -y) %>% mutate(across(ends_with("differentSourceID"), as.numeric)), scale = TRUE) # ~20 s @ 621k rows, 9900X
  factoextra::fviz_eig(predictorPca, ncp = 20)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(1, 2), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 0102 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(3, 4), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 0304 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(5, 6), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 0506 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(7, 8), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 0708 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(9, 10), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 0910 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(11, 12), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetop PCA 621k axes 1112 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)

  # MCA  
  predictorFamd = FactoMineR::FAMD(s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(-treetop, -tile, -id, -uniqueID, -sourceID, -x, -y), graph = FALSE, ncp = 6)
  factoextra::fviz_famd_var(predictorFamd, col.ind = "cos2", gradient.cols = c("blue", "orange", "red"), axes = c(5, 6), labelsize = 2, repel = TRUE)
  
  # LDA
  ldaData = s4268maxima %>% filter(is.na(treetop) == FALSE) %>% slice_sample(n = 100000) %>% # ordination plots start to get slow with > 100k points
    select(-tile, -id, -uniqueID, -sourceID, -x, -y, # exclude non-predictors
           -dsmZ, -cmmZ, # exclude non-portable predictors
           -rangeMean, -rangeMeanNormalized, -neighborDistance2Mean, -neighborDistance3Mean, -neighborDistance4Mean, -neighborDistance5Mean, -starts_with("neighborDifferentSourceID"), -neighborDistance5Variance, -neighborProminence5Mean, -prominenceMean, -prominenceMeanNormalized, -rangeNeighbor5Normalized, -netProminence, -netProminenceNormalized, -starts_with("netProminenceNeighbor"), -netRange, -netRangeNormalized) %>% # exclude colinear variables
    mutate(across(where(is.numeric), scale))
  predictorLda = MASS::lda(treetop ~ ., ldaData, tol = 1E-4)
  ggord::ggord(predictorLda, ldaData$treetop, arrow = 0.2, axes = c(1, 2), direction = "both", ext = 0.99, force = 10, grp_title = NULL, max.overlaps = 150, repel = TRUE, size = 0.5, txt = 2.5, vec_ext = 100, veccol = "grey40", xlims = NA * c(-1, 1), ylims = NA * c(-1, 1)) +
  ggord::ggord(predictorLda, ldaData$treetop, arrow = 0.2, axes = c(3, 4), direction = "both", ext = 0.99, force = 10, grp_title = NULL, max.overlaps = 150, repel = TRUE, size = 0.5, txt = 2.5, vec_ext = 100, veccol = "grey40", xlims = NA * c(-1, 1), ylims = NA * c(-1, 1)) +
  patchwork::plot_annotation(theme = theme(plot.margin = margin())) +
  patchwork::plot_layout(guides = "collect")
  
  # rows   predictors   CPU    VSURF  cores selected  accuracy   tune  trees  threads   mtry  min node size  sample fraction
  # 407k   PCA+MCA+LDA  9900X         12    13
  #
  # 164k   48           5950X   3.0h  14    12        ~98%       2.0h  1000   12
  # 308k   48           5950X   6.3h  14    12        ~98%       5.4h  1000   12
  # 308k   78           5950X  15.3h  14    13        99.1%      7.4h   500   15        8     5
  # 461k   99           5950X   1.1d  14    13                   10h    500   14        7     9              0.547
  # 461k   101          5950X   2.2d  14    14                   10h    500   14        8     7              0.549
  library(forcats)
  library(VSURF)
  treetopVsurf = VSURF(treetop ~ ., s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)), ncores = treetopOptions$rangerThreads, parallel = TRUE, RFimplem = "ranger")
  saveRDS(treetopVsurf, "trees/segmentation/treetop vsurf s4268 407k PMCADAV13.Rds")
  #treetopVsurf = readRDS("trees/segmentation/treetop vsurf s4268 407k PMCADAV13.Rds")
  treetopVsurf$nums.varselect # threshold -> interpretation -> prediction
  treetopVsurf$mean.perf # fractional error rate
  treetopVsurf$overall.time
  treetopVsurf$comput.times
  
  plot(treetopVsurf)
  (variablesTreshold = attributes(treetopVsurf$terms[treetopVsurf$varselect.thres])$term.labels)
  (variablesInterpretation = attributes(treetopVsurf$terms[treetopVsurf$varselect.interp])$term.labels)
  (variablesPrediction = attributes(treetopVsurf$terms[treetopVsurf$varselect.pred])$term.labels)
  
  predictorImportance = tibble(predictor = as.character(attr(treetopVsurf$terms, "predvars"))[treetopVsurf$imp.mean.dec.ind + 2], importance = treetopVsurf$imp.mean.dec) %>% # offset as.character() by two since first element is "list" and second is classification
    mutate(importance = importance / sum(importance), selection = factor(if_else(predictor %in% variablesPrediction, "prediction", if_else(predictor %in% variablesInterpretation, "interpretation", if_else(predictor %in% variablesTreshold, "thresholding", "excluded"))), levels = c("prediction", "interpretation", "thresholding", "excluded")))
  ggplot() +
    geom_col(aes(x = importance, y = fct_reorder(predictor, importance), fill = selection), predictorImportance) +
    labs(x = "normalized variable importance", y = NULL, fill = "VSURF") +
    scale_fill_manual(values = c("forestgreen", "blue2", "darkviolet", "black"))
  ggsave("trees/segmentation/treetop vsurf s4268 407k PMCADAV13 importance.png", width = 14, height = 0.33 * nrow(predictorImportance), units = "cm", dpi = 150)
  
  library(tuneRanger)
  library(mlr)
  #predictorVariables = c("treetop", variablesPrediction)
  rangerTuneData = as.data.frame(s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables))) # placing inline to makeClassifTask() fails with "Must have length 1" for certain combinations of predictor variables (but works with most combinations?!)
  rangerTuneTask = makeClassifTask(data = rangerTuneData, target = "treetop")
  #estimateStart = Sys.time()
  #estimateTimeTuneRanger(rangerTuneTask, num.trees = 500, num.threads = 14, iters = 70)# 5950X: 1.5 min to estimate 1h47m @ 164k tops and 12 predictors
  #Sys.time() - estimateStart
  (rangerTuneStart = Sys.time())
  rangerTuning = tuneRanger(rangerTuneTask, measure = list(multiclass.brier), num.trees = 500, num.threads = 12, iters = 70,
                            build.final.model = FALSE)
  (rangerTuneTime = Sys.time() - rangerTuneStart)
  (rangerTuning)
  saveRDS(rangerTuning, "trees/segmentation/treetop ranger tuning s4268 407k PMCADAV13.Rds")
  detach("package:tuneRanger", unload = TRUE)
  detach("package:mlrMBO", unload = TRUE)
  detach("package:mlr", unload = TRUE)
  #rangerTuning = readRDS("trees/segmentation/treetop ranger tuning s4268 407k PMCADAV18.Rds")
}

# fitting
if (treetopOptions$fitRandomForest)
{
  # predictor variable selection
  # s04200w06810 only
  #predictorVariables = c("treetop", "prominenceStdDevNormalized", "height", "radius", "prominence1normalized", "prominence3normalized", "netProminenceNormalized", "prominence4normalized", "prominence5normalized", "netRange", "range1", "deltaCmm", "ring2max")
  # s04200w06810 and s04230w06810
  #predictorVariables = c("treetop", "height", "prominenceStdDevNormalized", "prominence2normalized", "radius", "prominence2", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "prominence4normalized", "netRange", "range2", "range1", "prominence1", "ring1min" )
  #predictorVariables = c("treetop", "height", "radius", "prominence2normalized", "prominence3normalized", "netProminenceNormalized", "prominence1normalized", "ring2variance", "prominence4normalized", "mean2delta", "ring1variance", "netRange", "mean5delta", "neighbor1distance", "deltaCmm")
  # controls
  #predictorVariables = c("treetop", "height", "radius")
  #rangerTuning = tibble(mtry = 2, minNodeSize = 756, sampleFraction = 0.282)
  #predictorVariables = c("treetop", "height", "radius", "netProminenceNormalized")
  #rangerTuning = tibble(mtry = 2, minNodeSize = 143, sampleFraction = 0.255)
  # s04200w06810, s04230w06810, and s04230w06810 VSURF
  #predictorVariables = c("treetop", "prominence2normalized", "height", "netProminenceNormalized", "radius", "ring2variance", "prominence3normalized", "prominence1normalized", "ring1variance", "range2", "prominence4normalized", "mean2deltaNormalized", "netRange", "mean4delta")
  #predictorVariables = c("treetop", "prominence2normalized", "height", "radius", "netProminenceNormalized", "ring2variance", "prominence3normalized", "prominence1normalized", "prominence4normalized", "mean2deltaNormalized", "cmmSlope3", "mean3delta", "rangeMean", "neighbor1prominence", "mean4delta") # legacy14
  #rangerTuning = tibble(mtry = 8, minNodeSize = 7, sampleFraction = 0.549)
  #predictorVariables = c("treetop", "height", "radius", "netProminenceNormalized", "deltaCmm", "cmmSlope3", "prominenceMean", "netProminence", "mean3delta", "rangeMeanNormalized", "netRange", "netRangeNormalized", "rangeVarianceNormalized", "prominenceNeighbor5Variance", "neighborDistance5Mean", "netProminenceNeighbor5", "ring2variance", "prominence1normalized", "prominence2normalized", "prominence3normalized", "prominence4normalized", "mean2deltaNormalized", "rangeMean", "neighbor1prominence", "mean4delta") # PCA+MCA+LDA+VSURF
  #rangerTuning = tibble(mtry = 13, minNodeSize = 4, sampleFraction = 0.629)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "prominenceMean", "netProminence", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "netRange", "mean2deltaNormalized", "netRangeNormalized", "rangeMeanNormalized", "deltaCmm") # PMCDAV18
  #rangerTuning = tibble(mtry = 11, minNodeSize = 7, sampleFraction = 0.577)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2deltaNormalized", "netRangeNormalized", "rangeMeanNormalized", "deltaCmm") # PMCADAV15
  #rangerTuning = tibble(mtry = 8, minNodeSize = 4, sampleFraction = 0.565)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2deltaNormalized", "netRangeNormalized", "deltaCmm") # PMCADAV14
  #rangerTuning = tibble(mtry = 8, minNodeSize = 7, sampleFraction = 0.591)
  #predictorVariables = c("treetop", "height", "radius", "prominence1", "prominence2", "prominence3", "prominence4", "netProminence", "ring2variance", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2delta", "netRange", "deltaCmm") # PMCADAV14un
  #rangerTuning = tibble(mtry = 10, minNodeSize = 5, sampleFraction = 0.533)
  predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "mean2deltaNormalized", "netRangeNormalized", "deltaCmm") # PMCADAV13
  rangerTuning = tibble(mtry = 7, minNodeSize = 11, sampleFraction = 0.549)
  #rangerTuning = tibble(mtry = rangerTuning$recommended.pars$mtry, minNodeSize = rangerTuning$recommended.pars$min.node.size, sampleFraction = rangerTuning$recommended.pars$sample.fraction)
  
  # variables which need to flow through cross validation for accuracy measurements but should not be used for predictors
  accuracyVariables = c("tile", "id", "uniqueID", "sourceID", "x", "y", "dsmZ", "cmmZ") # x and y coordinates are also used for merge point kNNs
  if (sum(accuracyVariables %in% predictorVariables) > 0)
  {
    stop("At least one accuracy assessment variable excluded from training data is indicated as a predictor variable. Should it be removed from the exclusion list?")
  }
  
  # cross validated accuracy estimation
  #                                                                                                                     cross validated accuracy
  # tiles                              maxima predictors        tune          mtry minNode sampling cross validation    treetop  overall
  # s04200+s04230w06810 + s4200w06840  407k   hR                 1.5h 5950X    2   756     0.282    2x25 @ 0.18h 9900X  0.861    0.950
  # s04200+s04230w06810 + s4200w06840  407k   hRnnp              2.0h 5950X    2   143     0.255    2x25 @ 0.20h 9900X  0.940    0.970
  # s04200+s04230w06810 + s4200w06840  407k   legacy14          10h   5950X    8     7     0.549    2x25 @ 0.77h 9900X  0.943    0.976
  # s04200+s04230w06810 + s4200w06840  407k   PMCADAV13          3.66h 9900X   7    11     0.650    2x25 @ 0.80h 9900X  0.943    0.975
  # s04200+s04230w06810 + s4200w06840  407k   PMCADAV14          3.91h 9900X   8     7     0.591    2x25 @ 0.79h 9900X  0.942    0.975
  # s04200+s04230w06810 + s4200w06840  407k   PMCADAV14un        2.99h 9900X  10     5     0.533    2x25 @ 0.58h 9900X  0.940    0.975
  # s04200+s04230w06810 + s4200w06840  407k   PMCADAV15          4.59h 9900X   8     5     0.564    2x25 @ 0.94h 9900X  0.943    0.975
  # s04200+s04230w06810 + s4200w06840  407k   PMCADAV18          5.00h 9900X  11     7     0.577    2x25 @ 1.21h 9900X  0.943    0.975
  # s04200+s04230w06810 + s4200w06840  407k   PCA+MCA+LDA+VSURF  7.89h 9900X  13     4     0.629    2x25 @ 1.61h 9900X  0.944    0.976
  handlers(global = TRUE)
  handlers("cli")
  s4268training = s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables), all_of(accuracyVariables)) # drop maxima in tiles neighboring training region
  (crossValidationStartTime = Sys.time())
  s4268crossValidation = fit_ranger_treetop(s4268training, s4268maxima, 
                                            mtry = rangerTuning$mtry, minNodeSize = rangerTuning$minNodeSize, sampleFraction = rangerTuning$sampleFraction, 
                                            folds = 2, repetitions = 25)
  crossValidationTime = Sys.time() - crossValidationStartTime
  s4268crossValidation %>% summarize(cvTime = crossValidationTime, treetopAccuracy = mean(treetopAccuracy), overallAccuracy = mean(overallAccuracy), treetopMiscountPct = 100 * mean((treetop - expectedTreetop) / expectedTreetop))
  saveRDS(s4268crossValidation, paste0("trees/segmentation/treetopRandomForest s4268 407k PMCADAV13 2x25 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))
  #s4268crossValidation = readRDS(paste0("trees/segmentation/treetopRandomForest s4268 407k PMCADAV13 2x25 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))

  s4268crossValidation$confusionSubmatrix[1]$Accuracy["merge"]
  ggplot() +
    geom_violin(aes(x = overallAccuracy, y = "overall", color = "overall"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    geom_violin(aes(x = treetopAccuracy, y = "treetop", color = "treetop"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    geom_violin(aes(x = nonTreetopAccuracy, y = "not treetop", color = "not treetop"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    coord_cartesian(xlim = c(0.9, 1)) +
    guides(color = "none") +
    labs(x = "accuracy", y = NULL) + 
    scale_y_discrete(limits = rev(c("overall", "treetop", "not treetop")))

  # fit random forest
  treetopRandomForest = ranger::ranger(treetop ~ ., data = s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                       mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                       sample.fraction = rangerTuning$sampleFraction,
                                       num.threads = treetopOptions$rangerThreads)
  treetopRandomForest
  saveRDS(treetopRandomForest, paste0("trees/segmentation/treetopRandomForest s4268 407k PMCADAV13 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))
  #treetopRandomForest = readRDS("trees/segmentation/treetopRandomForest s4268 407k legacy14 m8n25.Rds")

  # variable importance
  randomForestImportance = ranger::ranger(treetop ~ ., data = s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                          importance = "impurity_corrected",
                                          mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                          sample.fraction = rangerTuning$sampleFraction,
                                          num.threads = treetopOptions$rangerThreads)
  randomForestLocalImportance = ranger::ranger(treetop ~ ., data = s4268maxima %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                               importance = "permutation", local.importance = TRUE,
                                               mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                               sample.fraction = rangerTuning$sampleFraction,
                                               num.threads = treetopOptions$rangerThreads)
  localImportance = as_tibble(randomForestLocalImportance$variable.importance.local) %>% mutate(treetop = (s4268maxima %>% filter(is.na(treetop) == FALSE))$treetop) %>% group_by(treetop) %>%
    summarize(across(everything(), mean)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~100 * .x / max(across(where(is.numeric))))) %>%
    rename(any_of(c(NNP = "netProminenceNormalized"))) %>%
    relocate(treetop, any_of(c("height", "radius")))

  ggplot() +
    geom_raster(aes(x = predictor, y = treetop, fill = importance), localImportance %>% pivot_longer(!treetop, names_to = "predictor", values_to = "importance")) +
    scale_x_discrete(limits = names(localImportance %>% select(-treetop)), labels = NULL) +
    scale_y_discrete(limits = rev(c("yes", "merge", "noise", "maybe noise", "no"))) +
  ggplot() +
    geom_raster(aes(x = predictor, y = "all subclasses combined", fill = importance), tibble::as_tibble_row(randomForestImportance$variable.importance) %>% 
                                                                                        rename(any_of(c(NNP = "netProminenceNormalized"))) %>%
                                                                                        pivot_longer(everything(), names_to = "predictor", values_to = "importance") %>%
                                                                                        mutate(importance = 100 * importance / max(importance))) +
    scale_x_discrete(limits = names(localImportance %>% select(-treetop))) +
    plot_annotation(theme = theme(plot.margin = margin())) +
    plot_layout(nrow = 2, ncol = 1, guides = "collect") &
    coord_equal() &
    labs(x = NULL, y = NULL, fill = "variable\nimportance, %") &
    scale_fill_viridis_c(option = "plasma", limits = c(0, 100 + 1E-14)) &
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 10))
  ggsave(file.path(getwd(), paste0("trees/segmentation/treetopRandomForest s4268 407k PMCADAV13 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, " cubic subclass.png")), units = "cm", height = 8, width = 12, dpi = 200)
  
  print(tibble(variable = names(randomForestImportance$variable.importance), importance = randomForestImportance$variable.importance, relativePct = 100 * importance / max(importance)) %>% arrange(desc(importance)), n = 25)
  
  # leave one out cross validation
  s04200w06810data = s4268maxima %>% filter(tile == "s04200w06810")
  s04200w06840data = s4268maxima %>% filter(tile == "s04200w06840")
  s04230w06810data = s4268maxima %>% filter(tile == "s04230w06810")
  # omit s04200w06840: unweighted accuracy = 0.979, κ = 0.933
  s04200s04230w06810forest = readRDS("trees/segmentation/treetopRandomForest vsurf 99.13 s04200+s04230w06810.Rds")
  s04200w06840prediction = predict(s04200s04230w06810forest, s04200w06840data)
  # omit s04230w06810: unweighted accuracy = 0.990, κ = 0.967
  s04200w06810w06840forest = readRDS("trees/segmentation/treetopRandomForest vsurf 99.13 s04200w06810+w06840.Rds")
  s04230w06810prediction = predict(s04200w06810w06840forest, s04230w06810data)
  # omit s04200w06810: unweighted accuracy = 0.979, κ = 0.904, weighted results within 0.002%
  s04200w06840s04230w06810forest = readRDS("trees/segmentation/treetopRandomForest vsurf 99.13 s04200w06840+s04230w06810.Rds")
  s04200w06810prediction = predict(s04200w06840s04230w06810forest, s04200w06810data)
  
  (tileLooAccuracy = bind_rows(get_treetop_accuracy(s04200w06840prediction, s04200w06840data) %>% mutate(tile = "s04200w06840"),
                               get_treetop_accuracy(s04200w06810prediction, s04200w06810data) %>% mutate(tile = "s04200w06810"),
                               get_treetop_accuracy(s04230w06810prediction, s04230w06810data) %>% mutate(tile = "s04230w06810")) %>% relocate(tile))
  (tileLooError = bind_rows(get_treetop_error(s04200w06840prediction, s04200w06840data) %>% mutate(tile = "s04200w06840"),
                            get_treetop_error(s04200w06810prediction, s04200w06810data) %>% mutate(tile = "s04200w06810"),
                            get_treetop_error(s04230w06810prediction, s04230w06810data) %>% mutate(tile = "s04230w06810")) %>% relocate(tile))
  
  
  
  # sanity check for random forest's expected near perfect recall
  #confusionMatrix(randomForestFit)
  randomForestConfusionMatrix = confusionMatrix(predict(treetopRandomForest, s4268maxima %>% select(-treetop)), s4268maxima$treetop)
  randomForestConfusionMatrix$table
  randomForestConfusionMatrix$overall
}


## predict treetops, merge, and noise points
if (treetopOptions$includeInvestigatory)
{
  #randomForestFit = readRDS("trees/segmentation/treetopRandomForest vsurf 48.14 s04200+s04230w06810.Rds")
  #randomForestFit = readRDS("trees/segmentation/treetopRandomForest vsurf 78.13 s04200+s04230w06810.Rds")
  treetopRandomForest = readRDS("trees/segmentation/treetopRandomForest vsurf 103.14 s04200w06840 + s4200+s04230w06810.Rds")
  tileName = "s04230w06840"
  tileMaxima = get_treetop_eligible_maxima(tileName) %>% filter(is.na(cmmSlope3) == FALSE)
  tileMaxima$treetop = predict(treetopRandomForest$finalModel, get_treetop_predictors(tileMaxima))$predictions
  table(tileMaxima$treetop)
  #writeVector(vect(tileMaxima, crs = tileCrs, geom = c("x", "y")), file.path(candidateTreetopsPath, "rf", paste0(tileName, " ranger.gpkg")), layer = "treetops", insert = TRUE, overwrite = TRUE)
  
  # group merge points and define merge treetops
  # In general, relative to the initial random forest classification, number of non-treetop maxima declines, number of singleton 
  # treetops increases, number of merge points increases, and noise points remain unchanged.
  tileNeighborhood = get_treetop_eligible_neighborhood(tileName, tileMaxima) # ~20 s, 9900X
  neighborhoodMaxima = tileNeighborhood
  tileMergePoints = get_merge_points(tileMaxima, tileNeighborhood)
  tileMaxima$treetop[tileMergePoints$mergePointIndices] = "merge" # update classifications based on merge point clustering results
  tileMaxima$treetop[tileMergePoints$treetopIndices] = "yes"
  table(tileMaxima$treetop)
  
  #                                tile           no           yes          merge         noise        maybe noise
  # from random forest             s04230w06840   143899       14453        1114          42           0             before merging
  # v0: 25 DSM @ 4 + h^0.45        s04230w06840   143557       14411        1498          42           0
  # v1: 1/60 + 1/40                s04230w06840   143441       14356        1669          42           0
  
  # write tile's GeoPackage
  # layers are treetops (single maxima and from merges)
  tileTreetops = bind_rows(tileMaxima %>% filter(tileMaxima$treetop == "yes"), tileMergePoints$treetops) %>%
    mutate(maxima = replace_na(maxima, as.integer(1)), sourceIDs = replace_na(maxima, as.integer(1))) # treetops not from merges have one maxima and one source ID but these values are only set on data from get_merge_points()
  
  treetopsFilePath = file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg"))
  writeVector(vect(tileTreetops, crs = tileCrs, geom = c("x", "y")), treetopsFilePath, layer = "treetops", insert = TRUE, overwrite = TRUE)
  writeVector(vect(tileMaxima %>% filter(tileMaxima$treetop == "merge"), crs = tileCrs, geom = c("x", "y")), treetopsFilePath, layer = "merge points", insert = TRUE, overwrite = TRUE)
  
  noisePoints = vect(tileMaxima %>% filter(treetop == "noise"), crs = tileCrs, geom = c("x", "y"))
  writeVector(noisePoints, treetopsFilePath, layer = "noise points", insert = TRUE, overwrite = TRUE)
  maybeNoisePoints = vect(tileMaxima %>% filter(treetop == "maybe noise"), crs = tileCrs, geom = c("x", "y"))
  writeVector(maybeNoisePoints, treetopsFilePath, layer = "maybe noise points", insert = TRUE, overwrite = TRUE)
  
  # checks
  range(tileNeighborhood$x)
  range(tileMaxima$x)
  range(tileNeighborhood$y)
  range(tileMaxima$y)
}


## merge radius functions
if (treetopOptions$includeInvestigatory)
{
  mergeRadiiEnglish = crossing(height = seq(0, 200), scale = seq(0.1, 1, by = 0.1)) %>%
    mutate(radiusMultiply = scale*height,
           radiusPower = height^scale)
  ggplot() +
    geom_line(aes(x = radiusPower, y = height, group = scale, color = as.factor(scale), linetype = "power"), mergeRadiiEnglish) +
    geom_line(aes(x = radiusMultiply, y = height, group = scale, color = as.factor(scale), linetype = "multiplier"), mergeRadiiEnglish) +
    labs(x = "merge radius in vertically exaggerated space, ft", y = "height, ft", color = NULL, linetype = NULL) +
    scale_linetype_manual(breaks = c("power", "multiplier"), values = c("solid", "dashed"))
}


## diff two treetop identifications
if (treetopOptions$includeInvestigatory)
{
  diff_treetops = function(treetops1, treetops2)
  {
    return(tibble(treeID = treetops1$treeID != treetops2$treeID,
                  x = treetops1$x != treetops2$x,
                  y = treetops1$y != treetops2$y,
                  z = treetops1$z != treetops2$z,
                  height = treetops1$height != treetops2$height))
  }
  
  read_treetops = function(geoPackagePath)
  {
    geoPackage = read_sf(geoPackagePath) # terra::vect() doesn't support z
    coordinates = st_coordinates(geoPackage$geom) 
    return(tibble(treeID = geoPackage$treeID, x = coordinates[, "X"], y = coordinates[, "Y"], z = coordinates[, "Z"], height = geoPackage$height))
  }
  
  referenceTreetops = read_treetops("GIS/DOGAMI/2021 OLC Coos County/treetops DSM/s04020w06690.gpkg")
  newTreetops = read_treetops("GIS/DOGAMI/2021 OLC Coos County/treetops DSM/s04020w06690 new.gpkg")
  
  diff = diff_treetops(referenceTreetops, newTreetops)
  diff %>% filter((treeID == TRUE) | (x == TRUE) | (y == TRUE) | (z == TRUE) | (height == TRUE))
  
  referenceTreetops %>% arrange(desc(y), x, z)
  newTreetops %>% arrange(desc(y), x, z)
}


## distribution of merge points
if (treetopOptions$includeInvestigatory)
{
  get_knnx = function(data, query, k)
  {
    knn = get.knnx(data, query, k)
    knnTibble = tibble(index = knn$nn.index[, 1], distance1 = knn$nn.dist[, 1])
    if (k > 1)
    {
      knnTibble$distance2 = knn$nn.dist[, 2]
    }
    return(knnTibble)
  }
  
  localMaxima = vect("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta/local maxima/s04230w06810.gpkg", layer = "localMaxima")
  mergeTreetops = vect("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted/s04230w06810.gpkg", layer = "merge treetops")

  maximaKnn = get_knnx(geom(localMaxima)[, c("x", "y")], geom(mergeTreetops)[, c("x", "y")], k = 1)
  isMergeKnn = get_knnx(geom(mergeTreetops)[, c("x", "y")], geom(mergeTreetops)[, c("x", "y")], k = 2)
  
  mergeTreetopsGeometry = geom(mergeTreetops)[, c("x", "y")]
  mergeDistance = tibble(x = mergeTreetopsGeometry[, "x"], y = mergeTreetopsGeometry[, "y"], z = localMaxima$dsmZ[maximaKnn$index], sourceID = localMaxima$sourceID[maximaKnn$index], height = localMaxima$height[maximaKnn$index], mergeDistance = isMergeKnn$distance2, maximaDistance = maximaKnn$distance1) %>%
    mutate(x = x, y = y, z = z, height = height, mergeDistance = mergeDistance, maximaDistance = maximaDistance)
  
  ggplot() +
    geom_point(aes(x = mergeDistance, y = height), mergeDistance, alpha = 0.05, shape = 16) +
    labs(x = "distance to nearest merge point, m", y = "height, m")
  
  mergeDistance %>% filter(mergeDistance > 5) # merge distance will be long when other merge points are in a different tile
}


## comparison of tree height distributions between 2016 cruise data and 2021 LiDAR
# trees2016 from trees/height-diameter/setup.R
if (treetopOptions$includeInvestigatory)
{
  #trees2016 %>% group_by(StandID) %>% summarize(standArea = sum(standArea[1])) %>% summarize(measureArea = sprintf("%.4f", sum(standArea)))
  trees2016height = trees2016 %>% mutate(heightClass = if_else(TotalHt > 0, round(TotalHt, 0), if_else(Ht2 > 0, round(Ht2, 0), NA_real_))) %>%
    filter(is.na(heightClass) == FALSE) %>%
    group_by(StandID) %>%
    mutate(liveHeightMeasureTph = sum(isLive * measureTreeTphContribution)) %>%
    group_by(speciesGroup, heightClass) %>%
    summarize(tph = sum(standArea * tph / liveHeightMeasureTph * measureTreeTphContribution) / 15981.0273, .groups = "drop") # 15981.0273 ha = total area of stands measured
  
  dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops"
  treetopTiles = read_sf(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg"), "Elliott tile index")
  treetopTileNames = (treetopTiles %>% filter(bufferDistance == 0))$Tile_ID
  trees2021 = list()
  for (treetopFileName in treetopTileNames) # ~1.8 minutes to load all 663 tiles
  {
    treetopTileName = str_remove(treetopFileName, "\\.gpkg")
    tile = read_sf(file.path(dataSourcePath, str_c(treetopFileName, ".gpkg")))
    tileCoordinates = st_coordinates(tile$geom)
    tile = tibble(tileName = treetopTileName, treeID = tile$treeID, x = tileCoordinates[, "X"], y = tileCoordinates[, "Y"], elevation = tileCoordinates[, "Z"], height = tile$height)
    trees2021[[treetopTileName]] = tile
  }
  trees2021 = bind_rows(trees2021)
  
  segmentedArea = length(treetopTileNames) * 0.3048^2 * 3000^2 / 10000 # ha
  trees2021height = left_join(left_join(trees2021 %>% mutate(heightClass = round(0.3048 * height, 0)) %>% # ~4s
                                          group_by(heightClass) %>%
                                          summarize(segmentedTph2021 = n() / segmentedArea),
                                        trees2016height %>% group_by(heightClass) %>% summarize(tph = sum(tph)) %>% rename(tph2016 = tph),
                                        by = c("heightClass")),
                              # crude height growth estimate: concave down parabola
                              trees2016height %>% group_by(heightClass) %>% summarize(tph = sum(tph), .groups = "drop") %>% 
                                mutate(heightClass = round(heightClass + (2021 - 2015) * pmax(0.5 + 0.06 * heightClass - 0.0012 * heightClass^2, 0), 0)) %>% 
                                group_by(heightClass) %>% summarize(tph = sum(tph), .groups = "drop") %>% rename(grownTph2021 = tph),
                              by = c("heightClass"))
  #print(trees2021 %>% filter(height > 300), n = 75)
  #ggplot() + geom_line(aes(x = seq(0, 80), y = pmax(0.5 + 0.06 * seq(0, 80) - 0.0012 * seq(0, 80)^2, 0)))
  
  ggplot() +
    geom_col(aes(x = heightClass, y = tph, group = speciesGroup, fill = speciesGroup), trees2016height) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 55)) +
    labs(x = "2016 height, m", y = "trees per hectare", fill = NULL) +
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggplot() +
    geom_col(aes(x = heightClass, y = segmentedTph2021), trees2021height) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 55)) +
    labs(x = "2021 height, m", y = NULL) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggplot() +
    geom_line(aes(x = heightClass, y = 100 * segmentedTph2021 / grownTph2021), trees2021height, na.rm = TRUE) +
    geom_smooth(aes(x = heightClass, y = 100 * segmentedTph2021 / grownTph2021), trees2021height, alpha = 0.1, formula = y ~ x, method = "loess", na.rm = TRUE, span = 0.5) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = "2021 height, m", y = "estimated fraction of treetops detected, %") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  tibble(tph2016 = sum(trees2016height$tph), tph2021 = sum(trees2021height$tph)) %>% 
    mutate(total2016M = (33408.277 + 318.838) * tph2016 / 1E6,
           total2021M = (33408.277 + 318.838) * tph2021 / 1E6,
           segmentedPct = 100 * tph2021 / tph2016)
}

# break 2009 DSM from dsmDtmJob.R into tiles after conversion to EPSG:6557 alignment in QGIS to match 2021 flight
# https://gis.stackexchange.com/questions/441960/qgis-clipping-virtual-raster-into-tiles-using-grid-layer-polygons
if (treetopOptions$includeSetup)
{
  tiles = st_read(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg"), quiet = TRUE)
  dsm2009 = rast("D:/Elliott/GIS/DOGAMI/2009 OLC South Coast/DSM 2009 tiles/DSM.tif")
  
  tile = tiles %>% filter(Tile_ID == "s03840w07290")
  
  lapply(1:nrow(tiles), function(tileIndex) {
    tile = tiles[tileIndex,]
    dsmTile = crop(dsm2009, tile, ext = TRUE)
    writeRaster(dsmTile, file.path("D:/Elliott/GIS/DOGAMI/2009 OLC South Coast/DSM", paste0(tile$Tile_ID, ".tif")), datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL9"), overwrite = TRUE)
  })
  
  dsm2009tiles = file.path(list.files("D:/Elliott/GIS/DOGAMI/2009 OLC South Coast/DSM", pattern = "\\.tif$"))
  vrt(file.path("D:/Elliott/GIS/DOGAMI/2009 OLC South Coast/DSM", dsm2009tiles), "D:/Elliott/GIS/DOGAMI/2009 OLC South Coast/DSM/DSM.vrt", overwrite = TRUE, set_names = TRUE)
}

# ring DSM diagnostics: statistics
if (treetopOptions$includeInvestigatory)
{
  library(dplyr)
  library(ggplot2)
  library(terra)
  
  dsmTileDiagnostics = as_tibble(rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM with outlier rejection/ring diagnostics/s04230w06810.tif"))
  treetops = as_tibble(vect("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM with outlier rejection/ring diagnostics/s04230w06810.gpkg")) %>%
    rename(netProminence = `1`, rangeProminence = `2`, totalProminence = `3`, totalRange = `4`, radius2 = `5`) %>%
    mutate(netPromTotalRangeRatio = totalRange / netProminence)
  suspectTops = as_tibble(vect("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM with outlier rejection/ring diagnostics/suspect tops.gpkg")) %>%
    rename(netProminence = `1`, rangeProminence = `2`, totalProminence = `3`, totalRange = `4`, radius2 = `5`) %>%
    mutate(netPromTotalRangeRatio = totalRange / netProminence) %>%
    select(-netProminenceNormalized)
  
  ggplot() +
    geom_segment(aes(x = 0.02, xend = 0.02, y = 0, yend = 1.4), color = "grey70", linewidth = 0.3, linetype = "longdash") +
    geom_bin2d(aes(x = `net prominence normalized`, y = `total prominence normalized`), dsmTileDiagnostics, binwidth = c(0.03333, 0.01)) +
    labs(x = "net prominence, normalized", y = "total prominence, normalized", fill = "treetop\ncandidates", title = "(a) candidates") +
  ggplot() +
    geom_segment(aes(x = 0.02, xend = 0.02, y = 0, yend = 1.4), color = "grey70", linewidth = 0.3, linetype = "longdash") +
    geom_bin2d(aes(x = netProminence, y = totalProminence), treetops, binwidth = c(0.03333, 0.01)) +
    labs(x = "net prominence, normalized", y = NULL, fill = "treetop\ncandidates", title = "(b) accepted") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(widths = c(4.5, 3), guides = "collect") &
    scale_fill_viridis_c(limits = c(0, 220))
  
  suspectTops %>% filter(type == "branch")
  
  ggplot() +
    geom_histogram(aes(x = totalProminence), treetops %>% filter(totalProminence > -100))
}

# distribution of local maxima radii
if (treetopOptions$includeInvestigatory)
{
  localMaxima = as_tibble(rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM with outlier rejection/s04200w06810 local maxima.tif")) %>% rename(localMaximaRadius = `s04200w06810 local maxima`)
  localMaximaEdf = localMaxima %>% group_by(localMaximaRadius) %>% summarize(n = n()) %>% mutate(edf = n / sum(n), cdf = cumsum(edf))
  (pTreetop = (17785 + 434) / nrow(localMaxima))
  
  ggplot() +
    geom_segment(aes(x = 0, y = 1 - pTreetop, xend = 10, yend = 1 - pTreetop), color = "grey70", linetype = "longdash", linewidth = 0.3) +
    geom_step(aes(x = localMaximaRadius, y = cdf), localMaximaEdf, direction = "mid") +
    annotate("text", x = 10, y = 1 - pTreetop, label = "detected treetop fraction", color = "grey70", hjust = 1, vjust = -0.2, size = 3.0) +
    labs(x = "local maxima radius, 46 cm DSM cells", y = "cumulative probability") +
    scale_x_continuous(breaks = seq(0, 10), minor_breaks = NULL) +
    scale_y_continuous(labels = scales::label_percent())
}

# radius selection
if (treetopOptions$includeInvestigatory)
{
  localMaximaR = tibble(height = seq(0, 90), radius = floor(pmin(0.045 * height + 0.5, 4.0) / 0.46 + 0.5),
                        innerRadius = pmax(pmin(floor(radius / 2 + 0.5), 3), 1))
  
  ggplot() +
    geom_step(aes(x = radius, y = height, color = "maximum"), localMaximaR, direction = "hv") +
    geom_step(aes(x = innerRadius, y = height, color = "inner"), localMaximaR, direction = "hv") +
    coord_cartesian(xlim = c(0, 10)) +
    labs(x = "ring radius, 46 cm DSM cells", y = "treetop height, m", color = NULL) +
    scale_x_continuous(breaks = seq(0, 10), minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0, 90, by = 10))
}

# default Gaussian filter creation via Pascal's triangle
# https://dsp.stackexchange.com/questions/10057/gaussian-blur-standard-deviation-radius-and-kernel-size
if (treetopOptions$includeSetup)
{
  pascal3 = c(1, 2, 1) / sum(c(1, 2, 1))
  pascal3 * matrix(pascal3, nrow = length(pascal3), ncol = length(pascal3), byrow = TRUE)
  
  pascal5 = c(1, 4, 6, 4, 1) / sum(c(1, 4, 6, 4, 1))
  pascal5 * matrix(pascal5, nrow = length(pascal5), ncol = length(pascal5), byrow = TRUE)
  
  pascal7 = c(1, 6, 15, 20, 15, 6, 1) / sum(c(1, 6, 15, 20, 15, 6, 1))
  pascal7 * matrix(pascal7, nrow = length(pascal7), ncol = length(pascal7), byrow = TRUE)
}


## convert LiDAR scan times from adjusted standard GPS seconds to local time
if (treetopOptions$includeInvestigatory)
{
  library(lubridate)
  lasOrigin = make_datetime(2011, 9, 14, 1, 46, 25, "UTC")
  
  # Coos Bay sunrise + sunset, August 30 2021: 6:39 AM - 7:55 PM
  lasOrigin + dseconds(314405000) + dhours(-7) # 2021-08-30 5:29PM UTC-7
  lasOrigin + dseconds(314408100) + dhours(-7) # 2021-08-30 6:21PM UTC-7
  lasOrigin + dseconds(314560690) + dhours(-7) # 2021-09-01 5:04AM UTC-7
  
  scanTime = lasOrigin + dseconds(314406103.75467044) + dhours(-7)
  tz(scanTime) = "PST8PDT" # simply sets timezone, does not adjust time from UTC to local
  format(scanTime, "%Y-%m-%dT%H:%M:%S %Z")
  
  scanTime = lasOrigin + dseconds(314416843.39940876) + dhours(-7)
}


## merge treetops into single GeoPackage after Get-Treetops has processed all tiles
# Now handled by MergeTreetops cmdlet in Clouds.
if (treetopOptions$includeInvestigatory)
{
  treetopSourcePath = file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/treetops DSM ring")
  treetopFilePaths = file.path(treetopSourcePath, list.files(treetopSourcePath, "\\.gpkg$"))
  treetopLayers = list()
  for (treetopFileIndex in 1:length(treetopFilePaths)) # 2 min, 47 s with terra but terra drops Z coordinates
  {
   treetopLayers[[treetopFileIndex]] = read_sf(treetopFilePaths[treetopFileIndex])
  }
  
  treetops = data.table::rbindlist(treetopLayers) # https://github.com/r-spatial/sf/issues/2254
}


## sun positions from image centers
if (treetopOptions$includeSetup)
{
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(magrittr)
  library(readr)
  library(suntools)
  library(sf)

  imageColumnTypes = cols(lift = "i", station = "c", sequenceNumber = "i", localDate = "c", localTime = "c", year = "i", month = "i", day = "i", hour = "i", second = "i", .default = "d")
  imagePositions = read_csv(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/210319_Elliot_SF_Lifts_Sep-2021_FMK2.csv"), col_types = imageColumnTypes) %>%
    mutate(dateTime = make_datetime(year, month, day, hour, minute, second, tz = "America/Los_Angeles"))
  sunPositions = solarpos(st_as_sf(imagePositions[, c(21, 20, 16)], coords = c(1, 2, 3), crs = st_crs(4326)), imagePositions$dateTime)
  imagePositions %<>% mutate(sunElevation = sunPositions[, 2], sunAzimuth = sunPositions[, 1])

  ggplot(imagePositions) +
    geom_line(aes(x = hourOfDay, y = sunElevation, color = "elevation")) +
    geom_line(aes(x = hourOfDay, y = sunAzimuth, color = "azimuth")) +
    labs(x = "hour of day", y = "sun angle, degrees", color = "sun position")
  
  write_csv(imagePositions %>% select(-dateTime), "GIS/DOGAMI/2021 OLC Coos County/210319_Elliot_SF_Lifts_Sep-2021_FMK2 with sun position.csv")
}

 
## create review grid for tile
# After .gpkg creation import into QT Modeler and, in Edit -> Edit Style -> Display Style,
#  3D Display Style -> Show as Terrain Hugging
#  Polygon Style -> uncheck filled
if (treetopOptions$includeSetup)
{
  tileName = "s04200w06840"
  elliott2021tiles = vect("GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg", layer = "Elliott tile index")
  elliott2021tiles = subset(elliott2021tiles, elliott2021tiles$Tile_ID == tileName)
  tileReviewGrid = as.polygons(rast(elliott2021tiles, resolution = 3.28084 * 40))
  writeVector(tileReviewGrid, file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted", paste0(tileName, " review grid.gpkg")))
}


## tile checks
if (treetopOptions$includeInvestigatory)
{
  dsmThreading = tibble(tiles = seq(1, 20), readThreads = 1, workerThreads = pmin(pmin(tiles, pmax(25 - floor(1.5*tiles), 2))), 16 - 1)
  ggplot() +
    geom_line(aes(x = tiles, y = workerThreads), dsmThreading) +
    labs(x = "tiles", y = "worker threads")
  
  library(lidR)
  #tile = readLAS("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/tiles RGB+NIR/s03780w06390.las")
  tile = readLAS("C:/Users/westjoh/PhD/data/McDonald-Dunn/Stand 50603/scan 3/scan 3 RGB+class registered -5.las")
  tile2 = readLAS("C:/Users/westjoh/PhD/data/McDonald-Dunn/Stand 50603/scan 3/scan 3 RGB+class registered +5.las")
  highNoiseIndex = which(tile$Z > 2684736)
  head(sort(tile$Z, decreasing = TRUE))
  unique(tile$PointSourceID)
  
  tileFixup = readLAS("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/tiles testing/s03780w06390.las")
  tileDiff = tibble(x = sum(tile$X != tileFixup$X), y = sum(tile$Y != tileFixup$Y), z = sum(tile$Z != tileFixup$Z), intensity = sum(tile$Intensity != tileFixup$Intensity),
                    classification = sum(tile$Classification != tileFixup$Classification), gpstime = sum(tile$gpstime != tileFixup$gpstime), returnNumber = sum(tile$ReturnNumber != tileFixup$ReturnNumber),
                    numReturns = sum(tile$NumberOfReturns != tileFixup$NumberOfReturns), scan = sum(tile$ScanDirectionFlag != tileFixup$ScanDirectionFlag), edge = sum(tile$EdgeOfFlightline != tileFixup$EdgeOfFlightline),
                    key = sum(tile$Keypoint_flag != tileFixup$Keypoint_flag), withheld = sum(tile$Withheld_flag != tileFixup$Withheld_flag), overlap = sum(tile$Overlap_flag != tileFixup$Overlap_flag),
                    scanAngle = sum(tile$ScanAngle != tileFixup$ScanAngle), userData = sum(tile$UserData != tileFixup$UserData), pointSource = sum(tile$PointSourceID != tileFixup$PointSourceID),
                    r = sum(tile$R != tileFixup$R), g = sum(tile$G != tileFixup$G), b = sum(tile$B != tileFixup$B), nir = sum(tile$NIR != tileFixup$NIR))
  print(tileDiff, width = Inf)
  
  scan = readLAS(file.path(getwd(), "../data/McDonald-Dunn/Stand 50603/scan 1/scan 1 RGB+class.las"))
  registeredScan = readLAS(file.path(getwd(), "../data/McDonald-Dunn/Stand 50603/scan 1/scan 1 RGB+class registered.las"))
  registeredScan2 = readLAS(file.path(getwd(), "../data/McDonald-Dunn/Stand 50603/scan 2/scan 2 RGB+class registered.las"))

  scanDiff = tibble(x = sum(scan$X != registeredScan$X), y = sum(scan$Y != registeredScan$Y), z = sum(scan$Z != registeredScan$Z), intensity = sum(scan$Intensity != registeredScan$Intensity),
                    classification = sum(scan$Classification != registeredScan$Classification), gpstime = sum(scan$gpstime != registeredScan$gpstime), returnNumber = sum(scan$ReturnNumber != registeredScan$ReturnNumber),
                    numReturns = sum(scan$NumberOfReturns != registeredScan$NumberOfReturns), scanDir = sum(scan$ScanDirectionFlag != registeredScan$ScanDirectionFlag), edge = sum(scan$EdgeOfFlightline != registeredScan$EdgeOfFlightline),
                    key = sum(scan$Keypoint_flag != registeredScan$Keypoint_flag), withheld = sum(scan$Withheld_flag != registeredScan$Withheld_flag), overlap = sum(scan$Overlap_flag != registeredScan$Overlap_flag),
                    scanAngle = sum(scan$ScanAngle != registeredScan$ScanAngle), userData = sum(scan$UserData != registeredScan$UserData), pointSource = sum(scan$PointSourceID != registeredScan$PointSourceID),
                    r = sum(scan$R != registeredScan$R), g = sum(scan$G != registeredScan$G), b = sum(scan$B != registeredScan$B), nir = sum(scan$NIR != registeredScan$NIR))
  print(scanDiff, width = Inf)
}


## treetop SVM
if (treetopOptions$includeInvestigatory)
{
  #library(kernlab)
  #ggplot() +
  #  geom_point(aes(x = ringRange1, y = ringProminence1, color = isTreetop), ringDataSample, alpha = 0.05, shape = 16) +
  #  labs(x = "ring 1 range, ft", y = "ring 1 prominence, ft", color = "treetop") +
  #ggplot() +
  #  geom_point(aes(x = ringRange2, y = ringProminence2, color = isTreetop), ringDataSample, alpha = 0.05, shape = 16, na.rm = TRUE) +
  #  labs(x = "ring 2 range, ft", y = "ring 2 prominence, ft", color = "treetop") +
  #ggplot() +
  #  geom_point(aes(x = ringRange3, y = ringProminence3, color = isTreetop), ringDataSample, alpha = 0.05, shape = 16, na.rm = TRUE) +
  #  labs(x = "ring 3 range, ft", y = "ring 3 prominence, ft", color = "treetop") +
  #ggplot() +
  #  geom_point(aes(x = ringRange4, y = ringProminence4, color = isTreetop), ringDataSample, alpha = 0.05, shape = 16, na.rm = TRUE) +
  #  labs(x = "ring 4 range, ft", y = "ring 4 prominence, ft", color = "treetop") +
  #plot_annotation(theme = theme(plot.margin = margin())) +
  #plot_layout(guides = "collect") &
  #  guides(color = guide_legend(override.aes = list(alpha = 0.5)))
  #
  # predictors                        accuracy
  # netProminenceNormalized           0.979
  # height, netProminence             0.958
  # height, netProminenceNormalized   0.979
  # fitStart = Sys.time() # ~8 s, 
  # svmFitLinear = ksvm(isTreetop ~ ., data = ringDataSample %>% select(isTreetop, height, netProminenceNormalized))
  # (svmFitTimeLinear = Sys.time() - fitStart)
  # 1 - svmFitLinear@error
  # 
  # svmWeights = colSums(coef(svmFitLinear)[[1]] * (ringDataSample %>% select(height, netProminence))[unlist(alphaindex(svmFitLinear)), ])
  # svmB = b(svmFitLinear)
  # 
  # ggplot() +
  #   geom_abline(slope = -svmWeights[2] / svmWeights[1], intercept = svmB / svmWeights[1]) +
  #   geom_point(aes(x = netProminenceNormalized, y = height, color = isTreetop), ringDataSample, alpha = 0.05, shape = 16) +
  #   labs(x = "net prominence, normalized", y = "height, ft", color = "treetop")
  # 
  # 
  #library(e1071)
  #fitStart = Sys.time() # 
  #svmFitLinear = svm(isTreetop ~ ., data = ringDataSample %>% select(isTreetop, height, netProminence), kernel = "linear", degree = 3)
  #(svmFitTimeLinear = Sys.time() - fitStart)
  #1 - sum(ringDataSample$isTreetop != svmFitLinear$fitted) / nrow(ringDataSample)
  #
  #svmWeights = t(svmFitLinear$coefs) %*% svmFitLinear$SV
  #plot(svmFitLinear, ringDataSample, height ~ netProminenceNormalized, dataSymbol = 16)
  #
  #library(caret)
  #fitStart = Sys.time()
  #svmFitLinear = train(isTreetop ~ ., data = ringData %>% sample_n(25000), method = "svmLinear", # ~3 s @ 25k rows, ~9 s @ 50k, k x r = 2 x 1, no parallel, verboseIter has no effect
  #                     trControl = trainControl(method = "cv", number = 2)) # preProcess = c("center", "scale"), 
  #                     #tuneGrid = data.frame(C = c(3.8, 4.0, 4.2, 4.4, 4.6, 4.8))) # c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1) for polygon statistics
  #svmFitTimeLinear = Sys.time() - fitStart
  #getModelInfo(svmFitLinear)
  #svmWeights = svmFitLinear$finalModel@coef[[1]] %*% svmFitLinear$finalModel@xmatrix[[1]]
  #kernlab::plot(svmFitLinear$finalModel, ringData %>% sample_n(25000))
}
  

## modified normalized net prominence requirement
if (treetopOptions$includeInvestigatory)
{
  #netProminenceMod = tibble(totalProminence = seq(-5, 10, length.out = 100), mod = 0.02 / (1 + exp(-totalProminence)))
  #
  #ggplot() +
  #  geom_line(aes(x = totalProminence, y = mod), netProminenceMod) + 
  #  #coord_cartesian(ylim = c(0, 0.02)) +
  #  labs(x = "total prominence", y = "normalized net prominence threshold")
  
  
  ## tile diagnostics
  #treetopDiagnostics = vect(file.path(dataPath, "ring diagnostics", "s04200w06810 control.gpkg"))
  #
  #treetopDiagnosticsKnn = get.knnx(geom(treetopDiagnostics)[, c("x", "y")], geom(acceptedTreetops)[, c("x", "y")], k = 1)
  #treetopDiagnosticsKnn = tibble(distance = treetopDiagnosticsKnn$nn.dist[, 1], diagnosticIndex = treetopDiagnosticsKnn$nn.index[, 1]) %>%
  #  mutate(acceptedIndex = row_number())
  #
  #treetopDiagnosticsKnnInexact = treetopDiagnosticsKnn %>% filter(distance > 1E-6)
  #treetopDiagnosticsInexact = treetopDiagnostics[treetopDiagnosticsKnnInexact$diagnosticIndex] # sorts into accepted order
  #
  #acceptedNotInDiagnostics = subset(acceptedTreetops, treetopDiagnosticsKnn$distance > 1E-6)
  #acceptedNotInDiagnostics$inEqualHeightPatch = treetopDiagnosticsInexact$inEqualHeightPatch
  #
  #writeVector(acceptedNotInDiagnostics, file.path(dataPath, paste0(tileName, " accepted treetops.gpkg")), layer = "accepted not in diagnostics", insert = TRUE, overwrite = TRUE)
  #
  #ggplot() +
  #  geom_histogram(aes(x = distance), treetopDiagnosticsKnn %>% filter(distance > 1E-6), binwidth = 1) +
  #  labs(x = "accepted treetop to nearest recorded local maxima, m", y = "treetops")
}
