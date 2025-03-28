# some investigatory code here assumes segmentation/setup.R
library(dplyr)
library(FNN)
library(furrr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(progressr)
library(ranger)
library(sf)
library(stringr)
library(terra)
library(tidyr)
library(WeightedROC)

options(cli.progress_format_iterator = "{cli::pb_bar} {cli::pb_percent} | {cli::pb_current} of {cli::pb_total} increments completed, {cli::pb_elapsed_clock} elapsed, {prettyunits::pretty_sec(as.numeric(cli::pb_eta_raw))} remaining",
        future.globals.maxSize = 1 * 1024^3) # increase from default 500 MB to 1 GB for fit_ranger_treetop()
theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             legend.title = element_text(size = 10),
                             panel.border = element_blank(), 
                             plot.subtitle = element_text(size = 9),
                             plot.title = element_text(size = 10, vjust = 0.5),
                             plot.title.position = "plot"))
plotLetters = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)")

# ranger performance maxima
# Zen 3 + DDR4-3200: one thread per core (default of two threads per core is slower and bogs the UX)
# Zen 5 + DDR5-5600: one thread per core
treetopOptions = tibble(fitRandomForest = FALSE,
                        includeInvestigatory = FALSE,
                        includeSetup = FALSE,
                        folds = 2,
                        repetitions = 25,
                        neighborhoodBufferWidth = 25, # in CRS units, so feet
                        rangerThreads = 0.5 * future::availableCores(),
                        dsmCellSize = 1.5, # feet
                        tileSize = 3000) # ft

acceptedTreetopsDsmPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted"
acceptedTreetopsChmPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted chm"
acceptedTreetopsCmmPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted cmm"
dsmPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3"
localMaximaPathV3beta = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta/local maxima"
localMaximaPathV3 = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3/local maxima"
candidateTreetopsDsmPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops"
tileCrs = NULL

extend_confusion_matrix_to_string = function(classCounts)
{
  classCountMatrix = as.matrix(classCounts)
  classes = nrow(classCountMatrix)
  dimensionNames = list(predicted = c(rownames(classCountMatrix), "n", "user's accuracy"),
                        reference = c(colnames(classCountMatrix), "n", "producer's accuracy"))
  columwiseSums = colSums(classCountMatrix)
  rowwiseSums = rowSums(classCountMatrix)
  n = sum(rowwiseSums)
  #confuMatrix = matrix(nrow = classes + 2, ncol = classes + 2, dimnames = dimensionNames)
  #confuMatrix[1:classes, 1:classes] = classCountMatrix
  #confuMatrix[1:classes, "n"] = rowwiseSums
  #confuMatrix[1:classes, "producer's accuracy"] = diag(classCountMatrix) / confuMatrix[1:classes, "n"]
  #confuMatrix["n", 1:classes] = columwiseSums
  #confuMatrix["user's accuracy", 1:classes] = diag(classCountMatrix) / confuMatrix["n", 1:classes]
  #confuMatrix["n", "n"] = n
  #confuMatrix["user's accuracy", "producer's accuracy"] = sum(diag(classCountMatrix)) / n
  #confuMatrix[1:classes, 1:classes] = 100 * confuMatrix[1:classes, 1:classes] / n
  confusionMatrixAsString = matrix(data = "", nrow = classes + 2, ncol = classes + 2, dimnames = dimensionNames)
  # for now, assume percentage of area is the same as percentage of reference distribution
  confusionMatrixAsString[1:classes, 1:classes] = sprintf("%.2g", 100 * classCountMatrix / n)
  confusionMatrixAsString[1:classes, "n"] = sprintf("%.0f", rowwiseSums)
  confusionMatrixAsString[1:classes, "producer's accuracy"] = str_replace(sprintf("%.3f", diag(classCountMatrix) / rowwiseSums), "NaN", "")
  confusionMatrixAsString["n", 1:classes] = sprintf("%.0f", columwiseSums)
  confusionMatrixAsString["user's accuracy", 1:classes] = str_replace(sprintf("%.3f", diag(classCountMatrix) / columwiseSums), "NaN", "")
  confusionMatrixAsString["n", "n"] = sprintf("%.0f", n)
  confusionMatrixAsString["user's accuracy", "producer's accuracy"] = sprintf("%.3f", sum(diag(classCountMatrix)) / n)
  return(confusionMatrixAsString)
}

fit_ranger_treetop = function(trainingMaxima, neighborhoodMaxima = trainingMaxima, mtry, minNodeSize, sampleFraction, classWeights = NULL, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
{
  progressBar = progressor(steps = folds * repetitions)
  
  if (folds == 1)
  {
    return(bind_rows(future_map(1:repetitions, function(crossValidationRepetition)
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
      return(allFitMetrics %>% mutate(repetition = crossValidationRepetition, fold = 1))
    })))
  }
  
  # use random cross validation as the training dataset lacks the spatial extent to block by stands
  # ranger is expected to use all cores, so little advantage to future_map() rather than map() is assumed
  uniqueMergeClusterIDs = tibble(clusterID = unique(trainingMaxima$uniqueMergeClusterID))

  splitsAndFits = bind_rows(future_map(1:repetitions, function(crossValidationRepetition)
  {
    clusterFoldAssignments = uniqueMergeClusterIDs %>% mutate(fold = rep(1:folds, each = ceiling(n() / folds))[sample.int(n = folds * ceiling(n() / folds), size = n())])
    return(bind_rows(map(1:folds, function(crossValidationFold)
    {
      start = Sys.time()
      currentFoldIndices = which(clusterFoldAssignments$fold == crossValidationFold)
      trainingData = trainingMaxima %>% filter(uniqueMergeClusterID %in% uniqueMergeClusterIDs$clusterID[currentFoldIndices])
      
      randomForestFit = ranger(treetop ~ ., data = trainingData %>% select(-all_of(accuracyVariables)), mtry = mtry, min.node.size = minNodeSize, sample.fraction = sampleFraction, class.weights = classWeights, num.threads = treetopOptions$rangerThreads)
      
      otherFoldIndices = which(clusterFoldAssignments$fold != crossValidationFold)
      validationData = trainingMaxima %>% filter(uniqueMergeClusterID %in% uniqueMergeClusterIDs$clusterID[otherFoldIndices])
      
      validationPrediction = predict(randomForestFit, validationData %>% select(-treetop))$predictions
      tileMergePoints = get_merge_points(validationData, neighborhoodMaxima, treetopClassification = validationPrediction)
      validationPrediction[tileMergePoints$mergePointIndices] = "merge"
      validationPrediction[tileMergePoints$treetopIndices] = "yes"
      fitMetrics = get_treetop_accuracy(validationPrediction, validationData) %>%
        mutate(repetition = crossValidationRepetition,
               fold = crossValidationFold,
               nTrain = nrow(trainingData),
               fitTimeInS = as.numeric(difftime(Sys.time(), start, units = "secs")))
      
      progressBar()
      return(fitMetrics)
    })))
  }, .options = furrr_options(seed = TRUE)))
  return(splitsAndFits %>% relocate(repetition, fold, nTrain))
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
  #   Max displacement + misalignment noted: s04200w06840 3.5 m, s04200w06810 4.0 m
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
  # neighborhoodMaxima = tileNeighborhood
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
  # TODO: Exclude noise points from merge clustering as well as from the initiating kNN above? For now they're excluded from multipoint 
  # clusters post hoc as it's unclear whether connectivity through noise points during cluster formation is more of a feature or more of a bug.
  mergeClusterInitiatingPoints = bind_rows(mergePointKnn %>% filter(neighbors > 0),
                                           treetopMergeKnn %>% filter(neighbors > 0))
  mergePointUniqueIDs = unique(c(mergeClusterInitiatingPoints$uniqueID,
                                 na.omit(mergeClusterInitiatingPoints$neighbor1uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor2uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor3uniqueID),
                                 na.omit(mergeClusterInitiatingPoints$neighbor4uniqueID)))
  addedMergePointUniqueIDs = setdiff(mergePointUniqueIDs, mergeClusterInitiatingPoints$uniqueID)

  mergePoints = bind_rows(mergeClusterInitiatingPoints %>% mutate(isInitiating = TRUE), 
                          neighborhoodMaxima %>% filter(uniqueID %in% addedMergePointUniqueIDs) %>% mutate(isInitiating = FALSE)) %>%
    mutate(clusterPoints = 0)

  adjacencyMatrixFrame = pivot_longer(mergePoints %>% select(uniqueID, neighbor1uniqueID, neighbor2uniqueID, neighbor3uniqueID, neighbor4uniqueID) %>% mutate(neighbor0uniqueID = uniqueID) %>% rename(mergePointUniqueID = uniqueID),
                                      cols = c("neighbor0uniqueID", "neighbor1uniqueID", "neighbor2uniqueID", "neighbor3uniqueID", "neighbor4uniqueID"),
                                      names_pattern = "neighbor(.)uniqueID", names_to = "neighbor",
                                      values_to = "clusterMemberUniqueID") %>%
    filter(is.na(clusterMemberUniqueID) == FALSE) %>%
    mutate(row = match(mergePointUniqueID, mergePoints$uniqueID), column = match(clusterMemberUniqueID, mergePoints$uniqueID))
  adjacencyMatrix = Matrix::sparseMatrix(i = adjacencyMatrixFrame$row, j = adjacencyMatrixFrame$column, x = TRUE, dimnames = list(mergePoints$uniqueID, mergePoints$uniqueID))
  connectivityMatrix = Matrix::tcrossprod(adjacencyMatrix, boolArith = TRUE) > 0
  mergeClusterComponents = igraph::components(igraph::graph_from_adjacency_matrix(connectivityMatrix))
  
  tileNames = unique(tileMaxima$tile)
  if (mergeClusterComponents$no > 0)
  {
    mergePoints = left_join(mergePoints %>% filter(treetop %in% c("yes", "merge", "no")), # allow clusters to include random forest single tops and non-tops but exclude noise points
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
      # TODO: refine outlier removal?
      # If a noise point's included in a cluster it is typically due it to not being classified as such.
      filter((dsmZ >= (clusterMaxDsmZ - clusterThickness)) & (dsmZ <= clusterMaxDsmZ), # within height range
             clusterDistance < 2 * clusterMeanDistance + 1E-6) %>% # within distance of centroid, epsilon keeps single points from excluding themselves (clusterDistance = clusterMeanDistance = 0)
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
      summarize(treeID = min(if_else(tile == centroidTile[1], id, NA_integer_), na.rm = TRUE), # all clusters should have at least one merge point in tileMaxima
                treetop = factor("merge treetop"),
                mergePoints = n(),
                mergePointsOnTile = sum(tile == centroidTile[1]), # always 1 if calculated after tile is assigned
                tile = centroidTile[1],
                x = clusterCentroidX[1],
                y = clusterCentroidY[1],
                radius = max(radius), dsmZ = clusterCentroidDsmZ[1], cmmZ = mean(cmmZ), height = mean(height),
                # leave ring statistics unpopulated as they're currently used only from local maxima and are not well defined for treetops obtained from merge clusters
                # sourceID = if_else(length(unique(sourceID)) == 1, sourceID[1], as.integer(0)), # source ID variety in merge clusters not currently used
                # sourceIDs = length(unique(sourceID)), 
                .groups = "drop") %>%
      select(-mergeClusterNumber) # potentially useful to flow for investigation and debugging but not functionally needed by any caller
  } else {
    mergePointsInOnTileClusters = mergePoints %>% filter(FALSE) # should be an empty tibble
    treetopsFromMergePoints = tibble(tile = character(), treeID = integer(), x = numeric(), y = numeric(), 
                                     radius = numeric(), dsmZ = numeric(), cmmZ = numeric(), height = numeric(), 
                                     mergePoints = integer(), mergePointsOnTile = integer())
  }

  # identify merge points which need reclassification as treetops because they're in single point clusters
  # Treetops with neighbors in treetopMergeKnn but which do not aggregate into multi-point clusters do not need to be listed here
  # as they're already classified as treetops.
  # Singleton merge points which do not appear to merit conversion to treetops could be rejected here.
  uniqueIDsWithoutInitialNeighbors = (mergePointKnn %>% filter(tile %in% tileNames, neighbors == 0))$uniqueID # may have gotten neighbors through connectivity matrix
  uniqueIDsInSinglePointClusters = (mergePoints %>% filter(tile %in% tileNames, clusterPoints == 1, treetop == "merge"))$uniqueID # initiating merge points with in range neighbors rejected by cluster filtering
  uniqueIDsInMultipointClusters = (mergePoints %>% filter(tile %in% tileNames, clusterPoints > 1, treetop %in% c("yes", "merge", "no")))$uniqueID
  #tibble(withoutInitial = length(uniqueIDsWithoutInitialNeighbors), singlePoint = length(uniqueIDsInSinglePointClusters), multipoint = length(uniqueIDsInMultipointClusters), withoutInitialLinkedIntoMultipoint = length(intersect(uniqueIDsWithoutInitialNeighbors, uniqueIDsInMultipointClusters)))
  
  singleTreetops = tibble(uniqueID = setdiff(union(uniqueIDsWithoutInitialNeighbors, uniqueIDsInSinglePointClusters), uniqueIDsInMultipointClusters)) %>% 
    mutate(tileIndex = match(uniqueID, tileMaxima$uniqueID))

  # identify on tile local maxima which are merge points so that classifying callers update non-merge points that were clustered in and set cluster IDs
  # Contains NAs for merge points which are in neighborhoodMaxima but not in tileMaxima.
  mergePointsInMultipointClusters = mergePoints %>% filter(uniqueID %in% uniqueIDsInMultipointClusters) %>%
    mutate(tileIndex = match(uniqueID, tileMaxima$uniqueID)) %>%
    select(tile, id, x, y, mergeClusterNumber, tileIndex, treetop) # , clusterPoints) # if additional classes besides noise are excluded from final cluster definitions a few merge clusters will contain fewer points that clusterPoints

  singleTopAndMergePointTileIndices = intersect(singleTreetops$tileIndex, mergePointsInMultipointClusters$tileIndex)
  if (length(singleTopAndMergePointTileIndices) > 0)
  {
    stop(paste0("Internal consistency failure on tile ", tileMaxima$tile[1], ". ", length(singleTopAndMergePointTileIndices), " local maxima are both in merge clusters and identified as single tops."))
  }
  
  return(list(mergeTreetops = treetopsFromMergePoints, mergePoints = mergePointsInMultipointClusters, singleTreetops = singleTreetops))
}

get_merge_distance = function(tileMaxima, treetopMerge = FALSE)
{
  # TODO: include radius penalty for same source ID?
  return(treetopOptions$dsmCellSize * (1 + 1/60 * tileMaxima$height + if_else(tileMaxima$height < 125, 0, 1/40 * (tileMaxima$height - 125))))
  #if (treetopMerge)
  #{
  #  return(treetopOptions$dsmCellSize * (1 + 1/60 * tileMaxima$height + if_else(tileMaxima$height < 125, 0, 1/40 * (tileMaxima$height - 125))))
  #} else {
  #  verticalExaggeration = 25 # DSM multiplier
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
  confusionMatrix = caret::confusionMatrix(predictedBinary, expectedBinary, positive = "yes", mode = "everything")
  confusionSubmatrix = caret::confusionMatrix(predicted, validationData$treetop, mode = "everything")
  
  classAccuracy = diag(confusionMatrix$table) / rowSums(confusionMatrix$table)
  subclassAccuracy = diag(confusionSubmatrix$table) / rowSums(confusionSubmatrix$table)
  classCounts = table(predicted)
  expectedClassCounts = table(validationData$treetop)
  
  overallAccuracyByHeightClass = validationData %>% select(treetop, height) %>% 
    mutate(heightClass = round(height),
           predictedClass = predicted,
           expectedBinary = expectedBinary,
           predictedBinary = predictedBinary) %>%
    group_by(heightClass) %>%
    summarize(n = n(), 
              overallAccuracy = sum(treetop == predictedClass) / n(),
              treetopAccuracy = sum(expectedBinary == predictedBinary) / n(),
              .groups = "drop")
  
  # can't readily compute distance and height error on mismatches here
  # - with random k-fold cross validation the status, on average, of s/k nearby maxima is unknown (where s is the tile sampling fraction)
  # - with spatially blocked cross validation the window to use to exclude edge effects on the validation data must be known
  return(tibble(nValidate = nrow(validationData),
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
                confusionSubmatrix = list(confusionSubmatrix),
                overallAccuracyByHeight = list(overallAccuracyByHeightClass)))
}

# default to minimum height of 1 m
get_treetop_eligible_maxima = function(tileName, localMaximaLayer = "localMaximaDsm", minimumHeightInM = 1, acceptedTileName = NULL, localMaximaPath = localMaximaPathV3, acceptedTreetopsPath = acceptedTreetopsDsmPath, returnNullIfNoFile = FALSE)
{
  #tileName = "s04200w06840"
  #acceptedTileName = "s04200w06840 cmm"
  #acceptedTreetopsPath = acceptedTreetopsCmmPath
  #localMaximaLayer = "localMaximaCmm"
  #localMaximaPath = localMaximaPathV3
  
  localMaximaFilePath = file.path(localMaximaPath, paste0(tileName, ".gpkg"))
  if (file.exists(localMaximaFilePath) == FALSE)
  {
    if (returnNullIfNoFile)
    {
      return(NULL)
    } else {
      stop(paste0("Requested local maxima file '", localMaximaFilePath, "' does not exist."))
    }
  }
  localMaxima = st_read(localMaximaFilePath, layer = localMaximaLayer, quiet = TRUE)
  
  # use minimumHeight to exclude groundcover maxima and maxima from sensor noise or error
  # Height is still in English units at this point so convert minimumHeight from metric.
  # For now, also exclude treetop candidates with so few adjacent points ring 1 or ring 2 has no data since it's 1) it's difficult to tell if these are maxima, 2) it's unlikely they're actually maxima, and 3) these likely comprise < 0.01% of maximas.
  localMaximaCrsParameters = st_crs(localMaxima, parameters = TRUE)
  minimumHeightInCrsUnits = minimumHeightInM
  if (localMaximaCrsParameters$units_gdal == "foot")
  {
    minimumHeightInCrsUnits = 3.28084 * minimumHeightInM
  }
  
  localMaximaExcluded = localMaxima %>% filter((height < minimumHeightInCrsUnits) | is.na(localMaxima$ring1mean) | is.na(localMaxima$ring2mean))
  localMaxima %<>% filter((height >= minimumHeightInCrsUnits) & (is.na(ring1mean) == FALSE) & (is.na(ring2mean) == FALSE))
  
  # convert to metric
  if (localMaximaCrsParameters$units_gdal == "foot")
  {
    localMaxima$radius = 0.3048 * localMaxima$radius
    localMaxima$dsmZ = 0.3048 * localMaxima$dsmZ
    localMaxima$cmmZ = 0.3048 * localMaxima$cmmZ
    localMaxima$height = 0.3048 * localMaxima$height
    localMaxima$ring1max = 0.3048 * localMaxima$ring1max
    localMaxima$ring2max = 0.3048 * localMaxima$ring2max
    localMaxima$ring3max = 0.3048 * localMaxima$ring3max
    localMaxima$ring4max = 0.3048 * localMaxima$ring4max
    localMaxima$ring5max = 0.3048 * localMaxima$ring5max
    localMaxima$ring1mean = 0.3048 * localMaxima$ring1mean
    localMaxima$ring2mean = 0.3048 * localMaxima$ring2mean
    localMaxima$ring3mean = 0.3048 * localMaxima$ring3mean
    localMaxima$ring4mean = 0.3048 * localMaxima$ring4mean
    localMaxima$ring5mean = 0.3048 * localMaxima$ring5mean
    localMaxima$ring1min = 0.3048 * localMaxima$ring1min
    localMaxima$ring2min = 0.3048 * localMaxima$ring2min
    localMaxima$ring3min = 0.3048 * localMaxima$ring3min
    localMaxima$ring4min = 0.3048 * localMaxima$ring4min
    localMaxima$ring5min = 0.3048 * localMaxima$ring5min
    localMaxima$ring1variance = 0.3048^2 * localMaxima$ring1variance
    localMaxima$ring2variance = 0.3048^2 * localMaxima$ring2variance
    localMaxima$ring3variance = 0.3048^2 * localMaxima$ring3variance
    localMaxima$ring4variance = 0.3048^2 * localMaxima$ring4variance
    localMaxima$ring5variance = 0.3048^2 * localMaxima$ring5variance
  }
  
  # default every local maxima as a singleton, overridden below if an accepted tile with merge information is available
  localMaxima$mergeClusterID = localMaxima$id
  
  slopeAspect = rast(file.path(localMaximaPath, "../slopeAspect", paste0(tileName, ".tif")))
  localMaximaSlopeAspect = terra::extract(slopeAspect, localMaxima, method = "simple") # nearest neighbor
  localMaxima$dsmSlope = localMaximaSlopeAspect$dsmSlope
  localMaxima$cmmSlope3 = localMaximaSlopeAspect$cmmSlope3
  
  localMaximaCrs = crs(localMaxima) # capture CRS for callers needing to rebuild sf geometry later
  localMaximaXY = st_coordinates(localMaxima)[, c("X", "Y")]
  localMaxima = as_tibble(localMaxima) %>% 
    mutate(x = localMaximaXY[, "X"], 
           y = localMaximaXY[, "Y"], 
           uniqueID = 1000000 * as.integer(str_c(str_sub(tileName, 2, 6), str_sub(tileName, 8, 12))) + localMaxima$id,
           uniqueMergeClusterID = 1000000 * as.integer(str_c(str_sub(tileName, 2, 6), str_sub(tileName, 8, 12))) + localMaxima$mergeClusterID)
  attributes(localMaxima)$crs = localMaximaCrs
  
  if (is.null(acceptedTileName))
  {
    return(localMaxima)
  }
  
  acceptedTreetopsFilePath = file.path(acceptedTreetopsPath, paste0(acceptedTileName, ".gpkg"))
  if (file.exists(acceptedTreetopsFilePath) == FALSE)
  {
    stop(paste0("Accepted treetop GeoPackage '", acceptedTreetopsFilePath, "' does not exist."))
  }
  
  localMaxima$treetop = factor("no", levels = c("no", "yes", "merge", "noise", "maybe noise"))
  
  availableTruthLayers = st_layers(acceptedTreetopsFilePath)
  if ("treetops" %in% availableTruthLayers$name)
  {
    acceptedTreetops = st_read(acceptedTreetopsFilePath, layer = "treetops", quiet = TRUE)
    acceptedTreetopsBelowMinimumHeight = acceptedTreetops %>% filter(height < minimumHeightInCrsUnits)
    acceptedTreetops %<>% filter(height >= minimumHeightInCrsUnits)
    isTreetopKnn = get.knnx(localMaximaXY, st_coordinates(acceptedTreetops)[, c("X", "Y")], k = 1)
    isTreetopKnn = tibble(index = isTreetopKnn$nn.index[, 1], distance = isTreetopKnn$nn.dist[, 1]) %>% filter(distance < 0.1) # assume any top farther away is a merge top and thus not to be marked as a single top
    localMaxima$treetop[isTreetopKnn$index] = "yes"
  }
  
  hasMergePoints = "merge points" %in% availableTruthLayers$name
  hasMergeTreetops = "merge treetops" %in% availableTruthLayers$name # backwards compatibility
  if (hasMergePoints | hasMergeTreetops)
  {
    if (hasMergePoints & hasMergeTreetops)
    {
      stop(paste0("'", acceptedTreetopsFilePath, "' has a 'merge points' layer as well as a 'merge treetops' layer. Don't know which to use."))
    }

    # for now, assume dataset consistency checking ensures consistent merge clustering
    # Cross checking of merge points against treetops distanced from nearest local maxima can be implemented if needed.
    mergePoints = st_read(acceptedTreetopsFilePath, layer = if_else(hasMergePoints, "merge points", "merge treetops"), quiet = TRUE)
    if ((is.null(acceptedTreetopsBelowMinimumHeight) == FALSE) & (nrow(acceptedTreetopsBelowMinimumHeight) > 0))
    {
      treesBelowMinimumHeightByID = unique(acceptedTreetopsBelowMinimumHeight$treeID)
      mergePoints %<>% filter((clusterID %in% treesBelowMinimumHeightByID) == FALSE)
    }
    
    isMergeKnn = get.knnx(localMaximaXY, st_coordinates(mergePoints)[, c("X", "Y")], k = 1)
    isMergeKnn = tibble(index = isMergeKnn$nn.index[, 1], distance = isMergeKnn$nn.dist[, 1], clusterID = mergePoints$clusterID) 
    unsnappedMergePoints = sum(isMergeKnn$distance >= 0.1)
    if (unsnappedMergePoints > 0)
    {
      warning(paste0(unsnappedMergePoints, " merge points more than 0.1 ft from nearest retained local maxima. This may be due to a member of a merge cluster being below the specified minimum height of ", minimumHeightInM, " m."))
      #table(isMergeKnn$distance)
      #mergePoints[which(isMergeKnn$distance >= 0.1), ]
      #ggplot() +
      #  geom_histogram(aes(x = distance), isMergeKnn, binwidth = 1)
    }
    
    # drop merge points distanced from nearest local maxima
    # This is expected only to drop merge points on local maxima below minimumHeight.
    isMergeKnn %<>% filter(distance < 0.1)
    
    localMaxima$treetop[isMergeKnn$index] = "merge"
    localMaxima$mergeClusterID[isMergeKnn$index] = isMergeKnn$clusterID
    # BUGBUG: need a way for unique merge cluster IDs to extend across tiles instead of assigning multiple unique IDs to the cluster
    localMaxima$uniqueMergeClusterID[isMergeKnn$index] = 1000000 * as.integer(str_c(str_sub(tileName, 2, 6), str_sub(tileName, 8, 12))) + isMergeKnn$clusterID
  }
  if ("noise points" %in% availableTruthLayers$name)
  {
    noisePoints = st_read(acceptedTreetopsFilePath, layer = "noise points", quiet = TRUE)
    noisePointKnn = get.knnx(localMaximaXY, st_coordinates(noisePoints)[, c("X", "Y")], k = 1)
    noisePointKnn = tibble(index = noisePointKnn$nn.index[, 1], distance = noisePointKnn$nn.dist[, 1])
    unsnappedNoisePoints = sum(noisePointKnn$distance >= 0.1)
    if (unsnappedNoisePoints > 0)
    {
      warning(paste0(unsnappedNoisePoints, " noise points more than 0.1 ft from nearest retained local maxima."))
      #which(noisePointKnn$distance >= 0.1)
    }
    
    # mark local maxima as noise except where noise points lie over singleton tops or merge points
    localMaxima$treetop[noisePointKnn$index] = if_else(localMaxima$treetop[noisePointKnn$index] == "no", "noise", localMaxima$treetop[noisePointKnn$index])
  }
  if ("maybe noise points" %in% availableTruthLayers$name)
  {
    maybeNoisePoints = st_read(acceptedTreetopsFilePath, layer = "maybe noise points", quiet = TRUE)
    if (nrow(maybeNoisePoints) > 0)
    {
      maybeNoisePointsXY = st_coordinates(maybeNoisePoints)[, c("X", "Y")]
      if (nrow(maybeNoisePoints) == 1)
      {
        # matrix indexing collapses a single row to a vector, requiring conversion back to a matrix
        maybeNoisePointsXY = matrix(maybeNoisePointsXY, ncol = 2)
      }
      maybeNoisePointKnn = get.knnx(localMaximaXY, maybeNoisePointsXY, k = 1)
      maybeNoisePointKnn = tibble(index = maybeNoisePointKnn$nn.index[, 1], distance = maybeNoisePointKnn$nn.dist[, 1])
      unsnappedMaybeNoisePoints = sum(maybeNoisePointKnn$distance >= 0.1)
      if (unsnappedMaybeNoisePoints > 0)
      {
        warning(paste0(unsnappedMaybeNoisePoints, " maybe noise points more than 0.1 ft from nearest retained local maxima."))
      }
      
      localMaxima$treetop[maybeNoisePointKnn$index] = "maybe noise" # for now, override any other classification
    }
  }
  
  return(localMaxima)
}

get_treetop_eligible_neighborhood = function(tileName, tileMaxima, localMaximaPath = localMaximaPathV3)
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
    northwestMaxima = get_treetop_eligible_maxima(tileNameNorthwest, localMaximaPath = localMaximaPath) %>% 
      filter(x >= tileExtentBuffered$xMin, y <= tileExtentBuffered$yMax)
  }
  tileNameNorth = sprintf("s%05dw%05d", xTile, yTile + 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameNorth, ".gpkg"))))
  {
    northMaxima = get_treetop_eligible_maxima(tileNameNorth, localMaximaPath = localMaximaPath) %>% 
      filter(y <= tileExtentBuffered$yMax)
  }
  tileNameNortheast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile + 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameNortheast, ".gpkg"))))
  {
    northeastMaxima = get_treetop_eligible_maxima(tileNameNortheast, localMaximaPath = localMaximaPath) %>% 
      filter(x <= tileExtentBuffered$xMax, y <= tileExtentBuffered$yMax)
  }
  
  tileNameWest = sprintf("s%05dw%05d", xTile - 0.01 * treetopOptions$tileSize, yTile)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameWest, ".gpkg"))))
  {
    westMaxima = get_treetop_eligible_maxima(tileNameWest, localMaximaPath = localMaximaPath) %>% 
      filter(x >= tileExtentBuffered$xMin)
  }
  tileNameEast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameEast, ".gpkg"))))
  {
    eastMaxima = get_treetop_eligible_maxima(tileNameEast, localMaximaPath = localMaximaPath) %>% 
      filter(x <= tileExtentBuffered$xMax)
  }
  
  tileNameSouthwest = sprintf("s%05dw%05d", xTile - 0.01 * treetopOptions$tileSize, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSouthwest, ".gpkg"))))
  {
    southwestMaxima = get_treetop_eligible_maxima(tileNameSouthwest, localMaximaPath = localMaximaPath) %>% 
      filter(x >= tileExtentBuffered$xMin, y >= tileExtentBuffered$yMin)
  }
  tileNameSouth = sprintf("s%05dw%05d", xTile, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSouth, ".gpkg"))))
  {
    southMaxima = get_treetop_eligible_maxima(tileNameSouth, localMaximaPath = localMaximaPath) %>% 
      filter(y >= tileExtentBuffered$yMin)
  }
  tileNameSoutheast = sprintf("s%05dw%05d", xTile + 0.01 * treetopOptions$tileSize, yTile - 0.01 * treetopOptions$tileSize)
  if (file.exists(file.path(localMaximaPath, paste0(tileNameSoutheast, ".gpkg"))))
  {
    southeastMaxima = get_treetop_eligible_maxima(tileNameSoutheast, localMaximaPath = localMaximaPath) %>% 
      filter(x <= tileExtentBuffered$xMax, y >= tileExtentBuffered$yMin)
  }
  
  return(bind_rows(northwestMaxima, northMaxima, northeastMaxima, eastMaxima, tileMaxima, westMaxima, southwestMaxima, southMaxima, southeastMaxima))
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
           neighborDistance2mean = neighborKnn$meanDistance2, neighborDistance3mean = neighborKnn$meanDistance3, neighborDistance4mean = neighborKnn$meanDistance4, 
           neighborDistance5mean = neighborKnn$meanDistance5, neighborDistance10mean = neighborKnn$meanDistance10, neighborDistance20mean = neighborKnn$meanDistance20, neighborDistance50mean = neighborKnn$meanDistance50,
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
           neighborDistance5variance = 0.2 * ((neighbor1distance - neighborDistance5mean)^2 + (neighbor2distance - neighborDistance5mean)^2 + (neighbor3distance - neighborDistance5mean)^2 + (neighbor4distance - neighborDistance5mean)^2 + (neighbor5distance - neighborDistance5mean)^2),
           neighborProminence5mean = 0.2 * (neighbor1prominence + neighbor2prominence + neighbor3prominence + neighbor4prominence + neighbor5prominence),
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
           prominenceNeighbor5Variance = 0.2 * ((neighbor1prominence - neighborProminence5mean)^2 + (neighbor2prominence - neighborProminence5mean)^2 + (neighbor3prominence - neighborProminence5mean)^2 + (neighbor4prominence - neighborProminence5mean)^2 + (neighbor5prominence - neighborProminence5mean)^2),
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
           slope1 = 180 / pi * atan2(range1, 2 * treetopOptions$dsmCellSize),
           slope2 = 180 / pi * atan2(range2, 4 * treetopOptions$dsmCellSize),
           slope3 = 180 / pi * atan2(range3, 6 * treetopOptions$dsmCellSize),
           slope4 = 180 / pi * atan2(range4, 8 * treetopOptions$dsmCellSize),
           slope5 = 180 / pi * atan2(range5, 10 * treetopOptions$dsmCellSize),
           slope1normalized = slope1 / height,
           slope2normalized = slope2 / height,
           slope3normalized = slope3 / height,
           slope4normalized = slope4 / height,
           slope5normalized = slope5 / height,
           netProminence = prominence1 + prominence2 + prominence3 + prominence4 + prominence5,
           netProminenceNeighbor2 = neighbor1prominence + neighbor2prominence,
           netProminenceNeighbor3 = neighbor1prominence + neighbor2prominence + neighbor3prominence,
           netProminenceNeighbor5 = neighbor1prominence + neighbor2prominence + neighbor3prominence + neighbor4prominence + neighbor5prominence,
           netProminenceNormalized = netProminence / height,
           netProminenceNeighbor1normalized = (neighbor1prominence / neighbor1distance) / (height / neighbor1distance),
           netProminenceNeighbor2normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance)),
           netProminenceNeighbor3normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance)),
           netProminenceNeighbor4normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance + neighbor4prominence / neighbor4distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance + 1 / neighbor4distance)),
           netProminenceNeighbor5normalized = (neighbor1prominence / neighbor1distance + neighbor2prominence / neighbor2distance + neighbor3prominence / neighbor3distance + neighbor4prominence / neighbor4distance + neighbor5prominence / neighbor5distance) / (height * (1 / neighbor1distance + 1 / neighbor2distance + 1 / neighbor3distance + 1 / neighbor4distance + 1 / neighbor5distance)),
           netRange = range1 + range2 + range3 + range4 + range5,
           netRangeNormalized = netRange / height,
           varianceNormalized1 = ring1variance / height, 
           varianceNormalized2 = ring2variance / height, 
           varianceNormalized3 = ring3variance / height, 
           varianceNormalized4 = ring4variance / height, 
           varianceNormalized5 = ring5variance / height)

  return(treetopPredictors)
}

make_compound_crs = function(projectedEpsg, verticalEpsg)
{
  projectedCrs = st_crs(projectedEpsg)
  projectedCrsName = str_extract(projectedCrs$wkt, "PROJCRS\\[\"([A-Za-z0-9_\\[\\]\\(\\){}<=>\\.,:;\\+\\- #%&'*^/\\@|°]+)", group = 1) # <quoted Latin text> except <doublequote symbol>, https://docs.ogc.org/is/18-010r11/18-010r11.pdf
  verticalCrs = st_crs(verticalEpsg)
  verticalCrsName = str_extract(verticalCrs$wkt, "VERTCRS\\[\"([A-Za-z0-9_\\[\\]\\(\\){}<=>\\.,:;\\+\\- #%&'*^/\\@|°]+)", group = 1)
  
  compoundWkt = paste0("COMPOUNDCRS[\"", projectedCrsName, " + ", verticalCrsName, "\",\n    ", projectedCrs$wkt, ",\n    ", verticalCrs$wkt, "]")
  return(st_crs(compoundWkt))
}


## treetop dataset + surrounding context
# 15 tile load (2.6M maxima, single threaded): ~2m 35s 9900X + DDR5-5600
# number of predictor variables = ncol(s4268maximaDsm) - length(c("tile", "id", "uniqueID", "treetop", "x", "y"))
if (treetopOptions$includeSetup)
{
  loadStart = Sys.time()
  s4268maximaDsm = get_treetop_predictors(bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"), # truthed tiles
                                                    get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"), 
                                                    get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04170w06870", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"), # adjacent tiles for neighborhood calculations
                                                    get_treetop_eligible_maxima("s04200w06870", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"), # 15 warnings on z coordinate loss on read, one for each tile
                                                    get_treetop_eligible_maxima("s04230w06870", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"), # 2 additional warnings from s04200w06840 and s04230w06810 truth
                                                    get_treetop_eligible_maxima("s04170w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04230w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04260w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04170w06810", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04260w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04170w06780", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04200w06780", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04230w06780", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                                                    get_treetop_eligible_maxima("s04260w06780", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"))) %>% 
    filter(is.na(dsmSlope) == FALSE, 
           ((420000 - treetopOptions$neighborhoodBufferWidth) < x) & (x < (426000 + treetopOptions$neighborhoodBufferWidth)), ((681000 - treetopOptions$neighborhoodBufferWidth) < y) & (y < (687000 + treetopOptions$neighborhoodBufferWidth))) %>% # window neighboring tiles, all of s04230w06840 is still included
    select(-neighbor1sourceID, -neighbor2sourceID, -neighbor3sourceID, -neighbor4sourceID, -neighbor5sourceID, # drop most non-portable values (id, sourceID, x, y, dsmZ, and cmmZ are excluded from training by predictor variable selection)
           -ring1max, -ring2max, -ring3max, -ring4max, -ring5max, -ring1mean, -ring2mean, -ring3mean, -ring4mean, -ring5mean, -ring1min, -ring2min, -ring3min, -ring4min, -ring5min, # drop elevations
           -neighbor1dsmZ, -neighbor2dsmZ, -neighbor3dsmZ, -neighbor4dsmZ, -neighbor5dsmZ) 
  s4268maximaDsm %<>% filter(is.na(dsmSlope) == FALSE) # exclude 346 local maxima on accepted tiles whose local slopes are undefined (comment this out if DSM and CMM slope aren't being considered as predictors)
  Sys.time() - loadStart
  
  # drop unrelated predictors
  s4268maximaDsm %<>% select(-neighbor4differentSourceID, -neighbor5differentSourceID) # Boruta
  
  tibble(localMaximaInDataset = sum(is.na(s4268maximaDsm$treetop) == FALSE), incompleteCases = sum(is.na(s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(-dsmSlope)))) # check for incomplete cases in training data
  #print(bind_rows(colSums(is.na(s4268maximaDsm))) %>% pivot_longer(cols = everything()), n = 125) # NA treetops expected on tiles without accepted treetops
  #s4268maximaDsm %>% filter(neighbor1distance == 0) %>% reframe(tile = unique(tile)) # TODO: investigate zero distance neighbors on s04260w06840
  #colSums(is.na(s4268maximaDsm)))
  #(sum(s4268maximaDsm$neighbor1distance <= s4268maximaDsm$neighbor2distance) + sum(s4268maximaDsm$neighbor2distance <= s4268maximaDsm$neighbor3distance) + sum(s4268maximaDsm$neighbor3distance <= s4268maximaDsm$neighbor4distance) + sum(s4268maximaDsm$neighbor4distance <= s4268maximaDsm$neighbor5distance)) / nrow(s4268maximaDsm) # check neighbor distance sort; should be exactly 4
  
  # drop redundant predictors
  #s4268maximaDsm %<>% select(-mean1delta, -mean2delta, -mean3delta, -mean4delta, -mean5delta,
  #                           -neighbor1prominence, -neighbor2prominence, -neighbor3prominence, -neighbor4prominence, -neighbor5prominence,
  #                           # TODO: neighbor net prominence
  #                           -prominence1, -prominence2, -prominence3, -prominence4, -prominence5, -prominenceMean, -prominenceMeanNormalized, -prominenceVariance,
  #                           -range1, -range2, -range3, -range4, -range5, -rangeMean, -rangeMeanNormalized, 
  #                           -slope1normalized, -slope2normalized, -slope3normalized, -slope4normalized, -slope5normalized,
  #                           -varianceNormalized1, -varianceNormalized2, -varianceNormalized3, -varianceNormalized4, -varianceNormalized5)
  
  # DSM dataset without surrounding local maxima
  datasetDsmStart = Sys.time() # 32s
  treetopDataDsm = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>%
    select(tile, treetop, height, radius, dsmZ, mergeClusterID) %>%
    mutate(treetop = factor(treetop, levels = c("yes", "merge", "noise", "maybe noise", "no"))) %>%
    group_by(tile, mergeClusterID) %>%
    mutate(isTreetopRadius = (treetop == "yes") | ((treetop == "merge") & (dsmZ == max(dsmZ))),
           identicalElevationsInCluster = sum(isTreetopRadius),
           mergeClusterSize = if_else(identicalElevationsInCluster > 0, n(), 0)) %>%
    ungroup() %>%
    mutate(isTreetopRadius = factor(isTreetopRadius))
  Sys.time() - datasetDsmStart
  
  # CHM dataset
  datasetChmStart = Sys.time() # 36 s
  treetopDataChm = bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810 chm", acceptedTreetopsPath = acceptedTreetopsChmPath, localMaximaLayer = "localMaximaChm", localMaximaPath = localMaximaPathV3),
                             get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840 chm", acceptedTreetopsPath = acceptedTreetopsChmPath, localMaximaLayer = "localMaximaChm", localMaximaPath = localMaximaPathV3),
                             get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810 chm", acceptedTreetopsPath = acceptedTreetopsChmPath, localMaximaLayer = "localMaximaChm", localMaximaPath = localMaximaPathV3)) %>%
    mutate(treetop = factor(treetop, levels = c("yes", "merge", "noise", "maybe noise", "no"))) %>%
    group_by(tile, mergeClusterID) %>%
    mutate(isTreetopRadius = (treetop == "yes") | ((treetop == "merge") & (height == max(height))),
           identicalElevationsInCluster = sum(isTreetopRadius),
           mergeClusterSize = if_else(identicalElevationsInCluster > 0, n(), 0)) %>%
    ungroup() %>%
    mutate(isTreetopRadius = factor(isTreetopRadius))
  Sys.time() - datasetChmStart
  #table(treetopDataChm$treetop) # written to layers: 74560 tops (6166 from merges) + 7338 merge points + 2439 unmatched tops + 2299 noise points
  
  # CMM dataset
  datasetCmmStart = Sys.time() # 7 s
  treetopDataCmm = bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810 cmm", acceptedTreetopsPath = acceptedTreetopsCmmPath, localMaximaLayer = "localMaximaCmm", localMaximaPath = localMaximaPathV3),
                             get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840 cmm", acceptedTreetopsPath = acceptedTreetopsCmmPath, localMaximaLayer = "localMaximaCmm", localMaximaPath = localMaximaPathV3),
                             get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810 cmm", acceptedTreetopsPath = acceptedTreetopsCmmPath, localMaximaLayer = "localMaximaCmm", localMaximaPath = localMaximaPathV3)) %>%
    mutate(treetop = factor(treetop, levels = c("yes", "merge", "noise", "maybe noise", "no"))) %>%
    group_by(tile, mergeClusterID) %>%
    mutate(dsmHeight = height,
           height = dsmHeight + cmmZ - dsmZ,
           isTreetopRadius = (treetop == "yes") | ((treetop == "merge") & (cmmZ == max(cmmZ))),
           identicalElevationsInCluster = sum(isTreetopRadius),
           mergeClusterSize = if_else(identicalElevationsInCluster > 0, n(), 0)) %>%
    ungroup() %>%
    mutate(isTreetopRadius = factor(isTreetopRadius))
  Sys.time() - datasetCmmStart
  #table(treetopDataCmm$treetop) # written to layers: 61848 tops (5842 from merges) + 1746 merge points + 14930 unmatched tops + 381 noise points
  #tibble(radiusTops = sum(treetopDataCmm$isTreetopRadius == TRUE))
}

#tibble(chmSingletonTops = sum(treetopDataDsm$isTreetopRadiusChm & (treetopDataDsm$treetop == "yes")), # 70894 singleton tops
#       chmMergeTops = sum(treetopDataDsm$isTreetopRadiusChm & (treetopDataDsm$treetop == "merge")), # 6084 tops in merge clusters
#       dsmSingletonTops = sum(treetopDataDsm$isTreetopRadiusDsm & (treetopDataDsm$treetop == "yes")), # 70894 singleton tops
#       dsmMergeTops = sum(treetopDataDsm$isTreetopRadiusDsm & (treetopDataDsm$treetop == "merge"))) # 8261 tops in merge clusters
# number of distinct tops 
#treetopDataDsm %>% group_by(tile) %>% 
#  summarize(treetopsChm = sum(if_else(identicalElevationsInCluster > 0, 1 / mergeClusterSize, 0)), 
#            treetopsDsm = sum(if_else(identicalHeightsInCluster > 0, 1 / mergeClusterSize, 0)), 
#            `single top` = sum(treetop == "yes"), `merge point` = sum(treetop == "merge"), `residual noise` = sum(treetop == "noise"), other = sum(treetop == "no")) %>%
#  bind_rows(summarize(., across(where(is.numeric), sum)))

# variable selection
if (treetopOptions$includeInvestigatory)
{
  library(forcats)
  library(VSURF)
  #predictorVariables = names(s4268maximaDsm %>% select(-tile, -id, -mergeClusterID, -uniqueID, -uniqueMergeClusterID, -x, -y, -sourceID, -dsmZ, -cmmZ))
  #colSums(is.na(s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables))))
  
  treetopVsurf = VSURF(treetop ~ ., s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)), 
                       parallel = TRUE, ncores = treetopOptions$rangerThreads, RFimplem = "ranger")
  saveRDS(treetopVsurf, "trees/segmentation/treetops/treetop vsurf s4268 458k all.Rds")
  #treetopVsurf = readRDS("trees/segmentation/treetops/treetop vsurf s4268 458k all")
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
  ggsave("trees/segmentation/treetops/treetop vsurf s4268 458k all.png", width = 14, height = 0.33 * nrow(predictorImportance), units = "cm", dpi = 150)
  
  library(Boruta)
  borutaMaxima = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables))
  (borutaStart = Sys.time())
  treetopBoruta = Boruta(treetop ~ ., borutaMaxima, num.threads = treetopOptions$rangerThreads, doTrace = 2) # 1.5 days @ 88 iterations, 116 predictors @ 458k maxima, 9900X
  (borutaTime = Sys.time() - borutaStart)
  ggplot() +
    geom_violin(aes(x = importance, y = reorder(variable, importance, median)), as_tibble(treetopBoruta$ImpHistory) %>% pivot_longer(cols = everything(), names_to = "variable", values_to = "importance"), draw_quantiles = c(0.25, 0.50, 0.75), width = 6) +
    labs(y = NULL)
  ggsave("trees/segmentation/treetops/treetop boruta s4268 458k all.png", width = 22, height = 0.4 * ncol(treetopBoruta$ImpHistory), units = "cm", dpi = 100)
  saveRDS(treetopBoruta, "trees/segmentation/treetops/treetop boruta s4268 458k all.Rds")
  
  # correlation
  # ggcorplot() too slow to be useful
  s4268correlation = cor(s4268maximaDsm %>% select(-tile, -id, -mergeClusterID, -uniqueID, -uniqueMergeClusterID, -treetop, -x, -y)) # ~5 s, 9900X
  s4268correlationLongform = as_tibble(s4268correlation) %>% mutate(variable1 = rownames(s4268correlation)) %>% pivot_longer(cols = -variable1, names_to = "variable2", values_to = "correlation")
  ggplot() +
    geom_tile(aes(x = variable1, y = variable2, fill = correlation), s4268correlationLongform) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "correlation") +
    scico::scale_fill_scico(palette = "cork", limits = c(-1, 1)) +
    scale_y_discrete(limits = rev) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("trees/segmentation/treetops/treetop correlation s4268 458k all.png", width = 40, height = 40, units = "cm", dpi = 100)

  # PCA
  predictorPca = prcomp(~ ., s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(-tile, -id, -mergeClusterID, -uniqueID, -uniqueMergeClusterID, -sourceID, -treetop, -x, -y) %>% mutate(across(ends_with("differentSourceID"), as.numeric)), scale = TRUE) # ~20 s @ 458k rows, 9900X
  factoextra::fviz_eig(predictorPca, ncp = 20)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(1, 2), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 0102 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(3, 4), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 0304 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(5, 6), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 0506 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(7, 8), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 0708 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(9, 10), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 0910 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(11, 12), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/treetops/treetop PCA 458k axes 1112 v2.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)

  # MCA  
  predictorFamd = FactoMineR::FAMD(s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(-treetop, -tile, -id, -mergeClusterID, -uniqueID, -uniqueMergeClusterID, -sourceID, -x, -y), graph = FALSE, ncp = 6)
  factoextra::fviz_famd_var(predictorFamd, col.ind = "cos2", gradient.cols = c("blue", "orange", "red"), axes = c(5, 6), labelsize = 2, repel = TRUE)
  
  # LDA
  ldaData = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% slice_sample(n = 100000) %>% # ordination plots start to get slow with > 100k points
    select(-tile, -id, -uniqueID, -sourceID, -x, -y, # exclude non-predictors
           -dsmZ, -cmmZ, # exclude non-portable predictors
           -rangeMean, -rangeMeanNormalized, -neighborDistance2mean, -neighborDistance3mean, -neighborDistance4mean, -neighborDistance5mean, -starts_with("neighborDifferentSourceID"), -neighborDistance5variance, -neighborProminence5mean, -prominenceMean, -prominenceMeanNormalized, -rangeNeighbor5Normalized, -netProminence, -netProminenceNormalized, -starts_with("netProminenceNeighbor"), -netRange, -netRangeNormalized) %>% # exclude colinear variables
    mutate(across(where(is.numeric), scale))
  predictorLda = MASS::lda(treetop ~ ., ldaData, tol = 1E-4)
  ggord::ggord(predictorLda, ldaData$treetop, arrow = 0.2, axes = c(1, 2), direction = "both", ext = 0.99, force = 10, grp_title = NULL, max.overlaps = 150, repel = TRUE, size = 0.5, txt = 2.5, vec_ext = 100, veccol = "grey40", xlims = NA * c(-1, 1), ylims = NA * c(-1, 1)) +
  ggord::ggord(predictorLda, ldaData$treetop, arrow = 0.2, axes = c(3, 4), direction = "both", ext = 0.99, force = 10, grp_title = NULL, max.overlaps = 150, repel = TRUE, size = 0.5, txt = 2.5, vec_ext = 100, veccol = "grey40", xlims = NA * c(-1, 1), ylims = NA * c(-1, 1)) +
  patchwork::plot_annotation(theme = theme(plot.margin = margin())) +
  patchwork::plot_layout(guides = "collect")
}
if (treetopOptions$includeSetup)
{
  # rows   predictors       CPU    threads Boruta VSURF  selected  trees  tune 70  mtry  min node  sample fraction
  # 458k   all              9900X  12      36.22h 39.21h 27(81)    500    
  # 458k   VSURF Pde        9900X  12                    24        500    8.823h   9     3         0.6913  
  # 458k   hr               9900X  12                    2         500    50.71m   2     632       0.2199           
  library(tuneRanger)
  library(mlr)
  #predictorVariables = c("treetop", variablesPrediction)
  rangerTuneData = as.data.frame(s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables))) # placing inline to makeClassifTask() fails with "Must have length 1" for certain combinations of predictor variables (but works with most combinations?!)
  rangerTuneTask = makeClassifTask(data = rangerTuneData, target = "treetop")
  #estimateStart = Sys.time()
  #estimateTimeTuneRanger(rangerTuneTask, num.trees = 500, num.threads = 14, iters = 70)# 5950X: 1.5 min to estimate 1h47m @ 164k tops and 12 predictors
  #Sys.time() - estimateStart
  (rangerTuneStart = Sys.time())
  rangerTuning = tuneRanger(rangerTuneTask, measure = list(multiclass.brier), num.trees = 500, num.threads = 12, iters = 70,
                            build.final.model = FALSE)
  (rangerTuneTime = Sys.time() - rangerTuneStart)
  (rangerTuning)
  saveRDS(rangerTuning, "trees/segmentation/treetops/treetop tune s4268 458k VSURF Pde.Rds")
  detach("package:tuneRanger", unload = TRUE)
  detach("package:mlrMBO", unload = TRUE)
  detach("package:mlr", unload = TRUE)
  #rangerTuning = readRDS("trees/segmentation/treetops/treetop tune s4268 458k PMCADAV18.Rds")
}

# random forest fitting
if (treetopOptions$fitRandomForest)
{
  # predictor variable selection
  # controls
  #predictorVariables = c("treetop", "height", "radius")
  #rangerTuning = tibble(mtry = 2, minNodeSize = 632, sampleFraction = 0.2199076)
  #predictorVariables = c("treetop", "height", "radius", "netProminenceNormalized")
  #rangerTuning = tibble(mtry = 2, minNodeSize = 143, sampleFraction = 0.255)
  # VSURF selections for prediction and interpretation
  predictorVariables = c("treetop", "prominence2normalized", "netProminenceNormalized", "ring2variance", "height", "prominence3normalized", "cmmSlope3", "prominence1normalized", "radius", "ring1variance", "prominence4normalized", "ring3variance", "netProminenceNeighbor1normalized", "mean2deltaNormalized", "mean1delta", "prominence5normalized", "mean3deltaNormalized", "slope5normalized", "mean4delta", "rangeNeighbor3Normalized", "neighbor1distance", "neighborDistance20mean", "netProminenceNeighbor3normalized", "neighborDistance50mean", "ring4variance") # VSURF Pde
  rangerTuning = tibble(mtry = 9, minNodeSize = 3, sampleFraction = 0.6912568)
  #predictorVariables = c("treetop", "prominence2normalized", "netProminenceNormalized", "prominenceMeanNormalized", "ring2variance", "height", "prominence3normalized", "cmmSlope3", "prominence1normalized", "radius", "ring1variance", "prominence4normalized", "ring3variance", "netProminenceNeighbor1normalized", "neighbor1prominenceNormalized", "mean2deltaNormalized", "prominence5normalized", "mean3deltaNormalized", "deltaCmm", "slope5normalized", "prominenceVarianceNormalized", "slope4normalized", "slope2normalized", "mean1deltaNormalized", "prominence5", "slope3normalized", "slope1normalized", "netRangeNormalized", "rangeMeanNormalized", "netRange", "dsmSlope", "mean4deltaNormalized", "range2normalized", "range5normalized", "range1normalized", "range4normalized", "range3normalized", "netProminenceNeighbor2normalized", "rangeNeighbor3Normalized", "neighbor1distance", "neighbor2prominenceNormalized", "netProminenceNeighbor4normalized", "rangeNeighbor5Normalized", "neighborDistance20mean", "netProminenceNeighbor3normalized", "netProminenceNeighbor2", "neighborDistance3mean", "neighborDistance5mean", "mean5deltaNormalized", "neighborDistance50mean", "ring4variance", "neighborDistance10mean", "neighborDistance4mean", "netProminenceNeighbor5", "neighborProminence5mean", "netProminenceNeighbor3") # VSURF allInterpDedup
                         
  # variables which need to flow through cross validation for accuracy measurements but should not be used for predictors
  accuracyVariables = c("tile", "id", "uniqueID", "sourceID", "x", "y", "dsmZ", "cmmZ", "uniqueMergeClusterID") # x and y coordinates are also used for merge point kNNs
  if (sum(accuracyVariables %in% predictorVariables) > 0)
  {
    stop("At least one accuracy assessment variable excluded from training data is indicated as a predictor variable. Should it be removed from the exclusion list?")
  }
  
  # cross validated accuracy estimation
  #                                                                                              cross validated accuracy
  # rows   predictors   CPU    threads trees  mtry  min node  sample fraction  cross validation  treetop  overall
  # 458k   hr           9900X  2x12    500    2     632       0.2199           2x25 7.394m       0.838    0.926    worse than f(h); overfit?
  # 458k   VSURF Pde    9900X  12      500    9     3         0.6913           2x25 1.660h       0.945    0.965
  handlers(global = TRUE)
  handlers("cli")
  #plan(multisession, workers = 2) # two workers beneficial with short ranger prediction times, diminishing but still positive returns with three
  
  s4268training = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables), all_of(accuracyVariables)) # drop maxima in tiles neighboring training region
  (crossValidationStartTime = Sys.time())
  s4268crossValidation = fit_ranger_treetop(s4268training, s4268maximaDsm, 
                                            mtry = rangerTuning$mtry, minNodeSize = rangerTuning$minNodeSize, sampleFraction = rangerTuning$sampleFraction, 
                                            folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
  (crossValidationTime = Sys.time() - crossValidationStartTime)
  s4268crossValidation %>% summarize(n = mean(nTrain + nValidate), cvTime = crossValidationTime, treetopAccuracy = mean(treetopAccuracy), overallAccuracy = mean(overallAccuracy), treetopMiscountPct = 100 * mean((treetop - expectedTreetop) / expectedTreetop))
  saveRDS(s4268crossValidation, paste0("trees/segmentation/treetops/random forest s4268 458k VSURF Pde 2x25 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))
  #s4268crossValidation = readRDS(paste0("trees/segmentation/random forest s4268 458k VSURF Pde 2x25 m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))
  #print.noquote(extend_confusion_matrix_to_string(s4268crossValidation$confusionSubmatrix[[1]]$table))
  #s4268crossValidation %>% select(noise, expectedNoise)
  
  ggplot() +
    geom_violin(aes(x = overallAccuracy, y = "overall", color = "overall"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    geom_violin(aes(x = treetopAccuracy, y = "treetop", color = "treetop"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    geom_violin(aes(x = nonTreetopAccuracy, y = "not treetop", color = "not treetop"), s4268crossValidation, draw_quantiles = c(0.25, 0.5, 0.75), width = 0.6) +
    labs(x = "accuracy", y = NULL) + 
    scale_y_discrete(limits = rev(c("overall", "treetop", "not treetop"))) +
  ggplot() +
    geom_violin(aes(x = overallAccuracy, y = heightClass, color = "overall", group = heightClass), unnest(s4268crossValidation %>% select(repetition, fold, overallAccuracyByHeight), cols = overallAccuracyByHeight) %>% filter(n > 25)) +
    labs(x = "accuracy", y = "tree height, m") +
    scale_y_continuous(breaks = seq(0, 90, by = 10)) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    coord_cartesian(xlim = c(0.5, 1)) &
    guides(color = "none")

  range((unnest(s4268crossValidation %>% select(repetition, fold, overallAccuracyByHeight), cols = overallAccuracyByHeight))$n)
  
  # fit random forest
  rangerStartTime = Sys.time()
  treetopRandomForest = ranger::ranger(treetop ~ ., data = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                       mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                       sample.fraction = rangerTuning$sampleFraction,
                                       num.threads = treetopOptions$rangerThreads)
  treetopRandomForest
  Sys.time() - rangerStartTime
  saveRDS(treetopRandomForest, paste0("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".Rds"))
  #treetopRandomForest = readRDS("trees/segmentation/treetops/treetopRandomForest s4268 458 legacy14 m8n25.Rds")

  # variable importance
  randomForestImportance = ranger::ranger(treetop ~ ., data = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                          importance = "permutation",
                                          mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                          sample.fraction = rangerTuning$sampleFraction,
                                          num.threads = treetopOptions$rangerThreads)
  randomForestLocalImportance = ranger::ranger(treetop ~ ., data = s4268maximaDsm %>% filter(is.na(treetop) == FALSE) %>% select(all_of(predictorVariables)),
                                               importance = "permutation", local.importance = TRUE,
                                               mtry = rangerTuning$mtry, splitrule = 'gini', min.node.size = rangerTuning$minNodeSize, 
                                               sample.fraction = rangerTuning$sampleFraction,
                                               num.threads = treetopOptions$rangerThreads)
  
  predictorLabels = tibble(predictor = c("prominence2normalized", "prominence3normalized", "netProminenceNormalized", "radius", "height", "ring2variance", "prominence1normalized", "cmmSlope3", "prominence4normalized", "ring1variance", "ring3variance", "mean1delta",
                                         "prominence5normalized", "mean4delta", "mean2deltaNormalized", "mean3deltaNormalized", "slope5normalized", "netProminenceNeighbor1normalized", "ring4variance", "netProminenceNeighbor3normalized", "neighbor1distance", "neighborDistance20mean", "neighborDistance50mean", "rangeNeighbor3Normalized"),
                           label = c("ring 2 prominence, normalized", "ring 3 prominence, normalized", "net prominence, normalized", "dominance radiance", "height", "ring 2 variance", "ring 1 prominence, normalized", "canopy maxima model slope", "ring 4 prominence, normalized", "ring 1 variance", "ring 3 variance", "ring 1 mean, relative",
                                     "ring 5 prominence, normalized", "ring 4 mean, relative", "ring 2 mean, relative normalized", "ring 3 mean, relative normalized", "ring 5 slope", "neighbor 1 prominence, normalized", "ring 4 variance", "neighbor 3 prominence, normalized", "neighbor 1 distance", "neighborhood density, nearest 20 maxima", "neighborhood density, nearest 50 maxima", "neighbor 1, 2, 3 height range, normalized"))
  globalImportance = left_join(tibble::as_tibble_row(randomForestImportance$variable.importance) %>% # pivot_longer() since as_tibble_col() doesn't facilitate row names
                                  pivot_longer(everything(), names_to = "predictor", values_to = "importance"),
                               predictorLabels,
                               by = join_by(predictor))  %>%
    mutate(importance = 100 * importance / max(importance)) %>%
    arrange(desc(importance)) %>%
    mutate(predictor = factor(predictor, levels = predictor),
           label = factor(label, levels = label))
  localImportance = as_tibble(randomForestLocalImportance$variable.importance.local) %>% mutate(treetop = (s4268maximaDsm %>% filter(is.na(treetop) == FALSE))$treetop) %>% group_by(treetop) %>%
    summarize(across(everything(), mean)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~100 * .x / max(across(where(is.numeric))))) %>%
    pivot_longer(!treetop, names_to = "predictor", values_to = "importance") %>%
    mutate(predictor = factor(predictor, levels = levels(globalImportance$predictor)))

  saveRDS(globalImportance, paste0("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, " global importance.Rds"))  
  saveRDS(localImportance, paste0("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, " local importance.Rds"))
  
  ggplot() + # also Figure TBD in results.R
    geom_raster(aes(x = "global", y = label, fill = importance), globalImportance) +
    labs(x = NULL, y = NULL, title = plotLetters[1]) +
    scale_y_discrete(limits = rev) +
  ggplot() +
    geom_raster(aes(x = treetop, y = predictor, fill = importance), localImportance) +
    labs(x = NULL, y = NULL, title = plotLetters[2]) +
    scale_x_discrete(labels = c("single top", "merge point", "noise", "processing\nartifact", "other"), limits = c("yes", "merge", "noise", "maybe noise", "no")) +
    scale_y_discrete(labels = NULL, limits = rev) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 1, ncol = 2, guides = "collect") &
    coord_equal() &
    labs(x = NULL, y = NULL, fill = "relative\npermutation\nimportance, %") &
    scale_fill_viridis_c(option = "plasma", limits = c(0, 100 + 1E-14)) &
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 10))
  ggsave(file.path(getwd(), paste0("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m", rangerTuning$mtry, "n", rangerTuning$minNodeSize, ".png")), units = "cm", height = 16, width = 14, dpi = 200)
  
  print(tibble(variable = names(randomForestImportance$variable.importance), importance = randomForestImportance$variable.importance, relativePct = 100 * importance / max(importance)) %>% arrange(desc(importance)), n = 25)
}

if (treetopOptions$includeInvestigatory)
{
  # variable selection
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
  #predictorVariables = c("treetop", "height", "radius", "netProminenceNormalized", "deltaCmm", "cmmSlope3", "prominenceMean", "netProminence", "mean3delta", "rangeMeanNormalized", "netRange", "netRangeNormalized", "rangeVarianceNormalized", "prominenceNeighbor5Variance", "neighborDistance5mean", "netProminenceNeighbor5", "ring2variance", "prominence1normalized", "prominence2normalized", "prominence3normalized", "prominence4normalized", "mean2deltaNormalized", "rangeMean", "neighbor1prominence", "mean4delta") # PCA+MCA+LDA+VSURF
  #rangerTuning = tibble(mtry = 13, minNodeSize = 4, sampleFraction = 0.629)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "prominenceMean", "netProminence", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "netRange", "mean2deltaNormalized", "netRangeNormalized", "rangeMeanNormalized", "deltaCmm") # PMCDAV18
  #rangerTuning = tibble(mtry = 11, minNodeSize = 7, sampleFraction = 0.577)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2deltaNormalized", "netRangeNormalized", "rangeMeanNormalized", "deltaCmm") # PMCADAV15
  #rangerTuning = tibble(mtry = 8, minNodeSize = 4, sampleFraction = 0.565)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2deltaNormalized", "netRangeNormalized", "deltaCmm") # PMCADAV14
  #rangerTuning = tibble(mtry = 8, minNodeSize = 7, sampleFraction = 0.591)
  #predictorVariables = c("treetop", "height", "radius", "prominence1", "prominence2", "prominence3", "prominence4", "netProminence", "ring2variance", "cmmSlope3", "neighbor1prominence", "rangeMean", "mean2delta", "netRange", "deltaCmm") # PMCADAV14un
  #rangerTuning = tibble(mtry = 10, minNodeSize = 5, sampleFraction = 0.533)
  #predictorVariables = c("treetop", "radius", "prominence2normalized", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "height", "ring2variance", "prominence4normalized", "cmmSlope3", "neighbor1prominence", "mean2deltaNormalized", "netRangeNormalized", "deltaCmm") # PMCADAV13
  #rangerTuning = tibble(mtry = 7, minNodeSize = 11, sampleFraction = 0.549)
  #rangerTuning = tibble(mtry = rangerTuning$recommended.pars$mtry, minNodeSize = rangerTuning$recommended.pars$min.node.size, sampleFraction = rangerTuning$recommended.pars$sample.fraction)
  
  # leave one tile out cross validation
  s04200w06810data = s4268maximaDsm %>% filter(tile == "s04200w06810")
  s04200w06840data = s4268maximaDsm %>% filter(tile == "s04200w06840")
  s04230w06810data = s4268maximaDsm %>% filter(tile == "s04230w06810")
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
  #caret::confusionMatrix(randomForestFit)
  randomForestConfusionMatrix = caret::confusionMatrix(predict(treetopRandomForest, s4268maximaDsm %>% select(-treetop)), s4268maximaDsm$treetop)
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
  #writeVector(vect(tileMaxima, crs = tileCrs, geom = c("x", "y")), file.path(candidateTreetopsDsmPath, "rf", paste0(tileName, " ranger.gpkg")), layer = "treetops", insert = TRUE, overwrite = TRUE)
  
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
  tileTreetops = bind_rows(tileMaxima %>% filter(tileMaxima$treetop == "yes") %>% mutate(mergePoints = as.integer(1), mergePointOnTile = as.integer(1)),
                           tileMergePoints$treetops)
  
  treetopsFilePath = file.path(candidateTreetopsDsmPath, "rf", paste0(tileName, ".gpkg"))
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


## joins of local maxima with treetops and merge points added manually in QGIS
if (treetopOptions$includeSetup)
{
  maxAutomaticallyGeneratedMergeClusterSize = 9
  rebuildMergeClusters = FALSE
  recalculateMergeClusterSize = FALSE
  
  rowMins = function(x)
  {
    if (is.null(nrow(x))) # nrow returns NULL if x is a vector; in this case assume a single row sliced from a matrix
    {
      return(min(x))
    }
    
    return(apply(x, 1, min))
  }
  
  tileName = "s04230w06810" # s04200w06810, s04200w06840, s04230w06810
  acceptedTreetopsFilePath = file.path(acceptedTreetopsDsmPath, paste0(tileName, ".gpkg"))
  
  localMaximaTile = st_read(file.path(localMaximaPathV3beta, paste0(tileName, ".gpkg")), layer = "localMaxima", quiet = TRUE) # 2.5D so sf is needed
  localMaximaTileXY = st_coordinates(localMaximaTile)[, c("X", "Y")]

  # annotate merge points with local maxima IDs
  mergePointTile = st_read(acceptedTreetopsFilePath, layer = "merge treetops", quiet = TRUE)
  mergePointKnn = get.knnx(localMaximaTileXY, st_coordinates(mergePointTile)[, c("X", "Y")], k = 2)
  tibble(naID = sum(is.na(mergePointTile$id)), deltaID = sum(mergePointTile$id != localMaximaTile$id[mergePointKnn$nn.index[, 1]], na.rm = TRUE))
  mergePointTile$id = localMaximaTile$id[mergePointKnn$nn.index[, 1]]
  
  # check for misplaced merge points
  tibble(id = mergePointTile$id[mergePointKnn$nn.index[, 1]], distance = mergePointKnn$nn.dist[, 1]) %>% filter(distance > 1E-6)
  tibble(naID = sum(is.na(mergePointTile$id)))
  # check for duplicated merge points
  tibble(index = mergePointKnn$nn.index[, 2], distance = mergePointKnn$nn.dist[, 2]) %>% filter(distance < 1E-6)
  #mergePointTile$nearestNeighborDistance = mergePointKnn$nn.dist[, 1]
  
  # sync accepted treetops with local maxima and merge points
  acceptedTreetopTile = st_read(acceptedTreetopsFilePath, layer = "treetops", quiet = TRUE) # 2.5D so sf is needed
  acceptedTreetopCoordinates = st_coordinates(acceptedTreetopTile)[, c("X", "Y")]
  acceptedTreetopKnn = get.knnx(localMaximaTileXY, acceptedTreetopCoordinates, k = 1)
  acceptedTreetopTile$nearestNeighborDistance = acceptedTreetopKnn$nn.dist[, 1]
  
  # check for tops too near to each other (within two DSM cells) and for duplicated tops
  duplicatedTopsKnn = get.knnx(acceptedTreetopCoordinates, acceptedTreetopCoordinates, k = 2) # k = 2 since nearest neighbor is self
  duplicatedTopsKnn = tibble(index = duplicatedTopsKnn$nn.index[, 2], distance = duplicatedTopsKnn$nn.dist[, 2]) 
  duplicatedTopsKnn %>% filter(distance < (2 * 1.5 - 0.1)) %>% mutate(treeID = acceptedTreetopTile$treeID[index]) %>% relocate(treeID)
  
  if (rebuildMergeClusters)
  {
    if (recalculateMergeClusterSize)
    {
      acceptedMergeTopIndices = which(acceptedTreetopKnn$nn.dist[, 1] > 1E-6)
      #acceptedMergeKnn = get.knnx(localMaximaTileXY, st_coordinates(acceptedTreetopTile)[acceptedMergeTopIndices, c("X", "Y")], k = maxMergeClusterSize)
      #acceptedMergeKnn$neighborIDs = matrix(localMaximaTile$id[acceptedMergeKnn$nn.index], ncol = maxMergeClusterSize)
      acceptedMergeKnn = get.knnx(st_coordinates(mergePointTile)[, c("X", "Y")], st_coordinates(acceptedTreetopTile)[acceptedMergeTopIndices, c("X", "Y")], k = maxAutomaticallyGeneratedMergeClusterSize)
      acceptedMergeKnn$neighborIDs = matrix(mergePointTile$id[acceptedMergeKnn$nn.index], ncol = maxAutomaticallyGeneratedMergeClusterSize)
      
      normalizedNeighbor3kDistance = matrix(acceptedMergeKnn$nn.dist[, 3:maxAutomaticallyGeneratedMergeClusterSize] / rowMeans(acceptedMergeKnn$nn.dist[, 1:2]), ncol = maxAutomaticallyGeneratedMergeClusterSize - 2)
      mergeClusterSize = as.integer(2 + rowSums(normalizedNeighbor3kDistance < 1.4))
      
      acceptedTreetopTile$mergePoints = as.integer(1)
      acceptedTreetopTile$mergePoints[acceptedMergeTopIndices] = mergeClusterSize
      acceptedTreetopTile$mergePointsOnTile = acceptedTreetopTile$mergePoints
    } else {
      # check for invalid merge cluster sizes
      tibble(naMerge = sum(is.na(acceptedTreetopTile$mergePoints)), outOfRangeMerge = sum((acceptedTreetopTile$mergePoints < 1) | (acceptedTreetopTile$mergePoints > 20), na.rm = TRUE)) # sanity upper bound
      
      acceptedMergeTopIndices = which(acceptedTreetopTile$mergePoints > 1)
      mergeClusterID = acceptedTreetopTile$treeID[acceptedMergeTopIndices]
      mergeClusterSize = acceptedTreetopTile$mergePoints[acceptedMergeTopIndices]
      
      acceptedMergeKnn = get.knnx(st_coordinates(mergePointTile)[, c("X", "Y")], st_coordinates(acceptedTreetopTile)[acceptedMergeTopIndices, c("X", "Y")], k = max(mergeClusterSize))
      acceptedMergeKnn$neighborIDs = matrix(mergePointTile$id[acceptedMergeKnn$nn.index], ncol = max(mergeClusterSize))
      colnames(acceptedMergeKnn$nn.dist) = str_c("neighbor", seq(1, max(mergeClusterSize)))
      colnames(acceptedMergeKnn$neighborIDs) = str_c("neighbor", seq(1, max(mergeClusterSize)))

      sum(is.na(acceptedTreetopTile$mergePoints))
      table(mergeClusterSize)
      
      if (is.null(acceptedTreetopTile$mergePointsOnTile))
      {
        acceptedTreetopTile$mergePointsOnTile = acceptedTreetopTile$mergePoints
      }
    }
    
    mergeClusterSizeOnTile = acceptedTreetopTile$mergePointsOnTile[acceptedMergeTopIndices]
    mergeClusterIndices01 = which(mergeClusterSizeOnTile == 1)
    mergeClusterIndices02 = which(mergeClusterSizeOnTile == 2)
    mergeClusterIndices03 = which(mergeClusterSizeOnTile == 3)
    mergeClusterIndices04 = which(mergeClusterSizeOnTile == 4)
    mergeClusterIndices05 = which(mergeClusterSizeOnTile == 5)
    mergeClusterIndices06 = which(mergeClusterSizeOnTile == 6)
    mergeClusterIndices07 = which(mergeClusterSizeOnTile == 7)
    mergeClusterIndices08 = which(mergeClusterSizeOnTile == 8)
    mergeClusterIndices09 = which(mergeClusterSizeOnTile == 9)
    mergeClusterIndices10 = which(mergeClusterSizeOnTile == 10)
    minOnTileIDinMergeCluster = rep(NA, length(acceptedMergeTopIndices))
    if (length(mergeClusterIndices01) > 0)
    {
      minOnTileIDinMergeCluster[mergeClusterIndices01] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices01, 1:2])
    }
    minOnTileIDinMergeCluster[mergeClusterIndices02] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices02, 1:2])
    minOnTileIDinMergeCluster[mergeClusterIndices03] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices03, 1:3])
    minOnTileIDinMergeCluster[mergeClusterIndices04] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices04, 1:4])
    minOnTileIDinMergeCluster[mergeClusterIndices05] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices05, 1:5])
    minOnTileIDinMergeCluster[mergeClusterIndices06] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices06, 1:6])
    minOnTileIDinMergeCluster[mergeClusterIndices07] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices07, 1:7])
    if (length(mergeClusterIndices08) > 0)
    {
      minOnTileIDinMergeCluster[mergeClusterIndices08] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices08, 1:8])
    }
    if (length(mergeClusterIndices09) > 0)
    {
      minOnTileIDinMergeCluster[mergeClusterIndices09] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices09, 1:9])
    }
    if (length(mergeClusterIndices10) > 0)
    {
      minOnTileIDinMergeCluster[mergeClusterIndices10] = rowMins(acceptedMergeKnn$neighborIDs[mergeClusterIndices10, 1:10])
    }
    sum(is.na(minOnTileIDinMergeCluster))
    
    # check for suspiciously formed clusters
    # Current bug: cluster kNN looks for number of neighbors, not number of on tile neighbors
    # tile           edge cluster IDs           large cluster IDs
    # s04200w06840   99, 69568, 96369, 150712   10555, 13794, 18799, 106840, 112556, 120368, 124614, 132381, 140816, 145537, 149595, 151140, 153510, 153589, 153761        
    if (tileName == "s04200w06840")
    {
      excludedClusterIDs = c(10555, 13794, 18799, 106840, 112556, 120368, 124614, 132381, 140816, 145537, 149595, 151140, 153510, 153589, 153761) # s04200w06840 large clusters (edge clusters c(99, 69568, 96369, 150712))
      offTileIDs = c(154505, 154595)
    } else if (tileName == "s04200w06810") {
      excludedClusterIDs = c(410, 4948, 4977, 5947, 10938, 11918, 13015, 13087, 13166, 13924, 13925, 14013, 14016, 14806, 15065, 15243, 16613, 18175, 23575, 31851, 32207, 35310, 40527, 41011, 41471, 41830, 42755, 45546, 53199, 54532, 54982, 69388, 71085, 72892, 100220, 110327, 124651, 135388, 140352, 140544, 150869, 153971, 154010, 156117, 156752, 160291) # large clusters
      offTileIDs = c(68, 134341, 134484)
    } else if (tileName == "s04230w06810") {
      excludedClusterIDs = c(10502, 15322, 16041, 23269, 24432, 26093, 32584, 34481, 35229, 35534, 36055, 39790, 40318, 42755, 46360, 47106, 52243, 57901, 76449, 79431, 81699, 91841, 103176, 114516, 117815, 120271, 142139) # large clusters
      offTileIDs = c(60, 63, 38159, 43537)
    }

    mergeDistances = as_tibble(acceptedMergeKnn$nn.dist) %>% mutate(clusterID = minOnTileIDinMergeCluster, clusterSize = mergeClusterSizeOnTile) %>%
      pivot_longer(starts_with("neighbor"), names_prefix = "neighbor", names_to = "neighbor", values_to = "distance") %>%
      filter(neighbor <= clusterSize) %>%
      group_by(clusterID) %>%
      mutate(meanDistance = mean(distance))
    
    mergeClusters = as_tibble(acceptedMergeKnn$neighborIDs) %>% mutate(clusterID = minOnTileIDinMergeCluster, clusterSize = mergeClusterSizeOnTile) %>%
      pivot_longer(starts_with("neighbor"), names_prefix = "neighbor", names_to = "neighbor", values_to = "id") %>%
      filter(neighbor <= clusterSize)
    
    # overly large clusters likely indicate missing merge points or accepted tops with merge point counts set higher than the number of points in the cluster
    print(mergeDistances %>% filter(meanDistance > 4.3, (clusterID %in% excludedClusterIDs) == FALSE) %>% summarize(meanDistance = meanDistance[1], clusterSize = clusterSize[1]) %>% arrange(desc(meanDistance)), n = 80)
    #suspiciousID = 91841        
    #cbind(treeID = acceptedTreetopTile$treeID[acceptedMergeTopIndices], mergeClusterSize, minOnTileIDinMergeCluster, acceptedMergeKnn$neighborIDs)[which(minOnTileIDinMergeCluster == suspiciousID), ]
    
    # a local maxima should not appear in more than one merge cluster, maxima which do likely indicate missing merge points, duplicate merge points, or incorrect cluster counts
    mergeClusters %>% group_by(id) %>% summarize(numberOfClustersContainingLocalMaxima = n()) %>% filter(numberOfClustersContainingLocalMaxima > 1)
    #suspiciousID = 49575
    #mergeClusters %>% filter(id == suspiciousID) # find matching local maxima

    # set merge treetops' treeIDs to the lowest on tile value in each merge cluster and update their heights
    # Singleton tops are handled separately below.
    # BUGBUG: heights of off tile local maxima are not included. 
    heightByMergeClusterID = left_join(tibble(clusterID = minOnTileIDinMergeCluster), # alignment frame so heights match ordering of acceptedMergeTopIndices
                                       left_join(mergeClusters, localMaximaTile, by = join_by(id)) %>% group_by(clusterID) %>% summarize(height = mean(height)), # mean height by cluster
                                       by = join_by(clusterID))
    # sum(heightByMergeClusterID$clusterID != minOnTileIDinMergeCluster) # zero if heights are correctly aligned
    
    #acceptedTreetopTile$treeIDsync = acceptedTreetopTile$treeID
    #acceptedTreetopTile$treeIDsync[acceptedMergeTopIndices] = minOnTileIDinMergeCluster
    #table((acceptedTreetopTile %>% filter(treeID != treeIDsync))$mergePoints)
    #print(as_tibble(acceptedTreetopTile[acceptedMergeTopIndices, ] %>% mutate(neighborIDs = acceptedMergeKnn$neighborIDs) %>% filter(treeID != treeIDsync)), n = 100)
    tibble(newIDs = sum(is.na(acceptedTreetopTile$treeID[acceptedMergeTopIndices])), clusterIDchanges = sum(acceptedTreetopTile$treeID[acceptedMergeTopIndices] != minOnTileIDinMergeCluster, na.rm = TRUE), heightChanges = sum(acceptedTreetopTile$height[acceptedMergeTopIndices] != heightByMergeClusterID$height, na.rm = TRUE), heightAssigns = sum(is.na(acceptedTreetopTile$height[acceptedMergeTopIndices])), radiiClearance = sum(is.na(acceptedTreetopTile$radius[acceptedMergeTopIndices]) == FALSE))
    
    acceptedTreetopTile$treeID[acceptedMergeTopIndices] = minOnTileIDinMergeCluster
    acceptedTreetopTile$height[acceptedMergeTopIndices] = heightByMergeClusterID$height
    acceptedTreetopTile$radius[acceptedMergeTopIndices] = NA_real_ # clear any dominance radii as they are not well defined for merge clusters
    
    # set merge points' merge cluster IDs
    mergeClusters %>% group_by(id) %>% summarize(numberOfClustersContainingLocalMaxima = n()) %>% filter(numberOfClustersContainingLocalMaxima > 1) # if not empty if following left_join() adds rows!
    mergePointTileWithClusterIDs = left_join(mergePointTile %>% select(-clusterID), mergeClusters %>% select(-clusterSize, -neighbor), by = join_by(id))
    #mergeClusters %>% filter(id == mergePointTile[351, ]$id) # row 351 of x matches multiple rows in y - duplicate inclusion from edge tree 69568
    #mergePointTile %>% filter(id == mergeClusters[3396, ]$id) # row 1483 of y matches multiple rows in x - duplicate merge points, only reports first one
    if (nrow(mergePointTile) != nrow(mergePointTileWithClusterIDs))
    {
      stop("Some merge points appear in more than one cluster.")
      #mergeClusters %>% filter(id == 73573)
    }
    # NA clusterIDs likely indicate unneeded merge points for singleton tops or merge tops with too low merge point counts (and thus possibly also with missing merge points)
    # But NA IDs also indicate merge points belonging to an off tile merge top.
    # tile           merge points with off tile top
    # s04200w06840   154505, 154595
    tibble(idMismatches = sum(mergePointTile$id != mergePointTileWithClusterIDs$id), missingClusters = sum(is.na(mergePointTileWithClusterIDs$clusterID)) - length(offTileIDs))
    as_tibble(mergePointTileWithClusterIDs) %>% filter(is.na(clusterID)) %>% filter((id %in% offTileIDs) == FALSE)

    mergePointTile$clusterID = if_else(is.na(mergePointTileWithClusterIDs$clusterID), mergePointTile$clusterID, mergePointTileWithClusterIDs$clusterID) # don't revert off tile cluster IDs back to NA
    tibble(naIDs = sum(is.na(mergePointTile$clusterID)), offTileIDs = length(offTileIDs))
    mergePointTile %>% filter(is.na(clusterID))
  }

  acceptedSingletonTopIndices = which(acceptedTreetopTile$mergePoints == 1)
  
  # basic check for possible missing merge clusters
  # Anticipated to return potentially upwards of a hundred results for (re)review.
  missingMergeKnn = get.knnx(localMaximaTileXY, acceptedTreetopCoordinates[acceptedSingletonTopIndices, ], k = 2) # k = 2 since nearest neighbor is self
  missingMergeKnn = tibble(treeID = acceptedTreetopTile$treeID[acceptedSingletonTopIndices], index = missingMergeKnn$nn.index[, 2], height = acceptedTreetopTile$height[acceptedSingletonTopIndices], distance = missingMergeKnn$nn.dist[, 2], dsmZ = tibble::num(localMaximaTile$dsmZ[acceptedTreetopKnn$nn.index[acceptedSingletonTopIndices]], digits = 1), dsmZnearestNeighbor = tibble::num(localMaximaTile$dsmZ[index], digits = 1), deltaZabs = abs(dsmZ - dsmZnearestNeighbor))
  print(missingMergeKnn %>% filter(deltaZabs >= 0, deltaZabs <= 0.1, height > 12, distance <= 3) %>% arrange(treeID), n = 100)

  # set singleton treetops' treeIDs and sync their height and radius
  tibble(newIDs = sum(is.na(acceptedTreetopTile$treeID[acceptedSingletonTopIndices])), idChanges = sum(acceptedTreetopTile$treeID[acceptedSingletonTopIndices] != localMaximaTile$id[acceptedTreetopKnn$nn.index][acceptedSingletonTopIndices], na.rm = TRUE))
  acceptedTreetopTile$treeID[acceptedSingletonTopIndices] = localMaximaTile$id[acceptedTreetopKnn$nn.index][acceptedSingletonTopIndices] # or update only NA treeIDs if present?
  acceptedTreetopTile$height[acceptedSingletonTopIndices] = localMaximaTile$height[acceptedTreetopKnn$nn.index][acceptedSingletonTopIndices]
  acceptedTreetopTile$radius[acceptedSingletonTopIndices] = localMaximaTile$radius[acceptedTreetopKnn$nn.index][acceptedSingletonTopIndices]

  #st_write(acceptedTreetopTile, acceptedTreetopsFilePath, layer = "treetops", delete_dsn = FALSE, delete_layer = TRUE)
  #st_write(mergePointTile, acceptedTreetopsFilePath, layer = "merge treetops", delete_dsn = FALSE, delete_layer = TRUE)
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
  elliott2021tiles = st_read("GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg", layer = "Elliott tile index", quiet = TRUE)
  elliott2021tiles %<>% filter(Tile_ID == tileName)
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


## local maxima numbering stability between DSM versions
if (treetopOptions$includeInvestigatory)
{
  tile = "s04200w06840"
  localMaxima = st_read(file.path(localMaximaPathV3beta, paste0(tile, ".gpkg")), layer = "localMaxima", quiet = TRUE)
  localMaxima2 = st_read(file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3/local maxima", paste0(tile, ".gpkg")), layer = "localMaximaDsm", quiet = TRUE)
  
  tibble(currentMaxima = nrow(localMaxima), nextMaxima = nrow(localMaxima2))
  localMaximaKnn = get.knnx(st_coordinates(localMaxima2)[, c("X", "Y")], st_coordinates(localMaxima)[, c("X", "Y")], k = 1)
  localMaximaKnn = tibble(id = localMaxima$id, index = localMaximaKnn$nn.index[, 1], distance = localMaximaKnn$nn.dist[, 1])
  localMaximaKnn %>% filter(distance > 0)

  dsm = rast(file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta", paste0(tile, ".tif")))
  dsm2 = rast(file.path(dsmPath, paste0(tile, ".tif")))
  
  table(dsm2$dsm[,] - dsm$dsm[,])
  tibble(heightChange = sum(dsm2$dsm[,] != dsm$dsm[,], na.rm = TRUE), na = sum(is.na(dsm$dsm[,])), na2 = sum(is.na(dsm2$dsm[,])))
}


## dataset transfer from DSM to CHM and CMM
if (treetopOptions$includeSetup)
{
  bind_merge_points = function(...)
  {
    return(bind_rows(...) %>% group_by(tile, clusterID) %>% mutate(mergeClusterSize = n()) %>% ungroup())
  }
  
  # yields both singleton tops where enough merge points have been lost as well as tops from surviving merge clusters
  get_merge_tops = function(localMaxima, mergePoints, verify = TRUE)
  {
    #localMaxima = localMaxima46chm
    #mergePoints = chmMergePoints
    #mergePoints = chmMergePointsUnmatched
    
    # find tops from cluster centroids
    mergePointsXY = st_coordinates(mergePoints)
    mergeTops = as_tibble(st_drop_geometry(mergePoints)) %>% mutate(x = mergePointsXY[, "X"], y = mergePointsXY[, "Y"]) %>%
      group_by(tile, clusterID) %>% 
      summarize(across(.cols = any_of(c("id", "originatingClusterID")), .fns = min), height = mean(height), mergePoints = n(), x = mean(x), y = mean(y), radius = NA, .groups = "drop") %>%
      rename(treeID = id, any_of(c(originatingTreeID = "originatingClusterID")))
    if (verify & any(mergeTops$treeID != mergeTops$clusterID))
    {
      stop("Merge top generation failure. Some tree IDs differ from their merge cluster IDs, indicating dataset inconsistency.")
    } else {
      mergeTops %<>% select(-clusterID)
    }
    
    # set radius on single tops
    singleTopsIndices = which(mergeTops$mergePoints == 1)
    if (length(singleTopsIndices) > 0)
    {
      singleMergePointsXY = mergePointsXY[singleTopsIndices, c("X", "Y")]
      if (length(singleTopsIndices) == 1)
      {
        # someMatrix[one row index, ] collapses to a vector, so have to rematrix it
        singleMergePointsXY = matrix(singleMergePointsXY, ncol = 2, dimnames = list(c(), c("X", "Y")))
      }
      
      localMaximaXY = st_coordinates(localMaxima)[, c("X", "Y")]
      singleTopsKnn = get.knnx(localMaximaXY, singleMergePointsXY, k = 1)
      if (verify & any(singleTopsKnn$nn.dist > 1E-6)) # could use st_join() instead of get.knnx() but get.knnx() allows distance checking
      {
        stop("Merge point locations do not exactly match local maxima.")
      }
      
      matchedSingleTopsIndices = which(singleTopsKnn$nn.dist[, 1] < 1E-6)
      mergeTops$radius[singleTopsIndices[matchedSingleTopsIndices]] = localMaxima$radius[singleTopsKnn$nn.index[matchedSingleTopsIndices, 1]]
    }
    
    return(st_as_sf(mergeTops, coords = c("x", "y"), crs = crs(mergePoints), sf_column_name = "geom"))
  }
  
  # a given local maxima can be designated as 
  #  1) both a treetop and a merge point
  #  2) both a treetop and a noise point
  #  3) both a merge point and a noise point
  #  4) conceivably all three (though this did not occur on s04200w06810, s04200w06840, or s04230w06810)
  # so local maxima assignment is not exclusive among categories
  maxima_match_and_remove = function(localMaxima, singleTops, mergePoints, noisePoints, maybeNoisePoints, maxDistance = 1E-6) # or height if less than maxDistance
  {
    #localMaxima = localMaxima46chm
    #singleTops = st_transform(st_zm(acceptedTreetops46dsm), crs(mergePoints46dsm)) %>% filter(mergePoints == 1)
    #mergePoints = mergePoints46dsm
    #noisePoints = noise46dsm
    #maybeNoisePoints = maybeNoise46dsm
    #maxDistance = 1E-6
    #localMaxima = chmMatches0$localMaximaUnmatched
    #singleTops = chmMatches0$singleTopsUnmatched
    #mergePoints = chmMatches0$mergePointsUnmatched
    #noisePoints = chmMatches0$noisePointsUnmatched
    #maybeNoisePoints = chmMatches0$maybeNoiseUnmatched
    #maxDistance = 1 * sqrt(2) * treetopOptions$dsmCellSize + 0.01
    
    # assumes tops and merge points have been validated not to overlap
    if (any(singleTops$mergePoints != 1))
    {
      stop("Some rows of singleTops have mergePoints > 1, indicating the top is not actually single.")
    }
    
    # transfer single tops as best able given conflicts
    # Potentially simpler to slice singleTopsMatched from localMaxima but this isn't currently done as most localMaxima fields don't
    # need to flow.
    localMaximaXyh = bind_cols(st_coordinates(localMaxima)[, c("X", "Y")], H = localMaxima$height)
    singleTopsXyh = bind_cols(st_coordinates(singleTops)[, c("X", "Y")], H = singleTops$height)
    singleTopsKnn = get.knnx(localMaximaXyh, singleTopsXyh, k = 1)
    singleTopsMatchIndices = which(singleTopsKnn$nn.dist[, 1] <= pmin(singleTops$height, maxDistance))
    singleTopsUnmatchedIndices = setdiff(1:nrow(singleTops), singleTopsMatchIndices)
    
    singleTopMatching = tibble(topIndex = singleTopsMatchIndices,
                               maximaIndex = singleTopsKnn$nn.index[singleTopsMatchIndices, 1],
                               distance = singleTopsKnn$nn.dist[singleTopsMatchIndices, 1]) %>%
      group_by(maximaIndex) %>%
      slice_min(distance, n = 1, with_ties = FALSE) # when multiple tops match the same maxima, take only the closest top and leave the rest to potentially pair to another maxima in a subsequent iteration
      
    singleTopsMatched = st_set_geometry(singleTops[singleTopMatching$topIndex, ], st_geometry(st_as_sf(localMaximaXyh[singleTopMatching$maximaIndex, c("X", "Y")], coords = c("X", "Y"), crs = crs(singleTops)))) # st_set_geometry() drops CRS if not set on st_as_sf()
    singleTopsMatched$tile = localMaxima$tile[singleTopMatching$maximaIndex] # rarely changes but simplest to always update
    singleTopsMatched$originatingTreeID = singleTopsMatched$treeID # retain originating ID (DSM local maxima ID in current use cases)
    singleTopsMatched$treeID = localMaxima$id[singleTopMatching$maximaIndex] # update ID to local maxima ID in current surface (CHM or CMM)
    singleTopsMatched$height = localMaxima$height[singleTopMatching$maximaIndex] # update height and radius to new surface
    singleTopsMatched$radius = localMaxima$radius[singleTopMatching$maximaIndex]
    # flow through remaining fields: notes, mergePoints (1, by definition), mergePointsOnTile (also 1)
    singleTopsMatched$nearestNeighborDistance = NULL # TBD what to do here, for now field is retained by cleared to block flow of data that's no longer valid
  
    # transfer merge points as best able given conflicts
    mergePointsXyh = bind_cols(st_coordinates(mergePoints)[, c("X", "Y")], H = mergePoints$height)
    mergePointKnn = get.knnx(localMaximaXyh, mergePointsXyh)
    mergePointMatchIndices = which(mergePointKnn$nn.dist[, 1] <= pmin(mergePoints$height, maxDistance))
    mergePointsUnmatchedIndices = setdiff(1:nrow(mergePoints), mergePointMatchIndices)
    
    mergePointMatching = tibble(pointIndices = mergePointMatchIndices,
                                maximaIndex = mergePointKnn$nn.index[mergePointMatchIndices, 1],
                                distance = mergePointKnn$nn.dist[mergePointMatchIndices, 1]) %>%
      group_by(maximaIndex) %>%
      slice_min(distance, n = 1, with_ties = FALSE)
    
    mergePointsMatched = st_set_geometry(mergePoints[mergePointMatching$pointIndices, ], st_geometry(st_as_sf(localMaximaXyh[mergePointMatching$maximaIndex, c("X", "Y")], coords = c("X", "Y"), crs = crs(mergePoints))))
    mergePointsMatched$tile = localMaxima$tile[mergePointMatching$maximaIndex]
    mergePointsMatched$originatingID = mergePointsMatched$id
    mergePointsMatched$id = localMaxima$id[mergePointMatching$maximaIndex]
    mergePointsMatched$height = localMaxima$height[mergePointMatching$maximaIndex]
    # mergePointsMatched$clusterID is needed later for merge top regeneration and cannot be updated or cleared here
    mergePointsMatched$nearestNeighborDistance = NULL

    # transfer noise points as best able given conflicts
    # local maxima can be also be designated both treetops and noise points
    noisePointsXyh = bind_cols(st_coordinates(noisePoints)[, c("X", "Y")], H = noisePoints$height)
    noisePointKnn = get.knnx(localMaximaXyh, noisePointsXyh)
    noisePointMatchIndices = which(noisePointKnn$nn.dist[, 1] < pmin(noisePoints$height, maxDistance))
    noisePointUnmatchedIndices = setdiff(1:nrow(noisePoints), noisePointMatchIndices)
    
    noisePointMatching = tibble(pointIndices = noisePointMatchIndices,
                                maximaIndex = noisePointKnn$nn.index[noisePointMatchIndices, 1],
                                distance = noisePointKnn$nn.dist[noisePointMatchIndices, 1]) %>%
      group_by(maximaIndex) %>%
      slice_min(distance, n = 1, with_ties = FALSE)
    
    noisePointsMatched = st_set_geometry(noisePoints[noisePointMatching$pointIndices, ], st_geometry(st_as_sf(localMaximaXyh[noisePointMatching$maximaIndex, c("X", "Y")], coords = c("X", "Y"), crs = crs(noisePoints))))
    noisePointsMatched$originatingID = noisePointsMatched$id
    noisePointsMatched$id = localMaxima$id[noisePointMatching$maximaIndex]
    noisePointsMatched$height = localMaxima$height[noisePointMatching$maximaIndex]
    # flow notes through
    
    # transfer maybe noise points
    maybeNoisePointsXyh = bind_cols(st_coordinates(maybeNoisePoints)[, c("X", "Y")], H = maybeNoisePoints$height)
    maybeNoisePointKnn = get.knnx(localMaximaXyh, maybeNoisePointsXyh)
    maybeNoisePointMatchIndices = which(maybeNoisePointKnn$nn.dist[, 1] < pmin(maybeNoisePoints$height, maxDistance))
    maybeNoisePointUnmatchedIndices = setdiff(1:nrow(maybeNoisePoints), maybeNoisePointMatchIndices)
    
    maybeNoisePointMatching = tibble(pointIndices = maybeNoisePointMatchIndices,
                                     maximaIndex = maybeNoisePointKnn$nn.index[maybeNoisePointMatchIndices, 1],
                                     distance = maybeNoisePointKnn$nn.dist[maybeNoisePointMatchIndices, 1]) %>%
      group_by(maximaIndex) %>%
      slice_min(distance, n = 1, with_ties = FALSE)
    
    if (nrow(maybeNoisePointMatching) > 0)
    {
      maybeNoisePointsMatched = st_set_geometry(maybeNoisePoints[maybeNoisePointMatching$pointIndices, ], st_geometry(st_as_sf(localMaximaXyh[maybeNoisePointMatching$maximaIndex, c("X", "Y")], coords = c("X", "Y"), crs = crs(maybeNoisePoints))))
      maybeNoisePointsMatched$originatingID = maybeNoisePointsMatched$id
      maybeNoisePointsMatched$id = localMaxima$id[maybeNoisePointMatching$maximaIndex]
      maybeNoisePointsMatched$height = localMaxima$height[maybeNoisePointMatching$maximaIndex]
    # flow notes through
    } else {
      maybeNoisePointsMatched = maybeNoisePoints[c(), ]
    }

    localMaximaMatchIndices = c(singleTopsKnn$nn.index[singleTopsMatchIndices, 1], mergePointKnn$nn.index[mergePointMatchIndices, 1])
    localMaximaUnmatchedIndices = setdiff(1:nrow(localMaxima), localMaximaMatchIndices)
    return(list(localMaximaMatched = localMaxima[localMaximaMatchIndices, ],
                singleTopsMatched = singleTopsMatched,
                mergePointsMatched = mergePointsMatched,
                noisePointsMatched = noisePointsMatched,
                maybeNoiseMatched = maybeNoisePointsMatched,
                localMaximaUnmatched = localMaxima[localMaximaUnmatchedIndices, ],
                singleTopsUnmatched = singleTops[singleTopsUnmatchedIndices, ],
                mergePointsUnmatched = mergePoints[mergePointsUnmatchedIndices, ],
                noisePointsUnmatched = noisePoints[noisePointUnmatchedIndices, ],
                maybeNoiseUnmatched = maybeNoisePoints[maybeNoisePointUnmatchedIndices, ],
                stats = tibble(maxima = length(localMaximaMatchIndices), tops = nrow(singleTopsMatched), merge = nrow(mergePointsMatched), noise = nrow(noisePointsMatched), maybeNoise = nrow(maybeNoisePointsMatched),
                               maximaRemaining = length(localMaximaUnmatchedIndices), topsRemaining = length(singleTopsUnmatchedIndices), mergeRemaining = length(mergePointsUnmatchedIndices), noiseRemaining = length(noisePointUnmatchedIndices), maybeNoiseRemaining = length(maybeNoisePointUnmatchedIndices),
                               minDistance = min(c(singleTopsKnn$nn.dist[, 1], mergePointKnn$nn.dist[, 1])), maxDistance = maxDistance)))
  }
  
  # discard merge cluster remnants if part of a merge cluster was matched into a treetop
  remove_remnant_merge_points = function(mergePointsUnmatched, mergePointsMatched)
  {
    # would be one line with dplyr if cur_group() could refer to groups in mergePointsUnmatched instead of treetops
    # But it apparently it doesn't do that (errors out on NULL as mergePointsUnmatched is ungrouped).
    # nonRemnantPoints = mergePointsUnmatched %>% group_by(tile) %>% 
    #  filter((clusterID %in% (mergePointsMatched %>% filter(tile == cur_group()$tile))$originatingClusterID) == FALSE)
    filterFrame = left_join(mergePointsUnmatched %>% nest(.by = tile, .key = "mergePointsUnmatched"),
                            mergePointsMatched %>% nest(.by = tile, .key = "mergePointsMatched"),
                            by = join_by(tile)) %>%
      mutate(nonRemnantMergePoints = vector(mode = "list", length = n()))
    for (row in 1:nrow(filterFrame))
    {
      filterFrame$nonRemnantMergePoints[[row]] = filterFrame$mergePointsUnmatched[[row]] %>% filter((clusterID %in% filterFrame$mergePointsMatched[[row]]$originatingClusterID) == FALSE) %>%
        mutate(tile = filterFrame$tile[row]) # unnest() is not built for this, so restore tile and then bind_rows()
    }
    
    return(bind_rows(filterFrame$nonRemnantMergePoints))
  }
  
  verify_disjoint = function(treetops, mergePoints, unmatchedTreetops, noisePoints, maybeNoisePoints)
  {
    #treetops = chmTreetops
    #mergePoints = chmMergePoints
    #unmatchedTreetops = chmUnmatchedTreetops
    #noisePoints = chmNoisePoints
    #maybeNoisePoints = chmMaybeNoise
    
    treetopsXY = st_coordinates(treetops)[, c("X", "Y")]
    mergePointsXY = st_coordinates(mergePoints)[, c("X", "Y")]
    unmatchedTreetopsXY = st_coordinates(unmatchedTreetops)[, c("X", "Y")]
    noisePointsXY = st_coordinates(noisePoints)[, c("X", "Y")]
    maybeNoisePointsXY = st_coordinates(maybeNoisePoints)[, c("X", "Y")]
    
    allXY = rbind(treetopsXY, mergePointsXY, unmatchedTreetopsXY, noisePointsXY, maybeNoisePointsXY)
    allKnn = get.knn(allXY, k = 1)
    nondisjointIndices = which(allKnn$nn.dist < 0.4 * treetopOptions$dsmCellSize) # formal minimum is 0.5 cell size but allow tolerance for inexact merge top placement

    nondisjoint = tibble(allIndex = nondisjointIndices, # includes both members of a nondisjoint pair or all pairings in a cluster
                         source = if_else(allIndex <= nrow(treetops), "treetop", 
                                          if_else(allIndex <= (nrow(treetops) + nrow(mergePoints)), "merge point",
                                                  if_else(allIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops)), "unmatched top", 
                                                          if_else(allIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops) + nrow(noisePoints)), "noise point", "maybe noise")))),
                         sourceIndex = case_match(source, "treetop" ~ if_else(allIndex <= nrow(treetops), allIndex, NA_integer_),
                                                          "merge point" ~ if_else((allIndex > nrow(treetops)) & (allIndex <= (nrow(treetops) + nrow(mergePoints))), allIndex - nrow(treetops), NA_integer_),
                                                          "unmatched top" ~ if_else((allIndex > (nrow(treetops) + nrow(mergePoints))) & (allIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops))), allIndex - nrow(treetops) - nrow(mergePoints), NA_integer_),
                                                          "noise point" ~ if_else(allIndex > (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops)), allIndex - nrow(treetops) - nrow(mergePoints) - nrow(unmatchedTreetops), NA_integer_),
                                                          "maybe noise" ~ if_else(allIndex > (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops) + nrow(noisePoints)), allIndex - nrow(treetops) - nrow(mergePoints) - nrow(unmatchedTreetops) - nrow(noisePoints), NA_integer_)),
                         tile = case_match(source, "treetop" ~ treetops$tile[sourceIndex], "merge point" ~ mergePoints$tile[sourceIndex], "unmatched top" ~ unmatchedTreetops$tile[sourceIndex], "noise point" ~ noisePoints$tile[sourceIndex], "maybe noise" ~ maybeNoisePoints$tile[sourceIndex]),
                         id = case_match(source, "treetop" ~ treetops$treeID[sourceIndex], "merge point" ~ mergePoints$id[sourceIndex], "unmatched top" ~ unmatchedTreetops$originatingTreeID[sourceIndex], "noise point" ~ noisePoints$id[sourceIndex], "maybe noise" ~ maybeNoisePoints$id[sourceIndex]),
                         neighborAllIndex = allKnn$nn.index[nondisjointIndices, ], 
                         neighborSource = if_else(neighborAllIndex <= nrow(treetops), "treetop", 
                                              if_else(neighborAllIndex <= (nrow(treetops) + nrow(mergePoints)), "merge point",
                                                      if_else(neighborAllIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops)), "unmatched top", 
                                                              if_else(neighborAllIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops) + nrow(noisePoints)), "noise point", "maybe noise")))),
                         neighborIndex = case_match(neighborSource, "treetop" ~ if_else(neighborAllIndex <= nrow(treetops), neighborAllIndex, NA_integer_),
                                                                    "merge point" ~ if_else((neighborAllIndex > nrow(treetops)) & (neighborAllIndex <= (nrow(treetops) + nrow(mergePoints))), neighborAllIndex - nrow(treetops), NA_integer_),
                                                                    "unmatched top" ~ if_else((neighborAllIndex > (nrow(treetops) + nrow(mergePoints))) & (neighborAllIndex <= (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops))), neighborAllIndex - nrow(treetops) - nrow(mergePoints), NA_integer_),
                                                                    "noise point" ~ if_else(neighborAllIndex > (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops)), neighborAllIndex - nrow(treetops) - nrow(mergePoints) - nrow(unmatchedTreetops), NA_integer_),
                                                                    "maybe noise" ~ if_else(neighborAllIndex > (nrow(treetops) + nrow(mergePoints) + nrow(unmatchedTreetops) + nrow(noisePoints)), neighborAllIndex - nrow(treetops) - nrow(mergePoints) - nrow(unmatchedTreetops) - nrow(noisePoints), NA_integer_)),
                         neighborTile = case_match(source, "treetop" ~ treetops$tile[neighborIndex], "merge point" ~ mergePoints$tile[neighborIndex], "unmatched top" ~ unmatchedTreetops$tile[neighborIndex], "noise point" ~ noisePoints$tile[neighborIndex], "maybe noise" ~ maybeNoisePoints$tile[neighborIndex]),
                         neighborID = case_match(source, "treetop" ~ treetops$treeID[neighborIndex], "merge point" ~ mergePoints$id[neighborIndex], "unmatched top" ~ unmatchedTreetops$originatingTreeID[neighborIndex], "noise point" ~ noisePoints$id[neighborIndex], "maybe noise" ~ maybeNoisePoints$id[neighborIndex]),
                         height = case_match(source, "treetop" ~ treetops$height[sourceIndex], "merge point" ~ mergePoints$height[sourceIndex], "unmatched top" ~ unmatchedTreetops$height[sourceIndex], "noise point" ~ noisePoints$height[sourceIndex], "maybe noise" ~ maybeNoisePoints$height[neighborIndex]),
                         distance = allKnn$nn.dist[allIndex]) %>%
      relocate(tile, source, id, neighborTile, neighborSource, neighborID, height, distance) %>%
      filter((tile != neighborTile) | (id != neighborID)) # exclude duplicate rows from second, third, fourth... members of pairs and clusters
    
    conflictingNondisjoint = nondisjoint %>% filter(source == neighborSource)
    if (nrow(conflictingNondisjoint) > 0)
    {
      # stub handler for now based only on (near-)duplicate designations in the same class
      # Difficult to determine if treetops in close proximity to noise points or merge points are problematic as
      # treetops can directly overlap merge points and noise points can directly overlap both. Also, merge clusters shapes which don't
      # place their treetop directly over a merge point can still place a top quite close to a merge point.
      # Could also re-check for misplaced noise points here but it's assume that's done in dataset checking rather than in
      # dataset transfer.
      print(conflictingNondisjoint)
      stop("Conflicting non-disjoint points found. See above.")
    }
  }
  
  write_treetops = function(directoryPath, tileName, treetops, mergePoints, unmatchedTreetops, noisePoints, maybeNoisePoints, tileSuffix = " chm")
  {
    tilePath = file.path(directoryPath, paste0(tileName, tileSuffix, ".gpkg"))
    st_write(treetops %>% filter(tile == tileName), tilePath, layer = "treetops", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(mergePoints %>% filter(tile == tileName), tilePath, layer = "merge treetops", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(unmatchedTreetops %>% filter(tile == tileName), tilePath, layer = "unmatched treetops", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(noisePoints %>% filter(tile == tileName), tilePath, layer = "noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(maybeNoisePoints %>% filter(tile == tileName), tilePath, layer = "maybe noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
  }
  
  localMaxima46dsm = bind_rows(st_read(file.path(localMaximaPathV3beta, "s04200w06840.gpkg"), layer = "localMaxima", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(localMaximaPathV3beta, "s04200w06810.gpkg"), layer = "localMaxima", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(localMaximaPathV3beta, "s04230w06810.gpkg"), layer = "localMaxima", quiet = TRUE) %>% mutate(tile = "s04230w06810"))
  
  s04200w06810treetops = st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "treetops", quiet = TRUE) %>% mutate(tile = "s04200w06810")
  
  acceptedTreetops46dsm = bind_rows(st_zm(s04200w06810treetops, drop = FALSE, what = "Z"), # add z dimension since st breaks on layers with mixed xy and xyz geometries but leave as zero as z values are not used here
                                    st_transform(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "treetops", quiet = TRUE) %>% mutate(tile = "s04200w06840"), crs(s04200w06810treetops)),
                                    st_transform(st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "treetops", quiet = TRUE) %>% mutate(tile = "s04230w06810"), crs(s04200w06810treetops)))
  
  mergePoints46dsm = bind_rows(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "merge treetops", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "merge treetops", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "merge treetops", quiet = TRUE) %>% mutate(tile = "s04230w06810")) %>%
    mutate(height = left_join(., as_tibble(localMaxima46dsm), by = join_by(tile, id))$height)
  
  noise46dsm = bind_rows(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "noise points", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                         st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "noise points", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                         st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "noise points", quiet = TRUE) %>% mutate(tile = "s04230w06810")) %>%
    st_join(., st_transform(localMaxima46dsm %>% select(id, height), crs = crs(.)), join = st_nearest_feature) # use st_nearest_feature as the default of st_intersects has no numerical tolerance (st_snap() might help but is too slow to use

  maybeNoise46dsm = bind_rows(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "maybe noise points", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                              st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "maybe noise points", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                              st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "maybe noise points", quiet = TRUE) %>% mutate(tile = "s04230w06810")) %>%
    st_join(., st_transform(localMaxima46dsm %>% select(id, height), crs = crs(.)), join = st_nearest_feature)
  
  # transfer to CHM
  # 77049 treetops + 13600 merge points + 2435 noise points -> 74560 tops + 7338 merge points + 2439 unmatched tops + 2299 noise points
  # 3.2% treetop loss
  localMaxima46chm = bind_rows(st_read(file.path(localMaximaPathV3, "s04200w06840.gpkg"), layer = "localMaximaChm", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(localMaximaPathV3, "s04200w06810.gpkg"), layer = "localMaximaChm", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(localMaximaPathV3, "s04230w06810.gpkg"), layer = "localMaximaChm", quiet = TRUE) %>% mutate(tile = "s04230w06810"))
  
  chmMatches0 = maxima_match_and_remove(localMaxima46chm, st_transform(st_zm(acceptedTreetops46dsm), crs(mergePoints46dsm)) %>% filter(mergePoints == 1), mergePoints46dsm, noise46dsm, maybeNoise46dsm) # ~6s, due to QGIS's current lack of xyz support xyz treetops are flattened to xy points for merge point consistency
  chmMatches1 = maxima_match_and_remove(chmMatches0$localMaximaUnmatched, chmMatches0$singleTopsUnmatched, chmMatches0$mergePointsUnmatched, chmMatches0$noisePointsUnmatched, chmMatches0$maybeNoiseUnmatched, maxDistance = sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  chmMatches2 = maxima_match_and_remove(chmMatches1$localMaximaUnmatched, chmMatches1$singleTopsUnmatched, chmMatches1$mergePointsUnmatched, chmMatches1$noisePointsUnmatched, chmMatches1$maybeNoiseUnmatched, maxDistance = 2 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  chmMatches3 = maxima_match_and_remove(chmMatches2$localMaximaUnmatched, chmMatches2$singleTopsUnmatched, chmMatches2$mergePointsUnmatched, chmMatches2$noisePointsUnmatched, chmMatches2$maybeNoiseUnmatched, maxDistance = 3 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  chmMatches4 = maxima_match_and_remove(chmMatches3$localMaximaUnmatched, chmMatches3$singleTopsUnmatched, chmMatches3$mergePointsUnmatched, chmMatches3$noisePointsUnmatched, chmMatches3$maybeNoiseUnmatched, maxDistance = 4 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)

  bind_rows(tibble(ring = NA, maxima = 0, tops = 0, merge = 0, noise = 0, maybeNoise = 0, maximaRemaining = nrow(localMaxima46chm), topsRemaining = sum(acceptedTreetops46dsm$mergePoints == 1), mergeRemaining = nrow(mergePoints46dsm), noiseRemaining = nrow(noise46dsm), maybeNoiseRemaining = nrow(maybeNoise46dsm), minDistance = NA, maxDistance = NA),
            chmMatches0$stats %>% mutate(ring = 0),
            chmMatches1$stats %>% mutate(ring = 1),
            chmMatches2$stats %>% mutate(ring = 2),
            chmMatches3$stats %>% mutate(ring = 3),
            chmMatches4$stats %>% mutate(ring = 4)) %>%
    mutate(maximaTotal = maxima + maximaRemaining, topsTotal = tops + topsRemaining, mergeTotal = merge + mergeRemaining)
  
  chmMergePoints = bind_merge_points(chmMatches0$mergePointsMatched, chmMatches1$mergePointsMatched, chmMatches2$mergePointsMatched, chmMatches3$mergePointsMatched) %>%
    mutate(originatingClusterID = clusterID) %>%
    group_by(tile, clusterID) %>% 
    mutate(clusterID = min(id), mergeClusterSize = n()) %>% # update cluster sizes and IDs based on transference
    ungroup()
  chmMergeTops = get_merge_tops(localMaxima46chm, chmMergePoints)
  chmMergePointsUnmatched = remove_remnant_merge_points(chmMatches3$mergePointsUnmatched, chmMergePoints)
  chmMergePoints %<>% filter(mergeClusterSize > 1) # remove merge points which yielded single tops

  chmTreetops = bind_rows(chmMatches0$singleTopsMatched, chmMatches1$singleTopsMatched, chmMatches2$singleTopsMatched, chmMatches3$singleTopsMatched, chmMergeTops)
  chmUnmatchedTreetops = bind_rows(chmMatches3$singleTopsUnmatched, get_merge_tops(localMaxima46chm, chmMergePointsUnmatched, verify = FALSE)) %>% 
    rename(originatingTreeID = treeID)
  
  chmNoisePoints = bind_rows(chmMatches0$noisePointsMatched, chmMatches1$noisePointsMatched, chmMatches2$noisePointsMatched)
  chmMaybeNoise = bind_rows(chmMatches0$maybeNoiseMatched, chmMatches1$maybeNoisesMatched, chmMatches2$maybeNoiseMatched)
  
  verify_disjoint(chmTreetops, chmMergePoints, chmUnmatchedTreetops, chmNoisePoints, chmMaybeNoise)
  
  tibble(treetops = nrow(chmTreetops), mergePoints = nrow(chmMergePoints), unmatchedTops = nrow(chmUnmatchedTreetops), noisePoints = nrow(chmNoisePoints), maybeNoise = nrow(chmMaybeNoise), treetopLossPct = 100 * (1 - treetops / nrow(acceptedTreetops46dsm)), treetopsFromMerges = nrow(chmMergeTops))
  
  treetopsAcceptedChm = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted chm"
  write_treetops(treetopsAcceptedChm, tileName = "s04200w06810", chmTreetops, chmMergePoints, chmUnmatchedTreetops, chmNoisePoints, chmMaybeNoise)
  write_treetops(treetopsAcceptedChm, tileName = "s04200w06840", chmTreetops, chmMergePoints, chmUnmatchedTreetops, chmNoisePoints, chmMaybeNoise)
  write_treetops(treetopsAcceptedChm, tileName = "s04230w06810", chmTreetops, chmMergePoints, chmUnmatchedTreetops, chmNoisePoints, chmMaybeNoise)
  
  # transfer to CMM
  # 77049 treetops + 13600 merge points + 2435 noise points-> 61848 tops + 1746 merge points + 14930 unmatched tops + 381 noise points
  # 19.7% treetop loss
  localMaxima46cmm = bind_rows(st_read(file.path(localMaximaPathV3, "s04200w06840.gpkg"), layer = "localMaximaCmm", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(localMaximaPathV3, "s04200w06810.gpkg"), layer = "localMaximaCmm", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(localMaximaPathV3, "s04230w06810.gpkg"), layer = "localMaximaCmm", quiet = TRUE) %>% mutate(tile = "s04230w06810"))
  
  cmmMatches0 = maxima_match_and_remove(localMaxima46cmm, st_transform(st_zm(acceptedTreetops46dsm), crs(mergePoints46dsm)) %>% filter(mergePoints == 1), mergePoints46dsm, noise46dsm, maybeNoise46dsm) # a few seconds, due to QGIS's currently lack of xyz support xyz treetops are flattened to xy points for merge point consistency
  cmmMatches1 = maxima_match_and_remove(cmmMatches0$localMaximaUnmatched, cmmMatches0$singleTopsUnmatched, cmmMatches0$mergePointsUnmatched, cmmMatches0$noisePointsUnmatched, cmmMatches0$maybeNoiseUnmatched, maxDistance = sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  cmmMatches2 = maxima_match_and_remove(cmmMatches1$localMaximaUnmatched, cmmMatches1$singleTopsUnmatched, cmmMatches1$mergePointsUnmatched, cmmMatches1$noisePointsUnmatched, cmmMatches1$maybeNoiseUnmatched, maxDistance = 2 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  cmmMatches3 = maxima_match_and_remove(cmmMatches2$localMaximaUnmatched, cmmMatches2$singleTopsUnmatched, cmmMatches2$mergePointsUnmatched, cmmMatches2$noisePointsUnmatched, cmmMatches2$maybeNoiseUnmatched, maxDistance = 3 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  cmmMatches4 = maxima_match_and_remove(cmmMatches3$localMaximaUnmatched, cmmMatches3$singleTopsUnmatched, cmmMatches3$mergePointsUnmatched, cmmMatches3$noisePointsUnmatched, cmmMatches3$maybeNoiseUnmatched, maxDistance = 4 * sqrt(2) * treetopOptions$dsmCellSize + 0.01)
  
  bind_rows(tibble(ring = NA, maxima = 0, tops = 0, merge = 0, noise = 0, maybeNoise = 0, maximaRemaining = nrow(localMaxima46cmm), topsRemaining = sum(acceptedTreetops46dsm$mergePoints == 1), mergeRemaining = nrow(mergePoints46dsm), noiseRemaining = nrow(noise46dsm), maybeNoiseRemaining = nrow(maybeNoise46dsm), minDistance = NA, maxDistance = NA),
            cmmMatches0$stats %>% mutate(ring = 0),
            cmmMatches1$stats %>% mutate(ring = 1),
            cmmMatches2$stats %>% mutate(ring = 2),
            cmmMatches3$stats %>% mutate(ring = 3),
            cmmMatches4$stats %>% mutate(ring = 4)) %>%
    mutate(maximaTotal = maxima + maximaRemaining, topsTotal = tops + topsRemaining, mergeTotal = merge + mergeRemaining)
  
  cmmMergePoints = bind_merge_points(cmmMatches0$mergePointsMatched, cmmMatches1$mergePointsMatched, cmmMatches2$mergePointsMatched, cmmMatches3$mergePointsMatched) %>%
    mutate(originatingClusterID = clusterID) %>%
    group_by(tile, clusterID) %>% 
    mutate(clusterID = min(id), mergeClusterSize = n()) %>% 
    ungroup()
  cmmMergeTops = get_merge_tops(localMaxima46cmm, cmmMergePoints)
  cmmMergePointsUnmatched = remove_remnant_merge_points(cmmMatches3$mergePointsUnmatched, cmmMergePoints)
  cmmMergePoints %<>% filter(mergeClusterSize > 1)
  
  cmmTreetops = bind_rows(cmmMatches0$singleTopsMatched, cmmMatches1$singleTopsMatched, cmmMatches2$singleTopsMatched, cmmMatches3$singleTopsMatched, cmmMergeTops)
  cmmUnmatchedTreetops = bind_rows(cmmMatches3$singleTopsUnmatched, get_merge_tops(localMaxima46cmm, cmmMergePointsUnmatched, verify = FALSE)) %>% 
    rename(originatingTreeID = treeID)
  
  cmmNoisePoints = bind_rows(cmmMatches0$noisePointsMatched, cmmMatches1$noisePointsMatched, cmmMatches2$noisePointsMatched)
  cmmMaybeNoise = bind_rows(cmmMatches0$maybeNoiseMatched, cmmMatches1$maybeNoiseMatched, cmmMatches2$maybeNoiseMatched)
  
  verify_disjoint(cmmTreetops, cmmMergePoints, cmmUnmatchedTreetops, cmmNoisePoints, cmmMaybeNoise)
  
  tibble(treetops = nrow(cmmTreetops), mergePoints = nrow(cmmMergePoints), unmatchedTops = nrow(cmmUnmatchedTreetops), noisePoints = nrow(cmmNoisePoints), maybeNoise = nrow(cmmMaybeNoise), treetopLossPct = 100 * (1 - treetops / nrow(acceptedTreetops46dsm)), treetopsFromMerges = nrow(cmmMergeTops))
  
  treetopsAcceptedCmm = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted cmm"
  write_treetops(treetopsAcceptedCmm, tile = "s04200w06810", tileSuffix = " cmm", cmmTreetops, cmmMergePoints, cmmUnmatchedTreetops, cmmNoisePoints, cmmMaybeNoise)
  write_treetops(treetopsAcceptedCmm, tile = "s04200w06840", tileSuffix = " cmm", cmmTreetops, cmmMergePoints, cmmUnmatchedTreetops, cmmNoisePoints, cmmMaybeNoise)
  write_treetops(treetopsAcceptedCmm, tile = "s04230w06810", tileSuffix = " cmm", cmmTreetops, cmmMergePoints, cmmUnmatchedTreetops, cmmNoisePoints, cmmMaybeNoise)
  
  # ggplot() +
  #   geom_histogram(aes(x = 0.3048 * distance, fill = heightClass, group = heightClass), chmSingleTopKnn2 %>% mutate(heightClass = factor(10 * round(0.1 * 0.3048 * height), levels = rev(seq(0, 90, by = 10)))), binwidth = 0.1) +
  #   labs(x = "displacement, m", y = "single tops", fill = "height, m", title = paste(plotLetters[1], "DSM to CHM second iteration")) +
  # plot_annotation(theme = theme(plot.margin = margin())) +
  # plot_layout(guides = "collect") &
  #   coord_cartesian(xlim = c(0, 10), ylim = c(0, NA)) &
  #   scale_fill_manual(breaks = rev(seq(0, 90, by = 10)), values = rev(cetcolor::cet_pal(n = 10, name = "l11"))) & # linear_gow_65-90_c35_n256, https://github.com/coatless-rpkg/cetcolor/blob/main/R/cet_color_maps.R
  #   scale_x_continuous(breaks = seq(0, 10, by = round(2 * 0.3048 * 1.5, 2)))
}
