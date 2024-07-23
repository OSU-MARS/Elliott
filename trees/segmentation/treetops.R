# some investigatory code here assumes segmentation/setup.R
library(caret)
library(dplyr)
library(FNN)
library(ggplot2)
library(magrittr)
library(patchwork)
library(sf)
library(stringr)
library(terra)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             panel.border = element_blank(), 
                             plot.title = element_text(size = 10)))

treetopOptions = tibble(includeInvestigatory = FALSE)

localMaximaPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta/local maxima"
acceptedTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted"
candidateTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops"

get_treetop_eligible_maxima = function(tileName)
{
  localMaxima = vect(file.path(localMaximaPath, paste0(tileName, ".gpkg")))
  # for now, exclude treetop candidates with so few adjacent points ring 1 or ring 2 has no data since it's 1) it's difficult to tell if these are maxima, 2) it's unlikely they're actually maxima, and 3) these likely comprise < 0.01% of maximas
  localMaxima = subset(localMaxima, (is.na(localMaxima$ring1mean) == FALSE) & (is.na(localMaxima$ring2mean) == FALSE))
  return(localMaxima)
}

get_treetop_predictors = function(localMaxima, acceptedTileName = NULL)
{
  localMaximaXY = geom(localMaxima)[, c("x", "y")]
  
  neighborKnn = get.knn(localMaximaXY, k = 2)
  neighborKnn = tibble(index1 = neighborKnn$nn.index[, 1], distance1 = neighborKnn$nn.dist[, 1],
                       index2 = neighborKnn$nn.index[, 2], distance2 = neighborKnn$nn.dist[, 2])
  
  treetopPredictors = as_tibble(localMaxima) %>% 
    mutate(neighbor1sourceID = sourceID[neighborKnn$index1], neighbor1distance = neighborKnn$distance1, neighbor1dsmZ = dsmZ[neighborKnn$index1],
           neighbor2sourceID = sourceID[neighborKnn$index2], neighbor2distance = neighborKnn$distance2, neighbor2dsmZ = dsmZ[neighborKnn$index2],
           across(where(is.double), ~0.3048 * .x), # convert from feet to m
           deltaCmm = dsmZ - cmmZ,
           mean1delta = dsmZ - ring1mean, mean2delta = dsmZ - ring2mean, mean3delta = dsmZ - ring2mean, mean4delta = dsmZ - ring4mean, mean5delta = dsmZ - ring5mean,
           mean1deltaNormalized = mean1delta / height, mean2deltaNormalized = mean2delta / height, mean3deltaNormalized = mean3delta / height, mean4deltaNormalized = mean4delta / height, mean5deltaNormalized = mean5delta / height,
           neighbor1differentSourceID = sourceID != neighbor1sourceID, neighbor1prominence = dsmZ - neighbor1dsmZ, neighbor1prominenceNormalized = neighbor1prominence / height,
           neighbor2differentSourceID = sourceID != neighbor2sourceID, neighbor2prominence = dsmZ - neighbor2dsmZ, neighbor2prominenceNormalized = neighbor2prominence / height,
           prominence1 = dsmZ - ring1max, prominence2 = dsmZ - ring2max, prominence3 = dsmZ - ring3max, prominence4 = dsmZ - ring4max, prominence5 = dsmZ - ring5max,
           prominence1normalized = prominence1 / height, prominence2normalized = prominence2 / height, prominence3normalized = prominence3 / height, prominence4normalized = prominence4 / height, prominence5normalized = prominence5 / height,
           prominenceStdDev = sd(c_across(prominence1:prominence5)),
           prominenceStdDevNormalized = prominenceStdDev / height,
           range1 = ring1max - ring1min, range2 = ring2max - ring2min, range3 = ring3max - ring3min, range4 = ring4max - ring4min, range5 = ring5max - ring5min,
           range1normalized = range1 / height, range2normalized = range2 / height, range3normalized = range3 / height, range4normalized = range4 / height, range5normalized = range5 / height,
           rangeStdDev = sd(c_across(range1:range5)),
           rangeStdDevNormalized = rangeStdDev / height,
           netProminence = prominence1 + prominence2 + prominence3 + prominence4 + prominence5,
           netProminenceNormalized = netProminence / height,
           netRange = range1 + range2 + range3 + range4 + range5,
           netRangeNormalized = netRange / height,
           varianceNormalized1 = ring1variance / height, varianceNormalized2 = ring2variance / height, varianceNormalized3 = ring3variance / height, varianceNormalized4 = ring4variance / height, varianceNormalized5 = ring5variance / height)

  if (is.null(acceptedTileName))
  {
    return(treetopPredictors)
  }
  
  acceptedTreetopsFilePath = file.path(acceptedTreetopsPath, paste0(acceptedTileName, ".gpkg"))
  if (file.exists(acceptedTreetopsFilePath) == FALSE)
  {
    stop(paste0("Tile ", acceptedTileName, " lacks an accepted treetop GeoPackage at '", acceptedTreetopsFilePath, "'."))
  }
  
  treetopPredictors = treetopPredictors %>% mutate(treetop = factor("no", levels = c("no", "yes", "merge", "noise", "maybe noise")))
  
  availableTruthLayers = vector_layers(acceptedTreetopsFilePath)
  if ("treetops" %in% availableTruthLayers)
  {
    acceptedTreetops = vect(acceptedTreetopsFilePath, layer = "treetops")
    isTreetopKnn = get.knnx(localMaximaXY, geom(acceptedTreetops)[, c("x", "y")], k = 1)
    isTreetopKnn = tibble(index = isTreetopKnn$nn.index[, 1], distance = isTreetopKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    treetopPredictors$treetop[isTreetopKnn$index] = "yes"
  }
  if ("merge treetops" %in% availableTruthLayers)
  {
    mergePoints = vect(acceptedTreetopsFilePath, layer = "merge treetops")
    isMergeKnn = get.knnx(localMaximaXY, geom(mergePoints)[, c("x", "y")], k = 1)
    isMergeKnn = tibble(index = isMergeKnn$nn.index[, 1], distance = isMergeKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    treetopPredictors$treetop[isMergeKnn$index] = "merge"
  }
  if ("noise points" %in% availableTruthLayers)
  {
    noisePoints = vect(acceptedTreetopsFilePath, layer = "noise points")
    noisePointKnn = get.knnx(localMaximaXY, geom(noisePoints)[, c("x", "y")], k = 1)
    noisePointKnn = tibble(index = noisePointKnn$nn.index[, 1], distance = noisePointKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    treetopPredictors$treetop[noisePointKnn$index] = "noise"
  }
  if ("maybe noise points" %in% availableTruthLayers)
  {
    maybeNoisePoints = vect(acceptedTreetopsFilePath, layer = "maybe noise points")
    maybeNoisePointKnn = get.knnx(localMaximaXY, geom(maybeNoisePoints)[, c("x", "y")], k = 1)
    maybeNoisePointKnn = tibble(index = maybeNoisePointKnn$nn.index[, 1], distance = maybeNoisePointKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
    treetopPredictors$treetop[maybeNoisePointKnn$index] = "maybe noise"
  }

  return(treetopPredictors)
}
  
## treetop random forests
# tops   predictors   tune grid  cross validation   train time  accuracy   κ       loopback accuracy
# 164k   12           3x3        2x10               51m         0.979      0.900   1.000
# 308k   12           3x3        2x10               2.4h        0.980      0.924   1.000
# 308k   13           3x3        2x10               7.3h        0.980      0.922   1.000
maximaTraining = bind_rows(get_treetop_predictors(get_treetop_eligible_maxima("s04200w06810"), acceptedTileName = "s04200w06810"), 
                           get_treetop_predictors(get_treetop_eligible_maxima("s04230w06810"), acceptedTileName = "s04230w06810")) %>% 
  select(-tile, -id, -sourceID, -neighbor1sourceID, -neighbor2sourceID, # drop unique identifiers and non-portable values
         -dsmZ, -cmmZ, -ring1max, -ring2max, -ring3max, -ring4max, -ring5max, -ring1mean, -ring2mean, -ring3mean, -ring4mean, -ring5mean, -ring1min, -ring2min, -ring3min, -ring4min, -ring5min, -neighbor1dsmZ, -neighbor2dsmZ) # drop elevations
sum(is.na(maximaTraining)) # check for incomplete cases

# s04200w06810 only
#predictorVariables = c("treetop", "prominenceStdDevNormalized", "height", "radius", "prominence1normalized", "prominence3normalized", "netProminenceNormalized", "prominence4normalized", "prominence5normalized", "netRange", "range1", "deltaCmm", "ring2max")
# s04200w06810 and s04230w06810
#predictorVariables = c("treetop", "height", "prominenceStdDevNormalized", "prominence2normalized", "radius", "prominence2", "prominence3normalized", "prominence1normalized", "netProminenceNormalized", "prominence4normalized", "netRange", "range2", "range1", "prominence1", "ring1min" )
predictorVariables = c("treetop", "height", "radius", "prominence2normalized", "prominence3normalized", "netProminenceNormalized", "prominence1normalized", "ring2variance", "prominence4normalized", "mean2delta", "ring1variance", "netRange", "mean5delta", "neighbor1distance", "deltaCmm")

# tuneRanger(): s04200w06810 only @ 48->12 predictors: mtry = 6, min node size = 2
#               s04200w06810 and s04230w06810 @ 48->12: also mtry = 6, min node size = 2, ~2h train time
#                                               78->13: mtry = 8, min node size = 4, h train time 
fitStart = Sys.time() 
randomForestFit = train(treetop ~ ., data = maximaTraining %>% select(all_of(predictorVariables)), 
                        method = "ranger", # importance = "impurity_corrected", # see https://github.com/imbs-hl/ranger/issues/664 for impurity corrected
                        trControl = trainControl(method = "repeatedcv", number = 2, repeats = 10, verboseIter = TRUE),
                        tuneGrid = expand.grid(mtry = c(8, 9, 10),
                                               splitrule = "gini",
                                               min.node.size = c(3, 4, 5)))
(randomForestFitTime = Sys.time() - fitStart)
randomForestFit
saveRDS(randomForestFit, "trees/segmentation/treetopRandomForest vsurf 48.14 s04200+s04230w06810.Rds")
#varImp(randomForestFit)

randomForestConfusionMatrix = confusionMatrix(predict(randomForestFit), maximaTraining$treetop)
randomForestConfusionMatrix$table
randomForestConfusionMatrix$overall

if (treetopOptions$includeInvestigatory)
{
  # rows   predictors   VSURF  cores selected  accuracy   tune  trees  threads   mtry  min node size
  # 164k   48            3.0h  14    12        ~98%       2.0h  1000   12
  # 308k   48            6.3h  14    12        ~98%       5.4h  1000   12
  # 308k   78           15.3h  14    13        99.1%      7.4h   500   15        8     5
  library(forcats)
  library(VSURF)
  treetopVsurf = VSURF(treetop ~ ., maximaTraining, ncores = 14, parallel = TRUE, RFimplem = "ranger")
  saveRDS(treetopVsurf, "trees/segmentation/treetop vsurf 78 s04200+s04230w06810.Rds")
  treetopVsurf = readRDS("trees/segmentation/treetop vsurf 78 s04200+s04230w06810.Rds")
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
  ggsave("trees/segmentation/treetop vsurf 78 importance s04200+s04230w06810.png", width = 14, height = 35, units = "cm", dpi = 150)

  library(tuneRanger)
  library(mlr)
  rangerTuneTask = makeClassifTask(data = as.data.frame(maximaTraining %>% select(all_of(predictorVariables))), target = "treetop")
  estimateStart = Sys.time()
  estimateTimeTuneRanger(rangerTuneTask, num.trees = 500, num.threads = 14, iters = 70)# 1.5 min to estimate 1h47m @ 164k tops and 12 predictors
  Sys.time() - estimateStart
  tuneStart = Sys.time()
  rangerTuning = tuneRanger(rangerTuneTask, measure = list(multiclass.brier), num.trees = 500, num.threads = 14, iters = 70)
  Sys.time() - tuneStart
  (rangerTuning)
}


## predict treetops, merge, and noise points
#randomForestFit = readRDS("trees/segmentation/treetopRandomForest vsurf 48.14 s04200+s04230w06810.Rds")
randomForestFit = readRDS("trees/segmentation/treetopRandomForest vsurf 78.13 s04200+s04230w06810.Rds")
tileName = "s04230w06840"
tileMaxima = get_treetop_eligible_maxima(tileName)
tileMaxima$treetop = predict(randomForestFit$finalModel, get_treetop_predictors(tileMaxima))$predictions
table(tileMaxima$treetop)

# group merge points and define merge treetops
# TODO: reduce upper merge point radii, consider fusion of single merge points with nearby local maxima
# TODO: reach into adjacent tiles and pull merge points needed for edge treetop completion so that edge treetops aren't generated from incomplete merge sets
tileMergePoints = subset(tileMaxima, tileMaxima$treetop == "merge")
tileMergePointsXY = geom(tileMergePoints)[, c("x", "y")]
mergePointKnn = get.knn(tileMergePointsXY, k = 3) # could also use cutree(hclust(dist(geom()))) but that's O(N²)
mergePointKnn = tibble(tile = tileMergePoints$tile, id = tileMergePoints$id, sourceID = tileMergePoints$sourceID,
                       x = tileMergePointsXY[, "x"], y = tileMergePointsXY[, "y"], 
                       radius = tileMergePoints$radius, dsmZ = tileMergePoints$dsmZ, cmmZ = tileMergePoints$cmmZ, height = tileMergePoints$height,
                       mergeRadius = sqrt(tileMergePoints$height),
                       neighbor1index = mergePointKnn$nn.index[, 1], neighbor1distance = mergePointKnn$nn.dist[, 1],
                       neighbor2index = mergePointKnn$nn.index[, 2], neighbor2distance = mergePointKnn$nn.dist[, 2],
                       neighbor3index = mergePointKnn$nn.index[, 3], neighbor3distance = mergePointKnn$nn.dist[, 3],
                       neighbors = (neighbor1distance < mergeRadius) + (neighbor2distance < mergeRadius) + (neighbor3distance < mergeRadius),
                       centroidX = if_else(neighbors == 0, x, 
                                           if_else(neighbors == 1, 1/2 * (x + tileMergePointsXY[neighbor1index, "x"]), 
                                                   if_else(neighbors == 2, 1/3 * (x + tileMergePointsXY[neighbor1index, "x"] + tileMergePointsXY[neighbor2index, "x"]), 
                                                           1/4 * (x + tileMergePointsXY[neighbor1index, "x"] + tileMergePointsXY[neighbor2index, "x"] + tileMergePointsXY[neighbor3index, "x"])))),
                       centroidY = if_else(neighbors == 0, y, 
                                           if_else(neighbors == 1, 1/2 * (y + tileMergePointsXY[neighbor1index, "y"]), 
                                                   if_else(neighbors == 2, 1/3 * (y + tileMergePointsXY[neighbor1index, "y"] + tileMergePointsXY[neighbor2index, "y"]), 
                                                           1/4 * (y + tileMergePointsXY[neighbor1index, "y"] + tileMergePointsXY[neighbor2index, "y"] + tileMergePointsXY[neighbor3index, "y"]))))) %>% 
  group_by(centroidX, centroidY) %>%
  mutate(mergeID = as.integer(1000000) + cur_group_id()) %>%
  ungroup()
mergeTreetops = mergePointKnn %>% # leave ring statistics unpopulated as they're not well defined for merge points
  group_by(mergeID) %>%
  summarize(tile = str_c(unique(tile), sep = "+"), 
            id = mergeID[1], 
            sourceID = if_else(length(unique(sourceID)) == 1, sourceID[1], as.integer(0)), 
            x = centroidX[1],
            y = centroidY[1],
            radius = max(radius), dsmZ = mean(dsmZ), cmmZ = mean(cmmZ), height = mean(height),
            maxima = n(),
            sourceIDs = length(unique(sourceID)), .groups = "drop")
mergeTreetops = vect(mergeTreetops, geom = c("x", "y"), crs = crs(tileMergePoints))

tileTreetops = subset(tileMaxima, tileMaxima$treetop == "yes")
tileTreetops$maxima = as.integer(1)
tileTreetops$sourceIDs = as.integer(1)
writeVector(rbind(tileTreetops, mergeTreetops), file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "treetops", insert = TRUE, overwrite = TRUE)

writeVector(tileMergePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "merge points", insert = TRUE, overwrite = TRUE)

noisePoints = subset(tileMaxima, tileMaxima$treetop == "noise")
writeVector(noisePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "noise points", insert = TRUE, overwrite = TRUE)
maybeNoisePoints = subset(tileMaxima, tileMaxima$treetop == "maybe noise")
writeVector(maybeNoisePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "maybe noise points", insert = TRUE, overwrite = TRUE)


## diff two treetop identifications
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
    mutate(x = 0.3048 * x, y = 0.3048 * y, z = 0.3048 * z, height = 0.3048 * height, mergeDistance = 0.3048 * mergeDistance, maximaDistance = 0.3048 * maximaDistance)
  
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
  #  geom_histogram(aes(x = 0.3048 * distance), treetopDiagnosticsKnn %>% filter(distance > 1E-6), binwidth = 1) +
  #  labs(x = "accepted treetop to nearest recorded local maxima, m", y = "treetops")
}
