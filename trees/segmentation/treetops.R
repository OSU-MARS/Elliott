# some investigatory code here assumes segmentation/setup.R
library(caret)
library(dplyr)
library(FNN)
library(ggplot2)
library(magrittr)
library(patchwork)
library(sf)
library(terra)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             panel.border = element_blank(), 
                             plot.title = element_text(size = 10)))

treetopOptions = tibble(includeInvestigatory = FALSE)

localMaximaPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 beta/local maxima"
treetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops accepted"

get_treetop_training_data = function(tileName)
{
  localMaxima = vect(file.path(localMaximaPath, paste0(tileName, ".gpkg")))
  treetops = vect(file.path(treetopsPath, paste0(tileName, ".gpkg")), layer = "treetops")
  mergePoints = vect(file.path(treetopsPath, paste0(tileName, ".gpkg")), layer = "merge treetops")
  
  neighborKnn = get.knn(geom(localMaxima), k = 1)
  neighborKnn = tibble(index = neighborKnn$nn.index[, 1], distance = neighborKnn$nn.dist[, 1])
  
  isTreetopKnn = get.knnx(geom(localMaxima)[, c("x", "y")], geom(treetops)[, c("x", "y")], k = 1)
  isTreetopKnn = tibble(index = isTreetopKnn$nn.index[, 1], distance = isTreetopKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
  
  isMergeKnn = get.knnx(geom(localMaxima)[, c("x", "y")], geom(mergePoints)[, c("x", "y")], k = 1)
  isMergeKnn = tibble(index = isMergeKnn$nn.index[, 1], distance = isMergeKnn$nn.dist[, 1]) %>% filter(distance < 0.1)
  #mergeKnn %>% slice_max(distance, n = 10)
  #localMaxima[105196, ]
  
  treetopTraining = as_tibble(localMaxima) %>% 
    mutate(treetop = factor("no", levels = c("no", "yes", "merge")),
           neighborSourceID1 = sourceID[neighborKnn$index], neighborDistance1 = neighborKnn$distance, neighborDsmZ1 = dsmZ[neighborKnn$index],
           across(where(is.double), ~0.3048 * .x), # convert from feet to m
           neighborDifferentSourceID1 = sourceID != neighborSourceID1, neighborProminence1 = neighborDsmZ1 - dsmZ,
           deltaCmm = dsmZ - cmmZ,
           prominenceRing1 = dsmZ - maxRing1, prominenceRing2 = dsmZ - maxRing2, prominenceRing3 = dsmZ - maxRing3, prominenceRing4 = dsmZ - maxRing4, prominenceRing5 = dsmZ - maxRing5,
           rangeRing1 = maxRing1 - minRing1, rangeRing2 = maxRing2 - minRing2, rangeRing3 = maxRing3 - minRing3, rangeRing4 = maxRing4 - minRing4, rangeRing5 = maxRing5 - minRing5,
           neighborProminenceNormalized1 = neighborProminence1 / height,
           netProminence = prominenceRing1 + prominenceRing2 + prominenceRing3 + prominenceRing4 + prominenceRing5,
           netProminenceNormalized = netProminence / height,
           netRange = rangeRing1 + rangeRing2 + rangeRing3 + rangeRing4 + rangeRing5,
           netRangeNormalized = netRange / height,
           prominenceNormalized1 = prominenceRing1 / height, prominenceNormalized2 = prominenceRing2 / height, prominenceNormalized3 = prominenceRing3 / height, prominenceNormalized4 = prominenceRing4 / height, prominenceNormalized5 = prominenceRing5 / height,
           prominenceStdDev = sd(c_across(prominenceRing1:prominenceRing5)),
           prominenceStdDevNormalized = prominenceStdDev / height,
           rangeNormalized1 = rangeRing1 / height, rangeNormalized2 = rangeRing2 / height, rangeNormalized3 = rangeRing3 / height, rangeNormalized4 = rangeRing4 / height, rangeNormalized5 = rangeRing5 / height,
           rangeStdDev = sd(c_across(rangeRing1:rangeRing5)),
           rangeStdDevNormalized = rangeStdDev / height)
  treetopTraining$treetop[isTreetopKnn$index] = "yes"
  treetopTraining$treetop[isMergeKnn$index] = "merge"
  
  return(treetopTraining)
}
  
## treetop random forests
maximaTraining = get_treetop_training_data("s04200w06810") %>% # bind_rows(get_treetop_training_data("s04200w06810"), get_treetop_training_data("s04230w06810")) %>% 
  select(-tile, -id, -sourceID, -neighborSourceID1)

predictorVariables = c("treetop", "prominenceStdDevNormalized", "height", "radius", "prominenceNormalized1", "prominenceNormalized3", "netProminenceNormalized", "prominenceNormalized4", "prominenceNormalized5",
                       "netRange", "rangeRing1", "deltaCmm", "maxRing2")

fitStart = Sys.time() # ~1.6 minutes @ 25k rows, k = 10, 3x3 tune grid, 8.3 minutes @ 75 k rows, 33 minutes @ 154k, gini accuracy 99.4% κ = 0.972
randomForestFit = train(treetop ~ ., data = maximaTraining %>% select(predictorVariables), 
                        method = "ranger", # importance = "impurity_corrected", # see https://github.com/imbs-hl/ranger/issues/664 for impurity corrected
                        trControl = trainControl(method = "repeatedcv", number = 2, repeats = 10, verboseIter = TRUE),
                        tuneGrid = expand.grid(mtry = c(9, 10, 11),
                                               splitrule = "gini",
                                               min.node.size = c(2, 4, 8, 16)))
(randomForestFitTime = Sys.time() - fitStart)
randomForestFit
varImp(randomForestFit)

randomForestConfusionMatrix = confusionMatrix(randomForestPrediction$predictedClass, randomForestPrediction$actualClass)
randomForestConfusionMatrix$table
randomForestConfusionMatrix$overall


if (treetopOptions$includeInvestigatory)
{
  library(forcats)
  library(VSURF) # 3 h @ 48 predictors and 164 k rows
  treetopVsurf = VSURF(treetop ~ ., maximaTraining, ncores = 12, parallel = TRUE, RFimplem = "ranger")
  saveRDS(treetopVsurf, "trees/segmentation/treetop vsurf 49 s04200w06810.Rds")
  treetopVsurf$nums.varselect # threshold -> interpretation -> prediction
  treetopVsurf$mean.perf # fractional error rate
  treetopVsurf$overall.time
  treetopVsurf$comput.times
  
  plot(treetopVsurf)
  attributes(treetopVsurf$terms[treetopVsurf$varselect.interp])$term.labels
  attributes(treetopVsurf$terms[treetopVsurf$varselect.pred])$term.labels
  
  predictorImportance = tibble(responseVariable = "classification", predictor = as.character(attr(treetopVsurf$terms, "predvars"))[treetopVsurf$imp.mean.dec.ind + 2], importance = treetopVsurf$imp.mean.dec) %>% # offset as.character() by two since first element is "list" and second is classification
    mutate(importance = importance / sum(importance))
  ggplot() +
    geom_col(aes(x = importance, y = fct_reorder(predictor, importance), fill = responseVariable), predictorImportance) +
    guides(fill = "none") +
    labs(x = "normalized variable importance", y = NULL)  
  ggsave("trees/segmentation/treetop vsurf 49 importance s04200w06810.png", width = 11, height = 25, units = "cm", dpi = 150)

  library(tuneRanger)
  library(mlr)
  rangerTuneTask = makeClassifTask(data = as.data.frame(maximaTraining %>% select(all_of(predictorVariables))), target = "treetop")
  estimateStart = Sys.time()
  estimateTimeTuneRanger(rangerTuneTask, num.trees = 1000, num.threads = 12, iters = 70) # 1.5 min to estimate 1h47m
  Sys.time() - estimateStart
  tuneStart = Sys.time()
  rangerTuning = tuneRanger(rangerTuneTask, measure = list(multiclass.brier), num.trees = 1000, num.threads = 12, iters = 70)
  Sys.time() - tuneStart
}


## form accepted treetop set for tile
dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileName = "s04230w06810"
detectedTreetops = vect(file.path(dataPath, "DSM v3 beta", paste0(tileName, ".gpkg")), layer = "treetops")
acceptedTreetops = vect(file.path(dataPath, "treetops accepted", paste0(tileName, ".gpkg")), layer = "treetops")
acceptedTreetopKnn = get.knnx(geom(acceptedTreetops)[, c("x", "y")], geom(detectedTreetops)[, c("x", "y")], k = 1)
acceptedTreetopKnn = tibble(acceptedTreetopIndex = acceptedTreetopKnn$nn.index, distanceToAcceptedTreetop = acceptedTreetopKnn$nn.dist)
detectedTreetops$distanceToAcceptedTreetop = acceptedTreetopKnn$distanceToAcceptedTreetop
detectedTreetops$detectionType = if_else(detectedTreetops$distanceToAcceptedTreetop < 1E-3, "correct", "commission")

omissionKnn = get.knnx(geom(detectedTreetops)[, c("x", "y")], geom(acceptedTreetops)[, c("x", "y")], k = 1)
omissionKnn = tibble(detectedTreetopIndex = omissionKnn$nn.index[, 1], distanceToAcceptedTreetop = omissionKnn$nn.dist[, 1])
omittedTreetops = omissionKnn %>% filter(distanceToAcceptedTreetop > 1E-3)

# goodness of treetop detection
# tile          method            treetops correct omission commission accuracy treetopAccuracy omissionError commissionError kappa
# s04200w06810  initial npn 0.02  18406    17901   505      19         1.00     0.973           0.0274        0.00103         0.986
# s04200w06810  inner 2 npn 0.02  18402    17846   550      2018       0.999    0.970           0.0299        0.110           0.933
tileCells = (3000 / 1.5)^2
tibble(treetops = nrow(acceptedTreetops), correct = sum(detectedTreetops$detectionType == "correct"), omission = nrow(omittedTreetops), commission = sum(detectedTreetops$detectionType == "commission")) %>% 
  mutate(accuracy = (tileCells - omission - commission) / tileCells, treetopAccuracy = correct / treetops, omissionError = omission / treetops, commissionError = commission / treetops,
         kappa = 2 * (correct * (tileCells - treetops) - omission * commission) / ((correct + commission) * (commission + (tileCells - treetops)) + (correct + omission) * (omission + (tileCells - treetops))))

# local maxima to treetop conversion
localMaxima = vect(file.path(dataPath, "DSM v3 beta/local maxima", paste0(tileName, ".gpkg"))) # EPSG:6557
localMaxima$distanceToAcceptedTreetop = NA_real_
isTreetopKnn = get.knnx(geom(localMaxima)[, c("x", "y")], geom(acceptedTreetops)[, c("x", "y")], k = 1) # both EPSG:6557
isTreetopKnn = tibble(localMaximaIndex = isTreetopKnn$nn.index[, 1], distanceToAcceptedTreetop = isTreetopKnn$nn.dist[, 1])
localMaxima$distanceToAcceptedTreetop[isTreetopKnn$localMaximaIndex] = round(isTreetopKnn$distanceToAcceptedTreetop, 3)

localMaxima = as_tibble(localMaxima) %>%
  mutate(isTreetop = replace_na((distanceToAcceptedTreetop < 1E-3) | ((inSimilarElevationGroup == 1) & (distanceToAcceptedTreetop < 3)), FALSE),
         netProminence = netProminenceNormalized * height)

table(localMaxima$distanceToAcceptedTreetop)
geom(localMaxima)[which((localMaxima$distanceToAcceptedTreetop >= 0.612) & (localMaxima$distanceToAcceptedTreetop <= 0.693)), c("x", "y")]


#sum(localMaxima$distanceToAcceptedTreetop < 1E-3, na.rm = TRUE)
#sum(localMaxima$distanceToAcceptedTreetop >= 1E-3, na.rm = TRUE)
#localMaxima %>% summarize(maxima = n(), netProminence000 = sum(netProminenceNormalized < 0), netProminence001 = sum((netProminenceNormalized >= 0) & (netProminenceNormalized < 0.01)), netProminence002 = sum((netProminenceNormalized >= 0.01) & (netProminenceNormalized < 0.02)), netProminenceGt002 = sum(netProminenceNormalized >= 0.02))

localMaxima %<>% mutate(rings = 1 + is.na(maxRing2) + is.na(maxRing3) + is.na(maxRing4) + is.na(maxRing5),
                        netProminence1 = dsmZ - maxRing1, # netProminenceNormalized = (netProminence1 + 2 + ...) / height
                        netProminence2 = dsmZ - maxRing2,
                        netProminence3 = dsmZ - maxRing3,
                        netProminence4 = dsmZ - maxRing4,
                        netProminence5 = dsmZ - maxRing5,
                        absoluteProminence = if_else(is.na(netProminence1), 0, abs(netProminence1)) + if_else(is.na(netProminence2), 0, abs(netProminence2)) + if_else(is.na(netProminence3), 0, abs(netProminence3)) + if_else(is.na(netProminence4), 0, abs(netProminence4)) + if_else(is.na(netProminence5), 0, abs(netProminence5)),
                        outerProminence = if_else(is.na(netProminence3), 0, netProminence3) + if_else(is.na(netProminence4), 0, netProminence4) + if_else(is.na(netProminence5), 0, netProminence5),
                        outerProminenceNormalized = outerProminence / height)
localMaxima %<>% mutate(isTreetopPredicted = (radius >= pmax(0.5 * 3.28048 * pmin(0.045 * 0.3048 * height + 0.5, 4.0), 1.5)) & (netProminenceNormalized > 0.02))
tibble(accepted = nrow(acceptedTreetops), detected = nrow(detectedTreetops), matchedLocalMaxima = sum(localMaxima$isTreetop), predictedTreetops = sum(localMaxima$isTreetopPredicted), correct = sum(localMaxima$isTreetop & (localMaxima$isTreetop == localMaxima$isTreetopPredicted)))

ggplot() +
  geom_abline(color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = 0.02, color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_bin2d(aes(x = netProminenceNormalized, y = outerProminenceNormalized), localMaxima, binwidth = 0.01, na.rm = TRUE) + 
  labs(x = "net prominence, normalized", y = "outer prominence, normalized", fill = "n", title = "(a) maxima") +
ggplot() +
  geom_abline(color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_vline(xintercept = 0.02, color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_bin2d(aes(x = netProminenceNormalized, y = outerProminenceNormalized), localMaxima %>% filter(isTreetop), binwidth = 0.01, na.rm = TRUE) + 
  labs(x = "net prominence, normalized", y = NULL, fill = "n", title = "(b) accepted treetops") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) & # xlim = c(-2, 4), ylim = c(-2, 2.5))
  scale_fill_viridis_c(limits = c(1, 15000), trans = "log10")

ggplot() +
  geom_point(aes(x = netProminence, y = radius, color = isTreetop), localMaxima, alpha = 0.01, shape = 16) +
  labs(x = "net prominence, normalized", y = "radius, ft", color = "treetops")

localMaximaToTreetops = localMaxima %>% filter(radius > 3) %>%
  mutate(netProminenceBin = 0.01 * round(netProminenceNormalized / 0.01)) %>% group_by(netProminenceBin) %>%
  summarize(n = n(), pTreetop = sum(isTreetop) / n)

ggplot() +
  geom_segment(aes(x = 0.02, y = 0, xend = 0.02, yend = 1), color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_line(aes(x = netProminenceBin, y = pTreetop), localMaximaToTreetops) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.01)) +
  labs(x = "net prominence, normalized", y = "treetop probability") +
  scale_x_continuous(breaks = c(-10, -1, -0.1, 0, 0.1, 1, 10)) +
  scale_y_continuous(labels = scales::label_percent())

ggplot() +
  geom_segment(aes(x = 0.02, y = 0, xend = 0.02, yend = 3000), color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_line(aes(x = netProminenceBin, y = n), localMaximaToTreetops) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.01), y = scales::transform_pseudo_log()) +
  labs(x = "net prominence, normalized", y = "local maxima") +
  scale_x_continuous(breaks = c(-10, -1, -0.1, 0, 0.1, 1, 10)) +
  scale_y_continuous(breaks = c(0, 1, 10, 100, 1000, 10000))



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


# treetop SVM
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


## modified normalized net prominence requirement
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
