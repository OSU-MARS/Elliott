library(ranger)
library(stringr)
library(terra)

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
localMaximaPath = file.path(dataPath, "DSM v3 beta/local maxima")
tileNames = list.files(dsmPath, "\\.tif$")

treetopDestinationPath = file.path(dataPath, "treetops/rf")
treetopRandomForest = readRDS("trees/segmentation/treetopRandomForest vsurf 99.13 s04200+s04230w06810+s04230w06840.Rds")

# copy-paste from treetops.R
get_treetop_eligible_maxima = function(tileName, minimumHeight = 0.984252, acceptedTileName = NULL)
{
  localMaxima = vect(file.path(localMaximaPath, paste0(tileName, ".gpkg")))
  # use minimumHeight to exclude groundcover maxima and maxima from sensor noise or error
  # for now, also exclude treetop candidates with so few adjacent points ring 1 or ring 2 has no data since it's 1) it's difficult to tell if these are maxima, 2) it's unlikely they're actually maxima, and 3) these likely comprise < 0.01% of maximas
  localMaxima = subset(localMaxima, (localMaxima$height >= minimumHeight) & (is.na(localMaxima$ring1mean) == FALSE) & (is.na(localMaxima$ring2mean) == FALSE))
  
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
  localMaximaXY = geom(localMaxima)[, c("x", "y")]
  
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

get_treetop_predictors = function(treetopEligibleMaxima)
{
  localMaximaXY = geom(treetopEligibleMaxima)[, c("x", "y")]
  
  neighborKnn = get.knn(localMaximaXY, k = 50)
  neighborKnn = tibble(index1 = neighborKnn$nn.index[, 1], distance1 = neighborKnn$nn.dist[, 1],
                       index2 = neighborKnn$nn.index[, 2], distance2 = neighborKnn$nn.dist[, 2],
                       index3 = neighborKnn$nn.index[, 3], distance3 = neighborKnn$nn.dist[, 3],
                       index4 = neighborKnn$nn.index[, 4], distance4 = neighborKnn$nn.dist[, 4],
                       index5 = neighborKnn$nn.index[, 5], distance5 = neighborKnn$nn.dist[, 5],
                       meanDistance2 = rowMeans(neighborKnn$nn.dist[, 1:2]),
                       meanDistance3 = rowMeans(neighborKnn$nn.dist[, 1:3]),
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
           neighborDistance2Mean = neighborKnn$meanDistance2, neighborDistance3Mean = neighborKnn$meanDistance3, neighborDistance5Mean = neighborKnn$meanDistance5, neighborDistance10Mean = neighborKnn$meanDistance10, neighborDistance20Mean = neighborKnn$meanDistance20, neighborDistance50Mean = neighborKnn$meanDistance50,
           # convert all distances from feet to m
           across(where(is.double), ~0.3048 * .x),
           deltaCmm = dsmZ - cmmZ,
           mean1delta = dsmZ - ring1mean, 
           mean2delta = dsmZ - ring2mean, 
           mean3delta = dsmZ - ring2mean, 
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


for (tileName in tileNames)
{
  treetopFilePath = file.path(dataDestinationPath, tileName)
  if (file.exists(treetopFilePath))
  {
    next
  }
  cat(paste0(str_remove(tileName, "\\.tif"), "...\n"))
  
  # predict local maxima classifications
  tileMaxima = get_treetop_eligible_maxima(tileName)
  tileMaxima$treetop = predict(treetopRandomForest$finalModel, get_treetop_predictors(tileMaxima))$predictions
  
  # write tile geopackage
  # treetops layer, including merge tops
  tileTreetops = subset(tileMaxima, tileMaxima$treetop == "yes")
  tileTreetops$maxima = as.integer(1)
  tileTreetops$sourceIDs = as.integer(1)
  writeVector(rbind(tileTreetops, mergeTreetops), file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "treetops", insert = TRUE, overwrite = TRUE)
  # merge points
  writeVector(tileMergePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "merge points", insert = TRUE, overwrite = TRUE)
  # noise points
  noisePoints = subset(tileMaxima, tileMaxima$treetop == "noise")
  writeVector(noisePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "noise points", insert = TRUE, overwrite = TRUE)
  maybeNoisePoints = subset(tileMaxima, tileMaxima$treetop == "maybe noise")
  writeVector(maybeNoisePoints, file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg")), layer = "maybe noise points", insert = TRUE, overwrite = TRUE)
}
