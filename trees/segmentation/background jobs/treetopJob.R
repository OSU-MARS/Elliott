jobStartTime = Sys.time()
source("trees/segmentation/treetops.R")

localMaximaFileNames = list.files(localMaximaPath, "\\.gpkg$")

chunkIndex = 4 # chunks in progress: 1, 2
chunkSize = 128

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(localMaximaFileNames))
localMaximaFileNames = localMaximaFileNames[startIndex:endIndex]

treetopRandomForest = readRDS("trees/segmentation/treetopRandomForest vsurf 103.14 s04200w06840 + s4200+s04230w06810.Rds")

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(localMaximaFileNames), " tiles..."))
for (localMaximaFileName in localMaximaFileNames)
{
  tileName = tools::file_path_sans_ext(localMaximaFileName)
  treetopsFilePath = file.path(candidateTreetopsPath, "rf", paste0(tileName, ".gpkg"))
  if (file.exists(treetopsFilePath))
  {
    next
  }
  cat(paste0(tileName, "...\n"))
  
  # classify local maxima
  # See same code in treetops.R for comments.
  tileMaxima = get_treetop_eligible_maxima(tileName) %>% filter(is.na(cmmSlope3) == FALSE, is.na(ring4mean) == FALSE) # s03870w06600 and s04020w07230, at least, have ring4mean NAs that fail random forest prediction on ring4delta
  tileMaxima$treetop = predict(treetopRandomForest$finalModel, get_treetop_predictors(tileMaxima))$predictions

  tileNeighborhood = get_treetop_eligible_neighborhood(tileName, tileMaxima)
  neighborhoodMaxima = tileNeighborhood
  tileMergePoints = get_merge_points(tileMaxima, tileNeighborhood)
  tileMaxima$treetop[tileMergePoints$mergePointIndices] = "merge"
  tileMaxima$treetop[tileMergePoints$treetopIndices] = "yes"

  # write tile's treetop GeoPackage
  tileTreetops = bind_rows(tileMaxima %>% filter(tileMaxima$treetop == "yes"), tileMergePoints$treetops) %>%
    mutate(maxima = replace_na(maxima, as.integer(1)), sourceIDs = replace_na(maxima, as.integer(1)))
  writeVector(vect(tileTreetops, crs = tileCrs, geom = c("x", "y")), treetopsFilePath, layer = "treetops", insert = TRUE, overwrite = TRUE)
  
  writeVector(vect(tileMaxima %>% filter(tileMaxima$treetop == "merge"), crs = tileCrs, geom = c("x", "y")), treetopsFilePath, layer = "merge points", insert = TRUE, overwrite = TRUE)
  
  noisePoints = vect(tileMaxima %>% filter(treetop == "noise"), crs = tileCrs, geom = c("x", "y"))
  writeVector(noisePoints, treetopsFilePath, layer = "noise points", insert = TRUE, overwrite = TRUE)
  maybeNoisePoints = vect(tileMaxima %>% filter(treetop == "maybe noise"), crs = tileCrs, geom = c("x", "y"))
  writeVector(maybeNoisePoints, treetopsFilePath, layer = "maybe noise points", insert = TRUE, overwrite = TRUE)
}

warnings()
print(paste0("treetop classification ran for ", format(Sys.time() - jobStartTime), "."))