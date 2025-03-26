jobStartTime = Sys.time()
source("trees/segmentation/treetops.R")
treetopOptions$setTreetopStandIDs = FALSE

localMaximaFileNames = list.files(localMaximaChmCmmPath, "\\.gpkg$")

# TODO: investigate future_map() instead of chunking as 12 workers is probably fine
chunkIndex = 5 # 561 tiles -> chunks, ~2 GB DDR @ 5.5 GB/s peak per job
chunkSize = 128 # ~40 minutes/chunk with five concurrent jobs, 9900X

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(localMaximaFileNames))
localMaximaFileNames = localMaximaFileNames[startIndex:endIndex]

treetopRandomForest = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m9n3.Rds")
treetopsPathRandomForest = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/rf v1"

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(localMaximaFileNames), " tiles..."))
for (localMaximaFileName in localMaximaFileNames)
{
  tileName = tools::file_path_sans_ext(localMaximaFileName)
  treetopsFilePathRandomForest = file.path(treetopsPathRandomForest, paste0(tileName, ".gpkg"))
  if (file.exists(treetopsFilePathRandomForest))
  {
    next
  }
  cat(paste0(tileName, "...\n"))

  # load tile, also neighborhood if density predictors are used
  tileMaxima = get_treetop_eligible_maxima(tileName) %>% filter(is.na(cmmSlope3) == FALSE) #, is.na(ring4mean) == FALSE) # s03870w06600 and s04020w07230, at least, have ring4mean NAs that fail random forest prediction on ring4delta
  tileNeighborhood = get_treetop_eligible_neighborhood(tileName, tileMaxima)
  tileIndices = which(tileNeighborhood$tile == tileName)

  # classify local maxima
  tileMaxima$treetop = predict(treetopRandomForest, get_treetop_predictors(tileNeighborhood) %>% filter(tile == tileName), num.threads = treetopOptions$rangerThreads)$predictions
  tileNeighborhood$treetop = factor(NA, levels = levels(tileMaxima$treetop))
  tileNeighborhood$treetop[tileIndices] = tileMaxima$treetop

  tileMergePoints = get_merge_points(tileMaxima, tileNeighborhood)
  tileMaxima$treetop[tileMergePoints$mergePointIndices] = "merge"
  tileMaxima$treetop[tileMergePoints$treetopIndices] = "yes"

  # write tile's treetop GeoPackage
  tileTreetopPoints = st_as_sf(bind_rows(tileMaxima %>% filter(treetop != "no") %>% rename(treeID = id) %>% select(-mergeClusterID, -sourceID, -uniqueID, -uniqueMergeClusterID),
                                         tileMergePoints$treetops %>% select(-mergeClusterNumber, -sourceID, -sourceIDs) %>% mutate(treetop = factor("yes"))) %>%
                                 mutate(mergePoints = replace_na(mergePoints, as.integer(1))), 
                               coords = c("x", "y"), crs = attributes(tileMaxima)$crs, sf_column_name = "geom") # crs attribute set by get_treetop_eligible_maxima()

  st_write(tileTreetopPoints %>% filter(treetop == "yes"), treetopsFilePathRandomForest, layer = "treetops", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
  st_write(tileTreetopPoints %>% filter(treetop == "merge"), treetopsFilePathRandomForest, layer = "merge points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
  st_write(tileTreetopPoints %>% filter(treetop == "noise"), treetopsFilePathRandomForest, layer = "noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
  if (any(tileTreetopPoints$treetop == "maybe noise"))
  {
    # debatable whether to write an empty layer or omit the layer entirely, for now complete omission is used
    st_write(tileTreetopPoints %>% filter(treetop == "maybe noise"), treetopsFilePathRandomForest, layer = "maybe noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
  }
}

warnings()
cat(paste0("treetop clustered random forest classification ran for ", format(Sys.time() - jobStartTime), "."))


## join stand IDs to random forest treetops
if (treetopOptions$setTreetopStandIDs)
{
  handlers(global = TRUE)
  handlers("cli")
  plan(multisession, workers = 0.5 * future::availableCores()) # no gain for form selection with vfold_cv() but effective for best fit searches and classifying treetops in tiles
  
  stands2016 = st_transform(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", quiet = TRUE, layer = "unified stands 2016"),
                            make_compound_crs(6557, 8228)) %>% # keep in sync with same code in decision boundary.R
    select(standID2016)
  
  treetopTilePaths = list.files(treetopsPathRandomForest, "\\.gpkg$", full = TRUE)
  with_progress({
    future_map(treetopTilePaths, function(treetopTilePath)
    {
      tileTreetops = st_join(st_read(treetopTilePath, quiet = TRUE, layer = "treetops"), stands2016, left = TRUE)
      st_write(tileTreetops, treetopTilePath, delete_dsn = FALSE, delete_layer = TRUE)
    })
  }, .options = furrr_options(seed = TRUE))
}