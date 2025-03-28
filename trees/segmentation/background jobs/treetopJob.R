source("trees/segmentation/treetops.R")

handlers(global = TRUE)
handlers("cli")
plan(multisession, workers = 0.5 * future::availableCores()) # workers mostly run single threaded in sf, so fine to also use treetopOptions$rangerThreads = half cores

get_tile_by_height_class = function(tileName, tileTreetops, tileMergePoints, tileNoisePoints, tileMaybeNoisePoints = NULL)
{
  tileByHeightClass = full_join(full_join(st_drop_geometry(tileTreetops) %>% mutate(heightClassInM = round(height)) %>% group_by(standID2016, heightClassInM) %>%
                                            summarize(treetops = n(), .groups = "drop"),
                                          st_drop_geometry(tileMergePoints)  %>% mutate(heightClassInM = round(height)) %>% group_by(standID2016, heightClassInM) %>%
                                            summarize(mergePoints = n(), .groups = "drop"),
                                          by = join_by(standID2016, heightClassInM)),
                                st_drop_geometry(tileNoisePoints) %>% mutate(heightClassInM = round(height)) %>% group_by(standID2016, heightClassInM) %>%
                                  summarize(noisePoints = n(), .groups = "drop"),
                                by = join_by(standID2016, heightClassInM))
  if (is.null(tileMaybeNoisePoints) == FALSE)
  {
    tileByHeightClass %<>% full_join(st_drop_geometry(tileMaybeNoisePoints) %>% mutate(heightClassInM = round(height)) %>% group_by(standID2016, heightClassInM) %>%
                                       summarize(maybeNoisePoints = n(), .groups = "drop"),
                                     by = join_by(standID2016, heightClassInM))
  } else {
    tileByHeightClass %<>% mutate(maybeNoisePoints = 0)
  }
  
  return(tileByHeightClass %>% mutate(method = factor("DSM forest"), tile = tileName,
                                      treetops = replace_na(treetops, 0), mergePoints = replace_na(mergePoints, 0), noisePoints = replace_na(noisePoints, 0), maybeNoisePoints = replace_na(maybeNoisePoints, 0)) %>%
           relocate(method, tile, standID2016, heightClassInM))
}

localMaximaFileNames = list.files(localMaximaPathV3, "\\.gpkg$")

stands2016 = st_transform(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", quiet = TRUE, layer = "unified stands 2016"),
                          make_compound_crs(6557, 8228)) %>% # for DSM v3, keep in sync with same code in treetopJob.R
                          # st_crs(6557)) %>% # for runs against DSM v3 beta
  select(standID2016)
treetopRandomForest = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m9n3.Rds") 
forestTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/rf v1"

treetopStartTime = Sys.time()
with_progress({
  progressBar = progressor(steps = length(localMaximaFileNames))
  
  forestTreetops = bind_rows(future_map(localMaximaFileNames, function(localMaximaFileName) # 46.04m with 12 workers, 9900X
  {
    require(ranger) # ranger apparently doesn't flow to workers for some reason, so predict() fails to resolve without this
    
    tileName = tools::file_path_sans_ext(localMaximaFileName)
    tileForestTreetopsPath = file.path(forestTreetopsPath, paste0(tileName, ".gpkg"))
    if (file.exists(tileForestTreetopsPath))
    {
      return(NULL)
    }
    cat(paste0(tileName, "...\n"))
  
    # load tile, also neighborhood if density predictors are used
    tileMaxima = get_treetop_eligible_maxima(tileName) %>% filter(is.na(cmmSlope3) == FALSE) #, is.na(ring4mean) == FALSE) # s03870w06600 and s04020w07230, at least, have ring4mean NAs that fail random forest prediction on ring4delta
    tileNeighborhood = get_treetop_eligible_neighborhood(tileName, tileMaxima)
    tileIndices = which(tileNeighborhood$tile == tileName)
  
    # classify local maxima
    tileMaxima$treetop = predict(treetopRandomForest, get_treetop_predictors(tileNeighborhood) %>% filter(tile == tileName), num.threads = treetopOptions$rangerThreads)$predictions
    tileNeighborhood$treetop = factor(NA, levels = levels(tileMaxima$treetop)) # needed by get_merge_points()
    tileNeighborhood$treetop[tileIndices] = tileMaxima$treetop
  
    # cluster treetops and update local maxima classifications
    tileMergePoints = get_merge_points(tileMaxima, tileNeighborhood)
    tileMaxima$treetop[tileMergePoints$singleTreetops$tileIndex] = "yes"
    tileMaxima$treetop[tileMergePoints$mergePoints$tileIndex] = "merge"
    tileMaxima$mergeClusterNumber = NA_integer_
    tileMaxima$mergeClusterNumber[tileMergePoints$mergePoints$tileIndex] = tileMergePoints$mergePoints$mergeClusterNumber

    if (sum(tileMaxima$treetop == "merge") != nrow(tileMergePoints$mergePoints))
    {
      stop("Internal consistency failure. Expected merge point clustering and single treetop revisions to result in ", nrow(tileMergePoints$mergePoints), " merge points on tile ", tileName, " but ", sum(tileMaxima$treetop == "merge"), " merge points are defined after clustering.")
      #tileMaxima[setdiff(which(tileMaxima$treetop == "merge"), tileMergePoints$mergePoints$tileIndex), ] %>% select(-starts_with("ring"), -geom)
    }
    
    # write tile's treetop GeoPackage
    tileCrs = st_crs(attributes(tileMaxima)$crs)
    tileTreetopPoints = st_join(st_as_sf(bind_rows(tileMaxima %>% filter(treetop != "no") %>% rename(treeID = id) %>% select(-sourceID, -uniqueID, -uniqueMergeClusterID, -starts_with("ring"), -dsmSlope, -cmmSlope3),
                                                   tileMergePoints$mergeTreetops),
                                         coords = c("x", "y"), crs = tileCrs, sf_column_name = "geom"), # crs attribute set by get_treetop_eligible_maxima()
                                stands2016, left = TRUE)
    tileTreetopsToWrite = tileTreetopPoints %>% filter(treetop %in% c("yes", "merge treetop")) %>% 
      select(-treetop, -mergeClusterID, -mergeClusterNumber) %>% # no merge cluster ID update needed as it's the same as treeID
      mutate(mergePoints = replace_na(mergePoints, as.integer(1)), mergePointsOnTile = replace_na(mergePointsOnTile, as.integer(1))) %>%
      relocate(tile, treeID, height, radius, cmmZ, dsmZ, mergePoints, mergePointsOnTile, standID2016)
    
    tileByHeightClass = get_tile_by_height_class(tileName,
                                                 tileTreetopsToWrite, # 1 m height classes are metric, so get summary before English unit conversion (if applicable)
                                                 tileTreetopPoints %>% filter(treetop == "merge"), # merge, noise, and maybe noise layers aren't written with heights
                                                 tileTreetopPoints %>% filter(treetop == "noise"), 
                                                 tileTreetopPoints %>% filter(treetop == "maybe noise"))
    #tileByHeightClass %>% summarize(treetops = sum(treetops), mergePoints = sum(mergePoints), noisePoints = sum(noisePoints), maybeNoisePoints = sum(maybeNoisePoints))
    if (tileCrs$units_gdal == "foot")
    {
      # get_treetop_eligible_maxima() converts to metric, needs to converted back for correct write on English CRSes
      # Could also change CRSes to metric but convention is to flow input CRS.
      tileTreetopsToWrite$height = 3.28084 * tileTreetopsToWrite$height
      tileTreetopsToWrite$radius = round(3.28084 * tileTreetopsToWrite$radius, 3) # debatable but, for now, assume whole number preference in surface model cell size
      tileTreetopsToWrite$dsmZ = 3.28084 * tileTreetopsToWrite$dsmZ # debatable if DSM and CMM elevations need to flow, but they're flown for now
      tileTreetopsToWrite$cmmZ = 3.28084 * tileTreetopsToWrite$cmmZ
    }
    
    tileMergePointsToWrite = tileTreetopPoints %>% filter(treetop == "merge") %>% rename(id = treeID, clusterID = mergeClusterID) %>%
      group_by(mergeClusterNumber) %>%
      mutate(clusterID = min(id)) %>%
      ungroup() %>%
      select(id, clusterID) #, mergeClusterNumber)
    tileNoisePointsToWrite = tileTreetopPoints %>% filter(treetop == "noise") %>% rename(id = treeID) %>% select(id)
    tileMaybeNoisePointsToWrite = tileTreetopPoints %>% filter(treetop == "maybe noise") %>% select(treeID) %>% rename(id = treeID)
    
    st_write(tileTreetopsToWrite, tileForestTreetopsPath, layer = "treetops", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(tileMergePointsToWrite, tileForestTreetopsPath, layer = "merge points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(tileNoisePointsToWrite, tileForestTreetopsPath, layer = "noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    if (nrow(tileMaybeNoisePointsToWrite) > 0)
    {
      # debatable whether to write an empty layer or omit the layer entirely, for now complete omission is used
      st_write(tileMaybeNoisePointsToWrite, tileForestTreetopsPath, layer = "maybe noise points", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    }
    
    progressBar()
    return(tileByHeightClass)
  }, .options = furrr_options(seed = TRUE)))
})

warnings()
cat(paste0("treetop clustered random forest classification ran for ", format(Sys.time() - treetopStartTime), "."))

forestTreetops %>% 
  summarize(tiles = length(unique(tile)), stands = length(unique(standID2016)), maxHeightInM = max(heightClassInM), 
            elliottTreetops1m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000), 0, treetops)),
            elliottTreetops5m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000) | (heightClassInM < 5), 0, treetops)), 
            totalTreetops1m = sum(treetops), totalTreetops5m = sum(if_else(heightClassInM >= 5, treetops, 0)))
#writexl::write_xlsx(forestTreetops, file.path(forestTreetopsPath, "standsByHeightClass.xlsx")) # 7.5 MB
#forestTreetops %>% filter(heightClassInM > 85) %>% arrange(desc(heightClassInM))
