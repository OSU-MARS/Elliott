library(stars)
source("trees/segmentation/treetopDetection.R")

handlers(global = TRUE)
handlers("cli")
# workers mostly run single threaded in sf, so fine to also use treetopOptions$rangerThreads = half cores
# ~1.5 GB DDR/worker average, ~3 GB peak
options(future.globals.maxSize = 2 * 1024^3) # needed for physiographic rasters (slope, aspect, topographic shelter)
plan(multisession, workers = 0.5 * future::availableCores())

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

# accepted tiles don't need treetop prediction but do need physiographic variables attached
acceptedTreetopFileNames = setdiff(list.files(acceptedTreetopsDsmPath, "\\.gpkg$"), "s04050w06930.gpkg") # s04050w06930 not yet complete
localMaximaFileNames = setdiff(list.files(localMaximaPathV3, "\\.gpkg$"), acceptedTreetopFileNames)

stands2016 = st_transform(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", quiet = TRUE, layer = "unified stands 2016"),
                          make_compound_crs(6557, 8228)) %>% # for DSM v3, keep in sync with same code in treetopDetection.R
                          # st_crs(6557)) %>% # for runs against DSM v3
  select(standID2016)
treetopRandomForest = readRDS("trees/segmentation/treetops/random forest s4268 617k VSURF Pde m9n4.Rds") 
forestTreetopsPath = file.path(treetopOptions$dataPath, "treetops/rf v2")

elliottCompoundCrs = make_compound_crs(6557, 8228)
elliottILandResourceUnitBoundary = st_transform(st_read("iLand/gis/Elliott + Hakki 400 m buffer resource unit snapped.gpkg", quiet = TRUE), crs = elliottCompoundCrs)

#workarounds for https://github.com/r-spatial/stars/issues/777
#slope = st_transform(read_stars(file.path(treetopOptions$dataPath, "bare earth slope Gaussian 10 m EPSG6556.tif")), crs = st_crs(6557))
#write_stars(slope, file.path(treetopOptions$dataPath, "bare earth slope Gaussian 10 m EPSG6557.tif")) # broken in stars 0.7-3, workaround by reprojecting in QGIS
#aspect = st_as_stars(180/pi * atan2(terra::rast(file.path(treetopOptions$dataPath, "bare earth sin(aspect) Gaussian 10 m EPSG6557.tif")),
#                                    terra::rast(file.path(treetopOptions$dataPath, "bare earth cos(aspect) Gaussian 10 m EPSG6557.tif"))))
#write_stars(aspect, file.path(treetopOptions$dataPath, "bare earth aspect Gaussian 10 m EPSG6557.tif"))

slope = read_stars(file.path(treetopOptions$dataPath, "bare earth slope Gaussian 10 m EPSG6557.tif")) # slopes are in degrees
names(slope)[1] = "slope" # default raster band name is filename
aspect = read_stars(file.path(treetopOptions$dataPath, "bare earth aspect Gaussian 10 m EPSG6557.tif")) # aspects are in degrees
names(aspect)[1] = "aspect"
topographicShelter = read_stars(file.path(treetopOptions$dataPath, "topographic shelter Gaussian 10 m EPSG6557.tif")) # shelter index notionally in degrees
names(topographicShelter)[1] = "topographicShelter"

# do treetop prediction
cat(paste0("Checking for treetops on ", length(localMaximaFileNames), " tiles and classifying local maxima if not present...\n"))
treetopStartTime = Sys.time()
with_progress({
  progressBar = progressor(steps = length(localMaximaFileNames))
  
  forestTreetops = bind_rows(future_map(localMaximaFileNames, function(localMaximaFileName) # 46.04m with 12 workers, 9900X
  {
    require(igraph) # dynamically loaded when merge points are obtained
    require(ranger) # ranger apparently doesn't flow to workers for some reason, so predict() fails to resolve without this
    
    tileName = tools::file_path_sans_ext(localMaximaFileName)
    tileForestTreetopsPath = file.path(forestTreetopsPath, paste0(tileName, ".gpkg"))
    if (file.exists(tileForestTreetopsPath))
    {
      cat(paste0("Skipped ", tileName, " as a predicted treetops tile is already present...\n"))
      return(NULL)
    }
  
    # load tile, also neighborhood if density predictors are used
    tileMaxima = get_treetop_eligible_maxima(tileName) %>% filter(is.na(cmmSlope3) == FALSE) #, is.na(ring4mean) == FALSE) # s03870w06600 and s04020w07230, at least, have ring4mean NAs that fail random forest prediction on ring4delta
    tileNeighborhood = get_treetop_eligible_neighborhood(tileName, tileMaxima)
    tileIndices = which(tileNeighborhood$tile == tileName)
  
    # classify local maxima
    tileMaxima$treetop = predict(treetopRandomForest, get_treetop_predictors(tileNeighborhood) %>% filter(tile == tileName), num.threads = treetopOptions$rangerThreads)$predictions
    tileNeighborhood$treetop = factor(NA, levels = levels(tileMaxima$treetop)) # needed by get_merge_points()
    tileNeighborhood$treetop[tileIndices] = tileMaxima$treetop
  
    # cluster treetops and update local maxima classifications with revised cluster membership
    # number of output merge points = initial number of merge points - single treetops + number of merge cluster members added
    tileMergePoints = get_merge_points(tileMaxima, tileNeighborhood)
    tileMaxima$treetop[tileMergePoints$singleTreetops$tileIndex] = "yes"
    tileMaxima$treetop[tileMergePoints$mergePoints$tileIndex] = "merge"
    tileMaxima$treetop[tileMergePoints$ejectedTreetops$tileIndex] = "no"
    tileMaxima$mergeClusterNumber = NA_integer_
    tileMaxima$mergeClusterNumber[tileMergePoints$mergePoints$tileIndex] = tileMergePoints$mergePoints$mergeClusterNumber
    
    # debugging fork
    #if (FALSE)
    #{
    #  tileMaxima %>% filter(id == 175580) %>% select(id, uniqueID, treetop, sourceID, dsmZ, height, elevation, mergeClusterID)
    #  tileMergePoints$mergePoints %>% filter(id == 175580)
    #  tileMergePoints$mergeTreetops %>% filter(treeID == 175580)
    #  tileMergePoints$singleTreetops %>% filter(uniqueID == 366007020000000 + 175580)
    #  
    #  tileMaxima2 = tileMaxima # comment out updates to tileMaxima$treetop above
    #  tileMaxima2$treetop[tileMergePoints$singleTreetops$tileIndex] = "yes"
    #  tileMaxima2$treetop[tileMergePoints$mergePoints$tileIndex] = "merge"
    #  tileMaxima2$treetop[tileMergePoints$ejectedTreetops$tileIndex] = "no"
    #  tileMaxima2$mergeClusterNumber = NA_integer_
    #  tileMaxima2$mergeClusterNumber[tileMergePoints$mergePoints$tileIndex] = tileMergePoints$mergePoints$mergeClusterNumber
    #
    #  tibble(expectedMergePointsOut = nrow(tileMergePoints$mergePoints), uniqueMergeIndices = length(unique(tileMergePoints$mergePoints$tileIndex)), actualMergePointsOut = sum(tileMaxima2$treetop == "merge"), # should all be identical
    #         mergePointsIn = sum(tileMaxima$treetop == "merge"), singleTops = nrow(tileMergePoints$singleTreetops), ejectedTops = nrow(tileMergePoints$ejectedTreetops)) %>% 
    #    mutate(actualMergePointsAdded = actualMergePointsOut - mergePointsIn + singleTops)
    #
    #  st_write(tileMaxima, file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/debug", paste0(tileName, ".gpkg")), layer = "initial maxima classification", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    #  st_write(tileMaxima2, file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/debug", paste0(tileName, ".gpkg")), layer = "merge cluster reformed classification", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    #}
    
    if (sum(tileMaxima$treetop == "merge") != nrow(tileMergePoints$mergePoints))
    {
      stop("Internal consistency failure. Expected merge point clustering and single treetop revisions to result in ", nrow(tileMergePoints$mergePoints), " merge points on tile ", tileName, " but ", sum(tileMaxima$treetop == "merge"), " merge points are defined after clustering.")
      #tileMaxima[setdiff(which(tileMaxima$treetop == "merge"), tileMergePoints$mergePoints$tileIndex), ] %>% select(-starts_with("ring"), -geom)
    }
    
    # assemble treetops and treetop statistics for tile
    # get_treetop_eligible_maxima() standardizes to metric, need to convert back to English units for correct coordinates and attributes on English CRSes.
    # Could also change CRSes to metric but the implementation preference here is is to flow input CRS.
    tileCrs = st_crs(attr(tileMaxima, "crs"))
    tileSingleAndMergeTreetopsPlusNoise = st_join(st_as_sf(bind_rows(tileMaxima %>% filter(treetop != "no") %>% rename(treeID = id) %>% select(-sourceID, -uniqueID, -uniqueMergeClusterID, -starts_with("ring"), -dsmSlope, -cmmSlope3),
                                                                     tileMergePoints$mergeTreetops) %>%
                                                             mutate(x = if(tileCrs$units_gdal == "foot") { 3.28084 * x } else { x },
                                                                    y = if(tileCrs$units_gdal == "foot") { 3.28084 * y } else { y },
                                                                    elevation = if(tileCrs$units_gdal == "foot") { 3.28084 * elevation } else { elevation }, # DTM elevation
                                                                    height = if(tileCrs$units_gdal == "foot") { 3.28084 * height } else { height },
                                                                    radius = if(tileCrs$units_gdal == "foot") { round(3.28084 * radius, 3) } else { radius }, # debatable but, for now, assume whole number preference in surface model cell size
                                                                    dsmZ = if(tileCrs$units_gdal == "foot") { 3.28084 * dsmZ } else { dsmZ }, # debatable if DSM and CMM elevations need to flow, but they're flown for now
                                                                    cmmZ = if(tileCrs$units_gdal == "foot") { 3.28084 * cmmZ } else { cmmZ }),
                                                           coords = c("x", "y", "elevation"), crs = tileCrs, sf_column_name = "geom"), # crs attribute set by get_treetop_eligible_maxima()
                                                  stands2016, left = TRUE)
    tileTreetopsToWrite = tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop %in% c("yes", "merge treetop")) %>% 
      select(-treetop, -mergeClusterID, -mergeClusterNumber) %>% # no merge cluster ID update needed as it's the same as treeID
      mutate(mergePoints = replace_na(mergePoints, as.integer(1)), mergePointsOnTile = replace_na(mergePointsOnTile, as.integer(1))) %>%
      relocate(tile, treeID, height, radius, cmmZ, dsmZ, mergePoints, mergePointsOnTile, standID2016)
    tileTreetopsToWrite = st_filter(tileTreetopsToWrite, elliottILandResourceUnitBoundary) # clip to project area of interest
    if(nrow(tileTreetopsToWrite) == 0)
    {
      # tile is entirely outside the area of interest or has no trees within area of interest
      # Applies to Umpqua estuary tiles s03840w07290, s03870w07290, s03960w07260, and s04080w07170 at surrounding distance 1.
      cat(paste0(tileName, " does not have any treetops...\n")) # but it may have noise points and subsequent processing may expect .gpkg with an empty treetops layer to be present
      
      tileByHeightClass = tibble()
    } else {
      cat(paste0(tileName, "...\n"))
      
      # st_extract() is relatively slow on full rasters (~650 ms) compared to crop and extract (~50 ms)
      # Minor optimization but it got coded, so no reason to back it out and take the performance penalty.
      tileExtentInFeet = get_tile_extent_6557(tileName)
      tileBoundingBox = st_bbox(c(xmin = tileExtentInFeet$xMin, ymin = tileExtentInFeet$yMin, xmax = tileExtentInFeet$xMax, ymax = tileExtentInFeet$yMax), crs = st_crs(6557))
      tileTreetops2D = st_as_sf(as_tibble(st_coordinates(tileTreetopsToWrite)[, c("X", "Y")]), coords = c("X", "Y"), crs = get_projected_crs(tileCrs$wkt))
      
      slopeCrop = st_crop(slope, tileBoundingBox)
      physiographicExtract = st_extract(slopeCrop, tileTreetops2D)
      tileTreetopsToWrite$slope = physiographicExtract$slope
      
      aspectCrop = st_crop(aspect, tileBoundingBox)
      physiographicExtract = st_extract(aspectCrop, tileTreetops2D)
      tileTreetopsToWrite$aspect = physiographicExtract$aspect
      
      topographicShelterCrop = st_crop(topographicShelter, tileBoundingBox)
      physiographicExtract = st_extract(topographicShelterCrop, tileTreetops2D)
      tileTreetopsToWrite$topographicShelter = physiographicExtract$topographicShelter
      
      tileByHeightClass = get_tile_by_height_class(tileName,
                                                   tileTreetopsToWrite, # 1 m height classes are metric, so get summary before English unit conversion (if applicable)
                                                   tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "merge"), # merge, noise, and maybe noise layers aren't written with heights
                                                   tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "noise"), 
                                                   tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "maybe noise"))
      #tileByHeightClass %>% summarize(treetops = sum(treetops), mergePoints = sum(mergePoints), noisePoints = sum(noisePoints), maybeNoisePoints = sum(maybeNoisePoints))
    }

    # assemble merge, noise, and any maybe noise points
    # Currently merge points associated with treetops outside of the study area boundary are dropped. Noise and maybe noise points are
    # retained to tile edge, which is a bit arbitrary but it's presently unclear if clipping them to the study area boundary would be
    # more helpful or less helpful than not doing so.
    tileMergePointsToWrite = tileSingleAndMergeTreetopsPlusNoise %>% 
      filter(treetop == "merge") %>% 
      rename(id = treeID, clusterID = mergeClusterID) %>%
      group_by(mergeClusterNumber) %>%
      mutate(clusterID = min(id)) %>%
      filter(clusterID %in% tileTreetopsToWrite$treeID) %>%
      ungroup() %>%
      select(id, clusterID) #, mergeClusterNumber)
    tileNoisePointsToWrite = tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "noise") %>% rename(id = treeID) %>% select(id)
    tileMaybeNoisePointsToWrite = tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "maybe noise") %>% select(treeID) %>% rename(id = treeID)
    #tileTreetopsToWrite %>% filter(treeID %in% c(175058, 175309, 175580))
    #tileMergePointsToWrite %>% filter(id %in% c(175058, 175309, 175580))
    #tileSingleAndMergeTreetopsPlusNoise %>% filter(treetop == "merge") %>% select(treeID, mergeClusterID, mergeClusterNumber)
    
    # write tile's treetop GeoPackage
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

# attach physiographic variables to accepted treetops if needed
cat(paste0("Checking ", length(acceptedTreetopFileNames), " accepted treetop tiles for physiographic variables..."))
for(acceptedTreetopFileName in acceptedTreetopFileNames)
{
  tileName = tools::file_path_sans_ext(acceptedTreetopFileName)
  tileForestTreetopsPath = file.path(forestTreetopsPath, paste0(tileName, ".gpkg"))
  if (file.exists(tileForestTreetopsPath))
  {
    cat(paste0("Skipped ", tileName, " as an accepted treetops tile is already present...\n"))
    next
  }
  
  cat(paste0(tileName, "...\n"))
  # copy accepted tile into place as only treetops need to be updated
  # Merge and noise points are passthrough. For now it's assumed accepted treetops all lie within the area of interest, so cropping
  # is not needed.
  acceptedTreetopFileSourcePath = file.path(acceptedTreetopsDsmPath, acceptedTreetopFileName)
  file.copy(acceptedTreetopFileSourcePath, tileForestTreetopsPath)
  
  acceptedTreetops = st_read(tileForestTreetopsPath, layer = "treetops", quiet = TRUE)
  acceptedTreetopCoordinates = st_coordinates(acceptedTreetops)
  if (any(c("slope",  "aspect") %in% names(acceptedTreetops)))
  {
    stop(paste0("Accepted treetops tile ", tileName, " already has slope or aspect attached."))
  }

  # accepted treetops come from QGIS and thus may or may not have 2.5D points with elevations, though they have a compound CRS
  # Elevations thus need to be attached and a 2D version with only a projected CRS is required by st_extract().
  tileExtentInFeet = get_tile_extent_6557(tileName)
  tileBoundingBox = st_bbox(c(xmin = tileExtentInFeet$xMin, ymin = tileExtentInFeet$yMin, xmax = tileExtentInFeet$xMax, ymax = tileExtentInFeet$yMax), crs = st_crs(6557))
  acceptedTreetops2D = st_as_sf(as_tibble(acceptedTreetopCoordinates), coords = c("X", "Y"), crs = get_projected_crs(st_crs(acceptedTreetops)$wkt))

  if (ncol(acceptedTreetopCoordinates) == 2) # should any(is.na(acceptedTreetopCoordinates[, "Z"])) also be checked?
  {
    # s04200w06810, s04230w06840 have 2D geometry with compound CRSes (6557 + 8228)
    dtm = read_stars(file.path(dtmPath, paste0(tileName, ".tif")))
    names(dtm)[1] = "elevation"
    physiographicExtract = st_extract(dtm, acceptedTreetops2D)
    acceptedTreetopCoordinates = cbind(acceptedTreetopCoordinates, Z = physiographicExtract$elevation)
    acceptedTreetops = st_set_geometry(acceptedTreetops, st_geometry(st_as_sf(as_tibble(acceptedTreetopCoordinates), coords = c("X", "Y", "Z"), crs = st_crs(acceptedTreetops))))
    #range(acceptedTreetopCoordinates[, "Z"])
  }
  
  slopeCrop = st_crop(slope, tileBoundingBox)
  physiographicExtract = st_extract(slopeCrop, acceptedTreetops2D)
  acceptedTreetops$slope = physiographicExtract$slope
  
  aspectCrop = st_crop(aspect, tileBoundingBox)
  physiographicExtract = st_extract(aspectCrop, acceptedTreetops2D)
  acceptedTreetops$aspect = physiographicExtract$aspect
  
  topographicShelterCrop = st_crop(topographicShelter, tileBoundingBox)
  physiographicExtract = st_extract(topographicShelterCrop, acceptedTreetops2D)
  acceptedTreetops$topographicShelter = physiographicExtract$topographicShelter

  st_write(acceptedTreetops, tileForestTreetopsPath, layer = "treetops", append = FALSE, quiet = TRUE)
}

warnings()
cat(paste0("treetop clustered random forest classification ran for ", format(Sys.time() - treetopStartTime), " over ", length(localMaximaFileNames), " tiles."))

if (nrow(forestTreetops) > 0)
{
  forestTreetops %>% 
    summarize(tiles = length(unique(tile)), stands = length(unique(standID2016)), maxHeightInM = max(heightClassInM), 
              elliottTreetops1m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000), 0, treetops)),
              elliottTreetops5m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000) | (heightClassInM < 5), 0, treetops)), 
              totalTreetops1m = sum(treetops), totalTreetops5m = sum(if_else(heightClassInM >= 5, treetops, 0)))
  writexl::write_xlsx(forestTreetops, file.path(forestTreetopsPath, "standsByHeightClass.xlsx")) # 7.5 MB
  forestTreetops %>% filter(heightClassInM > 85) %>% arrange(desc(heightClassInM))
}