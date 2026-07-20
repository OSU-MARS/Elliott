library(arrow)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(sf)
library(stringr)
library(terra)
library(tidyr)
library(writexl)

treeOptions = tibble(rebuildDbhModels = FALSE,
                     rebuildTreeList = FALSE,
                     recalcDbh = TRUE,
                     dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County",
                     includeInvestigatory = FALSE)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), panel.border = element_blank(), plot.title = element_text(size = 10)))

# load treetop locations, elevations, and heights, merge other physiological predictor variables
# Stand level variables (top height, relative height, ABA, AAT) are calculated below based on trees' stand IDs.
if (treeOptions$rebuildTreeList)
{
  # manual setup in QGIS 3.34 since both sf and terra are both intractibly slow or fail
  # - sample raster values to attach DTM elevation to treetops
  # - rename the sampled column to elevation (layer properties -> fields -> edit)
  # - persist sampled layer to a GeoPackage to create a spatial index (creating a spatial index with the toolbox function either hangs or is very slow)
  # - select by location within the iLand simulation boundary (11.8 M trees)
  # - export selected in EPSG:6556 to treetops 400 m rf v2 (transitory).gpkg with layer name treetops merged 400 m
  
  # reproject trees and crop
  elliottILandResourceUnitBoundary = st_read("iLand/gis/Elliott + Hakki 400 m buffer resource unit snapped.gpkg", quiet = TRUE)
  mergedTreeReadStart = Sys.time() # 1.4 minute load + 1 minute transform, 9900X @ 2 GB (12.3 million trees)
  elliottTrees = st_transform(st_read(file.path(treeOptions$dataPath, "treetops/treetops merged rf v2 (transitory).gpkg"), quiet = TRUE), crs = st_crs(6556))
  Sys.time() - mergedTreeReadStart # 2.6 minutes, 9950X
  
  #st_write(elliottTrees, file.path(treeOptions$dataPath, "treetops/treetops merged rf v2 6556 (transitory).gpkg"))
  iLandTreeIntersectStart = Sys.time() # >50 minutes since runs single threaded, apparently without spatial indexing, 9900X (crop() in terra 1.7-55 appears computationally intractable)
  elliottTrees = st_intersection(elliottTrees, elliottILandResourceUnitBoundary)
  Sys.time() - iLandTreeIntersectStart
  st_write(elliottTrees, "treetops 400 m rf v2 (transitory).gpkg", layer = "treetops merged 400 m")
  
  # requires 72 GB DDR @ 11.8 M trees
  elliottTrees = st_read(file.path(treeOptions$dataPath, "treetops", "treetops 400 m rf v2 (transitory).gpkg"), layer = "treetops merged 400 m", quiet = TRUE) # ~35 s to load with terra::vect() but z is dropped, so 2.7 min with st_read()
  elliottTrees$elevation = 0.3048 * elliottTrees$elevation # CRS is metric from QGIS export but field values need conversion
  elliottTrees$height = 0.3048 * elliottTrees$height
  elliottTrees$radius = 0.3048 * elliottTrees$radius

  elliottStands2016 = st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", layer = "unified stands 2022 property boundary split", quiet = TRUE)
  standJoinStart = Sys.time() # 1.5 minutes, 9900X
  elliottTrees = st_join(elliottTrees, elliottStands2016 %>% select(standID2016, isBuffer, isExternalBoundarySplit))
  Sys.time() - standJoinStart
  # vector-vector extract() in terra 1.7-55 is computationally intractable
  #elliottStands2016 = vect("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", layer = "unified stands 2016") # EPSG:6556
  #elliottStands2016 = terra::extract(elliottStands2016[, "standID2016"], elliottTrees)
  
  tibble(elevation = range(elliottTrees$elevation, na.rm = TRUE), elevNA = sum(is.na(elliottTrees$elevation)),
         height = range(elliottTrees$height), radius = range(elliottTrees$radius))
  
  # join physiographic predictors
  # DTM.vrt for elevation uses 90 GB of memory before failing with long vectors not supported yet: ../include/Rinlinedfuns.h:537 in stars 0.6-8
  #library(stars)
  #slopeReadStart = Sys.time() # 50 s, 9900X
  #slope = st_transform(read_stars(file.path(treeOptions$dataPath, "bare earth slope Gaussian 10 m EPSG6557.tif")), st_crs(6556))
  #Sys.time() - slopeReadStart
  #
  #slopeExtractStart = Sys.time()
  #elliottTrees = st_extract(slope, elliottTrees) # fails with cannot allocate vector of size 7529715.4 Gb
  #Sys.time() - slopeExtractStart
  
  #elevation = rast(file.path(treeOptions$dataPath, "DTM", "DTM.vrt")) # terra 1-8.29 is intractably slow
  #names(elevation) = "elevation"
  slope = project(rast(file.path(treeOptions$dataPath, "bare earth slope Gaussian 10 m EPSG6557.tif")), crs("epsg:6556"), threads = TRUE) # 8.5 s 5950X
  names(slope) = "slope"
  aspect = project(180/pi * atan2(rast(file.path(treeOptions$dataPath, "bare earth sin(aspect) Gaussian 10 m EPSG6557.tif")), 
                                  rast(file.path(treeOptions$dataPath, "bare earth cos(aspect) Gaussian 10 m EPSG6557.tif"))), 
                   crs("epsg:6556"), 
                   threads = TRUE) # ~10 s to load and calculate aspect, ~12 s to reproject
  aspect = ifel(aspect >= 0, aspect, 360 + aspect) # dplyr::if_else() computationally intractable
  names(aspect) = "aspect"
  topographicShelter = rast(file.path(treeOptions$dataPath, "topgraphic shelter Gaussian 10 m.tif")) # already EPSG:6556
  names(topographicShelter) = "topographicShelterIndex"
    
  physiographicExtractStart = Sys.time() # 5.2 minutes 9900X, up to ~95 GB of DDR depending what R feels like
  elliottTrees = terra::extract(slope, elliottTrees, bind = TRUE)
  elliottTrees = terra::extract(aspect, elliottTrees, bind = TRUE)
  elliottTrees = terra::extract(topographicShelter, elliottTrees, bind = TRUE)
  Sys.time() - physiographicExtractStart
  
  # 2.3 GB on disk @ 11.8 M trees
  writeVector(elliottTrees, file.path(treeOptions$dataPath, "treetops", "treetops 400 m rf v2 predictors (transitory).gpkg"), layer = "treetops", overwrite = TRUE)
}

# read trees even if treeOptions$rebuildTreeList == TRUE to switch from terra to sf
elliottTreeReadStart = Sys.time() # 61 s, ~8 GB in memory
elliottTrees = st_read(file.path(treeOptions$dataPath, "treetops", "treetops 400 m rf v2 predictors (transitory).gpkg"), layer = "treetops", quiet = TRUE)
Sys.time() - elliottTreeReadStart

#elliottTreeReadStart = Sys.time() # 36 s, ~32 GB in memory
#elliottTrees = vect(file.path(treeOptions$dataPath, "treetops", "treetops 400 m rf v2 predictors (transitory).gpkg"), layer = "treetops")
#Sys.time() - elliottTreeReadStart


## assign species and predict DBH
# Hardwood-conifer-snag classification permits mainly differentiation of Douglas-fir and red alder. In 2015-16 Elliott cruise data,
#                 % of stems
# Douglas-fir     76.8
# red alder        9.6
# other hardwood   7.0
# other conifer    6.4
#
# so assigning only Douglas-fir and red alder is at least 86.4% correct, likely higher at whole forest scale as Douglas-fir and red alder 
# tend to be taller than other most conifers (western hemlock, western redcedar) and hardwoods (bigleaf maple, Oregon myrtle, cascara
# buckthorn, Pacific madrone) and are thus more likely to be the overstory trees most detectable in fixed-wing LiDAR. Implied species 
# accuracy is potentially >92% for conifers and 58% for hardwood.
#
# Parallel evaluation not viable due to chronic furrr 0.3.1 future_map() failures of the form MultisessionFuture (<none>) failed to call grmall() on cluster RichSOCKnode #1 (PID 21712 on localhost ‘localhost’). The reason reported was ‘error writing to connection’. Post-mortem diagnostic: No process exists with this PID, i.e. the localhost worker is no longer alive.
# Circumstantial evidence suggests worker processes may be exiting due to lack of thread safety in mgcv::predict.gam().
if (treeOptions$recalcDbh)
{
  ## most preferred models for initial DBH prediction: can't use ABA or AAT as DBH hasn't yet been predicted
  # height-diameter/setup.R + preferred species models saved from PSME.R and ALRU2.R
  stands2022 = as_tibble(st_drop_geometry(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", layer = "unified stands 2022 property boundary split", quiet = TRUE))) %>%
    # Elliott State Forest (ESF) stands beyond the Elliott State Research Forest (ESRF) boundary lack inventory data and thus have NA for isPlantation
    # DBHes predicted by models factorized on plantation status will thus be NA, so NA values have to be resolved somehow
    # 2025-05-07: ESF stands were reviewed in 2021 orthoimagery where available and plantations marked where obvious, making non-plantation the
    #  most accurate simple assumption even though it is likely to include some older plantations
    mutate(isPlantation = replace_na(isPlantation, 0))

  load("trees/height-diameter/data/ALRU2 preferred models.Rdata")
  load("trees/height-diameter/data/PSME preferred models.Rdata")
  rm(alruHeightFromDiameterPreferred, psmeHeightFromDiameterPreferred)
  
  # check snag classification threshold
  snagThreshold = 0.25
  snagFraction = (elliottTrees$BrownTree + elliottTrees$GreyTree) / (elliottTrees$Unclassified + elliottTrees$Bare + elliottTrees$BareShadow + elliottTrees$BrownTree + elliottTrees$GreyTree + elliottTrees$Conifer + elliottTrees$ConiferShadow + elliottTrees$ConiferDeepShadow + elliottTrees$Hardwood + elliottTrees$HardwoodShadow + elliottTrees$HardwoodDeepShadow)
  tibble(snagThreshold = snagThreshold, liveTreePercentage = 100 * sum(snagFraction < snagThreshold) / length(snagFraction))
  
  # predict diameters
  # Dataset limitations mean predictions are either as Douglas-fir or as red alder, see notes below.
  startTime = Sys.time() # 40 s 9900X @ 11.8 M trees (potentially ~3 minutes if R is slow), 5950X 25 s dplyr + prediction time (~24 seconds for three nonlinear iterations) @ 12.3 M trees
  elliottTreesMod = left_join(elliottTrees,
                              stands2022 %>% group_by(standID2016) %>% summarize(standArea = sum(standArea), isPlantation = any(isPlantation == 1)), # recombine portions of stands with isExternalBoundarySplit = 1
                              by = join_by("standID2016")) %>% # ~1 s for join
    # for now, simple conifer-hardwood separation
    mutate(classification = factor(if_else((BrownTree + GreyTree) / (Unclassified + Bare + BareShadow + BrownTree + GreyTree + Conifer + ConiferShadow + ConiferDeepShadow + Hardwood + HardwoodShadow + HardwoodDeepShadow) < snagThreshold,
                                           if_else((Conifer + 0.7 * ConiferShadow + 0.5 * ConiferDeepShadow) > (Hardwood + 0.7 * HardwoodShadow + 0.5 * HardwoodDeepShadow), "conifer", if_else(height < 53, "hardwood", "conifer")), "snag"), # block hardwood classification above 53 m as no hardwoods taller than 50 m are present in 2015-16 Elliott cruise data
                                     levels = c("hardwood", "conifer", "snag")),
           # for now, assume all trees lacking stand IDs are on plantations outside the current boundary of stand ID polygons
           # TODO: this is probably wrong and needs to be rechecked after outer boundary of stand polygons is snapped to the iLand simulation area
           standArea = if_else(is.na(standID2016), 41619.9 - 41095.8, standArea), # standArea will be NA if standID is NA
           isPlantation = if_else(is.na(standID2016), 1, isPlantation), # isPlantation will be NA if standID is NA
           standID2016 = replace_na(standID2016, 0)) %>% # set default stand ID
    rename(TotalHt = height) %>% # change height to TotalHt to integrate with DBH model fits
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% # put tallest detected trees first in each stand for top height calculation: currently no detection of snags or broken tops
    mutate(treeID = 1E6 * standID2016 + row_number(), # generate unique IDs within 2016 stands: permits 999,999 trees per stand (max segmented is about 432,000 after filtering) as stand ID is four positive digits (32 bit unsigned int max is 4294 967 296 => largest stand ID in .feather is 4294)
           measureTreeTphContribution = 1 / standArea, # since top height trees are the tallest in the stand assume all of them are segmented from LiDAR: expansion factor for top height calculation is thus 1/(delineated area of stand in hectares)
           topHeightTph = pmin(cumsum(if_else(is.na(TotalHt), 0, measureTreeTphContribution)), 100), # TPH total towards the H100 definition of top height, trees not measured for TotalHt are skipped
           topHeightWeight = pmax((topHeightTph - lag(topHeightTph, default = 0)) / measureTreeTphContribution, 0), # clamp remaining fraction to [0, 1] to get individual trees' contributions to the top height average
           topHeight = sum(topHeightWeight * TotalHt, na.rm = TRUE) / sum(topHeightWeight, na.rm = TRUE), # m, tallest 100 trees per hectare
           relativeHeight = TotalHt / topHeight) %>%
    select(-measureTreeTphContribution) %>%
    group_by(classification) %>%
    mutate(dbhBootstrap = case_when(cur_group()$classification %in% c("conifer", "snag") ~ mgcv::predict.gam(psmeDiameterFromHeightPreferred$gam, pick(everything())), # as no better option is available, assume all live conifers and all snags visible in aerial imagery and aerial LiDAR are Douglas-fir
                                    cur_group()$classification == "hardwood" ~ predict(alruDiameterFromHeightPreferred$ruark, pick(everything()))), # as no better option is available, assume all hardwoods are red alder
           # initial estimate of tree's basal area in m²
           # For now, treat basal area NAs as zero when totaling stand basal area. This is correct for trees shorter than breast height 
           # (where DBH is undefined) and prevents NAs from trees below DBH or with DBH prediction issues from propagating to entire stands.
           # 2025-05-07 dataset: 5 conifers + 103 broadleaves < breast height -> 1.6 M NA DBHes without BA zeroing
           basalArea = if_else(TotalHt >= 1.37, replace_na(pi/4 * (0.01 * dbhBootstrap)^2, 0), 0)) %>%
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% 
    # initial estimate basal area of stand in m²/ha with clamp to plausibility bound as bootstrap DBH estimates can be high, TODO: refine adjustment for undetected trees
    mutate(standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
           basalAreaAdjustmentFactor = if_else(standBasalAreaApprox <= 125, 1, 125 / standBasalAreaApprox),
           standBasalAreaApprox = basalAreaAdjustmentFactor * standBasalAreaApprox,
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(lag(basalArea, default = 0)) / standArea) %>%  # basal area of taller trees in m²/ha
    group_by(classification) %>%
    # avoid generalized GAMs due to tendency to out of range predictions
    mutate(dbh = case_when(cur_group()$classification %in% c("conifer", "snag") ~ predict(psmeDiameterFromHeightPreferred$ruarkAbatPhysio, pick(everything())),
                           cur_group()$classification == "hardwood" ~ predict(alruDiameterFromHeightPreferred$ruarkAbatPhysio, pick(everything())))) %>%
           # basal areas updated below
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% 
    mutate(standBasalAreaBootstrap = standBasalAreaApprox,
           tallerApproxBasalAreaBootstrap = tallerApproxBasalArea,
           basalArea = if_else(TotalHt >= 1.37, replace_na(pi/4 * (0.01 * dbh)^2, 0), 0),
           standBasalAreaApprox = 1 / 1 * sum(basalArea, na.rm = TRUE) / standArea,
           basalAreaAdjustmentFactor = if_else(standBasalAreaApprox <= 125, 1, 125 / standBasalAreaApprox),
           standBasalAreaApprox = basalAreaAdjustmentFactor * standBasalAreaApprox,
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(lag(basalArea, default = 0)) / standArea) %>%
    # further iteration results in small changes at most percentiles but drives the smallest <0.5% of trees to negative DBH
    #group_by(classification) %>%
    #mutate(dbh = case_when(cur_group()$classification %in% c("conifer", "snag") ~ predict(psmeDiameterFromHeightPreferred$, pick(everything())),
    #                       cur_group()$classification == "hardwood" ~ predict(alruDiameterFromHeightPreferred$, pick(everything()))) %>%
    #group_by(standID2016) %>%
    #mutate(standBasalAreaInitial = standBasalAreaApprox,
    #       tallerApproxBasalAreaInitial = tallerApproxBasalArea,
    #       if_else(TotalHt >= 1.37, replace_na(basalArea = pi/4 * (0.01 * dbh)^2, 0), 0),
    #       standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
    #       tallerApproxBasalArea = cumsum(lag(basalArea, default = 0)) / standArea) %>%
    #group_by(classification) %>%
    #mutate(dbh2 = case_when(cur_group()$classification %in% c("conifer", "snag") ~ predict(psmeDiameterFromHeightPreferred$, pick(everything())),
    #                        cur_group()$classification == "hardwood" ~ predict(alruDiameterFromHeightPreferred$, pick(everything()))) %>%
    #group_by(standID2016) %>%
    #mutate(standBasalArea2 = standBasalAreaApprox,
    #       tallerApproxBasalArea2 = tallerApproxBasalArea,
    #       if_else(TotalHt >= 1.37, replace_na(basalArea = pi/4 * (0.01 * dbh2)^2, 0), 0),
    #       standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
    #       tallerApproxBasalArea = cumsum(lag(basalArea, default = 0)) / standArea) %>%
    #group_by(classification) %>%
    #mutate(dbh = case_when(cur_group()$classification %in% c("conifer", "snag") ~ predict(psmeRuarkAbatPhysio, pick(everything())),
    #                       cur_group()$classification == "hardwood" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
    #                       cur_group()$classification == "TSHE" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    ungroup() %>%
    rename(height = TotalHt)
  Sys.time() - startTime

  # check DBH imputation: all trees taller than breast height should have DBH > 0
  # For now, trees less than breast height are retained as it's unclear if they should be filtered out.
  st_drop_geometry(elliottTreesMod) %>% group_by(classification) %>% summarize(trees = n(), inputNA = sum(is.na(height) | is.na(isPlantation) | is.na(slope) | is.na(aspect)), 
                                                                               subBreastHt = sum(height < 1.37), bootstrapDbhNA = sum(is.na(dbhBootstrap)), bootstrapBAna = sum(is.na(standBasalAreaBootstrap) | is.na(tallerApproxBasalAreaBootstrap)), dbhNegativeOrZero = sum(dbh <= 0, na.rm = TRUE), dbhNA = sum(is.na(dbh)), isOversize = sum(dbh > 300, na.rm = TRUE), pctValid = 100 * (trees - dbhNegativeOrZero - dbhNA - isOversize) / trees)
  st_drop_geometry(elliottTreesMod) %>% filter(standID2016 == 0) %>% summarize(treesInDefaultStand = n()) # only 1
  # setdiff(unique(elliottTrees$standID2016), unique(stands2022$standID2016)) # missing stand information is a common cause of NAs

  #saveRDS(elliottTreesMod, file = file.path(treeOptions$dataPath, "treetops", "trees rf v2.Rds")) # ~2 minutes, writes 1.1 GB
  startTime = Sys.time()
  st_write(elliottTreesMod %>% select(tile, standID2016, treeID, classification, height, dbh), dsn = file.path(treeOptions$dataPath, "treetops", "trees and snags rf v2.gpkg"), layer = "trees and snags 2021 rf v2") # minutes, writes 1.7 GB
  Sys.time() - startTime
} else {
  #elliottTreesMod = readRDS(file.path(treeOptions$dataPath, "treetops", "trees rf v2.Rds"))
  elliottTreesMod = st_read(file.path(treeOptions$dataPath, "treetops", "trees and snags rf v2+v1.gpkg"), layer = "trees and snags 2021 rf v2", quiet = TRUE)
}


## write trees for iLand
elliottTreesMod %>% summarize(naStand = sum(is.na(standID2016)), naTree = sum(is.na(treeID)), naSpecies = sum(is.na(species)), naHeight = sum(is.na(height)), naDbh = sum((height >= 1.37) & is.na(dbh))) # should all be zero for usable iLand tree list
elliottTreesMod %>% summarize(ruXbelow = sum(x <= 106000), ruXabove = sum(x >= 131700), ruYbelow = sum(y <= 194100), ryYabove = sum(y >= 222200)) # check against resource unit grid bounds (Elliott.xml <resourceUnitFile>), must be zero for usable iLand tree list

startTime = Sys.time()
elliottTreesArrow = arrow_table(st_drop_geometry(elliottTreesMod %>% mutate(x = st_coordinates(elliottTrees)[, "X"], # ~17 s conversion to Arrow
                                                                            y = st_coordinates(elliottTrees)[, "Y"])) %>%
                                  # for now, simple live tree-snag separation: iLand trees are only live trees, snags go to carbon and nitrogen pools
                                  # for now, trees less than iLand's shrub layer depth (4 m by default) are not converted to saplings
                                  # could remove remove ~2 M trees below iLand's definition of tree height as 4.0 m (TODO: translate these rows to iLand saplings)
                                  filter(height >= 1.37, classification != "snag") %>%
                                  filter(y > 194100) %>% # temporary workaround for <20 m spillover at southernmost edge of resource unit grid
                                  mutate(fiaCode = case_match(as.character(classification), "conifer" ~ 202, "hardwood" ~ 351), # for lack of a better option, species dub all conifers as Douglas-fir and all hardwoods as red alder ("TSHE" ~ 263)
                                         resourceUnitX = as.integer(x / 100), # dropped in final select
                                         resourceUnitY = as.integer(y / 100)) %>% # case_match() breaks on factors as of dplyr 1.1.4 (2023-12)
                                  arrange(resourceUnitY, resourceUnitX, classification, y, x) %>%
                                  rename(standID = standID2016) %>%
                                  select(standID, treeID, fiaCode, dbh, height, x, y),
                                schema = schema(standID = uint32(), treeID = uint32(), fiaCode = uint16(),
                                                dbh = float32(), height = float32(), x = float32(), y = float32()))
Sys.time() - startTime
write_feather(elliottTreesArrow, "iLand/init/ESRF trees 2025-05-06.feather", compression = "uncompressed") # 448 MB, leave uncompressed for considerably faster iLand startup

if (treeOptions$includeInvestiatory)
{
  elliottStandsMod = left_join(elliottTreesMod %>% group_by(standID2016) %>% 
                                 summarize(segmentedTrees = n(),
                                           segmentedTph = n() / standArea[1],
                                           segmentedTopHeight = topHeight[1],
                                           segmentedQmd = sqrt(standArea[1] * standBasalAreaApprox[1] / (pi/4 * 0.01^2 * n())), # basal area is in m²/ha so need to multiply by stand area since n() counts all trees in the stand
                                           standBasalAreaBootstrap = standBasalAreaBootstrap[1], 
                                           #standBasalAreaInitial = standBasalAreaInitial[1],
                                           standBasalAreaApprox = standBasalAreaApprox[1],
                                           topHeightTrees = sum(topHeightWeight > 0)),
                               stands2022,
                               by = "standID2016") %>%
    mutate(topHeightRatio = segmentedTopHeight / topHeight,
           standClass = factor(isPlantation + (standAge2016 <= 15), labels = c("natural regen", "pre-2001 plantation", "2001+ plantation"), levels = c(0, 1, 2)))
  
  segmentationByHeight = left_join(elliottTreesMod %>% mutate(relativeHeightClass = 0.05 * floor(relativeHeight / 0.05) + 0.5 * 0.05) %>%
                                     group_by(standID2016, relativeHeightClass) %>% 
                                     reframe(standAge2016 = standAge2016[1],
                                             topHeight = topHeight[1],
                                             standArea = standArea[1],
                                             isPlantation = isPlantation[1],
                                             segmentedTph = n() / standArea[1]),
                                   trees2016 %>% rename(standID2016 = StandID) %>%
                                     filter(is.na(measureTreeTphContribution) == FALSE) %>%
                                     mutate(relativeHeight = imputedHeight / topHeight, # no effect for height measure trees, creates relative height for DBH only measure trees to impute TPH by relative height class
                                            relativeHeightClass = 0.05 * floor(relativeHeight / 0.05) + 0.5 * 0.05) %>%
                                     group_by(standID2016, relativeHeightClass) %>% 
                                     reframe(groundTph = sum(meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution) / measurePlotsInStand),
                                   by = c("standID2016", "relativeHeightClass")) %>%
    mutate(groundTph = replace_na(groundTph, 0),
           segmentationPct = 100 * segmentedTph / groundTph)
  
  ggplot(segmentationByHeight %>% group_by(standID2016) %>%
           summarize(standArea = standArea[1], isPlantation = isPlantation[1], topHeight = topHeight[1], segmentedTph = sum(segmentedTph), groundTph = sum(groundTph)) %>%
           filter(groundTph != 0)) +
    geom_point(aes(x = topHeight, y = 100 * segmentedTph / groundTph, color = isPlantation, size = standArea), alpha = 0.3, shape = 16) +
    coord_cartesian(ylim = c(0, 150)) + # excludes ~6 outlying stands of 200-3000%
    labs(x = "top height, m", y = "trees segmented, %", color = "plantation", size = "stand area, ha") +
    theme(legend.spacing.y = unit(0.3, "line"))
  
  ggplot(segmentationByHeight %>% group_by(isPlantation, relativeHeightClass) %>%
           summarize(segmentedTph = sum(standArea * segmentedTph) / sum(standArea),
                     groundTph = sum(standArea * groundTph) / sum(standArea),
                     .groups = "drop") %>%
           group_by(relativeHeightClass) %>%
           mutate(estimatedTphFraction = segmentedTph / sum(segmentedTph)) %>%
           ungroup() %>%
           filter(groundTph != 0)) +
    geom_bar(aes(x = 100 * estimatedTphFraction * segmentedTph / groundTph, y = relativeHeightClass, fill = isPlantation), orientation = "y", stat = "identity") +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 2.5)) + # substantial quantization due to limited ground trees by relative height = 2
    labs(x = "segmented trees, % of ground estimate", y = "relative height", fill = "plantation")
  
  print(segmentationByHeight %>% group_by(isPlantation, relativeHeightClass) %>%
    summarize(segmentedTph = sum(standArea * segmentedTph) / sum(standArea),
              groundTph = sum(standArea * groundTph) / sum(standArea),
              .groups = "drop") %>%
    group_by(relativeHeightClass) %>%
    mutate(estimatedTphFraction = segmentedTph / sum(segmentedTph)), n = 300)
  
  ggplot(elliottStandsMod) + 
    geom_histogram(aes(x = standArea, fill = isPlantation), binwidth = 5) +
    labs(x = "stand area, h", y = "stands", fill = "plantation") +
    theme(legend.spacing.y = unit(0.3, "line"))
}

# check summaries and check plots for dataset alignment and DBH prediction
if (treetopOptions$includeInvestigatory)
{
  #elliottTrees %>% filter(is.na(standID2016)) %>% group_by(STD_ID) %>% summarize(n = n())
  #elliottTreesMod %>% filter(is.na(dbhBootstrap))
  #print(elliottTreesMod %>% select(standID2016, standArea, species, height, topHeightTph, topHeightWeight, topHeight), n = 750)
  #elliottStandsMod %>% filter(standID2016 == 12) %>% select(standID2016, tph, segmentedTph, topHeight, segmentedTopHeight)
  tibble(tph = cor(drop_na(elliottStandsMod %>% select(tph, segmentedTph)))[2,1], topHeight = cor(drop_na(elliottStandsMod %>% select(topHeight, segmentedQmd)))[2,1], ba = cor(drop_na(elliottStandsMod %>% select(standBasalAreaPerHectare, standBasalAreaApprox.y)))[2,1], qmd = cor(drop_na(elliottStandsMod %>% select(qmd, segmentedQmd)))[2,1])
  elliottTreesMod %>% group_by(classification) %>% summarize(trees = n(), minDbh = min(dbh), maxDbh = max(dbh), minHt = min(height), maxHt = max(height), maxRelHt = max(relativeHeight), minRelHt = min(relativeHeight), 
                                                             naDbh = sum(is.na(dbh)), underDbh = sum(dbh < 0.3), naHt = sum(is.na(height)), naTopHt = sum(is.na(topHeight)), underHt = sum(height < 1.37), overAba = sum(tallerApproxBasalArea > standBasalAreaApprox), naRelHt = sum(is.na(relativeHeight)), naX = sum(is.na(x)), naY = sum(is.na(y)))
  speciesLimits = get_species_limits(elliottTreesMod %>% rename(TotalHt = height, DBH = dbh) %>% mutate(speciesGroup = factor(if_else(species == "HW", "RA", species), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other")))) # get_species_limits() in height-diameter/setup.R
  elliottTreesMod %>% mutate(overDbh = dbh > speciesLimits$dbhMax, underTaper = (height / (0.01 * dbh)) < speciesLimits$heightDiameterRatioMin, overTaper = (height / (0.01 * dbh)) > speciesLimits$heightDiameterRatioMax, outOfRangeHt = height > speciesLimits$heightMax) %>% group_by(species) %>% summarize(n = n(), outOfRangeHt = sum(outOfRangeHt), overDbh = sum(overDbh), underTaper = sum(underTaper), overTaper = sum(overTaper))
  elliottTreesMod %>% reframe(quantiles = c(0, 0.005, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.995, 1), height = quantile(height, probs = quantiles), dbhBootstrap = quantile(dbhBootstrap, probs = quantiles), dbh = quantile(dbh, probs = quantiles))
  
  #elliottTreesMod %>% select(standID, treeID, species, height, dbhBootstrap, dbhInitial, dbh, basalArea, standBasalAreaApprox, tallerApproxBasalArea, standBasalAreaBootstrap, tallerApproxBasalAreaBootstrap)
  #print(elliottTreesMod %>% filter(dbhInitial < 0.03) %>% select(standID, treeID, species, height, dbhBootstrap, dbhInitial, standBasalAreaBootstrap, tallerApproxBasalAreaBootstrap, slope, relativeHeight), n = 350)
  
  # trees' height-diameter
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 75, yend = 75), color = "grey70", linetype = "longdash") +
    geom_segment(aes(x = 100, y = 50, xend = 190, yend = 95), color = "grey70", linetype = "longdash") +
    geom_bin2d(aes(x = dbhBootstrap, y = height), binwidth = c(2.5, 1), elliottTreesMod) +
    annotate("text", x = 60, y = 76, label = "H:D = 100", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
    annotate("text", x = 190, y = 96, label = "H:D = 50", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
    coord_cartesian(xlim = c(0, 300), ylim = c(0, 96)) +
    labs(x = "bootstrap DBH, cm", y = "segmented height, m", fill = "trees") +
    scale_fill_viridis_c(breaks = c(1, 100, 10000, 500000), labels = scales::label_comma(), limits = c(1, 500E3), trans = "log10") +
  #ggplot() +
  #  geom_segment(aes(x = 0, y = 0, xend = 75, yend = 75), color = "grey70", linetype = "longdash") +
  #  geom_segment(aes(x = 100, y = 50, xend = 190, yend = 95), color = "grey70", linetype = "longdash") +
  #  geom_bin2d(aes(x = dbhInitial, y = height), binwidth = c(2.5, 1), elliottTreesMod) +
  #  annotate("text", x = 60, y = 76, label = "H:D = 100", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
  #  annotate("text", x = 190, y = 96, label = "H:D = 50", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
  #  coord_cartesian(xlim = c(0, 300), ylim = c(0, 96)) +
  #  labs(x = "initial DBH, cm", y = NULL, fill = "trees") +
  #  scale_fill_viridis_c(breaks = c(1, 100, 10000, 500000), labels = scales::label_comma(), limits = c(1, 500E3), trans = "log10") +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 75, yend = 75), color = "grey70", linetype = "longdash") +
    geom_segment(aes(x = 100, y = 50, xend = 190, yend = 95), color = "grey70", linetype = "longdash") +
    geom_bin2d(aes(x = dbh, y = height), binwidth = c(2.5, 1), elliottTreesMod) +
    annotate("text", x = 60, y = 76, label = "H:D = 100", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
    annotate("text", x = 190, y = 96, label = "H:D = 50", color = "grey70", hjust = 0.5, size = 3, vjust = 0) +
    coord_cartesian(xlim = c(0, 300), ylim = c(0, 96)) +
    labs(x = "predicted DBH, cm", y = NULL, fill = "trees") +
    scale_fill_viridis_c(breaks = c(1, 100, 10000, 500000), labels = scales::label_comma(), limits = c(1, 500E3), trans = "log10") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 1, guides = "collect") &
    theme(legend.spacing.y = unit(0.3, "line"))
  
  # comparison of stand-level properties
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 3000, yend = 3000), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = tph, y = segmentedTph, color = standAge2016, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    coord_equal() +
    labs(x = "2015–16 trees per hectare", y = "2021 trees per hectare", color = "2016 age, years", shape = NULL, size = NULL) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 80, yend = 80), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = topHeight, y = segmentedTopHeight, color = standAge2016, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    coord_equal() +
    labs(x = bquote("2015–16 H"[100]*", m"), y = bquote("2021 H"[100]*", m"), color = "2016 age, years", shape = NULL, size = NULL) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 125, yend = 125), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = standBasalAreaPerHectare, y = standBasalAreaApprox.y, color = standAge2016, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    coord_equal() +
    labs(x = bquote("2015–16 basal area, m"^2*" ha"^-1), y = bquote("2021 basal area, m"^2*" ha"^-1), color = "2016 age, years", shape = NULL, size = NULL) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = qmd, y = segmentedQmd, color = standAge2016, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    coord_equal() +
    labs(x = "2015–16 QMD, cm", y = "2021 QMD, cm", color = "2016 age, years", shape = NULL, size = NULL) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(widths = c(1, 1), heights = c(1, 1), guides = "collect") &
    theme(legend.spacing.y = unit(0.3, "line")) &
    guides(color = guide_colorbar(order = 1), shape = guide_legend(order = 2, override.aes = list(alpha = 1)), size = guide_legend(order = 2)) &
    scale_color_viridis_c(end = 0.96) &
    scale_shape_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(16, 17, 18)) &
    scale_size_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(1.5, 1.5, 1.85))
  #ggsave("Presentation/2015-16 ground to 2021 LiDAR stand comparison.png", units = "cm", width = 20, height = 14, dpi = 150)
  
  # distribution of fraction of trees segmented  
  ggplot() +
    geom_histogram(aes(x = 100 * segmentedTph / tph, y = after_stat(100 * count / sum(count)), fill = isPlantation, weight = standArea), elliottStandsMod, binwidth = 5, na.rm = TRUE) +
    geom_line(aes(x = seq(0, 300), y = 5 * dgamma(0.01 * seq(0, 300), shape = 3.2, rate = 5)), color = "grey30") +
    coord_cartesian(xlim = c(0, 200)) +
    labs(x = "fraction of trees segmented, %", y = "fraction of stands, %", fill = "plantation") +
    theme(legend.spacing.y = unit(0.3, "line"))
  # relative height distribution of 2016 measure trees and 2021 segmented trees
  # total number of measure plots: stands2022 %>% group_by(isPlantation) %>% summarize(measurePlots = sum(measurePlotsInStand, na.rm = TRUE)) = 10,076 plots, 4192 natural regen, 5844 plantation
  #                                trees2016 %>% filter(is.na(DBH) == FALSE) %>% group_by(isPlantation) %>% summarize(measurePlots = n_distinct(PlotID))
  # total cruised area: sum((stands2022 %>% filter(is.na(standBasalAreaPerHectare) == FALSE))$standArea) = 15,981 ha, 7172 ha natural regen + 8809 ha plantation = 45% + 55%
  # total forest area: 33,403 ha, 16468 natural regen + 16935 ha plantation = 49% + 51%
  ggplot(elliottTreesMod %>% filter(isPlantation == FALSE)) +
    geom_histogram(aes(x = relativeHeight, fill = species, weight = 1 / sum((stands2022 %>% filter(isPlantation == FALSE))$standArea)), binwidth = 0.05) + # average across total segmented area: 33,403 ha
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 60)) +
    labs(x = NULL, y = "trees per hectare", title = "a) 2021 natural regen segmentation", fill = "LiDAR\nsegmentation") +
  ggplot(trees2016 %>% filter(isPlantation == FALSE, is.na(DBH) == FALSE) %>% mutate(relativeHeight = imputedHeight / topHeight)) +
    geom_histogram(aes(x = relativeHeight, fill = speciesGroup, weight = meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / n_distinct(PlotID)), binwidth = 0.05) + 
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 60)) +
    labs(x = NULL, y = NULL, title = "b) 2015–16 natural regen ground", fill = "ground") +
  ggplot(elliottTreesMod %>% filter(isPlantation)) +
    geom_histogram(aes(x = relativeHeight, fill = species, weight = 1 / sum((stands2022 %>% filter(isPlantation))$standArea)), binwidth = 0.05) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 60)) +
    labs(x = "relative height", y = "trees per hectare", title = "c) 2021 plantation segmentation", fill = "LiDAR\nsegmentation") +
  ggplot(trees2016 %>% filter(isPlantation, is.na(DBH) == FALSE) %>% mutate(relativeHeight = imputedHeight / topHeight)) +
    geom_histogram(aes(x = relativeHeight, fill = speciesGroup, weight = meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / n_distinct(PlotID)), binwidth = 0.05) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0, 60)) +
    labs(x = "relative height", y = NULL, title = "d) 2015–16 plantation ground", fill = "ground") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_fill_manual(breaks = c("DF", "HW", "RA", "WH", "BM", "OM", "RC", "other"), values = c("forestgreen", "red2", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65"))&
    theme(legend.spacing.y = unit(0.3, "line"))
  #ggsave("presentation/LiDAR-ground TPH distribution.png", height = 12, width = 16, units = "cm", dpi = 150)
    
  # comparison of top height by age
  library(gslnls)
  yearsToBreastHeight = 5
  kingSiteIndexModelGround = gsl_nls(topHeight ~ (standAge2016 - yearsToBreastHeight)^2 / (a + b * (standAge2016 - yearsToBreastHeight) + c * (standAge2016 - yearsToBreastHeight)^2), elliottStandsMod %>% filter(standAge2016 > 5, topHeight < 1.5 * standAge2016), start = list(a = 1, b = 1, c = 1))
  kingSiteIndexModelLidar = gsl_nls(segmentedTopHeight ~ (standAge2016 - yearsToBreastHeight)^2 / (a + b * (standAge2016 - yearsToBreastHeight) + c * (standAge2016 - yearsToBreastHeight)^2), elliottStandsMod %>% filter(standAge2016 > 5, topHeight < (10 + 1.1 * standAge2016)), start = list(a = 1, b = 1, c = 1))
  
  elliottStandsGisSiteIndex = left_join(elliottStandsMod,
                                        read_xlsx("GIS/Planning/Elliott Stand Data Feb2022.xlsx") %>% rename(standID2016 = StandID),
                                        by = "standID2016")
  
  ggplot() +
    geom_line(aes(x = seq(yearsToBreastHeight, 250), y = predict(kingSiteIndexModelGround, tibble(standAge2016 = seq(yearsToBreastHeight, 250))), linetype = "King's"), color = "grey30", linewidth = 0.5) +
    geom_point(aes(x = standAge2016,  y = topHeight, color = isPlantation, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    labs(x = NULL, y = bquote("ground measured H"[100]*" in 2016, m"), color = NULL, fill = NULL, linetype = NULL, shape = NULL, size = NULL) +
  ggplot() +
    geom_line(aes(x = seq(yearsToBreastHeight, 250), y = predict(kingSiteIndexModelLidar, tibble(standAge2016 = seq(yearsToBreastHeight, 250))), linetype = "King's"), color = "grey30", linewidth = 0.5) +
    geom_point(aes(x = standAge2016 + 5,  y = segmentedTopHeight, color = isPlantation, shape = standClass, size = standClass), elliottStandsMod, alpha = 0.3, na.rm = TRUE) +
    labs(x = NULL, y = bquote("LiDAR measured H"[100]*" in 2021, m"), color = NULL, fill = NULL, linetype = NULL, shape = NULL, size = NULL) +
  ggplot() +
    geom_smooth(aes(x = standAge2016, y = 0.3048 * Cruised_Si, linetype = "GAM"), elliottStandsGisSiteIndex %>% filter(Cruised_Si != 0), formula = y ~ s(x), method = "gam", alpha = 0.1, color = "grey30", linewidth = 0.5) +
    geom_point(aes(x = standAge2016, y = 0.3048 * Cruised_Si, color = isPlantation, shape = standClass, size = standClass), elliottStandsGisSiteIndex %>% filter(Cruised_Si != 0), alpha = 0.3) +
    labs(x = "stand age, years", y = "ground measured site index in 2016, m", color = NULL, fill = NULL, linetype = NULL, shape = NULL, size = NULL) +
  ggplot() +
    geom_smooth(aes(x = standAge2016, y = 0.3048 * if_else(ODSL_Site_ == 0, ODSL_Physi, pmin(ODSL_Physi, ODSL_Site_)), linetype = "GAM"), elliottStandsGisSiteIndex %>% filter(Cruised_Si == 0), formula = y ~ s(x), method = "gam", alpha = 0.1, color = "grey30", linewidth = 0.5) +
    geom_point(aes(x = standAge2016, y = 0.3048 * if_else(ODSL_Site_ == 0, ODSL_Physi, pmin(ODSL_Physi, ODSL_Site_)), color = isPlantation, shape = standClass, size = standClass), elliottStandsGisSiteIndex %>% filter(Cruised_Si == 0), alpha = 0.3) +
    labs(x = "stand age, years", y = "lower of ODSL site indices, m", color = NULL, fill = NULL, linetype = NULL, shape = NULL, size = NULL) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    coord_cartesian(ylim = c(0, 75)) &
    guides(color = guide_legend(order = 1), fill = guide_legend(order = 3), linetype = guide_legend(order = 3), shape = guide_legend(order = 2), size = guide_legend(order = 2)) &
    scale_color_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("forestgreen", "green2")) &
    scale_linetype_manual(breaks = c("King's", "GAM"), labels = c("King's curve", "GAM smooth"), values = c("solid", "longdash")) &
    scale_shape_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(16, 17, 18)) &
    scale_size_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(1.5, 1.5, 1.85)) &
    theme(legend.spacing.y = unit(0.4, "line"))
}

# cruised site index versus modeled
if (treetopOptions$includeInvestigatory)
{
  ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 50, yend = 50), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = 0.3048 * Cruised_Si, y = 0.3048 * if_else(ODSL_Site_ == 0, ODSL_Physi, pmin(ODSL_Physi, ODSL_Site_)), color = isPlantation, shape = standClass, size = standClass), elliottStandsGisSiteIndex %>% filter(Cruised_Si > 0), alpha = 0.3) +
  coord_equal(ylim = c(0, 50)) +
  labs(x = "ground measured site index in 2016, m", y = "lower of ODSL site indices, m", color = NULL, shape = NULL, size = NULL) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 3), linetype = guide_legend(order = 3), shape = guide_legend(order = 2), size = guide_legend(order = 2)) +
  scale_color_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("forestgreen", "green2")) +
  scale_linetype_manual(breaks = c("King's", "GAM"), labels = c("King's curve", "GAM smooth"), values = c("solid", "longdash")) +
  scale_shape_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(16, 17, 18)) +
  scale_size_manual(breaks = c("natural regen", "pre-2001 plantation", "2001+ plantation"), values = c(1.5, 1.5, 1.85)) +
  theme(legend.spacing.y = unit(0.4, "line"))
}

# breakout of segmented trees by elevation, slope, and topographic shelter
if (treetopOptions$includeInvestigatory)
{
  ggplot(elliottTreesMod) +
    geom_histogram(aes(x = elevation, fill = isPlantation, group = isPlantation), binwidth = 20) +
    labs(x = "elevation, m", y = "trees") +
    scale_y_continuous(labels = scales::label_comma()) +
  #ggplot(elliottTrees) +
  #  geom_histogram(aes(x = elevation)) +
  ggplot(elliottTreesMod) +
    geom_bin2d(aes(x = aspect, y = slope), binwidth = c(10, 1)) +
    labs(x = "aspect, °", y = "slope, °") +
    scale_x_continuous(breaks = seq(0, 360, by = 90)) +
    scale_fill_viridis_c(limits = c(1, 500E3), trans = "log10") +
  ggplot(elliottTreesMod) +
    geom_histogram(aes(x = topographicShelterIndex, fill = isPlantation, group = isPlantation), binwidth = 1) +
    labs(x = "topographic shelter index, °", y = "trees") +
    scale_y_continuous(labels = scales::label_comma()) +
  plot_layout(guides = "collect") &
    theme(legend.spacing.y = unit(0.3, "line"))
}


## nominal crown radius
if (treetopOptions$includeInvestigatory)
{
  library(gslnls)
  library(mgcv)
  library(quantreg)
  
  crownRadii = elliottTreesMod %>%
    select(standID2016, species, height, CanopyArea, isPlantation, elevation, slope, sinAspect, cosAspect, topographicShelterIndex) %>% 
    rename(projectedCrownArea = CanopyArea) %>%
    mutate(isPlantation = as.factor(isPlantation), # for bam()
           projectedCrownArea = 0.3048^2 * projectedCrownArea, # ft² to m²
           nominalCrownRadius = sqrt(projectedCrownArea / pi)) %>%
    filter(nominalCrownRadius <= height) # exclude physiologically implausible crowns, exclude NA crown areas
  # tibble(height = range(crownRadii$height), crownArea = range(crownRadii$projectedCrownArea))
  
  #crownRadiiLinear = lm(nominalCrownRadius ~ height, crownRadii) # adj R² = 0.651
  crownRadiiLinear = lm(nominalCrownRadius ~ height + I(height^2), crownRadii) # adj R² = 0.688
  #crownRadiiLinear = lm(nominalCrownRadius ~ height + I(height^2) + I(height^3), crownRadii) # adj R² = 0.690, asymmetric residuals as crown radius can't be negative but not clearly heterokedastic, inadequately asymptotic for small and large trees
  #crownRadiiLinear = lm(nominalCrownRadius ~ height + log(height), crownRadii) # adj R² = 0.682 but poorly behaved for height < 8 m
  #crownRadiiLinear = lm(nominalCrownRadius ~ height + I(height^2) + elevation + slope + sinAspect + cosAspect + topographicShelterIndex, crownRadii) # adj R² = 0.682, all predictors significant
  summary(crownRadiiLinear)
  crownRadiiLogistic = gsl_nls(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii)
  summary(crownRadiiLogistic)
  crownRadiiBam = bam(nominalCrownRadius ~ s(height, by = isPlantation, k = 6), data = crownRadii) # ~8 minutes 5950X, default k = 10 likely overfits
  summary(crownRadiiBam)
  k.check(crownRadiiBam)
  crownRadiiLogisticQ05 = nlrq(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii, tau = 0.5) # ~2.5 minutes, internal error on tau = c(0.5, 0.2)
  crownRadiiLogisticQ02 = nlrq(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii, tau = 0.2) # ~2.5 minutes
  crownRadiiLogisticQ01 = nlrq(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii, tau = 0.1) # ~2.5 minutes
  crownRadiiLogisticQ0025 = nlrq(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii, tau = 0.025) # ~2.5 minutes
  crownRadiiLogisticQ0010 = nlrq(nominalCrownRadius ~ SSlogis(height, asym, xmid, scale), crownRadii, tau = 0.01) # ~2.5 minutes
  #summary(crownRadiiLogisticQ05) # spins CPU for unclear duration, not interruptible (requires R session restart to abort)
  
  bind_cols(AIC(crownRadiiLinear, crownRadiiLogistic, crownRadiiBam) %>% mutate(AICn = AIC / nrow(crownRadii)),
            tibble(rmse = c(sqrt(mean(residuals(crownRadiiLinear)^2)), 
                            sqrt(mean(residuals(crownRadiiLogistic)^2)),
                            sqrt(mean(residuals(crownRadiiBam)^2))),
                   nse = 1 - c(sum(residuals(crownRadiiLinear)^2) / sum((crownRadii$height - mean(crownRadii$height))^2),
                               sum(residuals(crownRadiiLogistic)^2) / sum((crownRadii$height - mean(crownRadii$height))^2),
                               sum(residuals(crownRadiiBam)^2) / sum((crownRadii$height - mean(crownRadii$height))^2))))
  
  ggplot() +
    geom_bin2d(aes(x = height, y = nominalCrownRadius), crownRadii, binwidth = c(1, 0.2)) +
    #geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLinear, tibble(height = seq(1.5, 85))), color = "H + H²")) +
    #geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogistic, tibble(height = seq(1.5, 85))), color = "logistic(H)")) +
    #geom_line(aes(x = height, y = predictedCrownRadius, color = "bam(H)", group = isPlantation, linetype = isPlantation), crossing(height = seq(1.5, 85), isPlantation = factor(c(FALSE, TRUE))) %>% mutate(predictedCrownRadius = predict(crownRadiiBam, .))) +
    #geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogisticQ05, tibble(height = seq(1.5, 85))), color = "logistic(H, q = 0.5)")) +
    geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogisticQ02, tibble(height = seq(1.5, 85))), color = "logistic(H, q = 0.2)")) +
    geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogisticQ01, tibble(height = seq(1.5, 85))), color = "logistic(H, q = 0.1)")) +
    geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogisticQ0025, tibble(height = seq(1.5, 85))), color = "logistic(H, q = 0.025)")) +
    #geom_line(aes(x = seq(1.5, 85), y = predict(crownRadiiLogisticQ0010, tibble(height = seq(1.5, 85))), color = "logistic(H, q = 0.010)")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(5.7/(1 + exp((58 - seq(1.5, 85))/20)) / 0.4572, 0), 1), color = "manual 1")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(5.9/(1 + exp((56 - seq(1.5, 85))/17)) / 0.4572, 0), 1), color = "manual 2")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(5.9/(1 + exp((56 - seq(1.5, 85))/16)) / 0.4572, 0), 1), color = "manual 3")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(5.9/(1 + exp((56 - seq(1.5, 85))/15)) / 0.4572, 0), 1), color = "manual 4")) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "height, m", y = "nominal crown radius, m", color = "model", fill = "trees\nsegmented", linetype = "plantation") +
    scale_color_manual(breaks = c("H + H²", "logistic(H)", "logistic(H, q = 0.5)", "logistic(H, q = 0.2)", "logistic(H, q = 0.1)", "logistic(H, q = 0.025)", "logistic(H, q = 0.010)", "manual 1", "manual 2", "manual 3", "manual 4", "bam(H)"), values = c("cyan", "lawngreen", "green1", "green2", "green3", "green4", "darkgreen", "blue", "blue2", "blue3", "blue4","red")) +  scale_fill_viridis_c(trans = "log10") +
    scale_linetype_manual(breaks = c(FALSE, TRUE), values = c("solid", "longdash")) +
    scale_y_continuous(breaks = c(1, 5, 10, 15))
  #ggsave("trees/segmentation/figures/implied crown radius.png", height = 12, width = 16, units = "cm", dpi = 150)
  
  ggplot() +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(8.59/(1 + exp((58 - seq(1.5, 85))/19.42)) / 0.4572, 0), 1), color = "logistic radius")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(pmin(0.055*seq(1.5, 85) + 0.4, 5) / 0.4572, 0), 1), color = "linear radius")) +
    geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round((0.045*seq(1.5, 85) + 0.5) / 0.4572, 0), 1), color = "ring radius")) +
    #geom_line(aes(x = seq(1.5, 85), y = 0.4572 * pmax(round(6.0/(1 + exp((48 - seq(1.5, 85))/18)) / 0.4572, 0), 1), color = "alternate 3")) +
    labs(x = "height, m", y = "nominal crown radius, m", color = "model", fill = "trees\nsegmented", linetype = "plantation") +
    scale_color_discrete() +
    scale_linetype_manual(breaks = c(FALSE, TRUE), values = c("solid", "longdash")) +
    scale_y_continuous(breaks = c(1, 5, 10, 15))
  
  ggplot() +
    geom_hline(yintercept = 1, color = "grey70", linetype = "longdash") +
    geom_line(aes(x = height, y = 12.59/(1+exp((56.56-height)/22.34)) / 0.457, color = "q = 0.5"), tibble(height = seq(0, 100))) +
    geom_line(aes(x = height, y = 9.76/(1+exp((52.60-height)/19.37)) / 0.457, color = "q = 0.2"), tibble(height = seq(0, 100))) +
    geom_line(aes(x = height, y = 8.99/(1+exp((53.17-height)/18.79)) / 0.457, color = "q = 0.1"), tibble(height = seq(0, 100))) +
    geom_line(aes(x = height, y = 8.59/(1+exp((58.72-height)/19.42)) / 0.457, color = "q = 0.025"), tibble(height = seq(0, 100))) +
    guides(color = guide_legend(reverse = TRUE)) +
    labs(x = "height, m", y = "local maxima search radius, 45.7 cm raster cells", color = NULL)
  
  ggplot() +
    geom_bin2d(aes(x = height, y = residuals(crownRadiiBam)), crownRadii, binwidth = c(1, 0.25)) +
    scale_fill_viridis_c(trans = "log10")
}


if (treetopOptions$includeInvestigatory)
{
  # octants from Figure 3 of results dataset.R
  ggplot(distance) +
    geom_raster(aes(x = x, y = y, fill = as.factor(octant))) +
    geom_text(aes(x = x, y = y, label = octant), size = 3) +
    coord_equal() +
    labs(x = "x", y = "y", fill = "octant")
}

## stand distribution
if (treetopOptions$includeInvestigatory)
{
  stands2022 %>% filter(inElliottGis) %>% group_by(isPlantation) %>% 
    summarize(standArea = sum(standArea), .groups = "drop") %>% 
    mutate(standArea = if_else(isPlantation, standArea + 318.8, standArea), # add Hakki plantations
           pctArea = 100 * standArea / sum(standArea))
}