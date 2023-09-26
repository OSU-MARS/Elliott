library(arrow)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)
library(terra)
library(tidyr)
library(writexl)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), panel.border = element_blank()))


## read trees from .gpkg and predict their DBH
# QGIS preparation
#  1) spatially index ESRF_Trees2021 shapefile
#  2) GIS/fieldCalculations.py: join ESRF_Trees2021 shapefile with Elliott_CruiseStands_All_20160111 by location to add standID2016
#
# Single threaded prediction time: 3.9 hours total.
#
#                           Douglas-fir   red alder   western hemlock
# Chapman-Richards physio                 2.5 s
# GAM RelHt                                           30 min
# GAM RelHt physio          2.85 hours                30 min
# Ruark RelHt                                         4 s
# Ruark RelHt physio        4.7 s
#
# Parallel evaluation not viable due to chronic furrr 0.3.1 future_map() failures of the form MultisessionFuture (<none>) failed to call grmall() on cluster RichSOCKnode #1 (PID 21712 on localhost ‘localhost’). The reason reported was ‘error writing to connection’. Post-mortem diagnostic: No process exists with this PID, i.e. the localhost worker is no longer alive.
# Circumstantial evidence suggests worker processes may be exiting due to lack of thread safety in mgcv::predict.gam().
startTime = Sys.time()
elliottTrees = as_tibble(vect("GIS/Trees/2021 OLC Coos County/ESRF_Trees2021.gpkg", layer = "ESRF_Trees2021 StandID2016 physio")) # ~1.3 minutes @ ~2.6 GB, 10.3 million trees
Sys.time() - startTime

stands2022 = read_xlsx("GIS/Trees/2015-16 cruise with 2022 revisions.xlsx") # from height-diameter/setup.R plus manual revisions to stand 1672, 1807, and 2463 area for boundary shifts and slivers (2016 IDs)


# 26 s dplyr + prediction time (~25 seconds for three nonlinear iterations), Zen 3 4.7 GHz
if (recalcDbh)
{
  load("trees/height-diameter/data/segmentationDbhModels.Rdata")
  
  startTime = Sys.time()
  elliottTreesMod = left_join(elliottTrees %>% select(-StandID, -STD_ID), # drop other stand IDs to avoid errors
                              stands2022,
                              by = "standID2016") %>% # ~1 s for join
    rename(treeID = TreeID, TotalHt = Ht) %>% # change Ht to TotalHt to integrate with DBH model fits
    filter(treeID != "NA", TotalHt >= 1.5 * 3.2808) %>%  # remove total of 1271 rows missing TreeIDs, remove 273805 rows below iLand's definition of tree height as 4.0 m (TODO: translate these rows to iLand saplings?)
    mutate(species = factor(species, levels = c("DF", "WH", "HW")),
           TotalHt = 0.3048 * TotalHt,
           x = 0.3048 * x, # convert from English units to metric, assuming input is in EPSG:6557 (nudge trees off resource unit boundaries: - if_else(id %in% c(281406309, 342014985, 412900427), 0.06, 0))
           y = 0.3048 * y,
           elevation = 0.3048 * elevation, # slope, aspect, and topographic shelter index are already in degrees
           resourceUnitX = as.integer(x / 100), # dropped in final select
           resourceUnitY = as.integer(y / 100)) %>% 
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% # put tallest detected trees first in each stand for top height calculation: currently no detection of snags or broken tops
    mutate(treeID = 1E6 * standID2016 + row_number(), # generate unique IDs within 2016 stands: permits 999,999 trees per stand (max segmented is about 432,000 after filtering) as stand ID is four positive digits (32 bit unsigned int max is 4294 967 296 => largest stand ID in .feather is 4294)
           measureTreeTphContribution = 1 / standArea, # since top height trees are the tallest in the stand assume all of them are segmented from LiDAR: expansion factor for top height calculation is thus 1/(delineated area of stand in hectares)
           topHeightTph = pmin(cumsum(if_else(is.na(TotalHt), 0, measureTreeTphContribution)), 100), # TPH total towards the H100 definition of top height, trees not measured for TotalHt are skipped
           topHeightWeight = pmax((topHeightTph - lag(topHeightTph, default = 0)) / measureTreeTphContribution, 0), # clamp remaining fraction to [0, 1] to get individual trees' contributions to the top height average
           topHeight = sum(topHeightWeight * TotalHt, na.rm = TRUE) / sum(topHeightWeight, na.rm = TRUE), # m, tallest 100 trees per hectare
           relativeHeight = TotalHt / topHeight) %>%
    select(-measureTreeTphContribution) %>%
    group_by(species) %>%
    mutate(dbhBootstrap = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkRelHtPhysio, pick(everything())),
                                    cur_group()$species == "HW" ~ predict(alruChapmanRichardsPhysio, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
                                    cur_group()$species == "WH" ~ predict(tsheRuarkRelHt, pick(everything()))),
           # for now, skip GAM-based DBH prediction due to prediction of implausibly slender trees and trees with negative DBHes up to -45 m
           # Difficulty here is the range of GAM DBH underprediction is wide enough it's complex to justify any particular choice
           # of weighting scheme for blending GAM predictions with those from nonlinear regression in order to obtain an ensemble 
           # prediction. Finding an optimal combination predictor appears to be a research topic in its own right.
           #
           # since GAMs are fit with s(by = factor(isPlantation)) predict.gam() ends up calling factor(isPlantation) internally
           # This internal call fails when isPlantation happens to have only one level, workaround is to preemptively
           # convert isPlantation to a two level factor. Levels are specified as c(FALSE, TRUE) to match the order R uses
           # when factor() is called on a boolean variable.
           # Also, Douglas-fir GAM, at least, seems prone to failure with zero rank after access within future_map()
           #dbhGam = case_when(cur_group()$species == "DF" ~ predict(psmeGamRelHtPhysio, pick(everything()) %>% mutate(isPlantation = factor(isPlantation, levels = c(FALSE, TRUE)))), # 37211 DBHes < 3 mm: 0.39% physically impossible
           #                   cur_group()$species == "HW" ~ predict(alruGamRelHtPhysio, pick(everything()) %>% mutate(isPlantation = factor(isPlantation, levels = c(FALSE, TRUE)))), # for now, approximate all hardwoods as red alder
           #                   cur_group()$species == "WH" ~ predict(tsheGamRelHtPhysio, pick(everything()) %>% mutate(isPlantation = factor(isPlantation, levels = c(FALSE, TRUE)))))) %>%
           basalArea = pi/4 * (0.01 * dbhBootstrap)^2) %>% # initial estimate of tree's basal area in m²
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% 
    mutate(standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea, # initial estimate basal area of stand in m²/ha with clamp to plausibility bound as bootstrap DBH estimates can be high, TODO: refine adjustment for undetected trees
           basalAreaAdjustmentFactor = if_else(standBasalAreaApprox <= 125, 1, 125 / standBasalAreaApprox),
           standBasalAreaApprox = basalAreaAdjustmentFactor * standBasalAreaApprox,
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(lag(basalArea, default = 0)) / standArea) %>%  # basal area of taller trees in m²/ha
    group_by(species) %>%
    mutate(dbh = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkAbatPhysio, pick(everything())),
                           cur_group()$species == "HW" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
                           cur_group()$species == "WH" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    group_by(standID2016) %>%
    arrange(desc(TotalHt), .by_group = TRUE) %>% 
    mutate(standBasalAreaBootstrap = standBasalAreaApprox,
           tallerApproxBasalAreaBootstrap = tallerApproxBasalArea,
           basalArea = pi/4 * (0.01 * dbh)^2,
           standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
           basalAreaAdjustmentFactor = if_else(standBasalAreaApprox <= 125, 1, 125 / standBasalAreaApprox),
           standBasalAreaApprox = basalAreaAdjustmentFactor * standBasalAreaApprox,
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(lag(basalArea, default = 0)) / standArea) %>%
    # further iteration results in small changes at most percentiles but drives the smallest <0.5% of trees to negative DBH
    #group_by(species) %>%
    #mutate(dbh = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkAbatPhysio, pick(everything())),
    #                       cur_group()$species == "HW" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
    #                       cur_group()$species == "WH" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    #group_by(standID2016) %>%
    #mutate(standBasalAreaInitial = standBasalAreaApprox,
    #       tallerApproxBasalAreaInitial = tallerApproxBasalArea,
    #       basalArea = pi/4 * (0.01 * dbh)^2,
    #       standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
    #       tallerApproxBasalArea = cumsum(lag(basalArea, default = 0)) / standArea) %>%
    #group_by(species) %>%
    #mutate(dbh2 = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkAbatPhysio, pick(everything())),
    #                        cur_group()$species == "HW" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
    #                        cur_group()$species == "WH" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    #group_by(standID2016) %>%
    #mutate(standBasalArea2 = standBasalAreaApprox,
    #       tallerApproxBasalArea2 = tallerApproxBasalArea,
    #       basalArea = pi/4 * (0.01 * dbh2)^2,
    #       standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
    #       tallerApproxBasalArea = cumsum(lag(basalArea, default = 0)) / standArea) %>%
    #group_by(species) %>%
    #mutate(dbh = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkAbatPhysio, pick(everything())),
    #                       cur_group()$species == "HW" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
    #                       cur_group()$species == "WH" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    ungroup() %>%
    rename(height = TotalHt)
  Sys.time() - startTime
  save(file = "trees/height-diameter/data/segmentationDbh.Rdata", elliottTreesMod)
} else {
  load("trees/height-diameter/data/segmentationDbh.Rdata")
}

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

# check summaries and check plots for dataset alignment and DBH prediction
#elliottTrees %>% filter(is.na(standID2016)) %>% group_by(STD_ID) %>% summarize(n = n())
#elliottTreesMod %>% filter(is.na(dbhBootstrap))
#print(elliottTreesMod %>% select(standID2016, standArea, species, height, topHeightTph, topHeightWeight, topHeight), n = 750)
#elliottStandsMod %>% filter(standID2016 == 12) %>% select(standID2016, tph, segmentedTph, topHeight, segmentedTopHeight)
tibble(tph = cor(drop_na(elliottStandsMod %>% select(tph, segmentedTph)))[2,1], topHeight = cor(drop_na(elliottStandsMod %>% select(topHeight, segmentedQmd)))[2,1], ba = cor(drop_na(elliottStandsMod %>% select(standBasalAreaPerHectare, standBasalAreaApprox.y)))[2,1], qmd = cor(drop_na(elliottStandsMod %>% select(qmd, segmentedQmd)))[2,1])
elliottTreesMod %>% group_by(species) %>% summarize(trees = n(), minDbh = min(dbh), maxDbh = max(dbh), minHt = min(height), maxHt = max(height), maxRelHt = max(relativeHeight), minRelHt = min(relativeHeight), 
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

# cruised site index versus modeled
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


# breakout of segmented trees by elevation, slope, and topographic shelter
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


## write trees for iLand
elliottTreesArrow = arrow_table(elliottTreesMod %>%
                                  mutate(fiaCode = recode(species, "DF" = 202, "WH" = 263, "HW" = 351)) %>% # for now, approximate all hardwoods as red alder
                                  arrange(resourceUnitY, resourceUnitX, species, y, x) %>%
                                  rename(standID = standID2016) %>%
                                  select(standID, treeID, fiaCode, dbh, height, x, y),
                                schema = schema(standID = uint32(), treeID = uint32(), fiaCode = uint16(),
                                                dbh = float32(), height = float32(), x = float32(), y = float32()))
write_feather(elliottTreesArrow, "iLand/init/ESRF trees 2023-05.feather")


## transcode .csv files exported from QGIS to .feather
#tileCsvFileNames = list.files(path = "GIS/Trees", pattern = "^TSegD.*\\.csv")
#tileColumnTypes = cols(id = "i", species = "c", standID = "i", .default = "d")
#for (tileCsvFileName in tileCsvFileNames)
#{
#  featherFilePath = paste0("iLand/init/", str_replace(tileCsvFileName, ".csv", ".feather"))
#  if (file.exists(featherFilePath) == FALSE)
#  {
#    tile = read_csv(paste0("GIS/Trees/", tileCsvFileName), col_types = tileColumnTypes) %>%
#      mutate(resourceUnitX = as.integer(x / 100),
#             resourceUnitY = as.integer(y / 100)) %>%
#      arrange(resourceUnitY, resourceUnitX, species, y, x) %>%
#      select(-resourceUnitY, -resourceUnitX) %>%
#      rename(treeID = id) %>%
#      relocate(standID, treeID, species, height, dbh, x, y)
#    tileArrow = arrow_table(tile %>% mutate(fiaCode = recode(species, "psme" = 202)) %>% select(-species) %>% relocate(standID, treeID, fiaCode, dbh, height, x, y),
#                            schema = schema(standID = int32(), treeID = int32(), fiaCode = uint16(), 
#                                            dbh = float32(), height = float32(), x = float32(), y = float32()))
#    write_feather(tileArrow, featherFilePath)
#  }
#}


## most preferred models for initial DBH prediction: can't use ABA or AAT as DBH hasn't yet been predicted
# height-diameter/setup.R + species data from PSME.R, ALRU2.R, and TSHE.R
psmeGamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 8, folds = 1, repetitions = 1)
psmeRuarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012), folds = 1, repetitions = 1)
psmeRuarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.4, a6 = -0.05, a9 = 0.4, b1 = 0.83, b1p = -0.13, b2 = 0.0042, b2p = 0.008), folds = 1, repetitions = 1)
alruChapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 7, a1p = 16, a5 = 6.6, a6 = 0.7, a7 = 0.43, b1 = -0.033, b2 = 2.6, b2p = -1.16), folds = 1, repetitions = 1)
alruGamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 57, pc = gamConstraint), data = alru2016physio, constraint = alru2016gamConstraint, folds = 1, repetitions = 1)
alruRuarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a5 * sin(3.14159/180 * slope) + a9 * pmin(relativeHeight, 1.5))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016physio, start = list(a1 = 1.1, a2 = -0.01, a5 = 0.5, a9 = -0.4, b1 = 1.5, b1p = -0.30, b2 = -0.04, b2p = 0.028), folds = 1, repetitions = 1)
tsheGamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
tsheGamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 60, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint, nthreads = 8, folds = 1, repetitions = 1)
tsheRuarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.5, a3 = -0.0034, a5 = 0.5, b1 = 0.74, b2 = 0.014), folds = 1, repetitions = 1)
tsheRuarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01), folds = 1, repetitions = 1)
save(file = "trees/height-diameter/data/segmentationDbhModels.Rdata", psmeGamRelHtPhysio, psmeRuarkAbatPhysio, psmeRuarkRelHtPhysio, alruChapmanRichardsPhysio, alruGamRelHtPhysio, alruRuarkAbatPhysioRelHt, tsheGamRelHt, tsheGamRelHtPhysio, tsheRuarkAbatPhysio, tsheRuarkRelHt)


## nominal crown radius
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


## ring prominence
distance = crossing(x = seq(-10, 10), y = seq(-10, 10)) %>% 
  mutate(distance = sqrt(x^2 + y^2),
         ring = round(distance, 0)) %>%
  filter(distance < 10.5)
ggplot(distance) +
  geom_raster(aes(x = x, y = y, fill = as.factor(ring))) +
  geom_text(aes(x = x, y = y, label = ring), size = 3) +
  coord_equal() +
  labs(x = "x", y = "y", fill = "ring") +
  scale_fill_viridis_d()

write_xlsx(distance %>% arrange(ring), "trees/segmentation/rings.xlsx")
