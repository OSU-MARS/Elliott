library(arrow)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)
library(terra)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), panel.border = element_blank()))

load("trees/height-diameter/data/segmentationDbhModels.Rdata")


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
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(replace_na(lag(basalArea), 0)) / standArea) %>%  # basal area of taller trees in m²/ha
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
           tallerApproxBasalArea = basalAreaAdjustmentFactor * cumsum(replace_na(lag(basalArea), 0)) / standArea) %>%
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
    #       tallerApproxBasalArea = cumsum(replace_na(lag(basalArea), 0)) / standArea) %>%
    #group_by(species) %>%
    #mutate(dbh2 = case_when(cur_group()$species == "DF" ~ predict(psmeRuarkAbatPhysio, pick(everything())),
    #                        cur_group()$species == "HW" ~ predict(alruRuarkAbatPhysioRelHt, pick(everything())), # for now, approximate all hardwoods as red alder: all predictions physically possible
    #                        cur_group()$species == "WH" ~ predict(tsheRuarkAbatPhysio, pick(everything())))) %>%
    #group_by(standID2016) %>%
    #mutate(standBasalArea2 = standBasalAreaApprox,
    #       tallerApproxBasalArea2 = tallerApproxBasalArea,
    #       basalArea = pi/4 * (0.01 * dbh2)^2,
    #       standBasalAreaApprox = 1 / 1 * sum(basalArea) / standArea,
    #       tallerApproxBasalArea = cumsum(replace_na(lag(basalArea), 0)) / standArea) %>%
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

elliotStandsMod = left_join(elliottTreesMod %>% group_by(standID2016) %>% 
                              summarize(segmentedTph = n() / standArea[1],
                                        segmentedTopHeight = topHeight[1],
                                        segmentedQmd = sqrt(standArea[1] * standBasalAreaApprox[1] / (pi/4 * 0.01^2 * n())), # basal area is in m²/ha so need to multiply by stand area since n() counts all trees in the stand
                                        standBasalAreaBootstrap = standBasalAreaBootstrap[1], 
                                        #standBasalAreaInitial = standBasalAreaInitial[1],
                                        standBasalAreaApprox = standBasalAreaApprox[1]),
                            stands2022,
                            by = "standID2016")

# check summaries and check plots for dataset alignment and DBH prediction
#elliottTrees %>% filter(is.na(standID2016)) %>% group_by(STD_ID) %>% summarize(n = n())
#elliottTreesMod %>% filter(is.na(dbhBootstrap))
tibble(tph = cor(drop_na(elliotStandsMod %>% select(tph, segmentedTph)))[2,1], topHeight = cor(drop_na(elliotStandsMod %>% select(topHeight, segmentedQmd)))[2,1], ba = cor(drop_na(elliotStandsMod %>% select(standBasalAreaPerHectare, standBasalAreaApprox.y)))[2,1], qmd = cor(drop_na(elliotStandsMod %>% select(qmd, segmentedQmd)))[2,1])
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
  coord_cartesian(xlim = c(0, 700), ylim = c(0, 96)) +
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
plot_layout(nrow = 1, widths = c(7, 3), guides = "collect") &
  theme(legend.spacing.y = unit(0.3, "line"))

# comparison of stand-level properties
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 3000, yend = 3000), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = tph, y = segmentedTph, color = isPlantation), elliotStandsMod, alpha = 0.3, na.rm = TRUE, shape = 16) +
  coord_equal() +
  labs(x = "2015–16 trees per hectare", y = "trees segmented per hectare", color = "plantation") +
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 80, yend = 80), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = topHeight, y = segmentedTopHeight, color = isPlantation), elliotStandsMod, alpha = 0.3, na.rm = TRUE, shape = 16) +
  coord_equal() +
  labs(x = bquote("2015–16 H"[100]*", m"), y = bquote("segmented H"[100]*", m"), color = "plantation") +
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 125, yend = 125), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = standBasalAreaPerHectare, y = standBasalAreaApprox.y, color = isPlantation), elliotStandsMod, alpha = 0.3, na.rm = TRUE, shape = 16) +
  coord_equal() +
  labs(x = bquote("2015–16 basal area, m"^2*" ha"^-1), y = bquote("segmented basal area, m"^2*" ha"^-1), color = "plantation") +
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = qmd, y = segmentedQmd, color = isPlantation), elliotStandsMod, alpha = 0.3, na.rm = TRUE, shape = 16) +
  coord_equal() +
  labs(x = "2015–16 QMD, cm", y = "segmented QMD, cm", color = "plantation") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(widths = c(1, 1), heights = c(1, 1), guides = "collect")
  
ggplot() +
  geom_histogram(aes(x = 100 * segmentedTph / tph, y = after_stat(100 * count / sum(count)), fill = isPlantation), elliotStandsMod, binwidth = 5, na.rm = TRUE) +
  geom_line(aes(x = seq(0, 300), y = 5 * dgamma(0.01 * seq(0, 300), shape = 3.2, rate = 5)), color = "grey30") +
  coord_cartesian(xlim = c(0, 200)) +
  labs(x = "fraction of trees segmented, %", y = "fraction of stands, %", fill = "plantation") +
  theme(legend.spacing.y = unit(0.3, "line"))

# random stand selection for faster initial analysis
# 
elliotStandsModSubset = elliotStandsMod %>% group_by(isPlantation) %>% sample_n(10) %>% arrange(standID2016, .by_group = TRUE)
elliottTreesNaturalRegenSubset = elliottTreesMod %>% filter(standID2016 %in% c(18, 1051, 1059, 1102, 1172, 1319, 1391, 1399, 1645, 1984))
elliottTreesPlantationSubset = elliottTreesMod %>% filter(standID2016 %in% c(152, 222, 769, 979, 1222, 1483, 1555, 2373, 2518)) # 1470: 15.3 kTrees

ggplot(elliottTreesMod %>% filter(isPlantation == FALSE)) +
  geom_histogram(aes(x = relativeHeight, y = after_stat(100 * count/sum(count)), fill = species), binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 11)) +
  labs(x = NULL, y = "fraction of trees, %", title = "a) 2021 natural regen", fill = NULL) +
ggplot(trees2016 %>% filter(isPlantation == FALSE)) +
  geom_histogram(aes(x = relativeHeight, y = after_stat(100 * count/sum(count)), fill = speciesGroup), binwidth = 0.05, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 11)) +
  labs(x = NULL, y = NULL, title = "b) 2015–16 natural regen", fill = NULL) +
ggplot(elliottTreesMod %>% filter(isPlantation)) +
  geom_histogram(aes(x = relativeHeight, y = after_stat(100 * count/sum(count)), fill = species), binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 11)) +
  labs(x = "relative height", y = "fraction of trees, %", title = "c) 2021 plantation", fill = NULL) +
ggplot(trees2016 %>% filter(isPlantation)) +
  geom_histogram(aes(x = relativeHeight, y = after_stat(100 * count/sum(count)), fill = speciesGroup), binwidth = 0.05, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 11)) +
  labs(x = "relative height", y = NULL, title = "d) 2015–16 plantation", fill = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect")

ggplot(trees2016) +
  geom_histogram(aes(x = after_stat(100 * count / sum(count)), y = TotalHt, fill = isPlantation), binwidth = 1, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 83)) +
  labs(x = "fraction of trees, %", y = "tree height, m", fill = "plantation", title = "a) 2015–16") +
ggplot(elliottTreesMod) +
  geom_histogram(aes(x = after_stat(100 * count / sum(count)), y = height, fill = isPlantation), binwidth = 1) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 83)) +
  labs(x = "fraction of trees, %", y = NULL, fill = "plantation", title = "b) 2021") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect")

elliottTrees %>% reframe(quantiles = c(0, 0.005, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.995, 1), height = quantile(Ht, probs = quantiles, na.rm = TRUE))
trees2016 %>% reframe(quantiles = c(0, 0.005, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.995, 1), height = quantile(TotalHt, probs = quantiles, na.rm = TRUE), dbh = quantile(DBH, probs = quantiles, na.rm = TRUE))

#ggplot(elliottTreesMod) +
#  geom_bin2d(aes(x = dbhGam, y = height), binwidth = c(2.5, 1)) +
#  coord_cartesian(xlim = c(0, NA)) +
#  labs(x = "predicted DBH, cm", y = "segmented height, m", fill = "trees) +
#  scale_fill_viridis_c(labels = scales::label_comma(), limits = c(1, 500E3), trans = "log10") +
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

ggplot(elliottTreesMod) +
  geom_histogram(aes(y = height, fill = species), binwidth = 0.5) +
  labs(x = "trees segemented", y = "tree height, m", fill = NULL) +
  scale_fill_manual(breaks = c("DF", "WH", "HW"), values = c("forestgreen", "blue2", "red2")) +
  scale_x_continuous(labels = scales::label_comma())


## write trees for iLand
elliottTreesArrow = arrow_table(elliottTreesMod %>%
                                  mutate(fiaCode = recode(species, "DF" = 202, "WH" = 263, "HW" = 351)) %>% # for now, approximate all hardwoods as red alder
                                  select(standID, treeID, fiaCode, dbh, height, x, y) %>%
                                  arrange(resourceUnitY, resourceUnitX, species, y, x),
                                schema = schema(standID = uint32(), treeID = uint32(), fiaCode = uint16(),
                                                dbh = float32(), height = float32(), x = float32(), y = float32()))
write_feather(elliottTreesArrow, "iLand/init/ESRF_Trees2021.feather", compression = "uncompressed")


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
#    write_feather(tileArrow, featherFilePath, compression = "uncompressed")
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
