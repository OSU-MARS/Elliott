# assumes trees/area based/setup.R has been run

get_aba_cell_norms = function(abaCellTreeLists, stands2022, stands2021organon, trees2021lidarByHeightClass)
{
  #startTime = Sys.time()
  abaStands = left_join(abaCellTreeLists %>% mutate(heightClass = round(height)) %>% # ~4 s, match height classes in stands2021organon
                          group_by(stand, heightClass) %>%
                          summarize(treesConifer = sum(isConifer), treesHardwood = sum(isConifer == FALSE), .groups = "drop"),
                        stands2022 %>% select(stand, standArea),
                        by = join_by(stand)) %>%
    mutate(tphConifer = treesConifer / standArea, # trees per hectare
           tphHardwood = treesHardwood / standArea)
  #Sys.time() - startTime
  #colSums(is.na(abaStands))
  
  standsAbaCruise = full_join(stands2021organon %>% filter(stand %in% unique(abaStands$stand)) %>% select(stand, heightClass, tphConifer, tphHardwood) %>% rename(tphConiferCruise = tphConifer, tphHardwoodCruise = tphHardwood),
                              abaStands %>% filter(stand %in% unique(stands2021organon$stand)) %>% select(stand, standArea, heightClass, tphConifer, tphHardwood) %>% rename(tphConiferAba = tphConifer, tphHardwoodAba = tphHardwood), # exclude ABA data for uncruised stands
                              by = join_by(stand, heightClass)) %>%
    group_by(stand) %>%
    mutate(standArea = replace_na(standArea, max(standArea, na.rm = TRUE)), # if a height class exists only in cruise data the row will have NA for stand area and ABA tph
           tphConiferAba = replace_na(tphConiferAba, 0),
           tphHardwoodAba = replace_na(tphHardwoodAba, 0),
           tphConiferCruise = replace_na(tphConiferCruise, 0), # similarly, if a height class exists only in LiDAR data the cruise TPH will be NA
           tphHardwoodCruise = replace_na(tphHardwoodCruise, 0)) %>%
    ungroup()
  #colSums(is.na(standsAbaCruise))
  #standsAbaCruise %>% group_by(stand) %>% summarize(standArea = standArea[1], tphConiferCruise = sum(tphConiferCruise), tphHardwoodCruise = sum(tphHardwoodCruise), tphConiferAba = sum(tphConiferAba), tphHardwoodAba = sum(tphHardwoodAba))
  #ggplot() + geom_histogram(aes(x = tph), stands2021organon %>% group_by(stand) %>% summarize(tphConifer = sum(tphConifer), tphHardwood = sum(tphHardwood)) %>% mutate(tph = tphConifer + tphHardwood))
  
  trees2021lidarReference = trees2021lidarByHeightClass %>% filter(stand %in% unique(standsAbaCruise$stand)) %>%
    group_by(isConifer, heightClass) %>%
    summarize(tph = sum(standArea * tph) / sum(standArea), .groups = "drop")
  
  standsAbaCruiseError = standsAbaCruise %>% group_by(heightClass) %>% 
    summarize(tphConiferMad = sum(standArea * abs(tphConiferAba - tphConiferCruise)) / sum(standArea),
              tphConiferMrd = sum(standArea * abs(tphConiferAba - tphConiferCruise) / (pmax(tphConiferAba, tphConiferCruise) + pmax(tphHardwoodAba, tphHardwoodCruise)), na.rm = TRUE) / sum(standArea), # na.rm needed since both ABA and cruise TPH can be zero
              tphConiferRmse = sqrt(sum(standArea * (tphConiferAba - tphConiferCruise)^2) / sum(standArea)),
              tphHardwoodMad = sum(standArea * abs(tphHardwoodAba - tphHardwoodCruise)) / sum(standArea),
              tphHardwoodMrd = sum(standArea * abs(tphHardwoodAba - tphHardwoodCruise) / (pmax(tphConiferAba, tphConiferCruise) + pmax(tphHardwoodAba, tphHardwoodCruise)), na.rm = TRUE) / sum(standArea),
              tphHardwoodRmse = sqrt(sum(standArea * (tphHardwoodAba - tphHardwoodCruise)^2) / sum(standArea)),
              tphConiferAba = sum(standArea * tphConiferAba) / sum(standArea), # calculate totals last as this collapses vectors of stand data to scalars
              tphConiferCruise = sum(standArea * tphConiferCruise) / sum(standArea),
              tphHardwoodAba = sum(standArea * tphHardwoodAba) / sum(standArea),
              tphHardwoodCruise = sum(standArea * tphHardwoodCruise) / sum(standArea))
  #colSums(is.na(standsAbaCruiseError))
  #colSums(standsAbaCruiseError)
  #print(standsAbaCruise %>% filter(heightClass == 6) %>% select(-tphHardwoodCruise, -tphHardwoodAba) %>% mutate(error = tphConiferAba - tphConiferCruise, weightedError = standArea * abs(error) / sum(standArea), mae = sum(weightedError)) %>% relocate(stand, standArea), n = 200)
  return(list(standsAbaCruise = standsAbaCruise, standsAbaCruiseError = standsAbaCruiseError, trees2021lidarReference = trees2021lidarReference))
}

## cross-validated matching
# use vfold_cv() at stand level since group_vfold_cv() departs widely from balanced splits
# splitsAndFits = vfold_cv(plotHeightsScaled, v = 2, repeats = 1) %>% mutate(fit = future_map(splits, fitFunction))
splits = vfold_cv(stands2022 %>% filter(is.na(tph) == FALSE), v = 2, repeats = 2, strata = vegStrata) # cross validate only on cruised stands

split = splits[1, ]
trainingStands = analysis(split$splits[[1]])
validationStands = assessment(split$splits[[1]])

trainingPlotHeightsScaled = plotHeightsScaled %>% filter(stand %in% trainingStands$stand)
validationCells = abaCells %>% filter(stand1 %in% validationStands$stand)
validationCellsScaled = abaCellsScaled %>% filter(stand1 %in% validationStands$stand)

startTime = Sys.time() # ~0.8 s for 2-fold cross validation
abaCellPlots = get_aba_cell_plots(validationCells, validationCellsScaled, trainingPlotHeightsScaled, treeMatchBound = 8, lidarMetrics = c("pGround", "zQ10", "zQ20", "zQ30"))
Sys.time() - startTime
startTime = Sys.time() # ~8 s 2-fold
abaCellTrees = get_tree_lists(abaCellPlots, plotTreeCounts, trees2021lidar)
Sys.time() - startTime
startTime = Sys.time() # ~0.7 s
abaCellNorms = get_aba_cell_norms(abaCellTrees, stands2022, stands2021organon, trees2021lidarByHeightClass)
Sys.time() - startTime

xBreaks = c(0, 0.3, 1, 3, 10, 30, 100) # c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
xMinorBreaks = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2, 4, 5, 6, 7, 8, 9, 20, 40, 50, 60, 70, 80, 90) # c(0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90)
ggplot() +
  geom_col(aes(x = tph, y = heightClass, fill = speciesGroup, group = speciesGroup), abaCellNorms$trees2021lidarReference %>% mutate(speciesGroup = if_else(isConifer, "conifer", "hardwood")), orientation = "y") +
  #coord_cartesian(xlim = c(0, 90), ylim = c(0, 100)) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.1), xlim = c(0, 90), ylim = c(0, 100)) +
  labs(x = "mean trees per hectare", y = "tree height, m", fill = NULL, title = "a) LiDAR segmentation") +
  scale_x_continuous(breaks = xBreaks, minor_breaks = xMinorBreaks, labels = xBreaks) +
ggplot() +
  geom_col(aes(x = tphAba, y = heightClass, fill = speciesGroup, group = speciesGroup), abaCellNorms$standsAbaCruiseError %>% select(heightClass, tphConiferAba, tphHardwoodAba) %>% pivot_longer(cols = c("tphConiferAba", "tphHardwoodAba"), names_to = "speciesGroup", values_to = "tphAba") %>% mutate(speciesGroup = if_else(speciesGroup == "tphConiferAba", "conifer", "hardwood")), orientation = "y") +
  #coord_cartesian(xlim = c(0, 90), ylim = c(0, 100)) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.1), xlim = c(0, 90), ylim = c(0, 100)) +
  labs(x = "mean trees per hectare", y = "tree height, m", fill = NULL, title = "b) LiDAR plus missed trees") +
  scale_x_continuous(breaks = xBreaks, minor_breaks = xMinorBreaks, labels = xBreaks) +
ggplot() +
  geom_col(aes(x = tphCruise, y = heightClass, fill = speciesGroup, group = speciesGroup), abaCellNorms$standsAbaCruiseError %>% select(heightClass, tphConiferCruise, tphHardwoodCruise) %>% pivot_longer(cols = c("tphConiferCruise", "tphHardwoodCruise"), names_to = "speciesGroup", values_to = "tphCruise") %>% mutate(speciesGroup = if_else(speciesGroup == "tphConiferCruise", "conifer", "hardwood")), orientation = "y") +
  #coord_cartesian(xlim = c(0, 90), ylim = c(0, 100)) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.1), xlim = c(0, 90), ylim = c(0, 100)) +
  labs(x = "mean trees per hectare", y = "tree height, m", fill = NULL, title = "c) cruise estimate") +
  scale_x_continuous(breaks = xBreaks, minor_breaks = xMinorBreaks, labels = xBreaks) +
ggplot() +
  geom_col(aes(x = tphMad, y = heightClass, fill = speciesGroup, group = speciesGroup), abaCellNorms$standsAbaCruiseError %>% select(heightClass, tphConiferMad, tphHardwoodMad) %>% pivot_longer(cols = c("tphConiferMad", "tphHardwoodMad"), names_to = "speciesGroup", values_to = "tphMad") %>% mutate(speciesGroup = if_else(speciesGroup == "tphConiferMad", "conifer", "hardwood")), orientation = "y") +
  #coord_cartesian(xlim = c(0, 90), ylim = c(0, 100)) +
  coord_trans(x = scales::transform_pseudo_log(sigma = 0.1), xlim = c(0, 90), ylim = c(0, 100)) +
  labs(x = "trees per hectare", y = "tree height, m", fill = NULL, title = "d) mean absolute difference") +
  scale_x_continuous(breaks = xBreaks, minor_breaks = xMinorBreaks, labels = xBreaks) +
ggplot() +
  geom_col(aes(x = tphMrd, y = heightClass, fill = speciesGroup, group = speciesGroup), abaCellNorms$standsAbaCruiseError %>% select(heightClass, tphConiferMrd, tphHardwoodMrd) %>% pivot_longer(cols = c("tphConiferMrd", "tphHardwoodMrd"), names_to = "speciesGroup", values_to = "tphMrd") %>% mutate(speciesGroup = if_else(speciesGroup == "tphConiferMrd", "conifer", "hardwood")), orientation = "y") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 100)) +
  #coord_trans(x = scales::transform_pseudo_log(sigma = 0.01), ylim = c(0, 100)) +
  labs(x = "mean relative difference", y = NULL, fill = NULL, title = "e) relative difference") +
  #scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), minor_breaks = c(0.03, 0.04, 0.06, 0.07, 0.08, 0.09, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9), labels = scales::label_percent()) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, guides = "collect") &
  scale_fill_manual(breaks = c("conifer", "hardwood"), values = c("forestgreen", "green2")) &
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = 0.02, add = 0)) &
  theme(axis.title = element_text(size = 8), plot.title = element_text(size = 8))
#ggsave("trees/area based/figures/LiDAR+missed trees vs Organon grown ground cruise.png", height = 12, width = 26, units = "cm", dpi = 200)

# 11.4 M treetops detected in matched cells: 40,515 ha in matchable cells, 40,532 ha total ABA grid area
# Elliott (including Hakki Ridge): 33727 ha -> 9.5 M LiDAR treetops, 20.7 M imputed treetops
abaCellPlots %>% summarize(cells = nrow(abaCellPlots), totalArea = abaGrid$sizeHa * cells, metricsCells = sum(is.na(zMean) == FALSE), matchedCells = sum(is.na(plot1n) == FALSE), matchedArea = matchedCells * abaGrid$sizeHa, treetops = sum(n), elliottTreetops = sum((stand1 < 4000) * n), elliottAba = sum((stand1 < 4000) * plot1n, na.rm = TRUE)) %>% mutate(missedTrees = elliottTreetops - elliottAba)
abaCellNorms$standsAbaCruise %>% group_by(stand) %>% summarize(standArea = standArea[1], treesAba = sum(standArea * (tphConiferAba + tphHardwoodAba)), treesCruise = sum(standArea * (tphConiferCruise + tphHardwoodCruise))) %>% ungroup() %>%
  summarize(cruiseArea = sum(standArea), treesAba = sum(treesAba), treesCruise = sum(treesCruise)) %>% 
  mutate(elliottAba = 33408.3 / cruiseArea * treesAba, elliottCruise = 33408.3 / cruiseArea * treesCruise) # excluding 318.8 ha on Hakki Ridge

ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 125, yend = 125), color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_bin2d(aes(x = n, y = plot1n), abaCellPlots, binwidth = c(1, 1)) +
  coord_equal() +
  labs(x = bquote("detected treetops cell"^-1), y = bquote("plot estimate, trees cell"^-1), fill = "plots") +
  scale_fill_viridis_c(labels = scales::label_comma(), trans = "log10")

ggplot() +
  geom_histogram(aes(x = n - plot1n, y = after_stat(..count.. / sum(..count..)), fill = factor(round(isPlantation1)), group = round(isPlantation1)), abaCellPlots, binwidth = 2) +
  labs(x = "tree count mismatch", y = "probability", fill = NULL) +
ggplot() +
  geom_histogram(aes(x = speciesError1matched / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = factor(round(isPlantation1)), group = round(isPlantation1)), abaCellPlots, binwidth = 0.01) +
  labs(x = "species error, normalized", y = NULL, fill = NULL) +
ggplot() +
  geom_histogram(aes(x = heightError1matched / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = factor(round(isPlantation1)), group = round(isPlantation1)), abaCellPlots, binwidth = 1) +
  labs(x = "tree height MAE, m", y = NULL, fill = NULL) +
ggplot() +
ggplot() +
  geom_histogram(aes(x = speciesError1unmatched / (15 - treesMatched), y = after_stat(..count.. / sum(..count..)), fill = factor(round(isPlantation1)), group = round(isPlantation1)), abaCellPlots, binwidth = 0.01) +
  labs(x = "unmatched species error, normalized", y = NULL, fill = NULL) +
  plot_annotation(theme = theme(plot.margin = margin())) +
ggplot() +
  geom_histogram(aes(x = heightError1unmatched / (15 - treesMatched), y = after_stat(..count.. / sum(..count..)), fill = factor(round(isPlantation1)), group = round(isPlantation1)), abaCellPlots, binwidth = 1) +
  labs(x = "unmatched tree height MAE, m", y = NULL, fill = NULL) +
plot_layout(guides = "collect") &
  #coord_cartesian(ylim = c(0, 0.16)) &
  scale_fill_manual(breaks = c(0, 1, NA), labels = c("natural regeneration", "plantation", "unclassified"), values = c("forestgreen", "blue", "grey50")) &
  scale_y_continuous(labels = scales::percent)
  #scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), minor_breaks = c(0.003, 0.004, 0.006, 0.007, 0.008, 0.009, 0.03, 0.04, 0.06, 0.07, 0.08, 0.09), labels = scales::percent, trans = scales::pseudo_log_trans(sigma = 0.001)) # pseudolog doesn't stack bars correctly if there's more than one fill


# cell and match exploration
if (abaOptions$includeInvestigatory)
{
  abaCellPlots %>% group_by(stands) %>% summarize(n = n()) %>% mutate(pct = 100 * n / sum(n))
  
  tibble(plots = nrow(plotHeights), 
         plot1 = length(unique(abaCellPlots$plot1)), shannon1 = vegan::diversity(abaCellPlots$plot1), # 96% of available plots matched 
         plot2 = length(unique(abaCellPlots$plot2)), shannon2 = vegan::diversity(abaCellPlots$plot2), 
         plot3 = length(unique(abaCellPlots$plot3)), shannon3 = vegan::diversity(abaCellPlots$plot3))
  
  ggplot() +
    geom_histogram(aes(x = distance1 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
    labs(x = "normalized kNN distance\nto most similar plot", y = "probability", fill = "trees\nmatched") +
  ggplot() +
    geom_histogram(aes(x = distance2 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
    labs(x = "normalized kNN distance\nto second plot", y = NULL, fill = "trees\nmatched") +
  ggplot() +
    geom_histogram(aes(x = distance3 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
    labs(x = "normalized kNN distance\nto third plot", y = NULL, fill = "trees\nmatched") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.4)) &
    scale_y_continuous(labels = scales::percent)
  
  ggplot() +
    geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..)), fill = "cruise plots,\nunweighted by strata"), plotHeights, binwidth = 1) +
    geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..)), fill = "LiDAR"), abaCells, alpha = 0.5, binwidth = 1) +
    coord_cartesian(xlim = c(0, 30)) +
    labs(x = bquote("trees 20 x 20 m grid cell"^-1), y = "probability", fill = NULL) +
    scale_y_continuous(labels = scales::percent)
  
  left_join(abaCells %>% mutate(nClass = if_else(n < 16, n, 16)) %>% group_by(nClass) %>% summarize(cells = n()) %>% mutate(pctCells = 100 * cells / sum(cells)),
            plotHeights %>% mutate(nClass = if_else(n < 16, n, 16)) %>% group_by(nClass) %>% summarize(plots = n()) %>% mutate(pctPlots = 100 * plots / sum(plots)),
            by = "nClass")
  
  plots2021 = left_join(trees2021organon, plots2016, by = c("plot")) %>%
    # TODO: adjust fixed radius tree count? 1/200 acre fixed radius plots = 217.6 ft² = 20.215702 m² => expansion to 19.7866 * TreeCount
    #select(-PointID, -Dia1, -Ht1, -BHAge, -CrownRatio, -DefectMeasured, -STAND, -elevation, -slope, -aspect, -topographicShelterIndex, -x, -y, -heightDiameterRatio, -imputedBasalArea, -treeBasalAreaPerHectare, -treeBasalAreaPerHectareApprox, -plotsInStand, -standBasalAreaApprox, -measurePlotsInStand, -meanTreesPerBafPlot) %>%
    rename(plot = PlotID) %>%
    arrange(plot, desc(imputedHeight)) %>% 
    group_by(plot) %>%
    mutate(TreeID = row_number(), # renumber trees on plot in descending order of height
           minTreeID = 1, # add columns for cross join on tree ID; TODO: no longer needed, can be removed
           maxTreeID = n()) %>% 
    ungroup()
  #plots2021 %>% group_by(speciesGroup) %>% summarize(n = sum(TreeCount), minHt = min(imputedHeight, na.rm = TRUE), maxHt = max(imputedHeight, na.rm = TRUE), isNAdbh = sum(TreeCount * is.na(DBH)), isNAheight = sum(TreeCount * is.na(imputedHeight)))
  
  ggplot() +
    geom_histogram(aes(x = TotalHt, y = after_stat(..count.. / sum(..count..)), weight = TreeCount, fill = "plots"), plots2021, binwidth = 1) +
    geom_histogram(aes(x = height, y = after_stat(..count.. / sum(..count..)), fill = "cells"), trees2021lidar, alpha = 0.5, binwidth = 1) +
    coord_cartesian(xlim = c(0, 90)) +
    labs(x = "tree height, m", y = "probability", fill = NULL) +
    scale_y_continuous(labels = scales::percent)
  
  trees2016 %>% group_by(PlotID) %>% summarize(hasFixedRadiusPlot = any(SamplingMethod == "FIX"), .groups = "drop") %>% group_by(hasFixedRadiusPlot) %>% summarize(plots = n()) # 1924 of 18363 = 10.5% of plots have trees on nested fixed radius 
  #trees2016 %>% filter(TreeCount > 1, TotalHt > 4.9) %>% group_by(Species) %>% summarize(trees = sum(TreeCount)) # 973 trees with identical DBH, TotalHt, and TreeCount > 1
  
  ggplot() +
    geom_histogram(aes(x = nNonLidarTrees, y = after_stat(..count.. / sum(..count..))), trees2016 %>% filter(imputedHeight < 4.92) %>% group_by(PlotID) %>% summarize(nNonLidarTrees = sum(if_else(SamplingMethod == "FIX", 19.7866, 1) * TreeCount), .groups = "drop"), binwidth = 20) +
    labs(x = "non-LiDAR trees per 20 m ABA grid cell", y = "probability") +
    scale_y_continuous(labels = scales::percent)
}

# plot exploration
if (abaOptions$includeInvestigatory)
{
  ggplot() +
    geom_histogram(aes(x = abaGrid$sizeHa * liveExpansionFactor, y = after_stat(..count.. / sum(..count..)), fill = factor(if_else(species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), species, "other"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other"))), trees2021organon, binwidth = 0.1) +
    coord_cartesian(xlim = c(0, 5)) +
    labs(x = bquote("individual trees plot"^-1), y = "probability", fill = NULL) +
    scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), limits = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_y_continuous(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), labels = scales::percent, trans = scales::pseudo_log_trans(sigma = 0.01))
  
  ggplot() +
    geom_histogram(aes(x = plot1, y = after_stat(..count.. / sum(..count..)), fill = factor(treesMatched), group = treesMatched), abaCellPlots, binwidth = 100) +
    labs(x = "plot ID", y = "probability", fill = "trees\nmatched") +
    scale_y_continuous(labels = scales::percent)
}


# create plot data from 2016 measurements rather than Organon predictions of 2021 growth
if (abaOptions$use2016trees)
{
  load("trees/height-diameter/data/ACMA3 preferred models.Rdata")
  load("trees/height-diameter/data/ALRU2 preferred models.Rdata")
  load("trees/height-diameter/data/other preferred models.Rdata")
  load("trees/height-diameter/data/PSME preferred models.Rdata")
  load("trees/height-diameter/data/THPL preferred models.Rdata")
  load("trees/height-diameter/data/TSHE preferred models.Rdata")
  load("trees/height-diameter/data/UMCA preferred models.Rdata")
  rm(acmaDiameterFromHeightPreferred, alruDiameterFromHeightPreferred, otherDiameterFromHeightPreferred, psmeDiameterFromHeightPreferred, thplDiameterFromHeightPreferred, tsheDiameterFromHeightPreferred, umcaDiameterFromHeightPreferred)
  
  impute_height = function(treeMeasurements)
  {
    # can't pipe to case_match() outside of a mutate statement (dplyr 1.1.3)
    return(case_match(treeMeasurements$speciesGroup,
                      "DF" ~ predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh, treeMeasurements), # GAM BA+L physio?
                      "RA" ~ predict(alruHeightFromDiameterPreferred$sharmaPartonBalPhysio, treeMeasurements),
                      "WH" ~ predict(tsheHeightFromDiameterPreferred$gamBalRelDbh, treeMeasurements),
                      "BM" ~ predict(acmaHeightFromDiameterPreferred$sharmaParton, treeMeasurements),
                      "OM" ~ predict(umcaHeightFromDiameterPreferred$sharmaPartonPhysio, treeMeasurements),
                      "RC" ~ predict(thplHeightFromDiameterPreferred$sharmaPartonPhysio, treeMeasurements),
                      "other" ~ predict(otherHeightFromDiameterPreferred$gamBal, treeMeasurements)))
  }
  
  plots2021 = trees2016 # from height-diameter/setup.R
    filter(PlotType == "IP", is.na(elevation) == FALSE) %>% # exclude dropped plots (IB, CB), count plots, and measure plots missing coordinates as they also lack ABA statistics for kNN
    mutate(imputedHeight = if_else(is.na(TotalHt) & is.na(Ht2), impute_height(.), if_else(is.na(TotalHt), Ht2, TotalHt))) %>% # several seconds due to GAMs
    # Maybe adjust fixed radius tree count? 1/200 acre fixed radius plots = 217.6 ft² = 20.215702 m² => expansion to 19.7866 * TreeCount
    select(-PointID, -Dia1, -Ht1, -BHAge, -CrownRatio, -DefectMeasured, -STAND, -elevation, -slope, -aspect, -topographicShelterIndex, -x, -y, -heightDiameterRatio, -imputedBasalArea, -treeBasalAreaPerHectare, -treeBasalAreaPerHectareApprox, -plotsInStand, -standBasalAreaApprox, -measurePlotsInStand, -meanTreesPerBafPlot) %>%
    rename(plot = PlotID) %>%
    arrange(plot, desc(imputedHeight)) %>% 
    group_by(plot) %>%
    mutate(TreeID = row_number(), # renumber trees on plot in descending order of height
           minTreeID = 1, # add columns for cross join on tree ID; TODO: no longer needed, can be removed
           maxTreeID = n()) %>% 
    ungroup()
}


## correlations between predictors and plot measurements
if (abaOptions$includeInvestigatory)
{
  plotStatistics = plots2021 %>% group_by(plot) %>% 
    summarize(trees = max(TreeID), treeCount = sum(TreeCount), tph = sum(measureTreeTphContribution), 
              h1 = imputedHeight[1], h2 = if_else(n() > 1, imputedHeight[2], NA_real_), h3 = if_else(n() > 2, imputedHeight[3], NA_real_), .groups = "drop")
  plotCorrelations = cor(left_join(plotMetrics2021, plotStatistics, by = "plot"), use = "complete.obs")
  ggplot() +
    geom_tile(aes(x = var1, y = var2, fill = value), as_tibble(plotCorrelations) %>% mutate(var1 = rownames(plotCorrelations)) %>% pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
                mutate(var1 = as_factor(var1), var2 = as_factor(var2))) +
    labs(x = NULL, y = NULL, fill = NULL) +
    scale_fill_scico(palette = "cork", limits = c(-1, 1)) +
    scale_y_discrete(limits = rev) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

## specialized classification and regression random forests for tree count and tallest tree height
if (abaOptions$includeInvestigatory)
{
  library(caret)
  library(ranger)
  repeatedCrossValidation = trainControl(method = "repeatedcv", number = 2, repeats = 5, verboseIter = TRUE)
  
  treeCountData = left_join(plotMetrics2021, plotStatistics %>% select(plot, treeCount), by = "plot") %>% select(-plot)
  #treeCountForest = ranger(treeCount ~ ., treeCountData, classification = TRUE, num.threads = 12) # 32 s
  treeCountForest = train(treeCount ~ ., data = treeCountData %>% sample_n(1000) %>% mutate(treeCount = as_factor(treeCount)), method = "ranger", trControl = repeatedCrossValidation, 
                          tuneGrid = expand.grid(mtry = floor(sqrt(dim(treeCountData)[2])),
                                                 splitrule = "gini",
                                                 min.node.size = 1))
  
  treeCountData = left_join(plotMetrics2021, plotStatistics %>% select(plot, h1), by = "plot") %>% select(-plot)
  h1forest = train(h1 ~ ., data = treeCountData %>% sample_n(1000), method = "ranger", trControl = repeatedCrossValidation, 
                   tuneGrid = expand.grid(mtry = floor(sqrt(dim(treeCountData)[2])),
                                          splitrule = "variance",
                                          min.node.size = 1))
}


## random forest computational intractability of plot matching
if (abaOptions$includeInvestigatory)
{
  # Initial progress message is slow to appear on large fits and overestimates time by an order of magnitude.
  # n       threads   fit time     loopback accuracy, %
  #   100    8          0.1 s       0.9 - perfect recall of known plots
  #   200    8          0.5 s       1.9
  #   500    8          0.5 m       4.8
  #  1000    8          2.8 m       9.6
  # 10365   12        > 17 h
  library(caret)
  library(ranger)
  fitStart = Sys.time()
  plotForest = ranger(plot ~ ., plotMetrics, classification = TRUE, importance = "impurity", num.threads = 12) # 12 threads keep 14 cores at 100%
  Sys.time() - fitStart
  save(plotForest, file = file.path(getwd(), "trees/area based/plotRandomForestFit.Rdata"))
  
  sort(importance(plotForest), decreasing = TRUE)
  predictStart = Sys.time()
  plotPredictions = predict(plotForest, plotMetrics)
  Sys.time() - predictStart
  sum(plotPredictions$predictions == plotMetrics$plot)
  
  plotForest = train(plot ~ ., data = plotMetrics, method = "ranger", trControl = repeatedCrossValidation, 
                    tuneGrid = expand.grid(mtry = ceiling(sqrt(dim(plotMetrics)[2] - 1)),
                                           splitrule = 'gini',
                                           min.node.size = 2))
  
  fitStart = Sys.time()
  plotKknn = train(plot ~ ., data = plotMetrics %>% sample_n(2000), method = "kknn", tuneGrid = expand.grid(kmax = c(13, 15, 17), distance = c(2, 3, 4), kernel = "optimal")) # 18s @ k = 2, r = 5, n = 2000, 6 element grid
  Sys.time() - fitStart
  plotKnn = train(plot ~ ., data = plotMetrics %>% sample_n(2000), method = "knn", tuneGrid = data.frame(k = c(9, 11, 13, 15, 17, 19))) # 11s @ k = 2, r = 5, n = 2000, 6 k values
}


## 2009-2021 treetop matching
if (abaOptions$includeInvestigatory)
{
  library(dplyr)
  library(FNN)
  library(ggplot2)
  library(patchwork)
  library(sf)
  
  theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                               legend.title = element_text(size = 10),
                               panel.border = element_blank()))
  
  treetops2009 = st_read(file.path(getwd(), "GIS/DOGAMI/2009 OLC South Coast/treetops DSM ring/treetops clipped.gpkg"), quiet = TRUE) # slow
  treetops2021 = st_read(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/treetops DSM ring/treetops clipped.gpkg"), quiet = TRUE) # slow
  
  treetops2021neighbors2009 = get.knnx(st_coordinates(treetops2009)[, c("X", "Y")], st_coordinates(treetops2021)[, c("X", "Y")], k = 3) # ~15 s
  
  treetops2021match = st_drop_geometry(treetops2021) %>%
    mutate(x = st_coordinates(treetops2021)[, "X"], y = st_coordinates(treetops2021)[, "Y"], elevation = st_coordinates(treetops2021)[, "Z"],
           height = 0.3048 * height,
           distance1 = 0.3048 * treetops2021neighbors2009$nn.dist[, 1], height2009_1 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 1]],
           distance2 = 0.3048 * treetops2021neighbors2009$nn.dist[, 2], height2009_2 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 2]],
           distance3 = 0.3048 * treetops2021neighbors2009$nn.dist[, 3], height2009_3 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 3]],
           distance = if_else(height > height2009_1, distance1, if_else(distance2 < 0.3048 * 5, distance2, NA_real_)),
           height2009 = if_else(height > height2009_1, height2009_1, if_else(distance2 < 0.3048 * 5, height2009_2, NA_real_)),
           heightIncrement = height - height2009)
  
  treetops2021match %>% filter(distance < pmin(0.1 * height2009, 0.3048 * 5)) %>%
    mutate(classification = if_else(heightIncrement > 24, "> 2 m/year", if_else(heightIncrement > 12, "1–2 m/year", if_else(heightIncrement < -0.5 * height2009, "replaced", if_else(heightIncrement < -12 * 0.01, "broken", "-0.1–1 m/year"))))) %>% 
    group_by(classification) %>% summarize(n = n()) %>% mutate(pct = 100 * n / sum(n), pctAll = 100 * n / nrow(treetops2021)) %>%
    arrange(desc(n))
  
  ggplot(treetops2021match) +
    geom_histogram(aes(x = heightIncrement, y = after_stat(..count.. / sum(..count..))), binwidth = 1) +
    labs(x = "2009–2021 height growth, m", y = "trees") +
    scale_y_continuous(labels = scales::percent)
  
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), color = "grey70", linetype = "longdash") +
    geom_bin_2d(aes(x = height2009, y = height), treetops2021match %>% filter(height >= (height2009 - 3.2808 * 5), height <= (height2009 + 3.2808 * 30), distance < pmin(0.1 * height2009, 0.3048 * 5)), binwidth = c(1, 1)) +
    labs(x = "matched tree height in 2009, m", y = "tree height in 2021, m", fill = "trees\nmatched") +
    scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")
  ggplot() +
    geom_bin_2d(aes(x = height, y = 1/12 * (height - height2009)), treetops2021match %>% filter(height >= (height2009 - 3.2808 * 5), height <= (height2009 + 3.2808 * 30), distance < pmin(0.1 * height2009, 0.3048 * 5)), binwidth = c(1, 0.1)) +
    labs(x = "matched tree height in 2021, m", y = bquote("mean annual height growth, m year"^-1), fill = "trees\nmatched") +
    scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")
  
  
  ggplot(treetops2021match %>% filter(distance1 < 30)) +
    geom_histogram(aes(x = 0.3048 * distance1, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
    labs(x = "distance to nearest treetop, m", y = "") +
  ggplot(treetops2021match %>% filter(distance2 < 30)) +
    geom_histogram(aes(x = 0.3048 * distance2, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
    labs(x = "distance to second nearest treetop, m", y = NULL) +
  ggplot(treetops2021match %>% filter(distance3 < 30)) +
    geom_histogram(aes(x = 0.3048 * distance3, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
    labs(x = "distance to third nearest treetop, m", y = NULL) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() &
    coord_cartesian(xlim = c(0, 30), ylim = c(0, 10)) +
    scale_y_continuous(labels = scales::percent)
}


## LiDAR content checks
if (abaOptions$includeInvestigatory)
{
  library(lidR)
  library(dplyr)
  # Terra
  fishCreek = readLAS("E:/Elliott/GIS/biodiversity plots/FishCreek_Aug2022.las") # all first returns, ±35° scan angle, source ID 0, all flags unused, unclassified
  
  # Remove-Points versus lidR filtering
  s03480w06540lidR = readLAS("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/tiles RGB+NIR/interleave 1/s03480w06540.las")
  s03480w06540lidRf = filter_poi(s03480w06540lidR, (Classification != 7L) & (ReturnNumber != 18L) & (Withheld_flag == FALSE))
  s03480w06540rp = readLAS("F:/Elliott/GIS/DOGAMI/2021 OLC Coos County/tiles RGB+NIR/s03480w06540.las")

  headerDiff = tibble(metadata = (s03480w06540lidRf@header$`File Signature` != s03480w06540rp@header$`File Signature`) + (s03480w06540lidRf@header$`File Source ID` != s03480w06540rp@header$`File Source ID`) + (s03480w06540lidRf@header$`Project ID - GUID` != s03480w06540rp@header$`Project ID - GUID`) + (s03480w06540lidRf@header$`Version Major` != s03480w06540rp@header$`Version Major`) + (s03480w06540lidRf@header$`Version Minor` != s03480w06540rp@header$`Version Minor`) + (s03480w06540lidRf@header$`System Identifier` != s03480w06540rp@header$`System Identifier`) + (s03480w06540lidRf@header$`Generating Software` != s03480w06540rp@header$`Generating Software`) + (s03480w06540lidRf@header$`File Creation Day of Year` != s03480w06540rp@header$`File Creation Day of Year`) + (s03480w06540lidRf@header$`File Creation Year` != s03480w06540rp@header$`File Creation Year`) + (s03480w06540lidRf@header$`Global Encoding`$`GPS Time Type` != s03480w06540rp@header$`Global Encoding`$`GPS Time Type`) + (s03480w06540lidRf@header$`Global Encoding`$`Waveform Data Packets Internal` != s03480w06540rp@header$`Global Encoding`$`Waveform Data Packets Internal`) + (s03480w06540lidRf@header$`Global Encoding`$`Waveform Data Packets External` != s03480w06540rp@header$`Global Encoding`$`Waveform Data Packets External`) + (s03480w06540lidRf@header$`Global Encoding`$`Synthetic Return Numbers` != s03480w06540rp@header$`Global Encoding`$`Synthetic Return Numbers`) + (s03480w06540lidRf@header$`Global Encoding`$WKT != s03480w06540rp@header$`Global Encoding`$WKT) + (s03480w06540rp@header$`Global Encoding`$`Aggregate Model` != s03480w06540rp@header$`Global Encoding`$`Aggregate Model`),
                      structure = (s03480w06540lidRf@header$`Header Size` != s03480w06540rp@header$`Header Size`) + (s03480w06540lidRf@header$`Offset to point data` != s03480w06540rp@header$`Offset to point data`) + (s03480w06540lidRf@header$`Number of variable length records` != s03480w06540rp@header$`Number of variable length records`) + (s03480w06540lidRf@header$`Point Data Format ID` != s03480w06540rp@header$`Point Data Format ID`) + (s03480w06540lidRf@header$`Point Data Record Length` != s03480w06540rp@header$`Point Data Record Length`) + (s03480w06540lidRf@header$`Number of point records` != s03480w06540rp@header$`Number of point records`) + sum(s03480w06540lidRf@header$`Number of points by return` != s03480w06540rp@header$`Number of points by return`),
                      position = (s03480w06540lidRf@header$`X scale factor` != s03480w06540rp@header$`X scale factor`) + (s03480w06540lidRf@header$`Y scale factor` != s03480w06540rp@header$`Y scale factor`) + (s03480w06540lidRf@header$`Z scale factor` != s03480w06540rp@header$`Z scale factor`) + (s03480w06540lidRf@header$`X offset` != s03480w06540rp@header$`X offset`) + (s03480w06540lidRf@header$`Y offset` != s03480w06540rp@header$`Y offset`) + (s03480w06540lidRf@header$`Z offset` != s03480w06540rp@header$`Z offset`) + (s03480w06540lidRf@header$`Max X` != s03480w06540rp@header$`Max X`) + (s03480w06540lidRf@header$`Min X` != s03480w06540rp@header$`Min X`) + (s03480w06540lidRf@header$`Max Y` != s03480w06540rp@header$`Max Y`) + (s03480w06540lidRf@header$`Min Y` != s03480w06540rp@header$`Min Y`) + (s03480w06540lidRf@header$`Max Z` != s03480w06540rp@header$`Max Z`) + (s03480w06540lidRf@header$`Min Z` != s03480w06540rp@header$`Min Z`) + (s03480w06540lidRf@header$`WKT OGC CS`$reserved != s03480w06540rp@header$`WKT OGC CS`$reserved) + (s03480w06540lidRf@header$`WKT OGC CS`$`user ID` != s03480w06540rp@header$`WKT OGC CS`$`user ID`) + (s03480w06540lidRf@header$`WKT OGC CS`$`record ID` != s03480w06540rp@header$`WKT OGC CS`$`record ID`) + (s03480w06540lidRf@header$`WKT OGC CS`$`length after header` != s03480w06540rp@header$`WKT OGC CS`$`length after header`) + (s03480w06540lidRf@header$`WKT OGC CS`$description != s03480w06540rp@header$`WKT OGC CS`$description),
                      wkt = (s03480w06540lidRf@header$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM` != s03480w06540rp@header$`WKT OGC CS`$`WKT OGC COORDINATE SYSTEM`))
  pointDiff = tibble(x = sum(s03480w06540lidRf$X != s03480w06540rp$X), y = sum(s03480w06540lidRf$Y != s03480w06540rp$Y), z = sum(s03480w06540lidRf$Z != s03480w06540rp$Z), intensity = sum(s03480w06540lidRf$Intensity != s03480w06540rp$Intensity),
                     classification = sum(s03480w06540lidRf$Classification != s03480w06540rp$Classification), gpstime = sum(s03480w06540lidRf$gpstime != s03480w06540rp$gpstime), returnNumber = sum(s03480w06540lidRf$ReturnNumber != s03480w06540rp$ReturnNumber),
                     numReturns = sum(s03480w06540lidRf$NumberOfReturns != s03480w06540rp$NumberOfReturns), scanDir = sum(s03480w06540lidRf$ScanDirectionFlag != s03480w06540rp$ScanDirectionFlag), edge = sum(s03480w06540lidRf$EdgeOfFlightline != s03480w06540rp$EdgeOfFlightline),
                     key = sum(s03480w06540lidRf$Keypoint_flag != s03480w06540rp$Keypoint_flag), withheld = sum(s03480w06540lidRf$Withheld_flag != s03480w06540rp$Withheld_flag), overlap = sum(s03480w06540lidRf$Overlap_flag != s03480w06540rp$Overlap_flag),
                     scanAngle = sum(s03480w06540lidRf$ScanAngle != s03480w06540rp$ScanAngle), userData = sum(s03480w06540lidRf$UserData != s03480w06540rp$UserData), pointSource = sum(s03480w06540lidRf$PointSourceID != s03480w06540rp$PointSourceID),
                     r = sum(s03480w06540lidRf$R != s03480w06540rp$R), g = sum(s03480w06540lidRf$G != s03480w06540rp$G), b = sum(s03480w06540lidRf$B != s03480w06540rp$B), nir = sum(s03480w06540lidRf$NIR != s03480w06540rp$NIR))
  headerDiff
  print(pointDiff, width = Inf)
  tibble(zMax = s03480w06540lidR@header$`Max Z`, zMin = s03480w06540lidR@header$`Min Z`, zLidRfMax = s03480w06540lidRf@header$`Max Z`, zLidRfMin = s03480w06540lidRf@header$`Min Z`, zMaxF = s03480w06540rp@header$`Max Z`, zMinF = s03480w06540rp@header$`Min Z`)
}

## Organon growth intervals
if (abaOptions$includeInvestigatory)
{
  #library(readr)
  #organon2021csv = read_csv(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.csv"), col_types = cols(stand = "i", plot = "i", species = "c", tag = "i", year = "i", standAge = "i", .default = "d"))
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  organon2021 = read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE) %>%
    filter(year <= 2021)
  
  ggplot() + 
    geom_bin_2d(aes(x = height2021, y = 1/5 * (height2021 - height2016), weight = liveExpansionFactor2021), organon2021 %>% pivot_wider(id_cols = c(stand, plot, tag, species), names_from = year, names_sep = "", values_from = c(height, liveExpansionFactor)), binwidth = c(1, 0.1)) +
    coord_cartesian(ylim = c(-1.5, 5.75)) +
    labs(x = "predicted tree height in 2021, m", y = bquote("mean annual height growth, m year"^-1), fill = bquote("trees ha"^-1)) +
    scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")
}
