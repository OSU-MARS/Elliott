# assumes library()s from treetops.R setup
library(arrow)
library(readxl)

# TODO: recalc with Elliott State Forest + Hakki stands 2016.gpkg:unified stands 2022 property boundary split to get more exact Elliott area
draw_square = function(data, params, size) # workaround for https://github.com/tidyverse/ggplot2/issues/3669 from https://stackoverflow.com/questions/65812949/set-standard-legend-key-size-with-long-label-names-ggplot
{
  if (is.null(data$size)) 
  {
    data$size <- 0.5
  }
  lwd <- min(data$size, min(size) /4)
  grid::rectGrob(width  = unit(1, "snpc") - unit(lwd, "mm"),
                 height = unit(1, "snpc") - unit(lwd, "mm"),
                 gp = grid::gpar(col = data$colour %||% NA,
                                 fill = alpha(data$fill %||% "grey20", data$alpha),
                                 lty = data$linetype %||% 1,
                                 lwd = lwd * .pt,
                                 linejoin = params$linejoin %||% "mitre",
                                 lineend = if (identical(params$linejoin, "round")) "round" else "square"))
}

# data for figures
# 52 unified, on Elliott stands which are not in ODF.DSL stand definitions
elliottStands2022 = st_drop_geometry(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", layer = "unified stands 2022 property boundary split", quiet = TRUE)) %>%
  group_by(standID2016) %>%
  mutate(onForestArea = if_else((isBuffer == 0) & (isExternalBoundarySplit == 0), standArea, 0), # if isExternalBoundarySplit = 1 then isBuffer = 1 should always occur, but check both in case of GIS flagging error
         totalStandArea = sum(standArea)) %>% # should be the same as grossHa
  filter(isBuffer == 0, # exclude entirely off Elliott stands, should also exclude boundaried out portions of boundary crossing stands 
         isExternalBoundarySplit == 0) %>% # just in case
  summarize(standID2016 = standID2016[1], standAge2016 = standAge2016[1], isPlantation = isPlantation[1], isBuffer = isBuffer[1], isForested = isForested[1], vegetationLabel = vegetationLabel[1], siteIndexInM = siteIndexInM[1], # notes are dropped for now
            onForestArea = sum(standArea), onForestFraction = onForestArea / totalStandArea[1],
            netHa = onForestFraction * netHa[1], # inexact correction but better than no adjustment until/if netHa is updated in GIS
            .groups = "drop") %>%
  mutate(isPlantation = as.logical(isPlantation), isBuffer = as.logical(isBuffer), isForested = as.logical(isForested),
         vegStrata = if_else(vegetationLabel %in% c("1D1L", "1D2H", "1D2L", "1D3H", "1D4H", "1D5H", "DX1L", "DX2H", "DX2L", "DX34L", "DX3H", "DX4H", "DX5H"), vegetationLabel, "other"),
         vegStrata = if_else(is.na(vegStrata) == FALSE, vegStrata, if_else(is.na(siteIndexInM), "1D1L", if_else(standAge2016 < 30, "1D1H", if_else(standAge2016 < 55, "1D2H", "1D3H"))))) # impute recent clearcuts to 1D1L, Hakki to 1D1H or 1D3H depending on age
# vegLabels, grossHa, netHa, and siteIndices are NA in recent clearcuts and on Hakki where  unknown

forest = elliottStands2022 %>% group_by(isPlantation) %>% summarize(stands = n(), areaHa = sum(onForestArea), .groups = "drop") %>% mutate(areaPct = 100 * areaHa / sum(areaHa)) # 16435 + 17350 = 33785 ha, 0.17% over dissolved polygon area of 33727 ha

radius2021 = left_join(read_xlsx("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/radius/standsByHeightClass.xlsx"), # trees by tile, stand, and height class
                       elliottStands2022 %>% select(standID2016, onForestArea, onForestFraction, standAge2016, vegStrata, isPlantation, isBuffer),
                       by = join_by(standID2016)) %>%
  filter(isBuffer == 0) %>%
  mutate(treetops = treetops * onForestFraction) # basic adjustment for boundary crossing stands, TODO: differentiate on and off forest parts of stands in GIS and exclude off forest treetops when compiling standsByHeightClass

randomForest2021 = left_join(read_xlsx("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/rf v2/standsByHeightClass.xlsx"),
                             elliottStands2022 %>% select(standID2016, onForestArea, onForestFraction, standAge2016, vegStrata, isPlantation, isBuffer),
                             by = join_by(standID2016)) %>%
  filter(isBuffer == 0) %>%
  mutate(treetops = treetops * onForestFraction) # same as radius2021

#randomForest2021 %>% filter(treetops > 0) %>% slice_max(heightClassInM, with_ties = TRUE) %>% select(-vegStrata, -isPlantation, -isBuffer)


# data for figures
organon2021 = left_join(read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE) %>%
                          filter(year %in% c(2021, 2026)) %>%
                          group_by(stand, plot, tag) %>%
                          arrange(year, .by_group = TRUE) %>%
                          summarize(species = species[1], standAge = standAge[1] + 1, # model initialized in 2016 with start of growing season, convert to end of 2021 growing season (model year 2022) since LiDAR and imaging flights were August 30, September 13, and 15
                                    dbh = dbh[1] + 0.2 * (dbh[2] - dbh[1]), height = height[1] + 0.2 * (height[2] - height[1]), crownRatio = crownRatio[1] + 0.2 * (crownRatio[2] - crownRatio[1]),
                                    liveExpansionFactor = liveExpansionFactor[1] + 0.2 * (liveExpansionFactor[2] - liveExpansionFactor[1]), deadExpansionFactor = deadExpansionFactor[1] + 0.2 * (deadExpansionFactor[2] - deadExpansionFactor[1]),
                                    .groups = "drop"),
                        elliottStands2022 %>% select(standID2016, onForestArea, vegStrata, isPlantation, isBuffer),
                        by = join_by(stand == standID2016))
organonStands2021 = organon2021 %>% 
  group_by(stand) %>%
  mutate(#plotsInStand = length(unique(plot)),
         #standBasalAreaPerHectare = sum(liveExpansionFactor * 0.25*pi*(0.1*dbh)^2), # m²/ha
         treesPerHectare = sum(liveExpansionFactor),
         snagsPerHectare = sum(deadExpansionFactor)) %>%
  # top height by estimating H100 on each plot and then averaging all plots
  group_by(stand, plot) %>%
  # top height by estimating H100 across all plots
  # group_by(stand) %>%
  arrange(desc(height), .by_group = TRUE) %>% 
  mutate(topHeightTph = pmin(cumsum(liveExpansionFactor + deadExpansionFactor), 100),
         topHeightWeight = pmax((topHeightTph - lag(topHeightTph, default = 0)) / (liveExpansionFactor + deadExpansionFactor), 0)) %>%
  group_by(stand) %>%
  arrange(desc(height), .by_group = TRUE) %>% 
  mutate(topHeightInM = sum(topHeightWeight * (liveExpansionFactor + deadExpansionFactor) * height) / sum(topHeightWeight * (liveExpansionFactor + deadExpansionFactor))) %>%
         #relativeHeight = height / topHeightInM) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(stemsPerHectare = treesPerHectare + snagsPerHectare) %>%
  select(-plot, -tag, -species, -dbh, -crownRatio, -liveExpansionFactor, -deadExpansionFactor, -topHeightTph, -topHeightWeight)
#organonStands2021 %>% group_by(isPlantation) %>% summarize(stands = n())
cruise2016 = organonStands2021 %>% group_by(isPlantation) %>% summarize(stands = n(), areaHa = sum(onForestArea)) %>% mutate(areaPct = 100 * areaHa / sum(areaHa)) # 7193 + 8814 = 16007 ha

randomForestStands2021 = randomForest2021 %>% 
  group_by(standID2016) %>%
  arrange(desc(heightClassInM), .by_group = TRUE) %>%
  mutate(topHeightTreetops = pmin(treetops, 100 * onForestArea),
         topHeightWeight = pmax(topHeightTreetops - lag(topHeightTreetops, default = 0), 0),
         topHeightInM = sum(topHeightWeight * heightClassInM) / sum(topHeightWeight)) %>%
  summarize(treetops = sum(treetops), mergePoints = sum(mergePoints), noisePoints = sum(noisePoints), maybeNoisePoints = sum(maybeNoisePoints),
            onForestArea = onForestArea[1], standAge2016 = standAge2016[1], vegStrata = vegStrata[1], isPlantation = isPlantation[1], isBuffer = isBuffer[1],
            treetopsPerHectare = sum(treetops) / onForestArea[1],
            topHeightInM = topHeightInM[1], .groups = "drop") %>%
  mutate(isInGroundInventory = standID2016 %in% organonStands2021$stand)

combinedStands2021 = full_join(randomForestStands2021 %>% select(standID2016, isPlantation, treetopsPerHectare, topHeightInM, isInGroundInventory) %>% rename(treetopsPerHectareLidar = treetopsPerHectare, topHeightLidar = topHeightInM),
                               organonStands2021 %>% select(stand, topHeightInM, stemsPerHectare) %>% rename(stemsPerHectareGround = stemsPerHectare, topHeightGround = topHeightInM), 
                               by = join_by(standID2016 == stand))
#which((organonStands2021$stand %in% randomForestStands2021$standID2016) == FALSE) # 55  56  78 140 214 215 284 330 471 484
#combinedStands2021 %>% filter(isInGroundInventory) %>% summarize(topHeightCor = cor(topHeightLidar, topHeightGround))
#treeDetectionRegression = lm(treetopsPerHectareLidar ~ stemsPerHectareGround, combinedStands2021 %>% filter(isInGroundInventory))
#summary(treeDetectionRegression) # p << 0.001, adj R² = 0.21

# 286k natural regen + 758k plantation cruised stems -> 6.7 + 14.5 = 21.2 M stems (600 TPH), min height = 2.46 m
# TODO: investigate ODF vegetation strata weighting
organon2021forestTotal = organon2021 %>% group_by(isPlantation) %>% 
  summarize(method = "ground inventory", 
            cruisedStands = length(unique(stand)),
            stemsInCruisedArea = sum(onForestArea * (liveExpansionFactor + deadExpansionFactor)),
            stemsPerHectareSigma = sqrt(sum(onForestArea * var(liveExpansionFactor + deadExpansionFactor)) / sum(onForestArea))) %>% 
  mutate(cruisedAreaHa = cruise2016$areaHa, forestAreaHa = forest$areaHa, 
         stemsInForest = forestAreaHa / cruisedAreaHa * stemsInCruisedArea, 
         stemsPerHectare = stemsInForest / forestAreaHa) # or tph = stemsInCruisedArea / cruisedAreaHa
randomForest2021forestTotal = randomForest2021 %>% filter(heightClassInM >= 2) %>% group_by(isPlantation) %>% summarize(lidarTreetops = sum(treetops)) %>% mutate(forestAreaHa = forest$areaHa, lidarStemsPerHectare = lidarTreetops / forest$areaHa, lidarPct = 100 * lidarStemsPerHectare / c(392, 795)) # 3.1 + 6.6 = 9.7 M -> ~46% detection rate @ 180 treetops/ha natural regen + 362 treetops/ha plantation for 2+ m height classes

# TODO: are ODF veg classes reliable enough for strata?
organon2021densityByHeight = organon2021 %>% mutate(heightClassInM = round(height)) %>% 
  group_by(stand, heightClassInM) %>%
  summarize(isPlantation = isPlantation[1], onForestArea = onForestArea[1], stemsPerHectare = sum(liveExpansionFactor + deadExpansionFactor), .groups = "drop") %>%
  group_by(isPlantation, heightClassInM) %>%
  summarize(cruisedStands = if_else(isPlantation[1], cruise2016$stands[2], cruise2016$stands[1]), # modelAreaWithHeightClassHa = sum(onForestArea), 
            cruiseAreaHa = if_else(isPlantation[1], cruise2016$areaHa[2], cruise2016$areaHa[1]),
            stemsInModel = sum(onForestArea * stemsPerHectare), 
            stemsInModelSigma = sqrt(cruisedStands / (cruisedStands - 1) * (sum((onForestArea * stemsPerHectare - stemsInModel / cruisedStands)^2)) + (cruisedStands - n()) * (0 - stemsInModel / cruisedStands)^2), # calculate standard deviation accounting for stands which aren't in group and therefore have zero stems per hectare
            stemsPerHectare = stemsInModel / cruiseAreaHa, 
            stemsPerHectareSigma = stemsInModelSigma / cruiseAreaHa,
            stemsPerHectare025 = stemsPerHectare + qt(0.025, cruisedStands - 1) * stemsPerHectareSigma,
            stemsPerHectare975 = stemsPerHectare + qt(0.975, cruisedStands - 1) * stemsPerHectareSigma,
            .groups = "drop")

organon2021densityByHeight %>% group_by(isPlantation) %>% summarize(stemsInModel = sum(stemsInModel)) %>%
  mutate(modeledAreaHa = cruise2016$areaHa, stemsPerHectare = stemsInModel / modeledAreaHa, forestAreaHa = forest$areaHa, stemsInForest = stemsInModel * forestAreaHa / modeledAreaHa) %>%
  mutate(method = "ground inventory")

radius2021densityByHeight = radius2021 %>%
  group_by(method, standID2016, heightClassInM) %>%
  summarize(isPlantation = isPlantation[1], onForestArea = onForestArea[1], treetops = sum(treetops), .groups = "drop") %>%
  group_by(method, isPlantation, heightClassInM) %>%
  summarize(forestAreaHa = if_else(isPlantation[1], forest$areaHa[2], forest$areaHa[1]), forestAreaWithHeightClassHa = sum(onForestArea), treetops = sum(treetops), treetopsPerHectare = treetops / forestAreaHa, .groups = "drop")

randomForest2021densityByHeight = randomForest2021 %>% # mutate(heightClassInM = heightClassSize * round(heightClassInM / heightClassSize)) %>% 
  group_by(standID2016, heightClassInM) %>%
  summarize(isPlantation = isPlantation[1], onForestArea = onForestArea[1], treetops = sum(treetops), .groups = "drop") %>%
  group_by(isPlantation, heightClassInM) %>%
  summarize(forestAreaHa = if_else(isPlantation[1], forest$areaHa[2], forest$areaHa[1]), forestAreaWithHeightClassHa = sum(onForestArea), treetops = sum(treetops), treetopsPerHectare = treetops / forestAreaHa, .groups = "drop") %>%
  mutate(method = "DSM forest")

randomForest2021densityByHeight %>% group_by(isPlantation) %>% summarize(lidarTreetops = sum(treetops), treetopsPerHectare = sum(treetopsPerHectare), inventoryAreaHa = max(forestAreaWithHeightClassHa)) %>% 
  mutate(forestAreaHa = forest$areaHa, forestTreetops = forestAreaHa * treetopsPerHectare, areaDriftPct = 100 * (forestAreaHa / inventoryAreaHa - 1), treetopsDriftPct = 100 * (forestTreetops / lidarTreetops - 1)) # 3.1 + 6.9 = 10.0 M

combined2021densityByHeight = full_join(randomForest2021densityByHeight, 
                                        organon2021densityByHeight, 
                                        by = join_by(isPlantation, heightClassInM))


## Figure 11: whole forest inventory summary
stems2021 = bind_rows(randomForest2021, 
                      radius2021, 
                      organon2021forestTotal %>% 
                        mutate(stemsPerHectare025 = stemsPerHectare + qt(0.025, cruisedStands - 1) * stemsPerHectareSigma,
                               stemsPerHectare975 = stemsPerHectare + qt(0.975, cruisedStands - 1) * stemsPerHectareSigma,
                               stems025 = stemsPerHectare025 * forestAreaHa,
                               stems975 = stemsPerHectare975 * forestAreaHa) %>%
                        select(method, isPlantation, stemsInForest, stems025, stems975, stemsPerHectare025, stemsPerHectare975) %>% 
                        rename(treetops = stemsInForest)) %>% 
  group_by(method, isPlantation) %>%
  summarize(stems = sum(treetops), 
            stems025 = stems025[1], stems975 = stems975[1], stemsPerHectare025 = stemsPerHectare025[1], stemsPerHectare975 = stemsPerHectare975[1],
            .groups = "drop_last") %>%
  mutate(method = factor(method, levels = c("DSM forest", "DSM radius", "CHM radius", "CMM radius", "ground inventory")),
         isPlantation = factor(isPlantation, levels = c(TRUE, FALSE)),
         forestAreaHa = forest$areaHa,
         stemsPerHectare = stems / forestAreaHa,
         methodLabel = forcats::fct_recode(factor(method, levels = rev(c("DSM forest", "DSM radius", "CHM radius", "CMM radius", "ground inventory"))), `Organon grown\nground inventory` = "ground inventory"),
         stemsLabelX = 0.5 * stems + lag(stems, n = 1, default = 0)) %>%
  group_by(isPlantation) %>%
  mutate(stemsPerHectareGroundPct = 100 * stemsPerHectare / max(stemsPerHectare))

stems2021 %>% group_by(method) %>% summarize(nrStems = stems[1], plantationStems = stems[2], stems = sum(stems), .groups = "drop") %>% mutate(stemsPct = 100 * stems / max(stems))

ggplot() +
  geom_col(aes(x = stems, y = methodLabel, fill = isPlantation, group = isPlantation), stems2021, key_glyph = draw_square) +
  #geom_errorbarh(aes(xmin = stems025 + if_else(isPlantation == TRUE, stems[1], 0), xmax = stems975 + if_else(isPlantation == TRUE, stems[1], 0), y = methodLabel, group = isPlantation), stems2021 %>% filter(method == "ground inventory"), color = "grey30", height = 0.33, linewidth = 0.3) +
  geom_text(aes(x = stemsLabelX, y = methodLabel, color = isPlantation, label = sprintf("%.2f", 1E-6 * stems)), stems2021, size = 2.7) +
  geom_text(aes(x = stems, y = methodLabel, label = sprintf(if_else(methodLabel != "DSM forest", "  %.2f", "  %.2f total"), 1E-6 * stems)), stems2021 %>% group_by(methodLabel) %>% summarize(stems = sum(stems)), hjust = 0, size = 2.7) +
  coord_cartesian(xlim = c(0, NA), clip = "off") +
  labs(x = "millions of treetops", y = NULL, color = NULL, fill = NULL, title = paste("                        ", plotLetters[1], "individual tree inventory")) +
  scale_x_continuous(labels = scales::label_number(scale = 1E-6)) +
  theme(plot.margin = margin(l = 5, r = 25), axis.title.x = element_text(vjust = 0.55)) +
ggplot() +
  geom_col(aes(x = stemsPerHectare, y = methodLabel, fill = isPlantation), stems2021 %>% filter(isPlantation == FALSE)) +
  #geom_errorbarh(aes(xmin = stemsPerHectare025, xmax = stemsPerHectare975, y = methodLabel), stems2021 %>% filter(isPlantation == FALSE), color = "grey30", height = 0.33, linewidth = 0.3) +
  geom_text(aes(x = 0.5 * stemsPerHectare, y = methodLabel, color = isPlantation, label = sprintf("%0.0f", stemsPerHectare)), stems2021 %>% filter(isPlantation == FALSE), size = 2.7) +
  geom_text(aes(x = stemsPerHectare, y = methodLabel, label = if_else(methodLabel == "Organon grown\nground inventory", "", sprintf(if_else(methodLabel != "DSM forest", "  %0.1f%%", "  %0.1f%% of ground"), stemsPerHectareGroundPct))), stems2021 %>% filter(isPlantation == FALSE), hjust = 0, size = 2.7) +
  coord_cartesian(xlim = c(0, 800)) +
  guides(fill = "none") +
  labs(x = bquote("mean treetops ha"^-1), y = NULL, fill = NULL, title = paste(plotLetters[2], "natural regeneration")) +
  scale_y_discrete(labels = NULL) +
ggplot() +
  geom_col(aes(x = stemsPerHectare, y = methodLabel, fill = isPlantation), stems2021 %>% filter(isPlantation == TRUE)) +
  #geom_errorbarh(aes(xmin = stemsPerHectare025, xmax = stemsPerHectare975, y = methodLabel), stems2021 %>% filter(isPlantation == TRUE), color = "grey30", height = 0.33, linewidth = 0.3) +
  geom_text(aes(x = 0.5 * stemsPerHectare, y = methodLabel, color = isPlantation, label = sprintf("%0.0f", stemsPerHectare)), stems2021 %>% filter(isPlantation == TRUE), size = 2.7) +
  geom_text(aes(x = stemsPerHectare, y = methodLabel, label = if_else(methodLabel == "Organon grown\nground inventory", "", sprintf(if_else(methodLabel != "DSM forest", "  %0.1f%%", "  %0.1f%% of ground"), stemsPerHectareGroundPct))), stems2021 %>% filter(isPlantation == TRUE), hjust = 0, size = 2.7) +
  coord_cartesian(xlim = c(0, 800), clip = "off") +
  guides(fill = "none") +
  labs(x = bquote("mean treetops ha"^-1), y = NULL, fill = NULL, title = paste(plotLetters[3], "plantation density")) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, widths = c(1, 0.75, 0.75), guides = "collect") &
  guides(color = "none") &
  scale_color_manual(breaks = c(FALSE, TRUE), values = c("white", "black")) &
  scale_fill_manual(breaks = c(FALSE, TRUE), labels = c("naturally\nregenerated\nstands", "plantations"), values = c("forestgreen", "#5588FF")) # https://davidmathlogic.com/colorblind/, https://www.nceas.ucsb.edu/sites/default/files/2022-06/Colorblind%20Safe%20Color%20Schemes.pdf, https://rdrr.io/cran/khroma/, ...
  # scale_fill_manual(breaks = c("DSM forest", "DSM radius", "CHM radius", "CMM radius"), values = c("#336887FF", "#8197A4FF", "#706A6BFF", "#B7AA9FFF")) # calecopal::casj
#ggsave("trees/segmentation/treetops/figures/Figure 11 forest inventory summary.png", height = 6, width = 20, units = "cm", dpi = figureDpi)


## Figure 12: whole forest inventory by height class
ggplot() +
  geom_ribbon(aes(xmin = stemsPerHectare025, xmax = stemsPerHectare975, y = heightClassInM, fill = "95% confidence"), combined2021densityByHeight %>% filter(isPlantation == FALSE), alpha = 0.15) +
  geom_segment(aes(x = 0, y = 4.5, xend = 45, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  geom_line(aes(x = treetopsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = method, group = method), radius2021densityByHeight %>% filter(isPlantation == FALSE, treetopsPerHectare > 0), orientation = "y") +
  geom_line(aes(x = stemsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = "Organon grown\nground inventory"), combined2021densityByHeight %>% filter(isPlantation == FALSE, stemsPerHectare > 0), orientation = "y") +
  geom_line(aes(x = treetopsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = "DSM forest"), combined2021densityByHeight %>% filter(isPlantation == FALSE, treetopsPerHectare > 0), orientation = "y") +
  labs(x = bquote("stems ha"^-1), y = "height above ground, m", alpha = NULL, color = NULL, fill = NULL, title = paste("    ", plotLetters[1], "natural regeneration")) +
ggplot() +
  geom_ribbon(aes(xmin = stemsPerHectare025, xmax = stemsPerHectare975, y = heightClassInM, fill = "95% confidence"), combined2021densityByHeight %>% filter(isPlantation), alpha = 0.15) +
  geom_segment(aes(x = 0, y = 4.5, xend = 45, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  geom_line(aes(x = treetopsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = method, group = method), radius2021densityByHeight %>% filter(isPlantation, treetopsPerHectare > 0), orientation = "y") +
  geom_line(aes(x = stemsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = "Organon grown\nground inventory"), combined2021densityByHeight %>% filter(isPlantation, stemsPerHectare > 0), orientation = "y") +
  geom_line(aes(x = treetopsPerHectare, y = heightClassInM, alpha = heightClassInM >= 5, color = "DSM forest"), combined2021densityByHeight %>% filter(isPlantation, treetopsPerHectare > 0), orientation = "y") +
  labs(x = bquote("stems ha"^-1), y = NULL, alpha = NULL, color = NULL, fill = NULL, title = paste(plotLetters[2], "plantation density")) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, guides = "collect") &
  coord_cartesian(xlim = c(0, 45), ylim = c(0, 105)) &
  #coord_trans(x = scales::pseudo_log_trans(), xlim = c(0, 45), ylim = c(0, 105)) &
  guides(alpha = "none", color = guide_legend(order = 1), fill = guide_legend(order = 2)) &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_color_manual(breaks = c("DSM forest", "DSM radius", "CHM radius", "CMM radius", "Organon grown\nground inventory"), values = c("dodgerblue", "grey50", "grey65", "grey80", "forestgreen")) &
  scale_fill_manual(breaks = c("95% confidence"), labels = c("95% confidence\ninterval"), values = c("forestgreen")) &
  #scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 10, 20, 30, 50), minor_breaks = c(4, 6, 7, 8, 9, 40)) &
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 1))
#ggsave("trees/segmentation/treetops/figures/Figure 12 forest by height class.png", height = 17, width = 17, units = "cm", dpi = figureDpi)
#ggsave("trees/segmentation/treetops/figures/Figure 12 forest by height class wide.png", height = 10, width = 13, units = "cm", dpi = figureDpi)


if (treetopOptions$includeInvestigatory)
{
  # check plots
  ggplot() +
    geom_histogram(aes(stemsPerHectare, fill = TRUE), organonStands2021, binwidth = 25) +
    coord_cartesian(xlim = c(0, 3500), ylim = c(0, 300)) +
    guides(fill = "none") +
    labs(x = bquote("ground stems ha"^-1), y = "stands", fill = NULL, title = paste(plotLetters[1], "ground inventory grown by Organon")) +
  ggplot() +
    geom_histogram(aes(treetopsPerHectare, fill = isInGroundInventory, group = isInGroundInventory), randomForestStands2021, binwidth = 25) +
    coord_cartesian(xlim = c(0, 2200), ylim = c(0, 300)) +
    labs(x = bquote("LiDAR treetops ha"^-1), y = "stands", fill = NULL, title = paste(plotLetters[2], "LiDAR detections")) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 3500, yend = 3500), color = "grey20", linewidth = 0.3) +
    geom_point(aes(x = stemsPerHectareGround, y = treetopsPerHectareLidar, color = TRUE, shape = isPlantation, size = isPlantation), combinedStands2021, alpha = 0.25) +
    coord_equal(xlim = c(0, 3500), ylim = c(0, NA)) +
    labs(x = bquote("ground stems ha"^-1), y = bquote("LiDAR treetops ha"^-1), color = NULL, shape = NULL, size = NULL, title = paste(plotLetters[3], "")) +
  ggplot() +
    geom_point(aes(x = standAge, y = topHeightInM, color = TRUE, shape = isPlantation, size = isPlantation), organonStands2021, alpha = 0.25) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 55)) +
    labs(x = "stand age, years", y = "top height, m", color = NULL, shape = NULL, size = NULL, title = paste(plotLetters[4], "heights grown from ground measurement")) +
    scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  ggplot() +
    geom_point(aes(x = standAge2016 + 5, y = topHeightInM, color = isInGroundInventory, shape = isPlantation, size = isPlantation), randomForestStands2021 %>% filter(isInGroundInventory), alpha = 0.25) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 55)) +
    labs(x = "stand age, years", y = "top height, m", color = NULL, shape = NULL, size = NULL, title = paste(plotLetters[5], "LiDAR height measurement")) +
    scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 55, yend = 55), color = "grey20", linewidth = 0.3) +
    geom_point(aes(x = topHeightGround, y = topHeightLidar, color = TRUE, shape = isPlantation, size = isPlantation), combinedStands2021, alpha = 0.25) +
    coord_equal(xlim = c(0, 55), ylim = c(0, 55)) +
    labs(x = "ground top height, m, adjusted", y = "LiDAR top height, m", color = NULL, shape = NULL, size = NULL, title = paste(plotLetters[6], "")) +
    scale_y_continuous(breaks = seq(0, 60, by = 10)) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 3, heights = c(1, 1), guides = "collect") &
    guides(color = "none", shape = guide_legend(override.aes = list(alpha = 0.8))) &
    scale_color_manual(breaks = c(TRUE, FALSE), labels = c("ground and LiDAR", "LiDAR only"), values = c("forestgreen", "blue")) &
    scale_fill_manual(breaks = c(TRUE, FALSE), labels = c("ground and LiDAR", "LiDAR only"), values = c("forestgreen", "blue")) &
    scale_shape_manual(breaks = c(FALSE, TRUE), labels = c("naturally regenerated", "plantation"), values = c(16, 18)) &
    scale_size_manual(breaks = c(FALSE, TRUE), labels = c("naturally regenerated", "plantation"), values = c(1, 1.7))
    
  ggplot() +
    geom_point(aes(x = dbh, y = height), organon2021, alpha = 0.2, shape = 16) +
    labs(x = "DBH, cm", y = "height, m")
}
