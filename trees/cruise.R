library(dplyr)
library(ggplot2)
library(patchwork)
library(terra)
library(readxl)
library(stringr)

# trees2016 and stands2022 from trees/height-diameter/setup.R
cruiseTph2016 = left_join(stands2022 %>% select(StandID, ODSL_VEG_L, area_ha, standAge2016, isPlantation) %>% rename(vegLabel = ODSL_VEG_L, standArea = area_ha),
                          trees2016 %>% group_by(StandID) %>% 
                            summarize(tph = tph[1], topHeight = topHeight[1], totalBasalArea = basalArea[1],
                                      psmeBA = sum((speciesGroup == "DF") * treeBasalAreaPerHectare) / plotsInStand[1], # m²/ha
                                      alruBA = sum((speciesGroup == "RA") * treeBasalAreaPerHectare) / plotsInStand[1],
                                      tsheBA = sum((speciesGroup == "WH") * treeBasalAreaPerHectare) / plotsInStand[1],
                                      acmaBA = sum((speciesGroup == "BM") * treeBasalAreaPerHectare) / plotsInStand[1],
                                      umcaBA = sum((speciesGroup == "OM") * treeBasalAreaPerHectare) / plotsInStand[1],
                                      thplBA = sum((speciesGroup == "RC") * treeBasalAreaPerHectare) / plotsInStand[1],
                                      otherBA = sum((speciesGroup == "other") * treeBasalAreaPerHectare) / plotsInStand[1]),
                          by = c("StandID")) %>%
  mutate(isCruised = (is.na(tph) == FALSE))

strata2016 = cruiseTph2016 %>% group_by(vegLabel) %>% 
  summarize(stands = n(), area = sum(standArea), 
            cruisedStands = sum(isCruised), cruisedArea = sum(standArea * isCruised), 
            cruisedTph = sum(standArea * tph, na.rm = TRUE) / sum(standArea * (is.na(tph) == FALSE)),
            meanAge2016 = sum(standArea * standAge2016) / sum(standArea), 
            meanTopHeight = sum(standArea * topHeight, na.rm = TRUE) / sum(standArea * (is.na(topHeight) == FALSE)),
            meanPsmeBA = sum(standArea * psmeBA, na.rm = TRUE) / sum(standArea * (is.na(psmeBA) == FALSE)),
            meanAlruBA = sum(standArea * alruBA, na.rm = TRUE) / sum(standArea * (is.na(alruBA) == FALSE)),
            meanTsheBA = sum(standArea * tsheBA, na.rm = TRUE) / sum(standArea * (is.na(tsheBA) == FALSE)),
            meanAcmaBA = sum(standArea * acmaBA, na.rm = TRUE) / sum(standArea * (is.na(acmaBA) == FALSE)),
            meanUmcaBA = sum(standArea * umcaBA, na.rm = TRUE) / sum(standArea * (is.na(umcaBA) == FALSE)),
            meanThplBA = sum(standArea * thplBA, na.rm = TRUE) / sum(standArea * (is.na(thplBA) == FALSE)),
            meanOtherBA = sum(standArea * otherBA, na.rm = TRUE) / sum(standArea * (is.na(otherBA) == FALSE)),
            plantationArea = sum(isPlantation * standArea)) %>%
  mutate(cruisedTrees = cruisedArea * cruisedTph, estimatedTrees = area * cruisedTph, vegCode = str_trunc(vegLabel, 2, ellipsis = ""), 
         standGroup = if_else(stands > 35, vegLabel, if_else(vegCode %in% c("DX", "1D"), "DX34L", "other")))
#strata2016 %>% summarize(groups = length(unique(standGroup)), vegLabel = "total", forestArea = sum(area), uncruisedArea = sum(is.na(tph) * area), tph = sum(area * tph, na.rm = TRUE) / sum((is.na(area) == FALSE) * area), treesInStrata = sum(trees, na.rm = TRUE), trees = 1 / (1 - uncruisedArea / forestArea) * treesInStrata)

strata2016grouped = strata2016 %>% group_by(standGroup) %>%
  summarize(stands = sum(stands), 
            cruisedStands = sum(cruisedStands), 
            cruisedTrees = sum(cruisedArea * cruisedTrees, na.rm = TRUE) / sum(cruisedArea), # some strata have zero cruisedArea => NA for top height, TPH, BA
            meanAge2016 = sum(area * meanAge2016) / sum(area), # no NA ages
            meanTopHeight = sum(cruisedArea * meanTopHeight, na.rm = TRUE) / sum(cruisedArea),
            meanPsmeBA = sum(cruisedArea * meanPsmeBA, na.rm = TRUE) / sum(cruisedArea),
            meanAlruBA = sum(cruisedArea * meanAlruBA, na.rm = TRUE) / sum(cruisedArea),
            meanTsheBA = sum(cruisedArea * meanTsheBA, na.rm = TRUE) / sum(cruisedArea),
            meanAcmaBA = sum(cruisedArea * meanAcmaBA, na.rm = TRUE) / sum(cruisedArea),
            meanUmcaBA = sum(cruisedArea * meanUmcaBA, na.rm = TRUE) / sum(cruisedArea),
            meanThplBA = sum(cruisedArea * meanThplBA, na.rm = TRUE) / sum(cruisedArea),
            meanOtherBA = sum(cruisedArea * meanOtherBA, na.rm = TRUE) / sum(cruisedArea),
            area = sum(area), 
            cruisedArea = sum(cruisedArea),
            plantationArea = sum(plantationArea)) %>% 
  mutate(cruisedTph = cruisedTrees / cruisedArea, 
         estimatedTrees = area * cruisedTph,
         meanTotalBA = meanPsmeBA + meanAlruBA + meanTsheBA + meanAcmaBA + meanUmcaBA + meanThplBA + meanOtherBA,
         standGroup = if_else(standGroup == "other", "all others", standGroup),
         vegCode = str_trunc(standGroup, 2, ellipsis = "")) %>%
  relocate(standGroup, stands, meanAge2016, estimatedTrees, area, cruisedArea, cruisedStands, cruisedTrees)


## Dooley and Fairweather Table 7
# ODF/ODSL vegetation codes = species + size + density
# species: 1 = single species, D = Douglas-fir, W = hemlock, H = hardwood, X = mixture, OT = other
# size: 1 = DBH < 20.3 cm, 2 = 35.6 cm, 3 = DBH < 50.8 cm, 4 = DBH < 76.2 cm, 5 = DBH ≥ 76.2 cm
# density: L = SDI56% < 30, H = SDI56% ≥ 30, SDI65 = stand density index of trees ≥ 14.2 cm DBH
ggplot() +
  geom_col(aes(x = stands, y = fct_reorder(vegLabel, area), alpha = "forest", fill = vegCode), strata2016) +
  geom_col(aes(x = cruisedStands, y = fct_reorder(vegLabel, area), alpha = "cruised", fill = vegCode), strata2016) +
  labs(x = "stands", y = NULL, alpha = NULL, fill = NULL) +
ggplot() +
  geom_col(aes(x = area, y = fct_reorder(vegLabel, area), alpha = "forest", fill = vegCode), strata2016) +
  geom_col(aes(x = cruisedArea, y = fct_reorder(vegLabel, area), alpha = "cruised", fill = vegCode), strata2016) +
  labs(x = "hectares", y = NULL, alpha = NULL, fill = NULL) +
  guides(alpha = "none", fill = "none") +
ggplot() +
  geom_col(aes(x = cruisedTph, y = fct_reorder(vegLabel, area), alpha = "cruised", fill = vegCode), strata2016) +
  labs(x = "trees per hectare", y = NULL, alpha = NULL, fill = NULL) +
  guides(alpha = "none", fill = "none") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_alpha_manual(breaks = c("forest", "cruised"), values = c(0.5, 1)) &
  scale_fill_manual(breaks = c("DX", "1D", "HX", "OT", "WX", "1H", "XC", "1W", "XG", "XS", "XW"), values = c("forestgreen", "forestgreen", "red", "grey50", "blue", "red", "brown4", "blue", "brown4", "brown4", "brown4"))

ggplot() +
  geom_col(aes(x = stands, y = fct_reorder(standGroup, area), alpha = "forest", fill = vegCode), strata2016grouped) +
  geom_col(aes(x = cruisedStands, y = fct_reorder(standGroup, area), alpha = "cruised", fill = vegCode), strata2016grouped) +
  labs(x = "stands", y = NULL, alpha = NULL, fill = NULL) +
ggplot() +
  geom_col(aes(x = area, y = fct_reorder(standGroup, area), alpha = "forest", fill = vegCode), strata2016grouped) +
  geom_col(aes(x = cruisedArea, y = fct_reorder(standGroup, area), alpha = "cruised", fill = vegCode), strata2016grouped) +
  labs(x = "hectares", y = NULL, alpha = NULL, fill = NULL) +
  guides(alpha = "none", fill = "none") +
ggplot() +
  geom_col(aes(x = cruisedTph, y = fct_reorder(standGroup, area), alpha = "cruised", fill = vegCode), strata2016grouped) +
  labs(x = "trees per hectare", y = NULL, alpha = NULL, fill = NULL) +
  guides(alpha = "none", fill = "none") +
plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
  scale_alpha_manual(breaks = c("forest", "cruised"), values = c(0.5, 1)) &
  scale_fill_manual(breaks = c("DX", "1D", "HX", "ot", "WX", "1H", "XC", "1W", "XG", "XS", "XW"), values = c("forestgreen", "forestgreen", "red", "grey50", "blue", "red", "brown4", "blue", "brown4", "brown4", "brown4"))


## Figure 2: Elliott data descriptor
strata2016annotations = strata2016grouped %>% 
  mutate(ageLabel = if_else(area == max(area), sprintf('bar(age) == %.0f~"years"', meanAge2016), if_else(area > 4000, sprintf("%.0f~years", meanAge2016), sprintf("%.0f", meanAge2016))),
         ageLabelX = if_else(area > 4200, 3280, area + 100),
         ageLabelColor = if_else(area > 4200, "white", "black"),
         topHeightLabel = if_else(area == max(area), sprintf('bar(H)[100] == "%.1f"~"m"', meanTopHeight), if_else(area > 4000, sprintf('"%.1f m"', meanTopHeight), sprintf('"%.1f"', meanTopHeight))),
         topHeightLabelColor = if_else(meanTotalBA > 10, "white", "black"),
         topHeightLabelX = if_else(meanTotalBA > 10, 1, meanTotalBA + 1))

ggplot() +
  geom_col(aes(x = area, y = fct_reorder(standGroup, area), fill = isPlantation, group = fct_rev(isPlantation)), orientation = "y", strata2016grouped %>% mutate(naturalRegenArea = area - plantationArea) %>% 
             select(standGroup, naturalRegenArea, plantationArea) %>%
             pivot_longer(cols = c("naturalRegenArea", "plantationArea"), names_to = "isPlantation", values_to = "area") %>%
             mutate(isPlantation = factor(isPlantation, labels = c("natural regeneration", "plantation"), levels = c("naturalRegenArea", "plantationArea")))) +
  geom_col(aes(x = cruisedArea, y = as.numeric(fct_reorder(standGroup, area)) - 0.5 * (0.9 - 0.35), fill = "area of cruised stands"), strata2016grouped, width = 0.35) +
  #geom_segment(aes(x = cruisedArea, xend = cruisedArea, y = as.numeric(fct_reorder(standGroup, area)) - 0.45, yend = as.numeric(fct_reorder(standGroup, area)) + 0.45, color = "area of cruised stands"), strata2016grouped, linewidth = 0.9) +
  geom_text(aes(x = ageLabelX, y = fct_reorder(standGroup, area), label = ageLabel), strata2016annotations, color = strata2016annotations$ageLabelColor, hjust = 0, parse = TRUE, size = 3) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  labs(x = "area, ha", y = "cruise strata", color = NULL, fill = NULL, title = "a) area and mean age") +
  #scale_color_manual(breaks = c("area of cruised stands"), values = c("cyan2")) +
  #scale_fill_manual(breaks = c("natural regeneration", "plantation"), values = c("grey10", "grey35")) +
  scale_fill_manual(breaks = c("natural regeneration", "plantation", "area of cruised stands"), values = c("grey10", "grey35", "purple2")) +
  scale_x_continuous(expand = c(0.012, 0)) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.02)) +
ggplot() +
  geom_bar(aes(y = fct_reorder(standGroup, area), fill = fct_rev(speciesGroup), weight = basalArea), strata2016grouped %>% select(standGroup, area, meanPsmeBA, meanAlruBA, meanTsheBA, meanAcmaBA, meanUmcaBA, meanThplBA, meanOtherBA) %>%
             pivot_longer(cols = c("meanPsmeBA", "meanAlruBA", "meanTsheBA", "meanAcmaBA", "meanUmcaBA", "meanThplBA", "meanOtherBA"), names_to = "speciesGroup", values_to = "basalArea") %>%
             mutate(speciesGroup = factor(speciesGroup, levels = c("meanPsmeBA", "meanAlruBA", "meanTsheBA", "meanAcmaBA", "meanUmcaBA", "meanThplBA", "meanOtherBA"))), na.rm = TRUE) +
  geom_text(aes(x = topHeightLabelX, y = fct_reorder(standGroup, area), label = topHeightLabel), strata2016annotations, color = strata2016annotations$topHeightLabelColor, hjust = 0, parse = TRUE, size = 3) +
  coord_cartesian(xlim = c(0, 80)) +
  labs(x = bquote("mean basal area, m"^2*" ha"^-1), y = NULL, fill = NULL, title = "b) density and mean top height") +
  scale_fill_manual(breaks = c("meanPsmeBA", "meanAlruBA", "meanTsheBA", "meanAcmaBA", "meanUmcaBA", "meanThplBA", "meanOtherBA"), labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
  scale_x_continuous(expand = c(0.008, 0)) +
  scale_y_discrete(labels = NULL) +
  theme(legend.key.height = unit(1, "line"), legend.key.width = unit(1, "line")) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, ncol = 2, widths = c(0.4, 0.6))
#ggsave("trees/height-diameter/figures/Figure 02 Elliott cruise strata.png", height = 9.5, width = 22, units = "cm", dpi = 250)


# stand age distribution in 2016
ggplot() +
  geom_histogram(aes(x = 2016 - standAge2016, fill = isPlantation, weight = area_ha), stands2022, binwidth = 1) +
  labs(x = "stand origin year", y = "area, ha", fill = NULL)
