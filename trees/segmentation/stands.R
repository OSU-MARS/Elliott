library(dplyr)
library(patchwork)
library(terra)
library(tibble)

# 2009 LiDAR tiles: 3324x4544 feet, EPSG:2994
# 2021 LiDAR tiles: 3000x3000 feet, from 384000,699000 EPSG:6557
stands2022gis = left_join(as_tibble(vect("GIS/Planning/Elliott Stand Data Feb2022 EPSG6556.gpkg")) %>% rename(standID2016 = StandID),
                          stands2022 %>% mutate(plotsInStand = replace_na(plotsInStand, 0)),
                          by = c("standID2016")) %>%
  mutate(strata = paste(if_else(plotsInStand > 0, "cruised", "unsampled"), if_else(isPlantation, "plantation", "natural regen")))

unique(stands2022gis$strata)

ggplot(stands2022gis) +
  geom_histogram(aes(x = standAge2016, fill = strata, weight = area_ha), binwidth = 2) +
  labs(x = "stand age, years", y = "area, ha", fill = NULL) +
ggplot(stands2022gis) +
  geom_histogram(aes(x = Elev_Mean, fill = strata, weight = area_ha), binwidth = 50) +
  labs(x = "mean elevation, m", y = "area, ha", fill = NULL) +
ggplot(stands2022gis) +
  geom_histogram(aes(x = SlopeMean, fill = strata, weight = area_ha), binwidth = 1) +
  labs(x = "mean slope, °", y = "area, ha", fill = NULL) +
ggplot(stands2022gis) +
  geom_histogram(aes(x = PrecipNorm, fill = strata, weight = area_ha), binwidth = 25) +
  labs(x = bquote("mean precipitation, mm year"^-1), y = "area, ha", fill = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_fill_manual(breaks = c("unsampled natural regen", "unsampled plantation", "cruised natural regen", "cruised plantation"), values = c("forestgreen", "darkviolet", "green2", "blue3"))

# effective plot sizes
plotRadii = bind_rows(tibble(method = "variable radius", sampleFactor = c(10, 13.61, 17.78, 20, 22.5, 25.15, 27.78, 32.14, 33.61, 40, 46.94, 52.57, 54.44, 62.5, 71.11, 80.27, 90)),
                      tibble(method = "fixed radius", sampleFactor = 200)) %>%
  mutate(limitingDistance = if_else(method == "fixed radius", sqrt(43560 / (pi * sampleFactor)), 100 / 2.54 * 1 / (3.28084 * 12 * sqrt(sampleFactor / 10890)))) # m, BAF = 10890 * (D/R)² => D/R = sqrt(BAF / 10890) => R = D / sqrt(BAF / 10890) inches

plotSizes = trees2016 %>% filter(PlotType == "IP") %>% # no DBH measurements on count plots
  group_by(PlotID) %>%
  summarize(isPlantation = isPlantation[1], minRadius = min(plotRadius), meanRadius = mean(plotRadius), maxRadius = max(plotRadius), nominalRadius = 2/3 * meanRadius) # 2/3 since average distance to tree in circular plot is a radius

ggplot() +
  geom_hline(yintercept = 2.54, color = "grey70", linetype = "longdash") +
  geom_violin(aes(x = isPlantation, y = plotRadius, color = isPlantation), trees2016 %>% filter(PlotType == "IP"), draw_quantiles = c(0.25, 0.50, 0.75)) +
  labs(x = NULL, y = "plot radius for tree, m") +
ggplot() +
  geom_hline(yintercept = c(2.54, 11.284), color = "grey70", linetype = "longdash") +
  geom_violin(aes(x = isPlantation, y = nominalRadius, color = isPlantation), plotSizes, draw_quantiles = c(0.25, 0.50, 0.75)) +
  labs(x = NULL, y = "nominal radius of plot, m") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout() &
  coord_cartesian(ylim = c(0, 33)) &
  guides(color = "none") &
  scale_x_discrete(breaks = c(FALSE, TRUE), labels = c("natural\nregen", "plantation"))
ggsave("trees/segmentation/figures/plot sizes.png", width = 15, height = 10, units = "cm", dpi = 120)
