library(dplyr)
library(ggplot2)
library(terra)
library(readxl)
library(stringr)

# trees2016 and stands2022 from trees/height-diameter/setup.R
cruiseTph2016 = left_join(trees2016 %>% group_by(StandID) %>% summarize(tph = tph[1], standArea = standArea[1]),
                          stands2022 %>% select(StandID, ODSL_VEG_L) %>% rename(vegLabel = ODSL_VEG_L),
                          by = c("StandID"))

cruiseStrata2016 = cruiseTph2016 %>% group_by(vegLabel) %>% summarize(area = sum(standArea), tph = sum(tph * standArea) / area, trees = area * tph)

strata2016 = left_join(stands2022 %>% group_by(ODSL_VEG_L) %>% summarize(area = sum(area_ha)) %>% rename(vegLabel = ODSL_VEG_L),
                       cruiseStrata2016 %>% select(vegLabel, tph),
                       by = c("vegLabel")) %>%
  mutate(trees = area * tph, vegCode = str_trunc(vegLabel, 2, ellipsis = ""), size = str_ex)

strata2016 %>% summarize(vegLabel = "total", forestArea = sum(area), uncruisedArea = sum(is.na(tph) * area), tph = sum(area * tph, na.rm = TRUE) / sum((is.na(area) == FALSE) * area), treesInStrata = sum(trees, na.rm = TRUE), trees = 1 / (1 - uncruisedArea / forestArea) * treesInStrata)

# Dooley and Fairweather Table 7
# ODF/ODSL vegetation codes = species + size + density
# species: 1 = single species, D = Douglas-fir, W = hemlock, H = hardwood, X = mixture, OT = other
# size: 1 = DBH < 20.3 cm, 2 = 35.6 cm, 3 = DBH < 50.8 cm, 4 = DBH < 76.2 cm, 5 = DBH ≥ 76.2 cm
# density: L = SDI56% < 30, H = SDI56% ≥ 30, SDI65 = stand density index of trees ≥ 14.2 cm DBH
ggplot() +
  geom_col(aes(x = area, y = fct_reorder(vegLabel, area), fill = vegCode), strata2016) +
  labs(x = "hectares", y = NULL, fill = NULL) +
ggplot() +
  geom_col(aes(x = tph, y = fct_reorder(vegLabel, area), fill = vegCode), strata2016) +
  labs(x = "trees per hectare", y = NULL, fill = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_fill_manual(breaks = c("DX", "1D", "HX", "OT", "WX", "1H", "XC", "1W", "XG", "XS", "XW"), values = c("forestgreen", "forestgreen", "red", "grey50", "blue", "red", "brown4", "blue", "brown4", "brown4", "brown4"))
