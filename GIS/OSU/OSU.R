library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(sf)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             axis.title = element_text(size = 9),
                             legend.title = element_text(size = 9),
                             panel.border = element_blank(),
                             plot.title = element_text(size = 9)))

bioTrees2023 = read_xlsx("GIS/OSU/ESRF 2023 Biodiversity Surveys - Field Datasheets - 7_CanopyPlots.xlsx", sheet = "ESRF 2023 Biodiversity Surveys") %>%
  mutate(`Alive?` = replace_na(`Alive?`, "no")) # one snag missing entry
bioPlots2023gis = st_read("GIS/OSU/ESRF 2023 Biodiversity Surveys 2.gpkg", quiet = TRUE)
stands2022 = st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", layer = "unified stands 2022 property boundary split", quiet = TRUE) %>%
  filter(isExternalBoundarySplit == 0) %>% group_by(standID2016) %>%
  summarize(standArea = sum(standArea), standAge2016 = standAge2016[1])
#bioPlots2025 = read_xlsx("GIS/OSU/2025 plot data.xlsx", sheet = "Scanned & Cruised BD plots")
#st_write(left_join(bioPlots2023gis, bioPlots2025 %>% select(StationID) %>% mutate(scanYear = 2025), 
#                   by = join_by(StationID)), 
#         "GIS/OSU/ESRF 2023 Biodiversity Surveys 2.gpkg", layer = "esrf2023_biodiversity_surveys", append = FALSE)

bioPlots2023 = left_join(left_join(bioTrees2023 %>% group_by(StationID) %>% 
                                     summarize(trees = n(), liveTrees = sum(`Alive?` == "yes"), douglasFir = sum(Species == "PSME"), liveBasalArea = sum(pi/4 * (0.01 * DBH * (`Alive?` == "yes"))^2), snags = sum(`Alive?` == "no")),
                                   st_drop_geometry(bioPlots2023gis) %>% select(StationID, standID2016, elevation, slope, aspect, scanYear), 
                                   by = join_by(StationID)),
                         stands2022 %>% select(standID2016, standArea, standAge2016), 
                         by = join_by(standID2016)) %>%
  mutate(subset = if_else(is.na(scanYear), if_else(douglasFir / trees >= 0.7, "2023 only", "<70% Douglas-fir"), "2025 scans"))

diversePlots = bioPlots2023 %>% filter(douglasFir / trees < 0.7)
scanPlots2025 = bioPlots2023 %>% filter(StationID %in% (bioPlots2023gis %>% filter(scanYear == 2025))$StationID)

as.data.frame(table(bioTrees2023$Species, dnn = "species", useNA = "ifany"), responseName = "records") %>% 
  mutate(percentage = 100 * records / sum(records)) %>%
  arrange(desc(records))

as.data.frame(table((bioTrees2023 %>% filter(StationID %in% diversePlots$StationID))$Species, dnn = "species", useNA = "ifany"), responseName = "records") %>% 
  mutate(percentage = 100 * records / sum(records)) %>%
  arrange(desc(records))

as.data.frame(table((bioTrees2023 %>% filter(StationID %in% scanPlots2025$StationID))$Species, dnn = "species", useNA = "ifany"), responseName = "records") %>% 
  mutate(percentage = 100 * records / sum(records)) %>%
  arrange(desc(records))

bioTrees2023 %>% filter(StationID %in% diversePlots$StationID) %>% group_by(Species, CanopyPosition) %>% summarize(n = n()) %>%
  pivot_wider(names_from = "CanopyPosition", values_from = "n") %>% mutate(notRecorded = `NA` + `N/A`) %>% select(-`NA`, -`N/A`) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  rename(species = Species, outside = O, dominant = D, codominant = C, intermediate = I, suppressed = S) %>%
  relocate(species, outside, dominant, codominant, intermediate, suppressed, notRecorded) %>%
  arrange(desc(dominant + codominant))

print(scanPlots2025 %>% select(StationID) %>% arrange(StationID), n = 50)

ggplot() +
  geom_point(aes(x = liveTrees, y = liveBasalArea, color = as.factor(snags)), bioPlots2023, alpha = 0.5, shape = 16) +
  guides(color = guide_legend(override.aes = list(alpha = 0.8))) +
  labs(x = "live trees", y = bquote("live basal area on plot, m"^2), color = "snags") +
  scale_color_viridis_d()

ggplot() +
  geom_histogram(aes(x = 100 * douglasFir / trees), bioPlots2023, binwidth = 5) +
  labs(x = "Douglas-fir, % stems", y = "plots")

# elevation, slope, and aspect histograms
# Elevation slow to load due to large .vrt.
elliott2022 = terra::vect("GIS/Planning/ESRF boundary 2022-04.gpkg", layer = "ESRF boundary with Hakki 2022-04")
elevation2021 = terra::crop(terra::rast("GIS/DOGAMI/2021 OLC Coos County/DTM/DTM.vrt"), elliott2022)
slope2021 = terra::crop(terra::rast("GIS/DOGAMI/2021 OLC Coos County/bare earth slope Gaussian 10 m EPSG6557.tif"), elliott2022)
aspect2021 = 180 / pi * terra::atan2(terra::crop(terra::rast("GIS/DOGAMI/2021 OLC Coos County/bare earth sin(aspect) Gaussian 10 m EPSG6557.tif"), elliott2022),
                                     terra::crop(terra::rast("GIS/DOGAMI/2021 OLC Coos County/bare earth cos(aspect) Gaussian 10 m EPSG6557.tif"), elliott2022))
aspect2021 = terra::ifel(aspect2021 < 0, 360 + aspect2021, aspect2021)

elevation2021histogram = terra::hist(elevation2021, breaks = seq(0, 13 * 3.28084 * 50, by = 3.28084 * 50), maxcell = 2.147E9, plot = FALSE) # constrained to 66% sample, 2^31 - 1 seems to be the largest maxsize without early data is too long errors but nearby values eventually fail with long vectors not supported yet: ../include/Rinlinedfuns.h:551
slope2021histogram = terra::hist(slope2021, breaks = seq(0, 90), maxcell = prod(dim(slope2021)), plot = FALSE) # specify all cells as terra 1.8-86's sampling draws from all cells, not just ones with data
aspect2021histogram = terra::hist(aspect2021, breaks = seq(0, 360), maxcell = prod(dim(aspect2021)), plot = FALSE)

ggplot() +
  geom_histogram(aes(x = 0.3048 * elevation, fill = subset, group = subset), bioPlots2023, binwidth = 50) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 75)) +
  labs(x = NULL, y = "biodiversity plots", fill = "plot subset") +
  scale_x_continuous(breaks = seq(0, 500, by = 250)) +
ggplot() +
  geom_histogram(aes(x = slope, fill = subset, group = subset), bioPlots2023, binwidth = 5) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 75)) +
  labs(x = NULL, y = NULL, fill = "plot subset") +
  scale_x_continuous(breaks = seq(0, 60, by = 20)) +
ggplot() +
  geom_histogram(aes(x = aspect, fill = subset, group = subset), bioPlots2023, binwidth = 30) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 75)) +
  labs(x = NULL, y = NULL, fill = "plot subset") +
  scale_x_continuous(breaks = seq(0, 360, by = 180)) +
ggplot() +
  geom_histogram(aes(x = standAge2016 + 4, fill = subset, group = subset), bioPlots2023, binwidth = 10) +
  coord_cartesian(xlim = c(0, 200), ylim = c(0, 75)) +
  labs(x = NULL, y = NULL, fill = "plot subset") +
  scale_x_continuous(breaks = seq(0, 200, by = 100)) +
ggplot() +
  geom_histogram(aes(x = mid, weight = density), tibble(mid = 0.3028 * elevation2021histogram$mids, density = elevation2021histogram$counts / sum(elevation2021histogram$counts)), binwidth = 50, fill = "grey80") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 0.31)) +
  labs(x = "elevation, m", y = "P(Elliott)") +
  scale_x_continuous(breaks = seq(0, 500, by = 250)) +
ggplot() +
  geom_histogram(aes(x = mid, weight = density), tibble(mid = slope2021histogram$mids, density = slope2021histogram$density) %>% filter(density > 0), binwidth = 5, fill = "grey80") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 0.31)) +
  labs(x = "slope, °", y = NULL) +
  scale_x_continuous(breaks = seq(0, 60, by = 20)) +
ggplot() +
  geom_histogram(aes(x = mid, weight = density), tibble(mid = aspect2021histogram$mids, density = aspect2021histogram$density), binwidth = 30, fill = "grey80") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 0.31)) +
  labs(x = "aspect, °", y = NULL) +
  scale_x_continuous(breaks = seq(0, 360, by = 180)) +
ggplot() +
  geom_histogram(aes(x = standAge2016 + 4, y = after_stat(count / sum(count)), weight = standArea), stands2022, binwidth = 10, fill = "grey80") + 
  coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.31)) +
  labs(x = "stand age, years", y = NULL) +
  scale_x_continuous(breaks = seq(0, 200, by = 100)) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 4, heights = c(0.7, 0.3), guides = "collect") &
  scale_fill_brewer() &
  scale_y_continuous(expand = c(0, 0))
ggsave("GIS/OSU/ESRF biodiversity plot distributions.png", width = 16.5, height = 10, dpi = 200, units = "cm")
