library(arrow)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)
library(stringr)
library(tidyr)

future::plan(future::multisession, workers = 6) # multicore not useful as future 1.34 considers it unstable and silently overrides to workers = 1 (future's docs say there's a warning but this is incorrect)
theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             panel.border = element_blank(),
                             plot.title = element_text(size = 10)))

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileStatOptions = tibble(getLayerCounts = FALSE,
                         tileArea = 3000^2 / (2.47105 * 43550)) # ha

elliottStands = readxl::read_xlsx("GIS/Planning/Elliott Stand Data Feb2022.xlsx", sheet = "Elliott Stand Data Feb2022") %>%
  mutate(standAge2016 = pmax(if_else((Age_2020 - 4) > (Age_2015 + 1), Age_2015 + 1, Age_2020 - 4), 0),
         standArea = 0.404686 * GrossAc,  # ac to ha
         isPlantation = standAge2016 < 70)
organonTreePredictions = read_feather("trees/Organon/Elliott tree lists 2016-2116.feather")
tileIndex = st_drop_geometry(st_read("GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg", as_tibble = TRUE, layer = "Elliott tile index", quiet = TRUE))


if (tileStatOptions$getLayerCounts)
{
  localMaximaFiles = list.files(file.path(dataPath, "DSM v3 beta/local maxima"), ".gpkg", full.names = TRUE)
  localMaximaStartTime = Sys.time() # 2.5 minutes for 561 layers, 9900X with six workers (14 minutes with purrr::map)
  localMaxima = bind_rows(furrr::future_map(localMaximaFiles, function(localMaximaFile) { 
      localMaximaTile = st_read(localMaximaFile, as_tibble = TRUE, quiet = TRUE) # assume layer = "localMaxima"
      return(tibble(tile = str_remove(basename(localMaximaFile), ".gpkg"), localMaxima = nrow(localMaximaTile)))
    })) %>% arrange(tile)
  localMaximaTime = Sys.time() - localMaximaStartTime
  saveRDS(localMaxima, "trees/segmentation/localMaximaTile summary.Rds")

  treetopsRandomForestFiles = list.files(file.path(dataPath, "treetops/rf"), ".gpkg", full.names = TRUE)
  treetopsStartTime = Sys.time() #  22 s for 561 layers, 9900X
  treetops = bind_rows(furrr::future_map(treetopsRandomForestFiles, function(treetopsRandomForestFile) { 
    treetopsTile = st_drop_geometry(st_read(treetopsRandomForestFile, as_tibble = TRUE, layer = "treetops", quiet = TRUE)) #%>%
      #filter(height < 100000) # work around DTM no data bug in CHM generation
    mergePointsTile = st_read(treetopsRandomForestFile, as_tibble = TRUE, layer = "merge points", quiet = TRUE)
    noiseTile = st_read(treetopsRandomForestFile, as_tibble = TRUE, layer = "noise points", quiet = TRUE)
    maybeNoiseTile = st_read(treetopsRandomForestFile, as_tibble = TRUE, layer = "maybe noise points", quiet = TRUE)
    
    heightHistogram = left_join(tibble(heightClass = seq(0, 100)),
                                treetopsTile %>% mutate(heightClass = round(0.3048 * height)) %>% group_by(heightClass) %>% summarize(n = n()),
                                by = join_by(heightClass)) %>%
      mutate(n = replace_na(n, 0)) %>%
      arrange(heightClass) %>% 
      pivot_wider(names_prefix = "height", names_from = "heightClass", values_from = "n")
    return(heightHistogram %>% 
             mutate(tile = str_remove(basename(treetopsRandomForestFile), ".gpkg"), 
                    treetops = nrow(treetopsTile), 
                    mergePoints = nrow(mergePointsTile), 
                    noise = nrow(noiseTile), 
                    maybeNoise = nrow(maybeNoiseTile),
                    maxHeight = max(treetopsTile$height)) %>%
             relocate(tile, treetops, mergePoints, noise, maybeNoise))
  }, .options = furrr::furrr_options(seed = TRUE))) %>% arrange(tile)
  treetopsTime = Sys.time() - treetopsStartTime
  
  # tiles affected CHM height bug
  treetops %>% filter(maxHeight > 400) %>% select(tile, maxHeight)
} else {
  localMaxima = readRDS("trees/segmentation/localMaximaTile summary.Rds")
}

tileCounts = left_join(left_join(localMaxima, treetops, by = join_by(tile)),
                       tileIndex %>% select(Tile_ID, bufferDistance, treetopQuad), by = join_by(tile == Tile_ID)) %>%
  mutate(bufferDistance = as.factor(bufferDistance))
heightHistogram = tileCounts %>% select(tile, bufferDistance, starts_with("height")) %>%
  pivot_longer(cols = starts_with("height"), names_prefix = "height", names_to = "heightClass", values_to = "trees") %>%
  mutate(heightClass = as.numeric(heightClass),
         category = if_else(heightClass > 7, "tree", if_else(heightClass > 4, "treeOrTallShrub", "maybeTree")),
         tileArea = if_else(bufferDistance == 0, 510 * tileStatOptions$tileArea, 51 * tileStatOptions$tileArea),
         treesPerHectare = trees / tileArea)
heightHistogram %>% group_by(bufferDistance, category) %>% summarise(trees = sum(trees)) %>% mutate(treesElliott = trees * 34726 / (510 * tileStatOptions$tileArea))

tileCounts %>% select(-starts_with("height"), -treetopQuad) %>% group_by(bufferDistance) %>% summarize(across(-tile, sum)) %>% 
  mutate(tileArea = if_else(bufferDistance == 0, 510 * tileStatOptions$tileArea, 51 * tileStatOptions$tileArea),
         treetopsPerHectare = treetops / tileArea)
tileCounts %>% filter(tile %in% c("s04200w06840", "s04200w06810", "s04230w06810"))

ggplot() +
  geom_histogram(aes(x = localMaxima, fill = bufferDistance, group = bufferDistance), tileCounts, binwidth = 10000) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 200)) +
  labs(x = "local maxima", y = "tiles", fill = "tiles") +
  scale_x_continuous(labels = scales::comma) +
ggplot() +
  geom_histogram(aes(x = treetops, fill = bufferDistance, group = bufferDistance), tileCounts, binwidth = 1000) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 200)) +
  labs(x = "treetops", y = NULL, fill = "tiles") +
  scale_x_continuous(labels = scales::comma) +
ggplot() +
  geom_histogram(aes(x = mergePoints, fill = bufferDistance, group = bufferDistance), tileCounts, binwidth = 50) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 200)) +
  labs(x = "merge points", y = "tiles", fill = "tiles") +
ggplot() +
  geom_histogram(aes(x = noise, fill = bufferDistance, group = bufferDistance), tileCounts, binwidth = 1) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 200)) +
  labs(x = "noise points", y = NULL, fill = "tiles") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_fill_manual(breaks = c(0, 400), labels = c("Elliott", "adjacent"), values = c("forestgreen", "grey30"))

cruiseHeightDistribution2021 = left_join(organonTreePredictions %>% filter(year == 2021) %>% mutate(heightClass = round(height)) %>% group_by(stand, heightClass) %>%
                                           summarize(treesPerHectare = sum(liveExpansionFactor) / length(unique(plot)),
                                                     snagsPerHectare = sum(deadExpansionFactor) / length(unique(plot)), .groups = "drop"),
                                         elliottStands %>% select(StandID, isPlantation, standArea), by = join_by(stand == StandID))
forestHeightDistribution = cruiseHeightDistribution2021 %>% group_by(isPlantation, heightClass) %>%
  summarize(treesPerHectare = sum(standArea * treesPerHectare) / sum(standArea),
            snagsPerHectare = sum(standArea * snagsPerHectare) / sum(standArea), .groups = "drop")
forestHeightDistribution %>% group_by(isPlantation) %>% summarize(treesPerHectare = sum(treesPerHectare), snagsPerHectare = sum(snagsPerHectare))

ggplot() +
  geom_col(aes(x = heightClass, y = trees, fill = bufferDistance, group = bufferDistance), heightHistogram) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 1.2E6)) +
  labs(x = "height class, m", y = "treetop candidates", fill = "tiles", title = "a) LiDAR detected trees") +
  scale_y_continuous(breaks = seq(0, 1.2E6, by = 3E5), labels = scales::comma) +
ggplot() +
  geom_col(aes(x = heightClass, y = treesPerHectare, fill = bufferDistance, group = bufferDistance), heightHistogram %>% filter(bufferDistance == 0)) +
  geom_line(aes(x = heightClass, y = treesPerHectare + snagsPerHectare, color = isPlantation, group = isPlantation), forestHeightDistribution) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 100)) +
  guides(fill = "none") +
  labs(x = "height class, m", y = "stems per hectare", color = "cruise projection", fill = "tiles", title = "b) Elliott density") +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
inset_element(ggplot() + 
                geom_col(aes(x = heightClass, y = treesPerHectare, fill = bufferDistance, group = bufferDistance), heightHistogram %>% filter(bufferDistance == 0)) +
                geom_line(aes(x = heightClass, y = treesPerHectare + snagsPerHectare, color = isPlantation, group = isPlantation), forestHeightDistribution) +
                coord_cartesian(xlim = c(25, 100), ylim = c(0, 8)) +
                guides(color = "none", fill = "none") +
                labs(x = NULL, y = NULL) + 
                scale_x_continuous(breaks = seq(25, 100, by = 25), expand = c(0, NA)) +
                scale_y_continuous(breaks = seq(0, 8, by = 2)) +
                theme(plot.margin = margin()),
              0.33, 0.48, 0.96, 0.98, align_to = "plot") +
ggplot() +
  geom_col(aes(x = heightClass, y = treesPerHectare, fill = bufferDistance, group = bufferDistance), heightHistogram %>% filter(bufferDistance != 0)) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 100)) +
  guides(fill = "none") +
  labs(x = "height class, m", y = "stems per hectare", fill = "tiles", title = "c) adjacent density") +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_color_manual(breaks = c(FALSE, TRUE), labels = c("natural regen", "plantation"), values = c("forestgreen", "blue")) &
  scale_fill_manual(breaks = c(0, 400), labels = c("Elliott", "adjacent"), values = c("forestgreen", "grey30"))

