# assumes segmentation/setup.R
library(sf)

diff_treetops = function(treetops1, treetops2)
{
  return(tibble(treeID = treetops1$treeID != treetops2$treeID,
                x = treetops1$x != treetops2$x,
                y = treetops1$y != treetops2$y,
                z = treetops1$z != treetops2$z,
                height = treetops1$height != treetops2$height))
}

read_treetops = function(geoPackagePath)
{
  geoPackage = read_sf(geoPackagePath) # terra::vect() doesn't support z
  coordinates = st_coordinates(geoPackage$geom) 
  return(tibble(treeID = geoPackage$treeID, x = coordinates[, "X"], y = coordinates[, "Y"], z = coordinates[, "Z"], height = geoPackage$height))
}


## comparison of tree height distributions
# trees2016 from trees/height-diameter/setup.R
#trees2016 %>% group_by(StandID) %>% summarize(standArea = sum(standArea[1])) %>% summarize(measureArea = sprintf("%.4f", sum(standArea)))
trees2016height = trees2016 %>% mutate(heightClass = if_else(TotalHt > 0, round(TotalHt, 0), if_else(Ht2 > 0, round(Ht2, 0), NA_real_))) %>%
  filter(is.na(heightClass) == FALSE) %>%
  group_by(StandID) %>%
  mutate(liveHeightMeasureTph = sum(isLive * measureTreeTphContribution)) %>%
  group_by(speciesGroup, heightClass) %>%
  summarize(tph = sum(standArea * tph / liveHeightMeasureTph * measureTreeTphContribution) / 15981.0273, .groups = "drop") # 15981.0273 ha = total area of stands measured

dataSourcePath = file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/treetops DSM ring")
treetopTiles = read_sf(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg"))
treetopTileNames = (treetopTiles %>% filter(bufferDistance == 0))$Tile_ID
trees2021 = list()
for (treetopFileName in treetopTileNames) # ~1.8 minutes to load all 663 tiles
{
  treetopTileName = str_remove(treetopFileName, "\\.gpkg")
  tile = read_sf(file.path(dataSourcePath, str_c(treetopFileName, ".gpkg")))
  tileCoordinates = st_coordinates(tile$geom)
  tile = tibble(tileName = treetopTileName, treeID = tile$treeID, x = tileCoordinates[, "X"], y = tileCoordinates[, "Y"], elevation = tileCoordinates[, "Z"], height = tile$height)
  trees2021[[treetopTileName]] = tile
}
trees2021 = bind_rows(trees2021)

segmentedArea = length(treetopTileNames) * 0.3048^2 * 3000^2 / 10000 # ha
trees2021height = left_join(left_join(trees2021 %>% mutate(heightClass = round(0.3048 * height, 0)) %>% # ~4s
                                        group_by(heightClass) %>%
                                        summarize(segmentedTph2021 = n() / segmentedArea),
                                      trees2016height %>% group_by(heightClass) %>% summarize(tph = sum(tph)) %>% rename(tph2016 = tph),
                                      by = c("heightClass")),
                            # crude height growth estimate: concave down parabola
                            trees2016height %>% group_by(heightClass) %>% summarize(tph = sum(tph), .groups = "drop") %>% 
                              mutate(heightClass = round(heightClass + (2021 - 2015) * pmax(0.5 + 0.06 * heightClass - 0.0012 * heightClass^2, 0), 0)) %>% 
                              group_by(heightClass) %>% summarize(tph = sum(tph), .groups = "drop") %>% rename(grownTph2021 = tph),
                            by = c("heightClass"))
#print(trees2021 %>% filter(height > 300), n = 75)
#ggplot() + geom_line(aes(x = seq(0, 80), y = pmax(0.5 + 0.06 * seq(0, 80) - 0.0012 * seq(0, 80)^2, 0)))

ggplot() +
  geom_col(aes(x = heightClass, y = tph, group = speciesGroup, fill = speciesGroup), trees2016height) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 55)) +
  labs(x = "2016 height, m", y = "trees per hectare", fill = NULL) +
  scale_fill_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
ggplot() +
  geom_col(aes(x = heightClass, y = segmentedTph2021), trees2021height) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 55)) +
  labs(x = "2021 height, m", y = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
ggplot() +
  geom_line(aes(x = heightClass, y = 100 * segmentedTph2021 / grownTph2021), trees2021height, na.rm = TRUE) +
  geom_smooth(aes(x = heightClass, y = 100 * segmentedTph2021 / grownTph2021), trees2021height, alpha = 0.1, formula = y ~ x, method = "loess", na.rm = TRUE, span = 0.5) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "2021 height, m", y = "estimated fraction of treetops detected, %") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_x_continuous(breaks = seq(0, 100, by = 10))

tibble(tph2016 = sum(trees2016height$tph), tph2021 = sum(trees2021height$tph)) %>% 
  mutate(total2016M = (33408.277 + 318.838) * tph2016 / 1E6,
         total2021M = (33408.277 + 318.838) * tph2021 / 1E6,
         segmentedPct = 100 * tph2021 / tph2016)


## diff two treetop identifications
referenceTreetops = read_treetops("GIS/DOGAMI/2021 OLC Coos County/treetops DSM/s04020w06690.gpkg")
newTreetops = read_treetops("GIS/DOGAMI/2021 OLC Coos County/treetops DSM/s04020w06690 new.gpkg")

diff = diff_treetops(referenceTreetops, newTreetops)
diff %>% filter((treeID == TRUE) | (x == TRUE) | (y == TRUE) | (z == TRUE) | (height == TRUE))

referenceTreetops %>% arrange(desc(y), x, z)
newTreetops %>% arrange(desc(y), x, z)

