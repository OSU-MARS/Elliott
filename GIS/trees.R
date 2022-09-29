library(arrow)
library(dplyr)
library(readr)
library(stringr)
library(terra)
library(tidyr)

## transcode .gpkg to .feather
startTime = proc.time()
elliottTrees = as_tibble(vect("GIS/Trees/ESRF_TreesD_H10Cr20h10A50MD7s with IDs.gpkg")) # ~2.5 minutes @ 5.7 GB, 12.6 million trees
(elapsedTime = proc.time() - startTime)

elliottTreesMod = elliottTrees %>% # 2-3 seconds
  filter(is.na(resourceUniTid) == FALSE, resourceUniTbufferDist <= 1000) %>% # TODO: convert trees with DBH < 5 cm to saplings?
  mutate(StandID = replace_na(StandID, 0), # reassign NULL stand IDs (signalling trees outside of the Elliott) to iLand's default stand
         species = factor(202, labels = c("psme"), levels = c(202)),
         treeID = 100000 * gridID + treeID, 
         x = 0.3048 * x , # convert from English units to metric, assuming input is in EPSG:6557 (nudge trees off resource unit boundaries: - if_else(id %in% c(281406309, 342014985, 412900427), 0.06, 0))
         y = 0.3048 * y, # for now, use x and y due to error in XTOP and YTOP ( - if_else(id %in% c(512801827, 530705841), 0.06, 0))
         Z = 0.3048 * Z, # use z as ZTOP contains high side outliers
         dbh = 0.1822 * exp(2.7042*(Z - 1.37)^0.2115), # Wyckoff height regression from Elliott Stand Data Feb2022.R: height in m -> DBH in cm
         resourceUnitX = as.integer(x / 100),
         resourceUnitY = as.integer(y / 100)) %>% 
  rename(standID = StandID, height = Z) %>%
  arrange(resourceUnitY, resourceUnitX, species, y, x) %>%
  select(standID, species, treeID, x, y, dbh, height)

elliottTreesMod %>% summarize(minDbh = min(dbh), maxDbh = max(dbh), minHeight = min(height), maxHeight = max(height), naDbh = sum(is.na(dbh)), naHeight = sum(is.na(height)), naX = sum(is.na(x)), naY = sum(is.na(y)), stands = length(unique(standID)))
#elliottTreesMod %>% filter((x %% 100.0 == 0.0) | (y %% 100.0 == 0.0)) # try to check for trees exactly on resource unit edges: ineffective due to rounding


elliottTreesArrow = arrow_table(elliottTreesMod %>% mutate(fiaCode = recode(species, "psme" = 202)) %>% select(-species) %>% relocate(standID, treeID, fiaCode, dbh, height, x, y),
                                schema = schema(standID = int32(), treeID = int32(), fiaCode = uint16(),
                                                dbh = float32(), height = float32(), x = float32(), y = float32()))
write_feather(elliottTreesArrow, "iLand/init/TreesD_H10Cr20h10A50MD7s 1 km buffer.feather", compression = "uncompressed")

elliottTrees %>% filter(resourceUniTid == 73606, treeID == 128) %>% select(StandID, treeID, XTOP, YTOP, ZTOP,resourceUniTid, resourceUniTbufferDist)
(elliottTreesMod %>% filter(id == 342014985))$x
(elliottTreesMod %>% filter(id == 412900427))$

elliottTreesMod %>% filter(id == 342000128) %>% select(standID, species, id, x, y, dbh, height, resourceUniTid, resourceUniTbufferDist)

## transcode .csv files exported from QGIS to .feather
tileCsvFileNames = list.files(path = "GIS/Trees", pattern = "^TSegD.*\\.csv")
tileColumnTypes = cols(id = "i", species = "c", standID = "i", .default = "d")
for (tileCsvFileName in tileCsvFileNames)
{
  featherFilePath = paste0("iLand/init/", str_replace(tileCsvFileName, ".csv", ".feather"))
  if (file.exists(featherFilePath) == FALSE)
  {
    tile = read_csv(paste0("GIS/Trees/", tileCsvFileName), col_types = tileColumnTypes) %>%
      mutate(resourceUnitX = as.integer(x / 100),
             resourceUnitY = as.integer(y / 100)) %>%
      arrange(resourceUnitY, resourceUnitX, species, y, x) %>%
      select(-resourceUnitY, -resourceUnitX) %>%
      rename(treeID = id) %>%
      relocate(standID, treeID, species, height, dbh, x, y)
    tileArrow = arrow_table(tile %>% mutate(fiaCode = recode(species, "psme" = 202)) %>% select(-species) %>% relocate(standID, treeID, fiaCode, dbh, height, x, y),
                            schema = schema(standID = int32(), treeID = int32(), fiaCode = uint16(), 
                                            dbh = float32(), height = float32(), x = float32(), y = float32()))
    write_feather(tileArrow, featherFilePath, compression = "uncompressed")
  }
}
