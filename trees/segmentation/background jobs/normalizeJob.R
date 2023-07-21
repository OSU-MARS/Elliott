library(lidR)
library(stringr)
library(terra)

chunkIndex = 1
chunkSize =  116 # see dsmDtmJob.R

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

dataPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileNames = list.files(file.path(dataPath, "points with noise"), pattern = "\\.laz$")
tileNames = tileNames[seq(chunkSize * (chunkIndex - 1) + 1, min(chunkSize * chunkIndex, length(tileNames)))]

for (tileName in tileNames)
{
  tileBaseName = str_remove(tileName, "_las\\.las$")
  cat(paste0(tileBaseName, "..."))
  
  tile = readLAS(file.path(dataPath, "points with noise", tileName))
  tileDtm = readRaster(file.path(dataPath, "DTM", paste0(tileBaseName, ".tif")))
  
  tile = normalize_height(tile, dtm = tileDtm)
  tile = filter_poi(tile, Classification != 18L, Z >= -10, Z <= 350)
  writeLAS(tile, file.path(dataPath, "points normalized", paste0(tileBaseName, ".laz")))

  tileDsm = readRaster(file.path(dataPath, "DSM", paste0(tileBaseName, ".tif")))
  tileChm = tileDsm - tileDtm
  chmHeightQuantiles = quantile(values(tileChm), probs = c(0.85, 0.95), na.rm = TRUE, names = FALSE) # point density dependent: 0.65, 0.95? for 2009 flight?
  names(chmHeightQuantiles) = c("q85", "q95")
  
  get_window_size = function(z)
  {
    return(z * log((chmHeightQuantiles[2] - chmHeightQuantiles[1])) / 25 + (chmHeightQuantiles[2] - chmHeightQuantiles[1]) / 10)
  }
  tileTreetops = locate_trees(chm, lmf(get_window_size, hmin = 3.2808 * 1.5), uniqueness = "incremental") # shortest trees in MBG dataset are 1.52 m (5.0 ft), incremental sequentially numbers the tree tops in the tile
  writeVector(tileTreetops, file.path(dataPath, "treetops normalized", paste0(tileBaseName, ".gpkg")))
}

print(paste0(length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), "."))
