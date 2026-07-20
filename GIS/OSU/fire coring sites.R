library(dplyr)
library(sf)

fireSites = st_read("GIS/OSU/fire coring sites.gpkg", layer = "fire coring sites", quiet = TRUE)
fireTiles = unique(fireSites$Tile_ID)

sourceRoot = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
destinationRoot = "//for-mars-viridi/Elliott/GIS/DOGAMI/2021 OLC Coos County"
sourcePaths = c("tiles RGB+NIR", "treetops/rf v2", "DSM v3/local maxima", "DSM v3/local maxima", "DSM v3/sourceID")
fileExtensions = c(".las", ".gpkg", ".gpkg", ".tif", ".tif")
for (sourceIndex in 1:length(sourcePaths))
{
  sourcePath = sourcePaths[sourceIndex]
  fileExtension = fileExtensions[sourceIndex]
  for (tileName in fireTiles)
  {
    fileName = paste0(tileName, fileExtension)
    cat(paste0(fileName, "...\n"))
    file.copy(file.path(sourceRoot, sourcePath, fileName), file.path(destinationRoot, sourcePath, fileName), overwrite = FALSE, copy.date = TRUE)
  }
}