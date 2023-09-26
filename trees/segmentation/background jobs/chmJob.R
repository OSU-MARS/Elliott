library(terra)

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataDestinationPath = dataSourcePath

dsmTilePath = file.path(dataSourcePath, "DSM")
dtmTilePath = file.path(dataSourcePath, "DTM")
chmTilePath = file.path(dataSourcePath, "CHM")

## simple extraction of CHM = DSM - DTM
# 5950X single thread: ~120 2000 x 2000 tiles/min
chmStartTime = Sys.time()

dsmTileNames = list.files(dsmTilePath, "\\.tif")
cat(paste0("Processing ", length(dsmTileNames), " tiles..."))
for (tileName in dsmTileNames)
{
  cat(paste0(tileName, "...\n"))
  dsmTile = rast(file.path(dsmTilePath, tileName))
  dtmTile = rast(file.path(dtmTilePath, tileName))
  chmTile = dsmTile - dtmTile
  writeRaster(chmTile, file.path(chmTilePath, tileName), overwrite = TRUE)
}

cat(paste0("Canopy height models for ", length(dsmTileNames), " tiles in ", format(Sys.time() - chmStartTime), ".\n"))
