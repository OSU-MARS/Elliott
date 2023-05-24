library(dplyr)
library(readxl)

# remove unneeded tiles
# Downloading { 43123+43124, D+E+F } from ftp://lidar.engr.oregonstate.edu/OREGON%20LIDAR%20CONSORTIUM%20PROJECT%20DATA/OLC%20SOUTH%20COAST%202009/POINTS/
# pulls 1128 tiles, 322 of which are on the Elliott and 460 of which are within 2 km of the Elliott
southCoastTiles = read_xlsx("GIS/DOGAMI/2009 OLC South Coast/OLC South Coast 2009 ESRF and neighboring tiles.xlsx")

#tilesDeleted = 0
for (index in 1:nrow(southCoastTiles))
{
  if (southCoastTiles$bufferDist[index] <= 2000)
  {
    next
  }
  
  tilePath = file.path(getwd(), "GIS/DOGAMI/2009 OLC South Coast/Points", paste0(southCoastTiles$name[index], ".laz"))
  if (file.exists(tilePath))
  {
    file.remove(tilePath)
    #print(paste("Delete", tilePath))
    #tilesDeleted = tilesDeleted + 1
  }
}

## lidR
#library(lidR)
#southCoastTiles = readLAScatalog("GIS/DOGAMI/2009 OLC South Coast/Points")
#southCoastBoundaries = catalog_boundaries(southCoastTiles, concavity = Inf, length_threshold = 100) # 2+ orders of magnitude slower than using lasboundary.exe