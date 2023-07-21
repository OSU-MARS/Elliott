## rename downsampled DTM tiles
# DOGAMI DTM resolution is three feet. Building a virtual raster of the DTM and resampling with gdalwarp (reproject in QGIS)
# allows DTM resolution to be matched with higher resolution DSMs obtained from LiDAR point clouds, simplifying local maxima 
# identification and segmentation implementations and reducing CHM estimation error due to 1:1 correspondence between DSM and
# DTM raster cells. Saving the resampled DTM as a .vrt in QGIS (right click -> export) with 3000x3000 foot tiles produces
# DTM tiles spatially aligned to the DSM tiles. However, the tiles are named DTM.0.tif, DTM.1.tif, ... which does not match
# the s<easting>w<northing> naming convention used for point cloud and DTM tiles. The tile layout follows
#
#      s<easting>w<northing> tile names            incrementing .vrt tile names
# ...                                         ... 
# ...  s03690w06390  s03720w06390  ...        ...  908  909 ...
#                    s03720w06360  ...                  939 ...
#
# Tile indexing (gdaltindex) creates a polygon vector layer whose location field contains each tile's filename. Joining by 
# location to the DSM tile index (with a within predicate) adds the Tile_ID, binding DTM.0.tif file names to s<easting>w<northing>.
# The .tif files then need to be renamed to s<easting>w<northing>.tif which (given a joined tile index exported from QGIS) this 
# script accomplishes with two lines of code (though PyQGIS could also be used).
library(dplyr)
library(readxl)
library(stringr)

# if renamed path is not absolute file is moved (or copied) to path relative to working directory
dtmTiles = read_xlsx("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DTM/tile index.xlsx") %>%
  mutate(renamedPath = file.path(dirname(location), paste0(Tile_ID, str_extract(location, "\\.(\\w)+$"))))
file.rename(dtmTiles$errorPath, dtmTiles$renamedPath) # should return only true, indicating successful renames
