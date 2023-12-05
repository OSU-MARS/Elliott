library(readxl)

sourcePath = "T:/Groups/ElliottStateResearchForest"
#destinationPath = file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County")
#destinationPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
destinationPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileIndex = read_xlsx("GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.xlsx")

# copy .las tiles from network share to local drive with support for semi-manual concurrent compression
# Assumed directory structure:
#   downloading - copy destination point
#   Points - location fully downloaded .las files are staged to to avoid .las to .laz conversion running on partially downloaded files
#   Pointz - end location of local .laz compressed versions of network share's .las files
# Typical download size is circa 1.5 TB, depending on tile index selected in GIS, which is about four hours to pull over an
# uncongested 1 Gb Ethernet link (~110 MB/s) and may be closer to 15 hours during academic year daytime activity. laszip64.exe 
# compresses about three tiles per minute so compression takes about as many core hours as an uncongested download requires wall clock
# time. Running laszip64.exe in parallel from multiple command windows on Points, Points2, ... can be helpful catch up if downloading
# has gotten ahead of compression, albeit primarily on NVMe and SSD destinations.
# 
# .las to .laz compression ratios are about 4.9x.
for (tileName in tileIndex$Tile_ID)
{
  #existingLazPath = file.path(destinationPath, "pointz", paste0(tileName, ".laz"))
  #if (file.exists(existingLazPath))
  #{
  #  cat(paste0(tileName, " already copied.\n"))
  #  next
  #}
  existingLasPath = file.path(destinationPath, "points", paste0(tileName, ".las"))
  if (file.exists(existingLasPath))
  {
    cat(paste0(tileName, " already copied.\n"))
    next
  }
  
  processedLasPath = file.path(sourcePath, "ProcessedData/Lidar/OLC Coos County 2021/NV5_Class12", paste0(tileName, "_las.las"))
  destinationLasPath = file.path(destinationPath, "downloading", paste0(tileName, ".las"))
  if (file.exists(processedLasPath))
  {
    cat(paste0(tileName, " from NV5_Class12...\n"))
    file.copy(processedLasPath, destinationLasPath)
  } else {
    rawLasPath = file.path(sourcePath, "RawData/ESRF_Initiation2021/NV5_Geospatial/LIDAR/Points/LAS", paste0(tileName, ".las"))
    if (file.exists(rawLasPath))
    {
      cat(paste0(tileName, " from NV5_Geospatial...\n"))
      file.copy(rawLasPath, destinationLasPath)
    } else {
      cat(paste0(tileName, " not found!\n"))
    }
  }
}

# copy DTM tiles from network share
for (tileName in tileIndex$Tile_ID)
{
  destinationDtmPath = file.path(destinationPath, "Bare_Earth", paste0("be_", tileName, ".tif"))
  if (file.exists(destinationDtmPath))
  {
    #cat(paste0(destinationDtmPath, " already copied.\n"))
    next
  }
  
  sourceDtmPath = file.path(sourcePath, "RawData/OLC Coos County 2021/NV5_Elliott/DTM", paste0("be_", tileName, ".tif"))
  if (file.exists(sourceDtmPath))
  {
    cat(paste0(tileName, " from NV5_Elliott...\n"))
    file.copy(sourceDtmPath, destinationDtmPath)
  } else {
    sourceDtmPath = file.path(sourcePath, "RawData/OLC Coos County 2021/NV5_Geospatial/LIDAR/Raster/Bare_Earth", paste0("be_", tileName, ".tif"))
    if (file.exists(sourceDtmPath))
    {
      cat(paste0(tileName, " from NV5_Geospatial...\n"))
      file.copy(sourceDtmPath, destinationDtmPath)
    } else {
      cat(paste0(tileName, " not found!\n")) 
    }
  }
}

## 2009 tiles
# archive downloaded tiles beyond 400 m buffer distance
#tiles2009 = list.files("GIS/DOGAMI/2009 OLC South Coast/pointz", "\\.laz$")
#tileIndex = read_xlsx("GIS/DOGAMI/2009 OLC South Coast/OLC South Coast 2009 ESRF and neighboring tiles.xlsx")
#tilesToMove = tileIndex %>% filter(bufferDist > 400, paste0(name, ".laz") %in% tiles2009) %>%
#  mutate(path = file.path(getwd(), "GIS/DOGAMI/2009 OLC South Coast/Pointz", paste0(name, ".laz")),
#         destination = file.path("E:/Elliott/GIS/DOGAMI/2009 South Coast/pointsz/500+", paste0(name, ".laz")))
#file.rename(tilesToMove$path, tilesToMove$destination)
