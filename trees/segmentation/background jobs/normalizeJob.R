library(lidR)
library(stringr)
library(terra)

read_tile = function(tilePath)
{
  # has classes 1 (unclassified), 2 (ground), and maybe 7 (low noise), 9 (water), 17 (bridge deck), 18 (high noise), 20 (ignored ground, https://www.usgs.gov/ngp-standards-and-specifications/lidar-base-specification-revision-history), or others
  tile = readLAS(tilePath)
  tileClasses = unique(tile$Classification)
  
  tileBaseName = str_remove(basename(tilePath), "\\.laz$")
  cat(paste0(tileBaseName, ": classes ", str_c(sort(tileClasses), collapse = " "), "...\n"))
  
  if ((7 %in% tileClasses) | (18 %in% tileClasses))
  {
    # exclude low and high noise points classified by NV5
    tile = filter_poi(tile, (Classification != 7L) & (Classification != 18L))
  }
  
  return(tile)
}

# 663 tiles -> 7 chunks @ 100 tiles/chunk, ~1.8 minutes/tile -> ~3.0 hours/job, five jobs = (89)~105(125) GB DDR, highest during disk IO pulses
chunkIndex = 1
chunkSize = 100

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataDestinationPath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileNames = list.files(file.path(dataSourcePath, "pointz"), pattern = "\\.laz$")
startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  tileBaseName = str_remove(tileName, "\\.laz$")
  normalizedTilePath = file.path(dataDestinationPath, "pointz normalized", paste0(tileBaseName, ".laz"))
  if (file.exists(normalizedTilePath))
  {
    next
  }
  
  # lasheight.exe -ignore_class 7 18 -replace_z is around 65 seconds total per tile (load + normalize + write)
  # lidR: 61 s = 28 s load + 17 s normalize + 16 s write
  tile = read_tile(file.path(dataSourcePath, "pointz", tileName)) # status update + low and high noise removal
  tileDtm = rast(file.path(dataSourcePath, "DTM", paste0(tileBaseName, ".tif")))
  tile = normalize_height(tile, tileDtm)
  writeLAS(tile, normalizedTilePath)

  #standardMetrics = pixel_metrics(tile, .stdmetrics, res = 60)
  #writeRaster(standardMetrics, file.path(dataDestinationPath, "stdmetrics 60 ft", paste0(tileBaseName, ".tif")), datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

print(paste0(length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), "."))
