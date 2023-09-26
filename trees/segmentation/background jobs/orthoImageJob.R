## orthoimage extraction from raytraced RGBIR point clouds
library(data.table)
library(lidR)
library(stringr)
library(terra)

#get_highest_point_rgbn = function(r, g, b, nir, z)
#{
#  highestPointIndex = which.max(z)
#  bands = list(R = r[highestPointIndex],
#               G = g[highestPointIndex],
#               B = b[highestPointIndex],
#               NIR = nir[highestPointIndex])
#  return(bands)
#}

# extension of https://stackoverflow.com/questions/66841801/creating-orthomosaic-from-las-point-cloud-in-r to include near infrared
metrics_rgbn = function(returnNumber, r, g, b, nir)
{
  firstReturns = which(returnNumber == 1)
  return(data.table(r = sqrt(mean(r[firstReturns]^2)), g = sqrt(mean(g[firstReturns]^2)), b = sqrt(mean(b[firstReturns]^2)), nir = sqrt(mean(nir[firstReturns]^2))))
}

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

# worker-chunk setup much the same as dsmDtmJob.R
# even chunk sizes for 510 total tiles across five workers = ceiling(512/5) = 103
# even chunk sizes for ~15 hour runs: ~2.8 chunks/hour * 15 = 45 => 12 chunks => 3 15 hour runs @ 4 workers/run
chunkIndex = 1 # night 1: chunks 1-5, tiles 1-200, ~64 GB DDR, night 2: chunks 6-13, tiles 201-510, (90)~110(125) GB DDR
chunkSize = 40

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

dataSourcePath = "E:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataDestinationPath = dataSourcePath
tileNames = list.files(file.path(dataSourcePath, "pointz RGBIR"), pattern = "\\.laz$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  tileBaseName = str_remove(tileName, "\\.laz$")
  orthoimageFilePath = file.path(dataDestinationPath, "orthoimage", paste0(tileBaseName, ".tif"))
  if (file.exists(orthoimageFilePath))
  {
    next
  }
  
  tile = read_tile(file.path(dataSourcePath, "pointz RGBIR", tileName)) # status update + low and high noise removal
  orthoimage = pixel_metrics(tile, ~metrics_rgbn(ReturnNumber, R, G, B, NIR), res = 1) # 1 foot orthoimage
  writeRaster(orthoimage, orthoimageFilePath, datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

cat(paste0(length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
