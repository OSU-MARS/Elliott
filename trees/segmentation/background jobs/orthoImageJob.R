## orthoimage extraction from raytraced RGBIR point clouds
library(lidR)
library(stringr)
library(terra)

# extension of https://stackoverflow.com/questions/66841801/creating-orthomosaic-from-las-point-cloud-in-r to include near infrared
get_highest_point_rgbn = function(r, g, b, nir, z)
{
  highestPointIndex = which.max(z)
  bands = list(R = r[highestPointIndex],
               G = g[highestPointIndex],
               B = b[highestPointIndex],
               NIR = nir[highestPointIndex])
  return(bands)
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

# worker-chunk setup pretty much the same as dsmDtmJob.R
# even chunk sizes for 512 total tiles across five workers = ceiling(512/5) = 103
# chunk  indices   5950X runtime (h:mm) => pull tiles forward from chunk 4 to balance runtime
# 1      1:103    
# 2      104:206   
# 3      207:309   
# 4      310:412   
# 5      413:512   
chunkIndex = 1
chunkSize =  103

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
tileNames = list.files(file.path(dataPath, "RGBIR raytrace"), pattern = "\\.laz$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  tile = read_tile(file.path(dataPath, "RGBIR raytrace", tileName))

  orthoStartTime = Sys.time()  
  orthoimage = pixel_metrics(tile, ~get_highest_point_rgbn(R, G, B, NIR, Z), res = 1)
  cat(format(Sys.time() - orthoStartTime))

  tileBaseName = str_remove(tileName, "\\.laz$")
  writeRaster(orthoimage, file.path(dataPath, "orthoimage", paste0(tileBaseName, ".tif")), datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

cat(paste0(length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
