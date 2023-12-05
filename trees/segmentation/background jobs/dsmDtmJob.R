## noise identification and DSM[+DTM] extraction
# Seven concurrent jobs fully utilize 128 GB DDR with 2021 tiles and induce some disk swapping. Six workers also regularly
# reach 128 GB and five workers will also approach 128 GB (20-30 MB per worker is typical). Five or six workers max suggested 
# (depending on other tasks performed while tile jobs are running running) unless more DDR is available.
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

# even chunk sizes for 663 total tiles across five workers = ceiling(663/5) = 133
# chunk  indices   5950X runtime (h:mm) => pull tiles forward from chunk 4 to balance runtime
# 1      1:133     6:40
# 2      134:266   7:57
# 3      267:399   8:46
# 4      400:533   9:41
# 5      534:663   8:30
chunkIndex = 1
chunkSize =  133

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2009 South Coast"
#dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataDestinationPath = dataSourcePath

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

tileNames = list.files(file.path(dataSourcePath, "pointz"), pattern = "\\.laz$")
startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  tile = read_tile(file.path(dataSourcePath, "pointz", tileName))
  tileBaseName = str_remove(tileName, "\\.laz$")

  # check for and then remove other high and low noise points
  #tile = classify_noise(tile, ivf(res = 3, n = 1))
  #writeLAS(tile, file.path(dataSourcePath, "points with noise", paste0(tileBaseName, ".laz")))
  #tile = filter_poi(tile, (Classification != 7L) & (Classification != 18L))
  # generate DSM
  tileDsm = rasterize_canopy(tile, res = 1.5, shape = bbox(tile)) # get tile's digital surface model at 1.5 foot (45.7 cm) resolution using default p2r() algorithm, shape is passed as lidR by default crops the DSM to be smaller than its tile when leading or trailing rows or columns in the .laz contain no points
  writeRaster(tileDsm, file.path(dataDestinationPath, "DSM", paste0(tileBaseName, ".tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
  
  # if needed, a DTM can also be generated
  # TIN meshing is fairly slow, though, and it's debatable whether increasing resolution from DOGAMI's three foot default is
  # useful for reasons other than interpolating DTM tiles to a resolution which matches DSM tiles. Interpolation's done more
  # efficiently with gdalwarp in QGIS (Warp reproject) over a DTM .vrt, then exporting the upsampled raster as a .vrt (tick 
  # save as vrt in QGIS's export raster dialog), and renaming the tiles (GIS/DOGAMI/2021 OLC Coos County\rename tiles.R).
  #tileDtm = rasterize_terrain(tile, res = 1.5)
  #writeRaster(tileDtm, file.path(dataDestinationPath, "DTM", paste0(tileBaseName, ".tif")), gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

#if (exists("generateVrt"))
#{
#  dsmTilePaths = list.files(file.path(dataDestinationPath, "DSM"), pattern = "\\.tif$")
#  vrt(file.path(dataDestinationPath, "DSM", dsmTilePaths), file.path(dataDestinationPath, "DSM", "DSM.vrt"))
#}

cat(paste0(length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
