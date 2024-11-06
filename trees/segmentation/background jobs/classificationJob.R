library(caret)
library(dplyr)
library(FNN)
library(magrittr)
library(ranger)
library(stringr)
library(terra)
library(tidyr)
library(tidyterra)

jobStartTime = Sys.time()

# even chunk sizes for 510 total tiles across eight workers = ceiling(510/8) = 64
chunkIndex = 1 # peak of 60 GB DDR per worker
chunkSize = 256

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dsmSourcePath = file.path(dataPath, "DSM v3 beta")
dtmSourcePath = file.path(dataPath, "DTM")
orthoimageSourcePath = file.path(dataPath, "orthoimage v3")
tileNames = list.files(orthoimageSourcePath, "\\.tif$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit vsurf 172.4 10 m cubic.Rdata"))
gridMetrics = rast(file.path(dataPath, "metrics", "grid metrics 10 m non-normalized v2.tif"))

dataDestinationPath = file.path(dataPath, "classification vsurf 172.4 10m cubic")

imageCenters = vect("GIS/DOGAMI/2021 OLC Coos County/image positions.gpkg") # terra drops z coordinates but they're redundant with the elevation field
imageCentersForKnn = tibble(xy = geom(imageCenters), z = imageCenters$`elevation, ft`) %>% # flatten to coordinates for kNN
  mutate(x = xy[, "x"], y = xy[, "y"]) %>% select(-xy) %>% relocate(x, y, z)
imageSunPositions = tibble(azimuth = imageCenters$sunAzimuth, elevation = imageCenters$sunElevation) # flatten for in loop lookup


## classify orthoimages plus DSM
# 5950X single thread: TBD 2000 x 2000 tiles/min
cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  classificationFilePath = file.path(dataDestinationPath, tileName)
  if (file.exists(classificationFilePath))
  {
    next
  }
  cat(paste0(str_remove(tileName, "\\.tif"), "...\n"))
  
  tileDsm = rast(file.path(dsmSourcePath, tileName))
  tileDtm = rast(file.path(dtmSourcePath, tileName))
  tileOrthoimage = c(rast(file.path(orthoimageSourcePath, tileName)), rast(file.path(orthoimageSourcePath, "nPoints", tileName)))
  tileGridMetrics = resample(gridMetrics, tileOrthoimage, method = "lanczos") %>% # ~18 s bilinear
    rename(intensityFirstGridReturn = intensityFirstReturn)
  
  tileOrthoimageCrs = crs(tileOrthoimage, describe = TRUE)
  tileDsmCrs = crs(tileDsm, describe = TRUE)
  tileDtmCrs = crs(tileDtm, describe = TRUE)
  tileGridMetricsCrs = crs(tileGridMetrics, describe = TRUE)
  if ((tileOrthoimageCrs$authority == tileDsmCrs$authority) & (tileOrthoimageCrs$code == tileDsmCrs$code))
  {
    crs(tileDsm) = crs(tileOrthoimage) # give DSM the orthoimage's vertical CRS
  } else {
    stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match DSM CRS ", tileGridMetricsCrs$authority, ":", tileDsmCrs$code, "."))
  }
  if ((tileOrthoimageCrs$authority == tileDtmCrs$authority) & (tileOrthoimageCrs$code == tileDtmCrs$code))
  {
    crs(tileDtm) = crs(tileOrthoimage) # give DTM the orthoimage's vertical CRS
  } else {
    stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match DTM CRS ", tileGridMetricsCrs$authority, ":", tileDsmCrs$code, "."))
  }
  if ((tileOrthoimageCrs$authority == tileGridMetricsCrs$authority) & (tileOrthoimageCrs$code == tileGridMetricsCrs$code))
  {
    crs(tileGridMetrics) = crs(tileOrthoimage) # give grid metrics the orthoimage's vertical CRS
  } else {
    stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match grid metrics CRS ", tileGridMetricsCrs$authority, ":", tileGridMetricsCrs$code, "."))
  }

  tileViewAngles = tibble(crds(tileDsm, df = TRUE, na.rm = FALSE), dsm = as.vector(tileDsm$dsm)) %>%
    mutate(dsm = if_else(is.na(dsm), as.vector(tileDtm), dsm))
  imageCenterIndex = knnx.index(imageCentersForKnn, tileViewAngles, k = 1)
  tileViewAngles %<>% mutate(sunAzimuth = imageSunPositions$azimuth[imageCenterIndex], 
                             sunElevation = imageSunPositions$elevation[imageCenterIndex], 
                             deltaX = imageCentersForKnn$x[imageCenterIndex] - x,
                             deltaY = imageCentersForKnn$y[imageCenterIndex] - y,
                             deltaZ = imageCentersForKnn$z[imageCenterIndex] - dsm,
                             viewAzimuth = -180/pi * (atan2(deltaY, deltaX) - pi/2),
                             viewAzimuth = if_else(viewAzimuth > 0, viewAzimuth, 360 + viewAzimuth),
                             viewZenithAngle = 180/pi * atan2(sqrt(deltaX^2 + deltaY^2), deltaZ))
  tileDsm$sunAzimuth = tileViewAngles$sunAzimuth
  tileDms$sunElevation = tileViewAngles$sunElevation
  tileDsm$viewAzimuth = tileViewAngles$viewAzimuth
  tileDsm$viewZenithAngle = tileViewAngles$viewZenithAngle
  
  tileStatistics = c(tileOrthoimage, tileDsm, tileGridMetrics)
  tileStatistics$cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))
  tileStatisticsTibble = as_tibble(tileStatistics) %>% # ~11 s
    filter(firstReturns > 0, is.na(intensityFirstGridReturn) == FALSE) %>% # exclude pixels without data and grid metrics cells without points
    #rename(nir = nearInfrared) %>%
    mutate(chm = replace_na(chm, 0), # if no LiDAR hits, assume DSM = DTM => CHM = 0
           intensitySecondReturn = replace_na(intensitySecondReturn, 0), # if no LiDAR hits, consider second return intensity to be zero
           luminosity = 0.299 * red + 0.587 * green + 0.114 * blue, # NTSC, ITU BT.610
           gNormalized = green / luminosity,
           viewElevation = 90 - viewZenithAngle,
           viewAzimuthSunRelative = sunAzimuth - viewAzimuth, # 0 = forward scatter, ±180 = backscatter, absolute value broken out separately to allow for asymmetric BRDF
           viewAzimuthSunRelative = if_else(viewAzimuthSunRelative > 180, 360 - viewAzimuthSunRelative, if_else(viewAzimuthSunRelative < -180, 360 + viewAzimuthSunRelative, viewAzimuthSunRelative)), # clamp to [0, ±180] to constrain training complexity
           viewAzimuthSunRelativeAbsolute = abs(viewAzimuthSunRelative))
           #gli = (2*green - blue - red) / (2*green + blue + red),
           #zMaxNormalized = zMax - zGroundMean,
           #zMeanNormalized = zMean - zGroundMean,
           #zQ05normalized = zQ05 - zGroundMean,
           #zQ10normalized = zQ10 - zGroundMean,
           #zQ15normalized = zQ15 - zGroundMean,
           #zQ20normalized = zQ20 - zGroundMean,
           #zQ25normalized = zQ25 - zGroundMean,
           #zQ30normalized = zQ30 - zGroundMean,
           #zQ35normalized = zQ35 - zGroundMean,
           #zQ40normalized = zQ40 - zGroundMean,
           #zQ45normalized = zQ45 - zGroundMean,
           #zQ50normalized = zQ50 - zGroundMean,
           #zQ55normalized = zQ55 - zGroundMean,
           #zQ60normalized = zQ60 - zGroundMean,
           #zQ65normalized = zQ65 - zGroundMean,
           #zQ70normalized = zQ70 - zGroundMean,
           #zQ75normalized = zQ75 - zGroundMean,
           #zQ80normalized = zQ80 - zGroundMean,
           #zQ85normalized = zQ85 - zGroundMean,
           #zQ90normalized = zQ90 - zGroundMean,
           #zQ95normalized = zQ95 - zGroundMean)
  #colSums(is.na(tileStatisticsTibble))
  
  tileClassification = predict(randomForestFit$finalModel, data = tileStatisticsTibble) # ~2.1 m/tile, predict(randomForestFit, newdata = ...) fails on internal matrix datatype error

  tileClassificationTibble = left_join(tibble(cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))),
                                       tibble(cellID = tileStatisticsTibble$cellID, classification = tileClassification$predictions),
                                       by = "cellID")
  
  classificationRaster = rast(tileStatistics, nlyrs = 1, names = "classification", vals = tileClassificationTibble$classification)
  writeRaster(classificationRaster, file.path(dataDestinationPath, tileName), datatype = "INT1U", overwrite = TRUE) #, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9")
}

cat(paste0("Orthoimage vegetation classification over ", length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
