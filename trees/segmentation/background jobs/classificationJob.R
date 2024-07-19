library(caret)
library(dplyr)
library(ranger)
library(stringr)
library(terra)
library(tidyr)
library(tidyterra)

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

# even chunk sizes for 510 total tiles across eight workers = ceiling(510/8) = 64
chunkIndex = 1 # peak of 60 GB DDR per worker
chunkSize = 256

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dsmSourcePath = file.path(dataPath, "DSM v3 beta")
orthoimageSourcePath = file.path(dataPath, "orthoimage v3")
tileNames = list.files(orthoimageSourcePath, "\\.tif$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit vsurf 133.7 10nn m lanczos.Rdata"))
gridMetrics = rast(file.path(dataPath, "metrics", "grid metrics 10 m non-normalized v2.tif"))

dataDestinationPath = file.path(dataPath, "classification vsurf 133.7 10nn m lanczos")


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
  
  #tileDsm = rast(file.path(dsmSourcePath, tileName))
  tileOrthoimage = c(rast(file.path(orthoimageSourcePath, tileName)), rast(file.path(orthoimageSourcePath, "nPoints", tileName)))
  tileGridMetrics = resample(gridMetrics, tileOrthoimage, method = "lanczos") %>% # ~18 s bilinear
    rename(intensityFirstGridReturn = intensityFirstReturn)
  
  tileOrthoimageCrs = crs(tileOrthoimage, describe = TRUE)
  #tileDsmCrs = crs(tileOrthoimage, describe = TRUE)
  tileGridMetricsCrs = crs(tileGridMetrics, describe = TRUE)
  #if ((tileOrthoimageCrs$authority == tileDsmCrs$authority) & (tileOrthoimageCrs$code == tileDsmCrs$code))
  #{
  #  crs(tileDsmCrs) = crs(tileOrthoimage) # give DSM the orthoimage's vertical CRS
  #} else {
  #  stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match DSM CRS ", tileGridMetricsCrs$authority, ":", tileDsmCrs$code, "."))
  #}
  if ((tileOrthoimageCrs$authority == tileGridMetricsCrs$authority) & (tileOrthoimageCrs$code == tileGridMetricsCrs$code))
  {
    crs(tileGridMetrics) = crs(tileOrthoimage) # give grid metrics the orthoimage's vertical CRS
  } else {
    stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match grid metrics CRS ", tileGridMetricsCrs$authority, ":", tileGridMetricsCrs$code, "."))
  }
  
  #tileStatistics = c(tileOrthoimage, tileDsm, tileGridMetrics)
  tileStatistics = c(tileOrthoimage, tileGridMetrics)
  tileStatistics$cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))
  tileStatisticsTibble = as_tibble(tileStatistics) %>% # ~11 s
    filter(firstReturns > 0, is.na(intensityFirstGridReturn) == FALSE) %>% # exclude pixels without data and grid metrics cells without points
    #rename(nir = nearInfrared) %>%
    mutate(intensitySecondReturn = replace_na(intensitySecondReturn, 0),
           gli = (2*green - blue - red) / (2*green + blue + red),
           zMaxNormalized = zMax - zGroundMean,
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
           zQ75normalized = zQ75 - zGroundMean,
           #zQ80normalized = zQ80 - zGroundMean,
           #zQ85normalized = zQ85 - zGroundMean,
           #zQ90normalized = zQ90 - zGroundMean,
           zQ95normalized = zQ95 - zGroundMean)
  #colSums(is.na(tileStatisticsTibble))
  
  tileClassification = predict(randomForestFit$finalModel, data = tileStatisticsTibble) # ~2.1 m/tile, predict(randomForestFit, newdata = ...) fails on internal matrix datatype error

  tileClassificationTibble = left_join(tibble(cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))),
                                       tibble(cellID = tileStatisticsTibble$cellID, classification = tileClassification$predictions),
                                       by = "cellID")
  
  classificationRaster = rast(tileStatistics, nlyrs = 1, names = "classification", vals = tileClassificationTibble$classification)
  writeRaster(classificationRaster, file.path(dataDestinationPath, tileName), datatype = "INT1U", overwrite = TRUE) #, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9")
}

cat(paste0("Orthoimage vegetation classification over ", length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
