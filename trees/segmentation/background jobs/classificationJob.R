library(caret)
library(dplyr)
library(ranger)
library(stringr)
library(terra)
library(tidyr)

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

# even chunk sizes for 510 total tiles across eight workers = ceiling(510/8) = 64
chunkIndex = 1 # peak of 60 GB DDR per worker
chunkSize = 256

dataPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataSourcePath = file.path(dataPath, "orthoimage v2")
tileNames = list.files(dataSourcePath, "\\.tif$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit ortho6 10 m non-normalized 9 lanczos.Rdata"))
gridMetrics = rast(file.path(dataPath, "metrics", "grid metrics 10 m non-normalized v2.tif"))

dataDestinationPath = file.path(dataPath, "classification ortho6 + 10 m non-normalized 9 lanczos")


## classify orthoimages
# 5950X single thread: TBD 3000 x 3000 tiles/min
cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  classificationFilePath = file.path(dataDestinationPath, tileName)
  if (file.exists(classificationFilePath))
  {
    next
  }
  cat(paste0(str_remove(tileName, "\\.tif"), "...\n"))
  
  tileOrthoimage = rast(file.path(dataSourcePath, tileName))
  tileGridMetrics = resample(gridMetrics, tileOrthoimage, method = "lanczos") # ~18 s bilinear
  
  tileOrthoimageCrs = crs(tileOrthoimage, describe = TRUE)
  tileGridMetricsCrs = crs(tileGridMetrics, describe = TRUE)
  if ((tileOrthoimageCrs$authority == tileGridMetricsCrs$authority) & (tileOrthoimageCrs$code == tileGridMetricsCrs$code))
  {
    crs(tileGridMetrics) = crs(tileOrthoimage) # give grid metrics the orthoimage's vertical CRS
  } else {
    stop(paste0("Orthoimage CRS ", tileOrthoimageCrs$authority, ":", tileOrthoimageCrs$code, " does not match grid metrics CRS ", tileGridMetricsCrs$authority, ":", tileGridMetricsCrs$code, "."))
  }
  
  tileStatistics = c(tileOrthoimage, tileGridMetrics)
  tileStatistics$cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))
  tileStatisticsTibble = as_tibble(tileStatistics) %>% # ~15 s
    filter(firstReturns > 0, is.na(intensityFirstReturn) == FALSE) %>% # exclude pixels without data and grid cells without points
    mutate(red = replace_na(red, 65535), # restore any RGB+NIR and intensity values which were max brightness (firstReturns > 0 indicates there is data but it gets lost due to collision with no data)
           green = replace_na(green, 65535),
           blue = replace_na(blue, 65535),
           nir = replace_na(nir, 65535),
           intensity = replace_na(intensity, 65535),
           intensitySecondReturn = if_else(secondReturns == 0, 0, replace_na(intensitySecondReturn, 65535)), # replace any second return no datas with 0 and fix up no data collisions
           #greenness = green / red,
           #redBlueRatio = red / blue,
           #nirBlueRatio = nir / blue,
           #nirGreenRatio = nir / green,
           #nirRedRatio = nir / red,
           #normalizedGreen = green / (nir + red + green),
           #normalizedNir = nir / (nir + red + green),
           #normalizedBlue	= blue / (nir + red + green),
           #coloration	= (red - blue) / red,
           #chlorophyllGreen = nir / green - 1,
           #chlorophyllVegetation = nir * as.numeric(red) / green^2,
           ndvi = (nir - red) / (nir + red),
           gndvi = (nir - green) / (nir + green),
           bndvi = (nir - blue) / (nir + blue),
           #sipi = (nir - blue) / (nir + red),
           #ndgb = (green - blue) / (green + blue),
           ndgr = (green - red) / (green + red),
           #wbi = (blue - red) / (blue + red),
           #mgrv = (green^2 - r^2) / (green^2 + r^2),
           #rgbv = (green^2 - red * as.numeric(blue)) / (green^2 + red * as.numeric(blue)), # convert explicitly to avoid integer overflow
           #gli = (2*green - blue - red) / (2*green + blue + red),
           #mexg = 1.62 * green - 0.884 * red - 0.311 * blue,
           #luminosity = 0.299 * red + 0.587 * green + 0.114 * blue, # NTSC, ITU BT.610
           #rNormalized = red / luminosity,
           #gNormalized = green / luminosity,
           #bNormalized = blue / luminosity,
           #nirNormalized = nir / luminosity,
           #luminosity709 = 0.2126 * red + 0.7152 * green + 0.0722 * blue,
           #gemi = (2 * (nir^2 - red^2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5),
           #tvi = sqrt(ndvi + 0.5),
           #evi = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1),
           #triangularVegIndex = 0.5 * (120 * (nir - green) - 200 * (red - green)),
           #mtvi2 = 1.5 * (2.5 * (nir - green) - 2.5 * (red - green)) / sqrt((2 * nir + 1)^2 - 6 * nir - 5 * sqrt(red) - 0.5),
           #eviBackup = 2.5 * (nir - red) / (nir + 2.4 * red + 1),
           #savi = 1.5 * (nir - red) / (nir + red + 0.5),
           #msavi = (2 * nir + 1 - sqrt(2 * (2 * nir + 1)^2 - 8 * (nir - red))) / 2,
           #osavi = 1.16 * (nir - red) / (nir + red + 0.16),
           #zNormalizedEntropy = replace_na(zNormalizedEntropy, 0), # keep synchronized with species.R
           #zMaxNormalized = zMax - zGroundMean,
           #zMeanNormalized = zMean - zGroundMean,
           zQ05normalized = zQ05 - zGroundMean,
           zQ10normalized = zQ10 - zGroundMean,
           #zQ15normalized = zQ15 - zGroundMean,
           #zQ20normalized = zQ20 - zGroundMean,
           zQ25normalized = zQ25 - zGroundMean,
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
           zQ95normalized = zQ95 - zGroundMean) %>%
    select(cellID, ndgr, ndvi, gndvi, bndvi, intensity, intensityFirstReturn, intensitySecondReturn, intensityQ10, intensityPground, pThirdReturn, pZaboveThreshold, zQ05normalized, zQ10normalized, zQ25normalized, zQ75normalized)
    # select(cellID, ndgr, ndvi, gndvi, bndvi, intensity, intensitySecondReturn) # imagery only
  colSums(is.na(tileStatisticsTibble))
  
  tileClassification = predict(randomForestFit$finalModel, data = tileStatisticsTibble) # ~2.1 m, predict(randomForestFit, newdata = ...) fails on internal matrix datatype error

  tileClassificationTibble = left_join(tibble(cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))),
                                       tibble(cellID = tileStatisticsTibble$cellID, classification = tileClassification$predictions),
                                       by = "cellID")
  
  classificationRaster = rast(tileStatistics, nlyrs = 1, names = "classification", vals = tileClassificationTibble$classification)
  writeRaster(classificationRaster, file.path(dataDestinationPath, tileName), datatype = "INT1U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

cat(paste0("Orthoimage vegetation classification over ", length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
