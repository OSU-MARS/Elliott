library(caret)
library(dplyr)
library(ranger)
library(stringr)
library(terra)
library(tidyr)

jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))

# even chunk sizes for 510 total tiles across eight workers = ceiling(510/8) = 64
chunkIndex = 2 # peak of 60 GB DDR per worker
chunkSize = 256

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
dataDestinationPath = dataSourcePath

tileNames = list.files(file.path(dataSourcePath, "orthoimage"), "\\.tif$")

startIndex = chunkSize * (chunkIndex - 1) + 1
endIndex = min(chunkSize * chunkIndex, length(tileNames))
tileNames = tileNames[startIndex:endIndex]

load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit 10 m non-normalized grid metrics 44 var.Rdata"))
gridMetrics = rast(file.path(dataSourcePath, "metrics", "grid metrics 10 m non-normalized.tif"))
names(gridMetrics)[45] = "intensityPground" # temporary workaround for C# typo

## classify orthoimages
# 5950X single thread: TBD 3000 x 3000 tiles/min
cat(paste0("Processing chunk ", chunkIndex, " (indices ", startIndex, ":", endIndex, ") with ", length(tileNames), " tiles..."))
for (tileName in tileNames)
{
  classificationFilePath = file.path(dataDestinationPath, "classification 10 m non-normalized 44 var", tileName)
  if (file.exists(classificationFilePath))
  {
    next
  }
  cat(paste0(str_remove(tileName, "\\.tif"), "...\n"))
  
  # extract training data for this polygon
  # For now, assume orthoimagery and CHM are in the same CRS and that orthoimagery is higher resolution.
  # For now, assume all CHM cells have Z values as this is the case for the current set of training polygons (October 2023, n = 365)
  tileOrthoimage = rast(file.path(dataSourcePath, "orthoimage", tileName))
  tileGridMetrics = resample(gridMetrics, tileOrthoimage) # ~18 s
  tileStatistics = c(tileOrthoimage, tileGridMetrics)
  tileStatistics$cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))
  tileStatisticsTibble = as_tibble(tileStatistics) %>% # ~15 s
    filter(is.na(n) == FALSE, is.na(r) == FALSE, is.na(g) == FALSE, is.na(b) == FALSE, is.na(nir) == FALSE, 
           is.na(zGroundMean) == FALSE, is.na(zSkew) == FALSE, is.na(zKurtosis) == FALSE, is.na(intensityKurtosis) == FALSE) %>% # ranger prediction requires complete cases
    mutate(zNormalizedEntropy = replace_na(zNormalizedEntropy, 0), # keep synchronized with species.R
           zMaxNormalized = zMax - zGroundMean,
           zMeanNormalized = zMean - zGroundMean,
           zQ05normalized = zQ05 - zGroundMean,
           zQ10normalized = zQ10 - zGroundMean,
           zQ15normalized = zQ15 - zGroundMean,
           zQ20normalized = zQ20 - zGroundMean,
           zQ25normalized = zQ25 - zGroundMean,
           zQ30normalized = zQ30 - zGroundMean,
           zQ35normalized = zQ35 - zGroundMean,
           zQ40normalized = zQ40 - zGroundMean,
           zQ45normalized = zQ45 - zGroundMean,
           zQ50normalized = zQ50 - zGroundMean,
           zQ55normalized = zQ55 - zGroundMean,
           zQ60normalized = zQ60 - zGroundMean,
           zQ65normalized = zQ65 - zGroundMean,
           zQ70normalized = zQ70 - zGroundMean,
           zQ75normalized = zQ75 - zGroundMean,
           zQ80normalized = zQ80 - zGroundMean,
           zQ85normalized = zQ85 - zGroundMean,
           zQ90normalized = zQ90 - zGroundMean,
           zQ95normalized = zQ95 - zGroundMean) %>%
    select(-zQ05, -zQ10, -zQ15, -zQ20, -zQ25, -zQ30, -zQ35, -zQ40, -zQ45, -zQ50, -zQ55, -zQ60, -zQ65, -zQ70, -zQ75, -zQ80, -zQ85, -zQ90, -zQ95, -zMax, -zMean, -zPcumulative10) %>% # 56 predictor variables
    select(-intensityMax, -intensitySkew, -pCumulativeZQ10, -pCumulativeZQ30, -pCumulativeZQ50, -pCumulativeZQ70, -pCumulativeZQ90, -pFifthReturn, -zGroundMean, -zPcumulative20, -zPcumulative30, -zNormalizedEntropy) # 44 predictor variables
    #select(-intensityKurtosis, -intensityStdDev, -pFourthReturn, -pZaboveZmean, -zPcumulative40, -zPcumulative50, -zPcumulative60, -zPcumulative70, -zPcumulative80, -zPcumulative90, -zQ05normalized, -zSkew) # 32 predictor variables
    #select(-intensityMean, -pFirstReturn, -pSecondReturn, -pThirdReturn, -zQ10normalized, -zQ75normalized, -zQ80normalized, -zQ85normalized, -zQ90normalized, -zQ95normalized, -zMaxNormalized, -zStdDev)
    #select(cellID, r, g, b, nir, pZaboveThreshold, intensityPground, pGround, zQ20normalized, zQ30normalized, zQ40normalized)
    #select(cellID, r, g, b, nir)
  
  tileClassification = predict(randomForestFit$finalModel, data = tileStatisticsTibble) # ~2.1 m, predict(randomForestFit, newdata = ...) fails on internal matrix datatype error

  tileClassificationTibble = left_join(tibble(cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))),
                                       tibble(cellID = tileStatisticsTibble$cellID, classification = tileClassification$predictions),
                                       by = "cellID")
  
  classificationRaster = rast(tileStatistics, nlyrs = 1, names = "classification", vals = tileClassificationTibble$classification)
  writeRaster(classificationRaster, file.path(dataDestinationPath, "classification 10 m non-normalized 44 var", tileName), datatype = "INT1U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
}

cat(paste0("Orthoimage vegetation classification over ", length(tileNames), " tiles in ", format(Sys.time() - jobStartTime), ".\n"))
