#install.packages(c("caret", "ranger", "sf", "terra", "tidyterra", "xgboost"))
#library(doParallel)
library(dplyr)
library(caret)
library(magrittr)
library(terra)
library(tidyr)


## assemble training data
trainingPolygons = vect(file.path(getwd(), "GIS/Trees/segmentation/SpeciesITC.gpkg"), layer = "training polygons + height") # manually designated ground truth polygons
#trainingPolygons = subset(trainingPolygons, trainingPolygons$classification != "shadow") # drop shadow as it's not a cover type

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
#chm = rast(file.path(dataSourcePath, "CHM", "CHM.vrt"))
orthoimagery = rast(file.path(dataSourcePath, "orthoimage v2 16", "orthoimage.vrt"))
orthoimageCellSize = max(res(orthoimagery))
gridMetrics = rast(file.path(dataSourcePath, "metrics", "grid metrics 10 m non-normalized.tif")) # TODO: once recalced, disambiguate intensity naming
names(gridMetrics)[45] = "intensityPground" # temporary workaround for C# typo

trainingData = list()
for (trainingPolygonIndex in 1:nrow(trainingPolygons))
{
  if (trainingPolygonIndex %in% seq(316, 320)) # skip the few training polygons outside of the 400 m grid metrics buffer
  {
    # terra::intersect() is slow and there doesn't appear to be a way avoiding no intersection errors from crop()
    next
  }
  if (trainingPolygonIndex %% 25 == 0)
  {
    cat(paste0(trainingPolygonIndex, "... "))
  }
  
  # extract training data for this polygon
  # For now, assume orthoimagery and CHM are in the same CRS and that orthoimagery is higher resolution.
  # For now, assume all CHM cells have Z values as this is the case for the current set of training polygons (October 2023, n = 365)
  trainingPolygon = trainingPolygons[trainingPolygonIndex]
  if (relate(gridMetrics, trainingPolygon, relation = "contains") == FALSE)
  {
    next # skip two bare ground training polygons south of the grid metrics area
  }
  polygonOrthoimage = crop(orthoimagery, trainingPolygon, touches = FALSE) # take only cells with centroids within the training polygon
  polygonGridMetrics = resample(crop(gridMetrics, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage) # defaults to bilinear for rasters whose first layer is numeric
  
  polygonOrthoimageCrs = crs(polygonOrthoimage, describe = TRUE)
  polygonGridMetricsCrs = crs(polygonGridMetrics, describe = TRUE)
  if ((polygonOrthoimageCrs$authority == polygonGridMetricsCrs$authority) & (polygonOrthoimageCrs$code == polygonGridMetricsCrs$code))
  {
    crs(polygonGridMetrics) = crs(polygonOrthoimage) # give grid metrics the orthoimage's vertical CRS
  }
  
  trainingStatistics = as_tibble(c(polygonOrthoimage, polygonGridMetrics)) %>% # stack rasters
    mutate(polygon = trainingPolygonIndex,
           classification = trainingPolygon$classification)
  #polygonChm = resample(crop(chm, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage)
  trainingData[[trainingPolygonIndex]] = trainingStatistics
}

trainingData = bind_rows(trainingData) %>%
  filter(firstReturns > 0, is.na(n) == FALSE) %>% # drops image pixels with no hits and unpopulated metrics grid cells
  filter(is.na(zGroundMean) == FALSE) %>% # TODO; for now, drop metrics cells without any ground hits
  select(-areaOfPointBoundingBox, -n, -intensityTotal) %>% # drop nonstructural grid metrics
  mutate(classification = as.factor(classification),
         intensitySecondReturn = if_else(intensitySecondReturn == 65535, 0, intensitySecondReturn), # replace NAs with 0
         greenness = green / red, # greenness, NDGR, MGRV @ 98+% correlation, 96% to WBI
         redBlueRatio = red / blue, # rNormalized, WBI, NDGR 98+% correlation, coloration, bNormalized 95%
         nirBlueRatio = nir / blue, # bNDVI, BAI 98%, GEMI 94%
         nirGreenRatio = nir / green, # gNDVI 97%
         nirRedRatio = nir / red, # SAVI 94%, eviBackup 96%
         normalizedGreen = green / (nir + red + green),
         normalizedNir = nir / (nir + red + green),
         normalizedBlue	= blue / (nir + red + green), # NDVI, gNDVI, TVI 95+%
         coloration	= (red - blue) / red, # NDGB, WBI, bNormalized 98+% correlation
         chlorophyllGreen = nir / green - 1, # gNDVI 97%
         chlorophyllVegetation = nir * as.numeric(red) / green^2, # bNDVI 96%, convert explicitly to avoid integer overflow
         ndvi = (nir - red) / (nir + red),
         gndvi = (nir - green) / (nir + green),
         bndvi = (nir - blue) / (nir + blue),
         sipi = (nir - blue) / (nir + red),
         #bai = (blue - nir) / (blue + nir), # bai = -bndvi
         ndgb = (green - blue) / (green + blue),
         ndgr = (green - red) / (green + red),
         wbi = (blue - red) / (blue + red),
         mgrv = (green^2 - red^2) / (green^2 + red^2),
         rgbv = (green^2 - red * as.numeric(blue)) / (green^2 + red * as.numeric(blue)), # convert explicitly to avoid integer overflow
         gli = (2*green - blue - red) / (2*green + blue + red),
         mexg = 1.62 * green - 0.884 * red - 0.311 * blue,
         luminosity = 0.299 * red + 0.587 * green + 0.114 * blue, # NTSC, ITU BT.610
         rNormalized = red / luminosity, # 99% correlation with redBlueRatio, -97% WBI, 94% luminosity
         gNormalized = green / luminosity,
         bNormalized = blue / luminosity,
         nirNormalized = nir / luminosity, # gNDVI 98%
         luminosity709 = 0.2126 * red + 0.7152 * green + 0.0722 * blue, # ITU BT.709
         gemi = (2 * (nir^2 - red^2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5),
         tvi = sqrt(ndvi + 0.5), # NDVI, MTVI2 99%
         evi = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1), # coefficients from MODIS Vegetation Index User's Guide
         triangularVegIndex = 0.5 * (120 * (nir - green) - 200 * (red - green)), # RGBV 93%
         mtvi2 = 1.5 * (2.5 * (nir - green) - 2.5 * (red - green)) / sqrt((2 * nir + 1)^2 - 6 * nir - 5 * sqrt(red) - 0.5),
         eviBackup = 2.5 * (nir - red) / (nir + 2.4 * red + 1),
         savi = 1.5 * (nir - red) / (nir + red + 0.5), # NDVI 100%, eviBackup, oSAVI 99% correlation
         msavi = (2 * nir + 1 - sqrt(2 * (2 * nir + 1)^2 - 8 * (nir - red))) / 2, # mSAVI 94%
         osavi = 1.16 * (nir - red) / (nir + red + 0.16),
         vari = (green - red) / (green + red - blue), # NaNs
         zNormalizedEntropy = replace_na(zNormalizedEntropy, 0), # might be better to replace NAs with 1?
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
  relocate(polygon, classification)

sum(colSums(is.na(trainingData))) # random forest requires complete cases
cor(as.matrix(trainingData %>% select(ndgr, ndvi, gndvi, bndvi, intensity, intensitySecondReturn)))
#cor(as.matrix(trainingData %>% select(-polygon, -classification)))[, "ndvi"]

filterVarImp(x = trainingData %>% select(-classification, -polygon), y = trainingData$classification) %>%
  mutate(max = pmax(bare, conifer, hardwood, shadow), mean = 0.25 * (bare + conifer + hardwood + shadow)) %>%
  #mutate(mean = 0.333 * (bare + conifer + hardwood)) %>%
  arrange(desc(max), desc(mean)) %>%
  mutate(var = row_number()) %>%
  relocate(var)

#crossValidationCluster = makePSOCKcluster(10) # ~14 GB DDR in R @ one core per cross validation repeat
#registerDoParallel(crossValidationCluster)

# drop non-normalized heights with low importance relative to normalized ones
# with shadow class
#trainingData %<>% select(-zQ05, -zQ10, -zQ15, -zQ20, -zQ25, -zQ30, -zQ35, -zQ40, -zQ45, -zQ50, -zQ55, -zQ60, -zQ65, -zQ70, -zQ75, -zQ80, -zQ85, -zQ90, -zQ95, -zKurtosis, -zMax, -zMean, -zPcumulative10) # mean or max variable importance < ~0.68: κ ~0.9997 with 56 predictor variables
#trainingData %<>% select(-intensityMax, -intensitySkew, -pCumulativeZQ10, -pCumulativeZQ30, -pCumulativeZQ50, -pCumulativeZQ70, -pCumulativeZQ90, -pFifthReturn, -zGroundMean, -zPcumulative20, -zPcumulative30, -zNormalizedEntropy) # 44 predictor variables: κ ~0.9999
#trainingData %<>% select(-intensityKurtosis, -intensityStdDev, -pFourthReturn, -pZaboveZmean, -zPcumulative40, -zPcumulative50, -zPcumulative60, -zPcumulative70, -zPcumulative80, -zPcumulative90, -zQ05normalized, -zSkew) # 32 predictor variables (all with max variable importance < 0.90 dropped), no shadow κ = 1
#trainingData %<>% select(-intensityMean, -pFirstReturn, -pSecondReturn, -pThirdReturn, -zQ10normalized, -zQ75normalized, -zQ80normalized, -zQ85normalized, -zQ90normalized, -zQ95normalized, -zMaxNormalized, -zStdDev) # 20 variables: κ ~0.9998 but degraded predictive performance on actual landscape
#trainingData %<>% select(-zMeanNormalized, -zQ15normalized, -zQ25normalized, -zQ35normalized, -zQ45normalized, -zQ50normalized, -zQ55normalized, -zQ60normalized, -zQ65normalized, -zQ70normalized) # 10 variables: κ ~0.9993 but also poor predictive performance
#trainingData %<>% select(-zQ20normalized, -zQ30normalized, -zQ40normalized) # 7 variables (RDB+NIR, intensityPground, pGround, pZaboveThreshold): κ ~0.986, importance red = 100, green = 41, nir = 39, pGround = 5.8, intensityPground = 2.9, blue = 1.8, pZaboveThreshold = 0.0
#trainingData %<>% select(-pZaboveThreshold) # κ ~0.973, red = 100, nir = 34, green = 20, pGround = 6.8, intensityPground = 6.0, blue = 0.0
#trainingData %<>% select(-intensityPground, -pGround) # reduce to RGB+NIR control: κ ~0.8913
#trainingData %<>% select(classification, polygon, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized, ndgb, savi, eviBackup, osavi, nirRedRatio, ndvi, tvi, mtvi2, msavi, normalizedNir, rgbv, triangularVegIndex, nirNormalized, gemi, mexg, nirGreenRatio, chlorophyllGreen, gndvi, normalizedBlue, vari, nirBlueRatio, bndvi, bai, normalizedGreen, gli, evi, chlorophyllVegetation) # 35 normalized vegetation indices: κ ~
#trainingData %<>% select(classification, polygon, luminosity, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized, ndgb, savi, eviBackup, osavi, nirRedRatio, ndvi, tvi, mtvi2, msavi, normalizedNir, rgbv, triangularVegIndex, nirNormalized) # luminosity + 20 normalized vegetation indices: not run due to slow fitting
#trainingData %<>% select(classification, polygon, luminosity, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized) # luminosity + 10 normalized vegetation indices: κ ~ 
#predictorVariables = c("classification", "polygon", "luminosity", "rNormalized", "greenness", "mgrv", "ndgr", "redBlueRatio") # luminosity + 5 normalized vegetation indices: κ ~0.872, 55 minutes to fit
#predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "msavi", "rgbv", "gemi", "mexg", "normalizedGreen", "gli", "evi") # κ ~0.905, 47 minutes to fit
#predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "ndgb", "ndgr", "wbi") # κ ~0.886, 58 minutes to fit
#predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "ndgb", "ndgr") # κ ~0.886, 60 minutes to fit
predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.898, 14 minutes to fit (8 without RGB+NIR no data)
#predictorVariables = c("classification", "ndgr", "ndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.896, 13 minutes to fit
#predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.893, 7 minutes to fit
#predictorVariables = c("classification", "ndvi", "gndvi", "bndvi") # κ ~0.884, 42 minutes to fit
#predictorVariables = c("classification", "rNormalized", "gNormalized", "bNormalized", "nirNormalized") # κ ~0.885, 58 minutes to fit

#trainingData %<>% select(-intensityMax, -pCumulativeZQ10, -pCumulativeZQ30, -pCumulativeZQ50, -pCumulativeZQ90, -zGroundMean, -zNormalizedEntropy, -zPcumulative20) # additional drops: increase κ from ~0.9997 to ~0.9999
# without shadow class
#trainingData %<>% select(-nir, -pCumulativeZQ10, -zQ05, -zQ10, -zQ15, -zQ20, -zQ25, -zQ30, -zQ35, -zQ40, -zQ45, -zQ50, -zQ55, -zQ60, -zQ65, -zQ70, -zQ75, -zQ80, -zQ85, -zQ90, -zQ95, -zGroundMean, -zKurtosis, -zMax, -zMean, -zNormalizedEntropy, -zPcumulative10)

repeatedCrossValidation = trainControl(method = "repeatedcv", number = 2, repeats = 50, verboseIter = TRUE)
#trainingDataSubset = sample_n(trainingData, 10000)
#trainingDataSubset2 = trainingDataSubset %>% select(-zNormalizedEntropy, -pCumulativeZQ90, -pCumulativeZQ50, -pFifthReturn, -intensitySkew, -intensityStdDev, -pCumulativeZQ30, -pCumulativeZQ70, -zPcumulative10, -intensityMax, -zKurtosis, -pCumulativeZQ10, -intensityTotal)

# random forest
# cells   method      cross validation   node sizes   cores   fit time    accuracy    κ       grid metrics
# 103k    ranger      2x50               5            1       27 m        0.999       0.999   from normalized point clouds
# 104k    ranger      2x50               5            1       41 m        0.999       0.999   from non-normalized clouds with heights relative to mean ground elevation added
# 104k    ranger      2x50               5            1       37 m        0.999       0.999   ibid with low importance non-normalized z statistics dropped
fitStart = Sys.time()
randomForestFit = train(classification ~ ., data = trainingData %>% select(predictorVariables), method = "ranger", trControl = repeatedCrossValidation, # importance = "impurity_corrected", 
                        tuneGrid = expand.grid(mtry = floor(sqrt(length(predictorVariables) - 1)),
                                               splitrule = 'gini',
                                               min.node.size = c(8, 10, 12)))
randomForestFitTime = Sys.time() - fitStart # ~25 minutes with shadow, ~15 minutes @ RGB+NIR+16 grid metrics, ~14 minutes @ RGB+NIR+6 grid metrics
#varImp(randomForestFit)
#save(randomForestFit, file = file.path(getwd(), "trees/segmentation/classificationRandomForestFit NDGR+rgbNDVI+intensity12.Rdata"))
load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit 20 m grid metrics.Rdata"))

# linear SVM
# κmax = 0.905 with polygon statistics and 10x10 cross validation
# κmax = 0.999 with all grid metrics
#
# cells   SVM      cross validation   C values   cores   fit time    accuracy    κ       cost scaling
#  10k    linear   10x1               5          1       16.2 s        
#  20k    linear   10x1               5          1       39.6 s                          1.2x
#  50k    linear   10x1               5          1       149 s                           1.5x
# 100k    linear   10x1               5          1       630 s                           2.1x
# 100k    radial   10x1               5          1       318 s
#
# 100k    linear   2x50               6          1       1.3 h        0.999       0.999
# 100k    radial   2x50               6          1       15.1 h       0.999       0.999
fitStart = Sys.time() # <10 s single 5950X core with 114 statistics on 365 polygons, with 60 statistics on 103 k raster cells at 10 cores
svmFitLinear = train(classification ~ ., data = trainingData %>% select(-polygon), method = "svmLinear", # linear kernel
                     trControl = repeatedCrossValidation, preProcess = c("center", "scale"), 
                     tuneGrid = data.frame(C = c(3.8, 4.0, 4.2, 4.4, 4.6, 4.8))) # c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1) for polygon statistics
svmFitTimeLinear = Sys.time() - fitStart

fitStart = Sys.time()
svmFitRadial = train(classification ~ ., data = trainingData %>% select(-polygon), method = "svmRadial", # radial basis functions
                     trControl = repeatedCrossValidation, preProcess = c("center", "scale"), 
                     tuneGrid = expand.grid(C = c(3.8, 4.0, 4.2, 4.4, 4.6, 4.8), sigma = c(0.029, 0.032, 0.035, 0.038)))
svmFitTimeRadial = Sys.time() - fitStart # κmax = 0.935 with polygon statistics and 10x10 cross validation

#save(svmFitLinear, svmFitRadial, file = file.path(getwd(), "trees/segmentation/classificationSvmFits 20 m grid metrics.Rdata"))
#load(file.path(getwd(), "trees/segmentation/classificationSvmFits 20 m grid metrics.Rdata"))
#stopCluster(crossValidationCluster)

# gradient boosting
# cells   method      cross validation   nrounds         cores   fit time    accuracy    κ       cost scaling
# 103k    linear      2x50               50, 100, 150    1       2.5 h       0.999       0.999
# 103k    tree        2x50               50, 100, 150    1       3.4 h       0.999       0.999
fitStart = Sys.time()
xgbLinearFit = train(classification ~., data = trainingData %>% select(-polygon), method = "xgbLinear", trControl = repeatedCrossValidation)
xgbLinearFitTime = Sys.time() - fitStart

fitStart = Sys.time()
xgbTreeFit = train(classification ~., data = trainingData %>% select(-polygon), method = "xgbTree", trControl = repeatedCrossValidation)
xgbTreeFitTime = Sys.time() - fitStart

#save(xgbLinearFit, xgbTreeFit, file = file.path(getwd(), "trees/segmentation/classificationXgboostFits 20 m grid metrics.Rdata"))
load(file.path(getwd(), "trees/segmentation/classificationXgboostFits 20 m grid metrics.Rdata"))

# neural network
# Caret 6.0-94 can fail fairly quickly with "invalid argument to unary operator" if training data is collinear with response 
# variable. If not, fitting still fails with the same error at the end of cross validation.
#fitStart = Sys.time()
#neuralNetworkFit =- train(classification ~ ., data = trainingData %>% select(-polygon), method = "nnet", 
#                          trControl = repeatedCrossValidation, preProcess = c("center", "scale"), verbose = FALSE, trace = FALSE,
#                          tuneGrid = expand.grid(size = c(2, 3), decay = seq(0.05, 1, length.out = 5)))
#neuralNetworkFitTime = Sys.time() - fitStart
#save(neuralNetworkFit, file = file.path(getwd(), "trees/segmentation/neuralNetworkFit 20 m grid metrics.Rdata"))


## forest level classification
orthoimageTiles = list.files(file.path(dataSourcePath, "orthoimage"), "\\.tif")
for (tileFileName in orthoimageTiles)
{
  cat(paste0(tileFileName, "...\n"))

  # extract training data for this polygon
  # For now, assume orthoimagery and CHM are in the same CRS and that orthoimagery is higher resolution.
  # For now, assume all CHM cells have Z values as this is the case for the current set of training polygons (October 2023, n = 365)
  tileOrthoimage = rast(file.path(dataSourcePath, "orthoimage", tileFileName))
  tileGridMetrics = resample(gridMetrics, tileOrthoimage) # ~18 s
  tileStatistics = c(tileOrthoimage, tileGridMetrics)
  tileStatistics$cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))
  tileStatisticsTibble = as_tibble(tileStatistics) %>% # ~15 s
    filter(is.na(red) == FALSE, is.na(green) == FALSE, is.na(n) == FALSE, is.na(zNormalizedEntropy) == FALSE) # ranger prediction requires complete cases

  predictStart = Sys.time()
  tileClassifications = predict(randomForestFit$finalModel, data = tileStatisticsTibble) # predict(randomForestFit, newdata = ...) fails on internal matrix datatype error
  predictTime = Sys.time() - predictStart # ~2.2 m/tile
  
  tileClassificationTibble = left_join(tibble(cellID = 1:(nrow(tileStatistics) * ncol(tileStatistics))),
                                       tibble(cellID = tileStatisticsTibble$cellID, species = tileClassifications$predictions),
                                       by = "cellID")
  
  speciesRaster = rast(tileStatistics, nlyrs = 1, names = "species", vals = tileClassificationTibble$species)
  writeRaster(speciesRaster, file.path(dataSourcePath, "species", tileFileName))
}


## one time setup
# Since gdalbuildvrt doesn't flow metadata (https://github.com/OSGeo/gdal/issues/3627) VRTs need manual editing
# after creation to add <VRTRasterBand/Description> elements naming raster bands. This is needed to set CHM band
# name to Z and orthoimage names to red, green, blue, nir, intensity, intensitySecondReturn, firstReturns, and 
# secondReturns.
# create .vrt for CHM after chmJob.R has processed all tiles
chmSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/CHM"
chmFilePaths = file.path(chmSourcePath, list.files(orthoimageSourcePath, "\\.tif$"))
vrt(chmFilePaths, file.path(chmSourcePath, "CHM.vrt"), overwrite = TRUE)

# create .vrt for orthoimages after orthoImageJob.R has processed all tiles
orthoimageSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2"
orthoimageFilePaths = file.path(orthoimageSourcePath, list.files(orthoimageSourcePath, "\\.tif$"))
vrt(orthoimageFilePaths, file.path(orthoimageSourcePath, "orthoimage.vrt"), overwrite = TRUE)

# create .vrt for species classifications after speciesJob.R has processed all tiles
speciesSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/classification 10 m non-normalized 44 var"
speciesFilePaths = file.path(speciesSourcePath, list.files(speciesSourcePath, "\\.tif$"))
vrt(speciesFilePaths, file.path(speciesSourcePath, "classification.vrt"), overwrite = TRUE)

# recode 32 bit signed orthoimages as 16 bit unsigned
orthoimageDestinationPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16"
for (orthoimageSourcePath in orthoimageFilePaths)
{
  tileName = basename(orthoimageSourcePath)
  orthoimageDestinationFilePath = file.path(orthoimageDestinationPath, tileName)
  if (file.exists(orthoimageDestinationFilePath) == FALSE)
  {
    orthoimage = rast(orthoimageSourcePath)
    writeRaster(orthoimage, orthoimageDestinationFilePath, datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
  }
}


## tile checking
orthoimageTile = rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16/s03480w06540.tif")
as_tibble(orthoimageTile) %>% filter(firstReturns > 0, is.na(green))
