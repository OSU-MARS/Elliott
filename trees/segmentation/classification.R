#install.packages(c("caret", "ggcorrplot", "ranger", "sf", "terra", "tidyterra", "xgboost"))
library(dplyr)
library(caret)
library(FNN)
library(forcats)
library(magrittr)
library(terra)
library(tidyterra)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), panel.border = element_blank()))

classificationOptions = tibble(rebuildTrainingData = FALSE, includeExploratory = FALSE)

## assemble training data
# Training data setup takes just over an hour for ~2000 polygons. Can't multithread as terra 1.7 lacks support.
if (classificationOptions$rebuildTrainingData)
{
  trainingPolygons = vect(file.path(getwd(), "GIS/Trees/segmentation/classification training polygons.gpkg"), layer = "training polygons + height") # manually designated ground truth polygons
  
  dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
  dsm = rast(file.path(dataSourcePath, "DSM v3 beta", "dsm cmm3 chm aerialMean.vrt"))
  orthoimage = rast(file.path(dataSourcePath, "orthoimage v3", "orthoimage all.vrt")) # red, green, blue, nearInfrared, intensityFirstReturn, intensitySecondReturn, firstReturns, secondReturns, scanAngleMeanAbsolute
  orthoimageCellSize = max(res(orthoimage))
  imageCenters = vect("GIS/DOGAMI/2021 OLC Coos County/image positions.gpkg") # terra drops z coordinates but they're redundant with the elevation field
  gridMetrics = rast(file.path(dataSourcePath, "metrics", "grid metrics 10 m non-normalized v2.tif")) %>%
    rename(intensityFirstGridReturn = intensityFirstReturn)

  dsmCrs = crs(dsm, describe = TRUE)
  orthoimageCrs = crs(orthoimage, describe = TRUE)
  gridMetricsCrs = crs(gridMetrics, describe = TRUE)
  imageCentersCrs = crs(imageCenters, describe = TRUE)
  # check CRSes for horizontal consistency
  if ((orthoimageCrs$authority != dsmCrs$authority) | (orthoimageCrs$code != dsmCrs$code) |
      (orthoimageCrs$authority != gridMetricsCrs$authority) | (gridMetricsCrs$code != dsmCrs$code) |
      (orthoimageCrs$authority != imageCentersCrs$authority) | (imageCentersCrs$code != dsmCrs$code))
  {
    stop(paste0("DSM, orthoimage, and grid metrics CRSes do not match for training polygon ", trainingPolygonIndex))
  }
  # fix up vertical CRSes
  if ((orthoimageCrs$authority == gridMetricsCrs$authority) & (orthoimageCrs$code == gridMetricsCrs$code))
  {
    crs(dsm) = crs(orthoimage) # give DSM the orthoimage's vertical CRS, TODO: remove once DSM fix propagates
  }
  if ((orthoimageCrs$authority == gridMetricsCrs$authority) & (orthoimageCrs$code == gridMetricsCrs$code))
  {
    crs(gridMetrics) = crs(orthoimage) # give grid metrics the orthoimage's vertical CRS
  }

  imageCentersForKnn = tibble(xy = geom(imageCenters), z = imageCenters$`elevation, ft`) %>% # flatten to coordinates for kNN
    mutate(x = xy[, "x"], y = xy[, "y"]) %>% select(-xy) %>% relocate(x, y, z)
  imageSunPositions = tibble(azimuth = imageCenters$sunAzimuth, elevation = imageCenters$sunElevation) # flatten for in loop lookup

  polygonCreationStart = Sys.time()
  cat(paste0(nrow(trainingPolygons), " training polygons: 0... "))
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
    # Assumes
    # - training polygons are small so the error of using constant sun and view angles for all cells in the polygon is slight
    # - imaging flight altitude is high enough the mean DSM elevation is an acceptably accurate representation of the polygon
    # - all layers are in the same CRS per checks above and that orthoimagery is of equal or higher resolution than the DSM
    trainingPolygon = trainingPolygons[trainingPolygonIndex]
    if (relate(gridMetrics, trainingPolygon, relation = "contains") == FALSE)
    {
      next # skip two bare ground training polygons south of the grid metrics area
    }
    polygonDsm = crop(dsm, trainingPolygon, touches = FALSE) # take only cells with centroids within the training polygon
    polygonOrthoimage = crop(orthoimage, trainingPolygon, touches = FALSE)
    polygonGridMetrics = resample(crop(gridMetrics, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage, method = "lanczos") # defaults to bilinear for rasters whose first layer is numeric
    
    polygonCentroid = matrix(c(geom(centroids(trainingPolygon))[, c("x", "y")], z = mean(values(polygonDsm$dsm, na.rm = TRUE))), nrow = 1) # x, y, z
    
    imageCenterIndex = knnx.index(imageCentersForKnn, polygonCentroid, k = 1)
    imagePositionDelta = imageCentersForKnn[imageCenterIndex, ] - polygonCentroid
    polygonViewAzimuth = -180/pi * (atan2(imagePositionDelta$y, imagePositionDelta$x) - pi/2) # -180/pi * (atan2(c(1, -1, -1, 1), c(1, 1, -1, -1)) - pi/2) to rotate to 0° being north and increasing clockwise 
    polygonViewAzimuth = if_else(polygonViewAzimuth > 0, polygonViewAzimuth, 360 + polygonViewAzimuth)
    polygonZenithAngle = 180/pi * atan2(sqrt(imagePositionDelta$x^2 + imagePositionDelta$y^2), imagePositionDelta$z)
    
    trainingStatistics = as_tibble(c(polygonOrthoimage, polygonDsm, polygonGridMetrics)) %>% # stack rasters
      mutate(polygon = trainingPolygonIndex,
             classification = trainingPolygon$classification,
             sunAzimuth = imageSunPositions$azimuth[imageCenterIndex],
             sunElevation = imageSunPositions$elevation[imageCenterIndex],
             viewAzimuth = polygonViewAzimuth,
             viewElevation = 90 - polygonZenithAngle)
    #polygonChm = resample(crop(chm, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage)
    trainingData[[trainingPolygonIndex]] = trainingStatistics
  }
  (Sys.time() - polygonCreationStart)
  
  trainingData = bind_rows(trainingData) %>%
    filter(firstReturns > 0) %>% # drops image pixels with no RGB[+NIR] data
    rename(nir = nearInfrared) %>%
    mutate(classification = as.factor(classification),
           intensitySecondReturn = replace_na(intensitySecondReturn, 0)) %>% # allow cells without second returns to be classified
    relocate(polygon, classification)
  
  # ranger requires complete cases but not all raster cells in all training polygons have all grid metrics or point hits
  # While data is written here for all training cells containing LiDAR points, it's useful to check what filters are required below.
  #print(trainingData %>% filter(is.na(dsm) == FALSE, is.na(zMean) == FALSE) %>% summarise(across(everything(), ~sum(is.na(.)))), width = Inf)
  #sum(is.na(trainingData))
  saveRDS(trainingData, "trees/segmentation/classification training data.Rds")
} else {
  trainingData = readRDS("trees/segmentation/classification training data.Rds") %>%
    filter(is.na(dsm) == FALSE, is.na(zMean) == FALSE) # only complete cases for now: exclude polygons not covered by grid metrics
}

# remove unique polygon identifiers and add derived predictors
trainingData = trainingData %>% 
  select(-polygon) %>%
  mutate(sunRelativeAzimuth = sunAzimuth - viewAzimuth, # 0 = forward scatter, ±180 = backscatter, absolute value broken out separately to allow for asymmetric BRDF
         sunRelativeAzimuth = if_else(sunRelativeAzimuth > 180, 360 - sunRelativeAzimuth, if_else(sunRelativeAzimuth < -180, 360 + sunRelativeAzimuth, sunRelativeAzimuth)), # clamp to [0, ±180] to constrain training complexity
         sunRelativeAbsoluteAzimuth = abs(sunRelativeAzimuth),
         sunRelativeAzimuthCosine = cos(pi/180 * sunRelativeAzimuth),
         sunRelativeAzimuthSine = sin(pi/180 * sunRelativeAzimuth),
         scanAngleCosine = cos(pi/180 * scanAngleMeanAbsolute),
         sunZenithAngleCosine = cos(pi/180 * (90 - sunElevation)),
         viewZenithAngleCosine = cos(pi/180 * (90 - viewElevation)),
         # vegetation indices
         ari = 1/green - 1/red, # WDRVI
         arvi2 = (nir - (red - 2 * (blue - red))) / (nir + (red - 2 * (blue - red))), # arvi(γ = 0) = ndvi, arvi(γ = 1) = bndvi
         atsavi = 1.22 * (nir - 1.22 * red - blue) / (1.22 * nir + red - 0.23567),
         #bai = (blue - nir) / (blue + nir), # bai = -bndvi
         bgi = blue / green,
         bri = blue / red,
         bndvi = (nir - blue) / (nir + blue), # bndvi = psnd
         brvi = (nir / (green + 0.1 * red) - blue / (red + 0.5 * green)) / (nir / (green + 0.1 * red) + blue / (red + 0.5 * green)),
         coloration	= (red - blue) / red, # NDGB, WBI, bNormalized 98+% correlation
         chlorophyllGreen = nir / green - 1, # gNDVI 97%
         chlorophyllVegetation = nir * as.numeric(red) / green^2, # bNDVI 96%, convert explicitly to avoid integer overflow
         evi = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1), # coefficients from MODIS Vegetation Index User's Guide
         eviBackup = 2.5 * (nir - red) / (nir + 2.4 * red + 1), # backup evi = evi2
         gari = (nir - (green - 1.7 * (blue - red))) / (nir + (green - 1.7 * (blue - red))), # also GRARI but η and λ not defined by Gitelson et al. 1996
         gdvi = ((nir/red) - 1) / ((nir/red) + 1),
         gdvi2 = ((nir/red)^2- 1) / ((nir/red)^2 + 1), # also GDVI3 and 4 but even 2 not recommended for forested areas
         gemi = (2 * (nir^2 - red^2) + 1.5 * nir + 0.5 * red) / (nir + red + 0.5),
         gli = (2*green - blue - red) / (2*green + blue + red),
         gndvi = (nir - green) / (nir + green),
         greenness = green / red, # greenness, NDGR, MGRV @ 98+% correlation, 96% to WBI
         # grvi = ndgr,
         luminosity = 0.299 * red + 0.587 * green + 0.114 * blue, # NTSC, ITU BT.610
         luminosity709 = 0.2126 * red + 0.7152 * green + 0.0722 * blue, # ITU BT.709
         mari = nir * (1/green - 1/red),
         mexg = 1.62 * green - 0.884 * red - 0.311 * blue,
         mcari = 1.5 * (2.5 * (nir - red) - 1.3 * (nir - green)) / sqrt((2 * nir + 1)^2 - (6 * nir - 5 * red) - 0.5),
         mgrv = (green^2 - red^2) / (green^2 + red^2),
         msavi = (2 * nir + 1 - sqrt(2 * (2 * nir + 1)^2 - 8 * (nir - red))) / 2, # mSAVI 94%
         msr = (nir/red - 1) / sqrt(nir/red + 1),
         mtvi2 = 1.5 * (2.5 * (nir - green) - 2.5 * (red - green)) / sqrt((2 * nir + 1)^2 - 6 * nir - 5 * sqrt(red) - 0.5),
         nirBlueRatio = nir / blue, # bNDVI, BAI 98%, GEMI 94%
         nirGreenRatio = nir / green, # nirGreenRatio = pbi, gNDVI 97%
         nirRedRatio = nir / red, # nirRedRatio = pssr, SAVI 94%, eviBackup 96%
         nirNormalized = nir / luminosity, # gNDVI 98%
         normalizedBlue	= blue / (nir + red + green), # NDVI, gNDVI, TVI 95+%
         normalizedGreen = green / (nir + red + green),
         normalizedNir = nir / (nir + red + green),
         ndgb = (green - blue) / (green + blue),
         ndgr = (green - red) / (green + red), # same as GRVI
         ndvi = (nir - red) / (nir + red),
         osavi = 1.16 * (nir - red) / (nir + red + 0.16),
         redBlueRatio = red / blue, # rNormalized, WBI, NDGR 98+% correlation, coloration, bNormalized 95%
         rdvi = (nir - red) / sqrt(nir + red),
         rgbv = (green^2 - red * as.numeric(blue)) / (green^2 + red * as.numeric(blue)), # convert explicitly to avoid integer overflow
         rndvi = (nir - red) / sqrt(nir + red), # sometimes mistaken for rdvi?
         savi = 1.5 * (nir - red) / (nir + red + 0.5), # NDVI 100%, eviBackup, oSAVI 99% correlation
         sipi = (nir - blue) / (nir + red),
         triangularVegIndex = 0.5 * (120 * (nir - green) - 200 * (red - green)), # RGBV 93%
         vari = (green - red) / (green + red - blue), # NaNs
         wbi = (blue - red) / (blue + red),
         # NDVI transforms and other second step image metrics
         ctvi = (ndvi + 0.5) * ndvi + 0.5,
         tvi = sqrt(ndvi + 0.5), # NDVI, MTVI2 99%
         bNormalized = blue / luminosity,
         gNormalized = green / luminosity,
         rNormalized = red / luminosity, # 99% correlation with redBlueRatio, -97% WBI, 94% luminosity
         # LiDAR statistics
         zNormalizedEntropy = replace_na(zNormalizedEntropy, 0), # might be better to replace NAs with 1?
         zMaxNormalized = zMax - zGroundMean,
         zMeanNormalized = zMean - zGroundMean,
         zQ05normalized = zQ05 - zGroundMean, # z quantiles mostly 90-99% correlated, lower quantiles more distinct
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
         zQ95normalized = zQ95 - zGroundMean)

if (classificationOptions$includeExploratory)
{
  library(VSURF)  # thresholding ~2 h trees + permutation + 10 h parallel + interpretation 10 h on 170k cell, 133 variables training set
  classificationVsurf = VSURF(classification ~ ., trainingData, ncores = 12, parallel = TRUE, RFimplem = "ranger")
  saveRDS(classificationVsurf, "trees/segmentation/classification vsurf 152.Rds")
  classificationVsurf$nums.varselect # threshold -> interpretation -> prediction
  classificationVsurf$mean.perf
  classificationVsurf$overall.time
  classificationVsurf$comput.times

  plot(classificationVsurf)
  attributes(classificationVsurf$terms[classificationVsurf$varselect.interp])$term.labels
  attributes(classificationVsurf$terms[classificationVsurf$varselect.pred])$term.labels

  predictorImportance = tibble(responseVariable = "classification", predictor = as.character(attr(classificationVsurf$terms, "predvars"))[classificationVsurf$imp.mean.dec.ind + 2], importance = classificationVsurf$imp.mean.dec) %>% # offset as.character() by two since first element is "list" and second is classification
    mutate(importance = importance / sum(importance))
  ggplot() +
    geom_col(aes(x = importance, y = fct_reorder(predictor, importance), fill = responseVariable), predictorImportance) +
    guides(fill = "none") +
    labs(x = "normalized variable importance", y = NULL)  
  ggsave("trees/segmentation/classification vsurf 152 importance.png", width = 10, height = 40, units = "cm", dpi = 150)
  # orthoimagery correlations
  #cor(as.matrix(trainingData %>% select(ndgr, ndvi, gndvi, bndvi, intensity, intensitySecondReturn)))
  # LiDAR correlations
  #cor(as.matrix(trainingData %>% select(zQ05normalized, zQ10normalized, zQ25normalized, zQ75normalized, zQ90normalized)))
  #cor(as.matrix(trainingData %>% select(intensityPground, intensityQ10, intensity, intensityFirstReturn, intensitySecondReturn)))
  #cor(as.matrix(trainingData %>% select(pFirstReturn, pSecondReturn, pThirdReturn, pFourthReturn, pFifthReturn, pZaboveZmean, pZaboveThreshold)))
  #cor(as.matrix(trainingData %>% select(pGround, pSecondReturn, pThirdReturn, pZaboveThreshold, zQ05normalized, zQ10normalized, zQ25normalized, zQ75normalized, intensityMean, intensityStdDev, zMean, zSkew, pZaboveZmean)))
  #predictorCorrelations = cor(as.matrix(trainingData %>% select(ndgr, ndvi, gndvi, bndvi, intensity, intensityFirstReturn, intensitySecondReturn, intensityQ10, intensityPground, pThirdReturn, pZaboveThreshold, zQ05normalized, zQ10normalized, zQ25normalized, zQ75normalized)))
  #ggcorrplot::ggcorrplot(predictorCorrelations)
  #cor(as.matrix(trainingData %>% select(-polygon, -classification)))[, "ndvi"]
  
  # imagery metric importance: { rNormalized, greenness, ndgr, mgrv, redBlueRatio, coloration, wbi, red, bNoramlized } 
  #                            { ndgb, luminosity, luminosity709, green, gNormalized, blue, mtvi2, nirRedRatio, ndvi, tvi, osavi, eviBackup, msavi, nir, rgbv, nirNormalized }
  #                            { gemi, sipi, mexg, nirGreenRatio, chlorophyllGreen, gndvi, normalizedBlue, vari, nirBlueRatio, bndvi, normalizedGreen, gli, secondReturns }
  #                            { evi, chlorophyllVegetation, intensitySecondReturn, intensity }
  #                    unused: { first returns }
  # LiDAR metric importance: { intensityPground, pZaboveThreshold, zQ40, 35, 45, zMeanNormalized, zQ55, 60, 70, 30, 75,  65, 80, pSecondReturn, zQ80, 25, 75, 85, 90 }
  #                          { pFirstReturn, zQ95, intensityMean, intensityQ30, zMaxNormalized, intensityQ20, intensityMeanBelowMeidanZ, intensityQ40, zQ10, intensityQ70, pThirdReturn, intensityQ50, 10, 80, zStdDev, 60, intensityFirstReturn }
  #                          { zQ05, intensityQ90, intensityStdDev, pFourthReturn, zSkew, pZaboveZmean, pFifthReturn, intensitySkew }
  #                  unused: { intensityMax, zNormalizedEntropy, zPc20, zGroundMean, zQ05, pCuZQ10, zQ10, 15, 20, 25, 30, 35, 40, 45, zPc10, zQ50, zMean, zZQ55, zQ60, zMax, zQ65, 70, 90, 85, 80, 75, 95, zKurtosis }
  filterVarImp(x = trainingData %>% select(-classification, -polygon), y = trainingData$classification) %>%
    mutate(max = pmax(bare, conifer, hardwood, shadow), mean = 0.25 * (bare + conifer + hardwood + shadow)) %>%
    #mutate(mean = 0.333 * (bare + conifer + hardwood)) %>%
    arrange(desc(max), desc(mean)) %>%
    mutate(var = row_number()) %>%
    relocate(var)

  # drop non-normalized heights with low importance relative to normalized ones
  predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.898, 15 minutes to fit with 8-10-12 min node sizes, 8 minutes with 6-8-10
  predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn", "pGround", "pSecondReturn", "pThirdReturn", "pZaboveThreshold", "zQ05normalized", "zQ10normalized", "zQ25normalized", "zQ75normalized", "intensityMean", "intensityStdDev", "zMean", "zSkew", "pZaboveZmean") # lanczos κ ~0.9998 + 13 minutes to fit, bilinear κ ~0.9998 + 14 minutes
  predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn", "pGround", "pSecondReturn", "pThirdReturn", "pZaboveThreshold", "zQ05normalized", "zQ10normalized", "zQ25normalized", "zQ75normalized", "intensityMean", "intensityStdDev") # lanczos κ ~0.9995 + 14 minutes to fit
  predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensity", "intensityFirstReturn", "intensitySecondReturn", "intensityQ10", "intensityPground", "pThirdReturn", "pZaboveThreshold", "zQ05normalized", "zQ10normalized", "zQ25normalized", "zQ75normalized") # lanczos κ ~0.9996 + 14 minutes to fit
  
  #predictorVariables = c("classification", "luminosity", "rNormalized", "greenness", "mgrv", "ndgr", "redBlueRatio") # luminosity + 5 normalized vegetation indices: κ ~0.872, 55 minutes to fit
  #predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "msavi", "rgbv", "gemi", "mexg", "normalizedGreen", "gli", "evi") # κ ~0.905, 47 minutes to fit
  #predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "ndgb", "ndgr", "wbi") # κ ~0.886, 58 minutes to fit
  #predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "ndgb", "ndgr") # κ ~0.886, 60 minutes to fit
  #predictorVariables = c("classification", "ndgr", "ndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.896, 13 minutes to fit
  #predictorVariables = c("classification", "ndvi", "gndvi", "bndvi", "intensity", "intensitySecondReturn") # κ ~0.893, 7 minutes to fit
  #predictorVariables = c("classification", "ndvi", "gndvi", "bndvi") # κ ~0.884, 42 minutes to fit
  #predictorVariables = c("classification", "rNormalized", "gNormalized", "bNormalized", "nirNormalized") # κ ~0.885, 58 minutes to fit
  #trainingData %<>% select(-zQ05, -zQ10, -zQ15, -zQ20, -zQ25, -zQ30, -zQ35, -zQ40, -zQ45, -zQ50, -zQ55, -zQ60, -zQ65, -zQ70, -zQ75, -zQ80, -zQ85, -zQ90, -zQ95, -zKurtosis, -zMax, -zMean, -zPcumulative10) # mean or max variable importance < ~0.68: κ ~0.9997 with 56 predictor variables
  #trainingData %<>% select(-intensityMax, -intensitySkew, -pCumulativeZQ10, -pCumulativeZQ30, -pCumulativeZQ50, -pCumulativeZQ70, -pCumulativeZQ90, -pFifthReturn, -zGroundMean, -zPcumulative20, -zPcumulative30, -zNormalizedEntropy) # 44 predictor variables: κ ~0.9999
  #trainingData %<>% select(-intensityKurtosis, -intensityStdDev, -pFourthReturn, -pZaboveZmean, -zPcumulative40, -zPcumulative50, -zPcumulative60, -zPcumulative70, -zPcumulative80, -zPcumulative90, -zQ05normalized, -zSkew) # 32 predictor variables (all with max variable importance < 0.90 dropped), no shadow κ = 1
  #trainingData %<>% select(-intensityMean, -pFirstReturn, -pSecondReturn, -pThirdReturn, -zQ10normalized, -zQ75normalized, -zQ80normalized, -zQ85normalized, -zQ90normalized, -zQ95normalized, -zMaxNormalized, -zStdDev) # 20 variables: κ ~0.9998 but degraded predictive performance on actual landscape
  #trainingData %<>% select(-zMeanNormalized, -zQ15normalized, -zQ25normalized, -zQ35normalized, -zQ45normalized, -zQ50normalized, -zQ55normalized, -zQ60normalized, -zQ65normalized, -zQ70normalized) # 10 variables: κ ~0.9993 but also poor predictive performance
  #trainingData %<>% select(-zQ20normalized, -zQ30normalized, -zQ40normalized) # 7 variables (RDB+NIR, intensityPground, pGround, pZaboveThreshold): κ ~0.986, importance red = 100, green = 41, nir = 39, pGround = 5.8, intensityPground = 2.9, blue = 1.8, pZaboveThreshold = 0.0
  #trainingData %<>% select(-pZaboveThreshold) # κ ~0.973, red = 100, nir = 34, green = 20, pGround = 6.8, intensityPground = 6.0, blue = 0.0
  #trainingData %<>% select(-intensityPground, -pGround) # reduce to RGB+NIR control: κ ~0.8913
  #trainingData %<>% select(classification, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized, ndgb, savi, eviBackup, osavi, nirRedRatio, ndvi, tvi, mtvi2, msavi, normalizedNir, rgbv, triangularVegIndex, nirNormalized, gemi, mexg, nirGreenRatio, chlorophyllGreen, gndvi, normalizedBlue, vari, nirBlueRatio, bndvi, bai, normalizedGreen, gli, evi, chlorophyllVegetation) # 35 normalized vegetation indices: κ ~
  #trainingData %<>% select(classification, luminosity, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized, ndgb, savi, eviBackup, osavi, nirRedRatio, ndvi, tvi, mtvi2, msavi, normalizedNir, rgbv, triangularVegIndex, nirNormalized) # luminosity + 20 normalized vegetation indices: not run due to slow fitting
  #trainingData %<>% select(classification, luminosity, rNormalized, greenness, mgrv, ndgr, redBlueRatio, coloration, wbi, bNormalized) # luminosity + 10 normalized vegetation indices: κ ~

  #trainingData %<>% select(-intensityMax, -pCumulativeZQ10, -pCumulativeZQ30, -pCumulativeZQ50, -pCumulativeZQ90, -zGroundMean, -zNormalizedEntropy, -zPcumulative20) # additional drops: increase κ from ~0.9997 to ~0.9999
  # without shadow class
  #trainingData %<>% select(-nir, -pCumulativeZQ10, -zQ05, -zQ10, -zQ15, -zQ20, -zQ25, -zQ30, -zQ35, -zQ40, -zQ45, -zQ50, -zQ55, -zQ60, -zQ65, -zQ70, -zQ75, -zQ80, -zQ85, -zQ90, -zQ95, -zGroundMean, -zKurtosis, -zMax, -zMean, -zNormalizedEntropy, -zPcumulative10)
}

repeatedCrossValidation = trainControl(method = "repeatedcv", number = 2, repeats = 25, verboseIter = TRUE)
#trainingDataSubset = sample_n(trainingData, 10000)
#trainingDataSubset2 = trainingDataSubset %>% select(-zNormalizedEntropy, -pCumulativeZQ90, -pCumulativeZQ50, -pFifthReturn, -intensitySkew, -intensityStdDev, -pCumulativeZQ30, -pCumulativeZQ70, -zPcumulative10, -intensityMax, -zKurtosis, -pCumulativeZQ10, -intensityTotal)

# random forest
# 1936 polygons
# cells   method              cross validation   node sizes   cores   fit time    accuracy    κ       grid metrics
# 169k    ranger + vsurf 7    2x25               4            1       25 m        0.944       0.933   10 m from non-normalized clouds with heights relative to mean ground elevation added
# 365 polygons: bare, conifer, hardwood, shadow
# 103k    ranger              2x50               5            1       27 m        0.999       0.999   from normalized point clouds
# 104k    ranger              2x50               5            1       41 m        0.999       0.999   from non-normalized clouds with heights relative to mean ground elevation added
# 104k    ranger              2x50               5            1       37 m        0.999       0.999   ibid with low importance non-normalized z statistics dropped
predictorVariables = c("classification", "zMaxNormalized", "zQ95normalized", "gli", "intensityMeanAboveMedianZ", "zQ75normalized", "intensitySkew", "intensityFirstGridReturn")
#cor(as.matrix(trainingData %>% select(all_of(predictorVariables[2:length(predictorVariables)]))))

fitStart = Sys.time() # 4.3 h @ 138 variables, 169k rows, 2x10 cross validation
randomForestFit = train(classification ~ ., data = trainingData %>% select(all_of(predictorVariables)), method = "ranger", trControl = repeatedCrossValidation, # importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = floor(sqrt(length(predictorVariables) - 1)),
                                               splitrule = 'gini',
                                               min.node.size = c(2, 4, 8, 16))) # @ 138->7
(randomForestFitTime = Sys.time() - fitStart)
save(randomForestFit, file = file.path(getwd(), "trees/segmentation/classificationRandomForestFit vsurf 138.7 10 m lanczos.Rdata"))
randomForestFit

randomForestPrediction = tibble(actual = trainingData$classification, predicted = predict(randomForestFit)) %>%
  mutate(actualClass = fct_collapse(actual, bare = c("bare", "bare shadow"),
                                            conifer = c("conifer", "conifer shadow", "conifer deep shadow"),
                                            hardwood = c("hardwood", "hardwood shadow", "hardwood deep shadow"),
                                            snag = c("brown tree", "grey tree")),
         predictedClass = fct_collapse(predicted, bare = c("bare", "bare shadow"),
                                                  conifer = c("conifer", "conifer shadow", "conifer deep shadow"),
                                                  hardwood = c("hardwood", "hardwood shadow", "hardwood deep shadow"),
                                                  snag = c("brown tree", "grey tree")))
randomForestConfusionMatrix = confusionMatrix(randomForestPrediction$predictedClass, randomForestPrediction$actualClass)
randomForestConfusionMatrix$table
randomForestConfusionMatrix$overall

randomForestImportance = varImp(randomForestFit) # only shows top 20
randomForestImportance = tibble(variable = rownames(randomForestImportance$importance), importance = randomForestImportance$importance$Overall) %>%
  arrange(desc(importance))
print(randomForestImportance, n = nrow(randomForestImportance))

#load(file.path(getwd(), "trees/segmentation/classificationRandomForestFit 20 m grid metrics.Rdata"))

ggplot() +
  geom_col(aes(x = importance, y = forcats::fct_reorder(variable, importance)), randomForestImportance) +
  labs(x = "normalized importance, %", y = NULL)

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


## tile checking
orthoimageTile = rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16/s03480w06540.tif")
as_tibble(orthoimageTile) %>% filter(firstReturns > 0, is.na(green))

tileIndex = as_tibble(vect("GIS/DOGAMI/2021 OLC Coos County/Elliott tile index.gpkg")) %>%
  mutate(hasRgbNir = as.logical(hasRgbNir))

rgbNirTiles = tibble(Tile_ID = str_remove(list.files("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/pointz RGBIR", "*.laz"), ".laz"))
#left_join(tileIndex, rgbNirTiles %>% mutate(isRgbNir = TRUE), by = join_by(Tile_ID)) %>% filter(hasRgbNir == FALSE, isRgbNir)

orthoimageTiles = str_remove(list.files("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16", "*.tif"), ".tif")
notActuallyOrthoImageTiles = str_c(setdiff(orthoimageTiles, rgbNirTiles$Tile_ID), ".tif") # have RGB but (presumably) aren't raytraced
file.rename(file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16", notActuallyOrthoImageTiles),
            file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/rgb", notActuallyOrthoImageTiles))


# recode 32 bit signed orthoimages as 16 bit unsigned
#orthoimageDestinationPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16"
#for (orthoimageSourcePath in orthoimageFilePaths)
#{
#  tileName = basename(orthoimageSourcePath)
#  orthoimageDestinationFilePath = file.path(orthoimageDestinationPath, tileName)
#  if (file.exists(orthoimageDestinationFilePath) == FALSE)
#  {
#    orthoimage = rast(orthoimageSourcePath)
#    writeRaster(orthoimage, orthoimageDestinationFilePath, datatype = "INT2U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"), overwrite = TRUE)
#  }
#}

# noise crop
lasTile = lidR::readLAS("E:/Elliott/GIS/DOGAMI/2021 OLC Coos County/points/s04230w06810.las") # EPSG::6557 + 8228 vertical
areaOfInterest = sf::st_read("GIS/Trees/segmentation/noise points.gpkg", layer = "unit test noise") # EPSG::6557
sf::st_crs(areaOfInterest) = crs(lasTile)
lasAoi = lidR::clip_roi(lasTile, areaOfInterest)
lidR::writeLAS(lasAoi[[1]], "GIS/Trees/segmentation/s04230w06810 main noise zone.las")


## virtual raster generation and manual editing
# Since gdalbuildvrt doesn't flow metadata (https://github.com/OSGeo/gdal/issues/3627) VRTs need manual editing
# after creation to add <VRTRasterBand/Description> elements naming raster bands. This is needed to set band
# name to Z and orthoimage names to red, green, blue, nir, intensity, intensitySecondReturn, firstReturns, and 
# secondReturns.
# create .vrt for DSM after Get-Dsm
library(terra)
dsmSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 alpha"
dsmFilePaths = file.path(dsmSourcePath, list.files(dsmSourcePath, "\\.tif$"))
vrt(dsmFilePaths, file.path(dsmSourcePath, "DSM dsm cmm3 chm terra.vrt"), options = c("b", 1, "b", 2, "b", 3), overwrite = TRUE, set_names = TRUE)

# create .vrt for orthoimages after orthoImageJob.R has processed all tiles
orthoimageSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16"
orthoimageFilePaths = file.path(orthoimageSourcePath, list.files(orthoimageSourcePath, "\\.tif$"))
vrt(orthoimageFilePaths, file.path(orthoimageSourcePath, "orthoimage.vrt"), overwrite = TRUE, set_names = TRUE)

# create .vrt for hardwood-conifer classification after classificationJob.R has processed all tiles
classificationSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/classification ortho6 + 10 m non-normalized 9 lanczos"
classificationFilePaths = file.path(classificationSourcePath, list.files(classificationSourcePath, "\\.tif$"))
vrt(classificationFilePaths, file.path(classificationSourcePath, "classification.vrt"), overwrite = TRUE, set_names = TRUE)
