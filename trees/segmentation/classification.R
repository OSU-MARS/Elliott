#install.packages(c("caret", "ggcorrplot", "ranger", "sf", "terra", "tidyterra", "xgboost"))
library(dplyr)
library(caret)
library(FNN)
library(forcats)
library(magrittr)
library(stringr)
library(terra)
library(tidyterra)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), panel.border = element_blank()))

classificationOptions = tibble(gridMetricsResamplingMethod = "cubic",
                               rangerThreads = 0.5 * future::availableCores(),
                               predictClasses = TRUE,
                               rebuildTrainingData = FALSE, 
                               includeExploratory = FALSE,
                               includeSetup = FALSE)

## assemble training data
# Training data setup takes ~4.30 hours for ~2000 polygons, single threaded on a 9900X. Can't multithread as terra 1.7 lacks support.
if (classificationOptions$rebuildTrainingData)
{
  trainingPolygons = vect(file.path(getwd(), "GIS/Trees/classification training polygons.gpkg"), layer = "training polygons + height") # manually designated ground truth polygons

  dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
  dsm = rast(file.path(dataSourcePath, "DSM v3 beta", "dsm cmm3 chm aerialMean density.vrt"))
  dsmSlopeAspect005 = rast(file.path(dataSourcePath, "DSM v3 beta", "slopeAspect", "slopeAspect.vrt"))
  dtmSlope100 = rast(file.path(dataSourcePath, "bare earth slope Gaussian 10 m EPSG6557.tif")) %>%
    rename(slope100 = `bare earth slope Gaussian 10 m EPSG6557`)
  dtmAspectCos100 = rast(file.path(dataSourcePath, "bare earth cos(aspect) Gaussian 10 m EPSG6557.tif")) %>%
    rename(aspect100cosine = `bare earth cos(aspect) Gaussian 10 m EPSG6557`)
  dtmAspectSin100 = rast(file.path(dataSourcePath, "bare earth sin(aspect) Gaussian 10 m EPSG6557.tif")) %>%
    rename(aspect100sine = `bare earth sin(aspect) Gaussian 10 m EPSG6557`)
  orthoimage = rast(file.path(dataSourcePath, "orthoimage v3", "orthoimage all.vrt")) # red, green, blue, nearInfrared, intensityFirstReturn, intensitySecondReturn, firstReturns, secondReturns, scanAngleMeanAbsolute
  imageCenters = vect("GIS/DOGAMI/2021 OLC Coos County/image positions.gpkg") # terra drops z coordinates but they're redundant with the elevation field
  gridMetrics018 = rast(file.path(dataSourcePath, "metrics", "1.8 m", "gridMetrics.vrt")) %>%
    rename_with(~str_c(., if_else(str_ends(., "\\d"), "_018", "018")))
  gridMetrics030 = rast(file.path(dataSourcePath, "metrics", "3.0 m", "gridMetrics.vrt")) %>%
    rename_with(~str_c(., if_else(str_ends(., "\\d"), "_030", "030")))
  gridMetrics046 = rast(file.path(dataSourcePath, "metrics", "4.6 m", "gridMetrics.vrt")) %>%
    rename_with(~str_c(., if_else(str_ends(., "\\d"), "_046", "046")))
  gridMetrics100 = rast(file.path(dataSourcePath, "metrics", "grid metrics 10 m non-normalized v2.tif")) %>%
    rename_with(~str_c(., if_else(str_ends(., "\\d"), "_100", "100")))

  #orthoimageCellSize = max(res(orthoimage))
  dtmCellSize100 = max(res(dtmSlope100))
  gridMetricsCellSize018 = max(res(gridMetrics018))
  gridMetricsCellSize030 = max(res(gridMetrics030))
  gridMetricsCellSize046 = max(res(gridMetrics046))
  gridMetricsCellSize100 = max(res(gridMetrics100))
  
  dsmCrs = crs(dsm, describe = TRUE)
  dsmSlopeAspect005Crs = crs(dsmSlopeAspect005, describe = TRUE)
  dtmSlope100Crs = crs(dtmSlope100, describe = TRUE)
  dtmAspectCos100Crs = crs(dtmAspectCos100, describe = TRUE)
  dtmAspectSin100Crs = crs(dtmAspectSin100, describe = TRUE)
  orthoimageCrs = crs(orthoimage, describe = TRUE)
  gridMetrics018Crs = crs(gridMetrics018, describe = TRUE)
  gridMetrics030Crs = crs(gridMetrics030, describe = TRUE)
  gridMetrics046Crs = crs(gridMetrics046, describe = TRUE)
  gridMetrics100Crs = crs(gridMetrics100, describe = TRUE)
  imageCentersCrs = crs(imageCenters, describe = TRUE)
  # check CRSes for horizontal consistency
  if ((orthoimageCrs$authority != dsmCrs$authority) | (orthoimageCrs$code != dsmCrs$code) |
      (orthoimageCrs$authority != dsmSlopeAspect005Crs$authority) | (orthoimageCrs$code != dsmSlopeAspect005Crs$code) |
      (orthoimageCrs$authority != dtmSlope100Crs$authority) | (orthoimageCrs$code != dtmSlope100Crs$code) |
      (orthoimageCrs$authority != dtmAspectCos100Crs$authority) | (orthoimageCrs$code != dtmAspectCos100Crs$code) |
      (orthoimageCrs$authority != dtmAspectSin100Crs$authority) | (orthoimageCrs$code != dtmAspectSin100Crs$code) |
      (orthoimageCrs$authority != gridMetrics018Crs$authority) | (gridMetrics018Crs$code != dsmCrs$code) |
      (orthoimageCrs$authority != gridMetrics030Crs$authority) | (gridMetrics030Crs$code != dsmCrs$code) |
      (orthoimageCrs$authority != gridMetrics046Crs$authority) | (gridMetrics046Crs$code != dsmCrs$code) |
      (orthoimageCrs$authority != gridMetrics100Crs$authority) | (gridMetrics100Crs$code != dsmCrs$code) |
      (orthoimageCrs$authority != imageCentersCrs$authority) | (imageCentersCrs$code != dsmCrs$code))
  {
    stop(paste0("DSM, orthoimage, and grid metrics CRSes do not match for training polygon ", trainingPolygonIndex))
  }
  # fix up vertical CRSes: terra doesn't know if vertical CRS is used or not so warns on horizontal-compound mismathc
  if ((orthoimageCrs$authority == dsmCrs$authority) & (orthoimageCrs$code == dsmCrs$code))
  {
    crs(dsm) = crs(orthoimage) # give DSM the orthoimage's vertical CRS, TODO: remove once DSM fix propagates
  }
  if ((orthoimageCrs$authority == dsmSlopeAspect005Crs$authority) & (orthoimageCrs$code == dsmSlopeAspect005Crs$code))
  {
    crs(dsmSlopeAspect005) = crs(orthoimage) # give DSM the orthoimage's vertical CRS, TODO: remove once DSM fix propagates
  }
  if ((orthoimageCrs$authority == dtmSlope100Crs$authority) & (orthoimageCrs$code == dtmSlope100Crs$code))
  {
    crs(dtmSlope100) = crs(orthoimage) # give slope the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == dtmAspectCos100Crs$authority) & (orthoimageCrs$code == dtmAspectCos100Crs$code))
  {
    crs(dtmAspectCos100) = crs(orthoimage) # give aspect the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == dtmAspectSin100Crs$authority) & (orthoimageCrs$code == dtmAspectSin100Crs$code))
  {
    crs(dtmAspectSin100) = crs(orthoimage) # give aspect the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == gridMetrics018Crs$authority) & (orthoimageCrs$code == gridMetrics018Crs$code))
  {
    crs(gridMetrics018) = crs(orthoimage) # give 1.8 m grid metrics the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == gridMetrics030Crs$authority) & (orthoimageCrs$code == gridMetrics030Crs$code))
  {
    crs(gridMetrics030) = crs(orthoimage) # give 3.0 m grid metrics the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == gridMetrics046Crs$authority) & (orthoimageCrs$code == gridMetrics046Crs$code))
  {
    crs(gridMetrics046) = crs(orthoimage) # give 4.6 m grid metrics the orthoimage's vertical CRS
  }
  if ((orthoimageCrs$authority == gridMetrics100Crs$authority) & (orthoimageCrs$code == gridMetrics100Crs$code))
  {
    crs(gridMetrics100) = crs(orthoimage) # give 10 m grid metrics the orthoimage's vertical CRS
  }

  imageCentersForKnn = tibble(xy = geom(imageCenters), z = imageCenters$`elevation, ft`) %>% # flatten to coordinates for kNN
    mutate(x = xy[, "x"], y = xy[, "y"]) %>% select(-xy) %>% relocate(x, y, z)
  imageSunPositions = tibble(azimuth = imageCenters$sunAzimuth, elevation = imageCenters$sunElevation) # flatten for in loop lookup

  polygonCreationStart = Sys.time()
  cat(paste0(nrow(trainingPolygons), " training polygons: 0... "))
  trainingData = vector(mode = "list", length = nrow(trainingPolygons))
  for (trainingPolygonIndex in 1:nrow(trainingPolygons))
  {
    # extract training data for this polygon
    # Assumes
    # - training polygons are small so the error of using constant sun and view angles for all cells in the polygon is slight
    # - imaging flight altitude is high enough the mean DSM elevation is an acceptably accurate representation of the polygon
    # - all layers are in the same CRS per checks above and that orthoimagery is of equal or higher resolution than the DSM
    trainingPolygon = trainingPolygons[trainingPolygonIndex]
    if (relate(gridMetrics100, trainingPolygon, relation = "contains") == FALSE)
    {
      next # skip two bare ground training polygons south of the grid metrics area
    }
    polygonOrthoimage = crop(orthoimage, trainingPolygon, touches = FALSE) # don't need to buffer as polygons are drawn on the orthoimage
    polygonDsm = crop(dsm, trainingPolygon, touches = FALSE) # take only cells with centroids within the training polygon, don't need to buffer as DSM and orthoimage have the same resolution and are aligned
    polygonDsmSlopeAspect = crop(dsmSlopeAspect005, trainingPolygon, touches = FALSE)
    
    # buffer to capture adjacent cells for resampling support window: bilinear needs 2x2, cubic and cubic spline are likely 4x4 and lanczos 6x6 (https://github.com/rspatial/terra/issues/1568)
    # resample() interpolates only over its y input, so outputs are matched to the orthoimage. If the buffering around the training 
    # polygon is not wide enough then the resampling method may return NAs.
    # Also, crop from input rasters as windowing may be unreliable (https://github.com/rspatial/terra/issues/1433).
    polygonGridMetrics018 = resample(crop(gridMetrics018, buffer(trainingPolygon, 5*gridMetricsCellSize018), touches = TRUE), polygonOrthoimage, method = classificationOptions$gridMetricsResamplingMethod) # defaults to bilinear for rasters whose first layer is numeric
    polygonGridMetrics030 = resample(crop(gridMetrics030, buffer(trainingPolygon, 5*gridMetricsCellSize030), touches = TRUE), polygonOrthoimage, method = classificationOptions$gridMetricsResamplingMethod)
    polygonGridMetrics046 = resample(crop(gridMetrics046, buffer(trainingPolygon, 5*gridMetricsCellSize046), touches = TRUE), polygonOrthoimage, method = classificationOptions$gridMetricsResamplingMethod)
    polygonGridMetrics100 = resample(crop(gridMetrics100, buffer(trainingPolygon, 5*gridMetricsCellSize100), touches = TRUE), polygonOrthoimage, method = classificationOptions$gridMetricsResamplingMethod)
    polygonsSlope100 = resample(crop(dtmSlope100, buffer(trainingPolygon, dtmCellSize100), touches = TRUE), polygonOrthoimage, method = "bilinear")
    polygonAspectCos100 = resample(crop(dtmAspectCos100, buffer(trainingPolygon, dtmCellSize100), touches = TRUE), polygonOrthoimage, method = "bilinear")
    polygonAspectSin100 = resample(crop(dtmAspectSin100, buffer(trainingPolygon, dtmCellSize100), touches = TRUE), polygonOrthoimage, method = "bilinear")
    
    polygonBrdf = tibble(crds(polygonDsm, df = TRUE, na.rm = FALSE), dsmZ = as.vector(polygonDsm$dsm)) %>%
      mutate(dsmZ = replace_na(dsmZ, min(dsmZ, na.rm = TRUE)))
    imageCenterIndex = knnx.index(imageCentersForKnn, polygonBrdf, k = 1)
    polygonBrdf %<>% mutate(sunAzimuth = imageSunPositions$azimuth[imageCenterIndex], 
                            sunZenithAngle = 90 - imageSunPositions$elevation[imageCenterIndex], 
                            deltaX = imageCentersForKnn$x[imageCenterIndex] - x,
                            deltaY = imageCentersForKnn$y[imageCenterIndex] - y,
                            deltaZ = imageCentersForKnn$z[imageCenterIndex] - dsmZ,
                            viewAzimuth = -180/pi * (atan2(deltaY, deltaX) - pi/2),
                            viewAzimuth = if_else(viewAzimuth > 0, viewAzimuth, 360 + viewAzimuth),
                            viewZenithAngle = 180/pi * atan2(sqrt(deltaX^2 + deltaY^2), deltaZ))
    polygonDsm$sunAzimuth = polygonBrdf$sunAzimuth
    polygonDsm$sunZenithAngle = polygonBrdf$sunZenithAngle
    polygonDsm$viewAzimuth = polygonBrdf$viewAzimuth
    polygonDsm$viewZenithAngle = polygonBrdf$viewZenithAngle

    trainingStatistics = as_tibble(c(polygonOrthoimage, polygonDsm, polygonDsmSlopeAspect, polygonGridMetrics018, polygonGridMetrics030, polygonGridMetrics046, polygonGridMetrics100, polygonsSlope100, polygonAspectCos100, polygonAspectSin100)) %>% # stack rasters
      mutate(polygon = trainingPolygonIndex,
             classification = trainingPolygon$classification)
    #polygonChm = resample(crop(chm, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage)
    trainingData[[trainingPolygonIndex]] = trainingStatistics
    
    if (trainingPolygonIndex %% 25 == 0)
    {
      cat(paste0(trainingPolygonIndex, "... "))
    }
  }
  (Sys.time() - polygonCreationStart)
  
  trainingData = bind_rows(trainingData) %>%
    rename(nir = nearInfrared) %>%
    mutate(classification = as.factor(classification),
           intensitySecondReturn = replace_na(intensitySecondReturn, 0)) %>% # allow cells without second returns to be classified
    relocate(polygon, classification)
  
  # while data is written here for all training cells, it's useful to check what filters are required below
  #print(trainingData %>% filter(is.na(intensityFirstReturn) == FALSE, is.na(dsm) == FALSE, is.na(zMean) == FALSE) %>% summarise(across(everything(), ~sum(is.na(.)))), width = Inf)
  #sum(is.na(trainingData))
  saveRDS(trainingData, paste0("trees/segmentation/classification training data ", nrow(trainingPolygons), " 0.5+1.8+3.0+4.6+10 m.Rds"))
} else {
  # ranger requires complete cases but not all raster cells in all training polygons have all grid metrics or point hits
  # At 1931 polygons (262,761 cells): 4929 cells without LiDAR grid metrics, 609 without RGB+NIR+I1, 233 without scan angle, 174 without DSM
  trainingData = readRDS("trees/segmentation/classification training data 1936 0.5+1.8+3.0+4.6+10 m.Rds") %>%
    filter(is.na(intensityFirstReturn) == FALSE, is.na(dsm) == FALSE, is.na(zMean100) == FALSE)
}

# remove unique polygon identifiers and add derived predictors
trainingData = trainingData %>% 
  filter(is.na(dsmSlope) == FALSE, is.na(cmmSlope3) == FALSE) %>%
  select(-polygon) %>%
  mutate(viewAzimuthSunRelative = sunAzimuth - viewAzimuth, # 0 = forward scatter, ±180 = backscatter, absolute value broken out separately to allow for asymmetric BRDF
         viewAzimuthSunRelative = if_else(viewAzimuthSunRelative > 180, 360 - viewAzimuthSunRelative, if_else(viewAzimuthSunRelative < -180, 360 + viewAzimuthSunRelative, viewAzimuthSunRelative)), # clamp to [0, ±180] to constrain training complexity
         viewAzimuthSunRelativeAbsolute = abs(viewAzimuthSunRelative),
         viewAzimuthSunRelativeCosine = cos(pi/180 * viewAzimuthSunRelative),
         viewAzimuthSunRelativeSine = sin(pi/180 * viewAzimuthSunRelative),
         scanAngleCosine = cos(pi/180 * scanAngleMeanAbsolute),
         scanAngleMeanNormalized = scanAngleMeanAbsolute / scanAngleCosine,
         sunZenithAngleCosine = cos(pi/180 * sunZenithAngle),
         viewZenithAngleCosine = cos(pi/180 * viewZenithAngle),
         dsmAspectSunRelative = sunAzimuth - dsmAspect,
         dsmAspectSunRelative = if_else(dsmAspectSunRelative > 180, 360 - dsmAspectSunRelative, if_else(dsmAspectSunRelative < -180, 360 + dsmAspectSunRelative, dsmAspectSunRelative)),
         dsmAspectSunRelativeAbsolute = abs(dsmAspectSunRelative),
         dsmAspectSunRelativeCosine = cos(pi/180 * dsmAspectSunRelative),
         dsmAspectSunRelativeSine = sin(pi/180 * dsmAspectSunRelative),
         cmmAspectSunRelative = sunAzimuth - cmmAspect3,
         cmmAspectSunRelative = if_else(cmmAspectSunRelative > 180, 360 - cmmAspectSunRelative, if_else(cmmAspectSunRelative < -180, 360 + cmmAspectSunRelative, cmmAspectSunRelative)),
         cmmAspectSunRelativeAbsolute = abs(cmmAspectSunRelative),
         cmmAspectSunRelativeCosine = cos(pi/180 * cmmAspectSunRelative),
         cmmAspectSunRelativeSine = sin(pi/180 * cmmAspectSunRelative),
         aspect100 = 180/pi * atan2(aspect100sine, aspect100cosine), # arctangent does not require inversion or rotation as aspect is being reconstructed from sin(aspect) and cos(aspect)
         aspect100 = if_else(aspect100 > 0, aspect100, 360 + aspect100),
         aspect100sunRelative = sunAzimuth - aspect100, 
         aspect100sunRelative = if_else(aspect100sunRelative > 180, 360 - aspect100sunRelative, if_else(aspect100sunRelative < -180, 360 + aspect100sunRelative, aspect100sunRelative)),
         aspect100sunRelativeAbsolute = abs(aspect100sunRelative),
         aspect100sunRelativeCosine = cos(pi/180 * aspect100sunRelative),
         aspect100sunRelativeSine = sin(pi/180 * aspect100sunRelative),
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
         ndgb = (green - blue) / (green + blue),
         ndgr = (green - red) / (green + red), # same as GRVI
         ndvi = (nir - red) / (nir + red),
         nirBlueRatio = nir / blue, # bNDVI, BAI 98%, GEMI 94%
         nirGreenRatio = nir / green, # nirGreenRatio = pbi, gNDVI 97%
         nirRedRatio = nir / red, # nirRedRatio = pssr, SAVI 94%, eviBackup 96%
         nirNormalized = nir / luminosity, # gNDVI 98%
         normalizedBlue	= blue / (nir + red + green), # NDVI, gNDVI, TVI 95+%
         normalizedGreen = green / (nir + red + green),
         normalizedNir = nir / (nir + red + green),
         osavi = 1.16 * (nir - red) / (nir + red + 0.16),
         redBlueRatio = red / blue, # rNormalized, WBI, NDGR 98+% correlation, coloration, bNormalized 95%
         rdvi = (nir - red) / sqrt(nir + red),
         rgbv = (green^2 - red * as.numeric(blue)) / (green^2 + red * as.numeric(blue)), # 87% with gNormalized, convert explicitly to avoid integer overflow
         rndvi = (nir - red) / sqrt(nir + red), # sometimes mistaken for rdvi?
         # TODO: estimate SAVI L ∈ [0, 1] = [all vegetation, no vegetation] from aerial and ground point counts?
         savi = (1 + 0.5) * (nir - red) / (nir + red + 0.5), # with moderate vegetation cover default L = 0.5, NDVI 100%, eviBackup, oSAVI 99% correlation
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
         aerialMeanNormalized = aerialMean - (dsm - chm), # DSM - CHM = DTM
         zNormalizedEntropy018 = replace_na(zNormalizedEntropy018, 0), # might be better to replace NAs with 1?
         zNormalizedEntropy030 = replace_na(zNormalizedEntropy030, 0),
         zNormalizedEntropy046 = replace_na(zNormalizedEntropy046, 0),
         zNormalizedEntropy100 = replace_na(zNormalizedEntropy100, 0),
         hMax018 = zMax018 - zGroundMean018,
         hMax030 = zMax030 - zGroundMean030,
         hMax046 = zMax046 - zGroundMean046,
         hMax100 = zMax100 - zGroundMean100,
         hMean018 = zMean018 - zGroundMean018,
         hMean030 = zMean030 - zGroundMean030,
         hMean046 = zMean046 - zGroundMean046,
         hMean100 = zMean100 - zGroundMean100,
         hQ10_018 = zQ10_018 - zGroundMean018,
         hQ20_018 = zQ20_018 - zGroundMean018,
         hQ30_018 = zQ30_018 - zGroundMean018,
         hQ40_018 = zQ40_018 - zGroundMean018,
         hQ50_018 = zQ50_018 - zGroundMean018,
         hQ60_018 = zQ60_018 - zGroundMean018,
         hQ70_018 = zQ70_018 - zGroundMean018,
         hQ80_018 = zQ80_018 - zGroundMean018,
         hQ90_018 = zQ90_018 - zGroundMean018,
         hQ10_030 = zQ10_030 - zGroundMean030,
         hQ20_030 = zQ20_030 - zGroundMean030,
         hQ30_030 = zQ30_030 - zGroundMean030,
         hQ40_030 = zQ40_030 - zGroundMean030,
         hQ50_030 = zQ50_030 - zGroundMean030,
         hQ60_030 = zQ60_030 - zGroundMean030,
         hQ70_030 = zQ70_030 - zGroundMean030,
         hQ80_030 = zQ80_030 - zGroundMean030,
         hQ90_030 = zQ90_030 - zGroundMean030,
         hQ10_046 = zQ10_046 - zGroundMean046,
         hQ20_046 = zQ20_046 - zGroundMean046,
         hQ30_046 = zQ30_046 - zGroundMean046,
         hQ40_046 = zQ40_046 - zGroundMean046,
         hQ50_046 = zQ50_046 - zGroundMean046,
         hQ60_046 = zQ60_046 - zGroundMean046,
         hQ70_046 = zQ70_046 - zGroundMean046,
         hQ80_046 = zQ80_046 - zGroundMean046,
         hQ90_046 = zQ90_046 - zGroundMean046,
         hQ05_100 = zQ05_100 - zGroundMean100, # z quantiles mostly 90-99% correlated, lower quantiles more distinct
         hQ10_100 = zQ10_100 - zGroundMean100,
         hQ15_100 = zQ15_100 - zGroundMean100,
         hQ20_100 = zQ20_100 - zGroundMean100,
         hQ25_100 = zQ25_100 - zGroundMean100,
         hQ30_100 = zQ30_100 - zGroundMean100, # unreferenced
         hQ35_100 = zQ35_100 - zGroundMean100, # unreferenced
         hQ40_100 = zQ40_100 - zGroundMean100, # unreferenced
         hQ45_100 = zQ45_100 - zGroundMean100, # unreferenced
         hQ50_100 = zQ50_100 - zGroundMean100, # unreferenced
         hQ55_100 = zQ55_100 - zGroundMean100, # unreferenced
         hQ60_100 = zQ60_100 - zGroundMean100,
         hQ65_100 = zQ65_100 - zGroundMean100, # unreferenced
         hQ70_100 = zQ70_100 - zGroundMean100, # unreferenced
         hQ75_100 = zQ75_100 - zGroundMean100,
         hQ80_100 = zQ80_100 - zGroundMean100,
         hQ85_100 = zQ85_100 - zGroundMean100,
         hQ90_100 = zQ90_100 - zGroundMean100,
         hQ95_100 = zQ95_100 - zGroundMean100)
# drop elevations
trainingData %<>% select(-aerialMean, -dsm, -cmm3, -zMax018, -zMax030, -zMax046, -zMax100, -zMean018, -zMean030, -zMean046, -zMean100, -zGroundMean018, -zGroundMean030, -zGroundMean046, -zGroundMean100, 
                         -zQ10_018, -zQ20_018, -zQ30_018, -zQ40_018, -zQ50_018, -zQ60_018, -zQ70_018, -zQ80_018, -zQ90_018,
                         -zQ10_030, -zQ20_030, -zQ30_030, -zQ40_030, -zQ50_030, -zQ60_030, -zQ70_030, -zQ80_030, -zQ90_030,
                         -zQ10_046, -zQ20_046, -zQ30_046, -zQ40_046, -zQ50_046, -zQ60_046, -zQ70_046, -zQ80_046, -zQ90_046,
                         -zQ05_100, -zQ10_100, -zQ15_100, -zQ20_100, -zQ25_100, -zQ30_100, -zQ35_100, -zQ40_100, -zQ45_100, -zQ50_100, -zQ55_100, -zQ60_100, -zQ65_100, -zQ70_100, -zQ75_100, -zQ80_100, -zQ85_100, -zQ90_100, -zQ95_100)
# drop non-BRDF angles
trainingData %<>% select(-sunAzimuth, -viewAzimuth, -dsmAspect, -cmmAspect3, -aspect100, -aspect100cosine, -aspect100sine)
# drop low importance variables
#trainingData %<>% select(-sipi, -zStdDev, -hMean100, -hQ05_100, -aerialMeanNormalized,
#                         -vari, -ari, -gemi, -hQ15_100, -hQ65_100, -rNormalized, -mexg, -normalizedGreen, -hQ60_100,
#                         -pZaboveZmean, -intensityQ10, -arvi2, -hQ20_100, -pFourthReturn, -rdvi, -rndvi, -evi,
#                         -blue, -intensityQ90, -intensityMean, -intensityStdDev,
#                         -red, -nir, -msavi, -luminosity, -green, -luminosity709, -atsavi, -nirGreenRatio, -chlorophyllGreen, -gndvi,
#                         -intensityQ20, -intensityMeanBelowMedianZ, -pSecondReturn, -pFifthReturn, -pFirstReturn, -intensityQ30, -intensityQ70, -intensityQ80, -intensityQ60, -intensityQ40, -intensityQ50,
#                         -pZaboveThreshold, -zSkew, -zNormalizedEntropy,
#                         -normalizedNir, -nirNormalized, -tvi, -nirRedRatio, -mtvi2, -ndvi, -gdvi2, -eviBackup, -msr, -gdvi, -savi, -ctvi, -osavi,
#                         -nAerial, -nGround, -firstReturns, -secondReturns, -intensitySecondReturn, -scanAngleMeanNormalized)
# drop BRDF trig
#trainingData %<>% select(-sunZenithAngleCosine,
#                         -cmmAspectSunRelativeCosine, -cmmAspectSunRelativeSine,
#                         -dsmAspectSunRelativeCosine, -dsmAspectSunRelativeSine, 
#                         -aspect100sunRelativeCosine, -aspect100sunRelativeSine,
#                         -scanAngleCosine,
#                         -viewAzimuthSunRelativeCosine, -viewAzimuthSunRelativeSine,
#                         -viewZenithAngleCosine)
# drop BRDF signed but leave sun and view elevations
#trainingData %<>% select(-cmmAspectSunRelative, -dsmAspectSunRelative, -aspect100sunRelative, -viewAzimuthSunRelative)
# drop BRDF absolute
#trainingData %<>% select(-cmmAspectSunRelativeAbsolute, -dsmAspectSunRelativeAbsolute, -aspect100sunRelativeAbsolute, -viewAzimuthSunRelativeAbsolute, -scanAngleMeanAbsolute)
# drop sun and view elevations
#trainingData %<>% select(-sunZenithAngle, -viewZenithAngle)
# gNormalized and gli are fairly similar; drop gli
#trainingData %<>% select(-gli)
# absolute view azimuth is preferred; drop signed azimuth
#trainingData %<>% select(-viewAzimuthSunRelative)

if (classificationOptions$predictClasses)
{
  # recode from subclasses to classes
  # levels: 1 = non-tree, 2 = snag, 3 = conifer, 4 = hardwood
  trainingData %<>% mutate(classification = fct_collapse(classification, `non-tree` = c("bare", "bare shadow"),
                                                                         conifer = c("conifer", "conifer shadow", "conifer deep shadow"),
                                                                         hardwood = c("hardwood", "hardwood shadow", "hardwood deep shadow"), 
                                                                         snag = c("brown tree", "grey tree")))
  #unique(trainingData$classification)
  
  # reference predictors
  # orthoimage 7: ndgr, ndvi, gndvi, bndvi, intensity, intensityFirstReturn, intensitySecondReturn
  # non-normalized 8: intensityQ10, intensityPground, pThirdReturn, pZaboveThreshold, hQ05_100, hQ10_100, hQ25_100, hQ75_100
  #predictorVariables = c("classification", "ndgr", "ndvi", "gndvi", "bndvi", "intensityFirstReturn", "intensityFirstGridReturn", "intensitySecondReturn", "intensityQ10", "intensityPground", "pThirdReturn", "pZaboveThreshold", "hQ05_100", "hQ10_100", "hQ25_100", "hQ75_100")
  #rangerTuning = tibble(minNodeSize = 2, mtry = 8, samplingFraction = 0.871)
  
  # PCA+MCA at 1.8, 3.0, 4.6 or 10 m
  predictorVariables = c("classification", "intensityQ40_018", "hMean018", "mari", "pGround018", "intensitySkew018", "msavi", "ndgb", "sunZenithAngleCosine", "intensityMax018", "hQ10_018", "cmmAspectSunRelative", "scanAngleCosine")
  rangerTuning = tibble(mtry = 8, minNodeSize = 2, samplingFraction = 0.894)
  #rangerTuning = tibble(mtry = 8, minNodeSize = 20, samplingFraction = 0.894)
  #predictorVariables = c("classification", "intensityQ40_030", "hMean030", "mari", "pGround030", "intensitySkew030", "msavi", "ndgb", "sunZenithAngleCosine", "intensityMax030", "hQ10_030", "cmmAspectSunRelative", "scanAngleCosine")
  #rangerTuning = tibble(mtry = 8, minNodeSize = 2, samplingFraction = 0.884)
  #predictorVariables = c("classification", "intensityQ40_046", "hMean046", "mari", "pGround046", "intensitySkew046", "msavi", "ndgb", "sunZenithAngleCosine", "intensityMax046", "hQ10_046", "cmmAspectSunRelative", "scanAngleCosine")
  #rangerTuning = tibble(mtry = 7, minNodeSize = 2, samplingFraction = 0.888)
  
  # PCA+MCA at 3.0, 4.6 and 10 m, selection could likely be refined (or intensityQ40_100 -> intensityQ30_046)
  #predictorVariables = c("classification", "ndvi", "mgrv", "pSecondReturn046", "red", "mari", "intensityQ40_100", "intensityQ80_030", "msavi", "nir", "zSkew030", "pThirdReturn030", "sunZenithAngleCosine", "scanAngleMeanAbsolute") # 12: "intensityQ10_046", "intensityStdDev030", "pGround100", "chm", 43: zStdDev100, hQ30_30
  #rangerTuning = tibble(minNodeSize = 2, mtry = 8, samplingFraction = 0.871)
                         
  # 125 -> 16 VSURF, class
  #predictorVariables = c("classification", "hMax100", "chm", "gNormalized", "intensitySkew", "hQ95_100", "intensityFirstGridReturn", "hQ75_100", "hQ85_100", "hQ90_100", "rgbv", "hQ80_100", "intensityQ10", "sunZenithAngle", "viewZenithAngle", "intensityMax", "viewAzimuthSunRelative")
  #rangerTuning = tibble(minNodeSize = 5, mtry = 2, samplingFraction = 0.881)
  #predictorVariables = c("classification", "hMax100", "chm", "gNormalized", "intensitySkew", "hQ95_100", "intensityFirstGridReturn", "hQ75_100", "hQ85_100", "hQ90_100", "rgbv", "hQ80_100", "intensityQ10", "intensityMax")
  #rangerTuning = tibble(minNodeSize = 5, mtry = 2, samplingFraction = 0.871)
  # PCA+MCA 25 subset -> 15 manual, class
  #predictorVariables = c("classification", "ndvi", "gndvi", "intensityFirstReturn", "intensityMean", "intensityQ10", "pGround", "normalizedGreen", "luminosity", "mari", "ari", "hMean100", "chm", "pSecondReturn", "intensityStdDev", "sipi") #, "intensityMeanAboveMedianZ", "msavi", "pThirdReturn", "brvi", "hQ10_100", "chlorophyllVegetation", "zSkew", "sunZenithAngle", "mexg")
  #rangerTuning = tibble(minNodeSize = 2, mtry = 9, samplingFraction = 0.877)
  # PCA+MCA 25 subset -> VSURF 16, class
  #predictorVariables = c("classification", "hMean100", "chm", "mari", "intensityMeanAboveMedianZ", "sipi", "ari", "intensityQ10", "sunZenithAngle", "intensityStdDev", "pThirdReturn", "hQ10_100", "intensityMean", "msavi", "mexg", "pGround", "zSkew")
  #rangerTuning = tibble(minNodeSize = 2, mtry = 6, samplingFraction = 0.864)
} else {
  # PCA+MCA at 1.8, 3.0, 4.6 or 10 m
  #ggplot() + geom_histogram(aes(x = sunZenithAngleCosine), trainingData, binwidth = 0.0001)
  predictorVariables = c("classification", "intensityQ40_018", "hMean018", "mari", "pGround018", "intensitySkew018", "msavi", "ndgb", "sunZenithAngleCosine", "intensityMax018", "hQ10_018", "cmmAspectSunRelative", "scanAngleCosine")
  rangerTuning = tibble(mtry = 8, minNodeSize = 27, samplingFraction = 0.699)
  
  # PCA+MCA at 3.0, 4.6 and 10 m
  #predictorVariables = c("classification", "ndvi", "mgrv", "pSecondReturn046", "red", "mari", "intensityQ40_100", "intensityQ80_030", "msavi", "nir", "zSkew030", "pThirdReturn030", "sunZenithAngleCosine", "scanAngleMeanAbsolute") # 12: "intensityQ10_046", "intensityStdDev030", "pGround100", "chm", 43: zStdDev100, hQ30_30
  #rangerTuning = tibble(minNodeSize = 23, mtry = 12, samplingFraction = 0.656)
  
  #predictorVariables = c("classification", "gNormalized", "chm", "viewZenithAngle", "viewAzimuthSunRelativeAbsolute", "sunZenithAngle")
  #rangerTuning = tibble(minNodeSize = 3, mtry = 2, samplingFraction = 0.883)
  #predictorVariables = c("classification", "gNormalized", "chm", "viewZenithAngle", "viewAzimuthSunRelativeAbsolute")
  #rangerTuning = tibble(minNodeSize = 3, mtry = 2, samplingFraction = 0.892)
}

if (classificationOptions$includeExploratory)
{
  # 9900X: tuneRanger() @ 12 cores and ~69 GB/s DDR bandwidth
  # rows   predictors                 VSURF    cores selected  tune   trees  threads   mtry  min node size  sample fraction
  # 263k   PCA+MCA 1.8 class                   12    12        2.47h  500    12        8     2              0.894
  # 263k   PCA+MCA 3.0 class                   12    12        1.97h  500    12        8     2              0.884
  # 263k   PCA+MCA 4.6 class                   12    12        1.85h  500    12        7     2              0.888
  # 263k   PCA+MCA 1.8 subclass                12    12        3.45h  500    12        8     27             0.699
  # 263k   PCA+MCA 3.0 subclass                12    12        3.59h  500    12        8     27             0.688
  # 263k   PCA+MCA 4.6 subclass                12    12        3.23h  500    12        8     27             0.689
  # 263k   7+8 reference class                 12    15        2.1h   500    12        8     2              0.871
  # 263k   172->125 class                      12    13        1.6h   500    12        5     2              0.871
  # 263k   25 PCA+MCA class manual             12    15        2.6h   500    12        9     2              0.877
  # 263k   25 PCA+MCA class VSURF     (4.4h)   12    16        2.4h   500    12        6     2              0.864
  #
  # 5950X
  # rows   predictors         VSURF  cores selected  tune  trees  threads   mtry  min node size  sample fraction
  # 308k   157                2.6d   15    15        2.3h  500    15        2     2
  # 263k   172                1.3d   15    13        2.2h  500    15        2     2              0.881
  # 263k   172 no trig        18h    15    4         1.8h  500    15        3     2              0.897
  # 263k   172 no trig - 2    18h    15    4         1.8h  500    15        3     2              0.897
  # 263k   172 trig           17h    15    4         1.8h  500    15        3     2              0.890
  # 263k   172 handpick 6                            5.3h  500    15        5     2              0.894
  # 263k   172->125 class     1.0d   15    17        3.9h  500    15        5     2              0.881
  # 263k   172->125 subclass  1.9d   16    4         1.8h  500    16        3     2              0.883
  #
  # variables from PCA+MCA
  vsurfVariables = c("classification", "ndvi", "gndvi", "intensityFirstReturn", "intensityMean", "intensityQ10", "pGround", "normalizedGreen", "luminosity", "mari", "ari", "hMean100", "chm", "pSecondReturn", "intensityStdDev", "sipi", "intensityMeanAboveMedianZ", "msavi", "pThirdReturn", "brvi", "hQ10_100", "chlorophyllVegetation", "zSkew", "sunZenithAngle", "mexg")

  library(VSURF)
  vsurfStartTime = Sys.time()
  classificationVsurf = VSURF(classification ~ ., trainingData %>% select(all_of(vsurfVariables)), ncores = classificationOptions$rangerThreads, parallel = TRUE, RFimplem = "ranger")
  saveRDS(classificationVsurf, "trees/segmentation/classification vsurf 172.25 classes.Rds")
  #classificationVsurf = readRDS("trees/segmentation/classification vsurf 172.25 classes.Rds")
  classificationVsurf$nums.varselect # threshold -> interpretation -> prediction
  classificationVsurf$mean.perf
  classificationVsurf$overall.time
  classificationVsurf$comput.times

  plot(classificationVsurf)
  (variablesThreshold = attributes(classificationVsurf$terms[classificationVsurf$varselect.thres])$term.labels)
  (variablesInterpretation = attributes(classificationVsurf$terms[classificationVsurf$varselect.interp])$term.labels)
  (variablesPrediction = attributes(classificationVsurf$terms[classificationVsurf$varselect.pred])$term.labels)
  
  predictorImportance = tibble(predictor = as.character(attr(classificationVsurf$terms, "predvars"))[classificationVsurf$imp.mean.dec.ind + 2], importance = classificationVsurf$imp.mean.dec) %>% # offset as.character() by two since first element is "list" and second is classification
    mutate(importance = importance / sum(importance), selection = factor(if_else(predictor %in% variablesPrediction, "prediction", if_else(predictor %in% variablesInterpretation, "interpretation", if_else(predictor %in% variablesThreshold, "thresholding", "excluded"))), levels = c("prediction", "interpretation", "thresholding", "excluded")))
  ggplot() +
    geom_col(aes(x = importance, y = fct_reorder(predictor, importance), fill = selection), predictorImportance) +
    labs(x = "normalized variable importance", y = NULL, fill = "VSURF") +
    scale_fill_manual(values = c("forestgreen", "blue2", "darkviolet", "black"))
  ggsave("trees/segmentation/classification vsurf 172.25 class importance.png", width = 14, height = 0.33 * nrow(predictorImportance), units = "cm", dpi = 150)

  
  # unclear if accurate over repeated runs: simplest mitigation is to restart R before every tuneRanger() call
  library(tuneRanger)
  library(mlr)
  #predictorVariables = c("classification", variablesPrediction) # makeClassifTask() breaks if not written separately
  # PCA+MCA 10 m predictorVariables = c("classification", "ndvi", "gndvi", "intensityFirstReturn", "intensityMean", "intensityQ10", "pGround", "normalizedGreen", "luminosity", "mari", "ari", "hMean100", "chm", "pSecondReturn", "intensityStdDev", "sipi", "intensityMeanAboveMedianZ", "msavi", "pThirdReturn", "brvi", "hQ10_100", "chlorophyllVegetation", "zSkew", "sunZenithAngle", "mexg")
  rangerTuneTask = makeClassifTask(data = as.data.frame(trainingData %>% select(all_of(predictorVariables))), target = "classification")
  #estimateStart = Sys.time()
  #estimateTimeTuneRanger(rangerTuneTask, num.trees = 500, num.threads = 15, iters = 70)
  #Sys.time() - estimateStart
  tuneStart = Sys.time()
  rangerTuning = tuneRanger(rangerTuneTask, measure = list(multiclass.brier), num.trees = 500, num.threads = 12, iters = 70) # passing classificationOptions$rangerThreads or a simple variable with the same value errors out
  Sys.time() - tuneStart
  (rangerTuning)
  
  # divide by zero -> infinity: arvi2, evi, vari, zSkew018
  #tibble(name = names(trainingData), infs = colSums(trainingData == Inf, na.rm = TRUE)) %>% filter(infs > 0)
  predictorPca = prcomp(~ ., trainingData %>% select(-classification, -arvi2, -evi, -vari, -zSkew018, -ends_with("030"), -ends_with("046"), -ends_with("100")), scale = TRUE) # ~10 s @ 100k rows, ~20 @ 263k
  factoextra::fviz_eig(predictorPca, ncp = 20)
  factoextra::fviz_pca_var(predictorPca, col.var = "cos2", axes = c(1, 2), labelsize = 2, repel = TRUE)
  ggsave("trees/segmentation/classification PCA 1936 axes 1.8 0102.png", width = 30, height = 30, units = "cm", bg = "white", dpi = 250)
  
  predictorFAMD = FactoMineR::FAMD(trainingData %>% select(-arvi2, -evi, -vari, -zSkew018, -ends_with("030"), -ends_with("046"), -ends_with("100")) %>% sample_n(100000), graph = FALSE)
  factoextra::fviz_famd_var(predictorFAMD, col.ind = "cos2", gradient.cols = c("blue", "orange", "red"), axes = c(5, 4), labelsize = 2, repel = TRUE)
  #corrplot::corrplot(factoextra::get_pca_var(predictorPca)$cos2, is.corr = FALSE)
  #ggcorrplot::ggcorrplot(factoextra::get_pca(predictorPca, element = "var")$cor)
  #factoextra::fviz_contrib(predictorPca, choice = "var", axes = 1:10) + theme(axis.text.x = element_text(size = 7))
  #ggsave("trees/segmentation/classification PCA 1936 contrib 1.8+3.0 0110.png", width = 50, height = 10, units = "cm", bg = "white", dpi = 250)
  
  # imagery metric importance: { rNormalized, greenness, ndgr, mgrv, redBlueRatio, coloration, wbi, red, bNoramlized } 
  #                            { ndgb, luminosity, luminosity709, green, gNormalized, blue, mtvi2, nirRedRatio, ndvi, tvi, osavi, eviBackup, msavi, nir, rgbv, nirNormalized }
  #                            { gemi, sipi, mexg, nirGreenRatio, chlorophyllGreen, gndvi, normalizedBlue, vari, nirBlueRatio, bndvi, normalizedGreen, gli, secondReturns }
  #                            { evi, chlorophyllVegetation, intensitySecondReturn, intensity }
  #                    unused: { first returns }
  # LiDAR metric importance: { intensityPground, pZaboveThreshold, zQ40, 35, 45, hMean100, zQ55, 60, 70, 30, 75,  65, 80, pSecondReturn, zQ80, 25, 75, 85, 90 }
  #                          { pFirstReturn, zQ95, intensityMean, intensityQ30, hMax100, intensityQ20, intensityMeanBelowMeidanZ, intensityQ40, zQ10, intensityQ70, pThirdReturn, intensityQ50, 10, 80, zStdDev, 60, intensityFirstReturn }
  #                          { zQ05, intensityQ90, intensityStdDev, pFourthReturn, zSkew, pZaboveZmean, pFifthReturn, intensitySkew }
  #                  unused: { intensityMax, zNormalizedEntropy, zPc20, zGroundMean100, zQ05, pCuZQ10, zQ10, 15, 20, 25, 30, 35, 40, 45, zPc10, zQ50, zMean, zZQ55, zQ60, zMax, zQ65, 70, 90, 85, 80, 75, 95, zKurtosis }
  filterVarImp(x = trainingData %>% select(-classification, -polygon), y = trainingData$classification) %>%
    mutate(max = pmax(bare, conifer, hardwood, shadow), mean = 0.25 * (bare + conifer + hardwood + shadow)) %>%
    #mutate(mean = 0.333 * (bare + conifer + hardwood)) %>%
    arrange(desc(max), desc(mean)) %>%
    mutate(var = row_number()) %>%
    relocate(var)

  # predictor correlations and components
  #cor(as.matrix(trainingData %>% select(ndgr, ndvi, gndvi, bndvi, intensity, intensitySecondReturn)))
  #cor(as.matrix(trainingData %>% select(hQ05_100, hQ10_100, hQ25_100, hQ75_100, hQ90_100)))
  #cor(as.matrix(trainingData %>% select(intensityPground, intensityQ10, intensity, intensityFirstReturn, intensitySecondReturn)))
  #cor(as.matrix(trainingData %>% select(pFirstReturn, pSecondReturn, pThirdReturn, pFourthReturn, pFifthReturn, pZaboveZmean, pZaboveThreshold)))
  #cor(as.matrix(trainingData %>% select(pGround, pSecondReturn, pThirdReturn, pZaboveThreshold, hQ05_100, hQ10_100, hQ25_100, hQ75_100, intensityMean, intensityStdDev, zMean, zSkew, pZaboveZmean)))
  #predictorCorrelations = cor(as.matrix(trainingData %>% select(all_of(predictorVariables), bNormalized, intensitySkew, slope100, hQ10_100, aspect100sunRelativeAbsolute) %>% select(-classification)))
  #predictorCorrelations = cor(as.matrix(trainingData %>% sample_n(20000) %>% mutate(classification = as.integer(classification))))
  ##ggcorrplot::ggcorrplot(predictorCorrelations)
  #predictorCorrelations[c("ndvi", "triangularVegIndex"), c("ndvi", "triangularVegIndex")]
  #predictorCorrelations[c("ndvi", "rndvi", "bndvi", "gndvi", "triangularVegIndex"), c("ndvi", "rndvi", "bndvi", "gndvi", "triangularVegIndex")]
  #tibble(variable = rownames(predictorCorrelations), classification = predictorCorrelations[, "classification"]) %>% 
  #  slice_max(abs(classification), n = 20)
  #cor(as.matrix(trainingData %>% select(-polygon, -classification)))[, "ndvi"]
}

# random forest: 1936 training polygons
# 9900X + DDR5-5600 or 6000: 12 core -> ~1.5 GB DDR @ 66 GB/s (subclass, 5 predictors) to 46 GB/s (class, 16 predictors)
# cells   method                             training                       cores   fit time    accuracy    grid metrics
# 263k    ranger 12 1.8 class                2x25 mtry 8 minNode 2 0.894    12      29.8m       0.976       1.8 m cubic
# 263k    ranger 12 1.8 class                2x25 mtry 8 minNode 20 0.894   12      23.6m       0.972       1.8 m cubic
# 263k    ranger 12 1.8 subclass             2x25 mtry 8 minNode 27 0.699   12      1.02h       0.874       1.8 m cubic
# 263k    ranger 12 3.0 class                2x25 mtry 8 minNode 2 0.884    12      21.8m       0.984       3.0 m cubic
# 263k    ranger 12 4.6 class                2x25 mtry 7 minNode 2 0.888    12      19.2m       0.988       4.6 m cubic (19.5m @ DDR5-5600 -> +1.5%)
# 263k    ranger + reference 7+8             2x25 mtry 8 minNode 2 0.871    12      21 m        0.994       10 m cubic from non-normalized clouds
# 263k    ranger + vsurf 172->5 subclass     2x25 mtry 3 minNode 2 0.881    12      13 m        0.996       10 m cubic from non-normalized clouds
# 263k    ranger PCA+MCA 25->15 class        2x25 mtry 9 minNode 2 0.877    12      18 m        0.987       10 m cubic from non-normalized clouds
# 263k    ranger PCA+MCA 25->16 class        2x25 mtry 6 minNode 2 0.864    12      17 m        0.996       10 m cubic from non-normalized clouds
# 263k    ranger + vsurf 172->14 class       2x25 mtry 5 minNode 2 0.871    12      8.3 m       0.993       10 m cubic from non-normalized clouds
# 263k    ranger + vsurf 172->16 class       2x25 mtry 5 minNode 2 0.883    12      7 m         0.999       10 m cubic from non-normalized clouds
#
# 5950X + DDR4-3200
# cells   method                             cross validation               cores   fit time   accuracy    κ       grid metrics
# 273k    ranger + vsurf 15                  2x10               3x3         16      1.3h       0.998       0.998   10 m cubic from non-normalized clouds
# 273k    ranger + vsurf 4                   2x25               1x1         16      35m        0.993       0.991   10 m cubic from non-normalized clouds
# 273k    ranger + vsurf 4 + sunE            2x25               1x1         16      28m        0.995       0.995   10 m cubic from non-normalized clouds
# 273k    ranger + vsurf 4 + sunE + preval   2x25               1x1         16      1.0h       0.996       0.995   10 m cubic from non-normalized clouds
# 273k    ranger + vsurf 4 + sunE + bNorm    2x25               1x1         16      55m        0.978       0.974   10 m cubic from non-normalized clouds
# 273k    ranger + vsurf 4 + handpick 6      2x25               1x1         16      1.3h       0.982       0.979   10 m cubic from non-normalized clouds
# PCA+MCA selected predictors
#predictorVariables = c("classification", "ndvi", "gndvi", "intensityFirstReturn", "intensityMean", "intensityQ10", "pGround", "normalizedGreen", "luminosity", "mari", "ari", "hMean100", "chm", "pSecondReturn", "intensityStdDev", "sipi", "intensityMeanAboveMedianZ", "msavi", "pThirdReturn", "brvi", "hQ10_100", "chlorophyllVegetation", "zSkew", "sunZenithAngle", "mexg")
# VSURF selected predictors without sun or view angles
# predictorVariables = c("classification", "hMax100", "hQ95_100", "gli", "intensityMeanAboveMedianZ", "hQ75_100", "intensitySkew", "intensityFirstGridReturn")
# VSURF selected predictors with sun and view angles
#predictorVariables = c("classification", "hMax100", "hQ95_100", "gli", "intensityMeanAboveMedianZ", "hQ75_100", "intensitySkew", "intensityFirstGridReturn", "hQ80_100", "viewAzimuthSunRelativeCosine", "viewAzimuthSunRelativeAbsolute", "viewZenithAngleCosine", "viewZenithAngle", "viewAzimuthSunRelative")
# VSURF selected predictors with sun, view, and prevailing slope angles + interpolation improvements
#predictorVariables = c("classification", "chm", "gNormalized", "gli", "hMax100", "viewZenithAngleCosine", "viewZenithAngle", "hQ95_100", "viewAzimuthSunRelativeCosine", "viewAzimuthSunRelativeAbsolute", "hQ85_100", "viewAzimuth", "intensityMeanAboveMedianZ", "viewAzimuthSunRelativeSine", "viewAzimuthSunRelative", "sunAzimuth")
#predictorVariables = c("classification", "chm", "gNormalized", "gli", "hMax100", "viewZenithAngleCosine", "viewZenithAngle", "hQ95_100", "viewAzimuthSunRelativeCosine", "viewAzimuthSunRelativeAbsolute", "hQ85_100", "viewAzimuth", "intensityFirstGridReturn", "viewAzimuthSunRelativeSine", "viewAzimuthSunRelative")
# VSURF selected predictors without BRDF trig
#predictorVariables = c("classification", "gNormalized", "chm", "viewZenithAngle", "viewAzimuthSunRelativeAbsolute")
# VSURF selected predictors with BRDF trig
#predictorVariables = c("classification", "chm", "gNormalized", "viewAzimuthSunRelativeCosine", "viewZenithAngleCosine")
# four predictor BRDF extensions
#predictorVariables = c("classification", "gNormalized", "chm", "viewZenithAngle", "viewAzimuthSunRelativeAbsolute", "sunZenithAngle", "slope100", "aspect100sunRelative")
#predictorCorrelations = cor(as.matrix(trainingData %>% select(all_of(predictorVariables[2:length(predictorVariables)]))))
#ggcorrplot::ggcorrplot(predictorCorrelations)

repeatedCrossValidation = trainControl(method = "repeatedcv", number = 2, repeats = 25, verboseIter = TRUE)

fitStart = Sys.time()
randomForestFit = train(classification ~ ., data = trainingData %>% select(all_of(predictorVariables)), method = "ranger", trControl = repeatedCrossValidation, # importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = rangerTuning$mtry,
                                               splitrule = 'gini',
                                               min.node.size = rangerTuning$minNodeSize),
                        sample.fraction = rangerTuning$samplingFraction,
                        num.threads = classificationOptions$rangerThreads)
(randomForestFitTime = Sys.time() - fitStart)
saveRDS(randomForestFit, file = file.path(getwd(), "trees/segmentation/classificationRandomForest PCA12 1.8m m8n27 cubic subclass.Rds"))
randomForestFit


# sanity check random forest recall
randomForestFit = readRDS(file = file.path(getwd(), "trees/segmentation/classificationRandomForest vsurf 172.13 class 10 m cubic.Rds"))
randomForestPrediction = tibble(actual = trainingData$classification, predicted = predict(randomForestFit$finalModel, trainingData, num.threads = classificationOptions$rangerThreads)$predictions) %>%
  mutate(actualClass = fct_collapse(actual, bare = c("bare", "bare shadow"),
                                            conifer = c("conifer", "conifer shadow", "conifer deep shadow"),
                                            hardwood = c("hardwood", "hardwood shadow", "hardwood deep shadow"),
                                            snag = c("brown tree", "grey tree")),
         predictedClass = fct_collapse(predicted, bare = c("bare", "bare shadow"),
                                                  conifer = c("conifer", "conifer shadow", "conifer deep shadow"),
                                                  hardwood = c("hardwood", "hardwood shadow", "hardwood deep shadow"),
                                                  snag = c("brown tree", "grey tree")))
randomForestConfusionMatrix = confusionMatrix(randomForestPrediction$predicted, randomForestPrediction$actual)
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

if (classificationOptions$includeExploratory)
{
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
}



## tile checking
if (classificationOptions$includeExploratory)
{
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
}


## virtual raster generation and manual editing
# Since gdalbuildvrt doesn't flow metadata (https://github.com/OSGeo/gdal/issues/3627) VRTs need manual editing
# after creation to add <VRTRasterBand/Description> elements naming raster bands. This is needed to set band
# name to Z and orthoimage names to red, green, blue, nir, intensity, intensitySecondReturn, firstReturns, and 
# secondReturns.
# create .vrt for DSM after Get-Dsm
if (classificationOptions$includeSetup)
{
  dsmSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/DSM v3 alpha"
  dsmFilePaths = file.path(dsmSourcePath, list.files(dsmSourcePath, "\\.tif$"))
  vrt(dsmFilePaths, file.path(dsmSourcePath, "DSM dsm cmm3 chm terra.vrt"), options = c("b", 1, "b", 2, "b", 3), overwrite = TRUE, set_names = TRUE)
  
  # create .vrt for orthoimages after orthoImageJob.R has processed all tiles
  orthoimageSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage v2 16"
  orthoimageFilePaths = file.path(orthoimageSourcePath, list.files(orthoimageSourcePath, "\\.tif$"))
  vrt(orthoimageFilePaths, file.path(orthoimageSourcePath, "orthoimage.vrt"), overwrite = TRUE, set_names = TRUE)
  
  # create .vrt for hardwood-conifer classification after classificationJob.R has processed all tiles
  classificationSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/classification ortho6 + 10 m non-normalized 9 cubic"
  classificationFilePaths = file.path(classificationSourcePath, list.files(classificationSourcePath, "\\.tif$"))
  vrt(classificationFilePaths, file.path(classificationSourcePath, "classification.vrt"), overwrite = TRUE, set_names = TRUE)
}