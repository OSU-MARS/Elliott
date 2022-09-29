## weather grid creation
# grid setup
# Vector -> Research -> Create Grid + Select By Location [on buffered Elliott boundary] + Invert Feature Selection
#  extent 4 km: 102000,138000,190000,226000 [EPSG:6556] (left, right, bottom, top)
#         1 km: 102000,136000,190000,226000
#        800 m: 102000,135600,190000,226000
#        400 m: 102000,135600,190400,226000
#        200 m: 102000,135400,190600,225800
#        100 m: 102000,135300,190700,225800
#           1°: -125,-116,41,50 (Oregon-California border 42°, Washington-Canada 49°)
# soils setup
# NASA ORNL: global 1 km soil thickness grids (https://doi.org/10.3334/ORNLDAAC/1304)
#   download from https://daac.ornl.gov/SOILS/guides/Global_Soil_Regolith_Sediment.html
# POLARIS: tiles download with POLARIS.R -> ESRF virtual raster builds -> crop ESRF van Genuchten parameter layers
# gSSURGO (preferred to SSURGO due to more convenient distribution format), https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_053628#value
#  add MUPOLYGON + component + muaggatt [+ Valu1 + ... ] from .gdb and join on mukey (mapunit key)
# 
# SSURGO: Web Soil Survey -> Download Soils Data -> Soil Survey Area (SSURGO) -> State = Oregon + County = { Coos, Douglas }
#   spatial/soilmu_a_*.shp + tabular/{ cfprod.txt -> cfprod.csv, comp.txt -> component.csv, muaggatt.txt -> muaggatt.csv }
#     .txt to .csv conversion in Excel: pipe separated data + paste in headers (search for SSURGO table name for metadata, e.g. https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/?cid=nrcs142p2_053627
#   QGIS: add .csvs as delimited text (no geometry) and join to .shp on mukey
#   transfer joined columns in field calculator (below)
#   Vector -> Data Management -> Merge Vector Layers + Coos and Douglas shapefiles + crop to Elliott


## weather grid field calculator
# following manual deletion of resource units (100 m grid) in Umpqua River
# latitude
y(transform(centroid($geometry), @layer_crs , 'EPSG:4326'))
# longitude
x(transform(centroid($geometry), @layer_crs , 'EPSG:4326'))
# grid cell identifiers
name = concat('N', lpad(format_number(0.001 * 0.5 * ("top" + "bottom"), 1, '0'), 4, '0'), ' E', lpad(format_number(0.001 * 0.5 * ("right" + "left"), 1, '0'), 4, '0'))
# grid cell indices, increasing in nothing and easting like EPSG:6556
grid_y = 8 + 1 - if((id % 8) = 0, 8, id % 8) # 4 km grid: 7 cells east-west by 8 cells north-south
grid_x = floor((id - 1) / 8) + 1 # 4 km

## POLARIS raster calculator: van Genuchten soil water retention curve parameters
# θr and θs are linear and can be therefore calculated as volume weighted arithmetic means
#1/200 * (5 * "thetaR mean ESRF 0-5 cm@1" + 10 * "thetaR mean ESRF 5-15 cm@1" + 15 * "thetaR mean ESRF 15-30 cm@1"
#+ 30 * "thetaR mean ESRF 30-60 cm@1" + 40 * "thetaR mean ESRF 60-100 cm@1" + 100 * "thetaR mean ESRF 100-200 cm@1")
#1/200 * (5 * "thetaS mean ESRF 0-5 cm@1" + 10 * "thetaS mean ESRF 5-15 cm@1" + 15 * "thetaS mean ESRF 15-30 cm@1"
#+ 30 * "thetaS mean ESRF 30-60 cm@1" + 40 * "thetaS mean ESRF 60-100 cm@1" + 100 * "thetaS mean ESRF 100-200 cm@1")

# .vrt clipping to ESRF 4 km
import processing
from qgis.core import *

elliott4km = QgsProject.instance().mapLayersByName('ESRF 4 km bounding box EPSG:6556')[0]
polarisRaster = QgsProject.instance().mapLayersByName('theta_s mean ESRF 100-200 cm')[0]

clipRasterParameters = { 'INPUT': polarisRaster, # processing.algorithmHelp('gdal:cliprasterbymasklayer')
                         'MASK': elliott4km, 
                         'TARGET_CRS': 'EPSG:6556',
                         'CROP_TO_CUTLINE': True,
                         'DATA_TYPE': 6, # loat32
                         'EXTRA': '-co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9',
                         'OUTPUT': 'C:/Users/westjoh/PhD/Elliott/GIS/POLARIS/thetaS/mean/thetaS mean ESRF 100-200 cm.tif' }
processing.runAndLoadResults("gdal:cliprasterbymasklayer", clipRasterParameters)

## weather grid Python
#from math import floor
from qgis.analysis import QgsZonalStatistics
import os

# elevation raster and progress update for zonal statistics on blocked :-( UI thread
#elliotPath = os.path.join(os.getenv('USERPROFILE'), 'PhD', 'Elliott', 'GIS')
#def progress_changed(progress):
#    iface.statusBarIface().showMessage("Zonal statistics: {:.2f}%...".format(progress))

# vector shapefile of weather cells
#weatherCells = QgsVectorLayer(os.path.join(elliotPath, 'weather', 'grid 4 km test.shp'), 'zonepolygons', "ogr") # open direct from disk rather than through project: creates separate file handle
#weatherCells = QgsProject.instance().mapLayersByName('iLand 4 km')[0]
#weatherCells = QgsProject.instance().mapLayersByName('iLand 800 m')[0]
#weatherCells = QgsProject.instance().mapLayersByName('iLand 400 m')[0]
#weatherCells = QgsProject.instance().mapLayersByName('iLand 200 m')[0]
weatherCells = QgsProject.instance().mapLayersByName('iLand 100 m')[0]
weatherCellProvider = weatherCells.dataProvider()

# names, latitudes, and longitudes of cells
if weatherCellProvider.fieldNameIndex("name") == -1:
    nameResult = weatherCellProvider.addAttributes([QgsField("name", QVariant.String)])
    nameResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("latitude") == -1:
    latitudeResult = weatherCellProvider.addAttributes([QgsField("latitude", QVariant.Double)])
    latitudeResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("longitude") == -1:
    longitudeResult = weatherCellProvider.addAttributes([QgsField("longitude", QVariant.Double)])
    longitudeResult = weatherCells.updateFields()

nameColumnIndex = weatherCellProvider.fieldNameIndex("name")
latitudeColumnIndex = weatherCellProvider.fieldNameIndex("latitude")
longitudeColumnIndex = weatherCellProvider.fieldNameIndex("longitude")
toWgs84 = QgsCoordinateTransform(weatherCells.crs(),
                                 QgsCoordinateReferenceSystem("EPSG:4326"), QgsProject.instance())

cellChanges = {}
for weatherCell in weatherCells.getFeatures():
    cellCentroid = weatherCell.geometry().centroid()
    cellCentroidAsPoint = cellCentroid.asPoint()
    cellName = 'N' + "{:0>5.1f}".format(0.001 * cellCentroidAsPoint.y()) + ' E' + "{:0>5.1f}".format(0.001 * cellCentroidAsPoint.x())
    #
    transformResult = cellCentroid.transform(toWgs84)
    cellCentroidAsPoint = cellCentroid.asPoint()
    cellChanges[weatherCell.id()] = { nameColumnIndex: cellName, latitudeColumnIndex: cellCentroidAsPoint.y(), longitudeColumnIndex: cellCentroidAsPoint.x() }

nameLatitudeLongitudeResult = weatherCellProvider.changeAttributeValues(cellChanges)

# weather cell names for resource units
if weatherCellProvider.fieldNameIndex("name4km") == -1:
    nameResult = weatherCellProvider.addAttributes([QgsField("name4km", QVariant.String)])
    nameResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("name800m") == -1:
    nameResult = weatherCellProvider.addAttributes([QgsField("name800m", QVariant.String)])
    nameResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("name400m") == -1:
    nameResult = weatherCellProvider.addAttributes([QgsField("name400m", QVariant.String)])
    nameResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("name200m") == -1:
    nameResult = weatherCellProvider.addAttributes([QgsField("name200m", QVariant.String)])
    nameResult = weatherCells.updateFields()

name4kmColumnIndex = weatherCellProvider.fieldNameIndex("name4km")
name800mColumnIndex = weatherCellProvider.fieldNameIndex("name800m")
name400mColumnIndex = weatherCellProvider.fieldNameIndex("name400m")
name200mColumnIndex = weatherCellProvider.fieldNameIndex("name200m")

cellChanges = {}
for weatherCell in weatherCells.getFeatures():
    cellCentroidAsPoint = weatherCell.geometry().centroid().asPoint()
    cellGridPositionX = cellCentroidAsPoint.x() - 102000
    cellGridPositionY = cellCentroidAsPoint.y() - 190000
    #
    weather4kmCellX = 4000 * math.floor(cellGridPositionX / 4000) + 102000 + 2000
    weather4kmCellY = 4000 * math.floor(cellGridPositionY / 4000) + 190000 + 2000
    weatherName4km = 'N' + "{:0>5.1f}".format(0.001 * weather4kmCellY) + ' E' + "{:0>5.1f}".format(0.001 * weather4kmCellX)
    #
    weather800mCellX = 800 * math.floor(cellGridPositionX / 800) + 102000 + 400
    weather800mCellY = 800 * math.floor(cellGridPositionY / 800) + 190000 + 400
    weatherName800m = 'N' + "{:0>5.1f}".format(0.001 * weather800mCellY) + ' E' + "{:0>5.1f}".format(0.001 * weather800mCellX)
    #
    weather400mCellX = 400 * math.floor(cellGridPositionX / 400) + 102000 + 200
    weather400mCellY = 400 * math.floor(cellGridPositionY / 400) + 190000 + 200
    weatherName400m = 'N' + "{:0>5.1f}".format(0.001 * weather400mCellY) + ' E' + "{:0>5.1f}".format(0.001 * weather400mCellX)
    #
    weather200mCellX = 200 * math.floor(cellGridPositionX / 200) + 102000 + 100
    weather200mCellY = 200 * math.floor(cellGridPositionY / 200) + 190000 + 100
    weatherName200m = 'N' + "{:0>5.1f}".format(0.001 * weather200mCellY) + ' E' + "{:0>5.1f}".format(0.001 * weather200mCellX)
    #
    cellChanges[weatherCell.id()] = { name4kmColumnIndex: weatherName4km, 
                                      name800mColumnIndex: weatherName800m,
                                      name400mColumnIndex: weatherName400m,
                                      name200mColumnIndex: weatherName200m }

weatherNameResult = weatherCellProvider.changeAttributeValues(cellChanges)

# get zonal statistics over cells
# https://api.qgis.org/api/3.22/classQgsZonalStatistics.html
# QgsZonalStatistics (QgsVectorLayer *polygonLayer, const QString &rasterFile, const QString &attributePrefix="", int rasterBand=1)
# Count = 1, Sum = 2, Mean = 4, Median = 8, StDev = 16, Min = 32, Max = 64, Range = 128, Minority = 256, Majority = 512, Variety = 1024, Variance = 2048
southCoastDem = QgsProject.instance().mapLayersByName('2009+2015 OLC Bare_Earth')[0]
elevationStatistics = QgsZonalStatistics(weatherCells, southCoastDem, 'elevation_', 1, QgsZonalStatistics.Min | QgsZonalStatistics.Mean | QgsZonalStatistics.Max)
elevationResult = elevationStatistics.calculateStatistics(None) # takes 3-4 minutes (Zen 3, 4.7 GHz) with 3 ft DEM over 4 km ESRF grid

# convert elevations to meters and "rename" fields
# QgsVectorDataProvider.renameAttributes() works with shapefiles and geopackages but nulls renamed columns with flat
# geobuffers. As a workaround, copy and delete the columns.
if weatherCellProvider.fieldNameIndex("elevationInMmin") == -1:
    elevationMinResult = weatherCellProvider.addAttributes([QgsField("elevationInMmin", QVariant.Double)])
    elevationMinResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("elevationInMmean") == -1:
    elevationMeanResult = weatherCellProvider.addAttributes([QgsField("elevationInMmean", QVariant.Double)])
    elevationMeanResult = weatherCells.updateFields()

if weatherCellProvider.fieldNameIndex("elevationInMmax") == -1:
    elevationMaxResult = weatherCellProvider.addAttributes([QgsField("elevationInMmax", QVariant.Double)])
    elevationMaxResult = weatherCells.updateFields()

elevationMinFrom = weatherCellProvider.fieldNameIndex("elevation_min")
elevationMeanFrom = weatherCellProvider.fieldNameIndex("elevation_mean")
elevationMaxFrom = weatherCellProvider.fieldNameIndex("elevation_max")

elevationMinTo = weatherCellProvider.fieldNameIndex("elevationInMmin")
elevationMeanTo = weatherCellProvider.fieldNameIndex("elevationInMmean")
elevationMaxTo = weatherCellProvider.fieldNameIndex("elevationInMmax")
cellChanges = {}
for weatherCell in weatherCells.getFeatures():
    minElevation = weatherCell[elevationMinFrom]
    if minElevation != NULL: minElevation *= 0.3048
    meanElevation = weatherCell[elevationMeanFrom]
    if meanElevation != NULL: meanElevation *= 0.3048
    maxElevation = weatherCell[elevationMaxFrom]
    if maxElevation != NULL: maxElevation *= 0.3048
    cellChanges[weatherCell.id()] = { elevationMinTo: minElevation, 
                                     elevationMeanTo: meanElevation, 
                                      elevationMaxTo: maxElevation }

elevationResult = weatherCellProvider.changeAttributeValues(cellChanges)

weatherCellProvider.deleteAttributes([ elevationMinFrom, elevationMeanFrom, elevationMaxFrom ])
weatherCells.updateFields()

#renameResult = weatherCellProvider.renameAttributes({ elevationMinIndex: 'elevationMin',
#                                                      elevationMeanIndex: 'elevationMean',
#                                                      elevationMaxIndex: 'elevationMax' })
#weatherCells.updateFields()

# soil depth by vector zonal statistics - replaced by ESRF soil depth raster
#gSSURGO = QgsProject.instance().mapLayersByName('ESRF + Hakki + 2 km soil depth in m 3.00 m EPSG6556')[0]
#soilDepthParameters = { 'INPUT': weatherCells.source(), 
#                        'JOIN': gSSURGO.source(),
#                        'JOIN_FIELDS': 'gSSURGO_OR — Valu1_rootznemc',
#                        'PREDICATE': 0, # intersects
#                        'SUMMARIES': [6],
#                        'OUTPUT': 'TEMPORARY_OUTPUT' }
#soilDepthResult = processing.run("qgis:joinbylocationsummary", soilDepthParameters)
#soilDepthLayer = soilDepthResult['OUTPUT']
#QgsProject.instance().addMapLayer(soilDepthLayer)

if weatherCellProvider.fieldNameIndex("soilDepthInCMmean") == -1:
    soilDepthResult = weatherCellProvider.addAttributes([QgsField("soilDepthInCMmean", QVariant.Double)])
    soilDepthResult = weatherCells.updateFields()

soilDepthIndex = weatherCellProvider.fieldNameIndex("soilDepthInCMmean")

cellChanges = {}
for weatherCell in soilDepthLayer.getFeatures():
    cellChanges[weatherCell.id() - 1] = { soilDepthIndex: weatherCell["gSSURGO_OR — Valu1_rootznemc_mean"] }

soilDepthResult = weatherCellProvider.changeAttributeValues(cellChanges)


# van Genuchten soil properties by POLARIS layer and rock fragments by SoilGrids layer
for layer in ["0-5", "5-15", "15-30", "30-60", "60-100", "100-200"]:
    fieldDepth = layer.replace('-', '_')
    #
    #vanGenuchtenAlpha = QgsProject.instance().mapLayersByName('alpha mean ESRF ' + layer + ' cm')[0]
    #alphaStatistics = QgsZonalStatistics(weatherCells, vanGenuchtenAlpha, 'vanGenuchtenLog10AlphaInCm' + fieldDepth, 1, QgsZonalStatistics.Mean)
    #alphaResult = alphaStatistics.calculateStatistics(None)
    #
    #vanGenuchtenN = QgsProject.instance().mapLayersByName('n mean ESRF ' + layer + ' cm')[0]
    #nStatistics = QgsZonalStatistics(weatherCells, vanGenuchtenN, 'vanGenuchtenN' + fieldDepth, 1, QgsZonalStatistics.Mean)
    #nResult = nStatistics.calculateStatistics(None)
    #
    #vanGenuchtenThetaR = QgsProject.instance().mapLayersByName('thetaR mean ESRF ' + layer + ' cm')[0]
    #thetaRstatistics = QgsZonalStatistics(weatherCells, vanGenuchtenThetaR, 'vanGenuchtenThetaR' + fieldDepth, 1, QgsZonalStatistics.Mean)
    #thetaRresult = thetaRstatistics.calculateStatistics(None)
    #
    #vanGenuchtenThetaS = QgsProject.instance().mapLayersByName('thetaS mean ESRF ' + layer + ' cm')[0]
    #thetaSstatistics = QgsZonalStatistics(weatherCells, vanGenuchtenThetaS, 'vanGenuchtenThetaS' + fieldDepth, 1, QgsZonalStatistics.Mean)
    #thetaSresult = thetaSstatistics.calculateStatistics(None)
    #
    coarseFragments = QgsProject.instance().mapLayersByName('coarse fragments ' + layer + ' cm')[0]
    fragmentStatistics = QgsZonalStatistics(weatherCells, coarseFragments, 'coarseFragments' + fieldDepth, 1, QgsZonalStatistics.Mean)
    fragmentResult = fragmentStatistics.calculateStatistics(None)

# filling of missing van Genuchten properties
# select resource units with NULL and run field calculator to fill
weatherCells = QgsProject.instance().mapLayersByName('iLand 100 m')[0]
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm0_5mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm0_5mean, filter := vanGenuchtenLog10AlphaInCm0_5mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm5_15mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm5_15mean, filter := vanGenuchtenLog10AlphaInCm5_15mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm15_30mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm15_30mean, filter := vanGenuchtenLog10AlphaInCm15_30mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm30_60mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm30_60mean, filter := vanGenuchtenLog10AlphaInCm30_60mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm60_100mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm60_100mean, filter := vanGenuchtenLog10AlphaInCm60_100mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenLog10AlphaInCm100_200mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenLog10AlphaInCm100_200mean, filter := vanGenuchtenLog10AlphaInCm100_200mean is not NULL, limit := 8, max_distance := 275))

weatherCells.selectByExpression("\"vanGenuchtenN0_5mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN0_5mean, filter := vanGenuchtenN0_5mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenN5_15mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN5_15mean, filter := vanGenuchtenN5_15mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenN15_30mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN15_30mean, filter := vanGenuchtenN15_30mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenN30_60mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN30_60mean, filter := vanGenuchtenN30_60mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenN60_100mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN60_100mean, filter := vanGenuchtenN60_100mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenN100_200mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenN100_200mean, filter := vanGenuchtenN100_200mean is not NULL, limit := 8, max_distance := 275))

weatherCells.selectByExpression("\"vanGenuchtenThetaR0_5mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR0_5mean, filter := vanGenuchtenThetaR0_5mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaR5_15mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR5_15mean, filter := vanGenuchtenThetaR5_15mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaR15_30mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR15_30mean, filter := vanGenuchtenThetaR15_30mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaR30_60mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR30_60mean, filter := vanGenuchtenThetaR30_60mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaR60_100mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR60_100mean, filter := vanGenuchtenThetaR60_100mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaR100_200mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaR100_200mean, filter := vanGenuchtenThetaR100_200mean is not NULL, limit := 8, max_distance := 275))

weatherCells.selectByExpression("\"vanGenuchtenThetaS0_5mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS0_5mean, filter := vanGenuchtenThetaS0_5mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaS5_15mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS5_15mean, filter := vanGenuchtenThetaS5_15mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaS15_30mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS15_30mean, filter := vanGenuchtenThetaS15_30mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaS30_60mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS30_60mean, filter := vanGenuchtenThetaS30_60mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaS60_100mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS60_100mean, filter := vanGenuchtenThetaS60_100mean is not NULL, limit := 8, max_distance := 275))
weatherCells.selectByExpression("\"vanGenuchtenThetaS100_200mean\" is NULL")
array_mean(overlay_nearest('grid_100_m_04f8d79a_9b89_4268_9973_9155a8b37489', vanGenuchtenThetaS100_200mean, filter := vanGenuchtenThetaS100_200mean is not NULL, limit := 8, max_distance := 275))

# export all weather grid resolutions .csv for manipulation in R
weatherCellLayerNames = ['iLand 4 km', 'iLand 800 m', 'iLand 400 m', 'iLand 200 m', 'iLand 100 m']
csvOptions = QgsVectorFileWriter.SaveVectorOptions()
csvOptions.driverName = 'CSV'
csvOptions.fileEncoding = 'UTF-8'

for layerName in weatherCellLayerNames:
    weatherCells = QgsProject.instance().mapLayersByName(layerName)[0]
    weatherCellProvider = weatherCells.dataProvider()
    csvFilePath = os.path.splitext(weatherCellProvider.dataSourceUri())[0] + '.csv'
    QgsVectorFileWriter.writeAsVectorFormatV3(weatherCells, csvFilePath, QgsCoordinateTransformContext(), csvOptions)

# release layers if opened from file rather than using project instances
#del weatherCells
#del southCoastDem

## single tree file joins
#   Coos_County_Tile_Index: gridID
#   Elliott Stand Data Feb2022 EPSG6556: standID
#   iLand 100 m: resourceUniT id, resourceUniT bufferDist
# use within as join predicate for bufferDist to avoid exact rounding issues with integer truncation to resource unit indices
# temporarily: create x = $x and y = $y fields as XTOP and YTOP aren't set correctly
# https://docs.qgis.org/3.22/en/docs/user_manual/processing_algs/qgis/vectorgeneral.html#join-attributes-by-location
coosTileIndex = QgsProject.instance().mapLayersByName('Coos_County_Tile_Index')[0]
elliottTrees = QgsProject.instance().mapLayersByName('ESRF_TreesD_H10Cr20h10A50MD7s')[0]
gridIDparameters = { 'INPUT': elliottTrees, 
                     'JOIN': coosTileIndex,
                     'JOIN_FIELDS': 'gridID',
                     'PREDICATE': 5, # within
                     'METHOD': 1, # first matching (one to one)
                     'OUTPUT': QgsProject.instance().readPath('.') + '/GIS/Trees/ESRF_TreesD_H10Cr20h10A50MD7s with grid ID.gpkg' }
gridIDresult = processing.run("qgis:joinattributesbylocation", gridIDparameters) # ~8 minutes
elliottTreesWithGridID = QgsVectorLayer(gridIDparameters['OUTPUT'], 'ESRF_TreesD_H10Cr20h10A50MD7s with grid ID')

elliottStands = QgsProject.instance().mapLayersByName('Elliott Stand Data Feb2022 EPSG6556')[0]
standIDparameters = { 'INPUT': elliottTreesWithGridID, 
                      'JOIN': elliottStands,
                      'JOIN_FIELDS': 'StandID',
                      'PREDICATE': 5, # within
                      'METHOD': 1, # first matching (one to one)
                      'OUTPUT': QgsProject.instance().readPath('.') + '/GIS/Trees/ESRF_TreesD_H10Cr20h10A50MD7s with stand ID.gpkg' }
standIDresult = processing.run("qgis:joinattributesbylocation", standIDparameters)
elliottTreesWithStandID = QgsVectorLayer(standIDresult['OUTPUT'], 'ESRF_TreesD_H10Cr20h10A50MD7s with stand ID')

resourceUnits = QgsProject.instance().mapLayersByName('iLand 100 m')[0]
resourceUnitParameters = { 'INPUT': elliottTreesWithStandID, 
                           'JOIN': resourceUnits,
                           'JOIN_FIELDS': ['id', 'bufferDist'],
                           'PREFIX': 'resourceUniT',
                           'PREDICATE': 5, # within
                           'METHOD': 1, # first matching (one to one)
                           'OUTPUT': QgsProject.instance().readPath('.') + '/GIS/Trees/ESRF_TreesD_H10Cr20h10A50MD7s with IDs.gpkg' }
resourceUnitResult = processing.run("qgis:joinattributesbylocation", resourceUnitParameters) # ~8 minutes
elliottTreesWithIDs = QgsVectorLayer(resourceUnitResult['OUTPUT'], 'ESRF_TreesD_H10Cr20h10A50MD7s with IDs')

QgsProject.instance().addMapLayer(elliottTreesWithIDs, False)
QgsProject.instance().layerTreeRoot().findGroup("TreesD_H10Cr20h10A50MD7s").addLayer(elliottTreesWithIDs)

## individual tree tile reprojection to EPSG:6556, field updates, and export to .csv
# Coos County 2021 LiDAR tiles are 3000 x 3000 feet in an 87 x 60 grid (north-south x east-west).
# In tree IDs, the 0-99,999 range is used as trees' LiDAR segmentation IDs, allowing 21474 tile IDs as prefixes in a 
# 32 bit integer. A semi-informative unique tree ID can therefore be formed as the nine digit number
#   <tile north-south grid index><tile east-west grid index><tree number in tile>
#tileNames = ['TSegD_H10Cr20h10A50MD7_s04140w06930', 'TSegD_H10Cr20h10A50MD7_s04110w07050', 'TSegD_H10Cr20h10A50MD7_s04110w07020', 'TSegD_H10Cr20h10A50MD7_s04110w06930']
#tileTreeIDprefixes = [ 4624, 5023, 4923, 4623 ]
#tileName = 'TSegD_H10Cr20h10A50MD7_s04110w07050'
#tileTreeIDprefix = 100000 * 5023
tileName = 'TSegD_H10Cr20h10A50MD7_s04110w07020'
tileTreeIDprefix = 100000 * 4923
tileInEpsg6557 = QgsProject.instance().mapLayersByName(tileName)[0]
tileInEpsg6556 = processing.run('native:reprojectlayer', {'INPUT': tileInEpsg6557, 'TARGET_CRS': 'EPSG:6556', 'OUTPUT': 'memory:Reprojected'})['OUTPUT']

# change treeID from a real to an integer, change units from feet to meters, remove unneeded fields, and change field names 
# from LiDAR processing to iLand nomenclature
# Coos County 2021 flight tiles are 83.61 ha. DBH is apparently already in cm.
tileDataProvider = tileInEpsg6556.dataProvider()
if tileDataProvider.fieldNameIndex("id") == -1:
    idResult = tileDataProvider.addAttributes([QgsField("id", QVariant.Int)])
    tileInEpsg6556.updateFields()

if tileDataProvider.fieldNameIndex("species") == -1:
    speciesResult = tileDataProvider.addAttributes([QgsField("species", QVariant.String)])
    tileInEpsg6556.updateFields()

idColumnIndex = tileDataProvider.fieldNameIndex("id") 
speciesColumnIndex = tileDataProvider.fieldNameIndex("species") 
xColumnIndex = tileDataProvider.fieldNameIndex("XTOP")
yColumnIndex = tileDataProvider.fieldNameIndex("YTOP")
heightColumnIndex = tileDataProvider.fieldNameIndex("ZTOP")
treeCoordinateUpdates = {}
for tree in tileInEpsg6556.getFeatures():
    treePosition = tree.geometry().asPoint()
    treeCoordinateUpdates[tree.id()] = { idColumnIndex: tileTreeIDprefix + int(tree["treeID"]),
                                         speciesColumnIndex: "psme", # for now, assume all trees are Douglas-firs as species classification has not been run
                                         xColumnIndex: treePosition.x(), 
                                         yColumnIndex: treePosition.y(), 
                                         heightColumnIndex: 0.3048 * tree[heightColumnIndex] }

tileDataProvider.changeAttributeValues(treeCoordinateUpdates)

tileFieldNames = [field.name() for field in tileInEpsg6556.fields()]
tileFieldsToRemove = []
for fieldIndex, fieldName in list(enumerate(tileFieldNames)):
    if fieldName not in ["id", "species", "XTOP", "YTOP", "dbh", "ZTOP"]:
        tileFieldsToRemove.append(fieldIndex)

tileDataProvider.deleteAttributes(tileFieldsToRemove)
tileDataProvider.renameAttributes({ tileDataProvider.fieldNameIndex("XTOP"): 'x',
                                    tileDataProvider.fieldNameIndex("YTOP"): 'y',
                                    tileDataProvider.fieldNameIndex("ZTOP"): 'height' })
tileInEpsg6556.updateFields() 

elliotStands = QgsProject.instance().mapLayersByName('Elliott Stand Data Feb2022 EPSG6556')[0] # EPSG:2992 version of layer has invalid geometry
standIDparameters = { 'INPUT': tileInEpsg6556, 
                      'JOIN': elliotStands,
                      'JOIN_FIELDS': 'StandID',
                      'PREDICATE': 0, # intersects
                      'METHOD': 1, # first matching (one to one)
                      'OUTPUT': 'TEMPORARY_OUTPUT' }
standIDresult = processing.run("qgis:joinattributesbylocation", standIDparameters)
tileWithStandID = standIDresult['OUTPUT']

tileDataProvider = tileWithStandID.dataProvider()
tileDataProvider.renameAttributes({ tileDataProvider.fieldNameIndex("StandID"): 'standID' })
tileWithStandID.updateFields() 

#QgsProject.instance().addMapLayer(tileWithStandID)

# write tile to .csv
csvOptions = QgsVectorFileWriter.SaveVectorOptions()
csvOptions.driverName = 'CSV'
csvOptions.fileEncoding = 'UTF-8'

tileCsvFilePath = os.path.splitext(tileInEpsg6557.source())[0] + '.csv'
QgsVectorFileWriter.writeAsVectorFormatV3(tileWithStandID, tileCsvFilePath, QgsCoordinateTransformContext(), csvOptions)
