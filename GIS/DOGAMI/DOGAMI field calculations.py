## plot properties
# horizon angles
#   r.horizon on 3 m DEM for azimuth 0, 45, ... 315° horizon angles at distance of 150 m = 50 grid cells
#   raster calculator: average all eight horizon angles to topographic shelter layer
# plot elevation, slope, aspect, and location
#   buffer plots by 15 m
#   zonal statistics: mean elevation, mean slope, mean cos(aspect), mean sin(aspect), mean topographic shelter index
#   join zonal stats to plots, copy elevation, slope, aspect, and topographicShelterIndex
#   aspect calculation: rotate arctangent from +x = 0 to north = 0, flip azimuth direction, and clear negative values
#      step1: aspect = 180/pi() * atan2("Zonal Statistics_sinAspectmean", "Zonal Statistics_cosAspectmean")
#        assumes cos(aspect) and sin(aspect) rasters are taken in azimuth coordinates (rather than +x = east), so
#
#                              N  cos(aspect) = 1
#                                 sin(aspect) = 0
#        cos(aspect) =  0   W     E  cos(aspect) = 0
#        sin(aspect) = -1            sin(aspect) = 1
#                              S  cos(aspect) = -1
#                                 sin(aspect) = 0
#
#      step 2: if(aspect < 0, 360 + aspect, aspect)
#        Step 1 produces north = 0, south = ±180, east = 90, and west = -90

## resample raster to smaller cell size: Raster -> Projections -> Warp (gdal_warp)

## aspect raster
cols: 11669
rows: 15454
nodata: -3.4028234663852886e+38
pi / 180: 0.01745329
if("bare earth aspect nearest degree 12 ft EPSG6557@1" >= 0, cos(0.01745329 * "bare earth aspect nearest degree 12 ft EPSG6557@1"), 0/0)
