#.libPaths(.libPaths()[2]) # remove N: for installing packages on C:
#install.packages(c("dplyr", "arrow", "ggplot2", "lidR", "lubridate", "magrittr", "nlstools" "patchwork", "readr", "readxl", "sf", "stringr, "tidyr", "writexl"))

library(arrow)
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(readr)
library(stringr)
library(terra)
library(tibble)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(size = 0.5),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             panel.border = element_blank()))

## translate weather grid .csvs exported from shapefiles in QGIS to ClimateNA location .csv files
# process in ClimateNA: multi-location -> { historical, future } time series + monthly all variables (needed for frost days) + year = 1901-most recent complete year
#                       select input file + specify output file
#                       Start TS
#                       repeat for future time series: climate future of interest (e.g. 13 GCM SSP370), present-2100, select input + specify output file
# ClimateNA bugs: expects Windows carriage return + line feed (readr::write_csv() defaults to Unix line feed only, resulting in ClimateNA processing no locations) 
#                 doesn't pick up file system changes made after being started
#                 if future climate series is selected but not double clicked so that directory chooser has closed it's not seen as selected
#                 pops path not found error and exits if output file not specified
#                 tries to make progress updates on UI thread blocked by output calculations
# output footprint:   81 4 km cells -> 5.5 MB/century
#                   1390 800 m cells -> 129 MB/century, ~5 minutes/century to generate
#                   5196 400 m cells -> 
#                  20440 200 m cells -> 
#                  79600 100 m cells -> 

#weatherCellSizes = c("4 km", "800 m", "400 m", "200 m", "100 m")
weatherCellSize = "100 m"
weatherColumnTypes = cols(name = "c", name4km = "c", name800m = "c", name400m = "c", name200m = "c", .default = "d")
weatherCells = read_csv(paste0("GIS/iLand/grid ", weatherCellSize, ".csv"), col_types = weatherColumnTypes) %>%
  rename(ID1 = id, ID2 = name, lat = latitude, long = longitude, el = elevationInMmean) %>%
  select(ID1, ID2, lat, long, el)
write_csv(weatherCells, paste0("../data/ClimateNA/Elliott ", weatherCellSize, ".csv"), eol = "\r\n")


## translate ClimateNA time series .csvs to iLand .csv or .feather files
# Since radiation columns (Rad*) are NA (-9999) in observed weather they aren't loaded, resulting in bind_rows() filling them with NA.
# tidyselect conflates construction and use, requiring inline col_select values. summarize_all() then collapses the overlapping years
# in favor of the observed values, replacing NA radiation with modeled values in years where it is available. It's assumed observed
# years are loaded first and that observed temperature and precipitation values therefore override GCM estimated ones.
#
# ClimateNA prefixes: Tmin, Tave, Tmax - mean minimum daily, mean, and mean maximum daily temperatures (°C)
#                     PPT, RAD - total monthly precipitation (mm) and mean daily solar radiation (MJ/m²-d)
#                     PAS - total monthly precipitation as snow (mm)
#
# uncompressed .feather footprints for monthly weather, 2011-2100: 36 MB/cell-century
#   cell size  cells (4 km buffer)  file size, MB
#   4 km           81                  2.6
#   800 m        1890                 60
#   400          5196                163
#   200        20,440                640
#   100        79,628               2370
climateNAcolumnTypes = cols(ID2 = "c", .default = "d")
for (weatherCellSize in c("100 m", "200 m", "400 m", "800 m", "4 km"))
{
  # read hindcast file first and then the projected file
  # This produces a mostly chronologically ordered tibble from (with ClimateNA default settings) 1901 through 2100 with an overlap
  # zone in the middle. Grouped by year and climate cell, groups have either one row (hindcast only or future prediction only), or
  # two rows. Since ClimateNA emits NA values for mean daily solar radiation in hindcast files those columns (Rad*) are not read 
  # from the hindcast file. This causes bind_rows() to fill them with NA.
  monthlyWeather = bind_rows(read_csv(file.path(getwd(), paste0("../data/ClimateNA/Elliott ", weatherCellSize, "_1901-2021M.csv")), col_select = c("Year", "ID2", starts_with("Tmax"), starts_with("Tmin"), starts_with("Tave"), starts_with("PPT"), starts_with("PAS"), starts_with("RH")), col_types = climateNAcolumnTypes) %>% mutate(source = "observed"),
                             read_csv(file.path(getwd(), paste0("../data/ClimateNA/Elliott ", weatherCellSize, "_13GCMs_ensemble_ssp370_2011-2100M.csv")), col_select = c("Year", "ID2", starts_with("Tmax"), starts_with("Tmin"), starts_with("Tave"), starts_with("PPT"), starts_with("RAD"), starts_with("PAS"), starts_with("RH")), col_types = climateNAcolumnTypes) %>% mutate(source = "13GCMssp370"))
  monthlyWeatherOverlap = monthlyWeather %>% filter(Year >= 2011, Year <= 2021) %>% 
    group_by(Year, ID2) %>% 
    summarize(across(.cols = !starts_with("Rad"), ~ .x[1]), # prefer hindcast (observed) value for all but solar radiaton columns
              across(.cols = starts_with("Rad"), ~ .x[2]),
              .groups = "drop")
  #  summarize_all(~first(na.omit(.x))) # deprecated and slow compared to across() even at 4 km and 800 m resolution, so only run for years where weather observations and predictions overlap
  #as.data.frame(monthlyWeather %>% filter(ID2 == "N192.0 E104.0", Year == 2011)) # diagnostics for checking merge correctness
  #as.data.frame(monthlyWeatherOverlap %>% filter(ID2 == "N192.0 E104.0", Year == 2011))
  
  # ClimateNA and iLand units mostly agree but some conversion is needed
  #                                              ClimateNA   iLand
  # monthly mean daily temperatures              °C          °C
  # monthly total precipitation (rain and snow)  mm          mm
  # relative humidity                            %           % -> vapor pressure deficit (debatable if should be pre-transformed in R)
  # solar radiation                              MJ/m²-day   MJ/m²
  monthlyWeather = bind_rows(monthlyWeatherOverlap, monthlyWeather %>% filter(Year >= 2022)) %>%
    select(-source) %>% 
    mutate(Rad01 = 31 * Rad01, # convert solar radiation from ClimateNA's mean MJ/m²-day to iLand's monthly total MJ/m²
           Rad02 = if_else(leap_year(Year), 29, 28) * Rad02,
           Rad03 = 31 * Rad03,
           Rad04 = 30 * Rad04,
           Rad05 = 31 * Rad05,
           Rad06 = 30 * Rad06,
           Rad07 = 31 * Rad07,
           Rad08 = 31 * Rad08,
           Rad09 = 30 * Rad09,
           Rad10 = 31 * Rad10,
           Rad11 = 30 * Rad11,
           Rad12 = 31 * Rad12) %>%
    arrange(ID2, Year) # group time series by ID2 to assist read performance rather than ClimateNA's grouping by year

  if (weatherCellSize == "4 km")
  {
    write_csv(monthlyWeather, paste0("iLand/database/weather ", weatherCellSize, " 2011-2100 13GCMssp370.csv"))
  }
  monthlyArrow = arrow_table(monthlyWeather, schema = schema(Year = int16(), ID2 = string(), 
                                                             Tmax01 = float32(), Tmax02 = float32(), Tmax03 = float32(), Tmax04 = float32(), Tmax05 = float32(), Tmax06 = float32(), Tmax07 = float32(), Tmax08 = float32(), Tmax09 = float32(), Tmax10 = float32(), Tmax11 = float32(), Tmax12 = float32(),
                                                             Tmin01 = float32(), Tmin02 = float32(), Tmin03 = float32(), Tmin04 = float32(), Tmin05 = float32(), Tmin06 = float32(), Tmin07 = float32(), Tmin08 = float32(), Tmin09 = float32(), Tmin10 = float32(), Tmin11 = float32(), Tmin12 = float32(),
                                                             Tave01 = float32(), Tave02 = float32(), Tave03 = float32(), Tave04 = float32(), Tave05 = float32(), Tave06 = float32(), Tave07 = float32(), Tave08 = float32(), Tave09 = float32(), Tave10 = float32(), Tave11 = float32(), Tave12 = float32(),
                                                             PPT01 = float32(), PPT02 = float32(), PPT03 = float32(), PPT04 = float32(), PPT05 = float32(), PPT06 = float32(), PPT07 = float32(), PPT08 = float32(), PPT09 = float32(), PPT10 = float32(), PPT11 = float32(), PPT12 = float32(),
                                                             PAS01 = float32(), PAS02 = float32(), PAS03 = float32(), PAS04 = float32(), PAS05 = float32(), PAS06 = float32(), PAS07 = float32(), PAS08 = float32(), PAS09 = float32(), PAS10 = float32(), PAS11 = float32(), PAS12 = float32(),
                                                             RH01 = float32(), RH02 = float32(), RH03 = float32(), RH04 = float32(), RH05 = float32(), RH06 = float32(), RH07 = float32(), RH08 = float32(), RH09 = float32(), RH10 = float32(), RH11 = float32(), RH12 = float32(),
                                                             Rad01 = float32(), Rad02 = float32(), Rad03 = float32(), Rad04 = float32(), Rad05 = float32(), Rad06 = float32(), Rad07 = float32(), Rad08 = float32(), Rad09 = float32(), Rad10 = float32(), Rad11 = float32(), Rad12 = float32()))
  write_feather(monthlyArrow, paste0("iLand/database/weather ", weatherCellSize, " 2011-2100 13GCMssp370.feather"), compression = "uncompressed") # work around lack of compression support in Arrow 8.0.0 C#
}

#library(fst)
#write_fst(monthlyWeather, paste0("iLand/database/weather ", weatherCellSize, " 2011-2100 13GCMssp370.fst"))

## translate resource unit .csv exported from QGIS to iLand resource unit .csv and .feather as a function of weather
# cell size
maxBufferDistance = 1000
maxBufferDistanceName = paste(0.001 * maxBufferDistance, "km")
resourceUnitColumnTypes = cols(name = "c", name4km = "c", name800m = "c", name400m = "c", name200m = "c", .default = "d")
for (weatherCellSize in c("4 km", "800 m", "400 m", "200 m", "100 m"))
{
  #start = Sys.time() # using geospatial file formats is slow as of terra 1.5-34 (https://github.com/rspatial/terra/issues/745)
  #resourceUnits = as_tibble(vect("GIS/iLand/grid 100 m.gpkg")) # 6.8 s total, 1.3 s vect() => 5.5 s as_tibble()
  #resourceUnits = as_tibble(vect("GIS/iLand/grid 100 m.fgb")) # 7.8 s total, 2.1 s vect() => 5.7 s as_tibble()
  #(elapsed = Sys.time() - start)
  resourceUnits = read_csv("GIS/iLand/grid 100 m.csv", col_types = resourceUnitColumnTypes) %>% # 0.17 seconds
    filter(bufferDist <= maxBufferDistance) %>%
    mutate(centerX = 0.5 * (left + right), # m
           centerY = 0.5 * (top + bottom), # m
           soilThickness0_5mean = pmax(0, pmin(5, soilDepthInCMmean)), # cm
           soilThickness5_15mean = pmax(0, pmin(10, soilDepthInCMmean - 5)), # cm
           soilThickness15_30mean = pmax(0, pmin(15, soilDepthInCMmean - 15)), # cm
           soilThickness30_60mean = pmax(0, pmin(30, soilDepthInCMmean - 30)), # cm
           soilThickness60_100mean = pmax(0, pmin(40, soilDepthInCMmean - 60)), # cm
           soilThickness100plusMean = pmax(0, soilDepthInCMmean - 100), # cm
           soilEffectiveThickness0_5mean = (1 - 0.001 * coarseFragments0_5mean) * soilThickness0_5mean, # cm
           soilEffectiveThickness5_15mean = (1 - 0.001 * coarseFragments5_15mean) * soilThickness5_15mean, # cm
           soilEffectiveThickness15_30mean = (1 - 0.001 * coarseFragments15_30mean) * soilThickness15_30mean, # cm
           soilEffectiveThickness30_60mean = (1 - 0.001 * coarseFragments30_60mean) * soilThickness30_60mean, # cm
           soilEffectiveThickness60_100mean = (1 - 0.001 * coarseFragments60_100mean) * soilThickness60_100mean, # cm
           soilEffectiveThickness100plusMean = (1 - 0.001 * coarseFragments100_200mean) * soilThickness100plusMean, # cm
           soilEffectiveDepthMean = soilEffectiveThickness0_5mean + soilEffectiveThickness5_15mean + soilEffectiveThickness15_30mean + soilEffectiveThickness30_60mean + soilEffectiveThickness60_100mean + soilEffectiveThickness100plusMean) %>% # cm
    #filter(centerX > 125200, centerX < 126200, centerY > 213900, centerY < 215800) %>% # optional: window resource units
    mutate(soilThetaR = 1 / soilEffectiveDepthMean * (soilEffectiveThickness0_5mean * vanGenuchtenThetaR0_5mean +
                                                      soilEffectiveThickness5_15mean * vanGenuchtenThetaR5_15mean + 
                                                      soilEffectiveThickness15_30mean * vanGenuchtenThetaR15_30mean + 
                                                      soilEffectiveThickness15_30mean * vanGenuchtenThetaR30_60mean + 
                                                      soilEffectiveThickness60_100mean * vanGenuchtenThetaR60_100mean + 
                                                      soilEffectiveThickness100plusMean * vanGenuchtenThetaR60_100mean), # m³/m³
           soilThetaS = soilThetaR + 1 / soilEffectiveDepthMean * (soilEffectiveThickness0_5mean * (vanGenuchtenThetaS0_5mean - vanGenuchtenThetaR0_5mean) +
                                                                   soilEffectiveThickness5_15mean * (vanGenuchtenThetaS5_15mean - vanGenuchtenThetaR5_15mean) +
                                                                   soilEffectiveThickness15_30mean * (vanGenuchtenThetaS15_30mean - vanGenuchtenThetaR15_30mean) +
                                                                   soilEffectiveThickness15_30mean * (vanGenuchtenThetaS30_60mean - vanGenuchtenThetaR30_60mean) +
                                                                   soilEffectiveThickness60_100mean * (vanGenuchtenThetaS60_100mean - vanGenuchtenThetaR60_100mean) +
                                                                   soilEffectiveThickness100plusMean * (vanGenuchtenThetaS100_200mean - vanGenuchtenThetaR100_200mean)), # m³/m³
           soilVanGenuchtenAlpha = 1 / soilEffectiveDepthMean * (soilEffectiveThickness0_5mean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm0_5mean +
                                                                 soilEffectiveThickness5_15mean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm5_15mean + 
                                                                 soilEffectiveThickness15_30mean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm15_30mean + 
                                                                 soilEffectiveThickness15_30mean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm30_60mean + 
                                                                 soilEffectiveThickness60_100mean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm60_100mean + 
                                                                 soilEffectiveThickness100plusMean * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm100_200mean), # kPa
           soilVanGenuchtenN = 1 / soilEffectiveDepthMean * (soilEffectiveThickness0_5mean * vanGenuchtenN0_5mean +
                                                             soilEffectiveThickness5_15mean * vanGenuchtenN5_15mean + 
                                                             soilEffectiveThickness15_30mean * vanGenuchtenN15_30mean + 
                                                             soilEffectiveThickness15_30mean * vanGenuchtenN30_60mean + 
                                                             soilEffectiveThickness60_100mean * vanGenuchtenN60_100mean + 
                                                             soilEffectiveThickness100plusMean * vanGenuchtenN100_200mean), # exponent
           weatherID = case_when(weatherCellSize == "4 km" ~ name4km, weatherCellSize == "800 m" ~ name800m, weatherCellSize == "400 m" ~ name400m, weatherCellSize == "200 m" ~ name200m, TRUE ~ name)) %>%
    rename(soilPlantAccessibleDepth = soilEffectiveDepthMean) %>%
    select(id, centerX, centerY, weatherID, soilPlantAccessibleDepth, soilThetaS, soilThetaR, soilVanGenuchtenAlpha, soilVanGenuchtenN) %>%
    mutate(# id = ,
           # speciesTable = ,
           # centerX = ,
           # centerY = ,
           # weatherID = ,
           # snagBranchRootC = ,
           # snagBranchRootCN = ,
           # snagCarbon = ,
           # snagCNRatio = ,
           # snagDecompositionRate = ,
           # snagHalfLife = ,
           # snagCount = ,
           # soilAnnualNitrogenDeposition = ,
           # soilAvailableNitrogen = ,
           # soilDepth = ,
           # soilEl = ,
           # soilEr = ,
           # soilLeaching = ,
           # soilHumificationRate = ,
           # soilOrganicC = ,
           # soilOrganicDecompositionRate = ,
           # soilOrganicN = ,
           # soilClayPercent = ,
           # soilSandPercent = ,
           # soilSiltPercent = ,
           # soilQh = ,
           # soilThetaR = ,
           # soilThetaS = ,
           # soilVanGenuchtenAlpha = ,
           # soilVanGenuchtenN = ,
           # soilYoungLabileC = ,
           # soilYoungLabileDecompositionRate = ,
           # soilYoungLabileN = ,
           # soilYoungRefractoryC = ,
           # soilYoungRefractoryDecompositionRate = ,
           # soilYoungRefractoryN = ,
           ) %>%
    arrange(centerY, centerX)
  if (weatherCellSize == "4 km")
  {
    write_csv(resourceUnits, paste0("iLand/gis/resource units ", maxBufferDistanceName, " buffer ", weatherCellSize, " weather.csv"))
  }
  resourceUnitsArrow = arrow_table(resourceUnits, schema = schema(id = int32(), 
                                                                  centerX = float32(), centerY = float32(), 
                                                                  weatherID = string(),
                                                                  soilPlantAccessibleDepth = float32(), soilThetaS = float32(), soilThetaR = float32(), soilVanGenuchtenAlpha = float32(), soilVanGenuchtenN = float32()))
  write_feather(resourceUnitsArrow, paste0("iLand/gis/resource units ", maxBufferDistanceName, " buffer ", weatherCellSize, " weather.feather"), compression = "uncompressed")
}

resourceUnitsWithPawc = resourceUnits %>% 
  mutate(plantAccessibleWater = 10 * soilPlantAccessibleDepth * (soilThetaS - soilThetaR))
pawcStats = as_tibble_row(quantile(resourceUnitsWithPawc$plantAccessibleWater, probs = c(0.05, 0.20, 0.50, 0.80, 0.95))) %>%
            mutate(mean = mean(resourceUnitsWithPawc$plantAccessibleWater))

ggplot() +
  geom_segment(aes(x = `5%`, y = 0, xend = `5%`, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = `20%`, y = 0, xend = `20%`, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = mean, y = 0, xend = mean, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = `50%`, y = 0, xend = `50%`, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = `80%`, y = 0, xend = `80%`, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = `95%`, y = 0, xend = `95%`, yend = 10000), pawcStats, color = "grey70", linetype = "longdash") +
  geom_histogram(aes(x = plantAccessibleWater), resourceUnitsWithPawc, binwidth = 20) +
  coord_cartesian(ylim = c(0, 5000)) +
  labs(x = "plant available soil water at field capacity, mm", y = "resource units")

# generate windowed resource units and matching windowed weather file for iLand unit tests
# Set resourceUnits and resourceUnitsArrow with loop body above after uncommenting the filter() statement for weather 
# windowing. The should find 190 resource units.
weatherCellSize = "200 m"
write_feather(resourceUnitsArrow, paste0("iLand/gis/unit test resource units ", weatherCellSize, " weather.feather"), compression = "uncompressed")
monthlyWeather = read_feather(paste0("iLand/database/weather ", weatherCellSize, " 2011-2100 13GCMssp370.feather")) %>%
  filter(ID2 %in% resourceUnits$weatherID)
monthlyArrow = arrow_table(monthlyWeather, schema = schema(Year = int16(), ID2 = string(), 
                                                           Tmax01 = float32(), Tmax02 = float32(), Tmax03 = float32(), Tmax04 = float32(), Tmax05 = float32(), Tmax06 = float32(), Tmax07 = float32(), Tmax08 = float32(), Tmax09 = float32(), Tmax10 = float32(), Tmax11 = float32(), Tmax12 = float32(),
                                                           Tmin01 = float32(), Tmin02 = float32(), Tmin03 = float32(), Tmin04 = float32(), Tmin05 = float32(), Tmin06 = float32(), Tmin07 = float32(), Tmin08 = float32(), Tmin09 = float32(), Tmin10 = float32(), Tmin11 = float32(), Tmin12 = float32(),
                                                           Tave01 = float32(), Tave02 = float32(), Tave03 = float32(), Tave04 = float32(), Tave05 = float32(), Tave06 = float32(), Tave07 = float32(), Tave08 = float32(), Tave09 = float32(), Tave10 = float32(), Tave11 = float32(), Tave12 = float32(),
                                                           PPT01 = float32(), PPT02 = float32(), PPT03 = float32(), PPT04 = float32(), PPT05 = float32(), PPT06 = float32(), PPT07 = float32(), PPT08 = float32(), PPT09 = float32(), PPT10 = float32(), PPT11 = float32(), PPT12 = float32(),
                                                           PAS01 = float32(), PAS02 = float32(), PAS03 = float32(), PAS04 = float32(), PAS05 = float32(), PAS06 = float32(), PAS07 = float32(), PAS08 = float32(), PAS09 = float32(), PAS10 = float32(), PAS11 = float32(), PAS12 = float32(),
                                                           RH01 = float32(), RH02 = float32(), RH03 = float32(), RH04 = float32(), RH05 = float32(), RH06 = float32(), RH07 = float32(), RH08 = float32(), RH09 = float32(), RH10 = float32(), RH11 = float32(), RH12 = float32(),
                                                           Rad01 = float32(), Rad02 = float32(), Rad03 = float32(), Rad04 = float32(), Rad05 = float32(), Rad06 = float32(), Rad07 = float32(), Rad08 = float32(), Rad09 = float32(), Rad10 = float32(), Rad11 = float32(), Rad12 = float32()))
write_feather(monthlyArrow, paste0("iLand/database/unit test weather ", weatherCellSize, " 2011-2100 13GCMssp370.feather"), compression = "uncompressed") # work around lack of compression support in Arrow 8.0.0 C#


resourceUnits %>% summarize(thetaRmin = min(soilThetaR), thetaRmax = max(soilThetaR),
                            thetaSmin = min(soilThetaS), thetaSmax = max(soilThetaS),
                            alphaMin = min(soilVanGenuchtenAlpha), alphaMax = max(soilVanGenuchtenAlpha),
                            nMin = min(soilVanGenuchtenN), nMax = max(soilVanGenuchtenN))
resourceUnits %>% filter(is.na(soilThetaR))

## CO₂ concentrations: Mauna Loa reference and SSP scenarios
# Mauna Loa: https://gml.noaa.gov/ccgg/trends/data.html
# NOAA Greenhouse Gas Marine Boundary Layer Reference: https://gml.noaa.gov/ccgg/mbl/data.php -> zonal -> custom -> 43.46-43.69 N, 1979-2020, open in Excel and 
# SSP: https://iiasa.ac.at/models-and-data/shared-socioeconomic-pathways-scenario-database for general data
#      https://greenhousegases.science.unimelb.edu.au/#!/ghg?mode=downloads&scenarioid=2 for CO₂ concentrations, monthly 0.5°
# maunaLoa = read_csv("Elliott/iLand/database/co2_mm_mlo.csv") # not needed if marine boundary layer used
marineBoundaryColumnTypes = cols(year = "i", month = "i", day = "i", decimal_date = "d", value = "d", uncertainty = "d") # uncertainty is standard deviation
marineBoundaryLayerReference = read_table("iLand/database/co2_GHGreference.2010931543_zonal.txt", col_types = marineBoundaryColumnTypes, skip = 74) %>%
  rename(co2 = value)
sspColumnTypes = cols(Var1 = "i", year = "i", month = "i", data = "d", datetime = "c", datenum = "i")
ssp370 = read_csv("iLand/database/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr-0p5x360deg_201501-250012.csv", col_types = sspColumnTypes) %>%
  filter(lat >= 43, lat <= 44) %>%
  select(year, month, lat, data) %>%
  pivot_wider(id_cols = c("year", "month"), names_prefix = "co2_", names_from = "lat", values_from = "data") %>%
  mutate(co2 = 0.5 * (co2_43.25 + co2_43.75), 
         decimal_date = year + (month - 0.5) / 12)

ssp370timeSeries = ssp370 %>% filter(year >= 2011, year <= 2100) %>% select(-co2_43.25, -co2_43.75, -decimal_date)
write_csv(ssp370timeSeries, "iLand/database/co2 ssp370.csv")

ssp370arrow = arrow_table(ssp370timeSeries, schema = schema(year = int16(), month = uint8(), co2 = float32()))
write_feather(ssp370arrow, "iLand/database/co2 ssp370.feather", compression = "uncompressed")

ggplot(marineBoundaryLayerReference) +
  # geom_ribbon(aes(x = decimal_date, ymin = value - 1.96 * uncertainty, ymax = value + 1.96 * uncertainty), alpha = 0.1, fill = "grey50") + # treated as normally distributed, see docs, but negligible for forest modeling since range is <±1 ppm
  geom_path(aes(x = decimal_date, y = co2)) +
  labs(x = "date", y = bquote("marine boundary layer atmospheric CO"[2]*" concentration, ppm"))
ggplot(ssp370) + 
  geom_path(aes(x = decimal_date, y = co2)) +
  labs(x = "date", y = bquote("SSP370 atmospheric CO"[2]*" concentration, ppm"))

print(ssp370, n = 100)

## ClimateNA outputs
weatherCellSize = "800 m"
monthlyWeather = read_feather(paste0("iLand/database/weather ", weatherCellSize, " 2011-2100 13GCMssp370.feather"), mmap = FALSE) %>%
  pivot_longer(cols = starts_with(c("Tmax", "Tmin", "Tave", "PPT", "PAS", "RH", "Rad")), names_pattern = "([a-zA-Z]+)(\\d+)", names_to = c(".value", "month")) %>%
  filter(Rad >= 0) %>% # exclude six tiles along coast, N221.6-N225.6 E102.4, where ClimateNA generates invalid solar radiation values
  mutate(month = as.integer(month),
         decimalYear = Year + (month - 0.5) / 12)

# monthly density
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = month, y = Tmax), binwidth = c(1, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean max temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = month, y = Tave), binwidth = c(1, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean mean temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = month, y = Tmin), binwidth = c(1, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean min temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = month, y = PPT), binwidth = c(1, 10)) +
  labs(x = "month", y = "monthly precipitation, mm", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = month, y = PAS), binwidth = c(1, 1)) +
  labs(x = "month", y = "monthly snow, mm SWE", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.position = "none") +
  ggplot(monthlyWeather) +
geom_bin_2d(aes(x = month, y = RH), binwidth = c(1, 1)) +
  coord_cartesian(ylim = c(0, 100)) +
  guides(fill = guide_colorbar(direction = "horizontal", title.position = "top")) +
  labs(x = "month", y = "mean relative humidity, %", fill = paste(weatherCellSize, "weather cell years")) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(1, 12, by = 2)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.02)) +
plot_layout(ncol = 3, nrow = 2)

# quarterly time series density
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = Tmax, fill = ..density..), binwidth = c(3/12, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean max temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = Tave, fill = ..density..), binwidth = c(3/12, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean mean temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = Tmin, fill = ..density..), binwidth = c(3/12, 0.5)) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean min temp, °C", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = PPT, fill = ..density..), binwidth = c(3/12, 10)) +
  labs(x = "year", y = "monthly precipitation, mm", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = PAS, fill = ..density..), binwidth = c(3/12, 1)) +
  labs(x = "year", y = "monthly snow, mm SWE", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
ggplot(monthlyWeather) +
  geom_bin_2d(aes(x = decimalYear, y = RH, fill = ..density..), binwidth = c(3/12, 1)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "year", y = "mean relative humidity, %", fill = NULL) +
  scale_fill_viridis_c(trans = "log10") +
  theme(legend.position = "none") +
plot_layout(ncol = 3, nrow = 2)
#ggsave(paste0("Presentation/Climate NA temp precip snow RH ", weatherCellSize, " monthly.png"), width = 3 * 10, height = 2 * 10, units = "cm", dpi = 200)
#ggsave(paste0("Presentation/Climate NA temp precip snow RH ", weatherCellSize, " monthly small.png"), width = 3 * 5, height = 2 * 5, units = "cm", dpi = 200)

ggplot(monthlyWeather) + # not particularly informative and takes several minutes to render but the binned version conveys even less information
  geom_path(aes(x = decimalYear, y = Rad, color = "radiation", group = ID2), alpha = 0.1) +
  labs(x = "year", y = "solar radiation, MJ m ", color = NULL) +
  scale_color_manual(breaks = c("radiation"), values = c("gold2")) +
  theme(legend.position = "none")

# annual trends
annualWeather = monthlyWeather %>% group_by(Year, ID2) %>%
  summarize(ID2 = ID2[1], Tmax = mean(Tmax), Tave = mean(Tave), Tmin = mean(Tmin), PPT = sum(PPT), PAS = sum(PAS), RH = mean(RH), Rad = sum(Rad), .groups = "drop")
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = Tmax, color = "Tmax", group = ID2), alpha = 0.1) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean max temp, °C", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = Tave, color = "Tave", group = ID2), alpha = 0.1) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean mean temp, °C", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = Tmin, color = "Tmin", group = ID2), alpha = 0.1) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = NULL, y = "daily mean min temp, °C", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = PPT, color = "precip", group = ID2), alpha = 0.1) +
  labs(x = "year", y = "annual precipitation, mm", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = PAS, color = "snow", group = ID2), alpha = 0.1) +
  labs(x = "year", y = "annual snow, mm SWE", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
ggplot(annualWeather) +
  geom_path(aes(x = Year, y = RH, color = "relative humidity", group = ID2), alpha = 0.1) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "year", y = "mean relative humidity, %", color = NULL) +
  scale_color_manual(breaks = c("Tmax", "Tave", "Tmean", "precip", "snow", "relative humidity"), values = c("red2", "orangered2", "darkorange2", "blue2", "cyan1", "grey50")) +
  theme(legend.position = "none") +
plot_layout(ncol = 3, nrow = 2)
#ggsave(paste0("Presentation/Climate NA temp precip snow RH ", weatherCellSize, " yearly.png"), width = 3 * 10, height = 2 * 10, units = "cm", dpi = 200)
#ggsave(paste0("Presentation/Climate NA temp precip snow RH ", weatherCellSize, " yearly small.png"), width = 3 * 5, height = 2 * 5, units = "cm", dpi = 200)

ggplot(annualWeather) +
  geom_path(aes(x = Year, y = Rad, color = "radiation", group = ID2), alpha = 0.1) +
  labs(x = "year", y = "annual solar radiation, MJ m⁻²", color = NULL) +
  scale_color_manual(breaks = c("radiation"), values = c("gold2")) +
  theme(legend.position = "none")

