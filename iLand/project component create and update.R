# R bits for maintenance of project files
#
# For project setup (input file generation)
# - Elliott/GIS/POLARIS.R: POLARIS 30 m soil grid -> iLand resource unit van Genuchten plant available water
#               weather.R: QGIS -> ClimateNA grid, ClimateNA monthly time series -> iLand .feather
#               trees.R: LiDAR treetops + classifications -> iLand .feather with { x, y, species, height, DBH }
# - Elliott/iLand/database/{CO2}: 
# - iLand/UnitTests/Elliott/database/species_param_pnw.sqlite: Pacific Northwest species parameterization (initial; Seidl et al. 2010)
# - iLand/UnitTests/Elliott/lip/*.feather: light intensity profiles for Pacific Northwest species (initial; Seidl et al. 2010)
