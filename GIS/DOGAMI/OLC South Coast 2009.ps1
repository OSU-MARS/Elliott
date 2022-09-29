$lastools = $env:USERPROFILE + "\PhD\tools\LAStools\bin"
$lasboundary = $lastools + "\lasboundary.exe"
$southCoast = $env:USERPROFILE + "\PhD\Elliott\GIS\DOGAMI\2009 OLC South Coast"

# index tiles since OLC doesn't distribute an index (tiling is distinct from PROCESSING_BINS .shp)
# Download tiles { 43123+43124, D+E+F } from ftp://lidar.engr.oregonstate.edu/OREGON%20LIDAR%20CONSORTIUM%20PROJECT%20DATA/OLC%20SOUTH%20COAST%202009/POINTS/
# Extract bounding boxes as shapefiles with LAStools' lasboundary.exe.
# Bulk add individual tiles' shapefiles in QGIS and merge (Vector -> Data Management Tools -> Merge Vector Layers...)
$inputFiles = """$southCoast\Points\43123*.laz""" # 
$outputDirectory = """$southCoast\Tiles"""
&$lasboundary -i $inputFiles -convex -cpu64 -cores 14 -odir $outputDirectory -oshp -v