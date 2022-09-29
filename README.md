### Overview
This repo contains scripts for manipulating data pertaining to the Elliott State Research Forest (previously the Elliot State Forest) on the Oregon coast. Since it doesn't contain the data (that's not ours to release, sorry) it's unlikely to be of much interest unless you're on one of the teams working on the Elliott. If you are on such a team, however, the usual benefits of revision control apply.

If you're looking LiDAR tiles covering the Elliott those are available from [DOGAMI](https://www.oregongeology.org/) with [GEO](https://www.oregon.gov/GEO/Pages/index.aspx) providing other state-level GIS data. Other open data sources—such as MODIS, POLARIS, and SSURGO—are linked in files which use them.

### Dependencies
[R](https://www.r-project.org/) is the primary tool used here, mainly via [RStudio](https://www.rstudio.com/) Desktop, and thus most code is in .R files. GIS processing is done mainly with [QGIS](https://qgis.org/), either QGIS proper or its associated installs of [GDAL](https://gdal.org/), [GRASS](https://grass.osgeo.org/). Python code (.py) is therefore likely intended to run in [PyQGIS](https://docs.qgis.org/3.22/en/docs/pyqgis_developer_cookbook/index.html).