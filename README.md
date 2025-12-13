### Overview
This repo contains scripts for manipulating data pertaining to the Elliott State Research Forest (previously the Elliott State Forest) on the Oregon coast. Since it doesn't contain the data (that's not ours to release, sorry) it's unlikely to be of much interest unless you're on one of the teams working on the Elliott. If you are on such a team, however, the usual benefits of revision control apply.

If you're looking for LiDAR tiles covering the Elliott those are available from [DOGAMI](https://www.oregongeology.org/) with [GEO](https://www.oregon.gov/GEO/Pages/index.aspx) providing other state-level GIS data. Other open data sources—such as MODIS, POLARIS, and SSURGO—are linked in files which use them.

If you're visiting for the source code of [West and Strimbu 2025](https://doi.org/10.1093/forestry/cpaf010) that is here. However,
starting from the [McDonald-Dunn](../McDonald-Dunn) code is suggested for establishing new height and diameter modeling efforts. It is essentially a v1.1 to the Elliott's 1.0.

### Dependencies
[R](https://www.r-project.org/) is the primary tool used here, mainly via [RStudio](https://www.rstudio.com/) Desktop, and thus most code is in .R files. GIS processing is done mainly with [sf](https://r-spatial.github.io/sf/) for vector data and [terra](https://rspatial.github.io/terra/) for rasters. [QGIS](https://qgis.org/), [GDAL](https://gdal.org/), and [GRASS](https://grass.osgeo.org/) are also used at times.