library(dplyr)
library(readr)
library(sf)
library(stringr)

## load shapefile and map unit descriptors
# shapefile mukey -> component cokey -> horizon
# https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_053631 -> SSURGO metadata
countyShapefile = st_read(file.path(getwd(), "GIS/SSURGO/Coos County (OR011)/spatial/soilmu_a_or011.shp"))

componentColumnNames = c("comppct_l", "comppct_r", "comppct_h", "compname", "compkind", "majcompflag", "otherph", "localphase", "slope_l", "slope_r", "slope_h", "slopelenusle_l", "slopelenusle_r", "slopelenusle_h", "runoff", "tfact", "wei", "weg", "erocl", "earthcovkind1", "earthcovkind2", "hydricon", "hydricrating", "drainagecl", "elev_l", "elev_r", "elev_h", "aspectccwise", "aspectrep", "aspectcwise", "geomdesc", "albedodry_l", "albedodry_r", "albedodry_h", "airtempa_l", "airtempa_r", "airtempa_h", "map_l", "map_r", "map_h", "reannualprecip_l", "reannualprecip_r", "reannualprecip_h", "ffd_l", "ffd_r", "ffd_h", "nirrcapcl", "nirrcapscl", "nirrcapunit", "irrcapcl", "irrcapscl", "irrcapunit", "cropprodindex", "constreeshrubgrp", "wndbrksuitgrp", "rsprod_l", "rsprod_r", "rsprod_h", "foragesuitgrpid", "wlgrain", "wlgrass", "wlherbaceous", "wlshrub", "wlconiferous", "wlhardwood", "wlwetplant", "wlshallowwat", "wlrangeland", "wlopenland", "wlwoodland", "wlwetland", "soilslippot", "frostact", "initsub_l", "initsub_r", "initsub_h", "totalsub_l", "totalsub_r", "totalsub_h", "hydgrp", "corcon", "corsteel", "taxclname", "taxorder", "taxsuborder", "taxgrtgroup", "taxsubgrp", "taxpartsize", "taxpartsizemod", "taxceactcl", "taxreaction", "taxtempcl", "taxmoistscl", "taxtempregime", "soiltaxedition", "castorieindex", "flecolcomnum", "flhe", "flphe", "flsoilleachpot", "flsoirunoffpot", "fltemik2use", "fltriumph2use", "indraingrp", "innitrateleachi", "misoimgmtgrp", "vasoimgtgrp", "mukey", "cokey")
componentColumnTypes = cols_only(comppct_r = "d", compname = "c", drainagecl = "c", hydgrp = "c", taxorder = "c", taxsuborder = "c", mukey = "i", cokey = "i")
countyComponents = read_delim("GIS/SSURGO/Coos County (OR011)/tabular/comp.txt", col_names = componentColumnNames, col_types = componentColumnTypes, delim = "|")

mapUnitAggregatedAttributeColumns = c("musym", "muname", "mustatus", "slopegraddcp", "slopegradwta", "brockdepmin", "wtdepannmin", "wtdepaprjunmin", "flodfreqdcd", "flodfreqmax", "pondfreqprs", "aws025wta", "aws050wta", "aws100wta", "aws150wta", "drclassdcd", "drclasswettest", "hydgrpdcd", "iccdcd", "iccdcdpct", "niccdcd", "niccdcdpct", "engdwobdcd", "engdwbdcd", "engdwbll", "engdwbml", "engstafdcd", "engstafll", "engstafml", "engsldcd", "engsldcp", "englrsdcd", "engcmssdcd", "engcmssmp", "urbrecptdcd", "urbrecptwta", "forpehrtdcp", "hydclprs", "awmmfpwwta", "mukey") # muagatt
componentColumnTypes = cols_only(brockdepmin = "i", wtdepannmin = "i", aws025wta = "d", aws050wta = "d", aws100wta = "d", aws150wta = "d", forpehrtdcp = "c", mukey = "i")
