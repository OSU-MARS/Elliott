setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/other.R")
warnings()