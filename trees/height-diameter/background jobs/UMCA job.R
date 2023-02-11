setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 1)

source("trees/height-diameter/UMCA.R")