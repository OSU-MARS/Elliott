# ~26 minutes with 10x10 cross validation and one worker(AMD Zen 3, 4.6 GHz)
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 1)

source("trees/height-diameter/UMCA.R")