# ~115 minutes with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
#  height fits ~14 minutes, DBH fits ~100 minutes
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/TSHE.R")
warnings()