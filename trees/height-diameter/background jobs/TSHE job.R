# ~4.5 hours with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
#  height: ~80 minutes fixed effects, ~60 minutes mixed 
#  DBH: ~55 minutes fixed, ~1.2 hours mixed
jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/TSHE.R")
warnings()
print(paste0("TSHE job ran for ", format(Sys.time() - jobStartTime), "."))