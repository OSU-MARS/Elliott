# about 16 hours with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
#  height: ~70 minutes fixed effect fits, 3.7 hours for mixed effects
#  DBH: 2.1 hours fixed, ~9.0 hours mixed effects
jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/ALRU2.R")
warnings()
print(paste0("ALRU job ran for ", format(Sys.time() - jobStartTime), "."))