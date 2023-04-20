# ~TBD minutes with 10x10 cross validation and two workers (AMD Zen 3, 4.6 GHz)
#  height: ~TBD fixed, ~20 minutes mixed
#  DBH: ~14 minutes fixed, ~30 mixed
jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/UMCA.R")
warnings()
print(paste0("UMCA job ran for ", format(Sys.time() - jobStartTime), "."))