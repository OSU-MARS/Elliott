# ~26 minutes with 10x10 cross validation and one worker(AMD Zen 3, 4.6 GHz)
# ~8 minutes height fits, ~TBD minutes DBH fits
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/UMCA.R")
warnings()