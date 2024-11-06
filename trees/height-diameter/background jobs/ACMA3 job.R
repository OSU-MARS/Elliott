# ~TBD minutes with 10x10 cross validation and one worker(AMD Zen 3, 4.6 GHz)
#  height: ~12 minutes fixed effect (gsl_nls, nlrob, gsl_nls default weight), ~33 minutes mixed effect
#  DBH: ~18 minutes fixed, ~TBD mixed effect
jobStartTime = Sys.time()
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2) # minimum two workers, https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/ACMA3.R")
warnings()
print(paste0("ACMA job ran for ", format(Sys.time() - jobStartTime), "."))