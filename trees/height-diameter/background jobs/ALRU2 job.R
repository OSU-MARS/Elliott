# about three hours with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
# ~70 minutes for height fixed effect fits, ~30 minutes for most DBH fits, and then an hour hour for DBH GAM ABA+T physio and RelHt physio.
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/ALRU2.R")
warnings()