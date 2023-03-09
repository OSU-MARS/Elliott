# about five hours with 10x10 cross validation and two workers (AMD Zen 3, 4.6 GHz)
# ~25 minutes for height fits, ~30 minutes for most DBH fits, and then four hours for DBH GAM ABA+T physio and RelHt physio.
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2)

source("trees/height-diameter/ALRU2.R")