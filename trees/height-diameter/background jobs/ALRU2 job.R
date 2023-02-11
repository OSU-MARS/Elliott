# four hours with 10x10 cross validation and two workers (AMD Zen 3, 4.6 GHz)
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2)

source("trees/height-diameter/ALRU2.R")