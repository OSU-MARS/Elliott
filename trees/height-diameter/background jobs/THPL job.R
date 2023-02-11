# change working directory back to project root
# RStudio defaults background jobs' working directory to the directory the .R file is in.
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 1) # THPL n is small enough single core is faster

source("trees/height-diameter/THPL.R")