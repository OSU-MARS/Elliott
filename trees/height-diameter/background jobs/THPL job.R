# ~25 minutes with 10x10 cross validation and one worker(AMD Zen 3, 4.6 GHz)
# 
# splits  models  stats file size, MB
# yes     no      985
# no      no      4
#
# change working directory back to project root
# RStudio defaults background jobs' working directory to the directory the .R file is in.
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 2) # THPL n is small enough single core is likely faster but https://github.com/HenrikBengtsson/globals/issues/87

source("trees/height-diameter/THPL.R")
warnings()