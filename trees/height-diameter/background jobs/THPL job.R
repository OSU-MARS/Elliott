# ~TBD minutes with 10x10 cross validation and two workers (AMD Zen 3, 4.6 GHz)
#  height fits: ~TBD fixed, ~27 minutes mixed
#  DBH fits: ~14 minutes fixed, ~30 mixed
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