# runtimes with 10x10 cross validation and eight workers (AMD Zen 3, 4.6 GHz)
# height: primary fixed effects 7.2 hours, nlrob() + gsl_nls() ~16 minutes, mixed: 15 hours
# DBH: 13.3 hours fixed, mixed 25 hours
#
# single worker GAM fitting time with 10x10 cross validation: approximately O(k²)
#    k     fit time, minutes
#    15      1.2
#    26      4.4
#    85     39
#   331   ~490 (8.2 hours)
#   455   ~950 (16 hours)
#
#                 height                                     DBH
#                 primary stats  nlrob() + gsl_nls() stats   primary stats  nlrob() + gsl_nls() stats
# splits  models  file size      file size                   file size      file size
# yes     no                                                 4.43 GB        5.74 GB
# no      no      1.89 MB        1.28 MB                     2.79 MB        1.23 MB
jobStartTime = Sys.time()
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 8) # increase worker count for mixed effects DBH GAMs

source("trees/height-diameter/PSME.R")
warnings()
print(paste0("PSME job ran for ", format(Sys.time() - jobStartTime), "."))