# runtimes with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
# height primary: 10.5 hours, mostly GAM BA+L physio
#                 14.5 hours GAM ABA+T RelHt physio
#        nlrob() + gsl_nls(): ~35 minutes
# DBH primary: 16.5 hours (mostly GAM ABA+T physio and RelHt physio)
#     nlrob() + gsl_nls(): ~21 minutes
#
# GAM fitting time with 10x10 cross validation: approximately O(k²)
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
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/PSME.R")
warnings()