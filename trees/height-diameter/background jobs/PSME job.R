# runtimes with 10x10 cross validation and four workers (AMD Zen 3, 4.6 GHz)
# height primary: >12 hours (mostly GAM BA+L physio)
#        gsl_nls(): ~35 minutes
# DBH primary: >12 hours (mostly GAM ABA+T physio)
#     gsl_nls(): ~21 minutes
#
# GAMs with 10x10 cross validation: approximately O(k²)
#    k     fit time, minutes
#    15      1.2
#    26      4.4
#    85     39
#   331   ~490 (8.2 hours)
#   455   ~950 (16 hours)
setwd(file.path(getwd(), "../../.."))
source("trees/height-diameter/setup.R")
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 4)

source("trees/height-diameter/PSME.R")