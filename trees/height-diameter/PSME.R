# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Douglas-fir height-diameter regression form sweep
# preferred forms: Sharma-Parton BA+L, Sharma-Parton, Sharma-Zhang BA+L, Chapman-Richards BA+L
# fits without isPlantation as a factor
#psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 74, b1 = -0.014, b2 = 1.14))
#psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = 0.49, a3 = -0.053, b1 = -0.011, b2 = 1.09))
#psmeHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 35, a2 = 0.70, a3 = 0.055, a9 = 38.32, b1 = -0.010, b2 = 0.92))
#psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -9.56, a2 = 0.0003, a2p = -0.011, b1 = -0.022, b2 = 1.50, b2p = -0.31)) # b1p not significant
#psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation + a3 * sin(pi/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 66.7, a1p = -9.66, a2 = 0.0003, a2p = -0.011, a3 = -2.06, b1 = -0.022, b2 = 1.51, b2p = -0.31)) # a3p not significant, little sensitivity to sin() or cos() of slope but sin() slightly more accurate
#psmeHeightFromDiameter$chapmanRichardsTopo = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * topographicShelterIndex) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, a2 = 0, b1 = -0.022, b2 = 1.51, b2p = -0.31)) # accuracy increases through cos(), sin(), linear but differences are limited
#psmeHeightFromDiameter$curtis = fit_gsl_nls(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b2, psme2016, start = list(a1 = 1.6, b2 = 0.24))
#psmeHeightFromDiameter$hossfeld = fit_gsl_nls(TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), psme2016, start = list(a1 = 99, b1 = 187, b2 = -1.2))
#psmeHeightFromDiameter$korf = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, a1p = 376, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = nls.control(maxiter = 500)) # b1p not significant
#psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls(TotalHt ~ 1.37 + a1*DBH / (a2 + DBH), psme2016, start = list(a1 = 138, a2 = 156))
#psmeHeightFromDiameter$prodan = fit_gsl_nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), psme2016, start = list(a1 = 0.008, a2 = 1.0, a3 = 2.7))
#psmeHeightFromDiameter$power = fit_gsl_nls(TotalHt ~ 1.37 + b0*DBH^b1, psme2016, start = list(b0 = 1.54, b1 = 0.77))
#psmeHeightFromDiameter$ratkowsky = fit_gsl_nls(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), psme2016, start = list(a1 = 93, b1 = -65, b2 = 15))
#psmeHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 58.3, d = 0.609, kU = 0.0134))
#psmeHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 62.6, Hap = -26.5, d = 0.609, kU = 0.0128, kUp = 0.0114))
#psmeHeightFromDiameter$sharmaParton = fit_nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 55.4, a2 = 0.082, b1 = -0.022, b2 = -0.17, b3 = 1.04))
#psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, psme2016, start = list(a1 = 39, a2 = 0.15, b1 = -0.018, b2 = -0.16, b3 = 1.01))
#psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 55, a2 = 0.07, b1 = -0.012, b2 = 0.03, b3 = 1.1))
#psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 66, a2 = 0.056, a3 = 0.006, b1 = -0.022, b2 = -0.13, b3 = 1.04))
#psmeHeightFromDiameter$sibbesen = fit_gsl_nls(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), psme2016, start = list(a1 = 0.034, b1 = 3.1, b2 = -0.15))
#psmeHeightFromDiameter$korf = fit_gsl_nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), psme2016, start = list(a1 = 172, b1 = -9.682, b2 = -0.4574))
#psmeHeightFromDiameter$weibull = fit_gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 69, b1 = -0.0075, b2 = 1.15))
#psmeHeightFromDiameter$weibullBal = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 73.7, a2 = 0.38, a3 = -0.007, b1 = -0.008, b2 = 1.09))
#psmeHeightFromDiameter$weibullBalRelHt = fit_nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation + b2 * pmin(relativeHeight, 1.25))*DBH^(b3 + b3p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a3p = 0.217, a4 = 41.9, a4p = 76.9, b1 = -0.013, b1p = 0.007, b2 = 0, b3 = 0.94, b3p = -0.096), control = nls.control(maxiter = 500)) # step factor with fit_nlrob() 
# fits with isPlantation
#psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 68.6, a1p = -3.2, a2 = 0.024, a2p = 0.96, a3 = 0.014, a3p = -0.15, b1 = -0.018, b1p = 0.0028, b2 = 1.29, b2p = -0.10))
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, TotalHt > 0) %>% # live Douglas-firs measured for height
  mutate(dbhWeight = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016gamConstraint = c(DBH = -1.2240/0.6566, TotalHt = 1.37, standBasalAreaPerHectare = median(psme2016$standBasalAreaPerHectare), basalAreaLarger = median(psme2016$basalAreaLarger), standBasalAreaApprox = median(psme2016$standBasalAreaApprox), tallerApproxBasalArea = median(psme2016$tallerApproxBasalArea), elevation = median(psme2016physio$elevation), slope = median(psme2016physio$slope), aspect = median(psme2016physio$aspect), topographicShelterIndex = median(psme2016physio$topographicShelterIndex), relativeHeight = median(psme2016$relativeHeight)) # point constraint for mgcv::s() where response variable is ignored, zero crossing of height from DBH from fit_lm(TotalHt ~ DBH, data = psme2016 %>% filter(DBH < 6))

psme2016defaultWeight = psme2016 %>%
  mutate(dbhWeight = pmin(DBH^-1, 1))
psme2016defaultWeightPhysio = psme2016defaultWeight %>% filter(is.na(elevation) == FALSE)
#psme2016natural = psme2016 %>% filter(isPlantation == FALSE)
#psme2016plantation = psme2016 %>% filter(isPlantation)
#psme2016plantationPhysio = psme2016physio %>% filter(isPlantation)

psmeHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))) # b1p not significant
psmeHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 72.9, a1p = -11.8, a2 = 0.087, a2p = 0.84, a3 = -0.0021, a3p = -0.073, b1 = -0.016, b2 = 1.26, b2p = -0.054)) # a3 not significant, step factor with b1p
psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 74.3, a2 = 0.096, a2p = 0.92, a3 = 0, a3p = 0, a4 = -0.015, a5 = -0.101, a6 = 0.793, a7 = 1.695, a8 = 0.183, b1 = -0.018, b1p = 0.005, b2 = 1.30, b2p = -0.154)) # a4 not significant
psmeHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a4 = 54.0, a4p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37))
psmeHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31)) # a4p not significant, a5p induces overfitting
psmeHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
psmeHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = psme2016gamConstraint), data = psme2016)
psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26, pc = psme2016gamConstraint), data = psme2016, nthreads = 4)
psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6) # long fitting time since k = 455, edf = 405
psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6)
psmeHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28))
psmeHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = nls.control(maxiter = 500)) # a1p parameter evaporation, b1p not significant
psmeHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016)
psmeHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30)) # b1p not significant
psmeHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), psme2016)
psmeHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
psmeHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
psmeHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
psmeHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
psmeHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22))
psmeHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 18.5, a1p = 11.3, b1 = 0.30, b1p = -0.14, b2 = -0.019, b2p = -0.011, b3 = 0.089, b4 = 1.49, b4p = -0.44)) # b3p not significant
psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51)) # b1p not significant
psmeHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 51.9, a4 = 0.0004, a5 = 0, a6 = 0.010, a7 = 0.003, a8 = 0.004, b1 = 0.056, b2 = -0.024, b3 = -0.012, b3p = -0.23, b4 = 1.56, b4p = -0.64)) # a7, b1p, b2p not significant, fit_nlrob() step factor with a1p
psmeHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 56.1, a1p = -23.1, b1 = 0.042, b1p = 0.117, b2 = -0.0247, b2p = -0.0131, b3 = -0.0217, b3p = -0.112, b4 = 1.476, b4p = -0.456)) # b1 not significant
psmeHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 53.5, a1p = -10.1, a2 = -0.00003, a2p = 0.012, b1 = 0.070, b2 = -0.027, b2p = -0.017, b3 = -0.064, b3p = -0.096, b4 = 1.28, b4p = -0.16)) # a2, b1p, b3 not significant
psmeHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.001, a1p = 0.168, b1 = 5.78, b1p = -3.54, b2 = -0.181, b2p = 0.051), control = list(maxiter = 50))
psmeHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16))
psmeHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64.1, a2 = -0.007, a2p = 1.05, a3 = 0.032, a3p = -0.193, b1 = -0.006, b1p = -0.001, b2 = 1.231, b2p = -0.070)) # a2 not significant
psmeHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a4 = 41.9, b1 = -0.013, b1p = 0.007, b2 = 0.94, b2p = -0.096)) # a3p, a4p not significant
#confint_nlrob(psmeHeightFromDiameter$sharmaPartonBalPhysio, level = 0.99)

psmeHeightFromDiameterNls = list(chapmanRichardsNls = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$m$getPars()))
psmeHeightFromDiameterNls$chapmanRichardsBalNls = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars())
psmeHeightFromDiameterNls$chapmanRichardsBalPhysioNls = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$m$getPars())
psmeHeightFromDiameterNls$chapmanRichardsPhysioNls = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$m$getPars())
psmeHeightFromDiameterNls$curtisNls = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = psmeHeightFromDiameter$curtis$m$getPars())
psmeHeightFromDiameterNls$hossfeldNls = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$hossfeld$m$getPars())
psmeHeightFromDiameterNls$korfNls = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$korf$m$getPars(), control = nls.control(maxiter = 500))
psmeHeightFromDiameterNls$michaelisMentenNls = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = psmeHeightFromDiameter$michaelisMenten$m$getPars())
psmeHeightFromDiameterNls$prodanNls = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = psmeHeightFromDiameter$prodan$m$getPars())
psmeHeightFromDiameterNls$powerNls = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = psmeHeightFromDiameter$power$m$getPars())
psmeHeightFromDiameterNls$ratkowskyNls = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$ratkowsky$m$getPars())
psmeHeightFromDiameterNls$richardsNls = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = psmeHeightFromDiameter$richards$m$getPars())
psmeHeightFromDiameterNls$sharmaPartonNls = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$m$getPars())
psmeHeightFromDiameterNls$sharmaPartonBalNls = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars())
psmeHeightFromDiameterNls$sharmaPartonBalPhysioNls = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$m$getPars())
psmeHeightFromDiameterNls$sharmaPartonPhysioNls = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonPhysio$m$getPars())
psmeHeightFromDiameterNls$sharmaZhangNls = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$m$getPars())
psmeHeightFromDiameterNls$sharmaZhangBalNls = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars())
psmeHeightFromDiameterNls$sibbesenNls = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$sibbesen$m$getPars(), control = list(maxiter = 50))
psmeHeightFromDiameterNls$weibullNls = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$m$getPars())
psmeHeightFromDiameterNls$weibullBalNls = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$m$getPars())

psmeHeightFromDiameterNlsDwt = list(chapmanRichardsNlsFwt = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$chapmanRichards$m$getPars()))
psmeHeightFromDiameterNlsDwt$chapmanRichardsBalNlsFwt = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars())
psmeHeightFromDiameterNlsDwt$chapmanRichardsBalPhysioNlsFwt = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$m$getPars())
psmeHeightFromDiameterNlsDwt$chapmanRichardsPhysioNlsFwt = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$m$getPars())
psmeHeightFromDiameterNlsDwt$curtisNlsFwt = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$curtis$m$getPars())
psmeHeightFromDiameterNlsDwt$hossfeldNlsFwt = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$hossfeld$m$getPars())
psmeHeightFromDiameterNlsDwt$korfNlsFwt = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$korf$m$getPars(), control = nls.control(maxiter = 500))
psmeHeightFromDiameterNlsDwt$michaelisMentenNlsFwt = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016defaultWeight, start = psmeHeightFromDiameter$michaelisMenten$m$getPars())
psmeHeightFromDiameterNlsDwt$prodanNlsFwt = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$prodan$m$getPars())
psmeHeightFromDiameterNlsDwt$powerNlsFwt = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$power$m$getPars())
psmeHeightFromDiameterNlsDwt$ratkowskyNlsFwt = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$ratkowsky$m$getPars())
psmeHeightFromDiameterNlsDwt$richardsNlsFwt = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$richards$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaPartonNlsFwt = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaParton$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaPartonBalNlsFwt = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaPartonBalPhysioNlsFwt = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaPartonPhysioNlsFwt = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$sharmaPartonPhysio$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaZhangNlsFwt = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaZhang$m$getPars())
psmeHeightFromDiameterNlsDwt$sharmaZhangBalNlsFwt = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars())
psmeHeightFromDiameterNlsDwt$sibbesenNlsFwt = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$sibbesen$m$getPars(), control = nls.control(maxiter = 100))
psmeHeightFromDiameterNlsDwt$weibullNlsFwt = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$weibull$m$getPars())
psmeHeightFromDiameterNlsDwt$weibullBalNlsFwt = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$weibullBal$m$getPars())

psmeHeightFromDiameterResults = bind_rows(bind_rows(lapply(psmeHeightFromDiameter, as_row)),
                                          bind_rows(lapply(psmeHeightFromDiameterNls, as_row, fitSet = "gsl_nls")),
                                          bind_rows(lapply(psmeHeightFromDiameterNlsDwt, as_row, fitSet = "gsl_nls", fixedWeight = -1))) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(psmeHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$curtis), color = "Curtis", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$korf), color = "Korf", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$linear), color = "linear", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$parabolic), color = "parabolic", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$power), color = "power", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$prodan), color = "Prodan", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$richards), color = "unified Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sibbesen), color = "Sibbesen", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibull), color = "Weibull", group = psme2016$isPlantation)) +
    annotate("text", x = 0, y = 85, label = "Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
  
  ggplot() +
    geom_point(aes(x = psme2016physio$elevation, y = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio)), alpha = 0.15, color = "grey25", shape = 16) +
    geom_smooth(aes(x = psme2016physio$elevation, y = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio)), alpha = 0.15, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    labs(x = "physiographic", y = "DBH error, cm")
  
  psmeHeightFromDiameter$Efficiency = psmeHeightFromDiameterResults %>% filter(fitting %in% c("nlrob", "gsl_nls"), (fitting == "nlrob") | (is.na(fixedWeight) == FALSE), str_detect(name, "RelHt") == FALSE) %>% 
    select(fitting, name, mae, mape, rmse, rmspe, aicN, nse, pearson) %>% 
    pivot_wider(names_from = fitting, values_from = c(mae, mape, rmse, rmspe, aicN, nse, pearson)) %>%
    mutate(deltaMae = mae_nlrob - mae_nls, deltaMape = mape_nlrob - mape_nls,
           deltaRmse = rmse_nlrob - rmse_nls, deltaRmspe = rmspe_nlrob - rmspe_nls,
           deltaAicN = aicN_nlrob - aicN_nls, deltaNse = nse_nlrob - nse_nls, deltaPearson = pearson_nlrob - pearson_nls)
  
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaMae), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, m", y = "number of regression forms") +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaRmse), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, m", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaAicN), bins = 30) +
    coord_cartesian(xlim = c(-3.75, -3.25), ylim = c(0, 7)) +
    labs(x = "normalized ΔAIC", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaNse), bins = 30) +
    coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(0, 7)) +
    labs(x = "change in model efficiency", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaPearson), bins = 30) +
    coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(0, 7)) +
    labs(x = "change in Pearson's R", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaMape), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, %", y = "number of regression forms") +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaRmspe), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, %", y = NULL) +
  plot_layout(ncol = 5, nrow = 2)
  
  ggplot(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, error = -residuals(psmeHeightFromDiameter$sharmaParton))) +
    geom_boxplot(aes(x = dbhClass, y = error, group = dbhClass)) +
    coord_cartesian(ylim = c(-40, 40)) +
    labs(x = NULL, y = "fit_nlrob() height residual, m") +
  ggplot(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, error = -residuals(psmeHeightFromDiameter$sharmaPartonNls))) +
    geom_boxplot(aes(x = dbhClass, y = error, group = dbhClass)) +
    coord_cartesian(ylim = c(-40, 40)) +
    labs(x = "DBH, cm", y = "nls() height residual, m") +
  plot_layout(ncol = 1, nrow = 2)
  
  ggplot(bind_rows(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, 
                          chapmanRichards = -residuals(psmeHeightFromDiameter$chapmanRichards),
                          chapmanRichardsBal = -residuals(psmeHeightFromDiameter$chapmanRichardsBal),
                          curtis = -residuals(psmeHeightFromDiameter$curtis),
                          gam = -residuals(psmeHeightFromDiameter$gam),
                          gamBal = -residuals(psmeHeightFromDiameter$gamBal),
                          hossfeld = -residuals(psmeHeightFromDiameter$hossfeld),
                          korf = -residuals(psmeHeightFromDiameter$korf),
                          linear = -residuals(psmeHeightFromDiameter$linear),
                          michaelisMenten = -residuals(psmeHeightFromDiameter$michaelisMenten),
                          parabolic = -residuals(psmeHeightFromDiameter$parabolic),
                          power = -residuals(psmeHeightFromDiameter$power),
                          prodan = -residuals(psmeHeightFromDiameter$prodan),
                          ratkowsky = -residuals(psmeHeightFromDiameter$ratkowsky),
                          richards = -residuals(psmeHeightFromDiameter$richards),
                          sharmaParton = -residuals(psmeHeightFromDiameter$sharmaParton), 
                          sharmaPartonBal = -residuals(psmeHeightFromDiameter$sharmaPartonBal),
                          sharmaZhang = -residuals(psmeHeightFromDiameter$sharmaZhang),
                          sharmaZhangBal = -residuals(psmeHeightFromDiameter$sharmaZhangBal),
                          sibbesen = -residuals(psmeHeightFromDiameter$sibbesen),
                          weibull = -residuals(psmeHeightFromDiameter$weibull),
                          weibullBal = -residuals(psmeHeightFromDiameter$weibullBal),
                          sharmaPartonNls = -residuals(psmeHeightFromDiameter$sharmaPartonNls)),
                   tibble(dbhClass = 2.5 * floor(psme2016physio$DBH / 2.5) + 0.5 * 2.5,
                          chapmanRichardsBalPhysio = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio),
                          chapmanRichardsPhysio = -residuals(psmeHeightFromDiameter$chapmanRichardsPhysio),
                          gamBalPhysio = -residuals(psmeHeightFromDiameter$gamBalPhysio),
                          gamPhysio = -residuals(psmeHeightFromDiameter$gamPhysio),
                          sharmaPartonBalPhysio = -residuals(psmeHeightFromDiameter$sharmaPartonBalPhysio),
                          sharmaPartonPhysio = -residuals(psmeHeightFromDiameter$sharmaPartonPhysio))) %>%
           pivot_longer(cols = !dbhClass, names_to = "name", values_to = "residual") %>%
           mutate(name = factor(name, levels = c("chapmanRichards", "chapmanRichardsBal", "chapmanRichardsBalPhysio", "chapmanRichardsPhysio", "curtis", "gam", "gamBal", "gamBalPhysio", "gamPhysio", "hossfeld", "korf", "linear", "michaelisMenten", "parabolic", "power", "prodan", "ratkowsky", "richards", "sharmaParton", "sharmaPartonBal", "sharmaPartonBalPhysio", "sharmaPartonPhysio", "sharmaZhang", "sharmaZhangBal", "sibbesen", "weibull", "weibullBal", "sharmaPartonNls"),
                                      labels = c("Chapman-Richards", "Chapman-Richards BA+L", "Chapman-Richards BA+L physio", "Chapman-Richards physio", "Curtis", "GAM", "GAM BAL", "GAM BAL physio", "GAM physio", "Hossfeld", "Korf", "linear", "Michaelis-Menten", "parabolic", "power", "Prodan", "Ratkowsky", "Richards", "Sharma-Parton", "Sharma-Parton BA+L", "Sharma-Parton BA+L physio", "Sharma-Parton physio", "Sharma-Zhang", "Sharma-Zhang BA+L", "Sibbesen", "Weibull", "Weibull BA+L", "Sharma-Parton nls()"))) %>%
           group_by(name, dbhClass) %>%
           summarize(n = n(), bias = mean(residual, na.rm = TRUE), 
                     fitting = if_else(name %in% c("GAM", "GAM BAL", "GAM BAL physio", "GAM"), "GAM", if_else(name %in% c("linear", "parabolic"), "linear", "nonlinear")),
                     .groups = "drop") %>%
           filter(n > 10, is.na(bias) == FALSE)) +
    geom_line(aes(x = dbhClass, y = bias, alpha = fitting, color = fitting, group = name)) +
    coord_cartesian(ylim = c(-5, 5)) +
    labs(x = "DBH, cm", y = "height bias, m", alpha = NULL, color = NULL) +
    scale_alpha_manual(breaks = c("GAM", "nonlinear", "linear"), values = c(1, 0.2, 0.3)) +
    scale_color_discrete(breaks = c("GAM", "nonlinear", "linear")) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))
}


## Douglas-fir height-diameter GNLS regressions
#psmeHeightFromDiameterGnls$chapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterGnls$chapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 0.001, tolerance = 0.01, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05, iterations at default msTol and tolerance
#psmeHeightFromDiameterGnls$sharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 500, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterGnls$sharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#psmeHeightFromDiameterGnls$sharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, msTol = 1E-5, tolerance = 1E-4, returnObject = FALSE)) # step having at nlsTol = 0.1
#psmeHeightFromDiameterGnls$sharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, nlsMaxIter = 50, msTol = 1E-7, tolerance = 1E-6, msVerbose = FALSE, returnObject = FALSE)) # step halving factor at nlsTol = 1 with plot correlation, step halving with nlsTol = 0.07
#psmeHeightFromDiameterGnls$weibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterGnls$weibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 5, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.2
#save(psmeHeightFromDiameterGnls, file = "trees/height-diameter/HtDia PSME GNLS.rdata")

load("trees/height-diameter/HtDia PSME GNLS.rdata")
psmeHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", psmeHeightFromDiameterGnls$chapmanRichards, psme2016)
psmeHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", psmeHeightFromDiameterGnls$chapmanRichardsBal, psme2016)
psmeHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", psmeHeightFromDiameterGnls$sharmaParton, psme2016)
psmeHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", psmeHeightFromDiameterGnls$sharmaPartonBal, psme2016)
psmeHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", psmeHeightFromDiameterGnls$sharmaZhang, psme2016)
psmeHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", psmeHeightFromDiameterGnls$sharmaZhangBal, psme2016)
psmeHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", psmeHeightFromDiameterGnls$weibull, psme2016)
psmeHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", psmeHeightFromDiameterGnls$weibullBal, psme2016)

psmeHeightFromDiameterResultsGnls = bind_rows(lapply(psmeHeightFromDiameterGnls, as_row)) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  psmeHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(psmeHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(psmeHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(psmeHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  
  ggplot() +
    geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsBal), color = "Chapman-Richards BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsPhysio), color = "Chapman-Richards physio", group = psme2016physio$isPlantation), alpha = 0.5) +
    geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsTopo), color = "Chapman-Richards topo", group = psme2016physio$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaPartonBal), color = "Sharma-Parton BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhangBal), color = "Sharma-Zhang BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibullBal), color = "Weibull BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterGnls$sharmaParton), color = "Sharma-Parton GNLS"), alpha = 0.35) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterGnls$weibull), color = "Weibull GNLS"), alpha = 0.35) + 
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibull), color = "Weibull", group = psme2016$isPlantation)) +
    annotate("text", x = 0, y = 85, label = "a) Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c("ElliottChapmanRichards", "ElliottWeibull", "ElliottChapmanRichardsBal", "ElliottWeibullBalNatural", "ElliottWeibullBalPlantation", "TemesgenWeibull"), labels = c("Chapman-Richards", "Weibull", "Chapman-Richards with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_color_manual(values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot(psme2016) +
    geom_point(aes(x = Elev_Mean, y = -residuals(psmeHeightFromDiameter$weibullBal)), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = Elev_Mean, y = -residuals(psmeHeightFromDiameter$weibullBal)), color = "red", formula = y ~ s(x, k = 10), method = "gam")
  
  ggplot(psme2016) +
    geom_point(aes(x = AspectCos, y = -residuals(psmeHeightFromDiameter$weibullBal)), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = AspectCos, y = -residuals(psmeHeightFromDiameter$weibullBal)), color = "red", formula = y ~ s(x, k = 10), method = "gam")
}


## Douglas-fir diameter-height regressions
# Chapman-Richards, Gomperz, Prodan, Ratkowsky, and Korf run opposite the desired curvature
# Huang has DBH as ln of height (https://doi.org/10.1093/forestry/cpz002, Eq. 1)
# Hossfeld is sensitive to going negative, extrapolating poorly
# Prodan sometimes produces negative DBH
# Sharma-Parton and Sharma-Zhang have wrong curvature, simply modified forms fail to start or fail to step
# Sibbesen and Korf are near identical
#psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999999)), psme2016, start = list(a1 = -65, b1 = 1/85, b2 = 0.75))
#psmeDiameterFromHeight$chapmanFormBalRelHt = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 109, a1p = -40.0, a2 = -1.42, a3 = 0.490, a3p = 0.120, a4 = -45.1, a4p = 36.5, b1 = 0.041, b2 = 0.763, b2p = -0.0092)) # a2p not significant
#psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016physio, start = list(a1 = -20.8, a1p = 5.59, a2 = 0.00218, a3 = -3.39, a4 = 0.209, a5 = 0.109, a6 = 0.101, b1 = 0.0176, b1p = 0.00776, b2 = 0.329)) # a5 not significant
#psmeDiameterFromHeight$michaelisMentenForm = fit_gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), psme2016, start = list(a1 = 108, a2 = 46.3, b1 = 0.77), weights = TotalHt^-2)
#psmeDiameterFromHeight$naslund = fit_gsl_nls(DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 3.7, a2 = -0.09), weights = TotalHt^-2) # diverges if a2 start value is less than ~-0.1
#psmeDiameterFromHeight$power = fit_gsl_nls(DBH ~ a1*(TotalHt - 1.37)^b1, psme2016, start = list(a1 = 1.26, b1 = 1.08))
#psmeDiameterFromHeight$prodan = fit_gsl_nls(DBH ~ (TotalHt - 1.37)^2 / (a0 + a1 * (TotalHt - 1.37) + a2 * (TotalHt - 1.37)^2), psme2016, start = list(a0 = -0.7, a1 = 0.8, a2 = -0.004)) # AIC 
#psmeDiameterFromHeight$prodan = get_dbh_error(psmeDiameterFromHeight$prodan, psme2016)
#psmeDiameterFromHeight$sharmaParton = fit_gsl_nls(DBH ~ a1*topHeight^a2*exp(b1*(tph/standBasalAreaPerHectare)^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 1, a2 = 1, b1 = -0.01, b2 = 1, b3 = 1))
#psmeDiameterFromHeight$sharmaParton = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 57.1, a1p = -25.5, a2 = -0.241, a2p = 0.060, b1 = 0.0261, b2 = -0.0535, b2p = 0.0906, b3 = 0.689, b3p = 0.0829))
#psmeDiameterFromHeight$sharmaParton = fit_gsl_nls(DBH ~ a1*topHeight^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, psme2016, start = list(a1 = 39.3, a2 = 0.127, b1 = 0.0180, b2 = -0.010, b3 = 0.824))
#psmeDiameterFromHeight$sharmaZhang = fit_gsl_nls(DBH ~ a1*standBasalAreaPerHectare^a2*exp(b1*tph^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.0005, b2 = 0.5, b3 = 1))
#psmeDiameterFromHeight$wykoff = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.081, a1p = -0.390, b1 = 1.576, b2 = 0.263, b2p = 0.0237))
#psmeDiameterFromHeight$wykoff = fit_gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 15, b1 = 0.1, b2 = 0.35)) # NaN or inf 
#psmeDiameterFromHeight$wykoff = fit_gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 36.6, b1 = 0.059, b2 = 0.77))
#psmeDiameterFromHeight$wykoff = fit_gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 63.0, b1 = 0.018, b2 = 0.87))
#psmeDiameterFromHeight$wykoff = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 17.8, a1p = -6.45, b1 = 0.216, b2 = 0.541, b2p = 0.048))
#psmeDiameterFromHeight$wykoffBal = fit_gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 121.7, a2 = -2.17, a3 = 1.04, b1 = 0.017, b2 = 0.85), control = list(maxiter = 500))
#psmeDiameterFromHeight$wykoffBalRelHt = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.340, a1p = -0.175, a2 = -0.0046, a2p = 0.00162, a3 = 0.00142, a3p = -0.00058, a4 = -0.132, a4p = 0.119, b1 = 3.040, b2 = 0.173, b2p = 0.0129)) # >500 iterations with b1p
#psmeDiameterFromHeight$wykoffRelHt = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.065, a1p = -0.386, a2 = 0.0071, b1 = 1.59, b2 = 0.262, b2p = 0.0234))
#psmeDiameterFromHeight$weibull = fit_gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -204, b1 = 0.011, b2 = 0.87))
#psmeDiameterFromHeight$weibull = fit_gsl_nls(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -374, a1p = 163, b1 = 0.009, b1p = 0.003, b2 = 0.82)) # b2p not significant
#cor(cbind(dbh = psme2016$DBH, bal = psme2016$basalAreaLarger, height = psme2016$TotalHt, relHt = psme2016$relativeHeight, topHt = psme2016$topHeight, deltaHt = psme2016$TotalHt - psme2016$topHeight, treesPerPlot = psme2016$treesPerPlot, tph = psme2016$tph, tallerTph = psme2016$tallerTph * psme2016$topHeight, tallerQuasiBA = psme2016$tallerQuasiBA))
psmeDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780)))
psmeDiameterFromHeight$chapmanFormAat = fit_nlrob("Chapman-Richards form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.4, a1p = -44.6, a2 = 0.0058, a3 = -0.0544, b1 = 0.00166, b1p = 0.017, b2 = 0.788, b2p = -0.056)) # a2, a2p, a3p not significant
psmeDiameterFromHeight$chapmanFormBal = fit_nlrob("Chapman-Richards form BA+L", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 135, a1p = -37.5, a2 = -1.2, b1 = 0.010, b1p = 0.002, b2 = 0.756, b2p = 0.064), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a3 * BA
psmeDiameterFromHeight$chapmanFormBalRelHt = fit_nlrob("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 139, a1p = -38.5, a2 = -1.2, a9 = -5.5, a9p = 0.22, b1 = 0.012, b1p = 0.003, b2 = 0.796, b2p = 0.066), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a9 * BA
psmeDiameterFromHeight$chapmanFormRelHt = fit_nlrob("Chapman-Richards form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 21.1, a1p = -6.56, a9 = 0.278, b1 = 0.170, b2 = 0.580, b2p = 0.0395)) # a4p not significant
psmeDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781)) # b2p not significant
psmeDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786)) # a2p, a3p, b2p not significant
psmeDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -13.8, a1p = -0.65, a5 = -3.58, a8 = 0.091, b1 = 0.019, b1p = 0.0062, b2 = 0.30), maxit = 20, control = nls.control(maxiter = 50)) # a4, a6, a7 not significant
psmeDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -1.67, a9 = -9.26, b1 = 0.020, b1p = 0.004, b2 = 0.004, b2p = 0.089), maxit = 20) # a1p not significant
psmeDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = psme2016gamConstraint), data = psme2016)
psmeDiameterFromHeight$gamAat = fit_gam("REML GAM AA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28, pc = psme2016gamConstraint), data = psme2016, nthreads = 4)
psmeDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM AA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 496, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6) # slow since minimum (k = 496, edf = 397), AIC 151550: 151552 without AAT, 151700 without ABA, 151783 without elevation, 151533 without slope, 151561 without sin(aspect), 151562 without cos(aspect), 151650 without topgraphic shelter
psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6)
psmeDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016)
psmeDiameterFromHeight$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08))
psmeDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018))
psmeDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016)
psmeDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108))
psmeDiameterFromHeight$powerAat = fit_nlrob("power AA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053)) # a3 not significant
psmeDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102)) # a7, a8 not significant
psmeDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068)) # a1p, a4p not significant
psmeDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096))
#psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 18.6)) # converges from red alder values but fails to reconverge (singular gradient or NaN-inf with nls())
psmeDiameterFromHeight$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 3.95, b1 = 0.681, b1p = -0.139, b2 = 0.097, b3 = -0.130, b3p = 0.157, b4 = 0.125, b4p = 0.0678), control = list(maxiter = 200)) # a1p NaN-inf (not significant?), singular gradient with all relative height forms attempted
psmeDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
psmeDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, a3p = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190)) # a2, a2p not significant
psmeDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = list(a1 = 3.812, a1p = -0.925, a4 = 0.00022, a6 = -0.0495, a7 = -0.0151, a8 = -0.00630, b1 = 0.520, b2 = 0.111, b2p = 0.0173)) # a5 not significant
psmeDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.90, a1p = -0.95, a9 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016)) # a4p not significant
psmeDiameterFromHeight$weibull = fit_nlrob("Weibull", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81)) # b2p not significant
#confint_nlrob(psmeDiameterFromHeight$sibbesenFormRelHt, level = 0.99)

psmeDiameterFromHeightNls = list(chapmanForm = fit_gsl_nls("Chapman-Richards form", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanForm$m$getPars()))
psmeDiameterFromHeightNls$chapmanFormAat = fit_gsl_nls("Chapman-Richards form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanFormAat$m$getPars())
psmeDiameterFromHeightNls$chapmanFormRelHt = fit_gsl_nls("Chapman-Richards form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = psmeDiameterFromHeight$chapmanFormRelHt$m$getPars())
psmeDiameterFromHeightNls$chapmanRichards = fit_gsl_nls("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichards$m$getPars())
psmeDiameterFromHeightNls$chapmanRichardsAat = fit_gsl_nls("Chapman-Richards AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsAat$m$getPars())
psmeDiameterFromHeightNls$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = psmeDiameterFromHeight$chapmanRichardsPhysio$m$getPars(), control = nls.control(maxiter = 75))
psmeDiameterFromHeightNls$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsRelHt$m$getPars())
psmeDiameterFromHeightNls$michaelisMentenForm = fit_gsl_nls("Michaelis-Menten form", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = psmeDiameterFromHeight$michaelisMentenForm$m$getPars())
psmeDiameterFromHeightNls$naslund = fit_gsl_nls("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$naslund$m$getPars())
psmeDiameterFromHeightNls$power = fit_gsl_nls("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$power$m$getPars())
psmeDiameterFromHeightNls$powerAat = fit_gsl_nls("power AA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$powerAat$m$getPars())
psmeDiameterFromHeightNls$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = psmeDiameterFromHeight$powerPhysio$m$getPars())
psmeDiameterFromHeightNls$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = psmeDiameterFromHeight$powerRelHt$m$getPars())
psmeDiameterFromHeightNls$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$ruark$m$getPars())
#psmeDiameterFromHeightNls$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = psmeDiameterFromHeight$schnute$m$getPars())
psmeDiameterFromHeightNls$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = psmeDiameterFromHeight$sharmaParton$m$getPars(), control = list(maxiter = 200))
psmeDiameterFromHeightNls$sibbesenForm = fit_gsl_nls("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenForm$m$getPars())
psmeDiameterFromHeightNls$sibbesenFormAat = fit_gsl_nls("Sibbesen form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormAat$m$getPars())
psmeDiameterFromHeightNls$sibbesenFormPhysio = fit_gsl_nls("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = psmeDiameterFromHeight$sibbesenFormPhysio$m$getPars())
psmeDiameterFromHeightNls$sibbesenFormRelHt = fit_gsl_nls("Sibbesen form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormRelHt$m$getPars())
psmeDiameterFromHeightNls$weibull = fit_gsl_nls("Weibull", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = psmeDiameterFromHeight$weibull$m$getPars())

psmeDiameterFromHeightResults = bind_rows(bind_rows(lapply(psmeDiameterFromHeight, as_row)),
                                          as_row(name = "Schnute"),
                                          bind_rows(lapply(psmeDiameterFromHeightNls, as_row, fitSet = "gsl_nls"))) %>% 
  mutate(responseVariable = "DBH", species = "PSME", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(psmeDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  psmeParameters = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_coefficients)),
                                       bind_rows(lapply(psmeHeightFromDiameterGnls, get_coefficients)),
                                       bind_rows(lapply(psmeHeightFromDiameterNls, get_coefficients)),
                                       bind_rows(lapply(psmeHeightFromDiameterNlsDwt, get_coefficients))) %>%
                               mutate(responseVariable = "height"),
                             bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_coefficients)),
                                       #get_coefficients(psmeDiameterFromHeight$schnute)
                                       bind_rows(lapply(psmeDiameterFromHeightNls, get_coefficients)) %>%
                               mutate(responseVariable = "DBH"))) %>%
    mutate(species = "PSME") %>%
    relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman form", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman form AA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$linear), y = TotalHt, color = "linear", group = isPlantation)) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$parabolic), y = TotalHt, color = "parabolic", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$schnute), y = TotalHt, color = "Schnute", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -65 * log(1 - pmin((1/85*(TotalHt - 1.37))^0.7, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = -155 * log(1 - pmin((0.00839*(TotalHt - 1.37))^0.970, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = 0.05*topHeight*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 3.5*sqrt(TotalHt - 1.37) / (1 - 0.1*sqrt(TotalHt - 1.37)), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = (-204*log(1 - pmin(0.011 * (TotalHt - 1.37), 0.9999)))^0.869, y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 90, label = "Douglas-fir, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot() +
    #geom_segment(aes(x = 0, y = 0, xend = 80, yend = 50), color = "grey70", linetype = "longdash") +
    #geom_segment(aes(x = 0, y = 0, xend = 80, yend = -50), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = psme2016$tallerTph, y = -residuals(psmeDiameterFromHeight$chapmanRH)), alpha = 0.05, shape = 16) +
    geom_smooth(aes(x = psme2016$tallerTph, y = -residuals(psmeDiameterFromHeight$chapmanRH)), color = "red", fill = "red") +
    labs(x = "height, m", y = "residual DBH error, cm")
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_line(aes(x = DBH, y = 1.37 + 300 * (1 * 0.1 * pmin(relativeHeight, 2)) * exp(-0.5*DBH^-1)), color = "red")
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = heightDiameterRatio, color = Species), alpha = 0.2, na.rm = TRUE, shape = 16) +
    labs(x = "DBH, cm", y = "height-diameter ratio", color = NULL) +
    scale_color_manual(breaks = c("DF", "WH"), values = c("green3", "blue2")) +
    scale_y_continuous(breaks = seq(0, 300, by = 40)) +
    theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.98))
  
  psmeDiameterFromHeight$Efficiency = psmeDiameterFromHeightResults %>% filter(fitting %in% c("nlrob", "nls"), str_detect(name, "BA\\+L") == FALSE) %>% 
    select(fitting, name, mae, mape, rmse, rmspe, aicN, nse, pearson) %>% 
    pivot_wider(names_from = fitting, values_from = c(mae, mape, rmse, rmspe, aicN, nse, pearson)) %>%
    mutate(deltaMae = mae_nlrob - mae_nls, deltaRmse = rmse_nlrob - rmse_nls, deltaAicN = aicN_nlrob - aicN_nls, deltaNse = nse_nlrob - nse_nls, deltaPearson = pearson_nlrob - pearson_nls)
  
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaMae), bins = 30) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 8)) +
    labs(x = "ΔMAE, m", y = "number of regression forms") +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaRmse), bins = 30) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 8)) +
    labs(x = "ΔRMSE, cm", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaAicN), bins = 30) +
    coord_cartesian(xlim = c(-6.5, -6), ylim = c(0, 8)) +
    labs(x = "normalized ΔAIC", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaNse), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 8)) +
    labs(x = "change in model efficiency", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaPearson), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 8)) +
    labs(x = "change in Pearson's R", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaMape), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, %", y = "number of regression forms") +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaRmspe), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, %", y = NULL) +
  plot_layout(ncol = 5, nrow = 2)
}


## Q-Q plots
if (includeInvestigatory)
{
  ggplot() + # symmetric t fits less well than skewed t, as expected
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'t(df = 7, '*alpha*' = 2.25)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
  ggplot() + # slow! and fits poorly
  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'pearsonIV(m = 1.72, '*nu*' = 0.656)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
}


## basal area from height
# essentially no difference between fit_gsl_nls() and fit_nlrob() fits
# Chapman-Richards has the wrong curvature
if (includeInvestigatory)
{
  #psmeBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1, psme2016, start = list(a1 = 1, a1p = 0, b1 = 0.00009, b1p = 0, b2 = 2.2, b2p = 0), weights = pmin(1/basalArea, 1E4)) # a1p not significant
  #psmeBasalAreaFromHeightPower = fit_gsl_nls(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 0.25 * pi * 0.01^2, a1p = 0, b1 = 2.4, b1p = 0), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
  psmeBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 0.689, a1p = -0.413, b1 = 0.0003, b1p = 0.0005, b2 = 1.91, b2p = -0.10), weights = pmin(1/basalArea, 1E4)) # a1p not significant
  psmeBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = 0.00005, b1 = 2.41, b1p = -0.248), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
  #confint_nlrob(psmeBasalAreaFromHeightKorf, level = 0.99, weights = 1/psme2016$basalArea)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(psmeBasalAreaFromHeightKorf), 100^2 * mean(-residuals(psmeBasalAreaFromHeightKorf)), mean(abs(residuals(psmeBasalAreaFromHeightKorf))), 1 - sum(residuals(psmeBasalAreaFromHeightKorf)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2)),
          "power", AIC(psmeBasalAreaFromHeightPower), 100^2 * mean(-residuals(psmeBasalAreaFromHeightPower)), mean(abs(residuals(psmeBasalAreaFromHeightPower))), 1 - sum(residuals(psmeBasalAreaFromHeightPower)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(psme2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "measured or imputed Douglas-fir height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
  
  ggplot(psme2016) +
    geom_point(aes(x = basalArea, y = -residuals(psmeBasalAreaFromHeightKorf), color = "Korf"), alpha = 0.1, shape = 16) +
    labs(x = "basal area, m²", y = "residual, m", color = NULL) +
    theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
  ggplot(psme2016) +
    geom_point(aes(x = basalArea, y = -residuals(psmeBasalAreaFromHeightPower), color = "power"), alpha = 0.1, shape = 16) +
    labs(x = "basal area, m²", y = "residual, m", color = NULL) +
    theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
  plot_layout(nrow = 1, ncol = 2)
    
  ggplot(psme2016) +
    geom_point(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.1, color = "grey25", shape = 16) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.10, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.95)) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.25, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.80)) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.35, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.50)) +
    geom_quantile(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), color = "darkviolet", formula = y ~ x) +
    labs(x = "Douglas-fir height, m", y = "basal area residual, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## exploratory plots
if (includeInvestigatory)
{
  ggplot() +
    geom_point(aes(x = psme2016natural$DBH, y = psme2016natural$TotalHt), alpha = 0.10, color = "green4", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = psme2016natural$DBH, y = psme2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "Douglas-fir naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = psme2016plantation$DBH, y = psme2016plantation$TotalHt), alpha = 0.10, color = "grey20", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = psme2016plantation$DBH, y = psme2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "Douglas-fir plantation height, m")
  #ggsave("Presentation/Douglas-fir height-diameter natural-plantation.png", width = 12.5, height = 11, units = "cm")
}