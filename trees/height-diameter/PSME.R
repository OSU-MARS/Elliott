# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Douglas-fir height-diameter regression form sweep
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, is.na(TotalHt) == FALSE) %>% # live Douglas-firs measured for height
  mutate(dbhWeight = pmin(TreeCount/(1.62*DBH^(0.73 - 0.041*isPlantation)), 5*TreeCount), # residuals.R::varianceForHeight$psme
         heightWeight = pmin(TreeCount/((5.661 - 5.655*isPlantation)*(TotalHt - 1.37)^(1.09 + 1.79*isPlantation)), 5*TreeCount)) # residuals.R:varianceForDbh$psme
psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016gamConstraint = c(DBH = -1.2240/0.6566, TotalHt = 1.37, standBasalAreaPerHectare = median(psme2016$standBasalAreaPerHectare), basalAreaLarger = median(psme2016$basalAreaLarger), standBasalAreaApprox = median(psme2016$standBasalAreaApprox), tallerApproxBasalArea = median(psme2016$tallerApproxBasalArea), elevation = median(psme2016physio$elevation), slope = median(psme2016physio$slope), aspect = median(psme2016physio$aspect), topographicShelterIndex = median(psme2016physio$topographicShelterIndex), relativeHeight = median(psme2016$relativeHeight), relativeDiameter = median(psme2016$relativeDiameter)) # point constraint for mgcv::s() where response variable is ignored, zero crossing of height from DBH from fit_lm(TotalHt ~ DBH, data = psme2016 %>% filter(DBH < 6))

psme2016defaultWeight = psme2016 %>% mutate(dbhWeight = pmin(TreeCount/DBH, 5*TreeCount),
                                            heightWeight = pmin(TreeCount/TotalHt, 5*TreeCount))
psme2016defaultWeightPhysio = psme2016defaultWeight %>% filter(is.na(elevation) == FALSE)

psmeOptions = tibble(fitHeightPrimary = TRUE,
                     fitHeightNlrobAndFixedWeight = fitHeightPrimary,
                     fitHeightGnls = FALSE,
                     fitHeightMixed = fitHeightPrimary,
                     fitDbhPrimary = TRUE,
                     fitDbhGslNlsAndGams = fitDbhPrimary,
                     fitDbhNlrobAndFixedWeight = fitDbhPrimary,
                     fitDbhMixed = fitDbhPrimary,
                     fitAbatRelHtPhysioGam = TRUE,
                     fitPhysioGams = TRUE)

if (psmeOptions$fitHeightPrimary)
{
  # linear regressions
  psmeHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016))
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), psme2016)
  # nonlinear regressions
  psmeHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))
  psmeHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65, a1p = -9, a2 = 0.07, a2p = 0.6, a3 = 0.07, a3p = -0.05, b1 = -0.016, b2 = 1.25, b2p = -0.07)) # a2 debatably significant
  psmeHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64, a1p = -23, a2 = 0, a2p = 0.48, a3 = 0.07, a3p = 0.09, a10 = -0.77, a10p = 2.1, b1 = -0.020, b2 = 1.44, b2p = -0.28)) # a2 not significant
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 63, a2 = 0.03, a2p = 0.67, a3 = 0.084, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, b1 = -0.021, b1p = 0.009, b2 = 1.5, b2p = -0.4)) # a2, a3p not significant
  psmeHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 61, a2 = 0, a2p = 0.72, a3 = 0.11, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, a10 = -0.4, a10p = 1.2, b1 = -0.023, b1p = 0.010, b2 = 1.54, b2p = -0.49)) # a2, a10 not significant
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31))
  psmeHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64, a1p = -22, a2 = 0, a2p = 0.5, a3 = -0.07, a3p = 0.08, a10 = -0.77, a10p = 2.1, b1 = -0.020, b2 = 1.43, b2p = -0.28)) # a2 not significant
  psmeHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 72, a1p = -18, a4 = -0.0045, a5 = -8.6, a6 = 0.7, a7 = 0.8, a8 = 0.23, a10 = -1.2, a10p = 1.6, b1 = -0.022, b2 = 1.5, b2p = -0.36))
  psmeHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
  psmeHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28))
  psmeHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = gsl_nls_control(maxiter = 500))
  psmeHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30))
  psmeHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
  psmeHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
  psmeHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
  psmeHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
  psmeHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 15, b1 = 0.35, b1p = -0.047, b2 = -0.023, b2p = -0.009, b3 = 0.002, b3p = -0.09, b4 = 1.52, b4p = -0.42))
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3 * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 6, a1p = 13, b1 = 0.6, b1p = -0.33, b2 = -0.025, b3 = 0.0, b4 = 1.53, b4p = -0.15)) # b2p not significant
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 20, a1p = -2.8, a4 = -0.0016, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.07, b1 = 0.30, b2 = -0.035, b2p = -0.007, b3 = -0.003, b3p = -0.07, b4 = 1.57, b4p = -0.51)) # b3 not significant
  psmeHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a10 = -0.35, a10p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40)) # b3 not significant
  psmeHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3 * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 4.2, a1p = 14.5, a10 = 0.16, a10p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66)) # a2, a2p not significant
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65))
  psmeHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 13.3, a10 = -0.2, a10p = 0.73, b1 = 0.40, b1p = -0.12, b2 = -0.020, b2p = -0.0233, b3 = 0.032, b3p = -0.20, b4 = 1.5, b4p = -0.39))
  psmeHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 15.5, a4 = -0.0016, a5 = -0.028, a7 = 0.17, a8 = 0.06, a10 = 0.48, a10p = -0.3, b1 = 0.33, b2 = -0.027, b3 = -0.043, b3p = -0.074, b4 = 1.52, b4p = -0.55)) # a6-a7 not reliably mutually significant
  psmeHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48))
  psmeHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 60, a2 = 0.026, a2p = 0.65, a3 = 0.09, a3p = -0.07, b1 = 0.05, b2 = -0.017, b3 = 0.02, b3p = -0.05, b4 = 1.35, b4p = -0.22)) # a2, b1 not significant
  psmeHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050))
  psmeHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1))
  psmeHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 60, a2 = 0, a2p = 0.68, a3 = 0.066, b1 = -0.0053, b1p = -0.0042, b2 = 1.3, b2p = -0.20)) # a2, a3p not significant
  
  # GAMs
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2)
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2)
  psmeHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 37, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint)
  psmeHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", TotalHt ~ s(DBH, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 40, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint)
  if (psmeOptions$fitPhysioGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2) # long fit time (k = 455, edf = 405) AIC 140007: 140237 without BAL, 140249 without BA, 140335 without elevation, 140504 without slope, 140129 without sin(aspect), 140126 without cos(aspect), 140979 without topographic shelter -> force drop of cos(aspect) as least disadvantageous option towards viable fitting times -> 140429 without BA, 140496 without BAL, 140448 without elevation, 140589 without slope, 140279 without sin(aspect), 141178 without topographic shelter
    psmeHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", TotalHt ~ s(DBH, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), topographicShelterIndex, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint) # long fit time (k = 496, edf = 360) -> force drop of BA
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2)
    psmeHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeDiameter, bs = "ts", k = 331, by = as.factor(isPlantation), pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint)

    psmeHeightFromDiameterGamBalPhysio = psmeHeightFromDiameter$gamBalPhysio
    psmeHeightFromDiameterGamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh
    psmeHeightFromDiameterGamPhysio = psmeHeightFromDiameter$gamPhysio
    psmeHeightFromDiameterGamRelDbhPhysio = psmeHeightFromDiameter$gamRelDbhPhysio
    save(file = "trees/height-diameter/data/PSME TotalHt primary GAMs.Rdata", psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamBalPhysioRelDbh, psmeHeightFromDiameterGamPhysio, psmeHeightFromDiameterGamRelDbhPhysio)
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamBalPhysioRelDbh, psmeHeightFromDiameterGamPhysio, psmeHeightFromDiameterGamRelDbhPhysio)
  } else {
    load("trees/height-diameter/data/PSME TotalHt primary GAMs.Rdata")
    psmeHeightFromDiameter$gamBalPhysio = psmeHeightFromDiameterGamBalPhysio
    psmeHeightFromDiameter$gamBalPhysioRelDbh = psmeHeightFromDiameter$gamBalPhysioRelDbh
    psmeHeightFromDiameter$gamPhysio = psmeHeightFromDiameterGamPhysio
    psmeHeightFromDiameter$gamRelDbhPhysio = psmeHeightFromDiameterGamRelDbhPhysio
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamPhysio)
  }

  #gamPrimary = fit_gam("REML GAM primary",
  #                     DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint) +
  #                           s(standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint) +
  #                           s(tallerApproxBasalArea, bs = "ts", by = as.factor(isPlantation), k = 18, pc = gamConstraint) +
  #                           s(slope, bs = "ts", k = 10, pc = gamConstraint) +
  #                           s(elevation, bs = "ts", k = 20, pc = gamConstraint) +
  #                           #s(sin(3.14159/180 * aspect), bs = "ts", k = 30, pc = gamConstraint) + # fit lacking plausible mechanism, significant but model efficiency increases slightly when dropped
  #                           #s(topographicShelterIndex, bs = "ts", k = 40, pc = gamConstraint) + # 35.9 EDF, fit lacking plausible mechanism, significant but model efficiency increases slightly when dropped
  #                           s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 18, pc = gamConstraint),
  #                     data = psme2016physio %>% mutate(relativeHeight = pmin(relativeHeight, 1.5)), constraint = psme2016gamConstraint, folds = 1, repetitions = 1)
  #par(mfrow = c(4, 3))
  #plot(gamPrimary, ylim = c(-8, 8))
  
  save(file = "trees/height-diameter/data/PSME TotalHt primary.rdata", psmeHeightFromDiameter)
}

if (psmeOptions$fitHeightNlrobAndFixedWeight)
{
  # nlrob()
  psmeHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))) # b1p not significant
  psmeHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65, a1p = -14.8, a2 = 0.04, a2p = 0.51, a3 = 0.069, a3p = -0.18, b1 = -0.017, b2 = 1.33, b2p = -0.16)) # a3 not significant, step factor with b1p
  psmeHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 62, a2 = 0.02, a2p = 0.7, a3 = 0.10, a4 = -0.006, a5 = -0.13, a6 = 0.9, a7 = 1.25, a8 = 0.30, b1 = -0.022, b1p = 0.008, b2 = 1.45, b2p = -0.41)) # a4 not significant
  psmeHeightFromDiameterNlrob$chapmanRichardsBalPhysioRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 61, a2 = 0, a2p = 0.75, a3 = 0.12, a4 = -0.007, a5 = -0.13, a6 = 0.98, a7 = 1.3, a8 = 0.3, a10 = -0.4, a10p = 1.2, b1 = -0.023, b1p = 0.010, b2 = 1.54, b2p = -0.49), control = nls.control(tol = 0.001)) # job step factor
  psmeHeightFromDiameterNlrob$chapmanRichardsBalRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64, a1p = -25, a2 = 0, a2p = 0.48, a3 = 0.07, a3p = 0.11, a10 = -0.8, a10p = 2.1, b1 = -0.021, b2 = 1.48, b2p = -0.32), control = nls.control(tol = 0.001)) # job step factor
  psmeHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -14.4, a4 = -0.0045, a5 = -8.09, a6 = 0.8, a7 = 0.9, a8 = 0.21, b1 = -0.022, b2 = 1.55, b2p = -0.34)) # a4p not significant, a5p induces overfitting
  psmeHeightFromDiameterNlrob$chapmanRichardsRelDbh = fit_nlrob("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64, a1p = -25, a2 = 0, a2p = 0.5, a3 = 0.07, a3p = 0.11, a10 = -0.77, a10p = 2.1, b1 = -0.020, b2 = 1.48, b2p = -0.28), control = nls.control(tol = 0.001)) # job step factor
  psmeHeightFromDiameterNlrob$chapmanRichardsRelDbhPhysio = fit_nlrob("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 72, a1p = -18, a4 = -0.0045, a5 = -8.6, a6 = 0.86, a7 = 0.96, a8 = 0.23, a10 = -1.1, a10p = 1.4, b1 = -0.022, b2 = 1.55, b2p = -0.36))
  psmeHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
  psmeHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28))
  psmeHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = nls.control(maxiter = 500)) # a1p parameter evaporation, b1p not significant
  psmeHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30)) # b1p not significant
  psmeHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
  psmeHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
  psmeHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
  psmeHeightFromDiameterNlrob$richardsW = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
  psmeHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 15, b1 = 0.36, b1p = -0.08, b2 = -0.022, b2p = -0.022, b3 = 0, b3p = -0.10, b4 = 1.5, b4p = -0.32), control = nls.control(maxiter = 100, tol = 0.001)) # job step factor
  psmeHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 22, a1p = -2.9, b1 = 0.27, b2 = -0.029, b3 = -0.06, b4 = 1.75, b4p = -0.60), control = nls.control(maxiter = 100, tol = 1E-4)) # job step factor
  psmeHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 22, a1p = -5.3, a4 = -0.0023, a5 = -0.04, a6 = 0.21, a7 = 0.21, a8 = 0.062, b1 = 0.273, b2 = -0.023, b2p = -0.016, b3 = -0.013, b3p = -0.08, b4 = 1.52, b4p = -0.39)) # b1p, b3 not significant
  psmeHeightFromDiameterNlrob$sharmaPartonBalPhysioRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 25, a1p = -9, a4 = -0.0024, a5 = -0.03, a6 = 0.28, a7 = 0.19, a8 = 0.06, a10 = -0.35, a10p = 1.05, b1 = 0.26, b2 = -0.021, b2p = -0.034, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.35))
  psmeHeightFromDiameterNlrob$sharmaPartonBalRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3 * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 3.8, a1p = 16, a10 = 0.16, a10p = 0.37, b1 = 0.63, b1p = -0.45, b2 = -0.035, b3 = -0.11, b4 = 1.9, b4p = -0.7))
  psmeHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 20, a4 = -0.002, a5 = -0.03, a6 = 0.19, a7 = 0.22, a8 = 0.07, b1 = 0.30, b2 = -0.023, b3 = 0.007, b3p = -0.12, b4 = 1.5, b4p = -0.50), control = nls.control(maxiter = 100, tol = 0.001)) # a1p, b1p, b2p not significant, b3 debatable, nlrob() job step factor and step factor with a1p
  psmeHeightFromDiameterNlrob$sharmaPartonRelDbh = fit_nlrob("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 13.3, a10 = -0.2, a10p = 0.9, b1 = 0.40, b1p = -0.14, b2 = -0.020, b2p = -0.044, b3 = 0.030, b3p = -0.21, b4 = 1.48, b4p = -0.35))
  psmeHeightFromDiameterNlrob$sharmaPartonRelDbhPhysio = fit_nlrob("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 15.5, a4 = -0.0021, a5 = -0.03, a7 = 0.19, a8 = 0.06, a10 = 0.7, a10p = -0.5, b1 = 0.30, b2 = -0.029, b3 = -0.08, b3p = -0.06, b4 = 1.52, b4p = -0.51))
  psmeHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 50, a1p = -25, b1 = 0.07, b1p = 0.11, b2 = -0.03, b2p = -0.01, b3 = 0, b3p = -0.06, b4 = 1.5, b4p = -0.48)) # b3 not significant
  psmeHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 49, a2 = 0.02, a2p = 0.46, a3 = 0, a3p = -0.055, b1 = 0.08, b2 = -0.017, b3 = 0.020, b3p = -0.05, b4 = 1.37, b4p = -0.28), control = nls.control(maxiter = 100, tol = 0.001)) # a1p, a2, a3, b1p, b2p, b3 not significant, job step factor
  psmeHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050))
  psmeHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1))
  psmeHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 58, a2 = 0, a2p = 0.72, a3 = 0.08, b1 = -0.005, b1p = -0.004, b2 = 1.3, b2p = -0.22), control = nls.control(maxiter = 100, tol = 1E-4)) # a2, a3 debatably significant, job step factor
  #psmeHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a9 = 54.0, a9p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37))
  #psmeHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.5)) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 3.0, a2 = 0.22, a2p = 0.95, a3 = 0.12, a9 = 53, b1 = -0.06, b1p = 0.04, b2 = 0.75, b2p = 0.1)) # a3p, a9p not significant
  #lapply(psmeHeightFromDiameter$sharmaPartonBalPhysio$fit, confint_nlrob)
  
  # gsl_nls(weights = default)
  psmeHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31)))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 72.9, a1p = -11.8, a2 = 0.087, a2p = 0.84, a3 = -0.0021, a3p = -0.073, b1 = -0.016, b2 = 1.26, b2p = -0.054))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 74.3, a2 = 0.096, a2p = 0.92, a3 = 0, a4 = -0.015, a5 = -0.101, a6 = 0.793, a7 = 1.695, a8 = 0.183, b1 = -0.018, b1p = 0.005, b2 = 1.30, b2p = -0.154))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 61, a2 = 0, a2p = 0.72, a3 = 0.11, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, a10 = -0.4, a10p = 1.2, b1 = -0.023, b1p = 0.010, b2 = 1.54, b2p = -0.49))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 64, a1p = -23, a2 = 0, a2p = 0.48, a3 = 0.07, a3p = 0.09, a10 = -0.77, a10p = 2.1, b1 = -0.020, b2 = 1.44, b2p = -0.28))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 64, a1p = -22, a2 = 0, a2p = 0.5, a3 = 0.07, a3p = 0.09, a10 = -0.77, a10p = 2.1, b1 = -0.020, b2 = 1.43, b2p = -0.28))
  psmeHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 72, a1p = -18, a4 = -0.0045, a5 = -8.6, a6 = 0.7, a7 = 0.8, a8 = 0.23, a10 = -1.2, a10p = 1.6, b1 = -0.022, b2 = 1.5, b2p = -0.36))
  psmeHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016defaultWeight, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
  psmeHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28))
  psmeHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = gsl_nls_control(maxiter = 500))
  psmeHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016defaultWeight, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30))
  psmeHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016defaultWeight, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
  psmeHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016defaultWeight, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
  psmeHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
  psmeHeightFromDiameterGslNlsDefault$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016defaultWeight, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
  psmeHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22))
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 20, a1p = -9, b1 = 0.3, b1p = 0.07, b2 = -0.023, b3 = 0.0, b4 = 1.53, b4p = -0.15))
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51))
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3 * DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 4.2, a1p = 14.5, a10 = 0.16, a10p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66)) # a2, a2p not significant
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a10 = -0.35, a10p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40)) # b3 not significant
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65))
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation) * DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 13.3, a10 = -0.2, a10p = 0.73, b1 = 0.40, b1p = -0.12, b2 = -0.020, b2p = -0.0233, b3 = 0.032, b3p = -0.20, b4 = 1.5, b4p = -0.39))
  psmeHeightFromDiameterGslNlsDefault$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 15.5, a4 = -0.0016, a5 = -0.028, a7 = 0.17, a8 = 0.06, a10 = 0.48, a10p = -0.3, b1 = 0.33, b2 = -0.027, b3 = -0.043, b3p = -0.074, b4 = 1.52, b4p = -0.55)) # a6-a7 not reliably mutually significant
  psmeHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48))
  psmeHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 54, a2 = 0.01, a2p = 0.57, a3 = 0, a3p = -0.06, b1 = 0.05, b2 = -0.017, b3 = -0.02, b3p = -0.045, b4 = 1.36, b4p = -0.25))
  psmeHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050))
  psmeHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1))
  psmeHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = list(a1 = 58, a2 = 0, a2p = 0.70, a3 = 0.084, b1 = -0.005, b1p = -0.003, b2 = 1.3, b2p = -0.19))

  save(file = "trees/height-diameter/data/PSME TotalHt nlrob gsl_nls default.rdata", psmeHeightFromDiameterNlrob, psmeHeightFromDiameterGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
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
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$richardsW), color = "unified Richards", group = psme2016$isPlantation)) +
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
                          richards = -residuals(psmeHeightFromDiameter$richardsW),
                          sharmaParton = -residuals(psmeHeightFromDiameter$sharmaParton), 
                          sharmaPartonBal = -residuals(psmeHeightFromDiameter$sharmaPartonBal),
                          sharmaZhang = -residuals(psmeHeightFromDiameter$sharmaZhang),
                          sharmaZhangBal = -residuals(psmeHeightFromDiameter$sharmaZhangBal),
                          sibbesen = -residuals(psmeHeightFromDiameter$sibbesen),
                          weibull = -residuals(psmeHeightFromDiameter$weibull),
                          weibullBal = -residuals(psmeHeightFromDiameter$weibullBal),
                          sharmaParton = -residuals(psmeHeightFromDiameter$sharmaPartonNls)),
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
# Typically fragile and prone to convergence failure.
if (psmeOptions$fitHeightGnls)
{
  psmeHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), control = gnlsControl(nlsTol = 0.1, msTol = 1E-4, tolerance = 1E-3, maxIter = 250, nlsMaxIter = 50)))
  psmeHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 70, a1p = -8, a2 = 0.1, a2p = 0.73, a3 = -0.02, a3p = -0.17, b1 = -0.016, b2 = 1.3, b2p = -0.06), control = gnlsControl(nlsTol = 0.1, msTol = 0.01, tolerance = 0.1))
  psmeHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22), control = gnlsControl(nlsTol = 0.1, maxIter = 500, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5))
  psmeHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 20, a1p = -9, b1 = 0.3, b1p = 0.07, b2 = -0.023, b2p = -0.02, b3 = 0.0, b4 = 1.53, b4p = -0.15), control = gnlsControl(nlsTol = 1, maxIter = 250, nlsMaxIter = 50, msTol = 0.01, tolerance = 0.1))
  #psmeHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 54, a2 = 0.03, a2p = 0.7, b1 = 0.05, b2 = -0.03, b3 = -0.07, b3p = -0.09, b4 = 1.55, b4p = -0.5), control = gnlsControl(nlsTol = 1, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, msTol = 1, tolerance = 1)) # always step halving
  psmeHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 65, a2 = 0.06, a2p = 0.8, b1 = 0.01, b2 = -0.025, b3 = -0.03, b3p = -0.08, b4 = 1.5, b4p = -0.35), control = gnlsControl(nlsTol = 0.5, maxIter = 250, nlsMaxIter = 50, msTol = 1E-7, tolerance = 1E-6))
  psmeHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4))
  psmeHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 62.5, a2 = 0.04, a2p = 0.85, a3 = 0.02, a3p = -0.2, b1 = -0.0045, b1p = -0.0026, b2 = 1.3, b2p = -0.15), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-3, tolerance = 0.01))
  save(file = "trees/height-diameter/data/PSME TotalHt GNLS.rdata", psmeHeightFromDiameterGnls)
}
if (htDiaOptions$includeInvestigatory)
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


## psme height-diameter nlme regressions
if (psmeOptions$fitHeightMixed)
{
  psmeHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, 
                                                                fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))))
  psmeHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, 
                                                            fixedFormula = a1 + a1p + a2 + a2p + a3 + a3p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 65, a1p = -9, a2 = 0.07, a2p = 0.6, a3 = -0.07, a3p = -0.05, b1 = -0.016, b2 = 1.25, b2p = -0.07)))
  psmeHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, 
                                                                  fixedFormula = a1 + a2 + a2p + a3 + a4 + a5 + a6 + a7 + a8 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 63, a2 = 0.03, a2p = 0.67, a3 = 0.084, a4 = -0.005, a5 = -0.14, a6 = 0.8, a7 = 1.1, a8 = 0.3, b1 = -0.021, b1p = 0.009, b2 = 1.5, b2p = -0.4)), control = nlmeControl(maxIter = 100, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # job step factor
  psmeHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, 
                                                               fixedFormula = a1 + a1p + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31)))
  psmeHeightFromDiameterMixed$curtis = fit_nlme("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation + a1r) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, 
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                                start = list(fixed = c(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156)))
  psmeHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28)))
  psmeHeightFromDiameterMixed$korf = fit_nlme("Korf", TotalHt ~ 1.37 + (a1 + a1r)*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, 
                                              fixedFormula = a1 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                              start = list(fixed = c(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # max iterations, job max iterations
  psmeHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation + a1r)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, 
                                                         fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30)))
  psmeHeightFromDiameterMixed$prodan = fit_nlme("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation + a3r), psme2016, 
                                                fixedFormula = a1 + a2 + a2p + a3 + a3p ~ 1, randomFormula = a3r ~ 1,
                                                start = list(fixed = c(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6)))
  psmeHeightFromDiameterMixed$power = fit_nlme("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*DBH^(b1 + b1p * isPlantation), psme2016, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14)))
  psmeHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, 
                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52)))
  psmeHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, 
                                                   fixedFormula = Ha + Hap + d + dp + kU + kUp ~ 1, randomFormula = Har ~ 1,
                                                   start = list(fixed = c(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126)), control = nlmeControl(tolerance = 0.001, pnlsTol = 0.01, msTol = 1E-4)) # nlminb() false convergence
  #psmeHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, 
  #                                                    fixedFormula = a1 + b1 + b1p + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                    start = list(fixed = c(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22)), control = nlmeControl(maxIter = 100, tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # step halving
  #psmeHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, 
  #                                                       fixedFormula = a1 + a1p + b1 + b1p + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                       start = list(fixed = c(a1 = 6, a1p = 13, b1 = 0.6, b1p = -0.33, b2 = -0.025, b3 = 0.0, b4 = 1.53, b4p = -0.15)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 1E-4)) # step halving
  #psmeHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, 
  #                                                             fixedFormula = a1 + a1p + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                             start = list(fixed = c(a1 = 20, a1p = -2.8, a4 = -0.0016, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.07, b1 = 0.30, b2 = -0.035, b2p = -0.007, b3 = -0.003, b3p = -0.07, b4 = 1.57, b4p = -0.51)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 1E-3)) # step halving
  #psmeHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1r + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, 
  #                                                          fixedFormula = a1 + a4 + a5 + a6 + a7 + a8 + b1 + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                          start = list(fixed = c(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65)), control = nlmeControl(tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # computationally singular
  #psmeHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, 
  #                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                   start = list(fixed = c(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving, singularity in backsolve
  psmeHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, 
                                                        fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 60, a2 = 0.026, a2p = 0.65, a3 = 0.09, a3p = -0.07, b1 = 0.05, b2 = -0.017, b3 = 0.02, b3p = -0.05, b4 = 1.35, b4p = -0.22)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # job max iterations
  #psmeHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, 
  #                                                fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
  #                                                start = list(fixed = c(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  psmeHeightFromDiameterMixed$weibull = fit_nlme("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1)))
  psmeHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, 
                                                    fixedFormula = a1 + a2+ a2p + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 60, a2 = 0, a2p = 0.68, a3 = 0.066, b1 = -0.0053, b1p = -0.0042, b2 = 1.3, b2p = -0.20)), control = nlmeControl(maxIter = 100, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # job step factor
  
  psmeHeightFromDiameterMixed$gamm = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15) + s(StandID, bs = "re"), data = psme2016, mixed = TRUE)
  psmeHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26) + s(StandID, bs = "re"), data = psme2016, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/PSME TotalHt mixed.Rdata", psmeHeightFromDiameterMixed)
}


## Douglas-fir diameter-height regressions
if (psmeOptions$fitDbhPrimary)
{
  # linear regressions
  psmeDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016))
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016)
  
  # nonlinear regressions
  if (psmeOptions$fitDbhGslNlsAndGams)
  {
    psmeDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 40, b1 = 0.03, b1p = -0.005, b2 = 0.61, b2p = 0.15)) # a1p not significant
    psmeDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07))
    psmeDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 40, a1p = -13, a9 = 5, b1 = 0.7, b2 = 0.7, b2p = 0.03))
    psmeDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -123, a1p = 53.1, b1 = 0.0085, b1p = 0.0041, b2 = 0.77))
    psmeDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79))
    psmeDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -12, a1p = -3.9, a5 = -2.2, a8 = 0.04, b1 = 0.020, b1p = 0.0054, b2 = 0.45))
    psmeDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -18, a9 = 0, b1 = 0.016, b1p = 0.0084, b2 = 0.08, b2p = 0.4), significant = FALSE)
    psmeDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^b1 / (a2 + a2p * isPlantation - (TotalHt - 1.37)^b1), psme2016, start = list(a1 = 150, a1p = -77, a2 = 50, a2p = -15, b1 = 0.72)) # b1p not significant
    psmeDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018))
    psmeDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108))
    #psmeDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053))
    #psmeDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102))
    #psmeDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068))
    psmeDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096))
    psmeDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013)) # a2p, a3p, b2 not significant, a2-a3 not mutually significant
    psmeDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012)) # a1p, a2p, a3p not significant
    psmeDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.5, a2 = -0.02, a9 = 0.2, b1 = 0.91, b1p = -0.16, b2 = 0, b2p = 0.009)) # a9p not significant
    psmeDiameterFromHeight$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.4, a2 = -0.013, a3 = -0.0037, a6 = -0.05, a9 = 0.1, b1 = 0.94, b1p = -0.16, b2 = 0.001, b2p = 0.01), significant = FALSE) # a1p, a2p, a3p, a9 not significant
    psmeDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008)) # a1p, a4, a5, a7, a8 not significant
    psmeDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.2, a9 = 0.37, b1 = 0.85, b1p = -0.11, b2 = 0.003, b2p = 0.007)) # a9p not significant
    psmeDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.2, a6 = -0.06, a9 = 0.35, b1 = 0.86, b1p = -0.11, b2 = 0.0033, b2p = 0.007))
    #psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1), 0.9999)), psme2016, start = list(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90))
    psmeDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13), control = nlmeControl(maxIter = 250))
    psmeDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
    psmeDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021)) # a3p not significant
    psmeDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077))  # a2p, a3p not significant
    psmeDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.01, a3 = -0.006, a9 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021), significant = FALSE) # a9, a9p not significant
    psmeDiameterFromHeight$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 3.1, a1p = -0.5, a2 = -0.01, a3 = -0.005, a6 = -0.06, a9 = 0.1, b1 = 0.65, b2 = 0.07), significant = FALSE) # a2p, a3p, a9 not significant
    psmeDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09))  # a4, a5, a7, a8, b2p not significant
    psmeDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.2, a1p = -0.65, a9 = 0.43, b1 = 0.62, b2 = 0.07, b2p = 0.04))
    psmeDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 3.0, a1p = -0.5, a6 = -0.05, a9 = 0.4, b1 = 0.62, b2 = 0.07))
    psmeDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81))
    # GAMs
    # individual term selection: TotalHt + ABA + AAT + RelHt by = isPlantation + slope + elevation + sin(aspect) + TSI + RelHt
    #   primary effects NSE 0.861
    psmeDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2)
    psmeDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2)
    psmeDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 13, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2)
  }
  
  if (psmeOptions$fitPhysioGams)
  {
    psmeDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2) # slow since minimum (k = 496, edf = 397, 2x2 cross validation = 18 minutes with Zen 3, 4.6 GHz), AIC 151550: 151552 without AAT, 151700 without ABA, 151783 without elevation, 151533 without slope, 151561 without sin(aspect), 151562 without cos(aspect), 151650 without topographic shelter -> drop slope -> AIC 151800 without AAT, 151920 without ABA, 151926 without elevation, 151742 without sin(aspect), 151693 without cos(aspect), 151902 without topographic shelter
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2)
    psmeDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2)

    psmeDiameterFromHeightGamAbatPhysio = psmeDiameterFromHeight$gamAbatPhysio
    psmeDiameterFromHeightGamPhysio = psmeDiameterFromHeight$gamPhysio
    psmeDiameterFromHeightGamRelHtPhysio = psmeDiameterFromHeight$gamRelHtPhysio
    save(file = "trees/height-diameter/data/PSME DBH primary GAMs.Rdata", psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio)
    rm(psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio)
  } else {
    load("trees/height-diameter/data/PSME DBH primary GAMs.Rdata")
    psmeDiameterFromHeight$gamAbatPhysio = psmeDiameterFromHeightGamAbatPhysio
    psmeDiameterFromHeight$gamPhysio = psmeDiameterFromHeightGamPhysio
    psmeDiameterFromHeight$gamRelHtPhysio = psmeDiameterFromHeightGamRelHtPhysio
    rm(psmeDiameterFromHeightGamAbatPhysio, psmeDiameterFromHeightGamPhysio, psmeDiameterFromHeightGamRelHtPhysio)
  }
  
  if (psmeOptions$fitAbatRelHtPhysioGam)
  {
    psmeDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, nthreads = 2) # slow due to k = 496, removing either ABA or AAT increases AIC
    psmeDiameterFromHeightGamAbatPhysioRelHt = psmeDiameterFromHeight$gamAbatPhysioRelHt
    save(file = "trees/height-diameter/data/PSME DBH primary GAM ABA+T RelHt physio.Rdata", psmeDiameterFromHeightGamAbatPhysioRelHt)
    rm(psmeDiameterFromHeightGamAbatPhysioRelHt)
  } else {
    load("trees/height-diameter/data/PSME DBH primary GAM ABA+T RelHt physio.Rdata")
    psmeDiameterFromHeight$gamAbatPhysioRelHt = psmeDiameterFromHeightGamAbatPhysioRelHt
    rm(psmeDiameterFromHeightGamAbatPhysioRelHt)
  }
  
  save(file = "trees/height-diameter/data/PSME DBH primary.rdata", psmeDiameterFromHeight)
}

if (htDiaOptions$includeInvestigatory)
{
  print(psmeDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman replace", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanReplaceAbat), y = TotalHt, color = "Chapman replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$linear), y = TotalHt, color = "linear", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$parabolic), y = TotalHt, color = "parabolic", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$schnute), y = TotalHt, color = "Schnute inverse", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen replace", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -65 * log(1 - pmin((1/85*(TotalHt - 1.37))^0.7, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = -155 * log(1 - pmin((0.00839*(TotalHt - 1.37))^0.970, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = 0.05*topHeight*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 3.5*sqrt(TotalHt - 1.37) / (1 - 0.1*sqrt(TotalHt - 1.37)), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = -1/0.003*log(1 - (1 - exp(-0.55))*(TotalHt^1.05 - 1.37^1.05)/(90^1.05 - 1.3^1.05)), y = TotalHt, color = "Schnute", group = isPlantation), alpha = 0.5) +
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

if (psmeOptions$fitDbhNlrobAndFixedWeight)
{
  psmeDiameterFromHeightNlrob = list(chapmanReplace = fit_nlrob("Chapman-Richards replace", DBH ~ a1*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 42, b1 = 0.028, b1p = -0.003, b2 = 0.66, b2p = 0.11)))
  psmeDiameterFromHeightNlrob$chapmanReplaceAbat = fit_nlrob("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 80, a1p = -40, a2 = -0.52, a3 = -0.08, b1 = 0.016, b1p = 0.09, b2 = 0.74, b2p = 0.04)) # a2, a2p, a3p not significant
  #psmeDiameterFromHeightNlrob$chapmanReplaceBal = fit_nlrob("Chapman-Richards replace BA+L", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 135, a1p = -37.5, a2 = -1.2, b1 = 0.010, b1p = 0.002, b2 = 0.756, b2p = 0.064), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a3 * BA
  #psmeDiameterFromHeightNlrob$chapmanReplaceBalRelHt = fit_nlrob("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 139, a1p = -38.5, a2 = -1.2, a9 = -5.5, a9p = 0.22, b1 = 0.012, b1p = 0.003, b2 = 0.796, b2p = 0.066), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a9 * BA
  psmeDiameterFromHeightNlrob$chapmanReplaceRelHt = fit_nlrob("Chapman-Richards replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 40, a1p = -10, a9 = 8, b1 = 0.08, b2 = 0.66, b2p = 0.02)) # a9p not significant
  psmeDiameterFromHeightNlrob$chapmanRichards = fit_nlrob("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -110, a1p = 40, b1 = 0.009, b1p = 0.0036, b2 = 0.77), control = nls.control(maxiter = 250)) # b2p not significant
  psmeDiameterFromHeightNlrob$chapmanRichardsAbat = fit_nlrob("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -160, a1p = 80, a2 = 1.2, a3 = 0.12, b1 = 0.007, b1p = 0.005, b2 = 0.79)) # a2p, a3p, b2p not significant
  psmeDiameterFromHeightNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -13, a1p = -45, a8 = 0.05, b1 = 0.0198, b1p = -0.0057, b2 = 0.74), control = nls.control(maxiter = 500)) # a4, a5, a6, a7 not significant
  psmeDiameterFromHeightNlrob$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -7.4, a9 = -4.6, b1 = 0.027, b2 = 0.09, b2p = 0.2)) # a1p, b1p not significant
  psmeDiameterFromHeightNlrob$michaelisMentenReplace = fit_nlrob("Michaelis-Menten replace", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^b1 / (a2 + a2p * isPlantation - (TotalHt - 1.37)^b1), psme2016, start = list(a1 = 140, a1p = -55, a2 = 50, a2p = -11, b1 = 0.74))
  psmeDiameterFromHeightNlrob$naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 4.95, a1p = -2.05, a2 = -0.085, a2p = -0.03), control = nls.control(maxiter = 250))
  psmeDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.2, a1p = 1.0, b1 = 1.1, b1p = -0.23))
  #psmeDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.5, a1p = 0.9, a2 = -0.003, a2p = -0.008, b1 = 1.05, b1p = -0.17))
  #psmeDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.2, a1p = 1.0, a6 = -0.03, b1 = 1.1, b1p = -0.22))
  #psmeDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 2.2, a9 = 0.6, b1 = 0.89, b1p = -0.09)) # a1p, a9p not significant
  psmeDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.42, b1 = 0.82, b1p = -0.08, b2 = 0.0069, b2p = 0.005))
  psmeDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.4, a2 = -0.027, b1 = 0.91, b1p = -0.15, b2 = 0.0015, b2p = 0.01))
  psmeDiameterFromHeightNlrob$ruarkAbatPhysio = fit_nlrob("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.7, a2 = -0.017, a3 = -0.0042, a6 = -0.04, b1 = 0.89, b1p = -0.19, b2 = 0.0025, b2p = 0.012))
  psmeDiameterFromHeightNlrob$ruarkAbatRelHt = fit_nlrob("Ruark ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.4, a2 = -0.02, a9 = 0.3, b1 = 0.9, b1p = -0.14, b2 = 0.0008, b2p = 0.008), control = nls.control(tol = 1E-4)) # job step factor
  psmeDiameterFromHeightNlrob$ruarkAbatRelHtPhysio = fit_nlrob("Ruark ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.6, a2 = -0.009, a3 = -0.004, a6 = -0.04, a9 = 0.04, b1 = 0.9, b1p = -0.12, b2 = 0.004, b2p = 0.007), significant = FALSE)
  psmeDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.5, a6 = -0.04, b1 = 0.82, b1p = -0.07, b2 = 0.007, b2p = 0.005))
  psmeDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.2, a9 = 0.5, b1 = 0.83, b1p = -0.11, b2 = 0.0035, b2p = 0.006))
  psmeDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio, start = list(a1 = 2.2, a6 = -0.06, a9 = 0.5, b1 = 0.83, b1p = -0.11, b2 = 0.004, b2p = 0.007))
  #psmeDiameterFromHeightNlrob$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 18.6)) # converges from red alder values but fails to reconverge (singular gradient or NaN-inf with nls())
  #psmeDiameterFromHeightNlrob$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b4 = 0.3, b4p = 0.13), control = nls.control(maxiter = 500)) # b3p not significant, a1p NaN-inf (not significant?), b4p significance debatable, singular gradient with all relative height forms attempted, >1000 iterations with nlrob()
  psmeDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 2.9, a1p = -0.2, b1 = 0.64, b2 = 0.08, b2p = 0), control = nls.control(tol = 1E-4)) # b2p not significant, job step factor
  psmeDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 2.63, a2 = -0.007, a3 = -0.005, b1 = 0.67, b2 = 0.08, b2p = -0.014))
  psmeDiameterFromHeightNlrob$sibbesenReplaceAbatPhysio = fit_nlrob("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.2, a1p = -0.5, a2 = -0.010, a3 = -0.0068, a6 = -0.06, b1 = 0.63, b2 = 0.088))
  psmeDiameterFromHeightNlrob$sibbesenReplaceAbatRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.4, a1p = -0.9, a2 = -0.014, a3 = -0.005, a9 = 0.2, b1 = 0.67, b2 = 0.063, b2p = 0.01), significant = FALSE)
  psmeDiameterFromHeightNlrob$sibbesenReplaceAbatRelHtPhysio = fit_nlrob("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.2, a1p = -0.5, a2 = -0.01, a3 = -0.006, a6 = -0.05, a9 = 0.09, b1 = 0.66, b2 = 0.07), significant = FALSE)
  psmeDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 2.95, a1p = -0.32, a6 = -0.06, b1 = 0.63, b2 = 0.08))
  psmeDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016, start = list(a1 = 3.0, a1p = -0.5, a9 = 0.6, b1 = 0.62, b2 = 0.07))
  psmeDiameterFromHeightNlrob$sibbesenReplaceRelHtPhysio = fit_nlrob("Sibbesen replace RelHt physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio, start = list(a1 = 3.1, a1p = -0.5, a6 = -0.07, a9 = 0.5, b1 = 0.61, b2 = 0.071))
  psmeDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -360, a1p = 140, b1 = 0.01, b1p = 0.003, b2 = 0.81), control = nls.control(maxiter = 250))
  #confint_nlrob(psmeDiameterFromHeightNlrob$sibbesenReplaceRelHt, level = 0.99)

  psmeDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 40, b1 = 0.028, b1p = -0.0045, b2 = 0.63, b2p = 0.19)))
  psmeDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016defaultWeight, start = list(a1 = 75, a1p = -40, a2 = -0.5, a3 = -0.009, b1 = 0.019, b1p = 0.010, b2 = 0.70, b2p = 0.056))
  psmeDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016defaultWeight, start = list(a1 = 40, a1p = -10, a9 =1.67, b1 = 0.1, b2 = 0.65, b2p = 0.04))
  psmeDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016defaultWeight, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781))
  psmeDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016defaultWeight, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786))
  psmeDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016defaultWeightPhysio, start = list(a1 = -13.8, a1p = -0.65, a5 = -3.58, a8 = 0.091, b1 = 0.019, b1p = 0.0062, b2 = 0.30))
  psmeDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016defaultWeight, start = list(a1 = -10, a9 = -1.3, b1 = 0.021, b1p = 0.0025, b2 = 0.03, b2p = 0.12))
  psmeDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), psme2016defaultWeight, start = list(a1 = 190, a2 = 115, b1 = 0.91)) # a1p, a2p not significant
  psmeDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016defaultWeight, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018))
  psmeDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016defaultWeight, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108))
  #psmeDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016defaultWeight, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053))
  #psmeDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016defaultWeightPhysio, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102))
  #psmeDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016defaultWeight, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068))
  psmeDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeight, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096))
  psmeDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeight, start = list(a1 = 2.6, a2 = -0.005, b1 = 0.88, b1p = -0.133, b2 = 0.003, b2p = 0.009))
  psmeDiameterFromHeightGslNlsDefault$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.7, a2 = -0.017, a3 = -0.0042, a6 = -0.04, b1 = 0.89, b1p = -0.19, b2 = 0.0025, b2p = 0.012))
  psmeDiameterFromHeightGslNlsDefault$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeight, start = list(a1 = 2.7, a2 = -0.019, a9 = 0.08, b1 = 0.88, b1p = -0.16, b2 = 0.002, b2p = 0.01))
  psmeDiameterFromHeightGslNlsDefault$ruarkAbatRelHtPhysio = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.6, a2 = -0.01, a3 = -0.0044, a6 = -0.04, a9 = 0, b1 = 0.90, b1p = -0.16, b2 = 0.003, b2p = 0.01), significant = FALSE)
  psmeDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008))
  psmeDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeight, start = list(a1 = 2.4, a9 = 0.1, b1 = 0.85, b1p = -0.11, b2 = 0.005, b2p = 0.007))
  psmeDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016defaultWeightPhysio, start = list(a1 = 2.4, a6 = -0.04, a9 = 0.11, b1 = 0.84, b1p = -0.10, b2 = 0.0052, b2p = 0.007))
  #psmeDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016defaultWeight, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 86)) # NaN-inf
  psmeDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016defaultWeight, start = list(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13), control = nlmeControl(maxiter = 250))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.2, a1p = -0.5, a2 = -0.010, a3 = -0.007, a6 = -0.06, b1 = 0.63, b2 = 0.088))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 4.6, a1p = -1.5, a2 = -0.014, a3 = -0.008, a9 = 0, b1 = 0.55, b2 = 0.094, b2p = 0.026), significant = FALSE) # a9, a9p not significant
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatRelHtPhysio = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.1, a1p = -0.5, a2 = -0.01, a3 = -0.006, a6 = -0.05, a9 = 0.08, b1 = 0.7, b2 = 0.07), significant = FALSE)
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.2, a1p = -0.24, a6 = -0.05, b1 = 0.55, b2 = 0.111))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = list(a1 = 3.8, a1p = -0.9, a9 = 0.18, b1 = 0.53, b2 = 0.10, b2p = 0.014))
  psmeDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016defaultWeightPhysio, start = list(a1 = 3.0, a1p = -0.4, a6 = -0.05, a9 = 0.2, b1 = 0.6, b2 = 0.08))
  psmeDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016defaultWeight, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81))
  
  save(file = "trees/height-diameter/data/PSME DBH nlrob gsl_nls default.Rdata", psmeDiameterFromHeightNlrob, psmeDiameterFromHeightGslNlsDefault)
}

if (psmeOptions$fitDbhMixed)
{
  psmeDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", DBH ~ (a1 + a1r) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, 
                                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 40, b1 = 0.03, b1p = -0.005, b2 = 0.61, b2p = 0.15)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 0.1, msTol = 0.001))) # max iterations
  psmeDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, 
                                                            fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 75, a1p = -35, a2 = 0.5, a3 = -0.07, b1 = 0.018, b1p = 0.01, b2 = 0.7, b2p = 0.07)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 0.1, msTol = 0.001)) # max iterations
  #psmeDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", DBH ~ (a1 + a1r + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, 
  #                                                           fixedFormula = a1 + a1p + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
  #                                                           start = list(fixed = c(a1 = 30, a1p = -10, a9 = 5, b1 = 0.12, b2 = 0.60, b2p = 0.035)), control = nlmeControl(tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # step halving
  #psmeDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", DBH ~ (a1 + a1r + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, 
  #                                                       fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                       start = list(fixed = c(a1 = -123, a1p = 53.1, b1 = 0.0085, b1p = 0.0041, b2 = 0.77)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations, false convergence
  psmeDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = -175, a1p = 100, a2 = 1.2, a3 = 0.14, b1 = 0.007, b1p = 0.0057, b2 = 0.79)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  #psmeDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", DBH ~ (a1 + a1r + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio,
  #                                                             fixedFormula = a1 + a1p + a5 + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                             start = list(fixed = c(a1 = -12, a1p = -3.9, a5 = -2.2, a8 = 0.04, b1 = 0.020, b1p = 0.0054, b2 = 0.45)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # false convergence
  #psmeDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, 
  #                                                            fixedFormula = a1 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
  #                                                            start = list(fixed = c(a1 = -9.0, a9 = -12.6, b1 = 0.013, b1p = 0.009, b2 = 0.004, b2p = 0.31)), control = nlmeControl(tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # step halving
  #control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.001, msTol = 1E-5)
  psmeDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", DBH ~ (a1 + a1r + a1p * isPlantation) * (TotalHt - 1.37)^b1 / (a2 + a2p * isPlantation - (TotalHt - 1.37)^b1), psme2016, 
                                                                fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = 150, a1p = -77, a2 = 50, a2p = -15, b1 = 0.72)), control = nlmeControl(maxIter = 500, tolerance = 0.1, pnlsTol = 1, msTol = 0.01))
  psmeDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", DBH ~ (a1 + a1r + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018)), control = nlmeControl(tolerance = 0.1, pnlsTol = 1, msTol = 0.001)) # false convergence
  psmeDiameterFromHeightMixed$power = fit_nlme("power", DBH ~ (a1 + a1r + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108)))
  #psmeDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", DBH ~ (a1 + a1r + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, 
  #                                                 fixedFormula = a1 + a1p + a2 + a2p + a3 + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
  #                                                 start = list(fixed = c(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053)))
  #psmeDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", DBH ~ (a1 + a1r + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio,
  #                                                   fixedFormula = a1 + a1p + a4 + a5 + a6 + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
  #                                                   start = list(fixed = c(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102)))
  #psmeDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, 
  #                                                  fixedFormula = a1 + a9 + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
  #                                                  start = list(fixed = c(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068)))
  psmeDiameterFromHeightMixed$ruark = fit_nlme("Ruark", DBH ~ (a1 + a1r + a1r) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, 
                                               fixedFormula = a1 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096)))
  psmeDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, 
                                                   fixedFormula = a1 + a2 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 2.5, a2 = -0.03, b1 = 0.92, b1p = -0.2, b2 = 0, b2p = 0.013)))
  psmeDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio,
                                                         fixedFormula = a1 + a2 + a3 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 2.5, a2 = -0.017, a3 = -0.004, a6 = -0.05, b1 = 0.93, b1p = -0.19, b2 = 0.002, b2p = 0.012)))
  psmeDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, 
                                                        fixedFormula = a1 + a2 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 2.5, a2 = -0.027, a9 = 0.3, b1 = 0.91, b1p = -0.2, b2 = 0, b2p = 0.013)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # false convergence, step factor
  psmeDiameterFromHeightMixed$ruarkAbatRelHtPhysio = fit_nlme("Ruark ABA+T RelHt physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio,
                                                              fixedFormula = a1 + a2 + a3 + a6 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 2.4, a2 = -0.017, a3 = -0.004, a6 = -0.05, a9 = 0.3, b1 = 0.92, b1p = -0.19, b2 = 0.001, b2p = 0.012)), significant = FALSE)
  psmeDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", DBH ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio,
                                                     fixedFormula = a1 + a6 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 2.5, a6 = -0.05, b1 = 0.84, b1p = -0.11, b2 = 0.005, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", DBH ~ (a1 + a1r + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, 
                                                    fixedFormula = a1 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 2.4, a9 = 0.45, b1 = 0.83, b1p = -0.13, b2 = 0.004, b2p = 0.008)))
  psmeDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", DBH ~ (a1 + a1r + a6 * cos(3.14159/180 * aspect) + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016physio,
                                                          fixedFormula = a1 + a6 + a9 + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                          start = list(fixed = c(a1 = 2.4, a6 = -0.05, a9 = 0.4, b1 = 0.83, b1p = -0.13, b2 = 0.0042, b2p = 0.008)))
  #psmeDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", DBH ~ -1/a1 * log(1 - pmin((1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1), 0.9999)), psme2016, 
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1,
  #                                               start = list(fixed = c(a1 = 0.003, a2 = 0.55, b1 = 1.05, Ha = 90)), control = nlmeControl(tolerance = 0.1, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #psmeDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", DBH ~ (a1 + a1r + a1r) * (TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, 
  #                                                    fixedFormula = a1 + b1 + b1p + b2 + b3 + b3p + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
  #                                                    start = list(fixed = c(a1 = 9, b1 = 0.4, b1p = -0.14, b2 = 0.04, b3 = -0.06, b3p = 0.11, b4 = 0.3, b4p = 0.13)), control = nlmeControl(maxIter = 250, tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # singularity in backsolve
  psmeDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", DBH ~ (a1 + a1r + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, 
                                                         fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, 
                                                             fixedFormula = a1 + a1p + a2 + a3 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, b1 = 0.61, b2 = 0.071, b2p = 0.021)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # step halving
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio,
                                                                   fixedFormula = a1 + a1p + a2 + a3 + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                   start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = -0.016, a3 = -0.016, a6 = -0.057, b1 = 0.66, b2 = 0.077)))
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, 
                                                                  fixedFormula = a1 + a1p + a2 + a3 + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 4.4, a1p = -1.6, a2 = -0.03, a3 = -0.006, a9 = 0.15, b1 = 0.61, b2 = 0.071, b2p = 0.021)), significant = FALSE)
  psmeDiameterFromHeightMixed$sibbesenReplaceAbatRelHtPhysio = fit_nlme("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1r + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio,
                                                                        fixedFormula = a1 + a1p + a2 + a3 + a6 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                        start = list(fixed = c(a1 = 3.2, a1p = -0.6, a2 = 0, a3 = 0, a6 = -0.06, a9 = 0.5, b1 = 0.60, b2 = 0.09)), significant = FALSE)
  psmeDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", DBH ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio,
                                                               fixedFormula = a1 + a1p + a6 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 3.0, a1p = -0.3, a6 = -0.06, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", DBH ~ (a1 + a1r + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, 
                                                              fixedFormula = a1 + a1p + a9 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 3.90, a1p = -0.95, a9 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # false convergence
  psmeDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", DBH ~ (a1 + a1r + a1p * isPlantation + a6 * cos(3.14159/180 * aspect) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), psme2016physio,
                                                                    fixedFormula = a1 + a1p + a6 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                    start = list(fixed = c(a1 = 3.2, a1p = -0.6, a6 = -0.06, a9 = 0.5, b1 = 0.60, b2 = 0.09)))
  psmeDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", DBH ~ ((a1 + a1r + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81)))
  
  psmeDiameterFromHeightMixed = list(gamm = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10) + s(StandID, bs = "re"), data = psme2016, mixed = TRUE))
  psmeDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28) + s(StandID, bs = "re"), data = psme2016, mixed = TRUE)
  psmeDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 13) + s(StandID, bs = "re"), data = psme2016, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/PSME DBH mixed.Rdata", psmeDiameterFromHeightMixed)
}


if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightMixed & psmeOptions$fitHeightNlrobAndFixedWeight & 
    psmeOptions$fitDbhPrimary & psmeOptions$fitDbhMixed & psmeOptions$fitDbhNlrobAndFixedWeight)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + nlrob + gsl_nls + DBH primary + nlrob = 120 GB DDR at 10x10 cross validation with models dropped
  # Therefore, batch formation height and diameter results.
  if (exists("psmeHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/PSME TotalHt primary.rdata") }
  if (exists("psmeHeightFromDiameterGslNlsDefault") == FALSE) { load("trees/height-diameter/data/PSME TotalHt nlrob gsl_nls default.rdata") }
  #if (exists("psmeHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/PSME TotalHt GNLS.rdata") }
  if (exists("psmeHeightFromDiameterMixed") == FALSE) { load("trees/height-diameter/data/PSME TotalHt mixed.Rdata") }
  if (exists("psmeDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/PSME DBH primary.rdata") }
  if (exists("psmeDiameterFromHeightGslNlsDefault") == FALSE) { load("trees/height-diameter/data/PSME DBH nlrob gsl_nls default.rdata") }
  if (exists("psmeDiameterFromHeightMixed") == FALSE) { load("trees/height-diameter/data/PSME DBH mixed.rdata") }
  
  psmeCoefficients = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_coefficients)),
                                         #bind_rows(lapply(psmeHeightFromDiameterGnls, get_model_coefficients)),
                                         bind_rows(lapply(psmeHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(psmeHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(psmeHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(psmeDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(psmeDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(psmeDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")
  psmeResults = bind_rows(bind_rows(bind_rows(lapply(psmeHeightFromDiameter, get_list_stats)),
                                    bind_rows(lapply(psmeHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(psmeHeightFromDiameterMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(psmeHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                                    #bind_rows(lapply(psmeHeightFromDiameterGnls, get_model_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(psmeDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary"),
                                    bind_rows(lapply(psmeDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    create_model_stats(name = "Schnute inverse", fitSet = "gsl_nls", fittingMethod = "gsl_nls", fixedWeight = -1),
                                    bind_rows(lapply(psmeDiameterFromHeightMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(psmeDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob")),
                                    create_model_stats(name = "Chapman-Richards replace BA+L", fitSet = "nlrob", fittingMethod = "nlrob"),
                                    create_model_stats(name = "Chapman-Richards replace BA+L RelHt", fitSet = "nlrob", fittingMethod = "nlrob"),
                                    create_model_stats(name = "Schnute inverse", fitSet = "nlrob", fittingMethod = "nlrob"),
                                    create_model_stats(name = "modified Sharma-Parton", fitSet = "nlrob", fittingMethod = "nlrob")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "PSME")

  check_plot_results(psmeResults)
  save(file = "trees/height-diameter/data/PSME results.Rdata", psmeCoefficients, psmeResults)
}


## preferred forms identified (results.R, Figure 8)
if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightNlrobAndFixedWeight & psmeOptions$fitDbhPrimary & psmeOptions$fitDbhNlrobAndFixedWeight)
{
  psmeHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, folds = 1, repetitions = 1))
  psmeHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 22, a1p = -7.6, a4 = -0.0020, a5 = -0.03, a6 = 0.14, a7 = 0.14, a8 = 0.06, a10 = -0.35, a10p = 0.79, b1 = 0.28, b2 = -0.021, b2p = -0.026, b3 = 0.02, b3p = -0.17, b4 = 1.53, b4p = -0.40), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3 * DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 4.2, a1p = 14.5, a10 = 0.16, a10p = 0.27, b1 = 0.63, b1p = -0.43, b2 = -0.025, b3 = -0.09, b4 = 1.73, b4p = -0.66), folds = 1, repetitions = 1)
  psmeHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050), folds = 1, repetitions = 1)
  #AIC(psmeHeightFromDiameterPreferred$prodan, psmeHeightFromDiameterPreferred$sibbesen)
  
  psmeDiameterFromHeightPreferred = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), folds = 1, repetitions = 1))
  psmeDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint), data = psme2016, constraint = psme2016gamConstraint, nthreads = 2, folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, folds = 1, repetitions = 1, nthreads = 4)
  psmeDiameterFromHeightPreferred$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 496, pc = gamConstraint), data = psme2016physio, constraint = psme2016gamConstraint, folds = 1, repetitions = 1, nthreads = 4)
  psmeDiameterFromHeightPreferred$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), folds = 1, repetitions = 1)
  psmeDiameterFromHeightPreferred$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), folds = 1, repetitions = 1)
  
  save(file = "trees/height-diameter/data/PSME preferred models.Rdata", psmeHeightFromDiameterPreferred, psmeDiameterFromHeightPreferred)
}


## Q-Q plots
if (htDiaOptions$includeInvestigatory)
{
  ggplot() + # symmetric t fits less well than skewed t, as expected
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanReplace), color = "Chapman-Richards form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenReplace), color = "Sibbesen form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanReplace), color = "Chapman-Richards form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenReplace), color = "Sibbesen form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'t(df = 7, '*alpha*' = 2.25)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
  ggplot() + # slow! and fits poorly
  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanReplace), color = "Chapman-Richards form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenReplace), color = "Sibbesen form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanReplace), color = "Chapman-Richards form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenReplace), color = "Sibbesen form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'pearsonIV(m = 1.72, '*nu*' = 0.656)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
}


## basal area from height
# essentially no difference between fit_gsl_nls() and fit_nlrob() fits
# Chapman-Richards has the wrong curvature
if (htDiaOptions$includeInvestigatory)
{
  psmeBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 0.5, b1 = 0.0008, b2 = 1.8, b2p = -0.07), weights = heightWeight^2) # a1p, b1p not significant
  psmeBasalAreaFromHeightPower = gsl_nls(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 0.0001, a1p = 0.0002, b1 = 2.2, b1p = -0.5), weights = heightWeight^2)
  #confint2(psmeBasalAreaFromHeightKorf, level = 0.99)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(psmeBasalAreaFromHeightKorf), 100^2 * mean(-residuals(psmeBasalAreaFromHeightKorf)), mean(abs(residuals(psmeBasalAreaFromHeightKorf))), 1 - sum(residuals(psmeBasalAreaFromHeightKorf)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2)),
          "power", AIC(psmeBasalAreaFromHeightPower), 100^2 * mean(-residuals(psmeBasalAreaFromHeightPower)), mean(abs(residuals(psmeBasalAreaFromHeightPower))), 1 - sum(residuals(psmeBasalAreaFromHeightPower)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(psme2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_line(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    #geom_line(aes(x = imputedHeight, y = 1*(exp((0.01 + 0*isPlantation)*(imputedHeight - 1.37)^(1 + 0*isPlantation)) - 1) - 0, color = "Korf")) +
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
if (htDiaOptions$includeInvestigatory)
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