# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## red alder height-diameter regression form sweep
# R options for nonlinear least squares
#   gslnls::fit_gsl_nls() - fast but fixed weights only, defaults to Levenberg-Marquadt
#   minpack.lm::nlsLM() - weights = wfcs() fails with index out of range, Levenberg-Marquadt only
#   nlme::gnls() - varPower() based weighting fragile andconvergence usually fails
#   robustbase::fit_nlrob() - defaults to iterative reweighted least squares with stats::nls()
#   stats::nls() - fixed weights only, defaults to Gauss-Newton, NL2SOL with algorithm = "port"
# Difference between stats::nls() with fixed weighting and robustbase::fit_nlrob() is relatively small but fit_nlrob() 
# does noticeably decrease asymptoticity of most regression forms (slight decreases occur in some cases).
#alruHeightFromDiameter$michaelisMenten = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH / (a2 + a2p * isPlantation + DBH), alru2016, start = list(a1 = 45.5, a1p = 18.7, a2 = 49.1, a2p = 20.5))
#alruHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, d = 1.17, kU = 0.0366))
#alruHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, kU = 0.0366, kUp = 0)) # confint2() step factor
#alruHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - d))))^(1/(1 - (d + dp*isPlantation))), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, dp = 0, kU = 0.0366, kUp = 0)) # confint2() NaN-inf
#alruHeightFromDiameter$sharmaParton = fit_nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 21.0, a2 = 0.045, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44))
#alruHeightFromDiameter$sharmaPartonBal = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 16.4, a1p = 6.7, a2 = 0.104, a2p = -0.017, b1 = -0.042, b1p = 0.003, b2 = 0.062, b2p = -0.120, b3 = 1.18, b3p = -0.133)) # a1p, a2p, b1p, b2, b2p not significant
#alruHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45)) # a1p, a2, a2p, a3, a4, a6, b1p, b2p, b3p not significant
#alruHeightFromDiameter$sharmaPartonPhysio = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46))
#alruHeightFromDiameter$sharmaZhang = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456))
#alruHeightFromDiameter$sharmaZhangBal = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370))
#alruHeightFromDiameter$sibbesen = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.0019, a1p = 0.163, b1 = 4.98, b1p = -2.72, b2 = -0.175, b2p = 0.0427))
#alruHeightFromDiameter$korf = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365))
#alruHeightFromDiameter$weibull = fit_nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16))
#alruHeightFromDiameter$weibullBal = fit_nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133))
alru2016 = trees2016 %>% filter(Species == "RA", isLiveUnbroken, TotalHt > 0) %>% # live red alders measured for height
  mutate(dbhWeight = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alru2016physio = alru2016 %>% filter(is.na(elevation) == FALSE)
alru2016gamConstraint = c(DBH = -1.5179/0.7812, TotalHt = 1.37, standBasalAreaPerHectare = median(alru2016$standBasalAreaPerHectare), basalAreaLarger = median(alru2016$basalAreaLarger), standBasalAreaApprox = median(alru2016$standBasalAreaApprox), tallerApproxBasalArea = median(alru2016$tallerApproxBasalArea), elevation = median(alru2016physio$elevation), slope = median(alru2016physio$slope), aspect = median(alru2016physio$aspect), topographicShelterIndex = median(alru2016physio$topographicShelterIndex), relativeHeight = median(alru2016$relativeHeight)) # point constraint for mgcv::s()
#alru2016natural = alru2016 %>% filter(isPlantation == FALSE)
#alru2016plantation = alru2016 %>% filter(isPlantation)
#alru2016plantationPhysio = alru2016physio %>% filter(isPlantation)

alruHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027)))
alruHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 28.5, a1p = 13.1, a2 = 0.212, a2p = 0.218, a3 = -0.173, b1 = -0.0404, b1p = 0.0198, b2 = 1.143, b2p = -0.122)) # a3p not significant
alruHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 30.1, a1p = 3.7, a2 = 0.05, a2p = 0.083, a3 = 0, a4 = -0.002, a5 = -0.28, a6 = 0.074, a7 = 0.26, a8 = 0.21, b1 = -0.052, b1p = 0.007, b2 = 1.24)) # a4, a6, a7 not significant
alruHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750)) # a2, a3, a3p, b2 not significant
alruHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 33.7, a1p = 4.73, a5 = -0.195, a8 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14)) # a4, a6, a7, b2p not significant
alruHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016, start = list(a1 = 1.6, b1 = 0.24))
alruHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 9, pc = alru2016gamConstraint), data = alru2016)
alruHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 18, pc = alru2016gamConstraint), data = alru2016)
alruHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4) # all predictors k = 451, edf < 270, AIC 18588: 18709 without BA, 18735 without BAL, 18505 without elevation, 18872 without slope, 18557 without sin(aspect), 18598 without cos(aspect), 18629 without topographic shelter -> eliminate elevation AIC 18505: 18355 without BA, 18423 without BAL, 18686 without slope, 18334 without sin(aspect), 18320 without cos(aspect), 18342 without topographic shelter -> eliminate cos(aspect)
alruHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4)
alruHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 30.7, a1p = 13.3, b1 = 49.1, b1p = 7.95, b2 = -1.25, b2p = 0.12))
alruHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^b2), alru2016, start = list(a1 = 223, a1p = 22.1, b1 = -6.04, b2 = -0.252)) # b1p, b2p not significant
alruHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), alru2016)
alruHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152)) # a2p not significant
alruHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH), alru2016) # isPlantation*DBH^2 not significant
alruHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3 + a3p* isPlantation), alru2016, start = list(a1 = 0.0249, a1p = -0.0054, a2 = 0.938, a3 = 0.761, a3p = -0.177)) # a2p not significant
alruHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047)) # a1p not significant
alruHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016, start = list(a1 = 33.6, a1p = 7.31, b1 = -22.3, b1p = -4.77, b2 = 5.31, b2p = 1.37))
alruHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361))
alruHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, alru2016, start = list(a1 = 21.3, b1 = 0.038, b1p = 0.068, b2 = -0.039, b3 = 0.084, b4 = 1.19)) # b2p, b3p, b4p not significant
alruHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016, start = list(a1 = 19.0, a1p = 2.9, b1 = 0.090, b2 = -0.064, b3 = -0.14, b4 = 1.18)) # b1p, b2p, b3p, b4p not significant
alruHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016physio, start = list(a1 = 18.1, a1p = 2.80, a5 = 0, a8 = 0.002, b1 = 0.097, b2 = -0.065, b3 = -0.131, b4 = 1.18)) # a3, a4, a6, a7, b1p, b2p, b3p, b4p not significant
alruHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, alru2016physio, start = list(a1 = 15.2, a1p = 2.54, a5 = -0.13, a8 = 0.10, b1 = 0.17, b2 = -0.068, b3 = -0.07, b4 = 1.26)) # a4, a6, a7, b1p, b2p, b3p, b4p not significant
alruHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), alru2016, start = list(a1 = 22.4, b1 = 0.029, b1p = 0.071, b2 = -0.040, b3 = -0.026, b3p = -0.053, b4 = 1.18, b4p = -0.123)) # a1p, b1, b2p not significant
alruHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, alru2016, start = list(a1 = 48.9, a1p = -21.2, a2 = 0.006, b1 = -0.178, b1p = 0.197, b2 = -0.042, b3 = -0.039, b4 = 1.07)) # a3p, b2p, b3p, b4p not significant
alruHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044))
alruHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083)) # b1p not significant
alruHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 30.4, a2 = 0.238, a3 = -0.215, a3p = 0.186, b1 = -0.0247, b2 = 1.11, b2p = -0.052)) # a1p, a2p, b1p not significant
alruHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 2.4, a2 = -0.018, a3 = 0.128, a3p = 0.128, a4 = 15.2, b1 = -0.095, b2 = 0.970, b2p = -0.219)) # a1p, a2p, a4p, b1p not significant
#confint_nlrob(alruHeightFromDiameter$sharmaPartonPhysio, level = 0.99)

alruHeightFromDiameterResults = bind_rows(lapply(alruHeightFromDiameter, as_row)) %>%
  mutate(responseVariable = "height", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(alruHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = alru2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$curtis), color = "Curtis", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$korf), color = "Korf", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$linear), color = "linear", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$parabolic), color = "parabolic", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$power), color = "power", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$prodan), color = "Prodan", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$richards), color = "unified Richards", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$sibbesen), color = "Sibbesen", group = alru2016$isPlantation)) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$weibull), color = "Weibull", group = alru2016$isPlantation)) +
    annotate("text", x = 0, y = 50, label = "red alder, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 50)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
  
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$chapmanRichards))), alpha = 0.1, color = "grey25", shape = 16) +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$chapmanRichards)), color = "Chapman-Richards GAM", fill = "Chapman-Richards GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$chapmanRichards)), color = "Chapman-Richards sqrt(DBH)", fill = "Chapman-Richards sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$chapmanRichards)), color = "Chapman-Richards DBH", fill = "Chapman-Richards DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
    labs(x = NULL, y = "|height error residual|, m", color = NULL, fill = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$michaelisMenten))), alpha = 0.1, color = "grey25", shape = 16) +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$michaelisMenten)), color = "Michaelis-Menten GAM", fill = "Michaelis-Menten GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$michaelisMenten)), color = "Michaelis-Menten sqrt(DBH)", fill = "Michaelis-Menten sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$michaelisMenten)), color = "Michaelis-Menten DBH", fill = "Michaelis-Menten DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
    labs(x = NULL, y = NULL, color = NULL, fill = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaParton))), alpha = 0.1, color = "grey25", shape = 16) +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaParton)), color = "Sharma-Parton GAM", fill = "Sharma-Parton GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaParton)), color = "Sharma-Parton sqrt(DBH)", fill = "Sharma-Parton sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaParton)), color = "Sharma-Parton DBH", fill = "Sharma-Parton DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
    labs(x = "DBH, cm", y = "|height error residual|, m", color = NULL, fill = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaZhang))), alpha = 0.1, color = "grey25", shape = 16) +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaZhang)), color = "Sharma-Zhang GAM", fill = "Sharma-Zhang GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaZhang)), color = "Sharma-Zhang sqrt(DBH)", fill = "Sharma-Zhang sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
    geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightFromDiameter$sharmaZhang)), color = "Sharma-Zhang DBH", fill = "Sharma-Zhang DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
    labs(x = "DBH, cm", y = NULL, color = NULL, fill = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
  plot_layout(nrow = 2, ncol = 2)
  
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = alruHeightFromDiameter$ratkowsky$w * alruHeightFromDiameter$ratkowsky$working.residuals), alpha = 0.1, color = "grey25", shape = 16) +
    labs(x = "DBH, cm", y = "nlrob weight")
}


## red alder height-diameter GNLS regressions
#alruHeightFromDiameterGnls$chapmanRichards = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameter$chapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving with default gnlsControl(tolerance = 1E-6, nlsTol = 0.001, msTol = 1E-7), ok at 1E-6, 0.01, 1E-7
#alruHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameter$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02, >250 iterations at nlsTol = 0.05
#alruHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameter$sharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.5, > 250+50 iterations at nlsTol = 1 or with tighter tolerances than nlsTol = 0.2, msTol = 1E-5, tolerance = 1E-4
#alruHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameter$sharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.1, maxiter at default msTol and tolerance
#alruHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameter$sharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#alruHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, alru2016, start = alruHeightFromDiameter$sharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
#alruHeightFromDiameterGnls$weibull = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameter$weibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msVerbose = FALSE, returnObject = FALSE)) # step factor at nlsTol = 0.001
#alruHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
#save(alruHeightFromDiameterGnls, file = "trees/height-diameter/HtDia ALRU2 GNLS.rdata")

load("trees/height-diameter/HtDia ALRU2 GNLS.rdata")
alruHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", alruHeightFromDiameterGnls$chapmanRichards, alru2016)
alruHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", alruHeightFromDiameterGnls$chapmanRichardsBal, alru2016)
alruHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", alruHeightFromDiameterGnls$sharmaParton, alru2016)
alruHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", alruHeightFromDiameterGnls$sharmaPartonBal, alru2016)
alruHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", alruHeightFromDiameterGnls$sharmaZhang, alru2016)
alruHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", alruHeightFromDiameterGnls$sharmaZhangBal, alru2016)
alruHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", alruHeightFromDiameterGnls$weibull, alru2016)
alruHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", alruHeightFromDiameterGnls$weibullBal, alru2016)

alruHeightFromDiameterResultsGnls = bind_rows(lapply(alruHeightFromDiameterGnls, as_row)) %>%
  mutate(responseVariable = "height", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  alruHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameter$chapmanRichardsGnls), color = "Chapman-Richards GNLS", group = alru2016$isPlantation)) +
    annotate("text", x = 0, y = 50, label = "a) red alder, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 115), ylim = c(0, 50)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## red alder diameter-height regressions
#alruDiameterFromHeight$chapmanRichards = fit_gsl_nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, start = list(a1 = -18.3, b1 = 0.0305, b2 = 0.427), control = list(maxiter = 50))
#alruDiameterFromHeight$chapmanRichards = fit_gsl_nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.9999)), alru2016, start = list(a1 = -15.8, b1 = 0.0304, b2 = 0.362), control = list(maxiter = 50))
#alruDiameterFromHeight$chapmanRichards = fit_gsl_nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, algorithm = "port",
#                                            lower = list(a1 = -150, b1 = 0.02, b2 = 0.1), 
#                                            start = list(a1 = -50, b1 = 0.03, b2 = 0.7), 
#                                            upper = list(a1 = -1, b1 = 0.05, b2 = 2), control = list(warnOnly = TRUE))
#alruDiameterFromHeight$chapmanRichards = nls_multstart(DBH ~ a1 * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100,
#                                                      start_lower = list(a1 = -20, b1 = -0.03, b2 = -1), 
#                                                      start_upper = list(a1 = -10, b1 = 0.03, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeight$chapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeight$chapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeight$chapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b1p = -0.02, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b1p = 0.00, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b1p = 0.01, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b1p = 0.02, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b1p not significant
#alruDiameterFromHeight$chapmanRichardsMod = fit_gsl_nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-1, 0.7))
#alruDiameterFromHeight$chapmanRichardsMod = fit_gsl_nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2, 0.7))
#alruDiameterFromHeight$chapmanRichardsMod = fit_gsl_nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2.6, 0.7))
#alruDiameterFromHeight$chapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.1415926/180 * slope) + a4 * cos(3.1415926/180 * aspect) + a5 * sin(3.1415926/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999)), alru2016, iter = 100,
#                                                            start_lower = list(a1 = 10, a1p = 5, a2 = -0.01, a3 = -2, a4 = -2, a5 = -0.5, a6 = -0.1, b1 = -0.1, b1p = -0.1, b2 = -2), # 3.14159, if used, is seen as a mutable parameter
#                                                            start_upper = list(a1 = 20, a1p = 15, a2 = 0.01, a3 = 1, a4 = 2, a5 = 0.5, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight$chapmanRichardsRelHt = nls_multstart(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100, 
#                                                           lower = c(a1 = -100, a2 = -10, b1 = 1/60, b2 = 0),
#                                                           start_lower = list(a1 = -35, a2 = -1, b1 = 0.025, b2 = 0.4),
#                                                           start_upper = list(a1 = -25, a2 = 1, b1 = 0.035, b2 = 0.8), 
#                                                           upper = c(a1 = 0.1, a2 = 100, b1 = 1/25, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightKorf = fit_gsl_nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), lower = list(a1 = 0.5, b1 = 1, b2 = 0.3), algorithm = "port")
#alruDiameterFromHeightKorf = fit_gsl_nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = pmin(TotalHt^-1.3, 0.7), control = list(maxiter = 200, tol = 0.001, warnOnly = TRUE))
#alruDiameterFromHeightKorf = gnls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = varPower(0.65, ~TotalHt), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruDiameterFromHeight$linear = fit_lm(DBH ~ exp(1*(TotalHt - 1.37)^0.4) + 0, alru2016)
#alruDiameterFromHeight$linear = fit_lm(DBH ~ I(TotalHt - 1.37) + 0, alru2016)
#alruDiameterFromHeight$michaelisMenten = fit_gsl_nls(DBH ~ (a1*(TotalHt - 1.37)/(a2 - (TotalHt - 1.37)))^b1, alru2016, start = list(a1 = 25, a2 = 50, b1 = 1), weights = TotalHt^-2, control = gsl_nls_control(scale = "levenberg")) # inversion of Michaelis-Menten collapses to linear
#alruDiameterFromHeight$naslund = fit_gsl_nls(DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.5, a2 = -0.12), weights = TotalHt^-2)
#alruDiameterFromHeight$naslund = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 4.6, a1p = -1.5, a2 = -0.11, a2p = -0.014), weights = TotalHt^-2)
#alruDiameterFromHeight$power = fit_gsl_nls(DBH ~ a1*(TotalHt - 1.37)^b1, alru2016, algorithm = "port",
#                                  lower = list(a1 = 0, b1 = 1.2), # b1 < 1 is asymptotic in diameter rather than height
#                                  start = list(a1 = 0.4, b1 = 1.5), 
#                                  upper = list(a1 = 10, b2 = 3))
#alruDiameterFromHeight$ratkowsky = fit_gsl_nls(DBH ~ a1 / log(a2 * (TotalHt - 1.37)) - a3, alru2016, start = list(a1 = 1, a2 = 1, a3 = 0)) # not numerically stable due to log underruns
#alruDiameterFromHeight$ruark = fit_gsl_nls(DBH ~ a1*(TotalHt - 1.37)^b1*exp(b2*(TotalHt - 1.37)), alru2016, start = list(a1 = 1.29, b1 = 1.31, b2 = -0.027))
#alruDiameterFromHeight$schnute = fit_gsl_nls(DBH ~ a1*log(1 - b1*(TotalHt^b2 - 1.37^b2) / (40^b2 - 1.37^b2)), alru2016, start = list(a1 = -150, b1 = 0.5, b2 = 1.5)) # diverges to concave up
#alruDiameterFromHeight$sharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, alru2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = 0.001, b2 = 0.1, b3 = 0.1), 
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 0.1, b2 = 1, b3 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight$sibbesen = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), alru2016, iter = 100, 
#                                               lower = c(a1 = 0.01, b1 = 0.2, b2 = 0.2),
#                                               start_lower = list(a1 = 2, b1 = 0.25, b2 = 0.25),
#                                               start_upper = list(a1 = 4, b1 = 0.30, b2 = 0.30), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight = fit_gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53)) # step factor
#alruDiameterFromHeight$weibull = fit_gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -1000, b1 = 0.002, b2 = 0.95)) # NaN-inf
#alruDiameterFromHeight$weibull = fit_gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, algorithm = "port",
#                                    lower = list(a1 = -1000, b1 = 0.001, b2 = 0.4), # b1 constraint prevents NaN-inf
#                                    start = list(a1 = -150, b1 = 0.02, b2 = 0.75)) # collapses to linear fit
#alruDiameterFromHeight$weibullForm = fit_gsl_nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53), control = list(maxiter = 500))
#alruDiameterFromHeight$wykoff = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                             lower = c(a1 = 1, b1 = 1/55, b2 = 0.2),
#                                             start_lower = list(a1 = 3, b1 = 0.2, b2 = 0.5), 
#                                             start_upper = list(a1 = 7, b1 = 0.4, b2 = 0.7), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight$wykoffAat = nls_multstart(DBH ~ (a1 + a2 * tallerApproxBasalArea) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a2 = -1, b1 = -0.05, b2 = -1), 
#                                                start_upper = list(a1 = 1, a2 = 1, b1 = 0.05, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight$wykoffAat = nls_multstart(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a1p = -100, a2 = -1, a2p = -1, a3 = -1, a3p = -1, b1 = -1, b2 = -1, b2p = -1), 
#                                                start_upper = list(a1 = 1, a1p = 1, a2 = 1, a2p = 1, a3 = 1, a3p = 1, b1 = 1, b2 = 2, b2p = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
# coefficients for physiologically plausible curvature but poor accuracy
#alruDiameterFromHeight$chapmanRichards = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -16.9, a1p = 2.676, b1 = 0.0307, b2 = 0.0304, b2p = 0.0723), control = list(maxiter = 210)) # b1p not significant, >200 iterations to converge to start point
#alruDiameterFromHeight$chapmanRichardsAat = fit_gsl_nls(DBH ~ (a1 + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -15.4, a2 = 0.234, a3 = 0.233, b1 = 0.0342, b2 = 0.263, b2p = 0.164), control = list(maxiter = 500)) # a1p, a2, a2p, a3, a3p not significant
#alruDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -10.5, a1p = 2.93, a2 = 0.0030, a3 = -18.3, a4 = -0.095, a5 = -0.42, a6 = 0.137, b1 = 0.0307, b2 = 0.342, b2p = 0.0683), control = list(maxiter = 200)) # a2, a4, a5, a6, b1p not significant
#alruDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -36.7, a1p = 30.7, a2 = 33.2, a2p = -40.8, b1 = 0.0307, b2 = 0.562, b2p = -0.360), control = list(maxiter = 200))
# coefficients for physically implausible accuracy
alruDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22))) # a1p not significant
alruDiameterFromHeight$chapmanFormAat = fit_nlrob("Chapman-Richards form AA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0)) # a1p not significant, all NaN-inf with a3 * standBasalAreaApprox
alruDiameterFromHeight$chapmanFormBal = fit_nlrob("Chapman-Richards form BA+L", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -32.1, a1p = 4.8, a2 = 1.72, a2p = -0.45, a3 = -1.46, a3p = 0.45, b1 = -0.040, b2 = 1.31)) # b2p not significant
alruDiameterFromHeight$chapmanFormBalRelHt = fit_nlrob("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a9 = 29.8, a9p = -14.8, b1 = -0.0302, b2 = 1.30)) # NaN-inf
alruDiameterFromHeight$chapmanFormRelHt = fit_nlrob("Chapman-Richards form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -54.2, a9 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19)) # a1p not significant
alruDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50)) # b1p not significant
alruDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), control = list(maxiter = 50)) # a2p not significant
alruDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a5 = 7.37, a6 = 0.41, a7 = 0.52, b1 = -0.031, b2 = 2.47, b2p = -1.26)) # a4, a8 not significant
alruDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016, start = list(a1 = 33.0, a1p = -7.77, a9 = -0.95, b1 = -0.043, b2 = 1.39)) # a2p, b1p not significant
alruDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = alru2016gamConstraint), data = alru2016)
alruDiameterFromHeight$gamAat = fit_gam("REML GAM AA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 14, pc = alru2016gamConstraint), data = alru2016)
alruDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM AA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 496, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 6) # with all predictors k = 496, edf < 220, AIC 24509: 24540 without AAT, 24578 without ABA, 24571 without elevation, 24555 without slope, 24563 without sin(aspect), 24576 without cos(aspect), 24556 without topographic shelter
alruDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4)
alruDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), alru2016)
alruDiameterFromHeight$michaelisMentenForm = fit_gsl_nls("Michaelis-Menten form", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), alru2016, start = list(a1 = 100, a2 = 100, b1 = 1), control = gsl_nls_control(maxiter = 80)) # collapses to linear with or without b1p
alruDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.8, a1p = -1.0, a2 = -0.12, a2p = -0.006))
alruDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), alru2016)
alruDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15))
alruDiameterFromHeight$powerAat = fit_nlrob("power AA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147))
alruDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^b1, alru2016physio, start = list(a1 = 2.38, a1p = -0.70, a5 = 1.26, a6 = 0.078, a7 = 0.0075, b1 = 0.88)) # a4, a8, b1p not significant
alruDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.55, a1p = -1.73, a9 = -0.98, a9p = 0.65, b1 = 0.86, b1p = 0.14))
alruDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.35, b1 = 1.37, b1p = -0.24, b2 = -0.033, b2p = 0.022)) # a1p not significant
alruDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), alru2016, start = list(a1 = 0.00006, a2 = 0.04, b1 = 0.94, Ha = 100)) # poorly conditioned, singular gradient with Levenberg, nls(), or fit_nlrob()
alruDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, alru2016, start = list(a1 = 29, b1 = 0.71, b2 = 0.0001, b3 = -0.65, b4 = 0.21), control = nls.control(maxiter = 500)) # nls() NaN-infinity even with parameters from nls_multstart(), fit_nlrob() NaN-inf
alruDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119))
alruDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12))  # a2 not significant
alruDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a5 = 0.231, a6 = 0.191, a7 = 0.018, b1 = 2.13, b2 = -0.163, b2p = 0.0197)) # a4, a8 not significant
alruDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.712, a1p = -0.341, a9 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348))
alruDiameterFromHeight$weibull = fit_gsl_nls("Weibull", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -1000, b1 = 0.002, b2 = 0.95), control = gsl_nls_control(scale = "levenberg")) # collapses to linear, singular gradient with a1p, b1p, b2p, fit_nlrob() NaN-inf
#confint_nlrob(alruDiameterFromHeight$sibbesenFormRelHt, level = 0.99)

alruDiameterFromHeightResults = bind_rows(lapply(alruDiameterFromHeight, as_row)) %>% 
  mutate(responseVariable = "DBH", species = "ALRU2", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(alruDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_smooth(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanForm), y = alru2016$TotalHt, color = "Chapman-Richards form", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichards), y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsAat), y = alru2016$TotalHt, color = "Chapman-Richards AA+T", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsPhysio), y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsRelHt), y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeightKorf), y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$linear), y = alru2016$TotalHt, color = "linear", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$michaelisMentenForm), y = alru2016$TotalHt, color = "Michaelis-Menten form", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$naslund), y = alru2016$TotalHt, color = "Näslund", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$parabolic), y = alru2016$TotalHt, color = "parabolic", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$power), y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$ruark), y = alru2016$TotalHt, color = "Ruark", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$schnute), y = alru2016$TotalHt, color = "Schnute", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$sibbesenForm), y = alru2016$TotalHt, color = "Sibbesen form", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$weibull), y = alru2016$TotalHt, color = "Weibull", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -30*log(1 - pmin((0.03*(alru2016$TotalHt - 1.37))^0.85, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards")) +
    #geom_line(aes(x = (-16.95 + 2.676 * alru2016$isPlantation)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^(0.304 + 0.0723 * alru2016$isPlantation), 0.9999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99999", group = alru2016$isPlantation)) +
    #geom_line(aes(x = (-15.886)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^0.304, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.999")) +
    #geom_line(aes(x = (-18.329)*log(1 - pmin((0.0305*(alru2016$TotalHt - 1.37))^0.427, 0.99)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99")) +
    #geom_line(aes(x = -30*log(1 - pmin(0.03*(alru2016$TotalHt - 1.37)^1.0, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
    #geom_line(aes(x = -150*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.79, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
    #geom_line(aes(x = -119*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.86, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 2")) +
    #geom_line(aes(x = -108*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.89, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 3")) +
    #geom_line(aes(x = 5*(exp(0.33*(alru2016$TotalHt - 1.37)^0.60) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
    #geom_line(aes(x = 39*(exp(0.1*(alru2016$TotalHt - 1.37)^0.66) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
    #geom_line(aes(x = 125.8*(exp(0.0182*(alru2016$TotalHt - 1.37)^0.880) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
    #geom_line(aes(x = -33 + 3*alru2016$relativeHeight)*log(1 - pmin((0.03*(alru2016$TotalHt - 1.37))^0.85, 0.999), y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation)) +
    #geom_line(aes(x = (0.016*(alru2016$TotalHt - 1.37))^0.027 / (1 - (0.016*(alru2016$TotalHt - 1.37))^0.027), y = alru2016$TotalHt, color = "Curtis")) +
    #geom_line(aes(x = 1*exp(0.7*(alru2016$TotalHt - 1.37)^0.5), y = alru2016$TotalHt, color = "Korf")) +
    #geom_line(aes(x = -1/0.08*log((40 - alru2016$TotalHt)/alru2016$TotalHt) + 30, y = alru2016$TotalHt, color = "logistic")) +
    #geom_line(aes(x = 100*(alru2016$TotalHt - 1.37)/(100 - (alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Michaelis-Menten")) +
    #geom_line(aes(x = 2.5*sqrt(alru2016$TotalHt - 1.37)/(1 - 0.13*sqrt(alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Näslund")) +
    #geom_line(aes(x = (25*(alru2016$TotalHt - 1.37)/(50 - (alru2016$TotalHt - 1.37)))^1, y = alru2016$TotalHt, color = "Michaelis-Menten")) +
    #geom_line(aes(x = 0.3*(alru2016$TotalHt - 1.37)^1.5, y = alru2016$TotalHt, color = "power")) +
    #geom_line(aes(x = 1*(alru2016$TotalHt - 1.37)^1.2, y = alru2016$TotalHt, color = "power")) +
    #geom_line(aes(x = pmin(pmax(-22.3/log(1/33.6*(alru2016$TotalHt - 1.37)) - 5.3, 0), 100), y = alru2016$TotalHt, color = "Ratkowsky")) +
    #geom_line(aes(x = 32.6 - 0.33/log(0.32*(alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Ratkowsky")) +
    #geom_line(aes(x = 30*(alru2016$TotalHt - 1.37) / (60 - (alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Ratkowsky form")) +
    #geom_line(aes(x = 2*(alru2016$TotalHt - 1.37)^0.5 / (1 - 0.12*(alru2016$TotalHt - 1.37)^0.5), y = alru2016$TotalHt, color = "Ratkowsky")) +
    #geom_line(aes(x = 100*(alru20 16$TotalHt - 1.37)^0.01*(exp(0.01*(alru2016$TotalHt - 1.37)) - 1), y = alru2016$TotalHt, color = "Ruark form")) +
    #geom_line(aes(x = -1/0.006*log(1 - (1 - exp(-0.5))*(alru2016$TotalHt^1.5 - 1.3^1.5) / (40^1.5 - 1.3^1.5)), y = alru2016$TotalHt, color = "Schnute")) +
    #geom_line(aes(x = 3*(alru2016$TotalHt - 1.37)^(0.30*(alru2016$TotalHt - 1.37)^0.30), y = alru2016$TotalHt, color = "Sibbesen form")) +
    #geom_line(aes(x = 5.3*(alru2016$TotalHt - 1.37)^(0.35*(alru2016$TotalHt - 1.37)^0.20), y = alru2016$TotalHt, color = "Sibbesen form")) +
    #geom_line(aes(x = (-250*log(1 - pmin(0.013*(alru2016$TotalHt - 1.37), 0.9999)))^0.75, y = alru2016$TotalHt, color = "Weibull (Yang))")) +
    #geom_line(aes(x = 20*(alru2016$TotalHt - 1.37)^0.5*(exp(0.003*(alru2016$tph/alru2016$topHeight)^0.4*(alru2016$TotalHt - 1.37)) - 1)^0.9, y = alru2016$TotalHt, color = "modified Sharma-Parton", group = alru2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = (-204*log(1 - pmin(0.011 * (alru2016$TotalHt - 1.37), 0.9999)))^0.869, y = alru2016$TotalHt, color = "Weibull", group = alru2016$isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 50, label = "red alder, diameter from height", hjust = 0, size = 3.5) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot(alru2016) +
    geom_point(aes(x = TotalHt, y = DBH), alpha = 0.10, color = "grey25", shape = 16) +
    geom_line(aes(x = TotalHt, y = -50*log(1 - pmin(0.004*(TotalHt - 1.37)^1.5, 0.9)), color = "Chapman-Richards")) +
    annotate("text", x = 0, y = 90, label = "red alder, diameter from height", hjust = 0, size = 3.5) +
    labs(x = "height, m", y = "DBH, cm", color = NULL) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot(alru2016) +
    geom_bin_2d(aes(x = DBH, y = TotalHt), binwidth = c(2.5, 1)) +
    #geom_line(aes(x = 1.246*exp(1*(TotalHt - 1.37)^0.4), y = TotalHt, color = "linear")) +
    #geom_line(aes(x = predict(alruDiameterFromHeightKorf), y = TotalHt, color = "Korf")) +
    geom_line(aes(x = 1.246*exp(1*(TotalHt - 1.37)^0.4), y = TotalHt, color = "Korf manual")) +
    geom_line(aes(x = 0.155*exp(2.708*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf unweighted")) +
    geom_line(aes(x = 0.146*exp(2.732*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf weighted")) +
    geom_line(aes(x = 0.019*exp(4.691*(TotalHt - 1.37)^0.15), y = TotalHt, color = "Korf weighted to step factor")) +
    geom_line(aes(x = 0.051*exp(3.724*(TotalHt - 1.37)^0.18), y = TotalHt, color = "Korf GNLS")) +
    labs(x = "DBH, cm", y = "height, m", color = "regression\nform", fill = "stems\nmeasured") +
    scale_fill_viridis_c(trans = "log10")
  
  ggplot() +
    geom_point(aes(x = alru2016$TotalHt, y = residuals(alruDiameterFromHeight$ruark)), alpha = 0.1, color = "grey25", shape = 16) +
    geom_smooth(aes(x = alru2016$TotalHt, y = residuals(alruDiameterFromHeight$ruark)), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    #geom_point(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * residuals(alruDiameterFromHeightKorf)), alpha = 0.1, color = "grey25", shape = 16) +
    #geom_smooth(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * residuals(alruDiameterFromHeightKorf)), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    labs(x = "height, m", y = "DBH error residual, cm")
  
  ggplot(alru2016) + 
    geom_histogram(aes(x = DBH, y = 100 * ..count../sum(..count..)), binwidth = 2.5) +
    labs(x = "DBH, cm", y = "percentage of stems measured") +
  ggplot(alru2016) + 
    geom_histogram(aes(x = 100 * ..count../sum(..count..), y = TotalHt), binwidth = 1) +
    labs(x = "percentage of stems measured", y = "height, m")
  
  ggplot(alru2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = DBH, y = 40 / (1 + exp(-0.08*(DBH - 30))), color = "logistic")) +
    geom_line(aes(x = DBH, y = 1.37 + (40 - 1.37) / (1 + exp(-0.08*(DBH - 30)))^1.1, color = "generalized logistic")) +
    #geom_line(aes(x = seq(0, 100, by = 0.1), y = 1.37 + 40 * (1 - exp(-0.01*seq(0, 100, by = 0.1))^1.5), color = "Korf")) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}


## collect model parameters
alruParameters = bind_rows(bind_rows(bind_rows(lapply(alruHeightFromDiameter, get_coefficients)),
                                     bind_rows(lapply(alruHeightFromDiameterGnls, get_coefficients))) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(lapply(alruDiameterFromHeight, get_coefficients)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "ALRU2") %>%
  relocate(responseVariable, species, fitting, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)

## robust regression
if (includeInvestigatory)
{
  # adaptive weighting by fit_nlrob()
  ggplot() +
    geom_histogram(aes(x = alruHeightFromDiameter$chapmanRichards$rweights), binwidth = 0.01) +
    labs(x = "fit_nlrob() adaptive weight reduction", y = "number of red alder heights imputed") +
    scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), minor_breaks = c(0.5, 3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 3000, 4000, 6000, 7000, 8000, 9000), trans = scales::pseudo_log_trans(base = 10)) +
  ggplot() +
    geom_histogram(aes(x = alruDiameterFromHeight$chapmanRichards$rweights), binwidth = 0.01) +
    labs(x = "fit_nlrob() adaptive weight reduction", y = "number of red alder diameters imputed") +
    scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), minor_breaks = c(0.5, 3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 3000, 4000, 6000, 7000, 8000, 9000), trans = scales::pseudo_log_trans(base = 10))
}


## basal area from height
if (includeInvestigatory)
{
  #alruBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), alru2016, start = list(a1 = 500, b1 = 0.0002, b2 = 1.9), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # fit_nlrob() step factor
  alruBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), alru2016, start = list(a1 = 682, b1 = 0.00004, b2 = 1.9, b2p = -0.11), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # a1p, b1p not significant, fit_nlrob() step factor
  alruBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^b1, alru2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.00006, b1 = 1.91), weights = pmin(1/basalArea, 1E4)) # b1p not significant
  #confint2(alruBasalAreaFromHeightKorf, level = 0.99)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(alruBasalAreaFromHeightKorf), 100^2 * mean(-residuals(alruBasalAreaFromHeightKorf)), mean(abs(residuals(alruBasalAreaFromHeightKorf))), 1 - sum(residuals(alruBasalAreaFromHeightKorf)^2) / sum((alru2016$basalArea - mean(alru2016$basalArea)^2)),
          "power", AIC(alruBasalAreaFromHeightPower), 100^2 * mean(-residuals(alruBasalAreaFromHeightPower)), mean(abs(residuals(alruBasalAreaFromHeightPower))), 1 - sum(residuals(alruBasalAreaFromHeightPower)^2) / sum((alru2016$basalArea - mean(alru2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(alru2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(alruBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(alruBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "red alder height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}