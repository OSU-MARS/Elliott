# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western redcedar height-diameter regression form sweep
#thplHeightFromDiameter$gamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint, select = TRUE, weights = dbhWeight) # bs = "ts" -> 367, gamma = 2 -> 367, k = 169 min vs 367 default, method = "REML" -> 367
#thplHeightFromDiameter$sharmaPartonBalPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b1p = -0.0143, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = thplHeightFromDiameterWeights)
thpl2016 = trees2016 %>% filter(Species == "RC", isLiveUnbroken, is.na(TotalHt) == FALSE) %>% # live western redcedars measured for height
  mutate(dbhWeight = pmin(TreeCount/(0.14*DBH^1.20), 5*TreeCount),
         heightWeight = pmin(TreeCount/(2.29*(TotalHt - 1.37)^1.45), 5*TreeCount))
thpl2016physio = thpl2016 %>% filter(is.na(elevation) == FALSE)
thpl2016gamConstraint = c(DBH = -1.2264/0.5099, TotalHt = 1.37, standBasalAreaPerHectare = median(thpl2016$standBasalAreaPerHectare), basalAreaLarger = median(thpl2016$basalAreaLarger), standBasalAreaApprox = median(thpl2016$standBasalAreaApprox), tallerApproxBasalArea = median(thpl2016$tallerApproxBasalArea), elevation = median(thpl2016physio$elevation), slope = median(thpl2016physio$slope), aspect = median(thpl2016physio$aspect), topographicShelterIndex = median(thpl2016physio$topographicShelterIndex), relativeHeight = median(thpl2016$relativeHeight), relativeDiameter = median(thpl2016$relativeDiameter)) # point constraint for mgcv::s()

thpl2016defaultWeight = thpl2016 %>% mutate(dbhWeight = pmin(TreeCount/DBH, 5*TreeCount),
                                            heightWeight = pmin(TreeCount/TotalHt, 5*TreeCount))
thpl2016defaultWeightPhysio = thpl2016defaultWeight %>% filter(is.na(elevation) == FALSE)

thplOptions = tibble(fitHeight = TRUE, 
                     fitHeightNlrob = fitHeight,
                     fitHeightGnls = FALSE,
                     fitHeightMixed = fitHeight,
                     fitDbh = TRUE,
                     fitDbhNlrob = fitDbh,
                     fitDbhMixed = fitDbh)

if (thplOptions$fitHeight)
{
  thplHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH, thpl2016)) # isPlantation*DBH not significant (p = 0.044)
  thplHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2), thpl2016) # isPlantation*DBH not quite significant (p = 0.106), isPlantation*DBH^2 not significant
  
  thplHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131)) # a1p, b1p, b2p not significant
  thplHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1)) # a3, a3p, b1p, b2p not significant
  thplHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), significant = FALSE) # a2, a3, a4, a5, a6, a7, a8p, b2p not significant
  thplHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 58, a1p = -16, a2 = 0, a2p = 0.4, a8 = 0.3, a10 = -1.3, b1 = -0.012, b1p = -0.003, b2 = 1.13), significant = FALSE) # a2, a10, a10p not significant
  thplHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 61, a1p = -9, a2 = -0.1, a2p = 0.6, a10 = -1.3, b1 = -0.012, b2 = 1.1), significant = FALSE) # a2, a10, a10p not significant
  thplHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 7, a1p = 5, a2 = 0.2, a2p = 0.24, a9 = 47, a9p = -27, b1 = -0.021, b2 = 0.8, b2p = 0.2)) # a2, a3, a3p, b1p not significant, job step factor with nlrob()
  thplHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15)) # a4, a5, a6, a7, a8p, b2p not significant
  thplHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 72, a10 = -3.2, b1 = -0.012, b2 = 1.09)) # a10p not significant
  thplHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), significant = FALSE) # a10, a10p not significant
  thplHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016, start = list(a1 = 0.560, b1 = 0.069)) # a1p, b1p not significant
  thplHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176)) # b2p not significant
  thplHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), thpl2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175)) # a1p, b1p, b2p not significant
  thplHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176)) # b1p not significant
  thplHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649)) # a2p, a3p not significant
  thplHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^b1, thpl2016, start = list(a1 = 0.542, b1 = 0.939)) # a1p, b1p not significant
  thplHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
  thplHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), thpl2016, start = list(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008)) # dp not significant, susceptible to NaN-inf
  thplHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09)) # a1p, b3p, b4p not significant
  thplHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 38, b1 = 0.1, b2 = -0.013, b3 = -0.1, b4 = 1.03)) # a1p, b1p, b2p, b3p, b4p not significant
  thplHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10)) # b1, b1p, a4, a5, a6, a7, b3p, b4p not significant
  thplHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 25, a1p = -6, a8 = 0.12, a10 = -0.7, b1 = 0.21, b2 = -0.008, b2p = -0.011, b3 = -0.01, b4 = 1.12), significant = FALSE) # a10, a10p not significant
  thplHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 39, a10 = -1.7, b1 = 0.12, b2 = -0.01, b3 = 0, b4 = 1.07), significant = FALSE) # a10, a10p not significant
  thplHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09)) # a4, a5, a5, a6, a7, b1p, b3p, b4p not significant
  thplHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = list(a1 = 21, a10 = 0, b1 = 0.25, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = 0, b4 = 1.12), significant = FALSE) # a10, a10p not significant
  thplHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 36, a8 = 0.18, a10 = -2, b1 = 0.13, b2 = -0.01, b3 = -0.03, b4 = 1.09), significant = FALSE) # a1p, a10, a10p, b2p not significant
  thplHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016, start = list(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097)) # b1, b1p, b2p, b3p not significant
  thplHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016, start = list(a1 = 45, a1p = -7, a2 = -0.1, a2p = 0.4, b1 = -0.05, b2 = -0.02, b3 = -0.078, b4 = 1.08)) # a2, b1, b1p, b3, b3p, b4p not significant
  thplHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), thpl2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078)) # a1p, b1p, b2p not significant
  thplHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141)) # b2p not significant
  thplHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131)) # a1p, a2, a3, b1p, b2p not significant
  thplHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778)) # a1p, a2, a3, a3p, b1p, b2p not significant
  
  if (thplOptions$fitHeightNlrob)
  {
    thplHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131)))
    thplHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1))
    thplHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # job step factor
    thplHeightFromDiameterNlrob$chapmanRichardsBalPhysioRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 58, a1p = -14, a2 = -0.11, a2p = 0.5, a8 = 0.3, a10 = -1.8, b1 = -0.012, b1p = -0.003, b2 = 1.14), significant = FALSE)
    thplHeightFromDiameterNlrob$chapmanRichardsBalRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 61, a1p = -8, a2 = -0.13, a2p = 0.6, a10 = -1.4, b1 = -0.012, b2 = 1.12), significant = FALSE)
    thplHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 0, a1p = 17, a2 = 0, a2p = 0.25, a3 = 0.02, a9 = 38, a9p = -28, b1 = -0.023, b2 = 0.4, b2p = 0.9), control = nls.control(tol = 0.01)) # job step factor
    thplHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15))
    thplHeightFromDiameterNlrob$chapmanRichardsRelDbh = fit_nlrob("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 66, a10 = -3.2, b1 = -0.010, b2 = 1.09))
    thplHeightFromDiameterNlrob$chapmanRichardsRelDbhPhysio = fit_nlrob("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), control = nls.control(tol = 1E-4), significant = FALSE)
    thplHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016, start = list(a1 = 0.560, b1 = 0.069))
    thplHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176))
    thplHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), thpl2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175))
    thplHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176))
    thplHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649))
    thplHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, thpl2016, start = list(a1 = 0.542, b1 = 0.939))
    thplHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
    thplHeightFromDiameterNlrob$richardsW = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), thpl2016, start = list(Ha = 43, Hap = -10, d = 0.9, kU = 0.012, kUp = 0.004), control = nls.control(tol = 0.001)) # job step factor
    thplHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09), control = nls.control(tol = 0.001)) # job step factor
    thplHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 44, b1 = 0.07, b2 = -0.013, b3 = -0.10, b4 = 1.03), control = nls.control(maxiter = 100, tol = 0.001)) # job step factor
    thplHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10), control = nls.control(tol = 1E-4)) # job step factor
    thplHeightFromDiameterNlrob$sharmaPartonBalPhysioRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 25, a1p = -8, a8 = 0.13, a10 = -0.9, b1 = 0.18, b2 = -0.011, b2p = -0.010, b3 = 0, b4 = 1.13), control = nls.control(tol = 0.001), significant = FALSE) # job step factor
    thplHeightFromDiameterNlrob$sharmaPartonBalRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 50, a10 = -3, b1 = 0.12, b2 = -0.01, b3 = 0, b4 = 1.07), control = nls.control(maxiter = 100, tol = 0.001), significant = FALSE) # step factor
    thplHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 28, a1p = -10, a8 = 0.13, b1 = 0.16, b2 = -0.011, b2p = -0.01, b3 = 0, b4 = 1.1), control = nls.control(tol = 0.01)) # b3 not significant, job step factor
    thplHeightFromDiameterNlrob$sharmaPartonRelDbh = fit_nlrob("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = list(a1 = 27, a10 = -0.7, b1 = 0.22, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = -0.003, b4 = 1.12), control = nls.control(tol = 0.001), significant = FALSE)
    thplHeightFromDiameterNlrob$sharmaPartonRelDbhPhysio = fit_nlrob("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 41, a8 = 0.2, a10 = -2, b1 = 0.13, b2 = -0.01, b3 = 0, b4 = 1.09), significant = FALSE)
    thplHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016, start = list(a1 = 36, a1p = -3.0, b1 = 0.1, b2 = -0.02, b3 = 0, b4 = 1.2, b4p = -0.2)) # b3 not significant
    thplHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016, start = list(a1 = 44, a1p = -7, a2 = -0.12, a2p = 0.45, b1 = 0.05, b2 = -0.017, b3 = -0.02, b4 = 1.1), control = nls.control(maxiter = 100, tol = 0.001)) # b3 not significant, job step factor
    thplHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), thpl2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078))
    thplHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141), control = nls.control(maxiter = 100, tol = 1E-4)) # job step factor
    thplHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131))
    thplHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778))
    #lapply(thplHeightFromDiameterNlrob$sharmaPartonPhysio$fit, confint_nlrob, level = 0.99)
  } else {
    thplHeightFromDiameterNlrob = list()
  }
  
  thplHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, thpl2016defaultWeight, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131)))
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, thpl2016defaultWeight, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1))
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016defaultWeightPhysio, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016defaultWeight, start = list(a1 = 7, a1p = 5, a2 = 0.2, a2p = 0.24, a3 = -0.03, a9 = 47, a9p = -27, b1 = -0.021, b2 = 0.8, b2p = 0.2))
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 58, a1p = -14, a2 = -0.07, a2p = 0.48, a8 = 0.27, a10 = -1.3, b1 = -0.012, b1p = -0.003, b2 = 1.13), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 61, a1p = -9, a2 = -0.1, a2p = 0.6, a10 = -1.4, b1 = -0.012, b2 = 1.11), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016defaultWeightPhysio, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15))
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, thpl2016defaultWeight, start = list(a1 = 74, a10 = -3.2, b1 = -0.011, b2 = 1.09))
  thplHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016defaultWeightPhysio, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016defaultWeight, start = list(a1 = 0.560, b1 = 0.069))
  thplHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016defaultWeight, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176))
  thplHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), thpl2016defaultWeight, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175))
  thplHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016defaultWeight, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176))
  thplHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016defaultWeight, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649))
  thplHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^b1, thpl2016defaultWeight, start = list(a1 = 0.542, b1 = 0.939))
  thplHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016defaultWeight, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
  thplHeightFromDiameterGslNlsDefault$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), thpl2016defaultWeight, start = list(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008))
  thplHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016defaultWeight, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09))
  thplHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016defaultWeight, start = list(a1 = 38, b1 = 0.12, b2 = -0.013, b3 = -0.1, b4 = 1.02))
  thplHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016defaultWeightPhysio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10))
  thplHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016defaultWeightPhysio, start = list(a1 = 23, a1p = -6, a8 = 0.12, a10 = -0.7, b1 = 0.21, b2 = -0.01, b2p = -0.010, b3 = -0.012, b4 = 1.14), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016defaultWeight, start = list(a1 = 35, a10 = -1.4, b1 = 0.15, b2 = -0.01, b3 = 0, b4 = 1.07), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016defaultWeightPhysio, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09))
  thplHeightFromDiameterGslNlsDefault$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016defaultWeight, start = list(a1 = 19, a10 = -0.3, b1 = 0.29, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = -0.03, b4 = 1.13), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016defaultWeightPhysio, start = list(a1 = 36, a8 = 0.18, a10 = 0, b1 = 0.2, b2 = -0.01, b3 = 0.03, b4 = 1.09), significant = FALSE)
  thplHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016defaultWeight, start = list(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097))
  thplHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016defaultWeight, start = list(a1 = 53.2, a1p = -8.857, a2 = -0.002, a2p = 0.10, b1 = -0.016, b2 = -0.025, b3 = -0.078, b4 = 1.126))
  thplHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), thpl2016defaultWeight, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078))
  thplHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016defaultWeight, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141))
  thplHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016defaultWeight, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131))
  thplHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), thpl2016defaultWeight, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778))
  
  thplHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint) # newton() step failure with family = scat, internal code errors with scat(theta = <fixed val>), see https://stats.stackexchange.com/questions/410515/how-different-are-restricted-cubic-splines-and-penalized-splines for discusson of thin plate versus other spline types
  thplHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint)
  thplHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 20, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint) # slope and elevation not supported, aspect not tested since insufficient data for full model
  thplHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 57, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint)
  thplHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 22, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint)
  thplHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 18, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint) # k reduces from 85 to 18 without aspect
  thplHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", TotalHt ~ s(DBH, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint)
  thplHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", TotalHt ~ s(DBH, elevation, slope, topographicShelterIndex, relativeDiameter, bs = "ts", k = 57, by = as.factor(isPlantation), pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint)

  save(file = "trees/height-diameter/data/THPL TotalHt.Rdata", thplHeightFromDiameter, thplHeightFromDiameterNlrob, thplHeightFromDiameterGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  print(thplHeightFromDiameterResults %>% select(-responseVariable, -species, -fixedWeight, -n, -power, -significant, -contains("NaturalRegen"), -contains("Plantation")), n = 30)
  ggplot() +
    geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = thpl2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = thpl2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$curtis), color = "Curtis", group = thpl2016$isPlantation)) +
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$gam), color = "GAM", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$korf), color = "Korf", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$linear), color = "linear", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$parabolic), color = "parabolic", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$power), color = "power", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$prodan), color = "Prodan", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$richardsW), color = "unified Richards", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sibbesen), color = "Sibbesen", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibull), color = "Weibull", group = thpl2016$isPlantation)) +
    annotate("text", x = 0, y = 65, label = "western redcedar, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 65)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
  
  # dbhClassSize = 50
  # errorByDbhClass = tibble(dbhClass = dbhClassSize*floor(thpl2016$DBH/dbhClassSize) + 0.5*dbhClassSize, fittedValue = predict(thplHeightFromDiameter$gam, thpl2016), height = thpl2016$TotalHt, residual = fittedValue - height) %>%
  #  #mutate(residual = residual - if_else(dbhClass == 50, -0.477/376, 0.107/95)) %>%
  #  group_by(dbhClass) %>%
  #  summarize(n = n(),
  #            totalHeight = sum(height),
  #            totalFitted = sum(fittedValue),
  #            meanBiasPerTree = sum(residual) / n,
  #            meanBiasPerTreePct = 100 * sum(residual/height) / n,
  #            minError = min(residual),
  #            meanError = mean(residual),
  #            maxError = max(residual),
  #            minPct = 100 * min(residual/height),
  #            meanPct = 100 * mean(residual/height),
  #            maxPct = 100 * max(residual/height),
  #            .groups = "drop") %>%
  #  filter(n >= 10)
  # errorByDbhClass
}


## western redcedar height-diameter GNLS regressions
if (thplOptions$fitHeightGnls)
{
  thplHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131), control = gnlsControl(nlsTol = 0.001))) # step halving at nlsTol = 1 with corSymm
  thplHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, thpl2016, start = thplHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving at nlsTol = 0.2 with corSymm
  thplHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = thplHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving at nlsTol = 0.2 with corSymm
  thplHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = thplHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving with plot correlation
  thplHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016, start = thplHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001)) # step halving with plot correlation
  thplHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016, start = thplHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving with plot correlation
  thplHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = thplHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # corSymm() viable but dropped
  thplHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = thplHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001)) # step halving at nlsTol = 1 with corSymm

  save(file = "trees/height-diameter/data/THPL TotalHt gnls.Rdata", thplHeightFromDiameterGnls)
}
if (htDiaOptions$includeInvestigatory)
{
  thplHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(thplHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(thplHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(thplHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  ggplot() +
    geom_point(aes(x = thpl2016natural$DBH, y = thpl2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = thpl2016natural$DBH, y = thpl2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "western redcedar naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = thpl2016plantation$DBH, y = thpl2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = thpl2016plantation$DBH, y = thpl2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "western redcedar plantation height, m")
  
  ggplot() +
    geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibullBal), color = "Weibull BA+L"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$power), color = "power")) +
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibull), color = "Weibull")) +
    annotate("text", x = 0, y = 85, label = "a) western redcedar, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


if (thplOptions$fitHeightMixed)
{
  thplHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1r)*(1 - exp(b1*DBH))^b2, thpl2016, 
                                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                start = list(fixed = c(a1 = 48.2, b1 = -0.015, b2 = 1.131)), control = nlmeControl(maxIter = 250)))
  thplHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, thpl2016, 
                                                            fixedFormula = a1 + a1p + a2 + a2p + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                            start = list(fixed = c(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1)), control = nlmeControl(maxIter = 500))
  thplHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, 
                                                                  fixedFormula = a1 + a1p + a2 + a2p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12)), control = nlmeControl(maxIter = 250), significant = FALSE)
  thplHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, 
                                                               fixedFormula = a1 + a1p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1, start = list(fixed = c(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15)))
  thplHeightFromDiameterMixed$curtis = fit_nlme("Curtis", TotalHt ~ 1.37 + (a1 + a1r) * DBH / (1 + DBH)^b1, thpl2016, 
                                                fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1, 
                                                start = list(fixed = c(a1 = 0.560, b1 = 0.069)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # max iterations in job
  thplHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176)), control = nlmeControl(maxIter = 250))
  thplHeightFromDiameterMixed$korf = fit_nlme("Korf", TotalHt ~ 1.37 + (a1 + a1r)*exp(b1*DBH^b2), thpl2016, 
                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                              start = list(fixed = c(a1 = 1825, b1 = -8.726, b2 = -0.175)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  thplHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, 
                                                         fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176)), control = nlmeControl(maxIter = 250)) # job >100 iterations
  thplHeightFromDiameterMixed$prodan = fit_nlme("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3 + a3r), thpl2016, 
                                                fixedFormula = a1 + a1p + a2 + a3 ~ 1, randomFormula = a3r ~ 1,
                                                start = list(fixed = c(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649)))
  thplHeightFromDiameterMixed$power = fit_nlme("power", TotalHt ~ 1.37 + (a1 + a1r)*DBH^b1, thpl2016, 
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 0.542, b1 = 0.939)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # job >500 iterations without relaxed tolerances
  thplHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, 
                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151)))
  thplHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), thpl2016, 
                                                   fixedFormula = Ha + Hap + d + kU + kUp ~ 1, randomFormula = Har ~ 1,
                                                   start = list(fixed = c(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008)))
  thplHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, 
                                                      fixedFormula = a1 + b1 + b1p + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                      start = list(fixed = c(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  thplHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, 
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 38, b1 = 0.1, b2 = -0.013, b3 = -0.1, b4 = 1.03)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve, step halving
  thplHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, 
                                                               fixedFormula = a1 + a1p + a8 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  thplHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, 
                                                            fixedFormula = a1 + a1p + a8 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 1E-3)) # singular precision matrix, step halving
  thplHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016, 
                                                     fixedFormula = a1 + a1p + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001))  # singularity in backsolve
  thplHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016, 
                                                        fixedFormula = a1 + a1p + a2 + a2p + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 45, a1p = -7, a2 = -0.1, a2p = 0.4, b1 = -0.05, b2 = -0.02, b3 = -0.078, b4 = 1.08)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #thplHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", TotalHt ~ 1.37 + a1*DBH^((b1 + b1r)*DBH^b2), thpl2016, 
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1,
  #                                                start = list(fixed = c(a1 = 0.302, b1 = 1.495, b2 = -0.078))) # a1r: step halving, singular precision
  thplHeightFromDiameterMixed$weibull = fit_nlme("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141)))
  thplHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016,
                                                    fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131)))

  thplHeightFromDiameterMixed$gamm = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8) + s(StandID, bs = "re"), data = thpl2016, mixed = TRUE)
  thplHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13) + s(StandID, bs = "re"), data = thpl2016, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/THPL TotalHt mixed.Rdata", thplHeightFromDiameterMixed)
}

 
## western redcedar diameter-height regressions
if (thplOptions$fitDbh)
{
  thplDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), thpl2016)) # isPlantation*(TotalHt - 1.37) not significant
  thplDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), thpl2016) # (TotalHt - 1.37)^2 not significant
  
  thplDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 200, b1 = 0.01, b2 = 0.95), control = gsl_nls_control(maxiter = 500, xtol = 1E-5)) # a1p, b1p, b2p not significant, a1-b1 parameter evaporation: singular gradient with nls(), no convergence from nls_multstart(), NaN-inf with nlrob()
  thplDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE) # NaN-inf with nls() and nlrob
  thplDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 200, a2 = -10, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 300), significant = FALSE) # step size with nls() and nlrob()
  thplDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * pmin(relativeHeight, 1.5)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 10, a2 = 0, a9 = 2.3, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250, xtol = 0.001), significant = FALSE) # a2, a3 not significant, a1-b1 parameter evaporation: nlrob() step factor with either a2 or a3
  thplDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8), control = gsl_nls_control(maxiter = 500)) # step size with nls(), >500 iterations with nlrob()
  thplDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -200, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250)) # a1p and b2p not significant, poor convergence with b1p, step factor with nlrob()
  thplDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE) # a1p, b1p not significant, step factor with nlrob()
  thplDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016physio, start = list(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 5E-5)) # no physiographic effects significant, a1-b1 parameter evaporation: step factor with nlrob()
  thplDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9), control = gsl_nls_control(maxiter = 500), significant = FALSE) # step factor with nlrob()
  thplDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), thpl2016, start = list(a1 = 519, a2 = 237, b1 = 1.00)) # a1p, a2p, b1p not significant, singular gradient with nlrob()
  thplDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024))
  thplDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08)) # no significant plantation effects
  #thplDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09)) # no significant plantation effects
  #thplDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE) # no significant physiographic effects
  #thplDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13)) # a1p and b1p not significant
  thplDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01)) # a1p, b1p, b2p not significant
  thplDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), significant = FALSE) # a2, a2p, a3, a3p, b1p, b2p not significant
  thplDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 2.9, a2 = -0.005, a4 = -0.001, b1 = 0.93, b2 = 0.006), significant = FALSE) # a2, a3 not significant, no AIC discrimination
  thplDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = -1, b1 = 0.9, b2 = 0), significant = FALSE) # a2, a3, a4, a9, b2 not significant, drop ABA on AIC
  thplDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.7, a3 = 0, a9 = 0, b1 = 0.95, b2 = 0.005), significant = FALSE) # a9, a9p, b2 not significant
  thplDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE) # a1p, a5, a6, a7, a8, b1p, b2p not significant
  thplDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE) # a9, a9p, b1p, b2, b2p not significant
  thplDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 3.2, a4 = 0, a9 = -1, b1 = 0.9, b2 = 0.01), significant = FALSE) # a4, a9 not significant
  #thplDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016, start = list(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30), control = gsl_nls_control(maxiter = 200)) # singular gradient with nlrob() and gsl_nls()
  thplDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, thpl2016, start = list(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1), control = gsl_nls_control(maxiter = 250, xtol = 0.025)) # a1-b2 evaporation, b1, b3 not significant, NaN-inf with nls() from nls_multstart() point, NaN-inf, singular gradient, or code syntax error with nlrob()
  thplDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12)) # no significant plantation effects
  thplDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029)) # no significant plantation effects
  thplDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1), significant = FALSE) # a2, a3, a8 not significant, drop ABA on AIC
  thplDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.3, a2 = 0, a8 = -0.017, a9 = 1.0, b1 = 0.7, b2 = 0), significant = FALSE) # a2, a3, a8, a9, b2 not significant, no a2-a3 AIC discrimination
  thplDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.5, a2 = 0, a9 = 0, a9p = 0, b1 = 0.6, b2 = 0.12), significant = FALSE) # a2, a9, a9p, b2 not significant
  thplDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE) # a1p, no physiographic effects significant
  thplDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  thplDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.6, a8 = -0.01, a9 = 0.7, b1 = 0.7, b2 = 0.1), significant = FALSE) # a9 not significant
  thplDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016, start = list(a1 = -300, b1 = 0.04, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 1E-4)) # a1p, b1p, b2p not significant, a1-b1 parameter evaporation: NaN-inf with nlrob()
  
  if (thplOptions$fitDbhNlrob)
  {
    thplDiameterFromHeightNlrob = list(naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024)))
    thplDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08))
    #thplDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09))
    #thplDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE)
    #thplDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13))
    thplDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01))
    thplDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # occasional job step factor
    thplDiameterFromHeightNlrob$ruarkAbatPhysio = fit_nlrob("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 2.9, a2 = 0, a4 = -0.001, b1 = 0.93, b2 = 0.006), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # step factor
    thplDiameterFromHeightNlrob$ruarkAbatPhysioRelHt = fit_nlrob("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = 2.3, b1 = 0.9, b2 = 0), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # job step factor
    thplDiameterFromHeightNlrob$ruarkAbatRelHt = fit_nlrob("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.4, a3 = 0, a9 = 3, b1 = 0.8, b2 = 0.01), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
    thplDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE)
    thplDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE)
    thplDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, start = list(a1 = 3.2, a4 = -0.001, a9 = 1, b1 = 0.8, b2 = 0), significant = FALSE) # a4, a9, b2 not significant
    thplDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12))
    thplDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.1, a2 = -0.004, b1 = 0.7, b2 = 0.1), control = nls.control(tol = 1E-4)) # job step factor
    thplDiameterFromHeightNlrob$sibbesenReplaceAbatPhysio = fit_nlrob("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1), significant = FALSE)
    thplDiameterFromHeightNlrob$sibbesenReplaceAbatPhysioRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.3, a2 = 0, a8 = -0.01, a9 = 0.5, b1 = 0.7, b2 = 0.1), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
    thplDiameterFromHeightNlrob$sibbesenReplaceAbatRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.1, a2 = 0, a9 = 0, a9p = 0.7, b1 = 0.6, b2 = 0.12), significant = FALSE)
    thplDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
    thplDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
    thplDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016, start = list(a1 = -250, b1 = 0.043, b2 = 0.58), control = nls.control(maxiter = 500))
    #confint_nlrob(thplDiameterFromHeight$sibbesenReplacePhysio, level = 0.99, weights = pmin(thpl2016$TotalHt^if_else(thpl2016$isPlantation, -1.7, -1.6), 0.5))
  } else {
    thplDiameterFromHeightNlrob = list()
  }
  
  thplDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016defaultWeight, start = list(a1 = 200, b1 = 0.01, b2 = 0.95), control = gsl_nls_control(maxiter = 250, xtol = 1E-5)))
  thplDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016defaultWeight, start = list(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016defaultWeight, start = list(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8), control = gsl_nls_control(maxiter = 500))
  thplDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016defaultWeight, start = list(a1 = -200, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250))
  thplDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016defaultWeight, start = list(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016defaultWeightPhysio, start = list(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 5E-5))
  thplDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016defaultWeight, start = list(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), thpl2016defaultWeight, start = list(a1 = 519, a2 = 237, b1 = 1.00))
  thplDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016defaultWeight, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024))
  thplDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016defaultWeight, start = list(a1 = 1.93, b1 = 1.08))
  #thplDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, thpl2016defaultWeight, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09))
  #thplDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016defaultWeightPhysio, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE)
  #thplDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016defaultWeight, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13))
  thplDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeight, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01))
  thplDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeight, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeightPhysio, start = list(a1 = 1.6, a2 = -0.01, a4 = -0.0006, b1 = 1.2, b2 = -0.009), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeightPhysio, start = list(a1 = 1.6, a3 = -0.003, a4 = -0.0006, a9 = 0.4, b1 = 1.27, b2 = -0.01), significant = FALSE) 
  thplDiameterFromHeightGslNlsDefault$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeight, start = list(a1 = 1.3, a3 = -0.003, a9 = 0.25, b1 = 1.3, b2 = -0.008), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeightPhysio, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeight, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016defaultWeightPhysio, start = list(a1 = 1.6, a4 = -0.0005, a9 = -0.4, b1 = 1.2, b2 = -0.01), significant = FALSE) # a4, a9 not significant
  #thplDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016defaultWeight, start = list(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30), control = gsl_nls_control(maxiter = 200))
  thplDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, thpl2016defaultWeight, start = list(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1), control = gsl_nls_control(maxiter = 250, xtol = 0.025))
  thplDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeight, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12))
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeight, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029))
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeightPhysio, start = list(a1 = 1.5, a2 = -0.009, a8 = -0.005, b1 = 1.2, b2 = -0.04), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeightPhysio, start = list(a1 = 1.41, a2 = -0.009, a8 = -0.005, a9 = 0, b1 = 1.4, b2 = -0.05), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeight, start = list(a1 = 1.5, a2 = -0.008, a9 = 0, a9p = 0, b1 = 1.4, b2 = 0), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeightPhysio, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeight, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  thplDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016defaultWeightPhysio, start = list(a1 = 1.4, a8 = 0, a9 = 0.3, b1 = 1.3, b2 = -0.035), significant = FALSE) # a8, a9 not significant
  thplDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016defaultWeight, start = list(a1 = -300, b1 = 0.04, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 1E-4))
  
  # individual term selection: TotalHt by = isPlantation only, AAT retained by AIC but not significant (p = 0.38)
  thplDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint) # newton() step failure with scat()
  thplDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint)
  thplDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint)
  thplDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", DBH ~ s(TotalHt, standBasalAreaApprox, topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 22, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint) # drop ABA and elevation on AIC
  thplDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint) # drop elevation and topographic shelter on AIC
  thplDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint)
  thplDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 57, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint) # drop elevation and aspect on AIC

  save(file = "trees/height-diameter/data/THPL DBH.Rdata", thplDiameterFromHeight, thplDiameterFromHeightNlrob, thplDiameterFromHeightGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  print(thplDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(thpl2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanReplaceAbat), y = TotalHt, color = "Chapman-Richards replace approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$schnute), y = TotalHt, color = "Schnute inverse", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen replace", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -100 * log(1 - pmin(0.015*(TotalHt - 1.37)^1.0, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
    #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.45, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -1/0.0003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute inverse"), alpha = 0.5) +
    geom_line(aes(x = 30*topHeight^0.5*(exp(0.01 * (tph/standBasalAreaPerHectare)^0.25*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "modified Sharma-Parton"), alpha = 0.5) +
    annotate("text", x = 0, y = 62, label = "western redcedar, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}

if (thplOptions$fitDbhMixed)
{
  thplDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", DBH ~ (a1 + a1r)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, 
                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                               start = list(fixed = c(a1 = 200, b1 = 0.01, b2 = 0.95)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001))) # singularity in backsolve, max iterations
  thplDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, 
                                                            fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                            start = list(fixed = c(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singularity in backsolve
  #thplDiameterFromHeightMixed$chapmanReplaceBal = fit_nlme("Chapman-Richards replace BA+L", DBH ~ (a1 + a1r + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, 
  #                                                         fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                         start = list(fixed = c(a1 = 200, a2 = -10, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 300, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving
  #thplDiameterFromHeightMixed$chapmanReplaceBalRelHt = fit_nlme("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1r + a2 * basalAreaLarger + a9 * pmin(relativeHeight, 1.5)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, 
  #                                                              fixedFormula = a1 + a2 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                              start = list(fixed = c(a1 = 10, a2 = 0, a9 = 2.3, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 250, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singularity in backsolve
  #thplDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", DBH ~ (a1 + a1r + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, 
  #                                                           fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                           start = list(fixed = c(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  thplDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", DBH ~ (a1 + a1r)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = -200, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  thplDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                             start = list(fixed = c(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations, step halving
  #thplDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", DBH ~ (a1 + a1r + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016physio, 
  #                                                             fixedFormula = a1 + a1p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                             start = list(fixed = c(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations, step halving
  thplDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                              start = list(fixed = c(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving, singularity in backsolve
  thplDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), thpl2016, 
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1, 
                                                                start = list(fixed = c(a1 = 519, a2 = 237, b1 = 1.00)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations, step halving
  thplDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", DBH ~ (a1 + a1r + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1, 
                                                 start = list(fixed = c(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024)))
  thplDiameterFromHeightMixed$power = fit_nlme("power", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1, thpl2016, 
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1, 
                                               start = list(fixed = c(a1 = 1.93, b1 = 1.08)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # job max iterations
  #thplDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, thpl2016, 
  #                                                 fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                 start = list(fixed = c(a1 = 1.94, a2 = -0.00051, b1 = 1.09)))
  #thplDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, 
  #                                                   fixedFormula = a1 + a8 + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                   start = list(fixed = c(a1 = 2.26, a8 = -0.0060, b1 = 1.08)), significant = FALSE)
  #thplDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", DBH ~ (a1 + a1r + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, 
  #                                                  fixedFormula = a1 + a9 + a9p + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                  start = list(fixed = c(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13)))
  thplDiameterFromHeightMixed$ruark = fit_nlme("Ruark", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, 
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                               start = list(fixed = c(a1 = 2.8, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations
  thplDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, 
                                                   fixedFormula = a1 + a3 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                   start = list(fixed = c(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # job max iterations
  thplDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, 
                                                         fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = 2.9, a2 = -0.005, a4 = -0.001, b1 = 0.93, b2 = 0.006)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # job max iterations
  #thplDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, 
  #                                                            fixedFormula = a1 + a3 + a4 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                            start = list(fixed = c(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = -1, b1 = 0.9, b2 = 0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations
  thplDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, 
                                                        fixedFormula = a1 + a3 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                        start = list(fixed = c(a1 = 2.7, a3 = 0, a9 = 0, b1 = 0.95, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  thplDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", DBH ~ (a1 + a1r + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, 
                                                     fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                     start = list(fixed = c(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  #thplDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016, 
  #                                                  fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                  start = list(fixed = c(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  thplDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", DBH ~ (a1 + a1r + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), thpl2016physio, 
                                                          fixedFormula = a1 + a4 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                          start = list(fixed = c(a1 = 3.2, a4 = 0, a9 = -1, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations, false convergence
  #thplDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1)), thpl2016, 
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1, 
  #                                               start = list(fixed = c(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  #thplDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, thpl2016, 
  #                                                    fixedFormula = a1 + b1 + b2 + b4 ~ 1, randomFormula = a1r ~ 1, 
  #                                                    start = list(fixed = c(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1)), control = nlmeControl(maxIter = 250, tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # singularity in backsolve
  thplDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = 3.4, b1 = 0.8, b2 = 0.12)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations
  thplDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                             start = list(fixed = c(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5))
  thplDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, 
                                                                   fixedFormula = a1 + a2 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                   start = list(fixed = c(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  #thplDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, 
  #                                                                      fixedFormula = a1 + a2 + a8 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                                      start = list(fixed = c(a1 = 3.3, a2 = 0, a8 = -0.017, a9 = 1.0, b1 = 0.7, b2 = 0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations
  thplDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, 
                                                                  fixedFormula = a1 + a2 + a9 + a9p + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                  start = list(fixed = c(a1 = 3.5, a2 = 0, a9 = 0, a9p = 0, b1 = 0.6, b2 = 0.12)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5), significant = FALSE) # singular precision matrix
  thplDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, 
                                                               fixedFormula = a1 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                               start = list(fixed = c(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1)), control = nlmeControl(maxIter = 500), significant = FALSE)
  thplDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                              start = list(fixed = c(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  thplDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016physio, 
                                                                    fixedFormula = a1 + a8 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                    start = list(fixed = c(a1 = 3.0, a8 = -0.01, a9 = 0, b1 = 0.73, b2 = 0.07)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  thplDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", DBH ~ ((a1 + a1r)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016, 
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                 start = list(fixed = c(a1 = -300, b1 = 0.04, b2 = 0.55)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve

  thplDiameterFromHeightMixed$gamm = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9) + s(StandID, bs = "re"), data = thpl2016, mixed = TRUE)
  thplDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16) + s(StandID, bs = "re"), data = thpl2016, mixed = TRUE)
  thplDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9) + s(StandID, bs = "re"), data = thpl2016, mixed = TRUE)

  save(file = "trees/height-diameter/data/THPL DBH mixed.Rdata", thplDiameterFromHeightMixed)
}


## collect model results and parameters
if (thplOptions$fitHeight & thplOptions$fitHeightMixed & thplOptions$fitDbh & thplOptions$fitDbhMixed)
{
  if (exists("thplHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/THPL TotalHt.Rdata") }
  #if (exists("thplHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/THPL TotalHt gnls.Rdata") }
  if (exists("thplHeightFromDiameterMixed") == FALSE) { load("trees/height-diameter/data/THPL TotalHt mixed.Rdata") }
  if (exists("thplDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/THPL DBH.Rdata") }
  if (exists("thplDiameterFromHeightMixed") == FALSE) { load("trees/height-diameter/data/THPL DBH mixed.Rdata") }

  thplCoefficients = bind_rows(bind_rows(bind_rows(lapply(thplHeightFromDiameter, get_list_coefficients)),
                                         #bind_rows(lapply(thplHeightFromDiameterGnls, get_model_coefficients)),
                                         bind_rows(lapply(thplHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(thplHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(thplHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(thplDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(thplDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(thplDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(thplDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "THPL")
  thplResults = bind_rows(bind_rows(bind_rows(lapply(thplHeightFromDiameter, get_list_stats)),
                                    #bind_rows(lapply(thplHeightFromDiameterGnls, get_model_stats)),
                                    bind_rows(lapply(thplHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(thplHeightFromDiameterMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(thplHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(thplDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitSet = "primary", fittingMethod = "gsl_nls"),
                                    bind_rows(lapply(thplDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(thplDiameterFromHeightMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(thplDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "THPL")
  
  check_plot_results(thplResults)
  save(file = "trees/height-diameter/data/THPL results.Rdata", thplCoefficients, thplResults)
}


## preferred forms identified (results.R, Figure 10)
if (thplOptions$fitHeight & thplOptions$fitDbh)
{
  thplHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint, folds = 1, repetitions = 1))
  thplHeightFromDiameterPreferred$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 20, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint, folds = 1, repetitions = 1)
  thplHeightFromDiameterPreferred$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176), folds = 1, repetitions = 1)
  thplHeightFromDiameterPreferred$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176), folds = 1, repetitions = 1)
  thplHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649), folds = 1, repetitions = 1)
  #thplHeightFromDiameterPreferred$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 50.6, a1p = -15.8, b1 = 0.023, b2 = -0.014, b2p = -0.009, b3 = -0.069, b4 = 1.130), folds = 1, repetitions = 1)
  #thplHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10), folds = 1, repetitions = 1)
  thplHeightFromDiameterPreferred$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09), folds = 1, repetitions = 1)
  #thplHeightFromDiameterPreferred$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151), folds = 1, repetitions = 1)
  #AIC(thplHeightFromDiameterPreferred$hossfeld, thplHeightFromDiameterPreferred$michaelisMenten, thplHeightFromDiameterPreferred$prodan, thplHeightFromDiameterPreferred$ratkowsky)
  
  thplDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint, folds = 1, repetitions = 1))
  thplDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), thpl2016, folds = 1, repetitions = 1)
  thplDiameterFromHeightPreferred$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08), folds = 1, repetitions = 1)
  thplDiameterFromHeightPreferred$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = thpl2016, constraint = thpl2016gamConstraint, folds = 1, repetitions = 1)
  thplDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = thpl2016physio, constraint = thpl2016gamConstraint, folds = 1, repetitions = 1)
  
  save(file = "trees/height-diameter/data/THPL preferred models.Rdata", thplHeightFromDiameterPreferred, thplDiameterFromHeightPreferred)
}


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  thplBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), thpl2016, start = list(a1 = 90, b1 = 0.000003, b2 = 2.18), weights = heightWeight^2) # a1p, b1p, b2p not significant
  thplBasalAreaFromHeightPower = gsl_nls(basalArea ~ a1*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), thpl2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, b1 = 2.14, b1p = 0.34), weights = heightWeight^2) # a1p not significant
  #confint2(thplBasalAreaFromHeightPower, level = 0.99)

  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(thplBasalAreaFromHeightKorf), 100^2 * mean(residuals(thplBasalAreaFromHeightKorf)), mean(abs(residuals(thplBasalAreaFromHeightKorf))), 1 - sum(residuals(thplBasalAreaFromHeightKorf)^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2)),
          "power", AIC(thplBasalAreaFromHeightPower), 100^2 * mean(residuals(thplBasalAreaFromHeightPower)), mean(abs(residuals(thplBasalAreaFromHeightPower))), 1 - sum(residuals(thplBasalAreaFromHeightPower)^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(thpl2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(thplBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(thplBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "western redcedar height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## exploratory plots
if (htDiaOptions$includeInvestigatory)
{
  library(GGally)
  ggpairs(thpl2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger, relativeHeight, regeneration), 
          aes(alpha = 0.1, color = regeneration, shape = "16"),
          columnLabels = c("DBH, cm", "height, m", "BA, m² ha⁻¹", "BAL, m² ha⁻¹", "relative height, %", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    #scale_color_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) + # https://github.com/ggobi/ggally/issues/445
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  ggpairs(thpl2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, slope, elevation, topographicShelterIndex, regeneration), 
          aes(alpha = 0.1, color = if_else(thpl2016$isPlantation, "plantation", "natural regen"), shape = "16"), 
          columnLabels = c("DBH, cm", "height, m", "slope, °", "elevation, m", "TSI, °", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  scatterPlotMatrix::scatterPlotMatrix(thpl2016 %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger))
}
