# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western hemlock height-diameter regression form sweep
tshe2016 = trees2016 %>% filter(Species == "WH", isLiveUnbroken, is.na(TotalHt) == FALSE) %>% # live western hemlocks measured for height
  mutate(dbhWeight = pmin(TreeCount/(1.11*DBH^0.84), 5*TreeCount),
         heightWeight = pmin(TreeCount/(1.00*(TotalHt - 1.37)^1.52), 5*TreeCount))
tshe2016physio = tshe2016 %>% filter(is.na(elevation) == FALSE) # 12 trees without physiographic variables
tshe2016gamConstraint = c(DBH = -1.2994/0.6005, TotalHt = 1.37, standBasalAreaPerHectare = median(tshe2016$standBasalAreaPerHectare), basalAreaLarger = median(tshe2016$basalAreaLarger), standBasalAreaApprox = median(tshe2016$standBasalAreaApprox), tallerApproxBasalArea = median(tshe2016$tallerApproxBasalArea), elevation = median(tshe2016physio$elevation), slope = median(tshe2016physio$slope), aspect = median(tshe2016physio$aspect), topographicShelterIndex = median(tshe2016physio$topographicShelterIndex), relativeHeight = median(tshe2016$relativeHeight), relativeDiameter = median(tshe2016$relativeDiameter)) # point constraint for mgcv::s()

tshe2016defaultWeight = tshe2016 %>% mutate(dbhWeight = pmin(TreeCount/DBH, 5*TreeCount),
                                            heightWeight = pmin(TreeCount/TotalHt, 5*TreeCount))
tshe2016defaultWeightPhysio = tshe2016defaultWeight %>% filter(is.na(elevation) == FALSE)

tsheOptions = tibble(fitHeight = TRUE, 
                     fitHeightGnls = FALSE,
                     fitHeightMixed = FALSE,
                     fitDbh = FALSE,
                     fitDbhMixed = FALSE)

if (tsheOptions$fitHeight)
{
  tsheHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), tshe2016))
  tsheHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), tshe2016)
  
  tsheHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 58, a1p = -17, b1 = -0.02, b1p = -0.018, b2 = 1.2, b2p = 0.17))
  tsheHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 44, a2 = -0.14, a2p = 0.95, a3 = 0.21, a3p = -0.14, b1 = -0.019, b2 = 1.25)) # a1p, a2, b1p, b2p not significant
  tsheHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 51, a2 = -0.15, a2p = 0.9, a3 = 0.21, a3p = -0.14, a4 = -0.014, a5 = -0.12, b1 = -0.021, b2 = 1.27)) # a1p, a5, a6, a7, a8, b1p, b2p not significant
  tsheHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 51, a2 = -0.15, a2p = 0.9, a3 = 0.21, a3p = -0.14, a4 = -0.014, a5 = -0.12, a10 = 0, b1 = -0.021, b2 = 1.27), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 45, a2 = -0.14, a2p = 0.95, a3 = 0.21, a3p = -0.14, a10 = 0, b1 = -0.019, b2 = 1.25), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -3.7, a1p = 10, a2 = 0, a2p = 0.5, a3 = 0.048, a9 = 62, a9p = -28, b1 = -0.031, b2 = 0.10, b2p = 0.93)) # a2, a3p, b1p not significant
  tsheHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 60, a1p = -6.5, a4 = -0.016, a5 = -8.4, b1 = -0.027, b2 = 1.45, b2p = -0.28)) # a6, a7, a8 not significant, b1p+b2p not both significant
  tsheHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 65, a1p = -22, a10 = 0, b1 = -0.016, b1p = -0.017, b2 = 1.2, b2p = 0.16), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a10 * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 62, a1p = -6.3, a4 = -0.016, a5 = -8.8, a10 = -0.7, b1 = -0.027, b2 = 1.42, b2p = -0.24), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)) # b1 not significant
  tsheHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
  tsheHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016, start = list(a1 = 200, b1 = -7.2, b2 = -0.33)) # a1p, b1p, b2p not significant
  tsheHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)) # {a1p, a2p}+b1p not mutually significant
  tsheHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)) # a3p not significant
  tsheHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
  tsheHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)) # b2p not significant
  tsheHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)) # dp not significant
  tsheHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 11, b1 = 0.36, b2 = -0.031, b3 = -0.12, b3p = 0.13, b4 = 1.27)) # a1p, b1p, b2p, b4p not significant
  tsheHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 13, a1p = -2.9, b1 = 0.38, b2 = -0.021, b2p = -0.018, b3 = -0.07, b4 = 1.28)) # b1p, b3p, b4p not significant
  tsheHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 16, a1p = -2.8, b1 = 0.34, a4 = -0.0034, a5 = -0.035, b2 = -0.021, b2p = -0.021, b3 = -0.06, b4 = 1.29)) # a6, a7, a8, b1p, b3p, b4p not significant
  tsheHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 16, a1p = -2.8, a10 = 0, b1 = 0.34, a4 = -0.0034, a5 = -0.035, b2 = -0.021, b2p = -0.021, b3 = -0.06, b4 = 1.29), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 13, a1p = -2.9, a10 = 0, b1 = 0.38, b2 = -0.021, b2p = -0.022, b3 = -0.07, b4 = 1.28), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 16, a1p = -2.8, a4 = -0.003, a5 = -0.035, b1 = 0.36, b2 = -0.022, b2p = -0.022, b3 = -0.06, b4 = 1.28)) # a6, a7, a8, b1p, b3p, b4p not significant
  tsheHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 13, a10 = 1.2, a10p = -1.1, b1 = 0.28, b2 = -0.036, b3 = -0.12, b3p = 0.13, b4 = 1.27))
  tsheHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 16, a1p = -3.1, a4 = -0.003, a5 = -0.035, a10 = 0, b1 = 0.36, b2 = -0.022, b2p = -0.022, b3 = -0.08, b4 = 1.28), significant = FALSE) # a10, a10p not significant
  tsheHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)) # b1p, b2p, b4p not significant
  tsheHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = list(a1 = 35, a2 = -0.07, a2p = 0.6, b1 = 0.11, b1p = -0.045, b2 = -0.020, b3 = 0.017, b4 = 1.24)) # a1p, a3, b1, b2p, b3p, b4p not significant
  tsheHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
  tsheHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 60, a1p = -25, b1 = -0.01, b2 = 1.10, b2p = 0.023)) # b1p not significant
  tsheHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 44, a2 = -0.15, a2p = 0.9, a3 = 0.20, a3p = -0.16, b1 = -0.0083, b2 = 1.20)) # a1p, b1p, b2p not significant
  tsheHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 4, a2 = 0.06, a2p = 0.29, a3 = 0.084, a9 = 50, a9p = -17.5, b1 = -0.04, b2 = 0.99), control = gsl_nls_control(maxiter = 200)) # a1p, a2, a3p, b1p, b2p not significant

  tsheHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 58, a1p = -19, b1 = -0.018, b1p = -0.018, b2 = 1.23, b2p = -0.17)))
  tsheHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 43, a2 = -0.13, a2p = 0.9, a3 = 0.22, a3p = -0.14, b1 = -0.019, b2 = 1.27))
  tsheHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 50, a2 = -0.14, a2p = 0.84, a3 = 0.22, a3p = -0.13, a4 = -0.015, a5 = -0.10, b1 = -0.022, b2 = 1.29))
  tsheHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1, a1p = 5, a2 = 0, a2p = 0.4, a3 = 0.15, a9 = 35, a9p = -20, b1 = -0.035, b2 = 0.9, b2p = 0.3), control = nls.control(maxiter = 100, tol = 0.001)) # job step factor
  tsheHeightFromDiameterNlrob$chapmanRichardsBalPhysioRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 51, a2 = -0.15, a2p = 0.2, a3 = 0.21, a3p = -0.12, a4 = -0.014, a5 = -0.12, a10 = 0, b1 = -0.021, b2 = 1.27), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheHeightFromDiameterNlrob$chapmanRichardsBalRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 45, a2 = -0.14, a2p = 0.95, a3 = 0.21, a3p = -0.14, a10 = 0, b1 = -0.022, b2 = 1.29), significant = FALSE)
  tsheHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58, a1p = -6, a4 = -0.015, a5 = -8.8, b1 = -0.026, b2 = 1.46, b2p = -0.23), control = nls.control(maxiter = 100, tol = 1E-4)) # step factor
  tsheHeightFromDiameterNlrob$chapmanRichardsRelDbh = fit_nlrob("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 65, a1p = -22, a10 = -1.8, b1 = -0.016, b1p = -0.018, b2 = 1.2, b2p = 0.16), significant = FALSE)
  tsheHeightFromDiameterNlrob$chapmanRichardsRelDbhPhysio = fit_nlrob("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a10 * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 64, a1p = -4, a4 = -0.015, a5 = -9, a10 = -2, b1 = -0.025, b2 = 1.41, b2p = -0.20), significant = FALSE)
  tsheHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054))
  tsheHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
  tsheHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016, start = list(a1 = 200, b1 = -7.2, b2 = -0.33))
  tsheHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264))
  tsheHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93))
  tsheHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
  tsheHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1))
  tsheHeightFromDiameterNlrob$richardsW = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013))
  tsheHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 11, b1 = 0.36, b2 = -0.032, b3 = -0.13, b3p = 0.13, b4 = 1.30), control = nls.control(tol = 0.001)) # job step factor
  tsheHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 15, a1p = -2.7, b1 = 0.31, b2 = -0.021, b2p = -0.019, b3 = -0.06, b4 = 1.29), control = nls.control(tol = 1E-4)) # job step factor
  tsheHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 18, a1p = -2.7, b1 = 0.32, a4 = -0.004, a5 = -0.004, b2 = -0.024, b2p = -0.020, b3 = -0.07, b4 = 1.30))
  tsheHeightFromDiameterNlrob$sharmaPartonBalPhysioRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 18, a1p = -3, a10 = 0.1, b1 = 0.33, a4 = -0.0038, a5 = -0.039, b2 = -0.024, b2p = -0.023, b3 = -0.09, b4 = 1.29), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheHeightFromDiameterNlrob$sharmaPartonBalRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 13, a1p = -3, a10 = 0, b1 = 0.32, b2 = -0.023, b2p = -0.022, b3 = -0.06, b4 = 1.29), control = nls.control(tol = 1E-4), significant = FALSE)
  tsheHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 17, a1p = -2.9, a4 = -0.004, a5 = -0.04, b1 = 0.33, b2 = -0.024, b2p = -0.022, b3 = -0.07, b4 = 1.29))
  tsheHeightFromDiameterNlrob$sharmaPartonRelDbh = fit_nlrob("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 14, a10 = 1.1, a10p = -1.1, b1 = 0.26, b2 = -0.037, b3 = -0.21, b3p = 0.18, b4 = 1.30))
  tsheHeightFromDiameterNlrob$sharmaPartonRelDbhPhysio = fit_nlrob("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 16, a1p = -3, a4 = -0.003, a5 = -0.035, a10 = 0, b1 = 0.34, b2 = -0.024, b2p = -0.023, b3 = -0.08, b4 = 1.30), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 23, a1p = -5.5, b1 = 0.23, b2 = -0.031, b3 = -0.07, b3p = 0.09, b4 = 1.27))
  tsheHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = list(a1 = 33, a2 = -0.05, a2p = 0.65, b1 = 0.13, b1p = -0.04, b2 = -0.020, b3 = -0.02, b4 = 1.27), control = nls.control(maxiter = 100, tol = 1E-4)) # b3 not significant, step factor
  tsheHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
  tsheHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 60, a1p = -20, b1 = -0.01, b2 = 1.1, b2p = 0.24))
  tsheHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 41, a2 = -0.13, a2p = 0.90, a3 = -0.13, a3p = -0.14, b1 = -0.008, b2 = 1.20), control = nls.control(tol = 0.001)) # job step factor
  tsheHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 1, a2 = 0.04, a2p = 0.30, a3 = 0.20, a9 = 35, a9p = -10, b1 = -0.03, b2 = 1.0), control = nls.control(maxiter = 200, tol = 1E-4)) # job step factor
  #lapply(tsheHeightFromDiameterNlrob$sharmaPartonBalPhysio$fit, confint_nlrob, level = 0.99)
  
  tsheHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeight, start = list(a1 = 58, a1p = -20, b1 = -0.017, b1p = -0.017, b2 = 1.24, b2p = 0.16)))
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeight, start = list(a1 = 44, a2 = -0.13, a2p = 0.9, a3 = 0.21, a3p = -0.14, b1 = -0.020, b2 = 1.27))
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeightPhysio, start = list(a1 = 50, a2 = -0.17, a2p = 0.88, a3 = 0.2, a3p = -0.13, a4 = -0.015, a5 = -0.12, b1 = -0.022, b2 = 1.29))
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeight, start = list(a1 = -3.1, a1p = 10, a2 = -0.005, a2p = 0.50, a3 = 0.032, a9 = 63.0, a9p = -28, b1 = -0.029, b2 = 0.007, b2p = 0.97))
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeightPhysio, start = list(a1 = 51, a2 = -0.15, a2p = 0.8, a3 = 0.21, a3p = -0.14, a4 = -0.014, a5 = -0.12, a10 = 0, b1 = -0.021, b2 = 1.27), control = nls.control(tol = 0.1), significant = FALSE) # job step factor
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeight, start = list(a1 = 44, a2 = -0.12, a2p = 0.9, a3 = 0.21, a3p = -0.14, a10 = 0, b1 = -0.020, b2 = 1.26), significant = FALSE)
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeightPhysio, start = list(a1 = 59, a1p = -5.8, a4 = -0.015, a5 = -0, b1 = -0.026, b2 = 1.45, b2p = -0.23))
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeight, start = list(a1 = 65, a1p = -21, a10 = -1, b1 = -0.016, b1p = -0.015, b2 = 1.2, b2p = 0.14), significant = FALSE)
  tsheHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a10 * relativeDiameter) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeightPhysio, start = list(a1 = 63, a1p = -5, a4 = -0.016, a5 = -9, a10 = -1.3, b1 = -0.025, b2 = 1.42, b2p = -0.22), significant = FALSE)
  tsheHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054))
  tsheHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016defaultWeight, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
  tsheHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016defaultWeight, start = list(a1 = 200, b1 = -7.2, b2 = -0.33))
  tsheHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016defaultWeight, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264))
  tsheHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016defaultWeight, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93))
  tsheHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
  tsheHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016defaultWeight, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1))
  tsheHeightFromDiameterGslNlsDefault$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016defaultWeight, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013))
  tsheHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 12, b1 = 0.34, b2 = -0.032, b3 = -0.15, b3p = 0.129, b4 = 1.28))
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 15, a1p = -2.9, b1 = 0.34, b2 = -0.021, b2p = -0.019, b3 = -0.066, b4 = 1.28))
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 19, a1p = -4, b1 = 0.31, a4 = -0.004, a5 = -0.04, b2 = -0.023, b2p = -0.021, b3 = -0.06, b4 = 1.29))
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 17, a1p = -3.1, a10 = 0, b1 = 0.31, a4 = -0.0036, a5 = -0.041, b2 = -0.023, b2p = -0.022, b3 = -0.08, b4 = 1.30), significant = FALSE)
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 14, a1p = -3.1, a10 = 0, b1 = 0.34, b2 = -0.021, b2p = -0.022, b3 = -0.07, b4 = 1.30), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 17, a1p = -3.3, a4 = -0.004, a5 = -0.04, b1 = 0.32, b2 = -0.022, b2p = -0.022, b3 = -0.07, b4 = 1.28))
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + (a10 + a10p * isPlantation) * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 14, a10 = 1.2, a10p = -1.2, b1 = 0.27, b2 = -0.038, b3 = -0.21, b3p = 0.18, b4 = 1.29))
  tsheHeightFromDiameterGslNlsDefault$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 16, a1p = -3.1, a4 = -0.003, a5 = -0.035, a10 = 0.05, b1 = 0.36, b2 = -0.022, b2p = -0.022, b3 = -0.07, b4 = 1.30), control = nls.control(tol = 0.01), significant = FALSE) # job step factor
  tsheHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 25, a1p = -6.2, b1 = 0.19, b2 = -0.031, b3 = -0.08, b3p = 0.09, b4 = 1.27))
  tsheHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 36, a2 = -0.06, a2p = 0.68, b1 = 0.11, b1p = -0.043, b2 = -0.019, b3 = 0.02, b4 = 1.28))
  tsheHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016defaultWeight, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014))
  tsheHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), tshe2016defaultWeight, start = list(a1 = 59, a1p = -22, b1 = -0.01, b2 = 1.12, b2p = -0.22))
  tsheHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016defaultWeight, start = list(a1 = 42, a2 = -0.12, a2p = 0.84, a3 = 0.20, a3p = -0.13, b1 = -0.008, b2 = 1.20))
  tsheHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), tshe2016defaultWeight, start = list(a1 = 4.6, a2 = 0.05, a2p = 0.28, a3 = 0.08, a9 = 50, a9p = -17.5, b1 = -0.038, b2 = 0.98), control = gsl_nls_control(maxiter = 200))
  
  tsheHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 19, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 90, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint) # ~2 minutes to fit+evaluate with all predictors at minimum (Zen 3 @ 3.4 GHz, k = 495, edf < 290) -> eliminate cos(aspect) (k = 330, edf < 240) AIC 11023: 10733 without BA, 10887 without BAL, 10729 without elevation, 10714 without slope, 10679 without sin(aspect), 10669 without topographic shelter -> eliminate sin(aspect)
  tsheHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, topographicShelterIndex, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 28, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", TotalHt ~ s(DBH, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 18, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, relativeDiameter, bs = "ts", k = 331, by = as.factor(isPlantation), pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint)

  save(file = "trees/height-diameter/data/TSHE TotalHt.Rdata", tsheHeightFromDiameter, tsheHeightFromDiameterNlrob, tsheHeightFromDiameterGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  print(tsheHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = tshe2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = tshe2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$curtis), color = "Curtis", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$korf), color = "Korf", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$linear), color = "linear", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$parabolic), color = "parabolic", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$power), color = "power", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$prodan), color = "Prodan", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$richardsW), color = "unified Richards", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sibbesen), color = "Sibbesen", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibull), color = "Weibull", group = tshe2016$isPlantation)) +
    annotate("text", x = 0, y = 70, label = "western hemlock, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 70)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}


## western hemlock height-diameter GNLS regressions
if (tsheOptions$fitHeightGnls)
{
  tsheHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = tsheHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.002, maxIter = 250, nlsMaxIter = 50))) # corSymm viable but dropped, step halving at nlsTol = 0.001
  tsheHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = tsheHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50)) # corSymm marginally viable but dropped, step halving at nlsTol = 0.02
  tsheHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50)) # corSymm viable but dropped
  tsheHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # corSymm viable but dropped
  tsheHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50)) # corSymm viable but dropped
  tsheHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250)) # corSymm dropped, step halving at nlsTol = 0.01
  tsheHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = tsheHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-7, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50)) # corSymm dropped, >250+50 iterations at nlsTol = 0.1, ok at msTol = 1E-3 and tol = 1E-3
  tsheHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = tsheHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50)) # corSymm viable but dropped

  save(tsheHeightFromDiameterGnls, file = "trees/height-diameter/data/TSHE GNLS.rdata")
}
if (htDiaOptions$includeInvestigatory)
{
  tsheHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(tsheHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(tsheHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(tsheHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  
  ggplot() +
    geom_point(aes(x = tshe2016natural$DBH, y = tshe2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = tshe2016natural$DBH, y = tshe2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "western hemlock naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = tshe2016plantation$DBH, y = tshe2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = tshe2016plantation$DBH, y = tshe2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "western hemlock plantation height, m")
  
  ggplot() +
    geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibullBAL), color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016natural$DBH, y = predict(tsheHeightFromDiameter$weibullBALnatural), color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016plantation$DBH, y = predict(tsheHeightFromDiameter$weibullBALplantation), color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterBase), color = "base")) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibull), color = "ElliottWeibull")) +
    annotate("text", x = 0, y = 85, label = "a) western hemlock, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


if (tsheOptions$fitHeightMixed)
{
  tsheHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, 
                                                                fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = 58, a1p = -17, b1 = -0.02, b1p = -0.018, b2 = 1.2, b2p = 0.17))))
  tsheHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, 
                                                            fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 44, a2 = -0.14, a2p = 0.95, a3 = 0.21, a3p = -0.14, b1 = -0.019, b2 = 1.25)))
  tsheHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope) * (1 - exp(b1*DBH))^b2, tshe2016physio, 
                                                                  fixedFormula = a1 + a2 + a2p + a3 + a3p + a4 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 51, a2 = -0.15, a2p = 0.9, a3 = 0.21, a3p = -0.14, a4 = -0.014, a5 = -0.12, b1 = -0.021, b2 = 1.27)))
  tsheHeightFromDiameterMixed$chapmanRichardsBalRelHt = fit_nlme("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, 
                                                                 fixedFormula = a1 + a1p + a2 + a2p + a3 + a9 + a9p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                                 start = list(fixed = c(a1 = -3.7, a1p = 10, a2 = 0, a2p = 0.5, a3 = 0.048, a9 = 62, a9p = -28, b1 = -0.031, b2 = 0.10, b2p = 0.93)))
  tsheHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * sin(3.14159/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio,
                                                               fixedFormula = a1 + a1p + a4 + a5 + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 60, a1p = -6.5, a4 = -0.016, a5 = -8.4, b1 = -0.027, b2 = 1.45, b2p = -0.28)))
  tsheHeightFromDiameterMixed$curtis = fit_nlme("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, 
                                                fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                                start = list(fixed = c(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)))
  tsheHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047)))
  tsheHeightFromDiameterMixed$korf = fit_nlme("Korf", TotalHt ~ 1.37 + (a1 + a1r)*exp(b1*DBH^b2), tshe2016, 
                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                              start = list(fixed = c(a1 = 200, b1 = -7.2, b2 = -0.33)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # max iterations
  tsheHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, 
                                                         fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)))
  tsheHeightFromDiameterMixed$prodan = fit_nlme("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3r), tshe2016, 
                                                fixedFormula = a1 + a1p + a2 + a2p + a3 ~ 1, randomFormula = a3r ~ 1,
                                                start = list(fixed = c(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)), control = nlmeControl(maxIter = 250)) # job >100 iterations
  tsheHeightFromDiameterMixed$power = fit_nlme("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*DBH^(b1 + b1p * isPlantation), tshe2016, 
                                               fixedFormula = a1 + a1p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05)))
  tsheHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, 
                                                   fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)))
  tsheHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, 
                                                   fixedFormula = Ha + Hap + d + kU + kUp ~ 1, randomFormula = Har ~ 1,
                                                   start = list(fixed = c(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)))
  #tsheHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, 
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b3p + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                    start = list(fixed = c(a1 = 11, b1 = 0.36, b2 = -0.031, b3 = -0.12, b3p = 0.13, b4 = 1.27)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #tsheHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, 
  #                                                       fixedFormula = a1 + a1p + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                       start = list(fixed = c(a1 = 13, a1p = -2.9, b1 = 0.38, b2 = -0.021, b2p = -0.018, b3 = -0.07, b4 = 1.28)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singular precision, step halving
  #tsheHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio,
  #                                                             fixedFormula = a1 + a1p + a4 + a5 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                             start = list(fixed = c(a1 = 16, a1p = -2.8, a4 = -0.0034, a5 = -0.035, b1 = 0.34, b2 = -0.021, b2p = -0.021, b3 = -0.06, b4 = 1.29)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  #tsheHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a4 * elevation + a5 * slope)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, 
  #                                                          fixedFormula = a1 + a1p + a4 + a5 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                          start = list(fixed = c(a1 = 16, a1p = -2.8, a4 = -0.003, a5 = -0.035, b1 = 0.36, b2 = -0.022, b2p = -0.022, b3 = -0.06, b4 = 1.28)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singular precision
  #tsheHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, 
  #                                                   fixedFormula = a1 + a1p + b1 + b2 + b3 + b3p + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                   start = list(fixed = c(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #tsheHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, 
  #                                                      fixedFormula = a1 + a2 + a2p + b1 + b1p + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                      start = list(fixed = c(a1 = 35, a2 = -0.07, a2p = 0.6, b1 = 0.11, b1p = -0.045, b2 = -0.020, b3 = 0.017, b4 = 1.24)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #tsheHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, 
  #                                                fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
  #                                                start = list(fixed = c(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singular precision
  tsheHeightFromDiameterMixed$weibull = fit_nlme("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), tshe2016, 
                                                 fixedFormula = a1 + a1p + b1 + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 60, a1p = -25, b1 = -0.01, b2 = 1.10, b2p = 0.023)))
  tsheHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, 
                                                    fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 44, a2 = -0.15, a2p = 0.9, a3 = 0.20, a3p = -0.16, b1 = -0.0083, b2 = 1.20)))
  
  tsheHeightFromDiameterMixed$gamm = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8) + s(StandID, bs = "re"), data = tshe2016, mixed = TRUE)
  tsheHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 19) + s(StandID, bs = "re"), data = tshe2016, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/TSHE TotalHt mixed.Rdata", tsheHeightFromDiameterMixed)
}


## western hemlock diameter-height regressions
if (tsheOptions$fitDbh)
{
  tsheDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), tshe2016)) # isPlantation*(TotalHt - 1.37) not significant(p = 0.036)
  tsheDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016)
  
  tsheDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 40, b1 = 0.028, b2 = 0.77)) # no significant plantation effects
  tsheDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 45, a2 = -0.12, b1 = 0.026, b2 = 0.77))
  tsheDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.25, a2 = -0.01, b1 = 2.8, b2 = 0.22), control = gsl_nls_control(maxiter = 250), significant = FALSE) # a1p, a2, a3, b1p not significant, a3 + b2p step size
  tsheDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.3, a2 = -0.002, a9 = -0.08, b1 = 2.6, b2 = 0.22), control = gsl_nls_control(maxiter = 250), significant = FALSE) # a1, a2, a3, a9 not significant
  tsheDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.3, a9 = 0.06, b1 = 2.6, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE) # a9 not significant, potentially >500 iterations with nlrob()
  tsheDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -95, b1 = 0.027, b2 = 0.81), control = gsl_nls_control(maxiter = 250)) # a1p, b1p not significant
  tsheDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -100, a2 = 0.3, b1 = 0.03, b2 = 0.75)) # a1p, b1p not significant
  tsheDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -85, a5 = -20, b1 = 0.03, b2 = 0.8), significant = FALSE) # a1p, a4, a5, a6, a7, a8, b1p not significant
  tsheDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -75, a9 = -20, b1 = 0.033, b2 = 0.74), significant = FALSE) # a1p, a9, b1p, b2p not significant, a1+a9-b1 not mutually significant
  tsheDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83)) # a1p, a2p, b1p not significant
  tsheDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
  tsheDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04)) # no significant plantation effects
  #tsheDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047)) # a1p not significant
  #tsheDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04)) # a4, a6, a7, a8 not significant
  #tsheDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
  tsheDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018)) # a1p, b1p, b2p not significant
  tsheDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 2.9, a2 = -0.004, a2p = -0.01, b1 = 0.75, b2 = 0.01)) # a2, a3, a3p, b1p, b2p not significant
  tsheDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.5, a3 = -0.0034, a5 = 0.5, b1 = 0.74, b2 = 0.014)) # a2, a3p not significant
  tsheDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.5, a3 = -0.002, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012), significant = FALSE) # a2, a9 not significant
  tsheDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 2.9, a2 = 0, a2p = -0.08, a9 = 0.7, b1 = 0.72, b2 = 0.01), significant = FALSE) # a9, a9p not significant
  tsheDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013)) # a1p, a4, a5, a6, a7, a8, b1p, b2p not significant
  tsheDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01)) # a9p, b2p not significant
  tsheDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012))
  #tsheDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32)) # converges from red alder values but fails to reconverge (singular gradient), NaN-inf or singular gradient with fit_gsl_nls()
  tsheDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016, start = list(a1 = 30, b1 = 0.5, b2 = 0.003, b3 = 0.4, b4 = 0.9), control = gsl_nls_control(maxiter = 250, xtol = 0.03)) # a1-b1 evaporation
  tsheDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057)) # no significant plantation effects
  tsheDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.7, a2 = -0.009, b1 = 0.6, b2 = 0.12)) # no significant plantation effects
  tsheDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.8, a3 = -0.0045, a5 = 0.5, b1 = 0.5, b2 = 0.15)) # a2, a3p not significant
  tsheDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.9, a3 = -0.005, a5 = 0.5, a9 = 0, b1 = 0.52, b2 = 0.13), significant = FALSE) # a2, a9 not significant
  tsheDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 3.1, a2 = -0.007, a9 = 0, b1 = 0.52, b2 = 0.13), significant = FALSE) # a9, a9p not significant
  tsheDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.115, a5 = 0.486, b1 = 0.736, b2 = 0.0593)) # a1p, a4, a6, a7, a8, b1p, b2p not significant
  tsheDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056))
  tsheDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.6, a5 = 0.7, a9 = 0.4, b1 = 0.5, b2 = 0.13))
  tsheDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -225, b1 = 0.011, b2 = 0.82), control = gsl_nls_control(maxiter = 250)) # a1p, b1p, b2p not significant
  #lapply(tsheDiameterFromHeight$ruarkAbatRelHt$fit, confint2, level = 0.99)
  #lapply(tsheDiameterFromHeight$sibbesenReplaceAbat$fit, get_model_coefficients)
  
  tsheDiameterFromHeightNlrob = list(chapmanReplace = fit_nlrob("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 32, b1 = 0.034, b2 = 0.73)))
  tsheDiameterFromHeightNlrob$chapmanReplaceAbat = fit_nlrob("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 32, a2 = -0.07, b1 = 0.034, b2 = 0.73))
  #tsheDiameterFromHeightNlrob$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.4, a9 = 0.06, b1 = 2.2, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  tsheDiameterFromHeightNlrob$chapmanRichards = fit_nlrob("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -80, b1 = 0.035, b2 = 0.77), control = nls.control(maxiter = 500)) # a1p, b1p not significant
  tsheDiameterFromHeightNlrob$chapmanRichardsAbat = fit_nlrob("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -85, a2 = 0.4, b1 = 0.033, b2 = 0.78), control = nls.control(maxiter = 500), significant = FALSE) # a1p, a2, b1p not significant
  tsheDiameterFromHeightNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -65, a5 = -18, b1 = 0.034, b2 = 0.77), control = nls.control(maxiter = 500)) # a1p, a4, b1p, b2p not significant
  tsheDiameterFromHeightNlrob$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -60, a9 = -20, b1 = 0.045, b2 = 0.70), control = nls.control(maxiter = 500), significant = FALSE)
  tsheDiameterFromHeightNlrob$michaelisMentenReplace = fit_nlrob("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83))
  tsheDiameterFromHeightNlrob$naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
  tsheDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04))
  #tsheDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047))
  #tsheDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04))
  #tsheDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
  tsheDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018))
  tsheDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.0, a2 = -0.009, a2p = -0.07, b1 = 0.70, b2 = 0.015))
  tsheDiameterFromHeightNlrob$ruarkAbatPhysio = fit_nlrob("Ruark ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.5, a3 = -0.0034, a5 = 0.5, b1 = 0.74, b2 = 0.014))
  tsheDiameterFromHeightNlrob$ruarkAbatPhysioRelHt = fit_nlrob("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.5, a3 = 0, a5 = 0.65, a9 = 0.7, b1 = 0.68, b2 = 0.014), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # job step factor
  tsheDiameterFromHeightNlrob$ruarkAbatRelHt = fit_nlrob("Ruark ABA+T RelHt", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 2.9, a2 = 0, a2p = -0.07, a9 = 0.9, b1 = 0.66, b2 = 0.013), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013))
  tsheDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01))
  tsheDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, a9 = 0.5, b1 = 0.65, b2 = 0.015))
  tsheDiameterFromHeightNlrob$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016, start = list(a1 = 30, b1 = 0.5, b2 = 0.003, b3 = 0.4, b4 = 0.9), control = nls.control(maxiter = 250, tol = 1)) # step factor with tol < 1
  tsheDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057))
  tsheDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 3.1, a2 = -0.001, b1 = 0.49, b2 = 0.15))
  tsheDiameterFromHeightNlrob$sibbesenReplaceAbatPhysio = fit_nlrob("Sibbesen replace ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 3.0, a3 = -0.004, a5 = 0.5, b1 = 0.5, b2 = 0.15)) # a2, a3p not significant
  tsheDiameterFromHeightNlrob$sibbesenReplaceAbatPhysioRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 3.0, a3 = -0.005, a5 = 0.5, a9 = 0, b1 = 0.50, b2 = 0.16), significant = FALSE) # a2, a9 not significant
  tsheDiameterFromHeightNlrob$sibbesenReplaceAbatRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 3.1, a2 = 0, a9 = 0, b1 = 0.52, b2 = 0.14), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
  tsheDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.8, a5 = 0.55, b1 = 0.5, b2 = 0.15))
  tsheDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 3.0, a9 = 0.35, b1 = 0.49, b2 = 0.15))
  tsheDiameterFromHeightNlrob$sibbesenReplaceRelHtPhysio = fit_nlrob("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.6, a5 = 0.7, a9 = 0.45, b1 = 0.5, b2 = 0.14))
  tsheDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -110, b1 = 0.068, b2 = 0.54), control = nls.control(maxiter = 500))
  #lapply(tsheDiameterFromHeightNlrob$chapmanRichardsAbat$fit, confint_nlrob)
  #lapply(tsheDiameterFromHeightNlrob$chapmanRichardsAbat$fit, get_model_coefficients)
  
  tsheDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016defaultWeight, start = list(a1 = 32, b1 = 0.034, b2 = 0.73)))
  tsheDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016defaultWeight, start = list(a1 = 3.5, a2 = -0.4, b1 = 0.34, b2 = 0.73)) # a1-b1 separation
  tsheDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016defaultWeight, start = list(a1 = 0.4, a9 = 0.06, b1 = 1.2, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE) # a1+a9-b1 separation
  tsheDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -170, b1 = 0.01, b2 = 0.93), control = gsl_nls_control(maxiter = 250)) # a1p, b1p not significant
  tsheDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -200, a2 = 0.3, b1 = 0.01, b2 = 0.93)) # a1p, b1p not significant, a1-b1 separation
  tsheDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeightPhysio, start = list(a1 = -200, a5 = -45, b1 = 0.01, b2 = 0.93), significant = FALSE) # a1p, a4, a5, b1p not significant
  tsheDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -1.7, a9 = -5.8, b1 = 0.8, b2 = 0.1), significant = FALSE) # prone to a1-b1 separation
  tsheDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016defaultWeight, start = list(a1 = 153, a2 = 68, b1 = 0.83))
  tsheDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
  tsheDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016defaultWeight, start = list(a1 = 1.52, b1 = 1.04))
  #tsheDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047))
  #tsheDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016defaultWeightPhysio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04))
  #tsheDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016defaultWeight, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
  tsheDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018))
  tsheDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.2, a2 = -0.004, a2p = -0.01, b1 = 0.90, b2 = 0.01))
  tsheDiameterFromHeightGslNlsDefault$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 1.9, a3 = -0.0030, a5 = 0.35, b1 = 0.9, b2 = 0.006)) # a2, a3p not significant
  tsheDiameterFromHeightGslNlsDefault$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 1.9, a3 = -0.003, a5 = 0.35, a9 = 0, b1 = 0.9, b2 = 0.007), significant = FALSE) # a2, a9 not significant
  tsheDiameterFromHeightGslNlsDefault$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.3, a2 = 0, a2p = -0.08, a9 = 0, b1 = 0.9, b2 = 0.004), significant = FALSE)
  tsheDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013))
  tsheDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.0, a9 = 1.1, b1 = 0.90, b1p = -0.033, b2 = 0)) # b2 not significant
  tsheDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 1.8, a5 = 0.49, a9 = 0.4, b1 = 0.87, b2 = 0.01))
  #tsheDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016defaultWeight, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32))
  tsheDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016defaultWeight, start = list(a1 = 45, b1 = 0.25, b2 = 0.01, b3 = -0.03, b4 = 0.67), control = gsl_nls_control(maxiter = 250, xtol = 0.03))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.3, a2 = -0.01, b1 = 0.74, b2 = 0.05))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.4, a3 = -0.004, a5 = 0.4, b1 = 0.65, b2 = 0.08))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.3, a3 = -0.0035, a5 = 0.45, a9 = 0, b1 = 0.7, b2 = 0.07), significant = FALSE)
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.5, a2 = -0.007, a9 = 0.8, b1 = 0.73, b2 = 0.08), significant = FALSE)
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.8, a5 = 0.5, b1 = 0.5, b2 = 0.15))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.3, a9 = 0.4, b1 = 0.7, b2 = 0.06))
  tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.1, a5 = 0.6, a9 = 0.4, b1 = 0.75, b2 = 0.1))
  tsheDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016defaultWeight, start = list(a1 = -120, b1 = 0.06, b2 = 0.55), control = gsl_nls_control(maxiter = 250))
  
  # individual term selection: TotalHt + ABA + AAT by = isPlantation + slope, cos(aspect) + TSI retained by AIC but not significant (p > 0.41)
  tsheDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint) # k = 495, ef < 210 with all predictors AIC 13838: 13828 without AAT, 13826 without ABA, 13825 without elevation, 13793 without slope, 13848 without sin(aspect), 13896 without cos(aspect), 13772 without topographic shelter -> eliminate topographic shelter (k = 330, edf < 185) AIC 13772: 13654 without AAT, 13665 without ABA, 13649 without elevation, 13649 without slope, 13635 without sin(aspect), 13656 without cos(aspect) -> eliminate sin(aspect)
  tsheDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, cos(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 331, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint)
  tsheDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), constraint = tshe2016gamConstraint, data = tshe2016physio)
  tsheDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint)
  tsheDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 60, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint)

  save(file = "trees/height-diameter/data/TSHE DBH.Rdata", tsheDiameterFromHeight, tsheDiameterFromHeightNlrob, tsheDiameterFromHeightGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  print(tsheDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(tshe2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanReplaceAbat), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    ##geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -1/0.003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute inverse"), alpha = 0.5) +
    #geom_line(aes(x = 30*(TotalHt - 1.37)^0.5*(exp(0.003*(tph/topHeight)^0.4*(TotalHt - 1.37)) - 1)^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 81, label = "western hemlock, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}

if (tsheOptions$fitDbhMixed)
{
  tsheDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", DBH ~ (a1 + a1r) * (exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, 
                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 40, b1 = 0.028, b2 = 0.77)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4))) # max iterations
  tsheDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, 
                                                            fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 45, a2 = -0.12, b1 = 0.026, b2 = 0.77)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # max iterations
  #tsheDiameterFromHeightMixed$chapmanReplaceBal = fit_nlme("Chapman-Richards replace BA+L", DBH ~ (a1 + a1r + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, 
  #                                                         fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                         start = list(fixed = c(a1 = 0.25, a2 = -0.01, b1 = 2.8, b2 = 0.22)), control = nlmeControl(maxIter = 250, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singular precision
  #tsheDiameterFromHeightMixed$chapmanReplaceBalRelHt = fit_nlme("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1r + a2 * basalAreaLarger + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, 
  #                                                              fixedFormula = a1 + a2 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                              start = list(fixed = c(a1 = 0.3, a2 = -0.002, a9 = -0.08, b1 = 2.6, b2 = 0.22)), control = nlmeControl(maxIter = 250, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving
  #tsheDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, 
  #                                                           fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                           start = list(fixed = c(a1 = 0.3, a9 = 0.06, b1 = 2.6, b2 = 0.2)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving
  #tsheDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", DBH ~ (a1 + a1r) * log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, 
  #                                                       fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                       start = list(fixed = c(a1 = -95, b1 = 0.027, b2 = 0.81)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations
  tsheDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = -100, a2 = 0.3, b1 = 0.03, b2 = 0.75)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations
  tsheDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio,
                                                               fixedFormula = a1 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = -85, a5 = -20, b1 = 0.03, b2 = 0.8)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  #tsheDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, 
  #                                                            fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
  #                                                            start = list(fixed = c(a1 = -75, a9 = -20, b1 = 0.033, b2 = 0.74)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singularity in backsolve
  tsheDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, 
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1,
                                                                start = list(fixed = c(a1 = 153, a2 = 68, b1 = 0.83)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  tsheDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", DBH ~ (a1 + a1r + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013)))
  tsheDiameterFromHeightMixed$power = fit_nlme("power", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1, tshe2016, 
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 1.52, b1 = 1.04)))
  #tsheDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", DBH ~ (a1 + a1r + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, 
  #                                                 fixedFormula = a1 + a2 + a2p + b1 + b1p ~ 1, randomFormula = a1r ~ 1,
  #                                                 start = list(fixed = c(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047)))
  #tsheDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, 
  #                                                   fixedFormula = a1 + a5 + b1 ~ 1, randomFormula = a1r ~ 1,
  #                                                   start = list(fixed = c(a1 = 1.33, a5 = 0.284, b1 = 1.04)))
  #tsheDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, 
  #                                                  fixedFormula = a1 + a9 + b1 ~ 1, randomFormula = a1r ~ 1,
  #                                                  start = list(fixed = c(a1 = 1.52, a9 = 0.08, b1 = 1.02)))
  tsheDiameterFromHeightMixed$ruark = fit_nlme("Ruark", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, 
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 2.8, b1 = 0.7, b2 = 0.018)))
  tsheDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", DBH ~ (a1 + a1r + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, 
                                                   fixedFormula = a1 + a2 + a2p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                   start = list(fixed = c(a1 = 2.9, a2 = -0.004, a2p = -0.01, b1 = 0.75, b2 = 0.01)))
  tsheDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, 
                                                         fixedFormula = a1 + a3 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 2.5, a3 = -0.0034, a5 = 0.5, b1 = 0.74, b2 = 0.014)))
  tsheDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio,
                                                              fixedFormula = a1 + a3 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 2.5, a3 = -0.002, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012)), control = nlmeControl(maxIter = 250), significant = FALSE)
  tsheDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", DBH ~ (a1 + a1r + (a2 + a2p*isPlantation) * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, 
                                                        fixedFormula = a1 + a2 + a2p + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 2.9, a2 = 0, a2p = -0.08, a9 = 0.7, b1 = 0.72, b2 = 0.01)), significant = FALSE)
  tsheDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio,
                                                     fixedFormula = a1 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013)), control = nlmeControl(tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # false convergence
  tsheDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", DBH ~ (a1 + a1r + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, 
                                                    fixedFormula = a1 + a9 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01)))
  tsheDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio,
                                                          fixedFormula = a1 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                          start = list(fixed = c(a1 = 2.4, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012)))
  #tsheDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1)), tshe2016, 
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1,
  #                                               start = list(fixed = c(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  #tsheDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016, 
  #                                                    fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
  #                                                    start = list(fixed = c(a1 = 30, b1 = 0.5, b2 = 0.003, b3 = 0.4, b4 = 0.9)), control = nlmeControl(maxIter = 250, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  tsheDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 2.32, b1 = 0.750, b2 = 0.057)))
  tsheDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                             start = list(fixed = c(a1 = 2.7, a2 = -0.009, b1 = 0.6, b2 = 0.12)))
  tsheDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio,
                                                                   fixedFormula = a1 + a3 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                   start = list(fixed = c(a1 = 2.8, a3 = -0.0045, a5 = 0.5, b1 = 0.5, b2 = 0.15)))
  tsheDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio,
                                                                        fixedFormula = a1 + a3 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                        start = list(fixed = c(a1 = 2.9, a3 = -0.005, a5 = 0.5, a9 = 0, b1 = 0.52, b2 = 0.13)), significant = FALSE)
  tsheDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, 
                                                                  fixedFormula = a1 + a2 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 3.1, a2 = -0.007, a9 = 0, b1 = 0.52, b2 = 0.13)), significant = FALSE)
  tsheDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio,
                                                               fixedFormula = a1 + a5 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 2.115, a5 = 0.486, b1 = 0.736, b2 = 0.0593)))
  tsheDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                              start = list(fixed = c(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056)))
  tsheDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", DBH ~ (a1 + a1r + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio,
                                                                    fixedFormula = a1 + a5 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                    start = list(fixed = c(a1 = 2.6, a5 = 0.7, a9 = 0.4, b1 = 0.5, b2 = 0.13)))
  tsheDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", DBH ~ ((a1 + a1r) * log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, 
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = -225, b1 = 0.011, b2 = 0.82)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # job max iterations
  
  tsheDiameterFromHeightMixed$gamm = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8) + s(StandID, bs = "re"), data = tshe2016, mixed = TRUE)
  tsheDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16) + s(StandID, bs = "re"), data = tshe2016, mixed = TRUE)
  tsheDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15) + s(StandID, bs = "re"), data = tshe2016, mixed = TRUE)
  
  save(file = "trees/height-diameter/data/TSHE DBH mixed.Rdata", tsheDiameterFromHeightMixed)
}


## collect model parameters
if (tsheOptions$fitHeight & tsheOptions$fitHeightMixed & tsheOptions$fitDbh & tsheOptions$fitDbhMixed)
{
  if (exists("tsheHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/TSHE TotalHt.Rdata") }
  #if (exists("tsheHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/TSHE TotalHt GNLS.Rdata") }
  if (exists("tsheHeightFromDiameterMixed") == FALSE) { load("trees/height-diameter/data/TSHE TotalHt mixed.Rdata") }
  if (exists("tsheDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/TSHE DBH.Rdata") }
  if (exists("tsheDiameterFromHeightMixed") == FALSE) { load("trees/height-diameter/data/TSHE DBH mixed.Rdata") }
  
  tsheCoefficients = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_coefficients)),
                                         #bind_rows(lapply(tsheHeightFromDiameterGnls, get_model_coefficients)),
                                         bind_rows(lapply(tsheHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(tsheHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(tsheHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(tsheDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(tsheDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(tsheDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                         #get_model_coefficients(tsheDiameterFromHeight$schnute),
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "TSHE")
  tsheResults = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_stats)),
                                    bind_rows(lapply(tsheHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(tsheHeightFromDiameterMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(tsheHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                                    #bind_rows(lapply(tsheHeightFromDiameterGnls, get_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitSet = "primary", fittingMethod = "gsl_nls"),
                                    bind_rows(lapply(tsheDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(tsheDiameterFromHeightMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(tsheDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "TSHE")
  
  check_plot_results(tsheResults)
  save(file = "trees/height-diameter/data/TSHE results.Rdata", tsheCoefficients, tsheResults)
} else if (tsheOptions$fitHeight & tsheOptions$fitHeightMixed & tsheOptions$fitDbh & tsheOptions$fitDbhMixed)
{
  if (exists("tsheHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/TSHE TotalHt.Rdata") }
  if (exists("tsheDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/TSHE DBH.Rdata") }

  tsheCoefficients = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "TSHE")
  tsheResults = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "TSHE")
  
  check_plot_results(tsheResults)
  save(file = "trees/height-diameter/data/TSHE results.Rdata", tsheCoefficients, tsheResults)
}



## preferred forms identified (results.R, Figure 9)
if (tsheOptions$fitHeight & tsheOptions$fitDbh)
{
  tsheHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 58, a1p = -17, b1 = -0.02, b1p = -0.018, b2 = 1.2, b2p = 0.17), folds = 1, repetitions = 1))
  #tsheHeightFromDiameterPreferred$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 44, a2 = -0.14, a2p = 0.95, a3 = 0.21, a3p = -0.14, b1 = -0.019, b2 = 1.25), folds = 1, repetitions = 1))
  tsheHeightFromDiameterPreferred$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 51, a2 = -0.15, a2p = 0.9, a3 = 0.21, a3p = -0.14, a4 = -0.014, a5 = -0.12, b1 = -0.021, b2 = 1.27), folds = 1, repetitions = 1)
  tsheHeightFromDiameterPreferred$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  #tsheHeightFromDiameterPreferred$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 28, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  tsheHeightFromDiameterPreferred$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047), folds = 1, repetitions = 1)
  tsheHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93), folds = 1, repetitions = 1)
  tsheHeightFromDiameterPreferred$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014), folds = 1, repetitions = 1)
  #AIC(tsheHeightFromDiameterPreferred$hossfeld, tsheHeightFromDiameterPreferred$prodan)
  
  tsheDiameterFromHeightPreferred = list(chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a9 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912), folds = 1, repetitions = 1))
  tsheDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  #tsheDiameterFromHeightPreferred$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  #tsheDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = tshe2016physio, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  tsheDiameterFromHeightPreferred$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15, pc = gamConstraint), data = tshe2016, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  tsheDiameterFromHeightPreferred$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), tshe2016, folds = 1, repetitions = 1)
  tsheDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016, folds = 1, repetitions = 1)
  tsheDiameterFromHeightPreferred$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057), folds = 1, repetitions = 1)
  
  save(file = "trees/height-diameter/data/TSHE preferred models.Rdata", tsheHeightFromDiameterPreferred, tsheDiameterFromHeightPreferred)
}


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  tsheBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 5.5, b1 = 0.00003, b2 = 2.09, b2p = -0.060), weights = heightWeight^2) # a1p, b1p not significant, nlrob() step factor on b1p without a1p
  tsheBasalAreaFromHeightPower = gsl_nls(basalArea ~ a1*(imputedHeight - 1.37)^b1, tshe2016, start = list(a1 = 11/7 * 0.25 * pi * 0.01^2, b1 = 2.18), weights = heightWeight^2) # a1p, b1p not significant
  #confint2(tsheBasalAreaFromHeightPower, level = 0.99)
  #coefficients(tsheBasalAreaFromHeightPower)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(tsheBasalAreaFromHeightKorf), 100^2 * mean(residuals(tsheBasalAreaFromHeightKorf)), mean(abs(residuals(tsheBasalAreaFromHeightKorf))), 1 - sum(residuals(tsheBasalAreaFromHeightKorf)^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2)),
          "power", AIC(tsheBasalAreaFromHeightPower), 100^2 * mean(residuals(tsheBasalAreaFromHeightPower)), mean(abs(residuals(tsheBasalAreaFromHeightPower))), 1 - sum(residuals(tsheBasalAreaFromHeightPower)^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(tshe2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(tsheBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(tsheBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "western hemlock height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  tsheHeightGam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint) + 
                            s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 6, pc = gamConstraint) + 
                            s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint) + 
                            s(elevation, bs = "ts", k = 4, pc = gamConstraint) + 
                            s(slope, bs = "ts", k = 4, pc = gamConstraint) + 
                            #s(aspect, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                            #s(topographicShelterIndex, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                            s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 5, pc = gamConstraint), 
                          data = tshe2016physio, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  k.check(tsheHeightGam)
  summary(tsheHeightGam)
  par(mfrow = c(3, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(tsheHeightGam, scale = 0)
  
  tsheDbhGam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = gamConstraint) +
                         s(standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 6, pc = gamConstraint) +
                         s(tallerApproxBasalArea, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint) +
                         s(elevation, bs = "ts", k = 4, pc = gamConstraint) +
                         #s(slope, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                         #s(aspect, bs = "ts", k = 4, pc = gamConstraint) + # not significant
                         s(topographicShelterIndex, bs = "ts", k = 4, pc = gamConstraint) +
                         s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint), 
                       data = tshe2016physio, constraint = tshe2016gamConstraint, folds = 1, repetitions = 1)
  k.check(tsheDbhGam)
  summary(tsheDbhGam)
  plot.gam(tsheDbhGam, scale = 0)
}


## random forest regression
if (htDiaOptions$includeInvestigatory)
{
  startTime = Sys.time()
  tsheHeightForest = train(TotalHt ~ DBH + standBasalAreaPerHectare + basalAreaLarger + elevation + slope + aspect + topographicShelterIndex + relativeDiameter, data = tshe2016physio, method = "ranger", trControl = repeatedCrossValidation, 
                           importance = "impurity_corrected",
                           tuneGrid = expand.grid(mtry = c(4, 5, 6),
                                                  splitrule = "variance",
                                                  min.node.size = c(1, 2, 3)))
  Sys.time() - startTime
  tsheHeightForest
  varImp(tsheHeightForest)
  
  startTime = Sys.time()
  tsheDbhForest = train(DBH ~ TotalHt + standBasalAreaApprox + tallerApproxBasalArea + elevation + slope + aspect + topographicShelterIndex + relativeHeight, data = tshe2016physio, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(7, 8),
                                               splitrule = "variance",
                                               min.node.size = c(4, 5, 6)))
  Sys.time() - startTime
  tsheDbhForest
  varImp(tsheDbhForest)
}
