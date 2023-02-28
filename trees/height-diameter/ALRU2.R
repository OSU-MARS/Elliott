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
alru2016 = trees2016 %>% filter(Species == "RA", isLiveUnbroken, TotalHt > 0) %>% # live red alders measured for height
  mutate(dbhWeight = pmin(1/(1.53*DBH^(0.78 - 0.051*isPlantation)), 5),
         heightWeight = pmin(1/(26.9*(TotalHt - 1.37)^(0.64 - 0.23*isPlantation)), 5))
alru2016physio = alru2016 %>% filter(is.na(elevation) == FALSE)
alru2016gamConstraint = c(DBH = -1.5179/0.7812, TotalHt = 1.37, standBasalAreaPerHectare = median(alru2016$standBasalAreaPerHectare), basalAreaLarger = median(alru2016$basalAreaLarger), standBasalAreaApprox = median(alru2016$standBasalAreaApprox), tallerApproxBasalArea = median(alru2016$tallerApproxBasalArea), elevation = median(alru2016physio$elevation), slope = median(alru2016physio$slope), aspect = median(alru2016physio$aspect), topographicShelterIndex = median(alru2016physio$topographicShelterIndex), relativeHeight = median(alru2016$relativeHeight)) # point constraint for mgcv::s()

alru2016defaultWeight = alru2016 %>% mutate(dbhWeight = pmin(1/DBH, 5),
                                            heightWeight = pmin(1/TotalHt, 5))
alru2016defaultWeightPhysio = alru2016defaultWeight %>% filter(is.na(elevation) == FALSE)

alruOptions = list(fitHeight = FALSE, 
                   fitDbh = TRUE)

if (alruOptions$fitHeight)
{
  alruHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), alru2016))
  alruHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH), alru2016) # isPlantation*DBH^2 not significant
  
  alruHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027))
  alruHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 28.5, a1p = 13.1, a2 = 0.212, a2p = 0.218, a3 = -0.173, b1 = -0.0404, b1p = 0.0198, b2 = 1.143, b2p = -0.122)) # a3p not significant
  alruHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 30.1, a1p = 3.7, a2 = 0.05, a2p = 0.083, a3 = 0, a4 = -0.002, a5 = -0.28, a6 = 0.074, a7 = 0.26, a8 = 0.21, b1 = -0.052, b1p = 0.007, b2 = 1.24)) # a4, a6, a7 not significant, job step factor with nlrob()
  alruHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750)) # a2, a3, a3p, b2 not significant, often step factor with nlrob()
  alruHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 33.7, a1p = 4.73, a5 = -0.195, a8 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14)) # a4, a6, a7, b2p not significant
  alruHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016, start = list(a1 = 1.6, b1 = 0.24))
  alruHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 27, a1p = 7, b1 = 60, b2 = -1.4, b2p = 0.12)) # b1p-b2p not mutually significant
  alruHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^b2), alru2016, start = list(a1 = 223, a1p = 22.1, b1 = -6.04, b2 = -0.252)) # b1p, b2p not significant
  alruHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152)) # a2p not significant
  alruHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3), alru2016, start = list(a1 = 0.03, a1p = -0.005, a2 = 0.7, a3 = 2.0)) # a2p, a3p not significant
  alruHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047)) # a1p not significant
  alruHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016, start = list(a1 = 33.6, a1p = 7.31, b1 = -22.3, b1p = -4.77, b2 = 5.31, b2p = 1.37))
  alruHeightFromDiameter$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361))
  alruHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, alru2016, start = list(a1 = 21.3, b1 = 0.038, b1p = 0.068, b2 = -0.039, b3 = 0.084, b4 = 1.19)) # b2p, b3p, b4p not significant
  alruHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016, start = list(a1 = 19.0, a1p = 2.9, b1 = 0.090, b2 = -0.064, b3 = -0.14, b4 = 1.18)) # b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016physio, start = list(a1 = 18.1, a1p = 2.80, a5 = 0, a8 = 0.002, b1 = 0.097, b2 = -0.065, b3 = -0.131, b4 = 1.18)) # a3, a4, a6, a7, b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, alru2016physio, start = list(a1 = 15.2, a1p = 2.54, a5 = -0.13, a8 = 0.10, b1 = 0.17, b2 = -0.068, b3 = -0.07, b4 = 1.26)) # a4, a6, a7, b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), alru2016, start = list(a1 = 22.4, b1 = 0.029, b1p = 0.071, b2 = -0.040, b3 = -0.026, b3p = -0.053, b4 = 1.18, b4p = -0.123)) # a1p, b1, b2p not significant
  alruHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, alru2016, start = list(a1 = 24, a1p = -2, a2 = 0.1, b1 = 0.0, b1p = 0.05, b2 = -0.1, b3 = -0.14, b4 = 1.1)) # a3p, b1, b2p, b3p, b4p not significant, job step factor with nlrob()
  alruHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044))
  alruHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083)) # b1p not significant
  alruHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 30.4, a2 = 0.238, a3 = -0.215, a3p = 0.186, b1 = -0.0247, b2 = 1.11, b2p = -0.052)) # a1p, a2p, b1p not significant
  alruHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 2.4, a2 = -0.018, a3 = 0.128, a3p = 0.128, a4 = 15.2, b1 = -0.095, b2 = 0.970, b2p = -0.219)) # a1p, a2p, a4p, b1p not significant
  
  alruHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027)))
  #alruHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 30.1, a1p = 3.7, a2 = 0.05, a2p = 0.083, a3 = 0, a4 = -0.002, a5 = -0.28, a6 = 0.074, a7 = 0.26, a8 = 0.21, b1 = -0.052, b1p = 0.007, b2 = 1.24)) # a4, a6, a7 not significant, job step factor with nlrob()
  #alruHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750)) # a2, a3, a3p, b2 not significant, often step factor with nlrob()
  alruHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 28.5, a1p = 13.1, a2 = 0.212, a2p = 0.218, a3 = -0.173, b1 = -0.0404, b1p = 0.0198, b2 = 1.143, b2p = -0.122)) # a3p not significant
  alruHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 33.7, a1p = 4.73, a5 = -0.195, a8 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14)) # a4, a6, a7, b2p not significant
  alruHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016, start = list(a1 = 1.6, b1 = 0.24))
  alruHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 27, a1p = 7, b1 = 60, b2 = -1.4, b2p = 0.12)) # b1p-b2p not mutually significant
  alruHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^b2), alru2016, start = list(a1 = 223, a1p = 22.1, b1 = -6.04, b2 = -0.252)) # b1p, b2p not significant
  alruHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152)) # a2p not significant
  alruHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3), alru2016, start = list(a1 = 0.03, a1p = -0.005, a2 = 0.7, a3 = 2.0)) # a2p, a3p not significant
  alruHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047)) # a1p not significant
  alruHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016, start = list(a1 = 33.6, a1p = 7.31, b1 = -22.3, b1p = -4.77, b2 = 5.31, b2p = 1.37))
  alruHeightFromDiameterNlrob$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361))
  alruHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, alru2016, start = list(a1 = 21.3, b1 = 0.038, b1p = 0.068, b2 = -0.039, b3 = 0.084, b4 = 1.19)) # b2p, b3p, b4p not significant
  alruHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016, start = list(a1 = 19.0, a1p = 2.9, b1 = 0.090, b2 = -0.064, b3 = -0.14, b4 = 1.18)) # b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016physio, start = list(a1 = 18.1, a1p = 2.80, a5 = 0, a8 = 0.002, b1 = 0.097, b2 = -0.065, b3 = -0.131, b4 = 1.18)) # a3, a4, a6, a7, b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, alru2016physio, start = list(a1 = 15.2, a1p = 2.54, a5 = -0.13, a8 = 0.10, b1 = 0.17, b2 = -0.068, b3 = -0.07, b4 = 1.26)) # a4, a6, a7, b1p, b2p, b3p, b4p not significant
  alruHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), alru2016, start = list(a1 = 22.4, b1 = 0.029, b1p = 0.071, b2 = -0.040, b3 = -0.026, b3p = -0.053, b4 = 1.18, b4p = -0.123)) # a1p, b1, b2p not significant
  #alruHeightFromDiameterNlrob$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, alru2016, start = list(a1 = 24, a1p = -2, a2 = 0.1, b1 = 0.0, b1p = 0.05, b2 = -0.1, b3 = -0.14, b4 = 1.1)) # a3p, b1, b2p, b3p, b4p not significant, job step factor with nlrob()
  alruHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044))
  alruHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083)) # b1p not significant
  alruHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 30.4, a2 = 0.238, a3 = -0.215, a3p = 0.186, b1 = -0.0247, b2 = 1.11, b2p = -0.052)) # a1p, a2p, b1p not significant
  alruHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 2.4, a2 = -0.018, a3 = 0.128, a3p = 0.128, a4 = 15.2, b1 = -0.095, b2 = 0.970, b2p = -0.219)) # a1p, a2p, a4p, b1p not significant
  #lapply(alruHeightFromDiameter$sharmaPartonPhysio$fit, confint_nlrob)
  
  alruHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016defaultWeight, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027)))
  alruHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016defaultWeight, start = list(a1 = 28.5, a1p = 13.1, a2 = 0.212, a2p = 0.218, a3 = -0.173, b1 = -0.0404, b1p = 0.0198, b2 = 1.143, b2p = -0.122))
  alruHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016defaultWeightPhysio, start = list(a1 = 30.1, a1p = 3.7, a2 = 0.05, a2p = 0.083, a3 = 0, a4 = -0.002, a5 = -0.28, a6 = 0.074, a7 = 0.26, a8 = 0.21, b1 = -0.052, b1p = 0.007, b2 = 1.24))
  alruHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016defaultWeight, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750))
  alruHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016defaultWeightPhysio, start = list(a1 = 33.7, a1p = 4.73, a5 = -0.195, a8 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14))
  alruHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016defaultWeight, start = list(a1 = 1.6, b1 = 0.24))
  alruHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1*DBH^(b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 27, a1p = 7, b1 = 60, b2 = -1.4, b2p = 0.12))
  alruHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^b2), alru2016defaultWeight, start = list(a1 = 223, a1p = 22.1, b1 = -6.04, b2 = -0.252))
  alruHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016defaultWeight, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152))
  alruHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3), alru2016defaultWeight, start = list(a1 = 0.03, a1p = -0.005, a2 = 0.7, a3 = 2.0))
  alruHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016defaultWeight, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047))
  alruHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 33.6, a1p = 7.31, b1 = -22.3, b1p = -4.77, b2 = 5.31, b2p = 1.37))
  alruHeightFromDiameterGslNlsDefault$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016defaultWeight, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361))
  alruHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, alru2016defaultWeight, start = list(a1 = 21.3, b1 = 0.038, b1p = 0.068, b2 = -0.039, b3 = 0.084, b4 = 1.19))
  alruHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016defaultWeight, start = list(a1 = 19.0, a1p = 2.9, b1 = 0.090, b2 = -0.064, b3 = -0.14, b4 = 1.18))
  alruHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016defaultWeightPhysio, start = list(a1 = 18.1, a1p = 2.80, a5 = 0, a8 = 0.002, b1 = 0.097, b2 = -0.065, b3 = -0.131, b4 = 1.18))
  alruHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, alru2016defaultWeightPhysio, start = list(a1 = 15.2, a1p = 2.54, a5 = -0.13, a8 = 0.10, b1 = 0.17, b2 = -0.068, b3 = -0.07, b4 = 1.26))
  alruHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), alru2016defaultWeight, start = list(a1 = 22.4, b1 = 0.029, b1p = 0.071, b2 = -0.040, b3 = -0.026, b3p = -0.053, b4 = 1.18, b4p = -0.123))
  alruHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, alru2016defaultWeight, start = list(a1 = 24, a1p = -2, a2 = 0.1, b1 = 0.0, b1p = 0.05, b2 = -0.1, b3 = -0.14, b4 = 1.1))
  alruHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044))
  alruHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016defaultWeight, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083))
  alruHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016defaultWeight, start = list(a1 = 30.4, a2 = 0.238, a3 = -0.215, a3p = 0.186, b1 = -0.0247, b2 = 1.11, b2p = -0.052))
  alruHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016defaultWeight, start = list(a1 = 2.4, a2 = -0.018, a3 = 0.128, a3p = 0.128, a4 = 15.2, b1 = -0.095, b2 = 0.970, b2p = -0.219))
  
  alruHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 9, pc = alru2016gamConstraint), data = alru2016)
  alruHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 18, pc = alru2016gamConstraint), data = alru2016)
  alruHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4) # all predictors k = 451, edf < 270, AIC 18588: 18709 without BA, 18735 without BAL, 18505 without elevation, 18872 without slope, 18557 without sin(aspect), 18598 without cos(aspect), 18629 without topographic shelter -> eliminate elevation AIC 18505: 18355 without BA, 18423 without BAL, 18686 without slope, 18334 without sin(aspect), 18320 without cos(aspect), 18342 without topographic shelter -> eliminate cos(aspect)
  alruHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4)
  
  save(file = "trees/height-diameter/data/ALRU2 TotalHt.Rdata", alruHeightFromDiameter, alruHeightFromDiameterNlrob, alruHeightFromDiameterGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
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
if (htDiaOptions$fitGnls)
{
  alruHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1)) # corSymm marginally viable but dropped, step halving nlsTol = 0.01
  alruHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250), folds = 1, repetitions = 1) # corSymm marginally viable but dropped, step halving at nlsTol = 0.02, >250 iterations at nlsTol = 0.05
  alruHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, alru2016, start = alruHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # step halving at nlsTol = 0.5, > 250+50 iterations at nlsTol = 1 or with tighter tolerances than nlsTol = 0.2, msTol = 1E-5, tolerance = 1E-4
  alruHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016, start = alruHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # step halving at nlsTol = 0.1, maxiter at default msTol and tolerance
  alruHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), alru2016, start = alruHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1)  # corSymm marginally viable but dropped, step halving at nlsTol = 0.02
  alruHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, alru2016, start = alruHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1)  # corSymm viable but dropped, step halving at nlsTol = 0.01
  alruHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1)  # corSymm marginally viable but dropped, step factor at nlsTol = 0.02
  alruHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # corSymm viable but dropped, step halving at nlsTol = 0.01

  save(alruHeightFromDiameterGnls, file = "trees/height-diameter/data/ALRU2 TotalHt GNLS.rdata")
}
if (htDiaOptions$includeInvestigatory)
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
if (alruOptions$fitDbh)
{
  alruDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), alru2016))
  alruDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), alru2016)
  
  alruDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22)) # a1p not significant
  alruDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0)) # a1p not significant, all NaN-inf with a3 * standBasalAreaApprox
  alruDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -28, a2 = 1.9, a2p = -0.71, a3 = -1.67, a3p = 0.7, b1 = -0.033, b2 = 1.33)) # a1p, b2p not significant
  alruDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a9 = 29.8, a9p = -14.8, b1 = -0.0302, b2 = 1.30)) # job step factor with nlrob()
  alruDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -54.2, a9 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19)) # a1p not significant
  alruDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50)) # b1p not significant
  alruDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), control = list(maxiter = 50)) # a2p not significant
  alruDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a5 = 7.37, a6 = 0.41, a7 = 0.52, b1 = -0.031, b2 = 2.47, b2p = -1.26)) # a4, a8 not significant
  alruDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016, start = list(a1 = 33.0, a1p = -7.77, a9 = -0.95, b1 = -0.043, b2 = 1.39)) # a2p, b1p not significant
  alruDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), alru2016, start = list(a1 = 100, a2 = 100, b1 = 1)) # collapses to linear with or without b1p
  alruDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.8, a1p = -1.0, a2 = -0.12, a2p = -0.006))
  alruDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15))
  #alruDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147))
  #alruDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^b1, alru2016physio, start = list(a1 = 2.38, a1p = -0.70, a5 = 1.26, a6 = 0.078, a7 = 0.0075, b1 = 0.88)) # a4, a8, b1p not significant
  #alruDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.55, a1p = -1.73, a9 = -0.98, a9p = 0.65, b1 = 0.86, b1p = 0.14))
  alruDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.3, b1 = 1.5, b1p = -0.28, b2 = -0.045, b2p = 0.025)) # a1p not significant
  alruDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + (a3 + a3p*isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.2, a2 = -0.013, a2p = 0.016, a3 = 0.01, a3p = -0.12, b1 = 1.5, b1p = -0.27, b2 = -0.047, b2p = 0.031))
  alruDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016physio, start = list(a1 = 1.1, a5 = 0.5, b1 = 1.5, b1p = -0.27, b2 = -0.04, b2p = 0.0027)) # a1p, a4, a6, a7, a8 not significant
  alruDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.2, a9 = -0.1, b1 = 1.5, b1p = -0.26, b2 = -0.043, b2p = 0.025), significant = FALSE) # a9, a9p not significant
  alruDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), alru2016, start = list(a1 = 0.00006, a2 = 0.04, b1 = 0.94, Ha = 100)) # poorly conditioned, singular gradient with Levenberg, nls(), or fit_gsl_nls()
  alruDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, alru2016, start = list(a1 = 29, b1 = 0.6, b2 = 0.0001, b3 = -1.0, b4 = 0.19), control = gsl_nls_control(maxiter = 500)) # nls() NaN-infinity even with parameters from nls_multstart(), nlrob() NaN-inf
  alruDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119))
  alruDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12)) # a2 not significant
  alruDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a5 = 0.231, a6 = 0.191, a7 = 0.018, b1 = 2.13, b2 = -0.163, b2p = 0.0197)) # a4, a8 not significant
  alruDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.712, a1p = -0.341, a9 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348))
  alruDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -60, b1 = 0.1, b2 = 0.57)) # singular gradient with a1p, b1p, b2p
  
  alruDiameterFromHeightNlrob = list(chapmanReplace = fit_nlrob("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22))) # a1p not significant
  alruDiameterFromHeightNlrob$chapmanReplaceAbat = fit_nlrob("Chapman-Richards replace ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0)) # a1p not significant, all NaN-inf with a3 * standBasalAreaApprox
  alruDiameterFromHeightNlrob$chapmanReplaceBal = fit_nlrob("Chapman-Richards replace BA+L", DBH ~ (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -28, a2 = 1.9, a2p = -0.71, a3 = -1.67, a3p = 0.7, b1 = -0.033, b2 = 1.33)) # a1p, b2p not significant
  #alruDiameterFromHeightNlrob$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a9 = 29.8, a9p = -14.8, b1 = -0.0302, b2 = 1.30)) # job step factor with nlrob()
  alruDiameterFromHeightNlrob$chapmanReplaceRelHt = fit_nlrob("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -54.2, a9 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19)) # a1p not significant
  alruDiameterFromHeightNlrob$chapmanRichards = fit_nlrob("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50)) # b1p not significant
  alruDiameterFromHeightNlrob$chapmanRichardsAbat = fit_nlrob("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), control = list(maxiter = 50)) # a2p not significant
  alruDiameterFromHeightNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a5 = 7.37, a6 = 0.41, a7 = 0.52, b1 = -0.031, b2 = 2.47, b2p = -1.26)) # a4, a8 not significant
  alruDiameterFromHeightNlrob$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016, start = list(a1 = 33.0, a1p = -7.77, a9 = -0.95, b1 = -0.043, b2 = 1.39)) # a2p, b1p not significant
  #alruDiameterFromHeightNlrob$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), alru2016, start = list(a1 = 100, a2 = 100, b1 = 1)) # collapses to linear with or without b1p
  alruDiameterFromHeightNlrob$naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.8, a1p = -1.0, a2 = -0.12, a2p = -0.006), maxit = 80)
  alruDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15))
  #alruDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147))
  #alruDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^b1, alru2016physio, start = list(a1 = 2.38, a1p = -0.70, a5 = 1.26, a6 = 0.078, a7 = 0.0075, b1 = 0.88)) # a4, a8, b1p not significant
  #alruDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.55, a1p = -1.73, a9 = -0.98, a9p = 0.65, b1 = 0.86, b1p = 0.14))
  alruDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.35, b1 = 1.37, b1p = -0.24, b2 = -0.033, b2p = 0.022)) # a1p not significant
  alruDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + (a3 + a3p*isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.1, a2 = -0.013, a2p = 0.015, a3 = 0.01, a3p = -0.12, b1 = 1.5, b1p = -0.23, b2 = -0.047, b2p = 0.026)) # a3 not reliably significant
  alruDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016physio, start = list(a1 = 0.85, a5 = 0.55, b1 = 1.45, b1p = -0.27, b2 = -0.038, b2p = 0.0022))
  alruDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.1, a9 = -0.07, b1 = 1.5, b1p = -0.24, b2 = -0.04, b2p = 0.02), significant = FALSE)
  #alruDiameterFromHeightNlrob$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), alru2016, start = list(a1 = 0.00006, a2 = 0.04, b1 = 0.94, Ha = 100)) # poorly conditioned, singular gradient with Levenberg, nls(), or fit_nlrob()
  #alruDiameterFromHeightNlrob$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, alru2016, start = list(a1 = 29, b1 = 0.71, b2 = 0.0001, b3 = -0.65, b4 = 0.21), control = nls.control(maxiter = 500)) # nls() NaN-infinity even with parameters from nls_multstart(), fit_nlrob() NaN-inf
  alruDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119))
  alruDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12))  # a2 not significant
  alruDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a5 = 0.231, a6 = 0.191, a7 = 0.018, b1 = 2.13, b2 = -0.163, b2p = 0.0197)) # a4, a8 not significant
  alruDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.712, a1p = -0.341, a9 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348))
  alruDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -60, b1 = 0.1, b2 = 0.55)) # singular gradient with a1p, b1p, b2p
  #confint_nlrob(alruDiameterFromHeight$sibbesenReplaceRelHt, level = 0.99)
  
  alruDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016defaultWeight, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22)))
  alruDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016defaultWeight, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0))
  alruDiameterFromHeightGslNlsDefault$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016defaultWeight, start = list(a1 = -28, a2 = 1.9, a2p = -0.71, a3 = -1.67, a3p = 0.7, b1 = -0.033, b2 = 1.33))
  alruDiameterFromHeightGslNlsDefault$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016defaultWeight, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a9 = 29.8, a9p = -14.8, b1 = -0.0302, b2 = 1.30))
  alruDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016defaultWeight, start = list(a1 = -54.2, a9 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19))
  alruDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016defaultWeight, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50))
  alruDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016defaultWeight, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), control = list(maxiter = 50))
  alruDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016defaultWeightPhysio, start = list(a1 = 6.3, a1p = 30.9, a5 = 7.37, a6 = 0.41, a7 = 0.52, b1 = -0.031, b2 = 2.47, b2p = -1.26))
  alruDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016defaultWeight, start = list(a1 = 33.0, a1p = -7.77, a9 = -0.95, b1 = -0.043, b2 = 1.39))
  alruDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), alru2016defaultWeight, start = list(a1 = 100, a2 = 100, b1 = 1))
  alruDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016defaultWeight, start = list(a1 = 3.8, a1p = -1.0, a2 = -0.12, a2p = -0.006))
  alruDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016defaultWeight, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15))
  #alruDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016defaultWeight, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147))
  #alruDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^b1, alru2016defaultWeightPhysio, start = list(a1 = 2.38, a1p = -0.70, a5 = 1.26, a6 = 0.078, a7 = 0.0075, b1 = 0.88))
  #alruDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016defaultWeight, start = list(a1 = 3.55, a1p = -1.73, a9 = -0.98, a9p = 0.65, b1 = 0.86, b1p = 0.14))
  alruDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016defaultWeight, start = list(a1 = 1.3, b1 = 1.5, b1p = -0.28, b2 = -0.045, b2p = 0.025))
  alruDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea + (a3 + a3p*isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016defaultWeight, start = list(a1 = 1.2, a2 = -0.013, a2p = 0.016, a3 = 0.01, a3p = -0.12, b1 = 1.5, b1p = -0.27, b2 = -0.047, b2p = 0.031))
  alruDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016defaultWeightPhysio, start = list(a1 = 1.1, a5 = 0.5, b1 = 1.5, b1p = -0.27, b2 = -0.04, b2p = 0.0027))
  alruDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016defaultWeight, start = list(a1 = 1.2, a9 = -0.1, b1 = 1.5, b1p = -0.26, b2 = -0.043, b2p = 0.025), significant = FALSE)
  alruDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), alru2016defaultWeight, start = list(a1 = 0.00006, a2 = 0.04, b1 = 0.94, Ha = 100))
  alruDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, alru2016defaultWeight, start = list(a1 = 39, b1 = 0.63, b2 = 0.0001, b3 = -1.2, b4 = 0.19), control = gsl_nls_control(maxiter = 500))
  alruDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119))
  alruDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12))
  alruDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016defaultWeightPhysio, start = list(a1 = 0.716, a1p = -0.336, a5 = 0.231, a6 = 0.191, a7 = 0.018, b1 = 2.13, b2 = -0.163, b2p = 0.0197))
  alruDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016defaultWeight, start = list(a1 = 0.712, a1p = -0.341, a9 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348))
  alruDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016defaultWeight, start = list(a1 = -68, b1 = 0.1, b2 = 0.55))

  alruDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = alru2016gamConstraint), data = alru2016)
  alruDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 14, pc = alru2016gamConstraint), data = alru2016)
  alruDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 496, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 6) # with all predictors k = 496, edf < 220, AIC 24509: 24540 without AAT, 24578 without ABA, 24571 without elevation, 24555 without slope, 24563 without sin(aspect), 24576 without cos(aspect), 24556 without topographic shelter
  alruDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = alru2016gamConstraint), data = alru2016physio, nthreads = 4)
  
  save(file = "trees/height-diameter/data/ALRU2 DBH.Rdata", alruDiameterFromHeight, alruDiameterFromHeightNlrob, alruDiameterFromHeightGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  print(alruDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

  ggplot() +
    geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_smooth(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanReplace), y = alru2016$TotalHt, color = "Chapman-Richards replace", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichards), y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsAbat), y = alru2016$TotalHt, color = "Chapman-Richards ABA+T", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsPhysio), y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$chapmanRichardsRelHt), y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeightKorf), y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$linear), y = alru2016$TotalHt, color = "linear", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$michaelisMentenReplace), y = alru2016$TotalHt, color = "Michaelis-Menten replace", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$naslund), y = alru2016$TotalHt, color = "Näslund", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$parabolic), y = alru2016$TotalHt, color = "parabolic", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$power), y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$ruark), y = alru2016$TotalHt, color = "Ruark", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$schnute), y = alru2016$TotalHt, color = "Schnute inverse", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$sibbesenReplace), y = alru2016$TotalHt, color = "Sibbesen replace", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(alruDiameterFromHeight$weibull), y = alru2016$TotalHt, color = "Weibull", group = alru2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -30*log(1 - pmin((0.03*(alru2016$TotalHt - 1.37))^0.85, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards")) +
    #geom_line(aes(x = (-16.95 + 2.676 * alru2016$isPlantation)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^(0.304 + 0.0723 * alru2016$isPlantation), 0.9999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99999", group = alru2016$isPlantation)) +
    #geom_line(aes(x = (-15.886)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^0.304, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.999")) +
    #geom_line(aes(x = (-18.329)*log(1 - pmin((0.0305*(alru2016$TotalHt - 1.37))^0.427, 0.99)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99")) +
    #geom_line(aes(x = -30*log(1 - pmin(0.03*(alru2016$TotalHt - 1.37)^1.0, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
    #geom_line(aes(x = -150*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.79, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
    #geom_line(aes(x = -119*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.86, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 2")) +
    #geom_line(aes(x = -108*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.89, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 3")) +
    #geom_line(aes(x = 5*(exp(0.33*(alru2016$TotalHt - 1.37)^0.60) - 1), y = alru2016$TotalHt, color = "Chapman-Richards replace")) +
    #geom_line(aes(x = 39*(exp(0.1*(alru2016$TotalHt - 1.37)^0.66) - 1), y = alru2016$TotalHt, color = "Chapman-Richards replace")) +
    #geom_line(aes(x = 125.8*(exp(0.0182*(alru2016$TotalHt - 1.37)^0.880) - 1), y = alru2016$TotalHt, color = "Chapman-Richards replace")) +
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
    #geom_line(aes(x = 30*(alru2016$TotalHt - 1.37) / (60 - (alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Ratkowsky replace")) +
    #geom_line(aes(x = 2*(alru2016$TotalHt - 1.37)^0.5 / (1 - 0.12*(alru2016$TotalHt - 1.37)^0.5), y = alru2016$TotalHt, color = "Ratkowsky")) +
    #geom_line(aes(x = 100*(alru20 16$TotalHt - 1.37)^0.01*(exp(0.01*(alru2016$TotalHt - 1.37)) - 1), y = alru2016$TotalHt, color = "Ruark replace")) +
    #geom_line(aes(x = -1/0.006*log(1 - (1 - exp(-0.5))*(alru2016$TotalHt^1.5 - 1.3^1.5) / (40^1.5 - 1.3^1.5)), y = alru2016$TotalHt, color = "Schnute inverse")) +
    #geom_line(aes(x = 3*(alru2016$TotalHt - 1.37)^(0.30*(alru2016$TotalHt - 1.37)^0.30), y = alru2016$TotalHt, color = "Sibbesen replace")) +
    #geom_line(aes(x = 5.3*(alru2016$TotalHt - 1.37)^(0.35*(alru2016$TotalHt - 1.37)^0.20), y = alru2016$TotalHt, color = "Sibbesen replace")) +
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
if (alruOptions$fitHeight & alruOptions$fitDbh)
{
  if (exists("alruHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/ALRU2 TotalHt.rdata") }
  #if (exists("alruHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/ALRU2 TotalHt GNLS.rdata") }
  if (exists("alruDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/ALRU2 DBH.rdata") }
  
  alruCoefficients = bind_rows(bind_rows(bind_rows(lapply(alruHeightFromDiameter, get_list_coefficients)),
                                         bind_rows(lapply(alruHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(alruHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                         #bind_rows(lapply(alruHeightFromDiameterGnls, get_model_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(alruDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(alruDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(alruDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "ALRU2")
  alruResults = bind_rows(bind_rows(bind_rows(lapply(alruHeightFromDiameter, get_list_stats)),
                                    bind_rows(lapply(alruHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(alruHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                                    #bind_rows(lapply(alruHeightFromDiameterGnls, get_model_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(alruDiameterFromHeight, get_list_stats)),
                                    bind_rows(lapply(alruDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(alruDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "ALRU2")
    
  save(file = "trees/height-diameter/data/ALRU2 results.Rdata", alruCoefficients, alruResults)
}


## preferred forms identified (results.R, Figure 5)
if (alruOptions$fitHeight & alruOptions$fitDbh)
{
  alruHeightFromDiameterPreferred = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027), folds = 1, repetitions = 1))
  alruHeightFromDiameterPreferred$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 30.1, a1p = 3.7, a2 = 0.05, a2p = 0.083, a3 = 0, a4 = -0.002, a5 = -0.28, a6 = 0.074, a7 = 0.26, a8 = 0.21, b1 = -0.052, b1p = 0.007, b2 = 1.24), folds = 1, repetitions = 1)
  alruHeightFromDiameterPreferred$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152), folds = 1, repetitions = 1)
  alruHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, alru2016physio, start = list(a1 = 18.1, a1p = 2.80, a5 = 0, a8 = 0.002, b1 = 0.097, b2 = -0.065, b3 = -0.131, b4 = 1.18), folds = 1, repetitions = 1)
  alruHeightFromDiameterPreferred$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a5 * slope + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, alru2016physio, start = list(a1 = 15.2, a1p = 2.54, a5 = -0.13, a8 = 0.10, b1 = 0.17, b2 = -0.068, b3 = -0.07, b4 = 1.26), folds = 1, repetitions = 1)
  alruHeightFromDiameterPreferred$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083), folds = 1, repetitions = 1)
  
  alruDiameterFromHeightPreferred = list(chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a5 = 7.37, a6 = 0.41, a7 = 0.52, b1 = -0.031, b2 = 2.47, b2p = -1.26), folds = 1, repetitions = 1))
  alruDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = alru2016gamConstraint), data = alru2016, folds = 1, repetitions = 1)
  alruDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), alru2016, folds = 1, repetitions = 1)
  alruDiameterFromHeightPreferred$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a5 = 0.231, a6 = 0.191, a7 = 0.018, b1 = 2.13, b2 = -0.163, b2p = 0.0197), folds = 1, repetitions = 1)
  
  save(file = "trees/height-diameter/data/ALRU2 preferred models.Rdata", alruHeightFromDiameterPreferred, alruDiameterFromHeightPreferred)
}


## robust regression
if (htDiaOptions$includeInvestigatory)
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
if (htDiaOptions$includeInvestigatory)
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