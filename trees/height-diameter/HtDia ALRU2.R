# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## red alder height-diameter regression form sweep
# R options for nonlinear least squares
#   gslnls::gsl_nls() - fast but fixed weights only, defaults to Levenberg-Marquadt
#   minpack.lm::nlsLM() - weights = wfcs() fails with index out of range, Levenberg-Marquadt only
#   nlme::gnls() - varPower() based weighting fragile andconvergence usually fails
#   robustbase::nlrob() - defaults to iterative reweighted least squares with stats::nls()
#   stats::nls() - fixed weights only, defaults to Gauss-Newton, NL2SOL with algorithm = "port"
# Difference between stats::nls() with fixed weighting and robustbase::nlrob() is relatively small but nlrob() 
# does noticeably decrease asymptoticity of most regression forms (slight decreases occur in some cases).
#alruHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH / (a2 + a2p * isPlantation + DBH), alru2016, start = list(a1 = 45.5, a1p = 18.7, a2 = 49.1, a2p = 20.5), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, d = 1.17, kU = 0.0366), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, kU = 0.0366, kUp = 0), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # confint2() step factor
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - d))))^(1/(1 - (d + dp*isPlantation))), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, dp = 0, kU = 0.0366, kUp = 0), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # confint2() NaN-inf
#alruHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 21.0, a2 = 0.045, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 16.4, a1p = 6.7, a2 = 0.104, a2p = -0.017, b1 = -0.042, b1p = 0.003, b2 = 0.062, b2p = -0.120, b3 = 1.18, b3p = -0.133), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p, a2p, b1p, b2, b2p not significant
#alruHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p, a2, a2p, a3, a4, a6, b1p, b2p, b3p not significant
#alruHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.0019, a1p = 0.163, b1 = 4.98, b1p = -2.72, b2 = -0.175, b2p = 0.0427), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
#alruHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alru2016 = trees2016 %>% filter(Species == "RA", isLiveUnbroken, TotalHt > 0) # live red alders measured for height
alru2016natural = alru2016 %>% filter(isPlantation == FALSE)
alru2016physio = alru2016 %>% filter(is.na(elevation) == FALSE)
alru2016plantation = alru2016 %>% filter(isPlantation)
alru2016plantationPhysio = alru2016physio %>% filter(isPlantation)

alruHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 28.5, a1p = 13.1, a2 = 0.212, a2p = 0.218, a3 = -0.173, b1 = -0.0404, b1p = 0.0198, b2 = 1.143, b2p = -0.122), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a3p not significant
alruHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 35.1, a1p = 6.02, a2 = 0.040, a2p = 0.131, a3 = -0.0050, a4 = -0.241, a5 = 0.364, a6 = 0.102, a7 = 0, b1 = -0.040, b1p = -0.040, b2 = 1.09), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a5, a6 not significant
alruHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2, a3, a3p, b2 not significant
alruHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 33.7, a1p = 4.73, a2 = -0.00025, a3 = -0.195, a4 = 0.195, a5 = -0.113, a6 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a4, a5, b2p not significant
alruHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016, start = list(a1 = 1.6, b1 = 0.24), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterGam = gam(TotalHt ~ s(DBH, by = as.factor(isPlantation)), data = alru2016)
alruHeightFromDiameterGamBal = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation)), data = alru2016)
#alruHeightFromDiameterGamBalPhysio = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = alru2016physio)
#alruHeightFromDiameterGamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = alru2016physio)
alruHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 30.7, a1p = 13.3, b1 = 49.1, b1p = 7.95, b2 = -1.25, b2p = 0.12), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^b2), alru2016, start = list(a1 = 223, a1p = 22.1, b1 = -6.04, b2 = -0.252), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # b1p, b2p not significant
alruHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), alru2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2p not significant
alruHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH), alru2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # isPlantation*DBH^2 not significant
alruHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3 + a3p* isPlantation), alru2016, start = list(a1 = 0.0249, a1p = -0.0054, a2 = 0.938, a3 = 0.761, a3p = -0.177), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2p not significant
alruHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p not significant
alruHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016, start = list(a1 = 33.6, a1p = 7.31, b1 = -22.3, b1p = -4.77, b2 = 5.31, b2p = 1.37), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 21.3, a2 = 0.038, a2p = 0.068, b1 = -0.039, b2 = 0.084, b2p = -0.127, b3 = 1.19, b3p = -0.14), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # b1p not significant
alruHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 17.5, a1p = 5.3, a2 = 0.089, b1 = -0.039, b2 = 0.086, b2p = -0.144, b3 = 1.18, b3p = -0.123), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2p, b1p not significant
alruHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, alru2016physio, start = list(a1 = 18.0, a1p = 3.17, a2 = 0.095, a3 = -0.00008, a4 = -0.012, a5 = -2.41, a6 = 0.0017, b1 = -0.042, b2 = -0.048, b3 = 1.09), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2p, a3, a4, a6, b1p, b2p, b3p not significant
alruHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 16.1, a1p = 5.11, a2 = 0.10, a3 = -0.00008, a4 = -0.012, a5 = -0.0023, a6 = 0.0018, b1 = -0.039, b2 = 0.0081, b2p = -0.123, b3 = 1.19, b3p = -0.14), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a2p, a3, a4, a6, b1p not significant
alruHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 22.4, a2 = 0.029, a2p = 0.071, b1 = -0.040, b2 = -0.026, b2p = -0.053, b3 = 1.18, b3p = -0.123), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p, a2, b1p not significant
alruHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, alru2016, start = list(a1 = 48.9, a1p = -21.2, a2 = -0.178, a2p = 0.197, a3 = 0.006, b1 = -0.042, b2 = -0.039, b3 = 1.07), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a3p, b1p, b2p, b3p not significant
alruHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1))
alruHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 24.3, a1p = -5.66, b1 = -0.0269, b2 = 1.16, b2p = -0.083), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # b1p not significant
alruHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 30.4, a2 = 0.238, a3 = -0.215, a3p = 0.186, b1 = -0.0247, b2 = 1.11, b2p = -0.052), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p, a2p, b1p not significant
alruHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 2.4, a2 = -0.018, a3 = 0.128, a3p = 0.128, a4 = 15.2, b1 = -0.095, b2 = 0.970, b2p = -0.219), weights = pmin(DBH^if_else(isPlantation, -0.6, -0.7), 1)) # a1p, a2p, a4p, b1p not significant
#confint2(alruHeightFromDiameterSibbesen, level = 0.99)

alruHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", alruHeightFromDiameterChapmanRichards, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BA+L", alruHeightFromDiameterChapmanRichardsBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BA+L physio", alruHeightFromDiameterChapmanRichardsBalPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BA+L RelHt", alruHeightFromDiameterChapmanRichardsBalRelHt, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", alruHeightFromDiameterChapmanRichardsPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterCurtis = get_height_error("Curtis", alruHeightFromDiameterCurtis, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterGam = get_height_error("GCV GAM", alruHeightFromDiameterGam, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterGamBal = get_height_error("GCV GAM BA+L", alruHeightFromDiameterGamBal, alru2016, alru2016natural, alru2016plantation)
#alruHeightFromDiameterGamBalPhysio = get_height_error("GCV GAM BA+L physio", alruHeightFromDiameterGamBalPhysio, alru2016physio, alru2016natural, alru2016plantation) # slow
#alruHeightFromDiameterGamPhysio = get_height_error("GCV GAM physio", alruHeightFromDiameterGamPhysio, alru2016physio, alru2016natural, alru2016plantation)
#save(alruHeightFromDiameterGamBalPhysio, alruHeightFromDiameterGamPhysio, file = "trees/height-diameter/HtDia ALRU2 spline height.rdata")
load("trees/height-diameter/HtDia ALRU2 spline height.rdata")
alruHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", alruHeightFromDiameterHossfeld, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterKorf = get_height_error("Korf", alruHeightFromDiameterKorf, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterLinear = get_height_error("linear", alruHeightFromDiameterLinear, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", alruHeightFromDiameterMichaelisMenten, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterParabolic = get_height_error("parabolic", alruHeightFromDiameterParabolic, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterProdan = get_height_error("Prodan", alruHeightFromDiameterProdan, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterPower = get_height_error("power", alruHeightFromDiameterPower, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", alruHeightFromDiameterRatkowsky, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterRichards = get_height_error("unified Richards", alruHeightFromDiameterRichards, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", alruHeightFromDiameterSharmaParton, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BA+L", alruHeightFromDiameterSharmaPartonBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BA+L physio", alruHeightFromDiameterSharmaPartonBalPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", alruHeightFromDiameterSharmaPartonPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", alruHeightFromDiameterSharmaZhang, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BA+L", alruHeightFromDiameterSharmaZhangBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSibbesen = get_height_error("Sibbesen", alruHeightFromDiameterSibbesen, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibull = get_height_error("Weibull", alruHeightFromDiameterWeibull, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBal = get_height_error("Weibull BA+L", alruHeightFromDiameterWeibullBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BA+L RelHt", alruHeightFromDiameterWeibullBalRelHt, alru2016, alru2016natural, alru2016plantation)

alruHeightFromDiameterResults = bind_rows(as_row(alruHeightFromDiameterChapmanRichards),
                                          as_row(alruHeightFromDiameterChapmanRichardsBal),
                                          as_row(alruHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row(alruHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(alruHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(alruHeightFromDiameterCurtis),
                                          as_row(alruHeightFromDiameterGam),
                                          as_row(alruHeightFromDiameterGamBal),
                                          as_row(alruHeightFromDiameterGamBalPhysio),
                                          as_row(alruHeightFromDiameterGamPhysio),
                                          as_row(alruHeightFromDiameterHossfeld),
                                          as_row(alruHeightFromDiameterKorf),
                                          as_row(alruHeightFromDiameterLinear),
                                          as_row(alruHeightFromDiameterMichaelisMenten),
                                          as_row(alruHeightFromDiameterParabolic),
                                          as_row(alruHeightFromDiameterPower),
                                          as_row(alruHeightFromDiameterProdan),
                                          as_row(alruHeightFromDiameterRatkowsky),
                                          as_row(alruHeightFromDiameterRichards),
                                          as_row(alruHeightFromDiameterSharmaParton),
                                          as_row(alruHeightFromDiameterSharmaPartonBal),
                                          as_row(alruHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(alruHeightFromDiameterSharmaPartonPhysio),
                                          as_row(alruHeightFromDiameterSharmaZhang),
                                          as_row(alruHeightFromDiameterSharmaZhangBal),
                                          as_row(alruHeightFromDiameterSibbesen),
                                          as_row(alruHeightFromDiameterWeibull),
                                          as_row(alruHeightFromDiameterWeibullBal),
                                          as_row(alruHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(alruHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterKorf$fitted.values, color = "Korf", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterLinear$fitted.values, color = "linear", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterPower$fitted.values, color = "power", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterProdan$fitted.values, color = "Prodan", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = alru2016$isPlantation)) +
  annotate("text", x = 0, y = 50, label = "red alder, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterChapmanRichards$residuals)), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterChapmanRichards$residuals), color = "Chapman-Richards GAM", fill = "Chapman-Richards GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterChapmanRichards$residuals), color = "Chapman-Richards sqrt(DBH)", fill = "Chapman-Richards sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterChapmanRichards$residuals), color = "Chapman-Richards DBH", fill = "Chapman-Richards DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
  labs(x = NULL, y = "|height error residual|, m", color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterMichaelisMenten$residuals)), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterMichaelisMenten$residuals), color = "Michaelis-Menten GAM", fill = "Michaelis-Menten GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterMichaelisMenten$residuals), color = "Michaelis-Menten sqrt(DBH)", fill = "Michaelis-Menten sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterMichaelisMenten$residuals), color = "Michaelis-Menten DBH", fill = "Michaelis-Menten DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
  labs(x = NULL, y = NULL, color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaParton$residuals)), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaParton$residuals), color = "Sharma-Parton GAM", fill = "Sharma-Parton GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaParton$residuals), color = "Sharma-Parton sqrt(DBH)", fill = "Sharma-Parton sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaParton$residuals), color = "Sharma-Parton DBH", fill = "Sharma-Parton DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
  labs(x = "DBH, cm", y = "|height error residual|, m", color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaZhang$residuals)), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaZhang$residuals), color = "Sharma-Zhang GAM", fill = "Sharma-Zhang GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaZhang$residuals), color = "Sharma-Zhang sqrt(DBH)", fill = "Sharma-Zhang sqrt(DBH)"), alpha = 0.1, formula = y ~ I(sqrt(x)), method = "lm") +
  geom_smooth(aes(x = alru2016$DBH, y = abs(alruHeightFromDiameterSharmaZhang$residuals), color = "Sharma-Zhang DBH", fill = "Sharma-Zhang DBH"), alpha = 0.1, formula = y ~ x, method = "lm") +
  labs(x = "DBH, cm", y = NULL, color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 2, ncol = 2)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alruHeightFromDiameterRatkowsky$w * alruHeightFromDiameterRatkowsky$working.residuals), alpha = 0.1, color = "grey25", shape = 16) +
  labs(x = "DBH, cm", y = "nlrob weight")


## red alder height-diameter GNLS regressions
#alruHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving with default gnlsControl(tolerance = 1E-6, nlsTol = 0.001, msTol = 1E-7), ok at 1E-6, 0.01, 1E-7
#alruHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = alruHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02, >250 iterations at nlsTol = 0.05
#alruHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameterSharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.5, > 250+50 iterations at nlsTol = 1 or with tighter tolerances than nlsTol = 0.2, msTol = 1E-5, tolerance = 1E-4
#alruHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameterSharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.1, maxiter at default msTol and tolerance
#alruHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = alruHeightFromDiameterSharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#alruHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, alru2016, start = alruHeightFromDiameterSharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
#alruHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msVerbose = FALSE, returnObject = FALSE)) # step factor at nlsTol = 0.001
#alruHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), alru2016, start = alruHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.30, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
#save(alruHeightFromDiameterChapmanRichardsGnls, alruHeightFromDiameterChapmanRichardsBalGnls, alruHeightFromDiameterSharmaPartonGnls, alruHeightFromDiameterSharmaPartonBalGnls, alruHeightFromDiameterSharmaZhangGnls, alruHeightFromDiameterSharmaZhangBalGnls, alruHeightFromDiameterWeibullGnls, alruHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia ALRU2 GNLS.rdata")

load("trees/height-diameter/HtDia ALRU2 GNLS.rdata")
alruHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", alruHeightFromDiameterChapmanRichardsGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BA+L GNLS", alruHeightFromDiameterChapmanRichardsBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", alruHeightFromDiameterSharmaPartonGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BA+L GNLS", alruHeightFromDiameterSharmaPartonBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", alruHeightFromDiameterSharmaZhangGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BA+L GNLS", alruHeightFromDiameterSharmaZhangBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", alruHeightFromDiameterWeibullGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BA+L GNLS", alruHeightFromDiameterWeibullBalGnls, alru2016, alru2016natural, alru2016plantation)

alruHeightFromDiameterResultsGnls = bind_rows(as_row(alruHeightFromDiameterChapmanRichardsGnls),
                                              as_row(alruHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(alruHeightFromDiameterSharmaPartonGnls),
                                              as_row(alruHeightFromDiameterSharmaPartonBalGnls),
                                              as_row(alruHeightFromDiameterSharmaZhangGnls),
                                              as_row(alruHeightFromDiameterSharmaZhangBalGnls),
                                              as_row(alruHeightFromDiameterWeibullGnls),
                                              as_row(alruHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
alruHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichardsGnls$fitted.values, color = "Chapman-Richards GNLS", group = alru2016$isPlantation)) +
  annotate("text", x = 0, y = 50, label = "a) red alder, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 115), ylim = c(0, 50)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

## red alder diameter-height regressions
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, start = list(a1 = -18.3, b1 = 0.0305, b2 = 0.427), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 50))
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.9999)), alru2016, start = list(a1 = -15.8, b1 = 0.0304, b2 = 0.362), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 50))
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, algorithm = "port",
#                                            lower = list(a1 = -150, b1 = 0.02, b2 = 0.1), 
#                                            start = list(a1 = -50, b1 = 0.03, b2 = 0.7), 
#                                            upper = list(a1 = -1, b1 = 0.05, b2 = 2), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(warnOnly = TRUE))
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ a1 * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100,
#                                                      start_lower = list(a1 = -20, b1 = -0.03, b2 = -1), 
#                                                      start_upper = list(a1 = -10, b1 = 0.03, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b1p = -0.02, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b1p = 0.00, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b1p = 0.01, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b1p = 0.02, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b1p not significant
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-1, 0.7))
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2, 0.7))
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2.6, 0.7))
#alruDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.1415926/180 * slope) + a4 * cos(3.1415926/180 * aspect) + a5 * sin(3.1415926/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999)), alru2016, iter = 100,
#                                                            start_lower = list(a1 = 10, a1p = 5, a2 = -0.01, a3 = -2, a4 = -2, a5 = -0.5, a6 = -0.1, b1 = -0.1, b1p = -0.1, b2 = -2), # 3.14159, if used, is seen as a mutable parameter
#                                                            start_upper = list(a1 = 20, a1p = 15, a2 = 0.01, a3 = 1, a4 = 2, a5 = 0.5, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightChapmanRichardsRelHt = nls_multstart(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100, 
#                                                           lower = c(a1 = -100, a2 = -10, b1 = 1/60, b2 = 0),
#                                                           start_lower = list(a1 = -35, a2 = -1, b1 = 0.025, b2 = 0.4),
#                                                           start_upper = list(a1 = -25, a2 = 1, b1 = 0.035, b2 = 0.8), 
#                                                           upper = c(a1 = 0.1, a2 = 100, b1 = 1/25, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightKorf = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), lower = list(a1 = 0.5, b1 = 1, b2 = 0.3), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), algorithm = "port")
#alruDiameterFromHeightKorf = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = pmin(TotalHt^-1.3, 0.7), control = list(maxiter = 200, tol = 0.001, warnOnly = TRUE))
#alruDiameterFromHeightKorf = gnls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = varPower(0.65, ~TotalHt), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruDiameterFromHeightLinear = lm(DBH ~ exp(1*(TotalHt - 1.37)^0.4) + 0, alru2016, weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightLinear = lm(DBH ~ I(TotalHt - 1.37) + 0, alru2016, weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightMichaelisMenten = gsl_nls(DBH ~ (a1*(TotalHt - 1.37)/(a2 - (TotalHt - 1.37)))^b1, alru2016, start = list(a1 = 25, a2 = 50, b1 = 1), weights = TotalHt^-2, control = gsl_nls_control(scale = "levenberg")) # inversion of Michaelis-Menten collapses to linear
#alruDiameterFromHeightNaslund = nls(DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.5, a2 = -0.12), weights = TotalHt^-2)
#alruDiameterFromHeightNaslund = gsl_nls(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 4.6, a1p = -1.5, a2 = -0.11, a2p = -0.014), weights = TotalHt^-2)
#alruDiameterFromHeightPower = nls(DBH ~ a1*(TotalHt - 1.37)^b1, alru2016, algorithm = "port",
#                                  lower = list(a1 = 0, b1 = 1.2), # b1 < 1 is asymptotic in diameter rather than height
#                                  start = list(a1 = 0.4, b1 = 1.5), 
#                                  upper = list(a1 = 10, b2 = 3), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightRatkowsky = gsl_nls(DBH ~ a1 / log(a2 * (TotalHt - 1.37)) - a3, alru2016, start = list(a1 = 1, a2 = 1, a3 = 0)) # not numerically stable due to log underruns
#alruDiameterFromHeightRuark = nls(DBH ~ a1*(TotalHt - 1.37)^b1*exp(b2*(TotalHt - 1.37)), alru2016, start = list(a1 = 1.29, b1 = 1.31, b2 = -0.027), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightSchnute = nls(DBH ~ a1*log(1 - b1*(TotalHt^b2 - 1.37^b2) / (40^b2 - 1.37^b2)), alru2016, start = list(a1 = -150, b1 = 0.5, b2 = 1.5), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # diverges to concave up
#alruDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, alru2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = 0.001, b2 = 0.1, b3 = 0.1), 
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 0.1, b2 = 1, b3 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightSibbesen = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), alru2016, iter = 100, 
#                                               lower = c(a1 = 0.01, b1 = 0.2, b2 = 0.2),
#                                               start_lower = list(a1 = 2, b1 = 0.25, b2 = 0.25),
#                                               start_upper = list(a1 = 4, b1 = 0.30, b2 = 0.30), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeight = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # step factor
#alruDiameterFromHeightWeibull = nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -1000, b1 = 0.002, b2 = 0.95), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # NaN-inf
#alruDiameterFromHeightWeibull = nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, algorithm = "port",
#                                    lower = list(a1 = -1000, b1 = 0.001, b2 = 0.4), # b1 constraint prevents NaN-inf
#                                    start = list(a1 = -150, b1 = 0.02, b2 = 0.75), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # collapses to linear fit
#alruDiameterFromHeightWeibullForm = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 500))
#alruDiameterFromHeightWykoff = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                             lower = c(a1 = 1, b1 = 1/55, b2 = 0.2),
#                                             start_lower = list(a1 = 3, b1 = 0.2, b2 = 0.5), 
#                                             start_upper = list(a1 = 7, b1 = 0.4, b2 = 0.7), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightWykoffAat = nls_multstart(DBH ~ (a1 + a2 * tallerQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a2 = -1, b1 = -0.05, b2 = -1), 
#                                                start_upper = list(a1 = 1, a2 = 1, b1 = 0.05, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
#alruDiameterFromHeightWykoffAat = nls_multstart(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea + (a3 + a3p * isPlantation) * standQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a1p = -100, a2 = -1, a2p = -1, a3 = -1, a3p = -1, b1 = -1, b2 = -1, b2p = -1), 
#                                                start_upper = list(a1 = 1, a1p = 1, a2 = 1, a2p = 1, a3 = 1, a3p = 1, b1 = 1, b2 = 2, b2p = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
# coefficients for physiologically plausible curvature but poor accuracy
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -16.9, a1p = 2.676, b1 = 0.0307, b2 = 0.0304, b2p = 0.0723), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 210)) # b1p not significant, >200 iterations to converge to start point
#alruDiameterFromHeightChapmanRichardsAat = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -15.4, a2 = 0.234, a3 = 0.233, b1 = 0.0342, b2 = 0.263, b2p = 0.164), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 500)) # a1p, a2, a2p, a3, a3p not significant
#alruDiameterFromHeightChapmanRichardsPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -10.5, a1p = 2.93, a2 = 0.0030, a3 = -18.3, a4 = -0.095, a5 = -0.42, a6 = 0.137, b1 = 0.0307, b2 = 0.342, b2p = 0.0683), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 200)) # a2, a4, a5, a6, b1p not significant
#alruDiameterFromHeightChapmanRichardsRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -36.7, a1p = 30.7, a2 = 33.2, a2p = -40.8, b1 = 0.0307, b2 = 0.562, b2p = -0.360), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 200))
# coefficients for physically implausible accuracy
alruDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a1p not significant
alruDiameterFromHeightChapmanFormAat = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a1p not significant, all NaN-inf with a3 * standQuasiBasalArea
alruDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -32.1, a1p = 4.8, a2 = 1.72, a2p = -0.45, a3 = -1.46, a3p = 0.45, b1 = -0.040, b2 = 1.31), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b2p not significant
alruDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a4 = 29.8, a4p = -14.8, b1 = -0.0302, b2 = 1.30), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # NaN-inf
alruDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -54.2, a2 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a1p not significant
alruDiameterFromHeightChapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # b1p not significant
alruDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = list(maxiter = 50)) # a2p not significant
alruDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a2 = -0.0004, a3 = 7.37, a4 = 0.41, a5 = 0.52, a6 = 0.037, b1 = -0.031, b2 = 2.47, b2p = -1.26), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a2, a4, a6 not significant
alruDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016, start = list(a1 = 33.0, a1p = -7.77, a2 = -0.95, b1 = -0.043, b2 = 1.39), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a2p, b1p not significant
alruDiameterFromHeightGam = gam(DBH ~ s(TotalHt, by = as.factor(isPlantation)), data = alru2016)
alruDiameterFromHeightGamAat = gam(DBH ~ s(TotalHt, tallerQuasiBasalArea, standQuasiBasalArea, by = as.factor(isPlantation)), data = alru2016)
#alruDiameterFromHeightGamAatPhysio = gam(DBH ~ s(TotalHt, tallerQuasiBasalArea, standQuasiBasalArea, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = alru2016physio)
#alruDiameterFromHeightGamPhysio = gam(DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = alru2016physio)
alruDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), alru2016, weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightMichaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), alru2016, start = list(a1 = 100, a2 = 100, b1 = 1), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # collapses to linear with or without b1p
alruDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), alru2016, start = list(a1 = 3.8, a1p = -1.0, a2 = -0.12, a2p = -0.006), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), alru2016, weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation  + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, alru2016physio, start = list(a1 = 2.38, a1p = -0.70, a2 = -0.00002, a3 = 1.26, a4 = 0.078, a5 = 0.0075, a6 = 0.00017, b1 = 0.88), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a2, a6, b1p not significant
alruDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.55, a1p = -1.73, a2 = -0.98, a2p = 0.65, b1 = 0.86, b1p = 0.14), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.35, b1 = 1.37, b1p = -0.24, b2 = -0.033, b2p = 0.022), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a1p not significant
alruDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), alru2016, start = list(a1 = 0.00006, a2 = 0.04, b1 = 0.94, Ha = 100), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # poorly conditioned, singular gradient with Levenberg, nls(), or nlrob()
alruDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, alru2016, start = list(a1 = 29, a2 = 0.71, b1 = 0.0001, b2 = -0.65, b3 = 0.21), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = nls.control(maxiter = 500)) # nls() NaN-infinity even with parameters from nls_multstart(), nlrob() NaN-inf
alruDiameterFromHeightSibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))
alruDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5))  # a2 not significant
alruDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a2 = -0.00002, a3 = 0.231, a4 = 0.191, a5 = 0.018, a6 = 0.0007, b1 = 2.13, b2 = -0.163, b2p = 0.0197), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # a2, a6 not significant
alruDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.712, a1p = -0.341, a2 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5)) # neither a2 not significant
alruDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, start = list(a1 = -1000, b1 = 0.002, b2 = 0.95), weights = pmin(TotalHt^if_else(isPlantation, -0.8, -0.9), 0.5), control = gsl_nls_control(scale = "levenberg")) # collapses to linear, singular gradient with a1p, b1p, b2p, nlrob() NaN-inf
#confint(alruDiameterFromHeightRuark, level = 0.99)

alruDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", alruDiameterFromHeightChapmanForm, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AA+T", alruDiameterFromHeightChapmanFormAat, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BA+L", alruDiameterFromHeightChapmanFormBal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BA+L RelHt", alruDiameterFromHeightChapmanFormBalRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", alruDiameterFromHeightChapmanFormRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", alruDiameterFromHeightChapmanRichards, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AA+T", alruDiameterFromHeightChapmanRichardsAat, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", alruDiameterFromHeightChapmanRichardsPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", alruDiameterFromHeightChapmanRichardsRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightGam = get_dbh_error("GCV GAM", alruDiameterFromHeightGam, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightGamAat = get_dbh_error("GCV GAM AA+T", alruDiameterFromHeightGamAat, alru2016, alru2016natural, alru2016plantation)
#alruDiameterFromHeightGamAatPhysio = get_dbh_error("GCV GAM AA+T physio", alruDiameterFromHeightGamAatPhysio, alru2016physio, alru2016natural, alru2016plantation)
#alruDiameterFromHeightGamPhysio = get_dbh_error("GCV GAM physio", alruDiameterFromHeightGamPhysio, alru2016physio, alru2016natural, alru2016plantation)
#save(alruDiameterFromHeightGamAatPhysio, alruDiameterFromHeightGamPhysio, file = "trees/height-diameter/HtDia PSME spline DBH.rdata")
load("trees/height-diameter/HtDia PSME spline DBH.rdata")
alruDiameterFromHeightLinear = get_dbh_error("linear", alruDiameterFromHeightLinear, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", alruDiameterFromHeightMichaelisMentenForm, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightNaslund = get_dbh_error("Näslund", alruDiameterFromHeightNaslund, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightParabolic = get_dbh_error("parabolic", alruDiameterFromHeightParabolic, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPower = get_dbh_error("power", alruDiameterFromHeightPower, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPowerAat = get_dbh_error("power AA+T", alruDiameterFromHeightPowerAat, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPowerPhysio = get_dbh_error("power physio", alruDiameterFromHeightPowerPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", alruDiameterFromHeightPowerRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightRuark = get_dbh_error("Ruark", alruDiameterFromHeightRuark, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSchnute = get_dbh_error("Schnute", alruDiameterFromHeightSchnute, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", alruDiameterFromHeightSharmaParton, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", alruDiameterFromHeightSibbesenForm, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AA+T", alruDiameterFromHeightSibbesenFormAat, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", alruDiameterFromHeightSibbesenFormPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", alruDiameterFromHeightSibbesenFormRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightWeibull = get_dbh_error("Weibull", alruDiameterFromHeightWeibull, alru2016, alru2016natural, alru2016plantation)

alruDiameterFromHeightResults = bind_rows(as_row(alruDiameterFromHeightChapmanRichards),
                                          as_row(alruDiameterFromHeightChapmanRichardsAat),
                                          as_row(alruDiameterFromHeightChapmanRichardsPhysio),
                                          as_row(alruDiameterFromHeightChapmanRichardsRelHt),
                                          as_row(alruDiameterFromHeightChapmanForm),
                                          as_row(alruDiameterFromHeightChapmanFormAat),
                                          as_row(alruDiameterFromHeightChapmanFormBal),
                                          as_row(alruDiameterFromHeightChapmanFormBalRelHt),
                                          as_row(alruDiameterFromHeightChapmanFormRelHt),
                                          as_row(alruDiameterFromHeightGam),
                                          as_row(alruDiameterFromHeightGamAat),
                                          as_row(alruDiameterFromHeightGamAatPhysio),
                                          as_row(alruDiameterFromHeightGamPhysio),
                                          as_row(alruDiameterFromHeightLinear),
                                          as_row(alruDiameterFromHeightMichaelisMentenForm),
                                          as_row(alruDiameterFromHeightNaslund),
                                          as_row(alruDiameterFromHeightParabolic),
                                          as_row(alruDiameterFromHeightPower),
                                          as_row(alruDiameterFromHeightPowerAat),
                                          as_row(alruDiameterFromHeightPowerPhysio),
                                          as_row(alruDiameterFromHeightPowerRelHt),
                                          as_row(alruDiameterFromHeightRuark),
                                          as_row(alruDiameterFromHeightSchnute),
                                          as_row(alruDiameterFromHeightSharmaParton),
                                          as_row(alruDiameterFromHeightSibbesenForm),
                                          as_row(alruDiameterFromHeightSibbesenFormAat),
                                          as_row(alruDiameterFromHeightSibbesenFormPhysio),
                                          as_row(alruDiameterFromHeightSibbesenFormRelHt),
                                          as_row(alruDiameterFromHeightWeibull)) %>% 
  mutate(responseVariable = "DBH", species = "ALRU2", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(alruDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_smooth(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  #geom_line(aes(x = alruDiameterFromHeightChapmanForm$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards form", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichards$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsAat$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards AA+T", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsRelHt$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightKorf$fitted.values, y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightLinear$fitted.values, y = alru2016$TotalHt, color = "linear", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightMichaelisMentenForm$fitted.values, y = alru2016$TotalHt, color = "Michaelis-Menten form", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightNaslund$fitted.values, y = alru2016$TotalHt, color = "Näslund", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightParabolic$fitted.values, y = alru2016$TotalHt, color = "parabolic", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightPower$fitted.values, y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightRuark$fitted.values, y = alru2016$TotalHt, color = "Ruark", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightSchnute$fitted.values, y = alru2016$TotalHt, color = "Schnute", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightSibbesenForm$fitted.values, y = alru2016$TotalHt, color = "Sibbesen form", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightWeibull$fitted.values, y = alru2016$TotalHt, color = "Weibull", group = alru2016$isPlantation), alpha = 0.5) +
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
  #geom_line(aes(x = alruDiameterFromHeightKorf$fitted.values, y = TotalHt, color = "Korf")) +
  geom_line(aes(x = 1.246*exp(1*(TotalHt - 1.37)^0.4), y = TotalHt, color = "Korf manual")) +
  geom_line(aes(x = 0.155*exp(2.708*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf unweighted")) +
  geom_line(aes(x = 0.146*exp(2.732*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf weighted")) +
  geom_line(aes(x = 0.019*exp(4.691*(TotalHt - 1.37)^0.15), y = TotalHt, color = "Korf weighted to step factor")) +
  geom_line(aes(x = 0.051*exp(3.724*(TotalHt - 1.37)^0.18), y = TotalHt, color = "Korf GNLS")) +
  labs(x = "DBH, cm", y = "height, m", color = "regression\nform", fill = "stems\nmeasured") +
  scale_fill_viridis_c(trans = "log10")


ggplot() +
  geom_point(aes(x = alru2016$TotalHt, y = alruDiameterFromHeightRuark$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$TotalHt, y = alruDiameterFromHeightRuark$residuals), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  #geom_point(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  #geom_smooth(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
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

## collect model parameters
alruParameters = bind_rows(bind_rows(get_coefficients(alruHeightFromDiameterChapmanRichards),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(alruHeightFromDiameterCurtis),
                                     get_coefficients(alruHeightFromDiameterGam),
                                     get_coefficients(alruHeightFromDiameterGamBal),
                                     get_coefficients(alruHeightFromDiameterGamBalPhysio),
                                     get_coefficients(alruHeightFromDiameterGamPhysio),
                                     get_coefficients(alruHeightFromDiameterHossfeld),
                                     get_coefficients(alruHeightFromDiameterKorf),
                                     get_coefficients(alruHeightFromDiameterLinear),
                                     get_coefficients(alruHeightFromDiameterMichaelisMenten),
                                     get_coefficients(alruHeightFromDiameterParabolic),
                                     get_coefficients(alruHeightFromDiameterPower),
                                     get_coefficients(alruHeightFromDiameterProdan),
                                     get_coefficients(alruHeightFromDiameterRatkowsky),
                                     get_coefficients(alruHeightFromDiameterRichards),
                                     get_coefficients(alruHeightFromDiameterSharmaParton),
                                     get_coefficients(alruHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(alruHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(alruHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(alruHeightFromDiameterSharmaZhang),
                                     get_coefficients(alruHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(alruHeightFromDiameterSibbesen),
                                     get_coefficients(alruHeightFromDiameterWeibull),
                                     get_coefficients(alruHeightFromDiameterWeibullBal),
                                     get_coefficients(alruHeightFromDiameterWeibullBalRelHt),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(alruHeightFromDiameterChapmanRichardsBalGnls),
                                     get_coefficients(alruHeightFromDiameterSharmaPartonGnls),
                                     get_coefficients(alruHeightFromDiameterSharmaPartonBalGnls),
                                     get_coefficients(alruHeightFromDiameterSharmaZhangGnls),
                                     get_coefficients(alruHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(alruHeightFromDiameterWeibullGnls),
                                     get_coefficients(alruHeightFromDiameterWeibullBalGnls)) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(alruDiameterFromHeightChapmanRichards),
                                     get_coefficients(alruDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(alruDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(alruDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(alruDiameterFromHeightChapmanForm),
                                     get_coefficients(alruDiameterFromHeightChapmanFormAat),
                                     get_coefficients(alruDiameterFromHeightChapmanFormBal),
                                     get_coefficients(alruDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(alruDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(alruDiameterFromHeightGam),
                                     get_coefficients(alruDiameterFromHeightGamAat),
                                     get_coefficients(alruDiameterFromHeightGamAatPhysio),
                                     get_coefficients(alruDiameterFromHeightGamPhysio),
                                     get_coefficients(alruDiameterFromHeightLinear),
                                     get_coefficients(alruDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(alruDiameterFromHeightNaslund),
                                     get_coefficients(alruDiameterFromHeightParabolic),
                                     get_coefficients(alruDiameterFromHeightPower),
                                     get_coefficients(alruDiameterFromHeightPowerAat),
                                     get_coefficients(alruDiameterFromHeightPowerPhysio),
                                     get_coefficients(alruDiameterFromHeightPowerRelHt),
                                     get_coefficients(alruDiameterFromHeightRuark),
                                     get_coefficients(alruDiameterFromHeightSchnute),
                                     get_coefficients(alruDiameterFromHeightSharmaParton),
                                     get_coefficients(alruDiameterFromHeightSibbesenForm),
                                     get_coefficients(alruDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(alruDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(alruDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(alruDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "ALRU2",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, fitting, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

## robust regression
# adaptive weighting by nlrob()
ggplot() +
  geom_histogram(aes(x = alruHeightFromDiameterChapmanRichards$rweights), binwidth = 0.01) +
  labs(x = "nlrob() adaptive weight reduction", y = "number of red alder heights imputed") +
  scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), minor_breaks = c(0.5, 3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 3000, 4000, 6000, 7000, 8000, 9000), trans = scales::pseudo_log_trans(base = 10)) +
ggplot() +
  geom_histogram(aes(x = alruDiameterFromHeightChapmanRichards$rweights), binwidth = 0.01) +
  labs(x = "nlrob() adaptive weight reduction", y = "number of red alder diameters imputed") +
  scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), minor_breaks = c(0.5, 3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 3000, 4000, 6000, 7000, 8000, 9000), trans = scales::pseudo_log_trans(base = 10))



## basal area from height
#alruBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), alru2016, start = list(a1 = 500, b1 = 0.0002, b2 = 1.9), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # nlrob() step factor
alruBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), alru2016, start = list(a1 = 682, b1 = 0.00004, b2 = 1.9, b2p = -0.11), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # a1p, b1p not significant, nlrob() step factor
alruBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^b1, alru2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.00006, b1 = 1.91), weights = pmin(1/basalArea, 1E4)) # b1p not significant
#confint2(alruBasalAreaFromHeightKorf, level = 0.99)

alruBasalAreaFromHeightKorf$fitted.values = predict(alruBasalAreaFromHeightKorf, alru2016)
alruBasalAreaFromHeightKorf$residuals = alruBasalAreaFromHeightKorf$fitted.values - alru2016$basalArea
alruBasalAreaFromHeightPower$fitted.values = predict(alruBasalAreaFromHeightPower, alru2016)
alruBasalAreaFromHeightPower$residuals = alruBasalAreaFromHeightPower$fitted.values - alru2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(alruBasalAreaFromHeightKorf), 100^2 * mean(alruBasalAreaFromHeightKorf$residuals), mean(abs(alruBasalAreaFromHeightKorf$residuals)), 1 - sum(alruBasalAreaFromHeightKorf$residuals^2) / sum((alru2016$basalArea - mean(alru2016$basalArea)^2)),
        "power", AIC(alruBasalAreaFromHeightPower), 100^2 * mean(alruBasalAreaFromHeightPower$residuals), mean(abs(alruBasalAreaFromHeightPower$residuals)), 1 - sum(alruBasalAreaFromHeightPower$residuals^2) / sum((alru2016$basalArea - mean(alru2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(alru2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = alruBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = alruBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "red alder height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
