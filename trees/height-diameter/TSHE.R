# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western hemlock height-diameter regression form sweep
tshe2016 = trees2016 %>% filter(Species == "WH", isLiveUnbroken, TotalHt > 0) %>% # live western hemlocks measured for height
  mutate(dbhWeight = pmin(1/(1.32*DBH^0.80), 5),
         heightWeight = pmin(1/(1.07*(TotalHt - 1.37)^1.50), 5))
tshe2016physio = tshe2016 %>% filter(is.na(elevation) == FALSE)
tshe2016gamConstraint = c(DBH = -1.2994/0.6005, TotalHt = 1.37, standBasalAreaPerHectare = median(tshe2016$standBasalAreaPerHectare), basalAreaLarger = median(tshe2016$basalAreaLarger), standBasalAreaApprox = median(tshe2016$standBasalAreaApprox), tallerApproxBasalArea = median(tshe2016$tallerApproxBasalArea), elevation = median(tshe2016physio$elevation), slope = median(tshe2016physio$slope), aspect = median(tshe2016physio$aspect), topographicShelterIndex = median(tshe2016physio$topographicShelterIndex), relativeHeight = median(tshe2016$relativeHeight)) # point constraint for mgcv::s()

tshe2016defaultWeight = tshe2016 %>% mutate(dbhWeight = pmin(1/DBH, 5),
                                            heightWeight = pmin(1/TotalHt, 5))
tshe2016defaultWeightPhysio = tshe2016defaultWeight %>% filter(is.na(elevation) == FALSE)

tsheHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), tshe2016))
tsheHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), tshe2016)

tsheHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048))
tsheHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 60.5, a2 = 0.026, a2p = 0.96, a3 = 0, a3p = 0, a4 = -0.012, a5 = -0.19, a6 = 0.44, a7 = 0.80, a8 = 0.15, b1 = -0.022, b2 = 1.31)) # a1p, a6, b1p, b2p not significant
tsheHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a9 = 63.0, a9p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01)) # a3, a3p, b1p not significant
tsheHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58.5, a1p = -3.5, a4 = -0.015, a5 = -8.07, a7 = 0.718, a8 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12)) # a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)) # b1 not significant
tsheHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
tsheHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016, start = list(a1 = 200, b1 = -7.2, b2 = -0.33)) # a1p, b1p, b2p not significant
tsheHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)) # a3p not significant
tsheHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
tsheHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)) # b2p not significant
tsheHeightFromDiameter$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)) # dp not significant
tsheHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 30.7, b1 = 0.15, b1p = -0.068, b2 = -0.034, b3 = -0.22, b3p = 0.184, b4 = 1.26)) # a1p, b2p, b3p not significant
tsheHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 46.7, a1p = -12.6, b1 = 0.054, b2 = -0.021, b2p = -0.013, b3 = -0.054, b4 = 1.26)) # b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)) # b1p, b2p, b4p not significant
tsheHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.0005, a2p = 0.0191, b1 = 0.100, b1p = -0.207, b2 = -0.020, b3 = -0.029, b4 = 1.24)) # a3, b1, b2p, b3p, b4p not significant
tsheHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
tsheHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029))
tsheHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212)) # a1p, b1p, b2p not significant
tsheHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 16.0, a1p = -9.7, a2 = 0.018, a2p = 0.63, a3 = 0.28, a9 = 73.5, b1 = -0.015, b2 = 0.76), control = list(maxiter = 200)) # a3p, a4p, b1p, b2p not significant

tsheHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048)))
tsheHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 60.5, a2 = 0.026, a2p = 0.96, a3 = 0, a3p = 0, a4 = -0.012, a5 = -0.19, a6 = 0.44, a7 = 0.80, a8 = 0.15, b1 = -0.022, b2 = 1.31)) # a1p, a6, b1p, b2p not significant
tsheHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a9 = 63.0, a9p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01)) # a3, a3p, b1p not significant
tsheHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58.5, a1p = -3.5, a4 = -0.015, a5 = -8.07, a7 = 0.718, a8 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12)) # a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)) # b1 not significant
tsheHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
tsheHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016, start = list(a1 = 200, b1 = -7.2, b2 = -0.33)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)) # a3p not significant
tsheHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
tsheHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)) # b2p not significant
tsheHeightFromDiameterNlrob$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)) # dp not significant
tsheHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 30.7, b1 = 0.15, b1p = -0.068, b2 = -0.034, b3 = -0.22, b3p = 0.184, b4 = 1.26)) # a1p, b2p, b3p not significant
tsheHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 46.7, a1p = -12.6, b1 = 0.054, b2 = -0.021, b2p = -0.013, b3 = -0.054, b4 = 1.26)) # b1p, b3p, b4p not significant
tsheHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)) # b1p, b2p, b4p not significant
tsheHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.0005, a2p = 0.0191, b1 = 0.100, b1p = -0.207, b2 = -0.020, b3 = -0.029, b4 = 1.24)) # a3, b1, b2p, b3p, b4p not significant
tsheHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
tsheHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029))
tsheHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 16.0, a1p = -9.7, a2 = 0.018, a2p = 0.63, a3 = 0.28, a9 = 73.5, b1 = -0.015, b2 = 0.76), control = list(maxiter = 200)) # a3p, a4p, b1p, b2p not significant
#confint_nlrob(tsheHeightFromDiameterNlrob$sharmaPartonBalPhysio, level = 0.99)

tsheHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeight, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048)))
tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeight, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, tshe2016defaultWeightPhysio, start = list(a1 = 60.5, a2 = 0.026, a2p = 0.96, a3 = 0, a3p = 0, a4 = -0.012, a5 = -0.19, a6 = 0.44, a7 = 0.80, a8 = 0.15, b1 = -0.022, b2 = 1.31)) # a1p, a6, b1p, b2p not significant
tsheHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeight, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a9 = 63.0, a9p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01)) # a3, a3p, b1p not significant
tsheHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016defaultWeightPhysio, start = list(a1 = 58.5, a1p = -3.5, a4 = -0.015, a5 = -8.07, a7 = 0.718, a8 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12)) # a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)) # b1 not significant
tsheHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016defaultWeight, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
tsheHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016defaultWeight, start = list(a1 = 200, b1 = -7.2, b2 = -0.33)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016defaultWeight, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016defaultWeight, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)) # a3p not significant
tsheHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
tsheHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016defaultWeight, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)) # b2p not significant
tsheHeightFromDiameterGslNlsDefault$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016defaultWeight, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)) # dp not significant
tsheHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 30.7, b1 = 0.15, b1p = -0.068, b2 = -0.034, b3 = -0.22, b3p = 0.184, b4 = 1.26)) # a1p, b2p, b3p not significant
tsheHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 46.7, a1p = -12.6, b1 = 0.054, b2 = -0.021, b2p = -0.013, b3 = -0.054, b4 = 1.26)) # b1p, b3p, b4p not significant
tsheHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016defaultWeightPhysio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)) # b1p, b2p, b4p not significant
tsheHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016defaultWeight, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.0005, a2p = 0.0191, b1 = 0.100, b1p = -0.207, b2 = -0.020, b3 = -0.029, b4 = 1.24)) # a3, b1, b2p, b3p, b4p not significant
tsheHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016defaultWeight, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
tsheHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016defaultWeight, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029))
tsheHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016defaultWeight, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), tshe2016defaultWeight, start = list(a1 = 16.0, a1p = -9.7, a2 = 0.018, a2p = 0.63, a3 = 0.28, a9 = 73.5, b1 = -0.015, b2 = 0.76), control = list(maxiter = 200)) # a3p, a4p, b1p, b2p not significant

tsheHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = tshe2016gamConstraint), data = tshe2016)
tsheHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 19, pc = tshe2016gamConstraint), data = tshe2016)
tsheHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 90, pc = tshe2016gamConstraint), data = tshe2016physio) # ~2 minutes to fit+evaluate with all predictors at minimum (Zen 3 @ 3.4 GHz, k = 495, edf < 290) -> eliminate cos(aspect) (k = 330, edf < 240) AIC 11023: 10733 without BA, 10887 without BAL, 10729 without elevation, 10714 without slope, 10679 without sin(aspect), 10669 without topographic shelter -> eliminate sin(aspect)
tsheHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio)

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
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$richards), color = "unified Richards", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sibbesen), color = "Sibbesen", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibull), color = "Weibull", group = tshe2016$isPlantation)) +
    annotate("text", x = 0, y = 70, label = "western hemlock, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 70)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}


## western hemlock height-diameter GNLS regressions
if (htDiaOptions$fitGnls)
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


## western hemlock diameter-height regressions
tsheDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), tshe2016)) # isPlantation*(TotalHt - 1.37) not significant(p = 0.036)
tsheDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016)

tsheDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 32, b1 = 0.034, b2 = 0.73)) # no significant plantation effects
tsheDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 35, a2 = -0.04, b1 = 0.034, b2 = 0.73))
tsheDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 35, a2 = -0.044, b1 = 0.033, b2 = 0.74), control = gsl_nls_control(maxiter = 250), significant = FALSE) # a1p, a2, a3, b1p not significant, a3 + b2p step size
tsheDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.3, a2 = -0.003, a9 = -0.08, b1 = 2.6, b2 = 0.22), control = gsl_nls_control(maxiter = 250), significant = FALSE) # a1, a2, a3, a9 not significant
tsheDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.4, a9 = 0.06, b1 = 2.2, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE) # a9 not significant, potentially >500 iterations with nlrob()
tsheDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -95, b1 = 0.027, b2 = 0.81), control = gsl_nls_control(maxiter = 250)) # a1p, b1p not significant
tsheDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -100, a2 = 0.3, b1 = 0.03, b2 = 0.75), significant = FALSE) # a1p, a2, b1p not significant
tsheDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -85, a5 = -20, b1 = 0.03, b2 = 0.8), significant = FALSE) # a1p, a4, a5, a6, a7, a8, b1p not significant
tsheDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -60, a9 = -14, b1 = 0.04, b2 = 0.72), significant = FALSE) # a1p, a9, b1p, b2p not significant, a1+a9-b1 not mutually significant
tsheDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83)) # a1p, a2p, b1p not significant
tsheDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
tsheDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04)) # no significant plantation effects
#tsheDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047)) # a1p not significant
#tsheDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04)) # a4, a6, a7, a8 not significant
#tsheDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
tsheDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018)) # a1p, b1p, b2p not significant
tsheDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 2.9, a2 = -0.004, a2p = -0.01, b1 = 0.75, b2 = 0.01)) # a2, a3, a3p, b1p, b2p not significant
tsheDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013)) # a1p, a4, a5, a6, a7, a8, b1p, b2p not significant
tsheDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01)) # a9p, b2p not significant
tsheDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012))
#tsheDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32)) # converges from red alder values but fails to reconverge (singular gradient), NaN-inf or singular gradient with fit_gsl_nls()
#tsheDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016, start = list(a1 = 1.4, b1 = 1.0, b2 = 0.001, b3 = 0.3, b4 = 0.9), control = gsl_nls_control(maxiter = 125, xtol = 1E-6), significant = FALSE) # unreliable convergence and NaN-inf with fit_gsl_nls() due to parameter evaporation and collapse to power form, NaN-inf or singular gradient with nlrob(), singular gradient with nls() even with nls_multstart() parameters
tsheDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057)) # no significant plantation effects
tsheDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578)) # no significant plantation effects
tsheDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.115, a5 = 0.486, b1 = 0.736, b2 = 0.0593)) # a1p, a4, a6, a7, a8, b1p, b2p not significant
tsheDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056))
tsheDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.6, a5 = 0.7, a9 = 0.4, b1 = 0.5, b2 = 0.13))
tsheDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -225, b1 = 0.011, b2 = 0.82), control = gsl_nls_control(maxiter = 250)) # a1p, b1p, b2p not significant

tsheDiameterFromHeightNlrob = list(chapmanReplace = fit_nlrob("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 32, b1 = 0.034, b2 = 0.73)))
tsheDiameterFromHeightNlrob$chapmanReplaceAbat = fit_nlrob("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 35, a2 = -0.04, b1 = 0.034, b2 = 0.73))
tsheDiameterFromHeightNlrob$chapmanReplaceBal = fit_nlrob("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 35, a2 = -0.044, b1 = 0.033, b2 = 0.74), control = nls.control(maxiter = 250), significant = FALSE)
tsheDiameterFromHeightNlrob$chapmanReplaceBalRelHt = fit_nlrob("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.3, a2 = -0.003, a9 = -0.08, b1 = 2.6, b2 = 0.22), control = nls.control(maxiter = 250), significant = FALSE)
#tsheDiameterFromHeightNlrob$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 0.4, a9 = 0.06, b1 = 2.2, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE)
tsheDiameterFromHeightNlrob$chapmanRichards = fit_nlrob("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -80, b1 = 0.035, b2 = 0.77), control = nls.control(maxiter = 250)) # a1p, b1p not significant
tsheDiameterFromHeightNlrob$chapmanRichardsAbat = fit_nlrob("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -85, a2 = 0.15, b1 = 0.033, b2 = 0.78), significant = FALSE) # a1p, a2, b1p not significant
tsheDiameterFromHeightNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -70, a5 = -15, b1 = 0.033, b2 = 0.77)) # a1p, a4, b1p, b2p not significant
tsheDiameterFromHeightNlrob$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -60, a9 = -14, b1 = 0.04, b2 = 0.72), significant = FALSE)
tsheDiameterFromHeightNlrob$michaelisMentenReplace = fit_nlrob("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83))
tsheDiameterFromHeightNlrob$naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
tsheDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04))
#tsheDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047))
#tsheDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04))
#tsheDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
tsheDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018))
tsheDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016, start = list(a1 = 2.9, a2 = -0.004, a2p = -0.01, b1 = 0.75, b2 = 0.01))
tsheDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013))
tsheDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01))
tsheDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016physio, start = list(a1 = 2.4, a5 = 0.6, a9 = 0.5, b1 = 0.65, b2 = 0.015))
tsheDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057))
tsheDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578))
tsheDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.8, a5 = 0.55, b1 = 0.5, b2 = 0.15))
tsheDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056))
tsheDiameterFromHeightNlrob$sibbesenReplaceRelHtPhysio = fit_nlrob("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016physio, start = list(a1 = 2.6, a5 = 0.7, a9 = 0.45, b1 = 0.5, b2 = 0.14))
tsheDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -110, b1 = 0.068, b2 = 0.54), control = nls.control(maxiter = 500))
#lapply(tsheDiameterFromHeight$powerPhysio, confint_nlrob)

tsheDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016defaultWeight, start = list(a1 = 32, b1 = 0.034, b2 = 0.73)))
tsheDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016defaultWeight, start = list(a1 = 35, a2 = -0.04, b1 = 0.034, b2 = 0.73))
tsheDiameterFromHeightGslNlsDefault$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016defaultWeight, start = list(a1 = 35, a2 = -0.044, b1 = 0.033, b2 = 0.74), control = gsl_nls_control(maxiter = 250), significant = FALSE)
tsheDiameterFromHeightGslNlsDefault$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016defaultWeight, start = list(a1 = 0.3, a2 = -0.003, a9 = -0.08, b1 = 2.6, b2 = 0.22), control = gsl_nls_control(maxiter = 250), significant = FALSE)
tsheDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016defaultWeight, start = list(a1 = 0.4, a9 = 0.06, b1 = 2.2, b2 = 0.2), control = gsl_nls_control(maxiter = 500), significant = FALSE)
tsheDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -170, b1 = 0.01, b2 = 0.93), control = gsl_nls_control(maxiter = 250)) # a1p, b1p not significant
tsheDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -200, a2 = 0.3, b1 = 0.01, b2 = 0.93)) # a1p, b1p not significant
tsheDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeightPhysio, start = list(a1 = -200, a5 = -45, b1 = 0.01, b2 = 0.93), significant = FALSE) # a1p, a4, a5, b1p not significant
tsheDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), tshe2016defaultWeight, start = list(a1 = -60, a9 = -14, b1 = 0.04, b2 = 0.72), significant = FALSE)
tsheDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016defaultWeight, start = list(a1 = 153, a2 = 68, b1 = 0.83))
tsheDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
tsheDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016defaultWeight, start = list(a1 = 1.52, b1 = 1.04))
#tsheDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016defaultWeight, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047))
#tsheDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016defaultWeightPhysio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04))
#tsheDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016defaultWeight, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
tsheDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.8, b1 = 0.7, b2 = 0.018))
tsheDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + (a2 + a2p*isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.9, a2 = -0.004, a2p = -0.01, b1 = 0.75, b2 = 0.01))
tsheDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 2.4, a5 = 0.6, b1 = 0.75, b2 = 0.013))
tsheDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9*relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp(b2 * (TotalHt - 1.37)), tshe2016defaultWeight, start = list(a1 = 2.7, a9 = 1.0, b1 = 0.74, b1p = -0.033, b2 = 0.01))
tsheDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9*relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2*(TotalHt - 1.37)), tshe2016defaultWeightPhysio, start = list(a1 = 2.4, a5 = 0.65, a9 = 0.5, b1 = 0.7, b2 = 0.012))
#tsheDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016defaultWeight, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32))
#tsheDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016defaultWeight, start = list(a1 = 1.4, b1 = 1.0, b2 = 0.001, b3 = 0.3, b4 = 0.9), control = gsl_nls_control(maxiter = 125, xtol = 1E-6), significant = FALSE)
tsheDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057))
tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578))
tsheDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.8, a5 = 0.5, b1 = 0.5, b2 = 0.15))
tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeight, start = list(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056))
tsheDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016defaultWeightPhysio, start = list(a1 = 2.6, a5 = 0.7, a9 = 0.4, b1 = 0.5, b2 = 0.13))
tsheDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016defaultWeight, start = list(a1 = -120, b1 = 0.06, b2 = 0.55), control = gsl_nls_control(maxiter = 250))

tsheDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = tshe2016gamConstraint), data = tshe2016)
tsheDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = tshe2016gamConstraint), data = tshe2016)
tsheDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio, nthreads = 4) # k = 495, ef < 210 with all predictors AIC 13838: 13828 without AAT, 13826 without ABA, 13825 without elevation, 13793 without slope, 13848 without sin(aspect), 13896 without cos(aspect), 13772 without topographic shelter -> eliminate topographic shelter (k = 330, edf < 185) AIC 13772: 13654 without AAT, 13665 without ABA, 13649 without elevation, 13649 without slope, 13635 without sin(aspect), 13656 without cos(aspect) -> eliminate sin(aspect)
tsheDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio)
tsheDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 15, pc = tshe2016gamConstraint), data = tshe2016, nthreads = 4)
tsheDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 60, pc = tshe2016gamConstraint), data = tshe2016physio, nthreads = 6)

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
    geom_line(aes(x = -1/0.0005*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute inverse"), alpha = 0.5) +
    geom_line(aes(x = 1.3*(TotalHt - 1.37)^1*exp(0.003*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 81, label = "western hemlock, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## collect model parameters
#if (exists("tsheHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/TSHE GNLS.rdata") }
tsheCoefficients = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_coefficients)),
                                       bind_rows(lapply(tsheHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                       bind_rows(lapply(tsheHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                       #bind_rows(lapply(tsheHeightFromDiameterGnls, get_model_coefficients)) %>%
                               mutate(responseVariable = "height"),
                             bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_coefficients)),
                                       bind_rows(lapply(tsheDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                       bind_rows(lapply(tsheDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                       #get_model_coefficients(tsheDiameterFromHeight$schnute),
                                       #get_model_coefficients(tsheDiameterFromHeight$sharmaParton),
                               mutate(responseVariable = "DBH")) %>%
  mutate(species = "TSHE")
tsheResults = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_stats)),
                                  bind_rows(lapply(tsheHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                  bind_rows(lapply(tsheHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                                  #bind_rows(lapply(tsheHeightFromDiameterGnls, get_model_stats))) %>%
                          mutate(responseVariable = "height"),
                        bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_stats)),
                                  get_model_stats(name = "Schnute inverse", fitting = "gsl_nls"),
                                  get_model_stats(name = "modified Sharma-Parton", fitting = "gsl_nls"),
                                  bind_rows(lapply(tsheDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                  bind_rows(lapply(tsheDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                          mutate(responseVariable = "DBH")) %>%
  mutate(species = "TSHE")

save(file = "trees/height-diameter/data/TSHE results.Rdata", tsheCoefficients, tsheResults)
if (htDiaOptions$folds * htDiaOptions$repetitions <= htDiaOptions$retainModelThreshold)
{
  save(file = "trees/height-diameter/data/TSHE models and stats.Rdata", 
       tsheHeightFromDiameter, tsheDiameterFromHeightGslNlsDefault, tsheHeightFromDiameterNlrob,
       tsheDiameterFromHeight, tsheDiameterFromHeightGslNlsDefault, tsheDiameterFromHeightNlrob)
} else {
  save(file = "trees/height-diameter/data/TSHE stats.Rdata", 
       tsheHeightFromDiameter, tsheDiameterFromHeightGslNlsDefault, tsheHeightFromDiameterNlrob,
       tsheDiameterFromHeight, tsheDiameterFromHeightGslNlsDefault, tsheDiameterFromHeightNlrob)
}


## preferred forms identified (results.R, Figure 5)
tsheHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1))
tsheHeightFromDiameterPreferred$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26), folds = 1, repetitions = 1)

tsheDiameterFromHeightPreferred = list(chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a9 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912), folds = 1, repetitions = 1))
tsheDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1)
tsheDiameterFromHeightPreferred$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1)
tsheDiameterFromHeightPreferred$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057), folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/TSHE preferred models.Rdata", tsheHeightFromDiameterPreferred, tsheDiameterFromHeightPreferred)


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  #tsheBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ a1*exp(b1*(imputedHeight - 1.37)^b2) - 1, tshe2016, start = list(a1 = 1, b1 = 0.00009, b2 = 2.2), weights = pmin(1/basalArea, 1E4))
  tsheBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 5.5, b1 = 0.00003, b2 = 2.09, b2p = -0.060), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor on b1p without a1p
  tsheBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 11/7 * 0.25 * pi * 0.01^2, a1p = -0.00005, b1 = 2.18, b1p = -0.178), weights = pmin(1/basalArea, 1E4))
  #confint_nlrob(tsheBasalAreaFromHeightKorf, level = 0.99, weights = pmin(1/tshe2016$basalArea, 1E4))
  
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
