# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## bigleaf maple height-diameter regression form sweep
# preferred forms: Sharma-Parton BA+L, Ratkowsky, Sharma-Parton, Hossfeld, Weibull, Chapman-Richards
#acmaHeightFromDiameterChapmanRichardsPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 31.9, a2 = -0.0044, a3 = -8.22, a4 = 0.198, a5 = 1.387, a6 = 0.014, b1 = -0.031, b2 = 1.02, b2p = -0.108), weights = pmin(DBH^-0.8, 1)) # a1p, a2, a4, b1p not significant
#acmaHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + a1*DBH / (a2 + a2p * isPlantation + DBH), acma2016, start = list(a1 = 39.5, a2 = 47.8, a2p = -9.11), weights = pmin(DBH^-0.8, 1))
#acmaHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), acma2016, start = list(Ha = 23.9, Hap = -1.1, d = 0.67, dp = -0.031, kU = 0.023, kUp = 0.008), weights = pmin(DBH^-0.8, 1))
acma2016 = trees2016 %>% filter(Species == "BM", isLiveUnbroken, TotalHt > 0) # live bigleaf maples measured for height
acma2016natural = acma2016 %>% filter(isPlantation == FALSE)
acma2016physio = acma2016 %>% filter(is.na(elevation) == FALSE)
acma2016plantation = acma2016 %>% filter(isPlantation)
acma2016plantationPhysio = acma2016physio %>% filter(isPlantation)

acmaHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^b2, acma2016, start = list(a1 = 26.8, a1p = 4.20, b1 = -0.026, b2 = 0.927), weights = pmin(DBH^-0.8, 1)) # b1p, b2p not significant
acmaHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), weights = pmin(DBH^-0.8, 1))
acmaHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016physio, start = list(a1 = 33.2, a2 = 0.156, a3 = 0, a3p = 0, a4 = -0.009, a5 = -0.130, a6 = 1.872, a7 = 0.276, a8 = 0, b1 = -0.024, b2 = 0.988, b2p = -0.129), weights = pmin(DBH^-0.8, 1)) # a1p, a2p, a3, a4, a7, a8, b1p not significant
acmaHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = -2.4, a1p = 1.51, a2 = -0.023, a2p = 0.159, a3 = 0.045, a4 = 59.4, a4p = -25.9, b1 = -0.033, b2 = 0.018, b2p = 0.406), weights = pmin(DBH^-0.8, 1)) # a2, a3p not significant, NaN-inf with b1p
acmaHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 32.5, a4 = -0.0052, a5 = -8.37, a6 = 0.187, a7 = 1.557, a8 = 0.014, b1 = -0.031, b2 = 1.01, b2p = -0.1-4), weights = pmin(DBH^-0.8, 1)) # a1p, a2, a4, b1p not significant
acmaHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^b1, acma2016, start = list(a1 = 1.24, a1p = 0.23, b1 = 0.28), weights = pmin(DBH^-0.8, 1)) # b1p not significant
acmaHeightFromDiameterGam = gam(TotalHt ~ s(DBH, by = as.factor(isPlantation)), data = acma2016)
acmaHeightFromDiameterGamBal = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation)), data = acma2016)
#acmaHeightFromDiameterGamBalPhysio = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = acma2016physio) # insufficient data
acmaHeightFromDiameterGamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = acma2016physio)
acmaHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1 * DBH^b2), acma2016, start = list(a1 = 40.1, a1p = 6.30, b1 = 44.9, b2 = -0.97), weights = pmin(DBH^-0.8, 1)) # b1p, b2p not significant
acmaHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), acma2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^-0.8, 1))
acmaHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), acma2016, offset = breastHeight, weights = pmin(DBH^-0.8, 1))
acmaHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), weights = pmin(DBH^-0.8, 1)) # a1p, b1p not significant
acmaHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), acma2016, offset = breastHeight, weights = pmin(DBH^-0.8, 1))
acmaHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1 * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), acma2016, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19), weights = pmin(DBH^-0.8, 1)) # a1p, a3p not significant
acmaHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, acma2016, start = list(a1 = 1.11, a1p = 0.21, b1 = 0.75), weights = pmin(DBH^-0.8, 1)) # b1p not significant
acmaHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1/(DBH + b2 + b2p * isPlantation)), acma2016, start = list(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62), weights = pmin(DBH^-0.8, 1)) # b1p not significant
acmaHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), acma2016, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064), weights = pmin(DBH^-0.8, 1)) # Hap, dp not significant
acmaHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 5.86, a1p = 39.2, b1 = 0.36, b1p = -0.502, b2 = -0.080, b3 = -0.398, b4 = 0.979, b4p = -0.116), weights = pmin(DBH^-0.8, 1)) # b2p, b3p not significant
acmaHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 46.0, b1 = -0.12, b2 = -0.054, b3 = -0.375, b4 = 0.960, b4p = -0.102), weights = pmin(DBH^-0.8, 1)) # a1p, b1p, b2p, b3p not significant
acmaHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), acma2016physio, start = list(a1 = 15.9, a4 = -0.002, a5 = 0, a6 = 0.49, a7 = -0.088, a8 = 0.0019, b1 = 0.18, b1p = -0.056, b2 = -0.038, b3 = -0.096, b4 = 1.02, b4p = -0.11), weights = pmin(DBH^-0.8, 1)) # a1p, a4, a7, a8, b1, b2p, b3p not significant
acmaHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^(b4 + b4p * isPlantation), acma2016physio, start = list(a1 = 20.1, a4 = -0.002, a5 = 0, a6 = 0.54, a7 = 0.14, a8 = 0.028, b1 = 0.14, b1p = 0.035, b2 = -0.041, b3 = -0.094, b4 = 1.05, b4p = -0.14), weights = pmin(DBH^-0.8, 1)) # a4, a7, a8, b1, b2p, b3p not significant
acmaHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 17.1, a1p = 1.36, b1 = 0.112, b2 = -0.044, b3 = -0.053, b4 = 1.038, b4p = -0.118), weights = pmin(DBH^-0.8, 1)) # b1p, b2p, b3p not significant
acmaHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 22.0, a2 = 0.002, b1 = 0.050, b2 = -0.044, b3 = -0.064, b4 = 1.070, b4p = -0.197), weights = pmin(DBH^-0.8, 1)) # a1p, b1p, a2, a2p, b2p, b3p not significant
acmaHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^b2), acma2016, start = list(a1 = 0.752, a1p = 0.120, b1 = 1.180, b2 = -0.087), weights = pmin(DBH^-0.8, 1)) # b1p, b2p not significant
acmaHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^b2)), acma2016, start = list(a1 = 27.3, a1p = 4.27, b1 = -0.033, b2 = 0.94), weights = pmin(DBH^-0.8, 1)) # b1p, b2p not significant
acmaHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = 26.6, a2 = 0.095, a3 = 0.073, b1 = -0.0267, b1p = -0.0092, b2 = 0.930), weights = pmin(DBH^-0.8, 1)) # a1p, a2p, a3p, b2p not significant
acmaHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = -29.1, a2 = 1.07, a3 = 0.65, a4 = 175, b1 = -0.020, b1p = 0.0006, b2 = 0.42), weights = pmin(DBH^-0.8, 1), control = nls.control(maxiter = 50)) # a3, a4p, b1p not significant, step factor with a1p, a2p, a3p, b2p. b3p
#confint_nlrob(acmaHeightFromDiameterSharmaPartonBalPhysio, level = 0.99, weights = pmin(acma2016physio$DBH^-0.8, 1))

acmaHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", acmaHeightFromDiameterChapmanRichards, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BA+L", acmaHeightFromDiameterChapmanRichardsBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BA+L physio", acmaHeightFromDiameterChapmanRichardsBalPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BA+L RelHt", acmaHeightFromDiameterChapmanRichardsBalRelHt, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", acmaHeightFromDiameterChapmanRichardsPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterCurtis = get_height_error("Curtis", acmaHeightFromDiameterCurtis, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterGam = get_height_error("GCV GAM", acmaHeightFromDiameterGam, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterGamBal = get_height_error("GCV GAM BA+L", acmaHeightFromDiameterGamBal, acma2016, acma2016natural, acma2016plantation)
#acmaHeightFromDiameterGamBalPhysio = get_height_error("GCV GAM BA+L physio", acmaHeightFromDiameterGamBalPhysio, acma2016physio, acma2016natural, acma2016plantation) # slow
acmaHeightFromDiameterGamPhysio = get_height_error("GCV GAM physio", acmaHeightFromDiameterGamPhysio, acma2016physio, acma2016natural, acma2016plantation)
#save(acmaHeightFromDiameterGamBalPhysio, acmaHeightFromDiameterGamPhysio, file = "trees/height-diameter/HtDia ACMA3 spline height.rdata")
#load("trees/height-diameter/HtDia ACMA3 spline height.rdata")
acmaHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", acmaHeightFromDiameterHossfeld, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterKorf = get_height_error("Korf", acmaHeightFromDiameterKorf, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterLinear = get_height_error("linear", acmaHeightFromDiameterLinear, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", acmaHeightFromDiameterMichaelisMenten, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterParabolic = get_height_error("parabolic", acmaHeightFromDiameterParabolic, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterProdan = get_height_error("Prodan", acmaHeightFromDiameterProdan, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterPower = get_height_error("power", acmaHeightFromDiameterPower, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", acmaHeightFromDiameterRatkowsky, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterRichards = get_height_error("unified Richards", acmaHeightFromDiameterRichards, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", acmaHeightFromDiameterSharmaParton, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BA+L", acmaHeightFromDiameterSharmaPartonBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BA+L physio", acmaHeightFromDiameterSharmaPartonBalPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", acmaHeightFromDiameterSharmaPartonPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", acmaHeightFromDiameterSharmaZhang, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BA+L", acmaHeightFromDiameterSharmaZhangBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSibbesen = get_height_error("Sibbesen", acmaHeightFromDiameterSibbesen, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibull = get_height_error("Weibull", acmaHeightFromDiameterWeibull, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBal = get_height_error("Weibull BA+L", acmaHeightFromDiameterWeibullBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BA+L RelHt", acmaHeightFromDiameterWeibullBalRelHt, acma2016, acma2016natural, acma2016plantation)

acmaHeightFromDiameterResults = bind_rows(as_row(acmaHeightFromDiameterChapmanRichards),
                                          as_row(acmaHeightFromDiameterChapmanRichardsBal),
                                          as_row(acmaHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row(acmaHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(acmaHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(acmaHeightFromDiameterCurtis),
                                          as_row(acmaHeightFromDiameterGam),
                                          as_row(acmaHeightFromDiameterGamBal),
                                          #as_row(acmaHeightFromDiameterGamBalPhysio),
                                          as_row(acmaHeightFromDiameterGamPhysio),
                                          as_row(acmaHeightFromDiameterHossfeld),
                                          as_row(acmaHeightFromDiameterKorf),
                                          as_row(acmaHeightFromDiameterLinear),
                                          as_row(acmaHeightFromDiameterMichaelisMenten),
                                          as_row(acmaHeightFromDiameterParabolic),
                                          as_row(acmaHeightFromDiameterPower),
                                          as_row(acmaHeightFromDiameterProdan),
                                          as_row(acmaHeightFromDiameterRatkowsky),
                                          as_row(acmaHeightFromDiameterRichards),
                                          as_row(acmaHeightFromDiameterSharmaParton),
                                          as_row(acmaHeightFromDiameterSharmaPartonBal),
                                          as_row(acmaHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(acmaHeightFromDiameterSharmaPartonPhysio),
                                          as_row(acmaHeightFromDiameterSharmaZhang),
                                          as_row(acmaHeightFromDiameterSharmaZhangBal, significant = FALSE),
                                          as_row(acmaHeightFromDiameterSibbesen),
                                          as_row(acmaHeightFromDiameterWeibull),
                                          as_row(acmaHeightFromDiameterWeibullBal),
                                          as_row(acmaHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "ACMA3", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(acmaHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = acma2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = acma2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterKorf$fitted.values, color = "Korf", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterLinear$fitted.values, color = "linear", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterPower$fitted.values, color = "power", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterProdan$fitted.values, color = "Prodan", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = acma2016$isPlantation)) +
  annotate("text", x = 0, y = 40, label = "bigleaf maple, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## bigleaf maple height-diameter GNLS regressions
#acmaHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^b2, acma2016, start = acmaHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#acmaHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), acma2016, start = acmaHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.005, msTol = 1E-4, tolerance = 1E-3, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001, maxiter at default msTol and tolerance
##acmaHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = acmaHeightFromDiameterSharmaParton$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # foreign NaN-inf with plot correlation, step halving at nlsTOl = 1
#acmaHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = acmaHeightFromDiameterSharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 0.01, tolerance = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#acmaHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = acmaHeightFromDiameterSharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, msTol = 0.01, tolerance = 0.01, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.005
#acmaHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3*basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = acmaHeightFromDiameterSharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, msTol = 1E-4, tolerance = 1E-3, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.005
#acmaHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^b2)), acma2016, start = acmaHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, maxIter = 250, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE))
#acmaHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = acmaHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.01, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE))
#save(acmaHeightFromDiameterChapmanRichardsGnls, acmaHeightFromDiameterChapmanRichardsBalGnls, acmaHeightFromDiameterSharmaPartonBalGnls, acmaHeightFromDiameterSharmaPartonBalGnls, acmaHeightFromDiameterSharmaZhangGnls, acmaHeightFromDiameterSharmaZhangBalGnls, acmaHeightFromDiameterWeibullGnls, acmaHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia ACMA3 GNLS.rdata")

load("trees/height-diameter/HtDia ACMA3 GNLS.rdata")
acmaHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", acmaHeightFromDiameterChapmanRichardsGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BA+L GNLS", acmaHeightFromDiameterChapmanRichardsBalGnls, acma2016, acma2016natural, acma2016plantation)
#acmaHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", acmaHeightFromDiameterSharmaPartonGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BA+L GNLS", acmaHeightFromDiameterSharmaPartonBalGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", acmaHeightFromDiameterSharmaZhangGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BA+L GNLS", acmaHeightFromDiameterSharmaZhangBalGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", acmaHeightFromDiameterWeibullGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BA+L GNLS", acmaHeightFromDiameterWeibullBalGnls, acma2016, acma2016natural, acma2016plantation)

acmaHeightFromDiameterResultsGnls = bind_rows(as_row(acmaHeightFromDiameterChapmanRichardsGnls),
                                              as_row(acmaHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(name = "Sharma-Parton GNLS"),
                                              as_row(acmaHeightFromDiameterSharmaPartonBalGnls),
                                              as_row(acmaHeightFromDiameterSharmaZhangGnls),
                                              as_row(acmaHeightFromDiameterSharmaZhangBalGnls),
                                              as_row(acmaHeightFromDiameterWeibullGnls),
                                              as_row(acmaHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "ACMA3", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
acmaHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(acmaHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(acmaHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(acmaHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
#  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
#  select(-bal, -balN, -balP)

ggplot() +
  geom_point(aes(x = acma2016natural$DBH, y = acma2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = acma2016natural$DBH, y = acma2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "natural regeneration DBH, cm", y = "bigleaf maple naturally regenerated height, m") +
ggplot() +
  geom_point(aes(x = acma2016plantation$DBH, y = acma2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = acma2016plantation$DBH, y = acma2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "plantation DBH, cm", y = "bigleaf maple plantation height, m")

ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterWeibullBAL$fitted.values, color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = acma2016natural$DBH, y = acmaHeightFromDiameterWeibullBALnatural$fitted.values, color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = acma2016plantation$DBH, y = acmaHeightFromDiameterWeibullBALplantation$fitted.values, color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterBase$fitted.values, color = "base")) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterWeibull$fitted.values, color = "ElliottWeibull")) +
  annotate("text", x = 0, y = 85, label = "a) bigleaf maple, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## bigleaf maple diameter-height regressions
#acmaDiameterFromHeightChapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, iter = 10000,
#                                                  start_lower = list(a1 = -15, b1 = -10, b2 = -1), 
#                                                  start_upper = list(a1 = 150, b1 = 100, b2 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # NaN-inf
#acmaDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, acma2016, iter = 100,
#                                                   start_lower = list(a1 = -1, a2 = 0.011, b1 = -1, b2 = -1, b3 = -1), 
#                                                   start_upper = list(a1 = 100, a2 = 3, b1 = 10, b2 = 1, b3 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
#acmaDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 1, a2 = 1, a2p = 0, b1 = 0.02, b2 = 0.26, b2p = 0, b3 = 0.9, b3p = 0), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # signular gradient
acmaDiameterFromHeightChapmanForm = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 540000, b1 = 0.00001, b2 = 1.11), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls() at multiple nls_multstart() positions, NaN-inf from nlrob()
acmaDiameterFromHeightChapmanFormAat = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 46000, a2 = 159, b1 = 0.000013, b2 = 1.11), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls() and nlrob()
acmaDiameterFromHeightChapmanFormBal = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 23600, a2 = -1955, b1 = 0.00001, b2 = 1.07), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls(), step factor from nlrob()
acmaDiameterFromHeightChapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 261000, a2 = -7700, a3 = 4881, a9 = -216000, b1 = 0.00001, b2 = 1.06), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # step factor from nls()
acmaDiameterFromHeightChapmanFormRelHt = gsl_nls(DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 241000, a9 = -114000, b1 = 0.000007, b2 = 1.23), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # step factor from nls() and nlrob()
acmaDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), acma2016, start = list(a1 = 23.5, b1 = -0.022, b2 = 2.01, b2p = -0.17), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # b1p not significant
acmaDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), acma2016, start = list(a1 = 24.2, a2 = -0.005, b1 = -0.021, b2 = 2.02, b2p = -0.19), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = list(maxiter = 500))
acmaDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), acma2016physio, start = list(a1 = 3.30, a1p = 4.985, a4 = 0.00031, a5 = -0.172, a6 = -0.215, a7 = -0.839, a8 = 1.459, b1 = -0.0035, b1p = 0.00193, b2 = 1.798), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a4, a6, a7 not significant
acmaDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), acma2016, start = list(a1 = 36.1, a9 = -2.78, b1 = -0.028, b1p = 0.009, b2 = 1.63), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a1p, a2p, b2p not significant
acmaDiameterFromHeightGam = gam(DBH ~ s(TotalHt, by = as.factor(isPlantation)), data = acma2016)
acmaDiameterFromHeightGamAat = gam(DBH ~ s(TotalHt, tallerQuasiBasalArea, standQuasiBasalArea, by = as.factor(isPlantation)), data = acma2016)
#acmaDiameterFromHeightGamAatPhysio = gam(DBH ~ s(TotalHt, tallerQuasiBasalArea, standQuasiBasalArea, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = acma2016physio) # insufficient data
acmaDiameterFromHeightGamPhysio = gam(DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation)), data = acma2016physio)
acmaDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), acma2016, weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
acmaDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), acma2016, start = list(a1 = -116, a2 = -128, b1 = 1.50), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
acmaDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), acma2016, start = list(a1 = 5.1, a1p = -2.0, a2 = -0.12, a2p = -0.023), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
acmaDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), acma2016, weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # (TotalHt - 1.37)^2 not significant (p = 0.058)
acmaDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
acmaDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a2p not significant
acmaDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016physio, start = list(a1 = 4.38, a1p = -1.88, a4 = 0.00014, a5 = -1.61, a6 = -0.021, a7 = -0.112, a8 = 0.0096, b1 = 0.895, b1p = 0.189), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a4, a6, a7, a8 not significant
acmaDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, acma2016, start = list(a1 = 2.37, a1p = -0.887, a9 = -1.508, a9p = 1513, b1 = 1.123), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) 
acmaDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a1p not significant
acmaDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 161), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # singular gradient or step factor with nlrob()
acmaDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, acma2016, start = list(a1 = 0.0001, b1 = 3.8, b2 = 0.0174, b3 = 0.112, b4 = -2.4), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # NaN-inf from nls() at nls_multstart() positions, NaN-inf or code error with nlrob()
acmaDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # no significant plantation effects
acmaDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # a2 not significant
acmaDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), acma2016physio, start = list(a1 = 0.792, a1p = -0.272, a4 = 0.00003, a5 = -0.345, a6 = -0.0044, a7 = -0.020, a8 = 0.0021, b1 = 2.433, b2 = -0.164, b2p = 0.0277), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5)) # no physiographic parameters significant
acmaDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.496, a9 = -0.188, b1 = 2.31, b2 = -0.12), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5))
acmaDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, acma2016, start = list(a1 = -185000, b1 = 0.00001, b2 = 1.11), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -0.9), 0.5), control = gsl_nls_control(maxiter = 50)) # unfavorable to concave up curvature, NaN-inf with nlrob()
#confint_nlrob(acmaDiameterFromHeightSibbesenFormPhysio, level = 0.99, weights = pmin(acma2016physio$TotalHt^if_else(acma2016physio$isPlantation, -1.6, -0.9), 0.5))

acmaDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", acmaDiameterFromHeightChapmanForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AA+T", acmaDiameterFromHeightChapmanFormAat, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BA+L", acmaDiameterFromHeightChapmanFormBal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BA+L RelHt", acmaDiameterFromHeightChapmanFormBalRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", acmaDiameterFromHeightChapmanFormRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", acmaDiameterFromHeightChapmanRichards, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AA+T", acmaDiameterFromHeightChapmanRichardsAat, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", acmaDiameterFromHeightChapmanRichardsPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", acmaDiameterFromHeightChapmanRichardsRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightGam = get_dbh_error("GCV GAM", acmaDiameterFromHeightGam, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightGamAat = get_dbh_error("GCV GAM AA+T", acmaDiameterFromHeightGamAat, acma2016, acma2016natural, acma2016plantation)
#acmaDiameterFromHeightGamAatPhysio = get_dbh_error("GCV GAM AA+T physio", acmaDiameterFromHeightGamAatPhysio, acma2016physio, acma2016natural, acma2016plantation)
acmaDiameterFromHeightGamPhysio = get_dbh_error("GCV GAM physio", acmaDiameterFromHeightGamPhysio, acma2016physio, acma2016natural, acma2016plantation)
#save(acmaDiameterFromHeightGamAatPhysio, acmaDiameterFromHeightGamPhysio, file = "trees/height-diameter/HtDia ACMA3 spline DBH.rdata")
#load("trees/height-diameter/HtDia ACMA3 spline DBH.rdata")
acmaDiameterFromHeightLinear = get_dbh_error("linear", acmaDiameterFromHeightLinear, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", acmaDiameterFromHeightMichaelisMentenForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightNaslund = get_dbh_error("Näslund", acmaDiameterFromHeightNaslund, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightParabolic = get_dbh_error("parabolic", acmaDiameterFromHeightParabolic, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPower = get_dbh_error("power", acmaDiameterFromHeightPower, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerAat = get_dbh_error("power AA+T", acmaDiameterFromHeightPowerAat, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerPhysio = get_dbh_error("power physio", acmaDiameterFromHeightPowerPhysio, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", acmaDiameterFromHeightPowerRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightRuark = get_dbh_error("Ruark", acmaDiameterFromHeightRuark, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSchnute = get_dbh_error("Schnute", acmaDiameterFromHeightSchnute, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", acmaDiameterFromHeightSharmaParton, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", acmaDiameterFromHeightSibbesenForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AA+T", acmaDiameterFromHeightSibbesenFormAat, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", acmaDiameterFromHeightSibbesenFormPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", acmaDiameterFromHeightSibbesenFormRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightWeibull = get_dbh_error("Weibull", acmaDiameterFromHeightWeibull, acma2016, acma2016natural, acma2016plantation)

acmaDiameterFromHeightResults = bind_rows(as_row(acmaDiameterFromHeightChapmanRichards),
                                          as_row(acmaDiameterFromHeightChapmanRichardsAat),
                                          as_row(acmaDiameterFromHeightChapmanRichardsPhysio),
                                          as_row(acmaDiameterFromHeightChapmanRichardsRelHt),
                                          as_row(acmaDiameterFromHeightChapmanForm),
                                          as_row(acmaDiameterFromHeightChapmanFormAat),
                                          as_row(acmaDiameterFromHeightChapmanFormBal),
                                          as_row(acmaDiameterFromHeightChapmanFormBalRelHt),
                                          as_row(acmaDiameterFromHeightChapmanFormRelHt),
                                          as_row(acmaDiameterFromHeightLinear),
                                          as_row(acmaDiameterFromHeightGam),
                                          as_row(acmaDiameterFromHeightGamAat),
                                          #as_row(acmaDiameterFromHeightGamAatPhysio),
                                          as_row(acmaDiameterFromHeightGamPhysio),
                                          as_row(acmaDiameterFromHeightMichaelisMentenForm),
                                          as_row(acmaDiameterFromHeightNaslund),
                                          as_row(acmaDiameterFromHeightParabolic),
                                          as_row(acmaDiameterFromHeightPower),
                                          as_row(acmaDiameterFromHeightPowerAat),
                                          as_row(acmaDiameterFromHeightPowerPhysio),
                                          as_row(acmaDiameterFromHeightPowerRelHt),
                                          as_row(acmaDiameterFromHeightRuark),
                                          as_row(acmaDiameterFromHeightSchnute),
                                          as_row(acmaDiameterFromHeightSharmaParton),
                                          as_row(acmaDiameterFromHeightSibbesenForm),
                                          as_row(acmaDiameterFromHeightSibbesenFormAat),
                                          as_row(acmaDiameterFromHeightSibbesenFormPhysio, significant = FALSE),
                                          as_row(acmaDiameterFromHeightSibbesenFormRelHt),
                                          as_row(acmaDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "DBH", species = "ACMA3", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(acmaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(acma2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = acmaDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman-Richards approximate BA+L", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards BA+L", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards form"), na.rm = TRUE) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AA+T", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = -116*(TotalHt - 1.37)^1.5/(-128 - (TotalHt - 1.37)^1.5), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation), alpha = 0.5) +#geom_line(aes(x = acmaDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation)) +
  #geom_line(aes(x = -1/0.0003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute"), alpha = 0.5) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1.5*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^0.5*(TotalHt - 1.37)))^0.9, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = (-15 * log(1 - pmin(0.01*(TotalHt - 1.37), 0.999)))^-0.5, y = TotalHt, color = "Weibull"), na.rm = TRUE) +
  annotate("text", x = 0, y = 41, label = "bigleaf maple, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
acmaParameters = bind_rows(bind_rows(get_coefficients(acmaHeightFromDiameterChapmanRichards),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(acmaHeightFromDiameterCurtis),
                                     get_coefficients(acmaHeightFromDiameterGam),
                                     get_coefficients(acmaHeightFromDiameterGamBal),
                                     #get_coefficients(acmaHeightFromDiameterGamBalPhysio),
                                     get_coefficients(acmaHeightFromDiameterGamPhysio),
                                     get_coefficients(acmaHeightFromDiameterHossfeld),
                                     get_coefficients(acmaHeightFromDiameterKorf),
                                     get_coefficients(acmaHeightFromDiameterLinear),
                                     get_coefficients(acmaHeightFromDiameterMichaelisMenten),
                                     get_coefficients(acmaHeightFromDiameterParabolic),
                                     get_coefficients(acmaHeightFromDiameterPower),
                                     get_coefficients(acmaHeightFromDiameterProdan),
                                     get_coefficients(acmaHeightFromDiameterRatkowsky),
                                     get_coefficients(acmaHeightFromDiameterRichards),
                                     get_coefficients(acmaHeightFromDiameterSharmaParton),
                                     get_coefficients(acmaHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(acmaHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(acmaHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(acmaHeightFromDiameterSharmaZhang),
                                     get_coefficients(acmaHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(acmaHeightFromDiameterSibbesen),
                                     get_coefficients(acmaHeightFromDiameterWeibull),
                                     get_coefficients(acmaHeightFromDiameterWeibullBal),
                                     get_coefficients(acmaHeightFromDiameterWeibullBalRelHt),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(acmaHeightFromDiameterChapmanRichardsBalGnls),
                                     #get_coefficients(acmaHeightFromDiameterSharmaPartonGnls),
                                     get_coefficients(acmaHeightFromDiameterSharmaPartonBalGnls),
                                     get_coefficients(acmaHeightFromDiameterSharmaZhangGnls),
                                     get_coefficients(acmaHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(acmaHeightFromDiameterWeibullGnls),
                                     get_coefficients(acmaHeightFromDiameterWeibullBalGnls)) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(acmaDiameterFromHeightChapmanRichards),
                                     get_coefficients(acmaDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(acmaDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(acmaDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(acmaDiameterFromHeightChapmanForm),
                                     get_coefficients(acmaDiameterFromHeightChapmanFormAat),
                                     get_coefficients(acmaDiameterFromHeightChapmanFormBal),
                                     get_coefficients(acmaDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(acmaDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(acmaDiameterFromHeightGam),
                                     get_coefficients(acmaDiameterFromHeightGamAat),
                                     #get_coefficients(acmaDiameterFromHeightGamAatPhysio),
                                     get_coefficients(acmaDiameterFromHeightGamPhysio),
                                     get_coefficients(acmaDiameterFromHeightLinear),
                                     get_coefficients(acmaDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(acmaDiameterFromHeightNaslund),
                                     get_coefficients(acmaDiameterFromHeightParabolic),
                                     get_coefficients(acmaDiameterFromHeightPower),
                                     get_coefficients(acmaDiameterFromHeightPowerAat),
                                     get_coefficients(acmaDiameterFromHeightPowerPhysio),
                                     get_coefficients(acmaDiameterFromHeightPowerRelHt),
                                     get_coefficients(acmaDiameterFromHeightRuark),
                                     get_coefficients(acmaDiameterFromHeightSchnute),
                                     get_coefficients(acmaDiameterFromHeightSharmaParton),
                                     get_coefficients(acmaDiameterFromHeightSibbesenForm),
                                     get_coefficients(acmaDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(acmaDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(acmaDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(acmaDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "ACMA3",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
#acmaBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), acma2016, start = list(a1 = 1.9, b1 = 0.00007, b2 = 2.3), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # a1p, b1p not significant, nlrob() step factor
acmaBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), acma2016, start = list(a1 = 690, b1 = 0.00007, b2 = 2.2, b2p = -0.15), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # nlrob() step factor
acmaBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 2/7 * 0.25 * pi * 0.01^2, a1p = -0.0002, b1 = 2.05, b1p = 0.324), weights = pmin(1/basalArea, 1E4))
#confint2(acmaBasalAreaFromHeightKorf, level = 0.99)
#confint_nlrob(acmaBasalAreaFromHeightPower, level = 0.99)

acmaBasalAreaFromHeightKorf$fitted.values = predict(acmaBasalAreaFromHeightKorf, acma2016)
acmaBasalAreaFromHeightKorf$residuals = acmaBasalAreaFromHeightKorf$fitted.values - acma2016$basalArea
acmaBasalAreaFromHeightPower$fitted.values = predict(acmaBasalAreaFromHeightPower, acma2016)
acmaBasalAreaFromHeightPower$residuals = acmaBasalAreaFromHeightPower$fitted.values - acma2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(acmaBasalAreaFromHeightKorf), 100^2 * mean(acmaBasalAreaFromHeightKorf$residuals), mean(abs(acmaBasalAreaFromHeightKorf$residuals)), 1 - sum(acmaBasalAreaFromHeightKorf$residuals^2) / sum((acma2016$basalArea - mean(acma2016$basalArea)^2)),
        "power", AIC(acmaBasalAreaFromHeightPower), 100^2 * mean(acmaBasalAreaFromHeightPower$residuals), mean(abs(acmaBasalAreaFromHeightPower$residuals)), 1 - sum(acmaBasalAreaFromHeightPower$residuals^2) / sum((acma2016$basalArea - mean(acma2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(acma2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = acmaBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = acmaBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "bigleaf maple height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
