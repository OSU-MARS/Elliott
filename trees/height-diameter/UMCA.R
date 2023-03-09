# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Oregon myrtle (California bay) height-diameter regression form sweep
#umcaHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), umca2016, start = list(Ha = 14.4, Hap = -3.6, d = 2.78, dp = -0.99, kU = 0.054, kUp = 0.026))
umca2016 = trees2016 %>% filter(Species == "OM", isLiveUnbroken, TotalHt > 0) %>% # live Oregon myrtles measured for height
  mutate(dbhWeight = pmin(1/(1.27*DBH^0.72), 5),
         heightWeight = pmin(1/(3.48*(TotalHt - 1.37)^1.65), 5))
umca2016physio = umca2016 %>% filter(is.na(elevation) == FALSE)
umca2016gamConstraint = c(DBH = -0.8547/0.8790, TotalHt = 1.37, standBasalAreaPerHectare = median(umca2016$standBasalAreaPerHectare), basalAreaLarger = median(umca2016$basalAreaLarger), standBasalAreaApprox = median(umca2016$standBasalAreaApprox), tallerApproxBasalArea = median(umca2016$tallerApproxBasalArea), elevation = median(umca2016physio$elevation), slope = median(umca2016physio$slope), aspect = median(umca2016physio$aspect), topographicShelterIndex = median(umca2016physio$topographicShelterIndex), relativeHeight = median(umca2016$relativeHeight)) # point constraint for mgcv::s()

umca2016defaultWeight = umca2016 %>% mutate(dbhWeight = pmin(1/DBH, 5),
                                            heightWeight = pmin(1/TotalHt, 5))
umca2016defaultWeightPhysio = umca2016defaultWeight %>% filter(is.na(elevation) == FALSE)

umcaHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), umca2016))
umcaHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), umca2016)

umcaHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, umca2016, start = list(a1 = 18, b1 = -0.05, b2 = 1.06)) # a1p, b1p, b2p not significant
umcaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH))^b2, umca2016, start = list(a1 = 18, a2 = 0.23, b1 = -0.04, b2 = 1.0), significant = FALSE) # a2, a2p, a3, a3p, b1p, b2p not significant
umcaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016physio, start = list(a1 = 20, a2 = 0.03, a4 = -0.014, b1 = -0.06, b2 = 1.1), significant = FALSE) # a1p, a2, a2p, a3, a4, a5, a6, a7, a8, b1p, b2p not significant
umcaHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = -1.5, a3 = 0.11, a4 = 37, a4p = -14, b1 = -0.1, b2 = 0.389, b2p = 0.684)) # a1p, a2, a2p, a3p, b1p, b2p not significant, a1 debatable
umcaHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016physio, start = list(a1 = 20, a4 = -0.014, b1 = -0.05, b2 = 1.1), significant = FALSE) # a4, a5, a6, a7, a8, b1p, b2p not significant
umcaHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, umca2016, start = list(a1 = 1.6, b1 = 0.4)) # a1p, b1p not significant
umcaHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), umca2016, start = list(a1 = 21.4, a1p = -4.37, b1 = 43.8, b1p = -19.9, b2 = -1.27)) # b2p not significant
umcaHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*DBH^b2), umca2016, start = list(a1 = 49.8, b1 = -5.02, b1p = 0.404, b2 = -0.386)) # a1p, b2p not significant
umcaHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27)) # b1p not significant
umcaHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016, start = list(a1 = 0.037, a1p = 0.020, a2 = 1.109, a2p = -0.472, a3 = 1.242)) # a3p not significant
umcaHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, umca2016, start = list(a1 = 0.821, a1p = 0.206, b1 = 0.810)) # b1p not significant
umcaHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(DBH + b2)), umca2016, start = list(a1 = 21.6, b1 = -16.5, b1p = 1.974, b2 = 3.629)) # a1p, b2p not significant
umcaHeightFromDiameter$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap * isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016, start = list(Ha = 14.6, Hap = -4.146, d = 2.198, kU = 0.048, kUp = 0.045)) # dp not significant, job step factor with nlrob()
umcaHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, umca2016, start = list(a1 = 8, b1 = 0.2, b2 = -0.045, b2p = -0.01, b3 = 0.085, b4 = 1.2)) # a1p, b1p, b3p, b4p not significant, job step factor with nlrob()
umcaHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 14, b1 = 0.05, b2 = -0.06, b3 = -0.01, b4 = 1.2, b4p = -0.22)) # a1p, b1, b1p, b2p, b3p not significant
umcaHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 15, a4 = -0.01, b1 = 0.1, b2 = -0.06, b3 = -0.06, b4 = 1.2, b4p = -0.22)) # a5, a6, a7, a8, b1p, b2p, b3p not significant, job step factor with nlrob()
umcaHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 13, a4 = -0.01, b1 = 0.08, b2 = -0.05, b3 = 0.0, b4 = 1.2, b4p = -0.23)) # a5, a6, a7, a8, b1p, b2p, b3, b3p not significant, job step factor with nlrob
umcaHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 13, b1 = 0.1, b2 = -0.06, b3 = 0.125, b4 = 1.1, b4p = -0.2)) # a1p, b1, b1p, b2p, b3p not significant, job step factor with nlrob()
umcaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, umca2016, start = list(a1 = 10, a2 = 0.03, b1 = 0.04, b2 = -0.04, b3 = 0.07, b3p = 0.05, b4 = 1.1), significant = FALSE) # a1p, a2, a2p, b1, b1p, b2p, b4p not significant, job step factor with nlrob()
umcaHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.3, a1p = 0.1, b1 = 2.0, b2 = -0.16, b2p = -0.03)) # b1p not significant
umcaHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), umca2016, start = list(a1 = 18, b1 = -0.04, b2 = 1.0)) # a1p, b1p, b2p not significant
umcaHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH^b2)), umca2016, start = list(a1 = 18, a2 = 0.02, b1 = -0.044, b2 = 1.1), significant = FALSE) # a1p, a2, a2p, a3, a3p, b1p, b2p not significant
umcaHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a9*pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), umca2016, start = list(a1 = -2.5, a2 = 0.09, a9 = 55, b1 = -0.5, b1p = 0.3, b2 = 0.4), control = nls.control(maxiter = 500)) # a1p, a2p, a3, a3p, a4p, b2p not significant

umcaHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, umca2016, start = list(a1 = 18, b1 = -0.05, b2 = 1.06)))
#umcaHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH))^b2, umca2016, start = list(a1 = 18, a2 = 0.23, b1 = -0.04, b2 = 1.0), significant = FALSE)
umcaHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016physio, start = list(a1 = 20, a2 = 0.03, a4 = -0.014, b1 = -0.06, b2 = 1.1), significant = FALSE)
#umcaHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = -1.5, a3 = 0.11, a4 = 37, a4p = -14, b1 = -0.1, b2 = 0.389, b2p = 0.684))
umcaHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016physio, start = list(a1 = 20, a4 = -0.014, b1 = -0.05, b2 = 1.1), significant = FALSE)
umcaHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, umca2016, start = list(a1 = 1.6, b1 = 0.4)) # a1p, b1p not significant
umcaHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), umca2016, start = list(a1 = 21.4, a1p = -4.37, b1 = 43.8, b1p = -19.9, b2 = -1.27))
umcaHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*DBH^b2), umca2016, start = list(a1 = 49.8, b1 = -5.02, b1p = 0.404, b2 = -0.386))
umcaHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27))
umcaHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016, start = list(a1 = 0.037, a1p = 0.020, a2 = 1.109, a2p = -0.472, a3 = 1.242))
umcaHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, umca2016, start = list(a1 = 0.821, a1p = 0.206, b1 = 0.810))
umcaHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(DBH + b2)), umca2016, start = list(a1 = 21.6, b1 = -16.5, b1p = 1.974, b2 = 3.629))
#umcaHeightFromDiameterNlrob$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap * isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016, start = list(Ha = 14.6, Hap = -4.146, d = 2.198, kU = 0.048, kUp = 0.045))
#umcaHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, umca2016, start = list(a1 = 8, b1 = 0.2, b2 = -0.045, b2p = -0.01, b3 = 0.085, b4 = 1.2))
umcaHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 14, b1 = 0.05, b2 = -0.06, b3 = -0.01, b4 = 1.2, b4p = -0.22))
#umcaHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 15, a4 = -0.01, b1 = 0.1, b2 = -0.06, b3 = -0.06, b4 = 1.2, b4p = -0.22))
#umcaHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 13, a4 = -0.01, b1 = 0.08, b2 = -0.05, b3 = 0.0, b4 = 1.2, b4p = -0.23))
#umcaHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 13, b1 = 0.1, b2 = -0.06, b3 = 0.125, b4 = 1.1, b4p = -0.2))
#umcaHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, umca2016, start = list(a1 = 10, a2 = 0.03, b1 = 0.04, b2 = -0.04, b3 = 0.07, b3p = 0.05, b4 = 1.1), significant = FALSE)
umcaHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.3, a1p = 0.1, b1 = 2.0, b2 = -0.16, b2p = -0.03))
umcaHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), umca2016, start = list(a1 = 18, b1 = -0.04, b2 = 1.0))
umcaHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH^b2)), umca2016, start = list(a1 = 18, a2 = 0.02, b1 = -0.044, b2 = 1.1), significant = FALSE)
umcaHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a9*pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), umca2016, start = list(a1 = -2.5, a2 = 0.09, a9 = 55, b1 = -0.5, b1p = 0.3, b2 = 0.4), control = nls.control(maxiter = 500))
#confint_nlrob(umcaHeightFromDiameterNlrob$sharmaZhangBal, level = 0.99, weights = pmin(umca2016$DBH^-0.8, 1))

umcaHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, umca2016defaultWeight, start = list(a1 = 18, b1 = -0.05, b2 = 1.06)))
umcaHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH))^b2, umca2016defaultWeight, start = list(a1 = 18, a2 = 0.23, b1 = -0.04, b2 = 1.0), significant = FALSE)
umcaHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016defaultWeightPhysio, start = list(a1 = 20, a2 = 0.03, a4 = -0.014, b1 = -0.06, b2 = 1.1), significant = FALSE)
umcaHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016defaultWeight, start = list(a1 = -1.5, a3 = 0.11, a4 = 37, a4p = -14, b1 = -0.1, b2 = 0.389, b2p = 0.684))
umcaHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation) * (1 - exp(b1*DBH))^b2, umca2016defaultWeightPhysio, start = list(a1 = 20, a4 = -0.014, b1 = -0.05, b2 = 1.1), significant = FALSE)
umcaHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, umca2016defaultWeight, start = list(a1 = 1.6, b1 = 0.4))
umcaHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), umca2016defaultWeight, start = list(a1 = 21.4, a1p = -4.37, b1 = 43.8, b1p = -19.9, b2 = -1.27))
umcaHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*DBH^b2), umca2016defaultWeight, start = list(a1 = 49.8, b1 = -5.02, b1p = 0.404, b2 = -0.386))
umcaHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016defaultWeight, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27))
umcaHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016defaultWeight, start = list(a1 = 0.037, a1p = 0.020, a2 = 1.109, a2p = -0.472, a3 = 1.242))
umcaHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, umca2016defaultWeight, start = list(a1 = 0.821, a1p = 0.206, b1 = 0.810))
umcaHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(DBH + b2)), umca2016defaultWeight, start = list(a1 = 21.6, b1 = -16.5, b1p = 1.974, b2 = 3.629))
umcaHeightFromDiameterGslNlsDefault$richards = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap * isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016defaultWeight, start = list(Ha = 14.6, Hap = -4.146, d = 2.198, kU = 0.048, kUp = 0.045))
umcaHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, umca2016defaultWeight, start = list(a1 = 8, b1 = 0.2, b2 = -0.045, b2p = -0.01, b3 = 0.085, b4 = 1.2))
umcaHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016defaultWeight, start = list(a1 = 14, b1 = 0.05, b2 = -0.06, b3 = -0.01, b4 = 1.2, b4p = -0.22))
umcaHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016defaultWeightPhysio, start = list(a1 = 15, a4 = -0.01, b1 = 0.1, b2 = -0.06, b3 = -0.06, b4 = 1.2, b4p = -0.22))
umcaHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^(b4 + b4p * isPlantation), umca2016defaultWeightPhysio, start = list(a1 = 13, a4 = -0.01, b1 = 0.08, b2 = -0.05, b3 = 0.0, b4 = 1.2, b4p = -0.23))
umcaHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), umca2016defaultWeight, start = list(a1 = 13, b1 = 0.1, b2 = -0.06, b3 = 0.125, b4 = 1.1, b4p = -0.2))
umcaHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, umca2016defaultWeight, start = list(a1 = 10, a2 = 0.03, b1 = 0.04, b2 = -0.04, b3 = 0.07, b3p = 0.05, b4 = 1.1), significant = FALSE)
umcaHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016defaultWeight, start = list(a1 = 0.3, a1p = 0.1, b1 = 2.0, b2 = -0.16, b2p = -0.03))
umcaHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), umca2016defaultWeight, start = list(a1 = 18, b1 = -0.04, b2 = 1.0))
umcaHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH^b2)), umca2016defaultWeight, start = list(a1 = 18, a2 = 0.02, b1 = -0.044, b2 = 1.1), significant = FALSE)
umcaHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a9*pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), umca2016defaultWeight, start = list(a1 = -2.5, a2 = 0.09, a9 = 55, b1 = -0.5, b1p = 0.3, b2 = 0.4), control = nls.control(maxiter = 500))

umcaHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 11, pc = umca2016gamConstraint), data = umca2016)
umcaHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15, pc = umca2016gamConstraint), data = umca2016)
umcaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, basalAreaLarger, elevation, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 36, pc = umca2016gamConstraint), data = umca2016physio)
umcaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 35, pc = umca2016gamConstraint), data = umca2016physio)

if (htDiaOptions$includeInvestigatory)
{
  print(umcaHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

  ggplot() +
    geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = umca2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = umca2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$curtis), color = "Curtis", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$korf), color = "Korf", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$linear), color = "linear", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$parabolic), color = "parabolic", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$power), color = "power", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$prodan), color = "Prodan", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$richards), color = "unified Richards", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$sibbesen), color = "Sibbesen", group = umca2016$isPlantation)) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$weibull), color = "Weibull", group = umca2016$isPlantation)) +
    annotate("text", x = 0, y = 35, label = "Oregon myrtle, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 35)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}

  
## Oregon myrtle height-diameter GNLS regressions
if (htDiaOptions$fitGnls)
{
  umcaHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = umcaHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1)) # step halving with plot correlation
  umcaHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, umca2016, start = umcaHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # corSymm viable but dropped
  umcaHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, umca2016, start = umcaHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # step halving with plot correlation
  umcaHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = umcaHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # step halving with plot correlation
  umcaHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = umcaHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # step halving with plot correlation
  umcaHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, umca2016, start = umcaHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.005, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # step halving with plot correlation or nlsTol = 0.001
  umcaHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = umcaHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # step halving with plot correlation
  umcaHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = umcaHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # corSymm viable but dropped

  save(file = "trees/height-diameter/data/UMCA GNLS.rdata", umcaHeightFromDiameterGnls)
}
if (htDiaOptions$includeInvestigatory)
{
  umcaHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(umcaHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(umcaHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(umcaHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  
  ggplot() +
    geom_point(aes(x = umca2016natural$DBH, y = umca2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = umca2016natural$DBH, y = umca2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "Oregon myrtle naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = umca2016plantation$DBH, y = umca2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = umca2016plantation$DBH, y = umca2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "Oregon myrtle plantation height, m")
  
  ggplot() +
    geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$weibullBAL), color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = umca2016natural$DBH, y = predict(umcaHeightFromDiameter$weibullBALnatural), color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = umca2016plantation$DBH, y = predict(umcaHeightFromDiameter$weibullBALplantation), color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$Base), color = "base")) +
    geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameter$weibull), color = "ElliottWeibull")) +
    annotate("text", x = 0, y = 85, label = "a) Oregon myrtle, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}

  
## Oregon myrtle diameter-height regressions
umcaDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), umca2016))
umcaDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), umca2016, significant = FALSE) # collapses to linear since (TotalHt - 1.37)^2 and isPlantation*(TotalHt - 1.37)^2 not significant

umcaDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 20, b1 = 0.1, b2 = 0.6)) # subject to a1-b1 evaporation, NaN-inf with nls() at multiple nls_multstart() points, NaN-inf with nlrob(), job NaN-inf with gsl_nls()
umcaDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 20, a2 = 0, b1 = 0.1, b2 = 0.6), significant = FALSE) # a2 not significant, a1-b1 parameter evaporation, NaN-inf with nls(), step factor with nlrob(), step factor with default gsl_nls tolerance
umcaDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 1, a2 = -0.005, b1 = 1.9, b2 = 0.3)) # NaN-inf with nls(), step factor with nlrob()
umcaDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 0.4, a2 = -0.012, a3 = 0.01, a9 = -0.3, b1 = 2.7, b2 = 0.23)) # step factor with nls() and nlrob()
umcaDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 0.9, a9 = -0.6, b1 = 1.8, b2 = 0.3)) # step factor with nls(), step factor with nlrob()
umcaDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 34, b1 = -0.047, b2 = 1.43, b2p = -0.15), control = gsl_nls_control(maxiter = 250, xtol = 1E-5)) # a1-b1 parameter evaporation, b2p debatable
umcaDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 35, a2 = -3, b1 = -0.01, b2 = 0.9, b2p = -0.2), control = gsl_nls_control(maxiter = 250, xtol = 1E-5), significant = FALSE) # a2 not significant, a1-b1 parameter evaporation
umcaDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), umca2016physio, start = list(a1 = 58.6, a5 = -52.7, b1 = -0.0637, b2 = 0.9), control = gsl_nls_control(maxiter = 250, xtol = 1E-5), significant = FALSE) # a1p, a4, a5, a6, a7, a8, b1p, b2p not significant, a1+a5-b1 parameter evaporation
umcaDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), umca2016, start = list(a1 = 35, a9 = -15, b1 = -0.03, b2 = 0.9), control = gsl_nls_control(maxiter = 250, xtol = 1E-5), significant = FALSE) # a1p, a9, b1p, b2p not significant, a1-b1 parameter evaporation
umcaDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), umca2016, start = list(a1 = 100, a2 = 100, b1 = 1)) # collapses to linear
umcaDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + a2*sqrt(TotalHt - 1.37)), umca2016, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.14)) # a2p not significant, a1p debatable
umcaDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332))
#umcaDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696)) # a1p, a2p not significant
#umcaDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016physio, start = list(a1 = 4, a1p = -0.95, a4 = 0.002, a5 = -2.8, a6 = 0.4, a8 = 0.05, b1 = 0.94)) # a1p, a2, a7, b1p not significant, a4 debatable
#umcaDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 2.37, a1p = -0.887, a9 = -1.508, a9p = 1513, b1 = 1.123))
umcaDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016, start = list(a1 = 3.95, b1 = 0.45, b2 = 0.06)) # a1p, b1p, b2p not significant
umcaDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016, start = list(a1 = 4.0, a2 = 0, b1 = 0.55, b2 = 0.045), significant = FALSE) # a1p, a2, a2p, a3, a3p, b1p, b2p not significant
umcaDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016physio, start = list(a1 = 6, a5 = 0, b1 = 0.6, b2 = 0.04), significant = FALSE) # a1p, a4, a5, a6, a7, a8, b1p, b2p not significant
umcaDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 4, a9 = -2, b1 = 0.7, b1p = -0.5, b2 = 0.04, b2p = 0.08), significant = FALSE) # a9, a9p not significant
umcaDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016physio, start = list(a1 = 6, a5 = -2, a9 = -2, b1 = 0.6, b2 = 0.05)) # a5 significant
umcaDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), umca2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 177)) # NaN-inf with nlrob()
#umcaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, umca2016, start = list(a1 = 27, b1 = 0.62, b2 = 0.03, b3 = -3.1, b4 = 0.17), control = gsl_nls_control(maxiter = 250, xtol = 1E-3)) # a2p, b2p, b3p not significant, NaN-inf with nls() at nls_multistart() point, NaN-inf or code error with nlrob(), NaN-inf with gls_nls()
umcaDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15)) # no significant plantation effects
umcaDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156)) # a2 not significant
umcaDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, start = list(a1 = 3.3, a4 = 0.002, b1 = 0.48, b2 = 0.21), significant = FALSE) # a1p, a4, a5, a6, a7, a8 not significant, step factor with b1p
umcaDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.496, a9 = -0.188, b1 = 2.31, b2 = -0.12), significant = FALSE) # a9 not significant
umcaDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, start = list(a1 = 4.2, a4 = 0, a9 = -3, b1 = 0.53, b2 = 0.20)) # a4 not significant
umcaDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, umca2016, start = list(a1 = -200, b1 = 0.027, b2 = 0.8), control = gsl_nls_control(maxiter = 500)) # NaN-inf with nlrob()
#confint_nlrob(umcaDiameterFromHeight$weibull, level = 0.99)

umcaDiameterFromHeightNlrob = list(naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + a2*sqrt(TotalHt - 1.37)), umca2016, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.14)))
umcaDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332))
#umcaDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696))
#umcaDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016physio, start = list(a1 = 4, a1p = -0.95, a4 = 0.002, a5 = -2.8, a6 = 0.4, a8 = 0.05, b1 = 0.94))
#umcaDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 2.37, a1p = -0.887, a9 = -1.508, a9p = 1513, b1 = 1.123))
umcaDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016, start = list(a1 = 3.95, b1 = 0.45, b2 = 0.06))
umcaDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016, start = list(a1 = 4.0, a2 = 0, b1 = 0.55, b2 = 0.045), significant = FALSE)
umcaDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016physio, start = list(a1 = 6, a5 = 0, b1 = 0.6, b2 = 0.04), significant = FALSE)
umcaDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 4, a9 = -2, b1 = 0.7, b1p = -0.5, b2 = 0.04, b2p = 0.08), significant = FALSE)
umcaDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016physio, start = list(a1 = 4.5, a5 = 0, a9 = -3, b1 = 0.45, b2 = 0.08), significant = FALSE)
#umcaDiameterFromHeightNlrob$schnute = fit_nlrob("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), umca2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 177)) # NaN-inf
umcaDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15))
umcaDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156))
umcaDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, start = list(a1 = 3.3, a4 = 0.002, b1 = 0.48, b2 = 0.21), significant = FALSE)
umcaDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.496, a9 = -0.188, b1 = 2.31, b2 = -0.12), significant = FALSE)
#umcaDiameterFromHeightNlrob$sibbesenReplaceRelHtPhysio = fit_nlrob("Sibbesen replace RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, maxit = 500, start = list(a1 = 4.2, a4 = 0, a9 = -2.8, b1 = 0.53, b2 = 0.18)) # >500 steps to converge

umcaDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016defaultWeight, start = list(a1 = 20, b1 = 0.1, b2 = 0.6), control = gsl_nls_control(maxiter = 250, xtol = 0.002))) # a1-b1 evaporation
umcaDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016defaultWeight, start = list(a1 = 20, a2 = 0, b1 = 0.1, b2 = 0.6), control = gsl_nls_control(maxiter = 250, xtol = 0.002), significant = FALSE) # a1-b1 evaporation
umcaDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016defaultWeight, start = list(a1 = 0.9, a9 = -0.6, b1 = 1.8, b2 = 0.3), control = gsl_nls_control(maxiter = 250, xtol = 0.002)) # a1 and b1 both tend to zero
umcaDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016defaultWeight, start = list(a1 = 34, b1 = -0.047, b2 = 1.43, b2p = -0.15), control = gsl_nls_control(maxiter = 500, xtol = 0.002))
umcaDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016defaultWeight, start = list(a1 = 35, a2 = -3, b1 = -0.01, b2 = 0.9, b2p = -0.2), control = gsl_nls_control(maxiter = 500, xtol = 0.001), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), umca2016defaultWeightPhysio, start = list(a1 = 58.6, a5 = -52.7, b1 = -0.0637, b2 = 0.9), control = gsl_nls_control(maxiter = 500, xtol = 0.001), significant = FALSE) # a1-b1 evaporation
umcaDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), umca2016defaultWeight, start = list(a1 = 35, a9 = -15, b1 = -0.03, b2 = 0.9), control = gsl_nls_control(maxiter = 500, xtol = 0.001), significant = FALSE) # a1-b1 evaporation
umcaDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), umca2016defaultWeight, start = list(a1 = 100, a2 = 100, b1 = 1))
umcaDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + a2*sqrt(TotalHt - 1.37)), umca2016defaultWeight, start = list(a1 = 3.9, a1p = -1.0, a2 = -0.14))
umcaDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016defaultWeight, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332))
#umcaDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016defaultWeight, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696))
#umcaDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016defaultWeightPhysio, start = list(a1 = 4, a1p = -0.95, a4 = 0.002, a5 = -2.8, a6 = 0.4, a8 = 0.05, b1 = 0.94))
#umcaDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016defaultWeight, start = list(a1 = 2.37, a1p = -0.887, a9 = -1.508, a9p = 1513, b1 = 1.123))
umcaDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016defaultWeight, start = list(a1 = 3.95, b1 = 0.45, b2 = 0.06))
umcaDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016defaultWeight, start = list(a1 = 4.0, a2 = 0, b1 = 0.55, b2 = 0.045), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016defaultWeightPhysio, start = list(a1 = 6, a5 = 0, b1 = 0.6, b2 = 0.04), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016defaultWeight, start = list(a1 = 4, a9 = -2, b1 = 0.7, b1p = -0.5, b2 = 0.04, b2p = 0.08), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope) + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), umca2016defaultWeightPhysio, start = list(a1 = 3.5, a5 = -1.8, a9 = -1.4, b1 = 1.2, b2 = 0.0), significant = FALSE) # a5, a9, b2 not significant
umcaDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), umca2016defaultWeight, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 177))
#umcaDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, umca2016defaultWeight, start = list(a1 = 27, b1 = 0.62, b2 = 0.03, b3 = -3.1, b4 = 0.17), control = gsl_nls_control(maxiter = 250, xtol = 1E-3))
umcaDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016defaultWeight, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15))
umcaDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016defaultWeight, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156))
umcaDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016defaultWeightPhysio, start = list(a1 = 3.3, a4 = 0.002, b1 = 0.48, b2 = 0.21), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016defaultWeight, start = list(a1 = 0.496, a9 = -0.188, b1 = 2.31, b2 = -0.12), significant = FALSE)
umcaDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016defaultWeightPhysio, start = list(a1 = 4.2, a4 = 0, a9 = -2.8, b1 = 0.53, b2 = 0.18)) # a4 not significant
umcaDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, umca2016defaultWeight, start = list(a1 = -200, b1 = 0.027, b2 = 0.8), control = gsl_nls_control(maxiter = 500))

umcaDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = umca2016gamConstraint), data = umca2016)
umcaDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 14, pc = umca2016gamConstraint), data = umca2016)
umcaDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, standBasalAreaApprox, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 28, pc = umca2016gamConstraint), data = umca2016physio) # >5 minute Zen 3 3.4 GHz fit time with ABA+T, elevation, slope, sin(aspect), and TSI: elevation dropped based on physio AICs, then TSI+aspect
umcaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 60, pc = umca2016gamConstraint), data = umca2016physio) # 72.5% deviance explained (AIC 6089): 75.7% (5995) without elevation, 71.3% (6111) without slope, 72.4% (6043) without aspect, 72.7% (6034) without topographic shelter
umcaDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 10, pc = umca2016gamConstraint), data = umca2016, nthreads = 4)
umcaDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 57, pc = umca2016gamConstraint), data = umca2016physio, nthreads = 6)

if (htDiaOptions$includeInvestigatory)
{
  print(umcaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

  ggplot(umca2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanReplaceAbat), y = TotalHt, color = "Chapman-Richards replace approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$schnute), y = TotalHt, color = "Schnute inverse", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen replace", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 1*(TotalHt - 1.37)^1.5, y = TotalHt, color = "power"), alpha = 0.5) +
    geom_line(aes(x = 1*(TotalHt - 1.37)^1.5*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #annotate("text", x = 0, y = 37, label = "Oregon myrtle, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## collect model parameters
#if (exists("umcaHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/UMCA GNLS.rdata") }
umcaCoefficients = bind_rows(bind_rows(bind_rows(lapply(umcaHeightFromDiameter, get_list_coefficients)),
                                       bind_rows(lapply(umcaHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                       bind_rows(lapply(umcaHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                       #bind_rows(lapply(umcaHeightFromDiameterGnls, get_model_coefficients))) %>%
                               mutate(responseVariable = "height"),
                             bind_rows(bind_rows(lapply(umcaDiameterFromHeight, get_list_coefficients)),
                                       bind_rows(lapply(umcaDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                       bind_rows(lapply(umcaDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                               mutate(responseVariable = "DBH")) %>%
  mutate(species = "UMCA")
umcaResults = bind_rows(bind_rows(bind_rows(lapply(umcaHeightFromDiameter, get_list_stats)),
                                  bind_rows(lapply(umcaHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                  bind_rows(lapply(umcaHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                                  #bind_rows(lapply(umcaHeightFromDiameterGnls, get_model_stats))) %>%
                          mutate(responseVariable = "height"),
                        bind_rows(bind_rows(lapply(umcaDiameterFromHeight, get_list_stats)),
                                  get_model_stats(name = "modified Sharma-Parton", fitting = "gsl_nls"),
                                  bind_rows(lapply(umcaDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                  bind_rows(lapply(umcaDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                          mutate(responseVariable = "DBH")) %>%
  mutate(species = "UMCA")

save(file = "trees/height-diameter/data/UMCA results.Rdata", umcaCoefficients, umcaResults)
if (htDiaOptions$folds * htDiaOptions$repetitions <= htDiaOptions$retainModelThreshold)
{
  save(file = "trees/height-diameter/data/UMCA models and stats.Rdata", 
       umcaHeightFromDiameter, umcaHeightFromDiameterGslNlsDefault, umcaHeightFromDiameterNlrob,
       umcaDiameterFromHeight, umcaDiameterFromHeightGslNlsDefault, umcaDiameterFromHeightNlrob)
} else {
  save(file = "trees/height-diameter/data/UMCA stats.Rdata", 
       umcaHeightFromDiameter, umcaHeightFromDiameterGslNlsDefault, umcaHeightFromDiameterNlrob,
       umcaDiameterFromHeight, umcaDiameterFromHeightGslNlsDefault, umcaDiameterFromHeightNlrob)
}


## preferred forms identified (results.R, Figure 7)
# Oregon myrtle     Prodan            Chapman-Richards physio        REML GAM                Chapman-Richards physio
#                   Sibbesen          Chapman-Richards BA+L physio   power                   Sibbesen form Physio
#                   REML GAM                                         parabolic
umcaHeightFromDiameterPreferred = list(chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation + a7 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 19, a2 = 0.05, a4 = -0.015, a7 = -1.2, b1 = -0.05, b2 = 1.15, b2p = -0.19), folds = 1, repetitions = 1))
umcaHeightFromDiameterPreferred$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a6 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 20, a4 = -0.014, a6 = -1.2, b1 = -0.05, b2 = 1.15, b2p = -0.16), folds = 1, repetitions = 1)
umcaHeightFromDiameterPreferred$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 11, pc = umca2016gamConstraint), data = umca2016, folds = 1, repetitions = 1)
umcaHeightFromDiameterPreferred$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016, start = list(a1 = 0.040, a1p = 0.017, a2 = 1.01, a2p = -0.40, a3 = 1.46), folds = 1, repetitions = 1)
umcaHeightFromDiameterPreferred$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.30, a1p = 0.13, b1 = 1.96, b2 = -0.165, b2p = -0.034), folds = 1, repetitions = 1)

umcaDiameterFromHeightPreferred = list(chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), umca2016physio, start = list(a1 = 60, a5 = -52.7, b1 = -0.01, b2 = 0.9), control = gsl_nls_control(maxiter = 250, xtol = 1E-5), folds = 1, repetitions = 1))
umcaDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = umca2016gamConstraint), data = umca2016, folds = 1, repetitions = 1)
umcaDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), umca2016, significant = FALSE, folds = 1, repetitions = 1)
umcaDiameterFromHeightPreferred$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332), folds = 1, repetitions = 1)
umcaDiameterFromHeightPreferred$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, start = list(a1 = 1.53, a1p = -0.38, a4 = 0.001, a5 = -1.15, a6 = 0.14, a8 = 0.02, b1 = 1.95, b2 = -0.15), folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/UMCA preferred models.Rdata", umcaHeightFromDiameterPreferred, umcaDiameterFromHeightPreferred)


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  #umcaBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), umca2016, start = list(a1 = 0.3, b1 = 0.0006, b2 = 2.1), weights = pmin(1/basalArea, 1E4)) # step factor with nlrob()
  umcaBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), umca2016, start = list(a1 = 1.36, b1 = 0.0002, b2 = 2.06, b2p = -0.27), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor with nlrob()
  umcaBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^b1, umca2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.0002, b1 = 2.00), weights = pmin(1/basalArea, 1E4)) # b1p not significant
  #confint2(umcaBasalAreaFromHeightKorf, level = 0.99)
  #confint_nlrob(umcaBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/umca2016$basalArea, 1E4))
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(umcaBasalAreaFromHeightKorf), 100^2 * mean(residuals(umcaBasalAreaFromHeightKorf)), mean(abs(residuals(umcaBasalAreaFromHeightKorf))), 1 - sum(residuals(umcaBasalAreaFromHeightKorf)^2) / sum((umca2016$basalArea - mean(umca2016$basalArea)^2)),
          "power", AIC(umcaBasalAreaFromHeightPower), 100^2 * mean(residuals(umcaBasalAreaFromHeightPower)), mean(abs(residuals(umcaBasalAreaFromHeightPower))), 1 - sum(residuals(umcaBasalAreaFromHeightPower)^2) / sum((umca2016$basalArea - mean(umca2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(umca2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(umcaBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(umcaBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "Oregon myrtle height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}