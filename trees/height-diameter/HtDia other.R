# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## minority species height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Parton, Sharma-Zhang, Chapman-Richards BAL
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 43.0, a1p = 13.5, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), control = list(maxiter = 50)) # a1p not significant
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 43.0, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875, b2p = 0), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), control = list(maxiter = 500)) # > 500 iterations
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 89.3, a2 = 1.166, a2p = 0, a3 = -0.039, a3p = 0, b1 = -0.00544, b2 = 0.873), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a2p, a3p not significant
#otherHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 58.1, a1p = -47.4, a2 = 0.162, a2p = 0.032, b1 = -0.021, b1p = -0.222, b2 = -0.292, b2p = 0.036, b3 = 0.818, b3p = 0.165), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1))
#otherHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.714, a1p = 0.088, b1 = 1.172, b2 = -0.074, b2p = -0.0040), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b2p not significant
other2016 = trees2016 %>% filter((Species %in% c("DF", "RA", "WH", "BM", "OM", "RC")) == FALSE, isLiveUnbroken, TotalHt > 0) # live western redcedars measured for height
other2016natural = other2016 %>% filter(isPlantation == FALSE)
other2016physio = other2016 %>% filter(is.na(elevation) == FALSE)
other2016plantation = other2016 %>% filter(isPlantation)
other2016plantationPhysio = other2016physio %>% filter(isPlantation)
#otherConifer2016 = other2016 %>% filter(Species %in% c("XX", "CX", "SS", "PC", "PY", "GF", "LP"))
#otherHardwood2016 = other2016 %>% filter(Species %in% c("CA", "HX", "CH", "PM", "GC", "PD", "TO", "WI", "OA", "WO"))

otherHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016, start = list(a1 = 24.7, b1 = -0.033, b1p = -0.003, b2 = 1.000), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 81.9, a2 = 1.97, a3 = -1.23, b1 = -0.006, b2 = 0.881), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), control = list(maxiter = 50)) # a1p, a2p, a2, a3p, a3, b1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, other2016physio, start = list(a1 = 28.5, a2 = 0.286, a3 = -0.011, a4 = -0.096, a5 = -0.602, a6 = 0.661, a7 = -0.176, b1 = -0.022, b2 = 0.908), maxit = 50, weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2, a2p, a3, a4, a5, a6, b1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = -1.16, a2 = -0.05, a3 = 0.07, a3p = 0.09, a4 = 16.6, a4p = -6.7, b1 = -0.05, b2 = 0.28, b2p = 0.131), maxit = 50, weights = pmin(DBH^if_else(isPlantation, -1.1, -1.94), 1)) # a2, a2p, a3, b1p not significant, fails to converge with natural regen weight power of -1.95 or higher
otherHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016physio, start = list(a1 = 27.0, a2 = -0.015, a3 = -0.025, a4 = 0.276, a5 = -0.514, a6 = -0.122, b1 = -0.037, b1p = -0.006, b2 = 0.999), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b2p not significant
otherHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1*DBH / (1 + DBH)^b1, other2016, start = list(a1 = 1.086, b1 = 0.190), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b1p not significant
otherHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), other2016, start = list(a1 = 34.6, b1 = 40.2, b2 = -1.03), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), other2016, start = list(a1 = 381, b1 = -6.26, b2 = -0.20), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), control = list(maxiter = 50)) # a1p, b1p, b2p not significant
otherHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), other2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1))
otherHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 34.6, a2 = 40.2, b1 = 1.03), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2p, b1p not significant
otherHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), other2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1))
otherHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), other2016, start = list(a1 = 0.028, a2 = 1.08, a3 = 0.15), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2p, a3p not significant
otherHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + a1*DBH^b1, other2016, start = list(a1 = 0.99, b1 = 0.84), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b1p not significant
otherHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), other2016, start = list(a1 = 24.9, b1 = -17.8, b2 = 4.65), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), other2016, start = list(Ha = 12.9, d = 2.59, kU = 0.064, kUp = 0.008), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # Hap, dp not significant
otherHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, other2016, maxit = 20, start = list(a1 = 160, a2 = -0.45, b1 = -0.0065, b2 = -0.41, b3 = 0.87), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2p, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016, start = list(a1 = 90.7, a2 = -0.328, b1 = -0.027, b2 = -0.277, b3 = 0.834), maxit = 20, weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2p, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016physio, start = list(a1 = 40.5, a2 = -0.359, a3 = -0.0004, a4 = -0.015, a5 = -0.020, a6 = -0.0016, b1 = -0.017, b2 = -0.330, b3 = 0.819), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), maxit = 30, control = nls.control(maxiter = 50)) # a2, a3, a4, a5, a6, b1p, b2p, b3p not significant, step factor with nlrob()
otherHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare))^b2*DBH))^b3, other2016physio, start = list(a1 = 168, a2 = -0.47, a3 = -0.0004, a4 = -0.015, a5 = -0.016, a6 = -0.002, b1 = -0.008, b2 = -0.406, b3 = 0.873), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), maxit = 20) # a1p, a2p, a4, a5, a6, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 14.6, a2 = 0.19, b1 = -0.13, b2 = -0.27, b3 = 0.96, b3p = -0.09), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2p, b1p, b2p not significant
otherHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = list(a1 = 62.2, a2 = -0.06, a3 = 0.007, b1 = -0.181, b2 = -0.333, b2p = 0.044, b3 = 0.912), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1), maxit = 50, control = list(maxiter = 50)) # a1p, a2, a2p, a3, a3p, b1p, b3p not significant
otherHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), other2016, start = list(a1 = 0.742, b1 = 1.207, b2 = -0.088), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 24.7, b1 = -0.035, b2 = 0.96, b2p = 0.07), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, b2p not significant
otherHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 45.5, a2 = 0.95, a3 = -0.25, b1 = -0.018, b2 = 0.772, b2p = 0.124), weights = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1)) # a1p, a2, a2p, a3, a3p, b1p not significant
otherHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = -13.4, a2 = -0.484, a3 = 1.36, a4 = 144, b1 = -0.02, b2 = 0.47, b2p = -0.08), weights = pmin(DBH^if_else(isPlantation, -1.6, -1.9), 1), maxit = 50, control = list(maxiter = 50)) # a2, a3p, a4p, b2p not significant, a4 evaporates for plantation weight powers smaller than 1.6
#confint_nlrob(otherHeightFromDiameterWeibullBalRelHt, level = -0.99, weights = pmin(other2016$DBH^if_else(other2016$isPlantation, -1.0, -1.9), 1))

otherHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", otherHeightFromDiameterChapmanRichards, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BAL", otherHeightFromDiameterChapmanRichardsBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BAL physio", otherHeightFromDiameterChapmanRichardsBalPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BAL RelHt", otherHeightFromDiameterChapmanRichardsBalRelHt, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", otherHeightFromDiameterChapmanRichardsPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterCurtis = get_height_error("Curtis", otherHeightFromDiameterCurtis, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", otherHeightFromDiameterHossfeld, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterKorf = get_height_error("Korf", otherHeightFromDiameterKorf, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterLinear = get_height_error("linear", otherHeightFromDiameterLinear, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", otherHeightFromDiameterMichaelisMenten, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterParabolic = get_height_error("parabolic", otherHeightFromDiameterParabolic, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterProdan = get_height_error("Prodan", otherHeightFromDiameterProdan, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterPower = get_height_error("power", otherHeightFromDiameterPower, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", otherHeightFromDiameterRatkowsky, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterRichards = get_height_error("unified Richards", otherHeightFromDiameterRichards, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", otherHeightFromDiameterSharmaParton, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BAL", otherHeightFromDiameterSharmaPartonBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BAL physio", otherHeightFromDiameterSharmaPartonBalPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", otherHeightFromDiameterSharmaPartonPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", otherHeightFromDiameterSharmaZhang, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BAL", otherHeightFromDiameterSharmaZhangBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSibbesen = get_height_error("Sibbesen", otherHeightFromDiameterSibbesen, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibull = get_height_error("Weibull", otherHeightFromDiameterWeibull, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBal = get_height_error("Weibull BAL", otherHeightFromDiameterWeibullBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BAL RelHt", otherHeightFromDiameterWeibullBalRelHt, other2016, other2016natural, other2016plantation)

otherHeightFromDiameterResults = bind_rows(as_row(otherHeightFromDiameterChapmanRichards),
                                           as_row(otherHeightFromDiameterChapmanRichardsBal, significant = FALSE),
                                           as_row(otherHeightFromDiameterChapmanRichardsBalPhysio, significant = FALSE),
                                           as_row(otherHeightFromDiameterChapmanRichardsBalRelHt),
                                           as_row(otherHeightFromDiameterChapmanRichardsPhysio),
                                           as_row(otherHeightFromDiameterCurtis),
                                           as_row(otherHeightFromDiameterHossfeld),
                                           as_row(otherHeightFromDiameterKorf),
                                           as_row(otherHeightFromDiameterLinear),
                                           as_row(otherHeightFromDiameterMichaelisMenten),
                                           as_row(otherHeightFromDiameterParabolic),
                                           as_row(otherHeightFromDiameterPower),
                                           as_row(otherHeightFromDiameterProdan),
                                           as_row(otherHeightFromDiameterRatkowsky),
                                           as_row(otherHeightFromDiameterRichards),
                                           as_row(otherHeightFromDiameterSharmaParton),
                                           as_row(otherHeightFromDiameterSharmaPartonBal, significant = FALSE),
                                           as_row(otherHeightFromDiameterSharmaPartonBalPhysio, significant = FALSE),
                                           as_row(otherHeightFromDiameterSharmaPartonPhysio),
                                           as_row(otherHeightFromDiameterSharmaZhang),
                                           as_row(otherHeightFromDiameterSharmaZhangBal),
                                           as_row(otherHeightFromDiameterSibbesen),
                                           as_row(otherHeightFromDiameterWeibull),
                                           as_row(otherHeightFromDiameterWeibullBal, significant = FALSE),
                                           as_row(otherHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "other", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(otherHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = other2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = other2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterKorf$fitted.values, color = "Korf", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterLinear$fitted.values, color = "linear", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterPower$fitted.values, color = "power", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterProdan$fitted.values, color = "Prodan", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = other2016$isPlantation)) +
  annotate("text", x = 0, y = 70, label = "merged minority species, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 70)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## other height-diameter GNLS regressions
#otherHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3*basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = otherHeightFromDiameterSharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameterChapmanRichards$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, logLik() NaN without
#otherHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameterChapmanRichardsBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, other2016, start = otherHeightFromDiameterSharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016, start = otherHeightFromDiameterSharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, maxIter = 250, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001
#otherHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), other2016, start = otherHeightFromDiameterSharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001
##otherHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3*basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = otherHeightFromDiameterSharmaZhangBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, logLik() NaN without
#otherHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameterWeibullBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, with nlsTol = 0.001 without plot level correlation
#save(otherHeightFromDiameterChapmanRichardsGnls, otherHeightFromDiameterChapmanRichardsBalGnls, otherHeightFromDiameterSharmaPartonGnls, otherHeightFromDiameterSharmaPartonBalGnls, otherHeightFromDiameterSharmaZhangGnls, otherHeightFromDiameterWeibullGnls, otherHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia other GNLS.rdata")

load("trees/height-diameter/HtDia other GNLS.rdata")
#otherHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", otherHeightFromDiameterChapmanRichardsGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BAL GNLS", otherHeightFromDiameterChapmanRichardsBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", otherHeightFromDiameterSharmaPartonGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BAL GNLS", otherHeightFromDiameterSharmaPartonBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", otherHeightFromDiameterSharmaZhangGnls, other2016, other2016natural, other2016plantation)
#otherHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BAL GNLS", otherHeightFromDiameterSharmaZhangBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", otherHeightFromDiameterWeibullGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BAL GNLS", otherHeightFromDiameterWeibullBalGnls, other2016, other2016natural, other2016plantation)

otherHeightFromDiameterResultsGnls = bind_rows(as_row(name = "Chapman-Richards GNLS"),
                                               as_row(otherHeightFromDiameterChapmanRichardsBalGnls),
                                               as_row(otherHeightFromDiameterSharmaPartonGnls),
                                               as_row(otherHeightFromDiameterSharmaPartonBalGnls),
                                               as_row(otherHeightFromDiameterSharmaZhangGnls),
                                               as_row(name = "Sharma-Zhang BAL GNLS"),
                                               as_row(otherHeightFromDiameterWeibullGnls),
                                               as_row(otherHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "other", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
otherHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(otherHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(otherHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(otherHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
#  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
#  select(-bal, -balN, -balP)

ggplot() +
  geom_point(aes(x = other2016natural$DBH, y = other2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = other2016natural$DBH, y = other2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "natural regeneration DBH, cm", y = "minority species naturally regenerated height, m") +
  ggplot() +
  geom_point(aes(x = other2016plantation$DBH, y = other2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = other2016plantation$DBH, y = other2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "plantation DBH, cm", y = "minority species plantation height, m")

ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterWeibullBAL$fitted.values, color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = other2016natural$DBH, y = otherHeightFromDiameterWeibullBALnatural$fitted.values, color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = other2016plantation$DBH, y = otherHeightFromDiameterWeibullBALplantation$fitted.values, color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterBase$fitted.values, color = "base")) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterWeibull$fitted.values, color = "ElliottWeibull")) +
  annotate("text", x = 0, y = 85, label = "a) minority species, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

ggplot(otherConifer2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.20, color = "forestgreen", shape = 16) +
  annotate("text", x = 0, y = 70, label = "other conifers", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 225), ylim = c(0, 70)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.position = "none") +
ggplot(otherHardwood2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.15, color = "green3", shape = 16) +
  annotate("text", x = 0, y = 70, label = "other hardwoods", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 225), ylim = c(0, 70)) +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  theme(legend.position = "none")

## other species diameter-height regressions
#otherDiameterFromHeightChapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, iter = 1000,
#                                              start_lower = list(a1 = 0.1, b1 = -1, b2 = -0.5), 
#                                              start_upper = list(a1 = 150, b1 = 1, b2 = 0.5), modelweights = pmin(TotalHt^-3.0, 0.5))
#otherDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016, iter = 100,
#                                                             start_lower = list(a1 = 1, a1p = -100, a2 = -0.1, a3 = -100, a4 = -1, a5 = -1, a6 = -1, b1 = -0.1, b2 = -2), 
#                                                             start_upper = list(a1 = 200, a1p = 100, a2 = 0.1, a3 = 100, a4 = 10, a5 = 1, a6 = 1, b1 = 0.1, b2 = 2), modelweights = pmin(TotalHt^-3.0, 0.5))
#otherDiameterFromHeightPowerAat = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.0, a1p = -0.38, a2 = -0.013, b1 = 1.06), weights = pmin(TotalHt^-3.0, 0.5)) # a2p, b2p not significant, nlrob() failure to converge >500 steps
#otherDiameterFromHeightRuark = gsl_nls(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 1.56, a1p = 0.90, b1 = 1.05, b1p = -0.57, b2 = 0.0057, b2p = 0.053), weights = pmin(TotalHt^-3.0, 0.5)) # step factor with nlrob()
#otherDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, other2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.01, b2 = 0.26, b3 = 0.5), weights = pmin(TotalHt^-3.0, 0.5), control = nls.control(maxiter = 50)) # step factor
#otherDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -138, b1 = 0.0116, b2 = 1.00), weights = pmin(TotalHt^-3.0, 0.5)) # nlrob() > 20 steps
otherDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 3.4, b1 = 0.476, b2 = 0.224), weights = pmin(TotalHt^-2.4, 0.5), maxit = 50, control = nls.control(maxiter = 50)) # a1p, b1p, b2p not significant, NaN-inf with nls(), no convergence from nls_multstart()
otherDiameterFromHeightChapmanFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 6.85, a2 = 0.019, b1 = 0.13, b2 = 0.52), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5)) # NaN-inf
otherDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.490, a2 = -0.022, a3 = 0.018, b1 = 1.72, b2 = 0.27), maxit = 30, weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 500)) # NaN-inf with nls()
otherDiameterFromHeightChapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.65, a2 = -0.03, a3 = 0.03, a4 = 0.1, b1 = 1.5, b2 = 0.27), weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 50)) # step factor with nls(), step factor or singular gradient with nlrob()
otherDiameterFromHeightChapmanFormRelHt = gsl_nls(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.09, a2 = -0.098, b1 = 1.14, b2 = 0.38), weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 1500)) # step factor with nls()
otherDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -7.17, b1 = 0.319, b2 = 0.415, b2p = -0.0044), weights = pmin(TotalHt^-2.4, 0.5), maxit = 20, control = list(maxiter = 50)) # a1p not significant
otherDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -8.91, a2 = -0.0010, b1 = 0.254, b2 = 0.466, b2p = -0.032), maxit = 150, weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 500)) # a1p, a2p not significant
otherDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = list(a1 = -6.4, a1p = 0.88, a2 = 0.001, a3 = 0.658, a4 = -0.216, a5 = -0.199, a6 = 0.012, b1 = 0.439, b2 = 0.300), weights = pmin(TotalHt^-2.4, 0.5), maxit = 500, control = list(maxiter = 500)) # no physiographic effect significant, convergence fails with b1p, unreliable nlrob() convergence
otherDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -7.17, b1 = 0.319, b2 = 0.415, b2p = -0.0044), weights = pmin(TotalHt^-2.4, 0.5), maxit = 20, control = list(maxiter = 50)) # a1p not significant
otherDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -21.3, a2 = 1.5, b1 = 0.11, b2 = 0.60, b2p = 0.006), maxit = 20, weights = pmin(TotalHt^-2.4, 0.5)) # a1p, a2p not significant, nlrob() fails to converge from closer positions
otherDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), other2016, weights = pmin(TotalHt^-2.4, 0.5))
otherDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = list(a1 = 28.7, a2 = 13.0, b1 = 0.60), weights = pmin(TotalHt^-2.4, 0.5), maxit = 30, control = nls.control(maxiter = 50)) # a1p, a2p, b1p not significant
otherDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 0.37, a1p = -0.27, a2 = -0.14, a2p = -0.038), weights = pmin(TotalHt^-2.4, 0.5)) # converges poorly, yields negative values
otherDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2), other2016, weights = pmin(TotalHt^-2.4, 0.5)) # isPlantation*(TotalHt - 1.37)^2 not significant
otherDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 1.36, a1p = -0.48, b1 = 1.15, b1p = 0.11), weights = pmin(TotalHt^-2.4, 0.5))
otherDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.08, a2 = -0.0004, b1 = 0.81), weights = pmin(TotalHt^-2.4, 0.5), maxit = 500) # a2p, b2p not significant, nlrob() failure to converge >500 steps with a1p
otherDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = list(a1 = 1.81, a2 = -0.0007, a3 = -1.05, a4 = 0.034, a5 = -0.024, a6 = 0.00049, b1 = 1.21), weights = pmin(TotalHt^-2.4, 0.5)) # a1p, a4, a5, a6 not significant
otherDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 1.32, a2 = 0.40, b1 = 1.11, b1p = -0.055), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5))
otherDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 2.72, b1 = 0.20, b1p = -0.04, b2 = 0.09, b2p = 0.009), maxit = 40, weights = pmin(TotalHt^-2.4, 0.5)) # step factor with a1p
otherDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.20, Ha = 187), weights = pmin(TotalHt^-2.4, 0.5)) # NaN-inf with nlrob()
otherDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, other2016, start = list(a1 = 8500, a2 = -1.2, b1 = 0.02, b2 = -0.05, b3 = 2.0), weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 50)) # step size or singular gradient with nls(), NaN-inf with nlrob()
otherDiameterFromHeightSibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.67, a1p = 1.84, b1 = 1.70, b1p = -1.26, b2 = -0.063, b2p = 0.30), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5))
otherDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.58, a1p = 1.93, a2 = -0.0005, b1 = 1.93, b1p = -1.49, b2 = -0.085, b2p = 0.32), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5)) # a2 not significant
otherDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = list(a1 = 3.47, a1p = -0.43, a2 = -0.0004, a3 = -0.53, a4 = 0.003, a5 = 0.106, a6 = -0.001, b1 = 0.34, b2 = 0.28, b2p = 0.019), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5)) # a4, a5, a6 not significant
otherDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 2.84, a2 = -0.007, b1 = 0.407, b1p = -0.121, b2 = 0.260, b2p = 0.122), maxit = 50, weights = pmin(TotalHt^-2.4, 0.5)) # a1p, a2p not significant
otherDiameterFromHeightWeibull = nlrob(DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -153, a1p = 37.0, b1 = 0.09, b2 = 0.42), weights = pmin(TotalHt^-2.4, 0.5), control = nls.control(maxiter = 50)) # NaN-inf with b1p
#confint2(otherDiameterFromHeightWeibull, level = 0.99)

otherDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", otherDiameterFromHeightChapmanForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AAT", otherDiameterFromHeightChapmanFormAat, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BAL", otherDiameterFromHeightChapmanFormBal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BAL RelHt", otherDiameterFromHeightChapmanFormBalRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", otherDiameterFromHeightChapmanFormRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", otherDiameterFromHeightChapmanRichards, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AAT", otherDiameterFromHeightChapmanRichardsAat, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", otherDiameterFromHeightChapmanRichardsPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", otherDiameterFromHeightChapmanRichardsRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightLinear = get_dbh_error("linear", otherDiameterFromHeightLinear, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", otherDiameterFromHeightMichaelisMentenForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightNaslund = get_dbh_error("Näslund", otherDiameterFromHeightNaslund, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightParabolic = get_dbh_error("parabolic", otherDiameterFromHeightParabolic, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPower = get_dbh_error("power", otherDiameterFromHeightPower, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerAat = get_dbh_error("power AAT", otherDiameterFromHeightPowerAat, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerPhysio = get_dbh_error("power physio", otherDiameterFromHeightPowerPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", otherDiameterFromHeightPowerRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightRuark = get_dbh_error("Ruark", otherDiameterFromHeightRuark, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSchnute = get_dbh_error("Schnute", otherDiameterFromHeightSchnute, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", otherDiameterFromHeightSharmaParton, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", otherDiameterFromHeightSibbesenForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AAT", otherDiameterFromHeightSibbesenFormAat, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", otherDiameterFromHeightSibbesenFormPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", otherDiameterFromHeightSibbesenFormRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightWeibull = get_dbh_error("Weibull", otherDiameterFromHeightWeibull, other2016, other2016natural, other2016plantation)

otherDiameterFromHeightResults = bind_rows(as_row(otherDiameterFromHeightChapmanRichards),
                                           as_row(otherDiameterFromHeightChapmanRichardsAat),
                                           as_row(otherDiameterFromHeightChapmanRichardsPhysio),
                                           as_row(otherDiameterFromHeightChapmanRichardsRelHt),
                                           as_row(otherDiameterFromHeightChapmanForm),
                                           as_row(otherDiameterFromHeightChapmanFormAat),
                                           as_row(otherDiameterFromHeightChapmanFormBal),
                                           as_row(otherDiameterFromHeightChapmanFormBalRelHt),
                                           as_row(otherDiameterFromHeightChapmanFormRelHt),
                                           as_row(otherDiameterFromHeightLinear),
                                           as_row(otherDiameterFromHeightMichaelisMentenForm),
                                           as_row(otherDiameterFromHeightNaslund),
                                           as_row(otherDiameterFromHeightParabolic),
                                           as_row(otherDiameterFromHeightPower),
                                           as_row(otherDiameterFromHeightPowerAat),
                                           as_row(otherDiameterFromHeightPowerPhysio),
                                           as_row(otherDiameterFromHeightPowerRelHt),
                                           as_row(otherDiameterFromHeightRuark),
                                           as_row(otherDiameterFromHeightSchnute),
                                           as_row(otherDiameterFromHeightSharmaParton),
                                           as_row(otherDiameterFromHeightSibbesenForm),
                                           as_row(otherDiameterFromHeightSibbesenFormAat),
                                           as_row(otherDiameterFromHeightSibbesenFormPhysio),
                                           as_row(otherDiameterFromHeightSibbesenFormRelHt),
                                           as_row(otherDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(otherDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(other2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = otherDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = otherDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form aBAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 5*(TotalHt - 1.37)^0.5*(exp(0.01*(tph/topHeight)^0.26*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 200*(TotalHt - 1.37)^0.9/(100 - (TotalHt - 1.37)^0.9), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 90, label = "other species, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

#ggplot(other2016) +
#  geom_point(aes(x = TotalHt, y = abs(residuals(otherDiameterFromHeightMichaelisMentenForm)), color = "Michaelis-Menten form", group = isPlantation), alpha = 0.1, color = "grey25", shape= 16)

## other diameter-height GNLS regressions
#otherDiameterFromHeightNaslundGnls = gnls(DBH ~ a1* sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 2.36, a2 = -0.13), weights = varPower(1.4, ~TotalHt | isPlantation))
otherDiameterFromHeightChapmanFormGnls = gnls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeightChapmanForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()) 
otherDiameterFromHeightChapmanFormAatGnls = gnls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeightChapmanFormAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightChapmanFormBalGnls = gnls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeightChapmanFormBal$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanFormBalRelHtGnls = gnls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeightChapmanFormBalRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanFormRelHtGnls = gnls(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeightChapmanFormRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanRichardsGnls = gnls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeightChapmanRichards$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanRichardsAatGnls = gnls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeightChapmanRichardsAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanRichardsPhysioGnls = gnls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = otherDiameterFromHeightChapmanRichardsPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanRichardsGnls = gnls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeightChapmanRichards$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightChapmanRichardsRelHtGnls = gnls(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeightChapmanRichardsRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightMichaelisMentenFormGnls = gnls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = otherDiameterFromHeightMichaelisMentenForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightNaslundGnls = gnls(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 2.28, a1p = 0, a2 = -0.13, a2p = 0), weights = varPower(1.4, ~TotalHt | isPlantation)) # step halving with or without isPlantation
otherDiameterFromHeightPowerGnls = gnls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeightPower$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightPowerAatGnls = gnls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, other2016, start = otherDiameterFromHeightPowerAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightPowerPhysioGnls = gnls(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = otherDiameterFromHeightPowerPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightPowerRelHtGnls = gnls(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeightPowerRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightRuarkGnls = gnls(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = otherDiameterFromHeightRuark$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
#otherDiameterFromHeightSchnuteGnls = gnls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = otherDiameterFromHeightSchnute$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation)) # step halving
#otherDiameterFromHeightSharmaPartonGnls = gnls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, other2016, start = otherDiameterFromHeightSharmaParton$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()) # NaN
otherDiameterFromHeightSibbesenFormGnls = gnls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeightSibbesenForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightSibbesenFormAatGnls = gnls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeightSibbesenFormAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightSibbesenFormPhysioGnls = gnls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = otherDiameterFromHeightSibbesenFormPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightSibbesenFormRelHtGnls = gnls(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeightSibbesenFormRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightWeibullGnls = gnls(DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = otherDiameterFromHeightWeibull$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())

otherDiameterFromHeightChapmanFormGnls = get_dbh_error("Chapman-Richards form GNLS", otherDiameterFromHeightChapmanFormGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormAatGnls = get_dbh_error("Chapman-Richards form AAT GNLS", otherDiameterFromHeightChapmanFormAatGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBalGnls = get_dbh_error("Chapman-Richards form BAL GNLS", otherDiameterFromHeightChapmanFormBalGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBalRelHtGnls = get_dbh_error("Chapman-Richards form BAL RelHt GNLS", otherDiameterFromHeightChapmanFormBalRelHtGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormRelHtGnls = get_dbh_error("Chapman-Richards form RelHt GNLS", otherDiameterFromHeightChapmanFormRelHtGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsGnls = get_dbh_error("Chapman-Richards GNLS", otherDiameterFromHeightChapmanRichardsGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsAatGnls = get_dbh_error("Chapman-Richards AAT GNLS", otherDiameterFromHeightChapmanRichardsAatGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsPhysioGnls = get_dbh_error("Chapman-Richards physio GNLS", otherDiameterFromHeightChapmanRichardsPhysioGnls, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightChapmanRichardsRelHtGnls = get_dbh_error("Chapman-Richards RelHt GNLS", otherDiameterFromHeightChapmanRichardsRelHtGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightMichaelisMentenFormGnls = get_dbh_error("Michaelis-Menten form GNLS", otherDiameterFromHeightMichaelisMentenFormGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightNaslundGnls = get_dbh_error("Näslund GNLS", otherDiameterFromHeightNaslundGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerGnls = get_dbh_error("power GNLS", otherDiameterFromHeightPowerGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerAatGnls = get_dbh_error("power AAT GNLS", otherDiameterFromHeightPowerAatGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerPhysioGnls = get_dbh_error("power physio GNLS", otherDiameterFromHeightPowerPhysioGnls, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightPowerRelHtGnls = get_dbh_error("power RelHt GNLS", otherDiameterFromHeightPowerRelHtGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightRuarkGnls = get_dbh_error("Ruark GNLS", otherDiameterFromHeightRuarkGnls, other2016, other2016natural, other2016plantation)
#otherDiameterFromHeightSchnuteGnls = get_dbh_error("Schnute GNLS", otherDiameterFromHeightSchnuteGnls, other2016, other2016natural, other2016plantation)
#otherDiameterFromHeightSharmaPartonGnls = get_dbh_error("modified Sharma-Parton GNLS", otherDiameterFromHeightSharmaPartonGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormGnls = get_dbh_error("Sibbesen form GNLS", otherDiameterFromHeightSibbesenFormGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormAatGnls = get_dbh_error("Sibbesen form AAT GNLS", otherDiameterFromHeightSibbesenFormAatGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormPhysioGnls = get_dbh_error("Sibbesen form physio GNLS", otherDiameterFromHeightSibbesenFormPhysioGnls, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightSibbesenFormRelHtGnls = get_dbh_error("Sibbesen form RelHt GNLS", otherDiameterFromHeightSibbesenFormRelHtGnls, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightWeibullGnls = get_dbh_error("Weibull GNLS", otherDiameterFromHeightWeibullGnls, other2016, other2016natural, other2016plantation)

otherDiameterFromHeightResultsGnls = bind_rows(as_row(otherDiameterFromHeightChapmanFormGnls),
                                               as_row(otherDiameterFromHeightChapmanFormAatGnls),
                                               as_row(otherDiameterFromHeightChapmanFormBalGnls),
                                               as_row(otherDiameterFromHeightChapmanFormRelHtGnls),
                                               as_row(otherDiameterFromHeightChapmanRichardsGnls),
                                               as_row(otherDiameterFromHeightChapmanRichardsAatGnls),
                                               as_row(otherDiameterFromHeightChapmanRichardsPhysioGnls),
                                               as_row(otherDiameterFromHeightChapmanRichardsRelHtGnls),
                                               as_row(otherDiameterFromHeightMichaelisMentenFormGnls),
                                               as_row(otherDiameterFromHeightNaslundGnls),
                                               as_row(otherDiameterFromHeightPowerGnls),
                                               as_row(otherDiameterFromHeightPowerAatGnls),
                                               as_row(otherDiameterFromHeightPowerPhysioGnls),
                                               as_row(otherDiameterFromHeightPowerRelHtGnls),                                           as_row(otherDiameterFromHeightRuarkGnls),
                                               as_row(otherDiameterFromHeightRuarkGnls),
                                               #as_row(otherDiameterFromHeightSchnuteGnls),
                                               #as_row(otherDiameterFromHeightSharmaPartonGnls),
                                               as_row(otherDiameterFromHeightSibbesenFormGnls),
                                               as_row(otherDiameterFromHeightSibbesenFormAatGnls),
                                               as_row(otherDiameterFromHeightSibbesenFormPhysioGnls),
                                               as_row(otherDiameterFromHeightSibbesenFormRelHtGnls),
                                               as_row(otherDiameterFromHeightWeibullGnls)) %>%
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

ggplot() +
  geom_histogram(aes(x = power), otherDiameterFromHeightResultsGnls, binwidth = 0.1) +
  geom_segment(aes(x = mean(otherDiameterFromHeightResultsGnls$power), xend = mean(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash") +
  geom_segment(aes(x = median(otherDiameterFromHeightResultsGnls$power), xend = median(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash")
  
otherDiameterFromHeightResultsGnls %>% summarize(power = mean(power), powerPlantation = mean(powerPlantation))

## collect model parameters
otherParameters = bind_rows(bind_rows(get_coefficients(otherHeightFromDiameterChapmanRichards),
                                      get_coefficients(otherHeightFromDiameterChapmanRichardsBal),
                                      get_coefficients(otherHeightFromDiameterChapmanRichardsBalPhysio),
                                      get_coefficients(otherHeightFromDiameterChapmanRichardsBalRelHt),
                                      get_coefficients(otherHeightFromDiameterChapmanRichardsPhysio),
                                      get_coefficients(otherHeightFromDiameterCurtis),
                                      get_coefficients(otherHeightFromDiameterHossfeld),
                                      get_coefficients(otherHeightFromDiameterKorf),
                                      get_coefficients(otherHeightFromDiameterLinear),
                                      get_coefficients(otherHeightFromDiameterMichaelisMenten),
                                      get_coefficients(otherHeightFromDiameterParabolic),
                                      get_coefficients(otherHeightFromDiameterPower),
                                      get_coefficients(otherHeightFromDiameterProdan),
                                      get_coefficients(otherHeightFromDiameterRatkowsky),
                                      get_coefficients(otherHeightFromDiameterRichards),
                                      get_coefficients(otherHeightFromDiameterSharmaParton),
                                      get_coefficients(otherHeightFromDiameterSharmaPartonBal),
                                      get_coefficients(otherHeightFromDiameterSharmaPartonBalPhysio),
                                      get_coefficients(otherHeightFromDiameterSharmaPartonPhysio),
                                      get_coefficients(otherHeightFromDiameterSharmaZhang),
                                      get_coefficients(otherHeightFromDiameterSharmaZhangBal),
                                      get_coefficients(otherHeightFromDiameterSibbesen),
                                      get_coefficients(otherHeightFromDiameterWeibull),
                                      get_coefficients(otherHeightFromDiameterWeibullBal),
                                      get_coefficients(otherHeightFromDiameterWeibullBalRelHt),
                                      #get_coefficients(otherHeightFromDiameterChapmanRichardsGnls),
                                      get_coefficients(otherHeightFromDiameterChapmanRichardsBalGnls),
                                      get_coefficients(otherHeightFromDiameterSharmaPartonGnls),
                                      get_coefficients(otherHeightFromDiameterSharmaPartonBalGnls),
                                      get_coefficients(otherHeightFromDiameterSharmaZhangGnls),
                                      #get_coefficients(otherHeightFromDiameterSharmaZhangBalGnls),
                                      get_coefficients(otherHeightFromDiameterWeibullGnls),
                                      get_coefficients(otherHeightFromDiameterWeibullBalGnls)) %>%
                              mutate(responseVariable = "height"),
                            bind_rows(get_coefficients(otherDiameterFromHeightChapmanRichards),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsAat),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsPhysio),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsRelHt),
                                      get_coefficients(otherDiameterFromHeightChapmanForm),
                                      get_coefficients(otherDiameterFromHeightChapmanFormAat),
                                      get_coefficients(otherDiameterFromHeightChapmanFormBal),
                                      get_coefficients(otherDiameterFromHeightChapmanFormBalRelHt),
                                      get_coefficients(otherDiameterFromHeightChapmanFormRelHt),
                                      get_coefficients(otherDiameterFromHeightLinear),
                                      get_coefficients(otherDiameterFromHeightMichaelisMentenForm),
                                      get_coefficients(otherDiameterFromHeightNaslund),
                                      get_coefficients(otherDiameterFromHeightParabolic),
                                      get_coefficients(otherDiameterFromHeightPower),
                                      get_coefficients(otherDiameterFromHeightPowerAat),
                                      get_coefficients(otherDiameterFromHeightPowerPhysio),
                                      get_coefficients(otherDiameterFromHeightPowerRelHt),
                                      get_coefficients(otherDiameterFromHeightRuark),
                                      get_coefficients(otherDiameterFromHeightSchnute),
                                      get_coefficients(otherDiameterFromHeightSharmaParton),
                                      get_coefficients(otherDiameterFromHeightSibbesenForm),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormAat),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormPhysio),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormRelHt),
                                      get_coefficients(otherDiameterFromHeightWeibull),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsAatGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsPhysioGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanRichardsRelHtGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanFormGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanFormAatGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanFormBalGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanFormBalRelHtGnls),
                                      get_coefficients(otherDiameterFromHeightChapmanFormRelHtGnls),
                                      get_coefficients(otherDiameterFromHeightMichaelisMentenFormGnls),
                                      get_coefficients(otherDiameterFromHeightNaslundGnls),
                                      get_coefficients(otherDiameterFromHeightPowerGnls),
                                      get_coefficients(otherDiameterFromHeightPowerAatGnls),
                                      get_coefficients(otherDiameterFromHeightPowerPhysioGnls),
                                      get_coefficients(otherDiameterFromHeightPowerRelHtGnls),
                                      get_coefficients(otherDiameterFromHeightRuarkGnls),
                                      #get_coefficients(otherDiameterFromHeightSchnuteGnls),
                                      #get_coefficients(otherDiameterFromHeightSharmaPartonGnls),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormGnls),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormAatGnls),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormPhysioGnls),
                                      get_coefficients(otherDiameterFromHeightSibbesenFormRelHtGnls),
                                      get_coefficients(otherDiameterFromHeightWeibullGnls)) %>%
                              mutate(responseVariable = "DBH")) %>%
  mutate(species = "other",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
otherBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), other2016, start = list(a1 = 1.36, b1 = 0.0002, b2 = 2.06, b2p = -0.27), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 50)) # a1p, b1p not significant, step factor with nlrob()
otherBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p*isPlantation), other2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = -0.00001, b1 = 2.53, b1p = -0.435), maxit = 70, weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 50))
#confint2(otherBasalAreaFromHeightKorf, level = 0.99)
#confint_nlrob(otherBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/other2016$basalArea, 1E4))

otherBasalAreaFromHeightKorf$fitted.values = predict(otherBasalAreaFromHeightKorf, other2016)
otherBasalAreaFromHeightKorf$residuals = otherBasalAreaFromHeightKorf$fitted.values - other2016$basalArea
otherBasalAreaFromHeightPower$fitted.values = predict(otherBasalAreaFromHeightPower, other2016)
otherBasalAreaFromHeightPower$residuals = otherBasalAreaFromHeightPower$fitted.values - other2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(otherBasalAreaFromHeightKorf), 100^2 * mean(otherBasalAreaFromHeightKorf$residuals), mean(abs(otherBasalAreaFromHeightKorf$residuals)), 1 - sum(otherBasalAreaFromHeightKorf$residuals^2) / sum((other2016$basalArea - mean(other2016$basalArea)^2)),
        "power", AIC(otherBasalAreaFromHeightPower), 100^2 * mean(otherBasalAreaFromHeightPower$residuals), mean(abs(otherBasalAreaFromHeightPower$residuals)), 1 - sum(otherBasalAreaFromHeightPower$residuals^2) / sum((other2016$basalArea - mean(other2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(other2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = otherBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = otherBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "minority species height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
