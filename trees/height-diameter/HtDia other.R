# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## minority species height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Parton, Sharma-Zhang, Chapman-Richards BAL
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 43.0, a1p = 13.5, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875), weights = pmin(DBH^-2, 1), control = list(maxiter = 50)) # a1p not significant
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 43.0, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875, b2p = 0), weights = pmin(DBH^-2, 1), control = list(maxiter = 500)) # > 500 iterations
#otherHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 89.3, a2 = 1.166, a2p = 0, a3 = -0.039, a3p = 0, b1 = -0.00544, b2 = 0.873), weights = pmin(DBH^-2, 1)) # a2p, a3p not significant
#otherHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 58.1, a1p = -47.4, a2 = 0.162, a2p = 0.032, b1 = -0.021, b1p = -0.222, b2 = -0.292, b2p = 0.036, b3 = 0.818, b3p = 0.165), weights = pmin(DBH^-2, 1))
#otherHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.714, a1p = 0.088, b1 = 1.172, b2 = -0.074, b2p = -0.0040), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
other2016physio = other2016 %>% filter(is.na(elevation) == FALSE)
other2016plantationPhysio = other2016physio %>% filter(isPlantation)
otherHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016, start = list(a1 = 24.7, b1 = -0.033, b1p = -0.003, b2 = 1.000), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 81.9, a2 = 1.97, a3 = -1.23, b1 = -0.006, b2 = 0.881), weights = pmin(DBH^-2, 1), control = list(maxiter = 50)) # a1p, a2p, a2, a3p, a3, b1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^b2, other2016physio, start = list(a1 = 27.2, a2 = 0.103, a3 = -00.10, a4 = -0.085, a5 = 0.045, a6 = 0.768, b1 = -0.034, b2 = 0.953), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, a3, a4, a5, a6, b1p, b2p not significant
otherHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = -1.846, a2 = -0.007, a3 = 0.019, a3p = 0.579, a4 = 64.0, a4p = 31.3, b1 = -0.002, b2 = 0.035, b2p = 0.444), weights = pmin(DBH^-2, 1)) # a2, a2p, a3, b1p not significant
otherHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016physio, start = list(a1 = 27.0, a2 = -0.015, a3 = -0.025, a4 = 0.276, a5 = -0.514, a6 = -0.122, b1 = -0.037, b1p = -0.006, b2 = 0.999), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
otherHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1*DBH / (1 + DBH)^b1, other2016, start = list(a1 = 1.086, b1 = 0.190), weights = pmin(DBH^-2, 1)) # a1p, b1p not significant
otherHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), other2016, start = list(a1 = 34.6, b1 = 40.2, b2 = -1.03), weights = pmin(DBH^-2, 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), other2016, start = list(a1 = 381, b1 = -6.26, b2 = -0.20), weights = pmin(DBH^-2, 1), control = list(maxiter = 50)) # a1p, b1p, b2p not significant
otherHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), other2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
otherHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 34.6, a2 = 40.2, b1 = 1.03), weights = pmin(DBH^-2, 1)) # a1p, a2p, b1p not significant
otherHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), other2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
otherHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), other2016, start = list(a1 = 0.028, a2 = 1.08, a3 = 0.15), weights = pmin(DBH^-2, 1)) # a1p, a2p, a3p not significant
otherHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + a1*DBH^b1, other2016, start = list(a1 = 0.99, b1 = 0.84), weights = pmin(DBH^-2, 1)) # a1p, b1p not significant
otherHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), other2016, start = list(a1 = 24.9, b1 = -17.8, b2 = 4.65), weights = pmin(DBH^-2, 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), other2016, start = list(Ha = 12.9, d = 2.59, kU = 0.064, kUp = 0.008), weights = pmin(DBH^-2, 1)) # Hap, dp not significant
otherHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, other2016, maxit = 20, start = list(a1 = 160, a2 = -0.45, b1 = -0.0065, b2 = -0.41, b3 = 0.87), weights = pmin(DBH^-2, 1)) # a1p, a2p, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016, start = list(a1 = 264, a2 = -0.465, b1 = -0.029, b2 = -0.436, b3 = 0.834), weights = pmin(DBH^-2, 1)) # a1p, a2p, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaPartonBalPhysio = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016physio, start = list(a1 = 478, a2 = -0.108, a3 = -0.00004, a4 = -0.0169, a5 = -0.0132, a6 = -0.00165, b1 = -0.0016, b2 = -0.243, b3 = 0.808), weights = pmin(DBH^-2, 1), control = nls.control(maxiter = 50)) # a2, a3, a4, a5, a6, b1p, b2p, b3p not significant, step factor with nlrob()
otherHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare))^b2*DBH))^b3, other2016physio, start = list(a1 = 168, a2 = -0.47, a3 = -0.0004, a4 = -0.015, a5 = -0.016, a6 = -0.002, b1 = -0.008, b2 = -0.406, b3 = 0.873), weights = pmin(DBH^-2, 1), maxit = 20) # a1p, a2p, a4, a5, a6, b1p, b2p, b3p not significant
otherHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 14.6, a2 = 0.19, b1 = -0.13, b2 = -0.27, b3 = 0.96, b3p = -0.09), weights = pmin(DBH^-2, 1)) # a1p, a2p, b1p, b2p not significant
otherHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = list(a1 = 16.6, a2 = 0.098, a3 = 0.007, b1 = -0.181, b2 = -0.333, b2p = 0.044, b3 = 0.912), weights = pmin(DBH^-2, 1), control = list(maxiter = 50)) # a1p, a2, a2p, a3, a3p, b1p, b3p not significant
otherHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), other2016, start = list(a1 = 0.742, b1 = 1.207, b2 = -0.088), weights = pmin(DBH^-2, 1)) # a1p, b1p, b2p not significant
otherHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 24.7, b1 = -0.035, b2 = 0.96, b2p = 0.07), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
otherHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 30.9, a2 = 0.60, a3 = -0.210, b1 = -0.025, b2 = 0.857, b2p = 0.123), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, a3, a3p, b1p not significant
otherHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), other2016, start = list(a1 = -1.74, a1p = 2.36, a2 = -0.016, a2p = 0.393, a3 = 0.027, a4 = 57.1, b1 = -1.56, b1p = 1.47, b2 = 0.59), weights = pmin(DBH^-2, 1), control = list(maxiter = 50)) # a2, a3p, a4p, b2p not significant
#confint2(otherHeightFromDiameterWeibullBalRelHt, level = -0.99)

otherHeightFromDiameterChapmanRichards = get_height_error(otherHeightFromDiameterChapmanRichards, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBal = get_height_error(otherHeightFromDiameterChapmanRichardsBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(otherHeightFromDiameterChapmanRichardsBalPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(otherHeightFromDiameterChapmanRichardsBalRelHt, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsPhysio = get_height_error(otherHeightFromDiameterChapmanRichardsPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterCurtis = get_height_error(otherHeightFromDiameterCurtis, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterHossfeld = get_height_error(otherHeightFromDiameterHossfeld, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterKorf = get_height_error(otherHeightFromDiameterKorf, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterLinear = get_height_error(otherHeightFromDiameterLinear, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterMichaelisMenten = get_height_error(otherHeightFromDiameterMichaelisMenten, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterParabolic = get_height_error(otherHeightFromDiameterParabolic, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterProdan = get_height_error(otherHeightFromDiameterProdan, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterPower = get_height_error(otherHeightFromDiameterPower, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterRatkowsky = get_height_error(otherHeightFromDiameterRatkowsky, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterRichards = get_height_error(otherHeightFromDiameterRichards, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaParton = get_height_error(otherHeightFromDiameterSharmaParton, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBal = get_height_error(otherHeightFromDiameterSharmaPartonBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBalPhysio = get_height_error(otherHeightFromDiameterSharmaPartonBalPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterSharmaPartonPhysio = get_height_error(otherHeightFromDiameterSharmaPartonPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherHeightFromDiameterSharmaZhang = get_height_error(otherHeightFromDiameterSharmaZhang, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaZhangBal = get_height_error(otherHeightFromDiameterSharmaZhangBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSibbesen = get_height_error(otherHeightFromDiameterSibbesen, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibull = get_height_error(otherHeightFromDiameterWeibull, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBal = get_height_error(otherHeightFromDiameterWeibullBal, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBalRelHt = get_height_error(otherHeightFromDiameterWeibullBalRelHt, other2016, other2016natural, other2016plantation)

otherHeightFromDiameterResults = bind_rows(as_row("Chapman-Richards", otherHeightFromDiameterChapmanRichards),
                                           as_row("Chapman-Richards BAL", otherHeightFromDiameterChapmanRichardsBal, significant = FALSE),
                                           as_row("Chapman-Richards BAL physio", otherHeightFromDiameterChapmanRichardsBalPhysio, significant = FALSE),
                                           as_row("Chapman-Richards BAL RelHt", otherHeightFromDiameterChapmanRichardsBalRelHt),
                                           as_row("Chapman-Richards physio", otherHeightFromDiameterChapmanRichardsPhysio),
                                           as_row("Curtis", otherHeightFromDiameterCurtis),
                                           as_row("Hossfeld", otherHeightFromDiameterHossfeld),
                                           as_row("Korf", otherHeightFromDiameterKorf),
                                           as_row("linear", otherHeightFromDiameterLinear),
                                           as_row("generalized Michaelis-Menten", otherHeightFromDiameterMichaelisMenten),
                                           as_row("parabolic", otherHeightFromDiameterParabolic),
                                           as_row("power", otherHeightFromDiameterPower),
                                           as_row("Prodan", otherHeightFromDiameterProdan),
                                           as_row("Ratkowsky", otherHeightFromDiameterRatkowsky),
                                           as_row("unified Richards", otherHeightFromDiameterRichards),
                                           as_row("Sharma-Parton", otherHeightFromDiameterSharmaParton),
                                           as_row("Sharma-Parton BAL", otherHeightFromDiameterSharmaPartonBal, significant = FALSE),
                                           as_row("Sharma-Parton BAL physio", otherHeightFromDiameterSharmaPartonBalPhysio, significant = FALSE),
                                           as_row("Sharma-Parton physio", otherHeightFromDiameterSharmaPartonPhysio),
                                           as_row("Sharma-Zhang", otherHeightFromDiameterSharmaZhang),
                                           as_row("Sharma-Zhang BAL", otherHeightFromDiameterSharmaZhangBal),
                                           as_row("Sibbesen", otherHeightFromDiameterSibbesen),
                                           as_row("Weibull", otherHeightFromDiameterWeibull),
                                           as_row("Weibull BAL", otherHeightFromDiameterWeibullBal, significant = FALSE),
                                           as_row("Weibull BAL RelHt", otherHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic)) %>%
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
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = other2016$isPlantation)) +
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
#otherHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#otherHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(otherHeightFromDiameterChapmanRichardsGnls, otherHeightFromDiameterChapmanRichardsBalGnls, otherHeightFromDiameterSharmaPartonGnls, otherHeightFromDiameterSharmaPartonBalGnls, otherHeightFromDiameterSharmaZhangGnls, otherHeightFromDiameterSharmaZhangBalGnls, otherHeightFromDiameterWeibullGnls, otherHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia other GNLS.rdata")
load("trees/height-diameter/HtDia other GNLS.rdata")
otherHeightFromDiameterWeibullGnls = otherHeightFromDiameterWykoffGnls # temporary naming error fixup
otherHeightFromDiameterWeibullBalGnls = otherHeightFromDiameterWykoffBalGnls

otherHeightFromDiameterChapmanRichardsGnls = get_height_error(otherHeightFromDiameterChapmanRichardsGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterChapmanRichardsBalGnls = get_height_error(otherHeightFromDiameterChapmanRichardsBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonGnls = get_height_error(otherHeightFromDiameterSharmaPartonGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaPartonBalGnls = get_height_error(otherHeightFromDiameterSharmaPartonBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaZhangGnls = get_height_error(otherHeightFromDiameterSharmaZhangGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterSharmaZhangBalGnls = get_height_error(otherHeightFromDiameterSharmaZhangBalGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullGnls = get_height_error(otherHeightFromDiameterWeibullGnls, other2016, other2016natural, other2016plantation)
otherHeightFromDiameterWeibullBalGnls = get_height_error(otherHeightFromDiameterWeibullBalGnls, other2016, other2016natural, other2016plantation)

otherHeightFromDiameterResultsGnls = bind_rows(as_row("Chapman-Richards GNLS", otherHeightFromDiameterChapmanRichardsGnls),
                                               as_row("Chapman-Richards BAL GNLS", otherHeightFromDiameterChapmanRichardsBalGnls),
                                               as_row("Sharma-Parton GNLS", otherHeightFromDiameterSharmaPartonGnls),
                                               as_row("Sharma-Parton BAL GNLS", otherHeightFromDiameterSharmaPartonBalGnls),
                                               as_row("Sharma-Zhang GNLS", otherHeightFromDiameterSharmaZhangGnls),
                                               as_row("Sharma-Zhang BAL GNLS", otherHeightFromDiameterSharmaZhangBalGnls),
                                               as_row("Weibull GNLS", otherHeightFromDiameterWeibullGnls),
                                               as_row("Weibull BAL GNLS", otherHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic)) %>%
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
#                                              start_upper = list(a1 = 150, b1 = 1, b2 = 0.5), modelweights = pmin(TotalHt^-2, 0.5))
#otherDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999999)), other2016, iter = 100,
#                                                             start_lower = list(a1 = 1, a1p = -100, a2 = -0.1, a3 = -100, a4 = -1, a5 = -1, a6 = -1, b1 = -0.1, b2 = -2), 
#                                                             start_upper = list(a1 = 200, a1p = 100, a2 = 0.1, a3 = 100, a4 = 10, a5 = 1, a6 = 1, b1 = 0.1, b2 = 2), modelweights = pmin(TotalHt^-2, 0.5))
otherDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 7.88, b1 = 0.11, b2 = 0.54), weights = pmin(TotalHt^-2, 0.5), maxit = 20, control = nls.control(maxiter = 50)) # a1p, b1p, b2p not significant, NaN-inf with nls(), no convergence from nls_multstart()
otherDiameterFromHeightChapmanFormAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 6.85, a2 = 0.019, b1 = 0.13, b2 = 0.52), maxit = 20, weights = pmin(TotalHt^-2, 0.5)) # NaN-inf
otherDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.03, a2 = -0.0025, b1 = 1.20, b2 = 0.36), maxit = 50, weights = pmin(TotalHt^-2, 0.5), control = nls.control(maxiter = 50)) # NaN-inf with nls()
otherDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.51, a2 = -0.0184, a3 = 0.0140, a4 = -0.199, b1 = 1.75, b2 = 0.286), weights = pmin(TotalHt^-2, 0.5)) # step factor with nls()
otherDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.12, a2 = -0.458, b1 = 1.13, b2 = 0.384), weights = pmin(TotalHt^-2, 0.5), control = nls.control(maxiter = 50)) # step factor with nls()
otherDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), other2016, start = list(a1 = -7.8, b1 = 0.344, b2 = 0.344, b2p = -0.005), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50)) # a1p not significant
otherDiameterFromHeightChapmanRichardsAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), other2016, start = list(a1 = -6.6, a2 = 0.0004, b1 = 0.398, b2 = 0.314, b2p = -0.023), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p not significant
otherDiameterFromHeightChapmanRichardsPhysio = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999999)), other2016physio, start = list(a1 = -200, a1p = 17.5, a2 = 0.046, a3 = 98, a4 = -1.1, a5 = -3.9, a6 = -0.35, b1 = 0.014, b2 = 0.96), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50)) # no physiographic effect significant, convergence fails with b1p, nlrob() failure to converge
otherDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), other2016, start = list(a1 = -253, a2 = -196, b1 = 0.007, b2 = 0.94, b2p = -0.112), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p not significant
otherDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), other2016, weights = TotalHt^-2)
otherDiameterFromHeightMichaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = list(a1 = 139, a2 = 78, b1 = 0.92), weights = TotalHt^-2) # a1p, a2p, b1p not significant, >50 iterations with nlrob()
otherDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 0.37, a1p = -0.27, a2 = -0.14, a2p = -0.038), weights = TotalHt^-2) # converges poorly, yields negative values
otherDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2), other2016, weights = TotalHt^-2) # isPlantation*(TotalHt - 1.37)^2 not significant
otherDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 1.36, a1p = -0.48, b1 = 1.15, b1p = 0.11), weights = pmin(TotalHt^-2, 0.5))
otherDiameterFromHeightPowerAal = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.0, a1p = -0.38, a2 = -0.013, b1 = 1.06), weights = pmin(TotalHt^-2, 0.5)) # a2p, b2p not significant, nlrob() failure to converge
otherDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = list(a1 = 1.81, a2 = -0.0007, a3 = -1.05, a4 = 0.034, a5 = -0.024, a6 = 0.00049, b1 = 1.21), weights = pmin(TotalHt^-2, 0.5)) # a1p, a4, a5, a6 not significant
otherDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 1.32, a2 = 0.40, b1 = 1.11, b1p = -0.055), weights = pmin(TotalHt^-2, 0.5))
otherDiameterFromHeightRuark = gsl_nls(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 1.56, a1p = 0.90, b1 = 1.05, b1p = -0.57, b2 = 0.0057, b2p = 0.053), weights = pmin(TotalHt^-2, 0.5)) # step factor with nlrob()
otherDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.20, Ha = 187), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with nlrob()
otherDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, other2016, start = list(a1 = 8500, a2 = -1.2, b1 = 0.02, b2 = -0.05, b3 = 2.0), weights = pmin(TotalHt^-2, 0.5), control = nls.control(maxiter = 50)) # step size or singular gradient with nls(), NaN-inf with nlrob()
otherDiameterFromHeightSibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.67, a1p = 1.84, b1 = 1.70, b1p = -1.26, b2 = -0.063, b2p = 0.30), weights = pmin(TotalHt^-2, 0.5))
otherDiameterFromHeightSibbesenFormAal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.58, a1p = 1.93, a2 = -0.0005, b1 = 1.93, b1p = -1.49, b2 = -0.085, b2p = 0.32), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
otherDiameterFromHeightSibbesenFormPhysio = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = list(a1 = 4.4, a1p = -0.56, a2 = -0.0012, a3 = -1.82, a4 = 0.007, a5 = 0.09, a6 = 0.002, b1 = 0.48, b2 = 0.18, b2p = 0.028), weights = pmin(TotalHt^-2, 0.5)) # a4, a5, a6 not significant, nlrob() failure to converge
otherDiameterFromHeightSibbesenFormRelHt = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 3.0, a1p = -0.57, a2 = -0.41, b1 = 0.49, b2 = 0.19, b2p = 0.035), weights = pmin(TotalHt^-2, 0.5)) # a2p not significant, nlrob() failure to converge
otherDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -138, b1 = 0.0116, b2 = 1.00), weights = pmin(TotalHt^-2, 0.5)) # nlrob() > 20 steps
otherDiameterFromHeightWeibull = nlrob(DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -163, a1p = 21.6, b1 = 0.0110, b2 = 0.99), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with b1p
#confint2(otherDiameterFromHeightWeibull, level = 0.99)

otherDiameterFromHeightChapmanForm = get_dbh_error(otherDiameterFromHeightChapmanForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormAal = get_dbh_error(otherDiameterFromHeightChapmanFormAal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBal = get_dbh_error(otherDiameterFromHeightChapmanFormBal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(otherDiameterFromHeightChapmanFormBalRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanFormRelHt = get_dbh_error(otherDiameterFromHeightChapmanFormRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichards = get_dbh_error(otherDiameterFromHeightChapmanRichards, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsAal = get_dbh_error(otherDiameterFromHeightChapmanRichardsAal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(otherDiameterFromHeightChapmanRichardsPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(otherDiameterFromHeightChapmanRichardsRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightLinear = get_dbh_error(otherDiameterFromHeightLinear, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightMichaelisMentenForm = get_dbh_error(otherDiameterFromHeightMichaelisMentenForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightNaslund = get_dbh_error(otherDiameterFromHeightNaslund, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightParabolic = get_dbh_error(otherDiameterFromHeightParabolic, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPower = get_dbh_error(otherDiameterFromHeightPower, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerAal = get_dbh_error(otherDiameterFromHeightPowerAal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightPowerPhysio = get_dbh_error(otherDiameterFromHeightPowerPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightPowerRelHt = get_dbh_error(otherDiameterFromHeightPowerRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightRuark = get_dbh_error(otherDiameterFromHeightRuark, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSchnute = get_dbh_error(otherDiameterFromHeightSchnute, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSharmaParton = get_dbh_error(otherDiameterFromHeightSharmaParton, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenForm = get_dbh_error(otherDiameterFromHeightSibbesenForm, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormAal = get_dbh_error(otherDiameterFromHeightSibbesenFormAal, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightSibbesenFormPhysio = get_dbh_error(otherDiameterFromHeightSibbesenFormPhysio, other2016physio, other2016natural, other2016plantationPhysio)
otherDiameterFromHeightSibbesenFormRelHt = get_dbh_error(otherDiameterFromHeightSibbesenFormRelHt, other2016, other2016natural, other2016plantation)
otherDiameterFromHeightWeibull = get_dbh_error(otherDiameterFromHeightWeibull, other2016, other2016natural, other2016plantation)

otherDiameterFromHeightResults = bind_rows(as_row("Chapman-Richards", otherDiameterFromHeightChapmanRichards),
                                           as_row("Chapman-Richards AAL", otherDiameterFromHeightChapmanRichardsAal),
                                           as_row("Chapman-Richards physio", otherDiameterFromHeightChapmanRichardsPhysio),
                                           as_row("Chapman-Richards RelHt", otherDiameterFromHeightChapmanRichardsRelHt),
                                           as_row("Chapman-Richards form", otherDiameterFromHeightChapmanForm),
                                           as_row("Chapman-Richards form AAL", otherDiameterFromHeightChapmanFormAal),
                                           as_row("Chapman-Richards form BAL", otherDiameterFromHeightChapmanFormBal),
                                           as_row("Chapman-Richards form BAL RelHt", otherDiameterFromHeightChapmanFormBalRelHt),
                                           as_row("Chapman-Richards form RelHt", otherDiameterFromHeightChapmanFormRelHt),
                                           as_row("linear", otherDiameterFromHeightLinear),
                                           as_row("generalized Michaelis-Menten form", otherDiameterFromHeightMichaelisMentenForm),
                                           as_row("Näslund", otherDiameterFromHeightNaslund),
                                           as_row("parabolic", otherDiameterFromHeightParabolic),
                                           as_row("power", otherDiameterFromHeightPower),
                                           as_row("power AAL", otherDiameterFromHeightPowerAal),
                                           as_row("power physio", otherDiameterFromHeightPowerPhysio),
                                           as_row("power RelHt", otherDiameterFromHeightPowerRelHt), # not AIC supported
                                           as_row("Ruark", otherDiameterFromHeightRuark),
                                           as_row("Schnute", otherDiameterFromHeightSchnute),
                                           as_row("modified Sharma-Parton", otherDiameterFromHeightSharmaParton),
                                           as_row("Sibbesen form", otherDiameterFromHeightSibbesenForm),
                                           as_row("Sibbesen form AAL", otherDiameterFromHeightSibbesenFormAal),
                                           as_row("Sibbesen form physio", otherDiameterFromHeightSibbesenFormPhysio),
                                           as_row("Sibbesen form RelHt", otherDiameterFromHeightSibbesenFormRelHt),
                                           as_row("Weibull", otherDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "height", species = "other", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(otherDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(other2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = otherDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = otherDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = otherDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "generalized Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  geom_line(aes(x = otherDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation)) +
  geom_line(aes(x = otherDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  #geom_line(aes(x = otherDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form aBAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 90, label = "other species, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
otherParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", otherHeightFromDiameterChapmanRichards$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards BAL", otherHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards BAL physio", otherHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards BAL RelHt", otherHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards physio", otherHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                      bind_rows(c(method = "Curtis", otherHeightFromDiameterCurtis$m$getPars())),
                                      bind_rows(c(method = "Hossfeld", otherHeightFromDiameterHossfeld$m$getPars())),
                                      bind_rows(c(method = "Korf", otherHeightFromDiameterKorf$m$getPars())),
                                      bind_rows(c(method = "linear", otherHeightFromDiameterLinear$coefficients)),
                                      bind_rows(c(method = "generalized Michaelis-Menten", otherHeightFromDiameterMichaelisMenten$m$getPars())),
                                      bind_rows(c(method = "parabolic", otherHeightFromDiameterParabolic$coefficients)),
                                      bind_rows(c(method = "power", otherHeightFromDiameterPower$m$getPars())),
                                      bind_rows(c(method = "Prodan", otherHeightFromDiameterProdan$m$getPars())),
                                      bind_rows(c(method = "Ratkowsky", otherHeightFromDiameterRatkowsky$m$getPars())),
                                      bind_rows(c(method = "Richards", otherHeightFromDiameterRichards$m$getPars())),
                                      bind_rows(c(method = "Sharma-Parton", otherHeightFromDiameterSharmaParton$m$getPars())),
                                      bind_rows(c(method = "Sharma-Parton BAL", otherHeightFromDiameterSharmaPartonBal$m$getPars())),
                                      bind_rows(c(method = "Sharma-Parton BAL physio", otherHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                      bind_rows(c(method = "Sharma-Parton physio", otherHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                      bind_rows(c(method = "Sharma-Zhang", otherHeightFromDiameterSharmaZhang$m$getPars())),
                                      bind_rows(c(method = "Sharma-Zhang BAL", otherHeightFromDiameterSharmaZhangBal$m$getPars())),
                                      bind_rows(c(method = "Sibbesen", otherHeightFromDiameterSibbesen$m$getPars())),
                                      bind_rows(c(method = "Weibull", otherHeightFromDiameterWeibull$m$getPars())),
                                      bind_rows(c(method = "Weibull BAL", otherHeightFromDiameterWeibullBal$m$getPars())),
                                      bind_rows(c(method = "Weibull RelHt", otherHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards GNLS", otherHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                      bind_rows(c(method = "Chapman-Richards BAL GNLS", otherHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                      bind_rows(c(method = "Sharma-Parton GNLS", otherHeightFromDiameterSharmaPartonGnls$coefficients)),
                                      bind_rows(c(method = "Sharma-Parton BAL GNLS", otherHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                      bind_rows(c(method = "Sharma-Zhang GNLS", otherHeightFromDiameterSharmaZhangGnls$coefficients)),
                                      bind_rows(c(method = "Sharma-Zhang BAL GNLS", otherHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                      bind_rows(c(method = "Weibull GNLS", otherHeightFromDiameterWeibullGnls$coefficients)),
                                      bind_rows(c(method = "Weibull BAL GNLS", otherHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                              mutate(responseVariable = "DBH"),
                            bind_rows(bind_rows(c(method = "Chapman-Richards", otherDiameterFromHeightChapmanRichards$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards AAL", otherDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards physio", otherDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                      bind_rows(c(method = "Chapman-Richards RelHt", otherDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                      #bind_rows(c(method = "Chapman-Richards form", otherDiameterFromHeightChapmanForm$m$getPars())),
                                      #bind_rows(c(method = "Chapman-Richards form AAL", otherDiameterFromHeightChapmanFormAal$m$getPars())),
                                      #bind_rows(c(method = "Chapman-Richards form BAL", otherDiameterFromHeightChapmanFormBal$m$getPars())),
                                      #bind_rows(c(method = "Chapman-Richards form BAL RelHt", otherDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                      #bind_rows(c(method = "Chapman-Richards form RelHt", otherDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                      bind_rows(c(method = "linear", otherDiameterFromHeightLinear$coefficients)),
                                      bind_rows(c(method = "generalized Michaelis-Menten form", otherDiameterFromHeightMichaelisMentenForm$m$getPars())),
                                      bind_rows(c(method = "Näslund", otherDiameterFromHeightNaslund$m$getPars())),
                                      bind_rows(c(method = "parabolic", otherDiameterFromHeightParabolic$coefficients)),
                                      bind_rows(c(method = "power", otherDiameterFromHeightPower$m$getPars())),
                                      bind_rows(c(method = "power AAL", otherDiameterFromHeightPowerAal$m$getPars())),
                                      bind_rows(c(method = "power physio", otherDiameterFromHeightPowerPhysio$m$getPars())),
                                      bind_rows(c(method = "power RelHt", otherDiameterFromHeightPowerRelHt$m$getPars())),
                                      bind_rows(c(method = "Ruark", otherDiameterFromHeightRuark$m$getPars())),
                                      bind_rows(c(method = "Schnute", otherDiameterFromHeightSchnute$m$getPars())),
                                      #bind_rows(c(method = "modified Sharma-Parton", otherDiameterFromHeightSharmaParton$m$getPars())),
                                      bind_rows(c(method = "Sibbesen form", otherDiameterFromHeightSibbesenForm$m$getPars())),
                                      bind_rows(c(method = "Sibbesen form AAL", otherDiameterFromHeightSibbesenFormAal$m$getPars())),
                                      bind_rows(c(method = "Sibbesen form physio", otherDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                      bind_rows(c(method = "Sibbesen form RelHt", otherDiameterFromHeightSibbesenFormRelHt$m$getPars())),
                                      bind_rows(c(method = "Weibull", otherDiameterFromHeightWeibull$m$getPars()))) %>%
                              mutate(responseVariable = "height")) %>%
  mutate(species = "other",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

