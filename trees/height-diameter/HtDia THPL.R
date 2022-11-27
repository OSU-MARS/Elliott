# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western redcedar height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Parton, Sharma-Zhang, Chapman-Richards BAL
#thplHeightFromDiameterSharmaPartonBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b1p = -0.0143, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = pmin(DBH^-1.2, 1))
thpl2016 = trees2016 %>% filter(Species == "RC", isLiveUnbroken, TotalHt > 0) # live western redcedars measured for height
thpl2016natural = thpl2016 %>% filter(isPlantation == FALSE)
thpl2016physio = thpl2016 %>% filter(is.na(elevation) == FALSE)
thpl2016plantation = thpl2016 %>% filter(isPlantation)
thpl2016plantationPhysio = thpl2016physio %>% filter(isPlantation)

thplHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131), weights = pmin(DBH^-1.2, 1)) # a1p, b1p, b2p not significant
thplHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 51.6, a1p = -9.09, a2 = -0.141, a2p = 0.464, a3 = -0.026, b1 = -0.016, b2 = 1.139), weights = pmin(DBH^-1.2, 1)) # a3p, b1p, b2p not significant
thplHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 42.7, a2 = -0.159, a2p = 0.530, a3 = 0.001, a4 = 0.195, a5 = 0.113, a6 = -0.532, a7 = 0, b1 = -0.015, b1p = 0.003, b2 = 1.101), weights = pmin(DBH^-1.2, 1)) # a2, a3, a4, a5, a6, b2p not significant, a1p and b1p not mutually significant
thplHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 8.7, a1p = -0.810, a2 = 21.1, a2p = 0.239, a3 = -0.043, a4 = 64.0, a4p = -42.8, b1 = -0.021, b2 = 0.120, b2p = 0.916), weights = pmin(DBH^-1.2, 1)) # a2, a3, a3p, b1p not significant
thplHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 45.8, a1p = -16.8, a2 = -0.0004, a3 = 6.265, a4 = 0.187, a5 = 0.393, a6 = 0.229, b1 = -0.013, b1p = -0.008, b2 = 1.148), weights = pmin(DBH^-1.2, 1)) # a2, a3, a4, a5, b2p not significant
thplHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016, start = list(a1 = 0.560, b1 = 0.069), weights = pmin(DBH^-1.2, 1)) # a1p, b1p not significant
thplHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176), weights = pmin(DBH^-1.2, 1)) # b2p not significant
thplHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), thpl2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175), weights = pmin(DBH^-1.2, 1)) # a1p, b1p, b2p not significant
thplHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH, thpl2016, offset = breastHeight, weights = pmin(DBH^-1.2, 1)) # isPlantation*DBH not significant (p = 0.044)
thplHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176), weights = pmin(DBH^-1.2, 1)) # b1p not significant
thplHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2), thpl2016, offset = breastHeight, weights = pmin(DBH^-1.2, 1)) # isPlantation*DBH not quite significant (p = 0.106), isPlantation*DBH^2 not significant
thplHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649), weights = pmin(DBH^-1.2, 1)) # a2p, a3p not significant
thplHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + a1*DBH^b1, thpl2016, start = list(a1 = 0.542, b1 = 0.939), weights = pmin(DBH^-1.2, 1)) # a1p, b1p not significant
thplHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151), weights = pmin(DBH^-1.2, 1))
thplHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), thpl2016, start = list(Ha = 55.3, Hap = -30.1, d = 0.546, dp = 0.298, kU = 0.009, kUp = 0.017), weights = pmin(DBH^-1.2, 1))
thplHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, thpl2016, start = list(a1 = 38.0, a2 = 0.131, a2p = -0.135, b1 = -0.015, b1p = -0.011, b2 = -0.114, b3 = 1.09), weights = pmin(DBH^-1.2, 1)) # a1p, b2p, b3p not significant
thplHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, thpl2016, start = list(a1 = 50.6, a1p = -15.8, a2 = 0.023, b1 = -0.014, b1p = -0.009, b2 = -0.069, b3 = 1.130), weights = pmin(DBH^-1.2, 1)) # a2p, b2p, b3p not significant
thplHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, thpl2016physio, start = list(a1 = 41.6, a1p = -13.2, a2 = 0.047, a3 = -0.00001, a4 = 0.163, a5 = 0.012, a6 = 0.006, b1 = -0.014, b1p = -0.010, b2 = -0.056, b3 = 1.121), weights = pmin(DBH^-1.2, 1)) # a2, a2p, a3, a4, a5, b2p, b3p not significant
thplHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^b2*DBH))^b3, thpl2016physio, start = list(a1 = 39.5, a1p = -12.5, a2 = 0.055, a3 = -0.00002, a4 = 0.016, a5 = 0.013, a6 = 0.007, b1 = -0.014, b1p = -0.010, b2 = -0.053, b3 = 1.122), weights = pmin(DBH^-1.2, 1)) # a2, a2p, a3, a5, a5, a6, b2p, b3p not significant
thplHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 40.1, a1p = -4.259, a2 = 0.040, b1 = -0.042, b2 = -0.148, b3 = 1.190, b3p = -0.097), weights = pmin(DBH^-1.2, 1)) # a2, a2p, b1p, b2p not significant
thplHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, thpl2016, start = list(a1 = 53.2, a1p = -8.857, a2 = -0.016, a3 = -0.002, a3p = 0.010, b1 = -0.025, b2 = -0.078, b3 = 1.126), weights = pmin(DBH^-1.2, 1)) # a2, a2p, a3, ab1, b2p, b3p not significant
thplHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), thpl2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078), weights = pmin(DBH^-1.2, 1)) # a1p, b1p, b2p not significant
thplHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141), weights = pmin(DBH^-1.2, 1)) # b2p not significant
thplHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131), weights = pmin(DBH^-1.2, 1)) # a1p, a2, a3, b1p, b2p not significant
thplHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a3 = -0.102, a3p = 0.068, a4 = 46.6, a4p = -9.98, b1 = -0.019, b2 = 0.778), weights = pmin(DBH^-1.2, 1)) # a1p, a2, a3, b1p, b2p not significant
#confint2(thplHeightFromDiameterWeibullBalRelHt, level = 0.99)

thplHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", thplHeightFromDiameterChapmanRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BAL", thplHeightFromDiameterChapmanRichardsBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BAL physio", thplHeightFromDiameterChapmanRichardsBalPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BAL RelHt", thplHeightFromDiameterChapmanRichardsBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", thplHeightFromDiameterChapmanRichardsPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterCurtis = get_height_error("Curtis", thplHeightFromDiameterCurtis, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", thplHeightFromDiameterHossfeld, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterKorf = get_height_error("Korf", thplHeightFromDiameterKorf, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterLinear = get_height_error("linear", thplHeightFromDiameterLinear, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", thplHeightFromDiameterMichaelisMenten, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterParabolic = get_height_error("parabolic", thplHeightFromDiameterParabolic, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterProdan = get_height_error("Prodan", thplHeightFromDiameterProdan, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterPower = get_height_error("power", thplHeightFromDiameterPower, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", thplHeightFromDiameterRatkowsky, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterRichards = get_height_error("unified Richards", thplHeightFromDiameterRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", thplHeightFromDiameterSharmaParton, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BAL", thplHeightFromDiameterSharmaPartonBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BAL physio", thplHeightFromDiameterSharmaPartonBalPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", thplHeightFromDiameterSharmaPartonPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", thplHeightFromDiameterSharmaZhang, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BAL", thplHeightFromDiameterSharmaZhangBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSibbesen = get_height_error("Sibbesen", thplHeightFromDiameterSibbesen, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibull = get_height_error("Weibull", thplHeightFromDiameterWeibull, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBal = get_height_error("Weibull BAL", thplHeightFromDiameterWeibullBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BAL RelHt", thplHeightFromDiameterWeibullBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)

thplHeightFromDiameterResults = bind_rows(as_row(thplHeightFromDiameterChapmanRichards),
                                          as_row(thplHeightFromDiameterChapmanRichardsBal),
                                          as_row(thplHeightFromDiameterChapmanRichardsBalPhysio, significant = FALSE),
                                          as_row(thplHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(thplHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(thplHeightFromDiameterCurtis),
                                          as_row(thplHeightFromDiameterHossfeld),
                                          as_row(thplHeightFromDiameterKorf),
                                          as_row(thplHeightFromDiameterLinear),
                                          as_row(thplHeightFromDiameterMichaelisMenten),
                                          as_row(thplHeightFromDiameterParabolic),
                                          as_row(thplHeightFromDiameterPower),
                                          as_row(thplHeightFromDiameterProdan),
                                          as_row(thplHeightFromDiameterRatkowsky),
                                          as_row(thplHeightFromDiameterRichards),
                                          as_row(thplHeightFromDiameterSharmaParton),
                                          as_row(thplHeightFromDiameterSharmaPartonBal),
                                          as_row(thplHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(thplHeightFromDiameterSharmaPartonPhysio, significant = FALSE),
                                          as_row(thplHeightFromDiameterSharmaZhang),
                                          as_row(thplHeightFromDiameterSharmaZhangBal),
                                          as_row(thplHeightFromDiameterSibbesen),
                                          as_row(thplHeightFromDiameterWeibull),
                                          as_row(thplHeightFromDiameterWeibullBal),
                                          as_row(thplHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "THPL", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(thplHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = thpl2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = thpl2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterKorf$fitted.values, color = "Korf", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterLinear$fitted.values, color = "linear", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterPower$fitted.values, color = "power", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterProdan$fitted.values, color = "Prodan", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = thpl2016$isPlantation)) +
  annotate("text", x = 0, y = 65, label = "western redcedar, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 65)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## western redcedar height-diameter GNLS regressions
#thplHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, thpl2016, start = thplHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = TRUE)) # step halving with varPower(0.60)
#thplHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, thpl2016, start = thplHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameterSharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameterSharmaPartonBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), thpl2016, start = thplHeightFromDiameterSharmaZhang$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameterSharmaZhangBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = thplHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = thplHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#save(thplHeightFromDiameterChapmanRichardsGnls, thplHeightFromDiameterChapmanRichardsBalGnls, thplHeightFromDiameterSharmaPartonGnls, thplHeightFromDiameterSharmaPartonBalGnls, thplHeightFromDiameterSharmaZhangGnls, thplHeightFromDiameterSharmaZhangBalGnls, thplHeightFromDiameterWeibullGnls, thplHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia THPL GNLS.rdata")

load("trees/height-diameter/HtDia THPL GNLS.rdata")
thplHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", thplHeightFromDiameterChapmanRichardsGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BAL GNLS", thplHeightFromDiameterChapmanRichardsBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", thplHeightFromDiameterSharmaPartonGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BAL GNLS", thplHeightFromDiameterSharmaPartonBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", thplHeightFromDiameterSharmaZhangGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BAL GNLS", thplHeightFromDiameterSharmaZhangBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", thplHeightFromDiameterWeibullGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BAL GNLS", thplHeightFromDiameterWeibullBalGnls, thpl2016, thpl2016natural, thpl2016plantation)

thplHeightFromDiameterResultsGnls = bind_rows(as_row(thplHeightFromDiameterChapmanRichardsGnls),
                                              as_row(thplHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(thplHeightFromDiameterSharmaPartonGnls),
                                              as_row(thplHeightFromDiameterSharmaPartonGnls),
                                              as_row(thplHeightFromDiameterSharmaZhangGnls),
                                              as_row(thplHeightFromDiameterSharmaZhangBalGnls),
                                              as_row(thplHeightFromDiameterWeibullGnls),
                                              as_row(thplHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "THPL", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
thplHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(thplHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(thplHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(thplHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
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
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterWeibullBAL$fitted.values, color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = thpl2016natural$DBH, y = thplHeightFromDiameterWeibullBALnatural$fitted.values, color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = thpl2016plantation$DBH, y = thplHeightFromDiameterWeibullBALplantation$fitted.values, color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterBase$fitted.values, color = "base")) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterWeibull$fitted.values, color = "ElliottWeibull")) +
  annotate("text", x = 0, y = 85, label = "a) western redcedar, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

## western redcedar diameter-height regressions
#thplDiameterFromHeightChapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, iter = 1000,
#                                                  start_lower = list(a1 = -10, b1 = -1, b2 = -1), 
#                                                  start_upper = list(a1 = 50, b1 = 1, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
#thplDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a4 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 115, a2 = 0, a4 = 13.3, b1 = 0.020, b2 = 0.92), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = nls.control(maxiter = 500)) # step size with nls(), step factor with BAL+nlrob()
#thplDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 139, a2 = 10, a3 = -0.066, a4 = 11.6, b1 = 0.016, b2 = 0.944), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = gsl_nls_control(maxiter = 250)) # step factor with BAL+nlrob()
#thplDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, iter = 100,
#                                                      start_lower = list(a1 = 1, b1 = -0.5, b2 = 0.1), 
#                                                      start_upper = list(a1 = 300, b1 = 0.5, b2 = 2.5), modelweights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
#thplDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), thpl2016, start = list(a1 = 229.6, b1 = -0.0057, b2 = 1.25, b2p = 0), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
#thplDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, iter = 100,
#                                                            start_lower = list(a1 = -200, a1p = -100, a2 = -0.1, a3 = -1, a4 = -1, a5 = -1, a6 = -0.1, b1 = -0.1, b1p = -1, b2 = -1), 
#                                                            start_upper = list(a1 = 1, a1p = 100, a2 = 0.1, a3 = 1, a4 = 1, a5 = 1, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
#thplDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, thpl2016, iter = 100,
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = -0.5, b2 = -0.5, b3 = -1), 
#                                                   start_upper = list(a1 = 10, a2 = 1.5, b1 = 0.5, b2 = 0.5, b3 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
#thplDiameterFromHeightSchnute = nlrob(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016, start = list(a1 = 0.0003, a2 = 0.1, b1 = 1.5, Ha = 75), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # NaN-inf with nlrob()
thplDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 960, b1 = 0.003, b2 = 1.05), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p, b1p, b2p not significant, singular gradient with nls(), no convergence from nls_multstart()
thplDiameterFromHeightChapmanFormAat = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 960, a2 = -100, b1 = 0.001, b2 = 1.07), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = nls.control(maxiter = 50)) # NaN-inf with nls() and nlrob
thplDiameterFromHeightChapmanFormBal = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 913, a2 = -40, b1 = 0.0003, b2 = 1.04), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = nls.control(maxiter = 50)) # step size with nls() and nlrob()
thplDiameterFromHeightChapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 655, a2 = 0, a3 = 0, a4 = 2.3, b1 = 0.003, b2 = 1.04), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = gsl_nls_control(maxiter = 50)) # nlrob() step factor with either a2 or a3
thplDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 655, a2 = 2.3, b1 = 0.003, b2 = 1.04), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = nls.control(maxiter = 50)) # step size with nls(), >50 iterations with nlrob()
thplDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -763, b1 = 0.0027, b2 = 1.05), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p and b2p not significant, poor convergence with b1p
thplDiameterFromHeightChapmanRichardsAat = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -1200, a2 = 0, b1 = 0.0016, b2 = 1.06), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = nls.control(maxiter = 50)) # a1p, b1p not significant, step factor with nlrob()
thplDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016physio, start = list(a1 = -500, a1p = 200, a2 = 0.07, a3 = 9.8, a4 = -7.2, a5 = -1.3, a6 = 0.7, b1 = 0.0047, b1p = 0.007, b2 = 1.005), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = list(maxiter = 50)) # no physiographic effects significant
thplDiameterFromHeightChapmanRichardsRelHt = gsl_nls(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -350, a2 = -66, b1 = 0.006, b2 = 0.98), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = list(maxiter = 50)) # step factor with nlrob()
thplDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37), thpl2016, weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # isPlantation*(TotalHt - 1.37) not significant
thplDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), thpl2016, start = list(a1 = 519, a2 = 237, b1 = 1.00), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p, a2p, b1p not significant
thplDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
thplDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), thpl2016, weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # (TotalHt - 1.37)^2 not significant
thplDiameterFromHeightPower = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no significant plantation effects
thplDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no significant plantation effects
thplDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, start = list(a1 = 2.26, a2 = -0.00047, a3 = -0.156, a4 = -0.0018, a5 = 0.0208, a6 = -0.0060, b1 = 1.08), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no significant physiographic effects
thplDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.68, a2 = -0.11, a2p = 0.23, b1 = 1.13), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p and b1p not significant
thplDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1 * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.12, b1 = 1.01, b2 = 0.0038, b2p = 0.0025), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p, b1p not significant
thplDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016, start = list(a1 = 0.00005, a2 = 0.0074, b1 = 1.12, Ha = 52), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # singular gradient with nlrob()
thplDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, thpl2016, start = list(a1 = 65, a2 = 0.4, b1 = 0.006, b2 = 0.03, b3 = 0.7), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf with nls() from nls_multstart() point, NaN-inf, singular gradient, or code syntax error with nlrob()
thplDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.42, b1 = 1.29, b2 = -0.27), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no significant plantation effects
thplDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no significant plantation effects
thplDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), thpl2016physio, start = list(a1 = 2.32, a1p = -0.294, a2 = -0.00069, a3 = -0.116, a4 = 0.0374, a5 = 0.0047, a6 = -0.0056, b1 = 1.084, b2 = -0.0048, b2p = 0.0203), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # no physiographic effects significant
thplDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.34, a2 = 0.30, b1 = 1.32, b2 = -0.040), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
thplDiameterFromHeightWeibull = nlrob(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016, start = list(a1 = -330, b1 = 0.0065, b2 = 1.01), weights = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5)) # a1p, b1p, b2p not significant
#confint_nlrob(thplDiameterFromHeightChapmanRichardsPhysio, level = 0.99, weights = pmin(thpl2016$TotalHt^if_else(thpl2016$isPlantation, -1.7, -1.6), 0.5))

thplDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", thplDiameterFromHeightChapmanForm, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AAT", thplDiameterFromHeightChapmanFormAat, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BAL", thplDiameterFromHeightChapmanFormBal, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BAL RelHt", thplDiameterFromHeightChapmanFormBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", thplDiameterFromHeightChapmanFormRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", thplDiameterFromHeightChapmanRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AAT", thplDiameterFromHeightChapmanRichardsAat, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", thplDiameterFromHeightChapmanRichardsPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", thplDiameterFromHeightChapmanRichardsRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightLinear = get_dbh_error("linear", thplDiameterFromHeightLinear, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", thplDiameterFromHeightMichaelisMentenForm, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightNaslund = get_dbh_error("Näslund", thplDiameterFromHeightNaslund, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightParabolic = get_dbh_error("parabolic", thplDiameterFromHeightParabolic, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPower = get_dbh_error("power", thplDiameterFromHeightPower, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerAat = get_dbh_error("power AAT", thplDiameterFromHeightPowerAat, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerPhysio = get_dbh_error("power physio", thplDiameterFromHeightPowerPhysio, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", thplDiameterFromHeightPowerRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightRuark = get_dbh_error("Ruark", thplDiameterFromHeightRuark, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSchnute = get_dbh_error("Schnute", thplDiameterFromHeightSchnute, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", thplDiameterFromHeightSharmaParton, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", thplDiameterFromHeightSibbesenForm, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AAT", thplDiameterFromHeightSibbesenFormAat, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", thplDiameterFromHeightSibbesenFormPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", thplDiameterFromHeightSibbesenFormRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightWeibull = get_dbh_error("Weibull", thplDiameterFromHeightWeibull, thpl2016, thpl2016natural, thpl2016plantation)

thplDiameterFromHeightResults = bind_rows(as_row(thplDiameterFromHeightChapmanRichards),
                                          as_row(thplDiameterFromHeightChapmanRichardsAat, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanRichardsPhysio, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanRichardsRelHt, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanForm),
                                          as_row(thplDiameterFromHeightChapmanFormAat, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanFormBal, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanFormBalRelHt, significant = FALSE),
                                          as_row(thplDiameterFromHeightChapmanFormRelHt),
                                          as_row(thplDiameterFromHeightLinear),
                                          as_row(thplDiameterFromHeightMichaelisMentenForm),
                                          as_row(thplDiameterFromHeightNaslund),
                                          as_row(thplDiameterFromHeightParabolic),
                                          as_row(thplDiameterFromHeightPower),
                                          as_row(thplDiameterFromHeightPowerAat),
                                          as_row(thplDiameterFromHeightPowerPhysio),
                                          as_row(thplDiameterFromHeightPowerRelHt),
                                          as_row(thplDiameterFromHeightRuark),
                                          as_row(thplDiameterFromHeightSharmaParton),
                                          as_row(thplDiameterFromHeightSchnute),
                                          as_row(thplDiameterFromHeightSibbesenForm),
                                          as_row(thplDiameterFromHeightSibbesenFormAat),
                                          as_row(thplDiameterFromHeightSibbesenFormPhysio),
                                          as_row(thplDiameterFromHeightSibbesenFormRelHt),
                                          as_row(thplDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "DBH", species = "THPL", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(thplDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(thpl2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = thplDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -100 * log(1 - pmin(0.015*(TotalHt - 1.37)^1.0, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
  #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.45, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAT", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = -1/0.0003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute"), alpha = 0.5) +
  geom_line(aes(x = 30*topHeight^0.5*(exp(0.01 * (tph/standBasalAreaPerHectare)^0.25*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "modified Sharma-Parton"), alpha = 0.5) +
  annotate("text", x = 0, y = 62, label = "western redcedar, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
thplParameters = bind_rows(bind_rows(get_coefficients(thplHeightFromDiameterChapmanRichards),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(thplHeightFromDiameterCurtis),
                                     get_coefficients(thplHeightFromDiameterHossfeld),
                                     get_coefficients(thplHeightFromDiameterKorf),
                                     get_coefficients(thplHeightFromDiameterLinear),
                                     get_coefficients(thplHeightFromDiameterMichaelisMenten),
                                     get_coefficients(thplHeightFromDiameterParabolic),
                                     get_coefficients(thplHeightFromDiameterPower),
                                     get_coefficients(thplHeightFromDiameterProdan),
                                     get_coefficients(thplHeightFromDiameterRatkowsky),
                                     get_coefficients(thplHeightFromDiameterRichards),
                                     get_coefficients(thplHeightFromDiameterSharmaParton),
                                     get_coefficients(thplHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(thplHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(thplHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(thplHeightFromDiameterSharmaZhang),
                                     get_coefficients(thplHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(thplHeightFromDiameterSibbesen),
                                     get_coefficients(thplHeightFromDiameterWeibull),
                                     get_coefficients(thplHeightFromDiameterWeibullBal),
                                     get_coefficients(thplHeightFromDiameterWeibullBalRelHt),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(thplHeightFromDiameterChapmanRichardsBalGnls),
                                     get_coefficients(thplHeightFromDiameterSharmaPartonGnls),
                                     get_coefficients(thplHeightFromDiameterSharmaPartonBalGnls),
                                     get_coefficients(thplHeightFromDiameterSharmaZhangGnls),
                                     get_coefficients(thplHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(thplHeightFromDiameterWeibullGnls),
                                     get_coefficients(thplHeightFromDiameterWeibullBalGnls)) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(thplDiameterFromHeightChapmanRichards),
                                     get_coefficients(thplDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(thplDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(thplDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(thplDiameterFromHeightChapmanForm),
                                     get_coefficients(thplDiameterFromHeightChapmanFormAat),
                                     get_coefficients(thplDiameterFromHeightChapmanFormBal),
                                     get_coefficients(thplDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(thplDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(thplDiameterFromHeightLinear),
                                     get_coefficients(thplDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(thplDiameterFromHeightNaslund),
                                     get_coefficients(thplDiameterFromHeightParabolic),
                                     get_coefficients(thplDiameterFromHeightPower),
                                     get_coefficients(thplDiameterFromHeightPowerAat),
                                     get_coefficients(thplDiameterFromHeightPowerPhysio),
                                     get_coefficients(thplDiameterFromHeightPowerRelHt),
                                     get_coefficients(thplDiameterFromHeightRuark),
                                     get_coefficients(thplDiameterFromHeightSchnute),
                                     get_coefficients(thplDiameterFromHeightSharmaParton),
                                     get_coefficients(thplDiameterFromHeightSibbesenForm),
                                     get_coefficients(thplDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(thplDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(thplDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(thplDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "THPL",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
thplBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), thpl2016, start = list(a1 = 90, b1 = 0.000003, b2 = 2.18), weights = pmin(1/basalArea, 1E4)) # a1p, b1p, b2p not significant, step factor with nlrob()
thplBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), thpl2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.00002, b1 = 2.14, b1p = 0.34), weights = pmin(1/basalArea, 1E4))
#confint2(thplBasalAreaFromHeightKorf, level = 0.99)
#confint_nlrob(thplBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/thpl2016$basalArea, 1E4))

thplBasalAreaFromHeightKorf$fitted.values = predict(thplBasalAreaFromHeightKorf, thpl2016)
thplBasalAreaFromHeightKorf$residuals = thplBasalAreaFromHeightKorf$fitted.values - thpl2016$basalArea
thplBasalAreaFromHeightPower$fitted.values = predict(thplBasalAreaFromHeightPower, thpl2016)
thplBasalAreaFromHeightPower$residuals = thplBasalAreaFromHeightPower$fitted.values - thpl2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(thplBasalAreaFromHeightKorf), 100^2 * mean(thplBasalAreaFromHeightKorf$residuals), mean(abs(thplBasalAreaFromHeightKorf$residuals)), 1 - sum(thplBasalAreaFromHeightKorf$residuals^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2)),
        "power", AIC(thplBasalAreaFromHeightPower), 100^2 * mean(thplBasalAreaFromHeightPower$residuals), mean(abs(thplBasalAreaFromHeightPower$residuals)), 1 - sum(thplBasalAreaFromHeightPower$residuals^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(thpl2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = thplBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = thplBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "western redcedar height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
