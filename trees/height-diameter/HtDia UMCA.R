# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Oregon myrtle (California bay) height-diameter regression form sweep
# preferred forms: Sibbesen, Korf, Prodan, Ratkowsky
#umcaHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), umca2016, start = list(Ha = 14.4, Hap = -3.6, d = 2.78, dp = -0.99, kU = 0.054, kUp = 0.026), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterChapmanRichards = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 17.7, a1p = -1.52, b1 = -0.0572, b2 = 1.23, b2p = -0.25), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 17.6, a1p = -6.97, a2 = 0.042, a2p = -0.198, a3 = -0.0039, a3p = 0.178, b1 = -0.049, b1p = -0.0743, b2 = 1.16, b2p = 0.156), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterChapmanRichardsBalPhysio = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(pi/180 * aspect) + a6 * cos(pi/180 * aspect)) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 69.4, a2 = 0.077, a2p = 0.72, a3 = -0.008, a4 = -0.031, a5 = 0.548, a6 = 0.980, b1 = -0.021, b1p = 0.0075, b2 = 1.49, b2p = -0.370), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterChapmanRichardsBalRelHt = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 17.6, a1p = -6.97, a2 = 0.042, a2p = -0.198, a3 = -0.0039, a3p = 0.178, a4 = 41.38, a4p = -16.4, b1 = -0.049, b2 = 1.16, b2p = 0.156), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterChapmanRichardsPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 68.5, a1p = -13.4, a2 = -0.0045, a3 = -8.09, a4 = 0.783, a5 = 0.766, a6 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^-2, 1)) # a4p not significant, a5p induces overfitting
umcaHeightFromDiameterCurtis = nls(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, umca2016, start = list(a1 = 1.6, b1 = 0.24), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterHossfeld = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterKorf = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), umca2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.0, a1p = -4.56, a2 = 44.6, a2p = -21.2, b1 = 1.30), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), umca2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterProdan = nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), umca2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterPower = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.9, a1p = -1.6, b1 = 0.63, b1p = 0.18), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterRatkowsky = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), umca2016, start = list(a1 = 22.9, a1p = -5.65, b1 = -17.68, b1p = 7.73, b2 = 3.93, b2p = -1.82), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016, start = list(Ha = 13.8, d = 2.46, kU = 0.0514, kUp = 0.0158), weights = pmin(DBH^-2, 1)) # Hap, dp not significant
umcaHeightFromDiameterSharmaParton = nls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSharmaPartonBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(pi/180 * aspect) + a5 * cos(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSharmaPartonPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(pi/180 * aspect) + a5 * cos(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 26, a1p = -16, a2 = -0.093, a2p = 0.21, b1 = -0.00469, b1p = -0.108, b2 = 0.387, b2p = -0.435, b3 = 1.228, b3p = -0.101), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSharmaZhangBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 37.9, a1p = -32.6, a2 = -0.203, a2p = 0.527, a3 = 0.0025, a3p = -0.0112, b1 = -0.00380, b1p = -0.0644, b2 = 0.408, b2p = -0.329, b3 = 1.200, b3p = 0.0300), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.0019, a1p = 0.163, b1 = 4.98, b1p = -2.72, b2 = -0.175, b2p = 0.0427), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterWeibull = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), umca2016, start = list(a1 = 18.0, a1p = -4.02, b1 = -0.0332, b1p = -0.0264, b2 = 1.11, b2p = -0.0031), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterWeibullBal = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), umca2016, start = list(a1 = 15.7, a2 = 0.090, a2p = 0.0356, a3 = -0.0302, a3p = -0.0689, b1 = -0.030, b1p = -0.0196, b2 = 1.184, b2p = -0.061), weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterWeibullBalRelHt = nls(TotalHt ~ 1.37 + (a1  + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation))), umca2016, start = list(a1 = -1.16, a2 = 0.0277, a2p = -0.141, a3 = -0.0060, a3p = 0.223, a4 = 51.6, a4p = -30.8, b1 = -0.658, b1p = 0.500, b2 = 0.729, b2p = 0.571), weights = pmin(DBH^-2, 1)) # a1p not significant
#confint2(umcaHeightFromDiameterWeibullBalRelHt, level = 0.99)

umca2016physio = umca2016 %>% filter(is.na(elevation) == FALSE)
umca2016plantationPhysio = umca2016physio %>% filter(isPlantation)
umcaHeightFromDiameterChapmanRichards = get_height_error(umcaHeightFromDiameterChapmanRichards, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBal = get_height_error(umcaHeightFromDiameterChapmanRichardsBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(umcaHeightFromDiameterChapmanRichardsBalPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(umcaHeightFromDiameterChapmanRichardsBalRelHt, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsPhysio = get_height_error(umcaHeightFromDiameterChapmanRichardsPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterCurtis = get_height_error(umcaHeightFromDiameterCurtis, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterHossfeld = get_height_error(umcaHeightFromDiameterHossfeld, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterKorf = get_height_error(umcaHeightFromDiameterKorf, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterLinear = get_height_error(umcaHeightFromDiameterLinear, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterMichaelisMenten = get_height_error(umcaHeightFromDiameterMichaelisMenten, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterParabolic = get_height_error(umcaHeightFromDiameterParabolic, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterProdan = get_height_error(umcaHeightFromDiameterProdan, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterPower = get_height_error(umcaHeightFromDiameterPower, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterRatkowsky = get_height_error(umcaHeightFromDiameterRatkowsky, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterRichards = get_height_error(umcaHeightFromDiameterRichards, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaParton = get_height_error(umcaHeightFromDiameterSharmaParton, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaPartonBal = get_height_error(umcaHeightFromDiameterSharmaPartonBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaPartonBalPhysio = get_height_error(umcaHeightFromDiameterSharmaPartonBalPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterSharmaPartonPhysio = get_height_error(umcaHeightFromDiameterSharmaPartonPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterSharmaZhang = get_height_error(umcaHeightFromDiameterSharmaZhang, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaZhangBal = get_height_error(umcaHeightFromDiameterSharmaZhangBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSibbesen = get_height_error(umcaHeightFromDiameterSibbesen, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibull = get_height_error(umcaHeightFromDiameterWeibull, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBal = get_height_error(umcaHeightFromDiameterWeibullBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBalRelHt = get_height_error(umcaHeightFromDiameterWeibullBalRelHt, umca2016, umca2016natural, umca2016plantation)

umcaHeightFromDiameterResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(umcaHeightFromDiameterChapmanRichards),
                                        "Chapman-Richards BAL", !!!as_row(umcaHeightFromDiameterChapmanRichardsBal),
                                        "Chapman-Richards BAL physio", !!!as_row(umcaHeightFromDiameterChapmanRichardsBalPhysio),
                                        "Chapman-Richards BAL RelHt", !!!as_row(umcaHeightFromDiameterChapmanRichardsBalRelHt),
                                        "Chapman-Richards physio", !!!as_row(umcaHeightFromDiameterChapmanRichardsPhysio),
                                        "Curtis", !!!as_row(umcaHeightFromDiameterCurtis),                                        "Hossfeld", !!!as_row(umcaHeightFromDiameterHossfeld),
                                        "Hossfeld", !!!as_row(umcaHeightFromDiameterHossfeld),
                                        "Korf", !!!as_row(umcaHeightFromDiameterKorf),
                                        "linear", !!!as_row(umcaHeightFromDiameterLinear),
                                        "Michaelis-Menten", !!!as_row(umcaHeightFromDiameterMichaelisMenten),
                                        "parabolic", !!!as_row(umcaHeightFromDiameterParabolic),
                                        "power", !!!as_row(umcaHeightFromDiameterPower),
                                        "Prodan", !!!as_row(umcaHeightFromDiameterProdan),
                                        "Ratkowsky", !!!as_row(umcaHeightFromDiameterRatkowsky),
                                        "unified Richards", !!!as_row(umcaHeightFromDiameterRichards),
                                        "Sharma-Parton", !!!as_row(umcaHeightFromDiameterSharmaParton),
                                        "Sharma-Parton BAL", !!!as_row(umcaHeightFromDiameterSharmaPartonBal),
                                        "Sharma-Parton BAL physio", !!!as_row(umcaHeightFromDiameterSharmaPartonBalPhysio),
                                        "Sharma-Parton physio", !!!as_row(umcaHeightFromDiameterSharmaPartonPhysio),
                                        "Sharma-Zhang", !!!as_row(umcaHeightFromDiameterSharmaZhang),
                                        "Sharma-Zhang BAL", !!!as_row(umcaHeightFromDiameterSharmaZhangBal),
                                        "Sibbesen", !!!as_row(umcaHeightFromDiameterSibbesen),
                                        "Weibull", !!!as_row(umcaHeightFromDiameterWeibull),
                                        "Weibull BAL", !!!as_row(umcaHeightFromDiameterWeibullBal),
                                        "Weibull BAL RelHt", !!!as_row(umcaHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "UMCA", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(umcaHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = umca2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = umca2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterKorf$fitted.values, color = "Korf", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterLinear$fitted.values, color = "linear", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterPower$fitted.values, color = "power", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterProdan$fitted.values, color = "Prodan", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = umca2016$isPlantation)) +
  annotate("text", x = 0, y = 35, label = "Oregon myrtle, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 35)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## Oregon myrtle height-diameter GNLS regressions
##umcaHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#umcaHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#umcaHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
##umcaHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
##umcaHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
##umcaHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#umcaHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), umca2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#umcaHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), umca2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(umcaHeightFromDiameterChapmanRichardsBalGnls, umcaHeightFromDiameterSharmaPartonGnls, umcaHeightFromDiameterSharmaZhangBalGnls, umcaHeightFromDiameterWeibullGnls, umcaHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia UMCA GNLS.rdata")
load("trees/height-diameter/HtDia UMCA GNLS.rdata")
umcaHeightFromDiameterWeibullGnls = umcaHeightFromDiameterWykoffGnls # temporary naming error fixup
umcaHeightFromDiameterWeibullBalGnls = umcaHeightFromDiameterWykoffBalGnls

#umcaHeightFromDiameterChapmanRichardsGnls = get_height_error(umcaHeightFromDiameterChapmanRichardsGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBalGnls = get_height_error(umcaHeightFromDiameterChapmanRichardsBalGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaPartonGnls = get_height_error(umcaHeightFromDiameterSharmaPartonGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaPartonBalGnls = get_height_error(umcaHeightFromDiameterSharmaPartonBalGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaZhangGnls = get_height_error(umcaHeightFromDiameterSharmaZhangGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaZhangBalGnls = get_height_error(umcaHeightFromDiameterSharmaZhangBalGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullGnls = get_height_error(umcaHeightFromDiameterWeibullGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBalGnls = get_height_error(umcaHeightFromDiameterWeibullBalGnls, umca2016, umca2016natural, umca2016plantation)

umcaHeightFromDiameterResultsGnls = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                            "Chapman-Richards GNLS", !!!as_row(NULL),
                                            "Chapman-Richards BAL GNLS", !!!as_row(umcaHeightFromDiameterChapmanRichardsBalGnls),
                                            "Sharma-Parton GNLS", !!!as_row(NULL),
                                            "Sharma-Parton BAL GNLS", !!!as_row(NULL),
                                            "Sharma-Zhang GNLS", !!!as_row(NULL),
                                            "Sharma-Zhang BAL GNLS", !!!as_row(NULL),
                                            "Weibull GNLS", !!!as_row(umcaHeightFromDiameterWeibullGnls),
                                            "Weibull BAL GNLS", !!!as_row(umcaHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "UMCA", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
umcaHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(umcaHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(umcaHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(umcaHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
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
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterWeibullBAL$fitted.values, color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = umca2016natural$DBH, y = umcaHeightFromDiameterWeibullBALnatural$fitted.values, color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = umca2016plantation$DBH, y = umcaHeightFromDiameterWeibullBALplantation$fitted.values, color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterBase$fitted.values, color = "base")) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterWeibull$fitted.values, color = "ElliottWeibull")) +
  annotate("text", x = 0, y = 85, label = "a) Oregon myrtle, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

## Oregon myrtle diameter-height regressions
#umcaDiameterFromHeightChapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, iter = 100,
#                                                  start_lower = list(a1 = -10, b1 = -2, b2 = -1), 
#                                                  start_upper = list(a1 = 100, b1 = 2.5, b2 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#umcaDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, umca2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.01, a2 = -1, b1 = -10, b2 = -0.5, b3 = 0.2),
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 10, b2 = 0.5, b3 = 1.5), modelweights = pmin(TotalHt^-2, 0.5))
#umcaDiameterFromHeightChapmanForm = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 2271, b1 = 0.001, b2 = 1.002), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf at multiple nls_multstart() points
#umcaDiameterFromHeightChapmanFormAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 15, a2 = -0, b1 = 0.01, b2 = 0.5), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf
#umcaDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 15, a2 = -0, b1 = 0.01, b2 = 0.5), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf
#umcaDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 13, a2 = -22.0, a3 = 10.4, a4 = -467, b1 = 0.0017, b2 = 0.971), weights = pmin(TotalHt^-2, 0.5)) # step factor
#umcaDiameterFromHeightChapmanFormRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 500, a2 = -21, b1 = 0.004, b2 = 1.0), weights = pmin(TotalHt^-2, 0.5)) # step factor
umcaDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), umca2016, start = list(a1 = 33.9, b1 = -0.047, b2 = 1.43, b2p = -0.15), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightChapmanRichardsAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), umca2016, start = list(a1 = 34.4, a2 = 0.004, b1 = -0.047, b2 = 1.41, b2p = -0.13), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightChapmanRichardsPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), umca2016, start = list(a1 = 58.6, a1p = 166.3, a2 = 0.0133, a3 = -52.7, a4 = 5.394, a5 = -0.837, a6 = 0.602, b1 = -0.0637, b1p = 0.0584, b2 = 1.329), weights = pmin(TotalHt^-2, 0.5)) # a2, a3, a4, a5 not significant
umcaDiameterFromHeightChapmanRichardsRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), umca2016, start = list(a1 = 39.0, a2 = -2.62, b1 = -0.043, b2 = 1.40, b2p = -0.132), weights = pmin(TotalHt^-2, 0.5)) # a1p convergence questionable, b1p not significant
umcaDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), umca2016, weights = TotalHt^-2)
umcaDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), umca2016, weights = TotalHt^-2)
umcaDiameterFromHeightPower = nls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightPowerAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p not significant
umcaDiameterFromHeightPowerPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 5.12, a1p = -0.67, a2 = 0.0013, a3 = -4.15, a4 = 0.46, a5 = -0.11, a6 = 0.044, b1 = 0.93), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2, b1p not significant
umcaDiameterFromHeightPowerRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 2.37, a1p = -0.887, a2 = -1.508, a2p = 1513, b1 = 1.123), weights = pmin(TotalHt^-2, 0.5)) 
#umcaDiameterFromHeightSharmaParton = nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, umca2016, start = list(a1 = 94.1, a2 = 0.65, b1 = 0.0001, b2 = -0.83, b3 = 0.32), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf at nls_multistart() point
#umcaDiameterFromHeightSharmaParton = nls(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 1, a2 = 1, a2p = 0, b1 = 0.02, b2 = 0.26, b2p = 0, b3 = 0.9, b3p = 0), weights = pmin(TotalHt^-2, 0.5)) # signular gradient
umcaDiameterFromHeightSibbesenForm = nls(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
umcaDiameterFromHeightSibbesenFormAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
umcaDiameterFromHeightSibbesenFormPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 3.688, a1p = -0.999, a2 = 0.00089, a3 = -2.834, a4 = -0.327, a5 = -0.091, a6 = 0.028, b1 = 1.278, b2 = -0.0715, b2p = 0.0496), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5 not significant
umcaDiameterFromHeightSibbesenFormRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.496, a2 = -0.188, b1 = 2.31, b2 = -0.12), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
#confint2(umcaDiameterFromHeightPowerPhysio, level = 0.99)

#umcaDiameterFromHeightChapmanForm = get_dbh_error(umcaDiameterFromHeightChapmanForm, umca2016, umca2016natural, umca2016plantation)
#umcaDiameterFromHeightChapmanFormAal = get_dbh_error(umcaDiameterFromHeightChapmanFormAal, umca2016, umca2016natural, umca2016plantation)
#umcaDiameterFromHeightChapmanFormBal = get_dbh_error(umcaDiameterFromHeightChapmanFormBal, umca2016, umca2016natural, umca2016plantation)
#umcaDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(umcaDiameterFromHeightChapmanFormBalRelHt, umca2016, umca2016natural, umca2016plantation)
#umcaDiameterFromHeightChapmanFormRelHt = get_dbh_error(umcaDiameterFromHeightChapmanFormRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichards = get_dbh_error(umcaDiameterFromHeightChapmanRichards, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichardsAal = get_dbh_error(umcaDiameterFromHeightChapmanRichardsAal, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(umcaDiameterFromHeightChapmanRichardsPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(umcaDiameterFromHeightChapmanRichardsRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightLinear = get_dbh_error(umcaDiameterFromHeightLinear, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightParabolic = get_dbh_error(umcaDiameterFromHeightParabolic, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPower = get_dbh_error(umcaDiameterFromHeightPower, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerAal = get_dbh_error(umcaDiameterFromHeightPowerAal, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerPhysio = get_dbh_error(umcaDiameterFromHeightPowerPhysio, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerRelHt = get_dbh_error(umcaDiameterFromHeightPowerRelHt, umca2016, umca2016natural, umca2016plantation)
#umcaDiameterFromHeightSharmaParton = get_dbh_error(umcaDiameterFromHeightSharmaParton, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenForm = get_dbh_error(umcaDiameterFromHeightSibbesenForm, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenFormAal = get_dbh_error(umcaDiameterFromHeightSibbesenFormAal, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenFormPhysio = get_dbh_error(umcaDiameterFromHeightSibbesenFormPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaDiameterFromHeightSibbesenFormFormRelHt = get_dbh_error(umcaDiameterFromHeightSibbesenFormRelHt, umca2016, umca2016natural, umca2016plantation)

umcaDiameterFromHeightResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(umcaDiameterFromHeightChapmanRichards),
                                        "Chapman-Richards AAL", !!!as_row(umcaDiameterFromHeightChapmanRichardsAal), # not AIC supported
                                        "Chapman-Richards physio", !!!as_row(umcaDiameterFromHeightChapmanRichardsPhysio),
                                        "Chapman-Richards RelHt", !!!as_row(umcaDiameterFromHeightChapmanRichardsRelHt), # not AIC supported
                                        "linear", !!!as_row(umcaDiameterFromHeightLinear),
                                        "parabolic", !!!as_row(umcaDiameterFromHeightParabolic),
                                        "power", !!!as_row(umcaDiameterFromHeightPower),
                                        "power AAL", !!!as_row(umcaDiameterFromHeightPowerAal), # not AIC supported
                                        "power physio", !!!as_row(umcaDiameterFromHeightPowerPhysio),
                                        "power RelHt", !!!as_row(umcaDiameterFromHeightPowerRelHt),
                                        "modified Sharma-Parton", !!!as_row(NULL),
                                        "Sibbesen form", !!!as_row(umcaDiameterFromHeightSibbesenForm),
                                        "Sibbesen form AAL", !!!as_row(umcaDiameterFromHeightSibbesenFormAal),
                                        "Sibbesen form physio", !!!as_row(umcaDiameterFromHeightSibbesenFormPhysio),
                                        "Sibbesen form RelHt", !!!as_row(umcaDiameterFromHeightSibbesenFormFormRelHt),
                                        "Chapman-Richards form", !!!as_row(NULL),
                                        "Chapman-Richards form AAL", !!!as_row(NULL),
                                        "Chapman-Richards form BAL", !!!as_row(NULL),
                                        "Chapman-Richards form BAL RelHt", !!!as_row(NULL),
                                        "Chapman-Richards form RelHt", !!!as_row(NULL)) %>%
  mutate(responseVariable = "height", species = "UMCA", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
umcaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic)

ggplot(umca2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = umcaDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = umcaDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = umcaDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = umcaDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 90, label = "Oregon myrtle, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
umcaParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", umcaHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", umcaHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", umcaHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", umcaHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", umcaHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", umcaHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", umcaHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", umcaHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", umcaHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "Michaelis-Menten", umcaHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", umcaHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", umcaHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", umcaHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", umcaHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", umcaHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", umcaHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", umcaHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", umcaHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", umcaHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", umcaHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", umcaHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", umcaHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", umcaHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", umcaHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", umcaHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards GNLS", umcaHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     bind_rows(c(method = "Chapman-Richards BAL GNLS", umcaHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     #bind_rows(c(method = "Sharma-Parton GNLS", umcaHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     #bind_rows(c(method = "Sharma-Parton BAL GNLS", umcaHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     #bind_rows(c(method = "Sharma-Zhang GNLS", umcaHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     #bind_rows(c(method = "Sharma-Zhang BAL GNLS", umcaHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", umcaHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", umcaHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", umcaDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", umcaDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", umcaDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", umcaDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form", umcaDiameterFromHeightChapmanForm$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form AAL", umcaDiameterFromHeightChapmanFormAal$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form BAL", umcaDiameterFromHeightChapmanFormBal$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form BAL RelHt", umcaDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form RelHt", umcaDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                     bind_rows(c(method = "linear", umcaDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "parabolic", umcaDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", umcaDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", umcaDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", umcaDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", umcaDiameterFromHeightPowerRelHt$m$getPars())),
                                     #bind_rows(c(method = "modified Sharma-Parton", umcaDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", umcaDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", umcaDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", umcaDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", umcaDiameterFromHeightSibbesenFormRelHt$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "UMCA",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

