# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western redcedar height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Parton, Sharma-Zhang, Chapman-Richards BAL
#thplHeightFromDiameterSharmaPartonBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b1p = -0.0143, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = pmin(DBH^-2, 1))
thpl2016physio = thpl2016 %>% filter(is.na(elevation) == FALSE)
thpl2016plantationPhysio = thpl2016physio %>% filter(isPlantation)
thplHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = pmin(DBH^-2, 1)) # b1p not significant
thplHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect)) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), thpl2016physio, start = list(a1 = 69.4, a2 = 0.077, a2p = 0.72, a3 = -0.008, a4 = -0.031, a5 = 0.548, a6 = 0.980, b1 = -0.021, b1p = 0.0075, b2 = 1.49, b2p = -0.370), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 8.7, a1p = 11.1, a2 = 0.18, a2p = 0.42, a3 = -0.0087, a3p = 0.070, a4 = 54.0, a4p = -28.4, b1 = -0.021, b2 = 0.65, b2p = 0.37), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016physio, start = list(a1 = 68.5, a1p = -13.4, a2 = -0.0045, a3 = -8.09, a4 = 0.783, a5 = 0.766, a6 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^-2, 1)) # a4p not significant, a5p induces overfitting
thplHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016, start = list(a1 = 1.6, b1 = 0.24), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), thpl2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
thplHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 72.1, a1p = -20.4, a2 = 208, a2p = -75.5, b1 = 1.18), weights = pmin(DBH^-2, 1)) # b1p not significant
thplHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), thpl2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
thplHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), thpl2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), thpl2016, start = list(a1 = 2.9, a1p = -1.6, b1 = 0.63, b1p = 0.18), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), thpl2016, start = list(Ha = 55.7, Hap = -29.2, d = 0.517, dp = 0.318, kU = 0.008, kUp = 0.016), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = pmin(DBH^-2, 1)) # step factor with b1p
thplHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 18.2, a1p = 7.93, a2 = 0.281, a2p = -0.187, b1 = -0.322, b1p = 0.293, b2 = -0.505, b2p = 0.463, b3 = 1.133, b3p = 0.020), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 18.2, a1p = 7.93, a2 = 0.281, a2p = -0.187, a3 = 0, a3p = 0, b1 = -0.322, b1p = 0.293, b2 = -0.505, b2p = 0.463, b3 = 1.133, b3p = 0.020), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 0.0019, a1p = 0.163, b1 = 4.98, b1p = -2.72, b2 = -0.175, b2p = 0.0427), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), thpl2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), thpl2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), weights = pmin(DBH^-2, 1))
thplHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1  + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * (1 + (b2 + b2p * isPlantation) * pmin(relativeHeight, 1.25)) * DBH^b3)), thpl2016, start = list(a1 = -3.82, a1p = 87.0, a2 = 0.117, a2p = 1.86, a3 = -0.00099, a3p = -0.0684, a4 = 68.2, a4p = -53.9, b1 = -0.146, b1p = 0.140, b2 = -0.610, b2p = 2.25, b3 = 0.800), weights = pmin(DBH^-2, 1))

thplHeightFromDiameterChapmanRichards = get_height_error(thplHeightFromDiameterChapmanRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBal = get_height_error(thplHeightFromDiameterChapmanRichardsBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(thplHeightFromDiameterChapmanRichardsBalPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(thplHeightFromDiameterChapmanRichardsBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsPhysio = get_height_error(thplHeightFromDiameterChapmanRichardsPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterCurtis = get_height_error(thplHeightFromDiameterCurtis, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterHossfeld = get_height_error(thplHeightFromDiameterHossfeld, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterKorf = get_height_error(thplHeightFromDiameterKorf, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterLinear = get_height_error(thplHeightFromDiameterLinear, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterMichaelisMenten = get_height_error(thplHeightFromDiameterMichaelisMenten, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterParabolic = get_height_error(thplHeightFromDiameterParabolic, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterProdan = get_height_error(thplHeightFromDiameterProdan, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterPower = get_height_error(thplHeightFromDiameterPower, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterRatkowsky = get_height_error(thplHeightFromDiameterRatkowsky, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterRichards = get_height_error(thplHeightFromDiameterRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaParton = get_height_error(thplHeightFromDiameterSharmaParton, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBal = get_height_error(thplHeightFromDiameterSharmaPartonBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBalPhysio = get_height_error(thplHeightFromDiameterSharmaPartonBalPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterSharmaPartonPhysio = get_height_error(thplHeightFromDiameterSharmaPartonPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplHeightFromDiameterSharmaZhang = get_height_error(thplHeightFromDiameterSharmaZhang, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangBal = get_height_error(thplHeightFromDiameterSharmaZhangBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSibbesen = get_height_error(thplHeightFromDiameterSibbesen, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibull = get_height_error(thplHeightFromDiameterWeibull, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBal = get_height_error(thplHeightFromDiameterWeibullBal, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBalRelHt = get_height_error(thplHeightFromDiameterWeibullBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)

thplHeightFromDiameterResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(thplHeightFromDiameterChapmanRichards),
                                        "Chapman-Richards BAL", !!!as_row(thplHeightFromDiameterChapmanRichardsBal),
                                        "Chapman-Richards BAL physio", !!!as_row(thplHeightFromDiameterChapmanRichardsBalPhysio),
                                        "Chapman-Richards BAL RelHt", !!!as_row(thplHeightFromDiameterChapmanRichardsBalRelHt),
                                        "Chapman-Richards physio", !!!as_row(thplHeightFromDiameterChapmanRichardsPhysio),
                                        "Curtis", !!!as_row(thplHeightFromDiameterCurtis),
                                        "Hossfeld", !!!as_row(thplHeightFromDiameterHossfeld),
                                        "Korf", !!!as_row(thplHeightFromDiameterKorf),
                                        "linear", !!!as_row(thplHeightFromDiameterLinear),
                                        "generalized Michaelis-Menten", !!!as_row(thplHeightFromDiameterMichaelisMenten),
                                        "parabolic", !!!as_row(thplHeightFromDiameterParabolic),
                                        "power", !!!as_row(thplHeightFromDiameterPower),
                                        "Prodan", !!!as_row(thplHeightFromDiameterProdan),
                                        "Ratkowsky", !!!as_row(thplHeightFromDiameterRatkowsky),
                                        "unified Richards", !!!as_row(thplHeightFromDiameterRichards),
                                        "Sharma-Parton", !!!as_row(thplHeightFromDiameterSharmaParton),
                                        "Sharma-Parton BAL", !!!as_row(thplHeightFromDiameterSharmaPartonBal),
                                        "Sharma-Parton BAL physio", !!!as_row(thplHeightFromDiameterSharmaPartonBalPhysio),
                                        "Sharma-Parton physio", !!!as_row(thplHeightFromDiameterSharmaPartonPhysio),
                                        "Sharma-Zhang", !!!as_row(thplHeightFromDiameterSharmaZhang),
                                        "Sharma-Zhang BAL", !!!as_row(thplHeightFromDiameterSharmaZhangBal),
                                        "Sibbesen", !!!as_row(thplHeightFromDiameterSibbesen),
                                        "Weibull", !!!as_row(thplHeightFromDiameterWeibull),
                                        "Weibull BAL", !!!as_row(thplHeightFromDiameterWeibullBal),
                                        "Weibull BAL RelHt", !!!as_row(thplHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "THPL", deltaAic = aic - min(aic)) %>%
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
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = thpl2016$isPlantation)) +
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
#thplHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), thpl2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#thplHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), thpl2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(thplHeightFromDiameterChapmanRichardsGnls, thplHeightFromDiameterChapmanRichardsBalGnls, thplHeightFromDiameterSharmaPartonGnls, thplHeightFromDiameterSharmaPartonBalGnls, thplHeightFromDiameterSharmaZhangGnls, thplHeightFromDiameterSharmaZhangBalGnls, thplHeightFromDiameterWeibullGnls, thplHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia THPL GNLS.rdata")
load("trees/height-diameter/HtDia THPL GNLS.rdata")
thplHeightFromDiameterWeibullGnls = thplHeightFromDiameterWykoffGnls # temporary naming error fixup
thplHeightFromDiameterWeibullBalGnls = thplHeightFromDiameterWykoffBalGnls

thplHeightFromDiameterChapmanRichardsGnls = get_height_error(thplHeightFromDiameterChapmanRichardsGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterChapmanRichardsBalGnls = get_height_error(thplHeightFromDiameterChapmanRichardsBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonGnls = get_height_error(thplHeightFromDiameterSharmaPartonGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaPartonBalGnls = get_height_error(thplHeightFromDiameterSharmaPartonBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangGnls = get_height_error(thplHeightFromDiameterSharmaZhangGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterSharmaZhangBalGnls = get_height_error(thplHeightFromDiameterSharmaZhangBalGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullGnls = get_height_error(thplHeightFromDiameterWeibullGnls, thpl2016, thpl2016natural, thpl2016plantation)
thplHeightFromDiameterWeibullBalGnls = get_height_error(thplHeightFromDiameterWeibullBalGnls, thpl2016, thpl2016natural, thpl2016plantation)

thplHeightFromDiameterResultsGnls = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                            "Chapman-Richards GNLS", !!!as_row(thplHeightFromDiameterChapmanRichardsGnls),
                                            "Chapman-Richards BAL GNLS", !!!as_row(thplHeightFromDiameterChapmanRichardsBalGnls),
                                            "Sharma-Parton GNLS", !!!as_row(thplHeightFromDiameterSharmaPartonGnls),
                                            "Sharma-Parton BAL GNLS", !!!as_row(thplHeightFromDiameterSharmaPartonBalGnls),
                                            "Sharma-Zhang GNLS", !!!as_row(thplHeightFromDiameterSharmaZhangGnls),
                                            "Sharma-Zhang BAL GNLS", !!!as_row(thplHeightFromDiameterSharmaZhangBalGnls),
                                            "Weibull GNLS", !!!as_row(thplHeightFromDiameterWeibullGnls),
                                            "Weibull BAL GNLS", !!!as_row(thplHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "THPL", deltaAic = aic - min(aic)) %>%
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
#                                                  start_upper = list(a1 = 50, b1 = 1, b2 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#thplDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999999)), thpl2016, iter = 100,
#                                                      start_lower = list(a1 = 1, b1 = -0.5, b2 = 0.1), 
#                                                      start_upper = list(a1 = 300, b1 = 0.5, b2 = 2.5), modelweights = pmin(TotalHt^-2, 0.5))
#thplDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), thpl2016, start = list(a1 = 229.6, b1 = -0.0057, b2 = 1.25, b2p = 0), weights = pmin(TotalHt^-2, 0.5))
#thplDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), thpl2016, iter = 100,
#                                                            start_lower = list(a1 = -200, a1p = -100, a2 = -0.1, a3 = -1, a4 = -1, a5 = -1, a6 = -0.1, b1 = -0.1, b1p = -1, b2 = -1), 
#                                                            start_upper = list(a1 = 1, a1p = 100, a2 = 0.1, a3 = 1, a4 = 1, a5 = 1, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#thplDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, thpl2016, iter = 100,
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = -0.5, b2 = -0.5, b3 = -1), 
#                                                   start_upper = list(a1 = 10, a2 = 1.5, b1 = 0.5, b2 = 0.5, b3 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#thplDiameterFromHeightChapmanForm = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 15, b1 = 0.01, b2 = 0.45), weights = pmin(TotalHt^-2, 0.5)) # singular gradient, no convergence from nls_multstart()
#thplDiameterFromHeightChapmanFormAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 15, a2 = -0.0045, b1 = 0.010, b2 = 0.924), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf
#thplDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 1000, a2 = -7, b1 = 0.001, b2 = 1.0), weights = pmin(TotalHt^-2, 0.5), trace = TRUE) # step size
#thplDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 1396, a2 = -22.0, a3 = 10.4, a4 = -467, b1 = 0.0017, b2 = 0.971), weights = pmin(TotalHt^-2, 0.5)) # step size
#thplDiameterFromHeightChapmanFormRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 1000, a2 = 100, b1 = 0.00154, b2 = 1.0), weights = pmin(TotalHt^-2, 0.5), trace = TRUE) # step size
thplDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -763, b1 = 0.0027, b2 = 1.05), weights = pmin(TotalHt^-2, 0.5)) # a1p and b2p not significant, poor convergence with b1p
thplDiameterFromHeightChapmanRichardsAal = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -1000, a2 = 4.9, b1 = 0.001, b2 = 1.07), weights = pmin(TotalHt^-2, 0.5)) # a1p, b1p not significant, step factor with nlrob()
thplDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016physio, start = list(a1 = -500, a1p = 200, a2 = 0.07, a3 = 9.8, a4 = -7.2, a5 = -1.3, a6 = 0.7, b1 = 0.0047, b1p = 0.007, b2 = 1.005), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50)) # no physiographic effects significant
thplDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999999)), thpl2016, start = list(a1 = -350, a2 = -66, b1 = 0.006, b2 = 0.98), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50))
thplDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), thpl2016, weights = TotalHt^-2)
thplDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), thpl2016, weights = TotalHt^-2)
thplDiameterFromHeightPower = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
thplDiameterFromHeightPowerAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
thplDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, start = list(a1 = 2.26, a2 = -0.00047, a3 = -0.156, a4 = -0.0018, a5 = 0.0208, a6 = -0.0060, b1 = 1.08), weights = pmin(TotalHt^-2, 0.5)) # no significant physiographic effects
thplDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.68, a2 = -0.11, a2p = 0.23, b1 = 1.13), weights = pmin(TotalHt^-2, 0.5)) # a1p and b1p not significant
#thplDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, thpl2016, start = list(a1 = 0.91, a2 = 1.39, b1 = 0.087, b2 = -0.029, b3 = -0.12), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf from nls_multstart() point
#thplDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), thpl2016, start = list(a1 = 2.86, a2 = 0.790, a2p = -0.130, b1 = 0.216, b2 = -0.153, b2p = 0.184, b3 = 0.045, b3p = 0.0278), weights = pmin(TotalHt^-2, 0.5)) # step factor
thplDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.42, b1 = 1.29, b2 = -0.27), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
thplDiameterFromHeightSibbesenFormAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
thplDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), thpl2016physio, start = list(a1 = 2.32, a1p = -0.294, a2 = -0.00069, a3 = -0.116, a4 = 0.0374, a5 = 0.0047, a6 = -0.0056, b1 = 1.084, b2 = -0.0048, b2p = 0.0203), weights = pmin(TotalHt^-2, 0.5)) # no physiographic effects significant
thplDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.34, a2 = 0.30, b1 = 1.32, b2 = -0.040), weights = pmin(TotalHt^-2, 0.5))
#confint2(thplDiameterFromHeightPowerPhysio, level = 0.99)

#thplDiameterFromHeightChapmanForm = get_dbh_error(thplDiameterFromHeightChapmanForm, thpl2016, thpl2016natural, thpl2016plantation)
#thplDiameterFromHeightChapmanFormAal = get_dbh_error(thplDiameterFromHeightChapmanFormAal, thpl2016, thpl2016natural, thpl2016plantation)
#thplDiameterFromHeightChapmanFormBal = get_dbh_error(thplDiameterFromHeightChapmanFormBal, thpl2016, thpl2016natural, thpl2016plantation)
#thplDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(thplDiameterFromHeightChapmanFormBalRelHt, thpl2016, thpl2016natural, thpl2016plantation)
#thplDiameterFromHeightChapmanFormRelHt = get_dbh_error(thplDiameterFromHeightChapmanFormRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichards = get_dbh_error(thplDiameterFromHeightChapmanRichards, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichardsAal = get_dbh_error(thplDiameterFromHeightChapmanRichardsAal, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(thplDiameterFromHeightChapmanRichardsPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(thplDiameterFromHeightChapmanRichardsRelHt, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightLinear = get_dbh_error(thplDiameterFromHeightLinear, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightParabolic = get_dbh_error(thplDiameterFromHeightParabolic, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPower = get_dbh_error(thplDiameterFromHeightPower, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerAal = get_dbh_error(thplDiameterFromHeightPowerAal, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerPhysio = get_dbh_error(thplDiameterFromHeightPowerPhysio, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightPowerRelHt = get_dbh_error(thplDiameterFromHeightPowerRelHt, thpl2016, thpl2016natural, thpl2016plantation)
#thplDiameterFromHeightSharmaParton = get_dbh_error(thplDiameterFromHeightSharmaParton, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenForm = get_dbh_error(thplDiameterFromHeightSibbesenForm, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenFormAal = get_dbh_error(thplDiameterFromHeightSibbesenFormAal, thpl2016, thpl2016natural, thpl2016plantation)
thplDiameterFromHeightSibbesenFormPhysio = get_dbh_error(thplDiameterFromHeightSibbesenFormPhysio, thpl2016physio, thpl2016natural, thpl2016plantationPhysio)
thplDiameterFromHeightSibbesenFormRelHt = get_dbh_error(thplDiameterFromHeightSibbesenFormRelHt, thpl2016, thpl2016natural, thpl2016plantation)

thplDiameterFromHeightResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(thplDiameterFromHeightChapmanRichards),
                                        "Chapman-Richards AAL", !!!as_row(thplDiameterFromHeightChapmanRichardsAal),
                                        "Chapman-Richards physio", !!!as_row(thplDiameterFromHeightChapmanRichardsPhysio),
                                        "Chapman-Richards RelHt", !!!as_row(thplDiameterFromHeightChapmanRichardsRelHt),
                                        "linear", !!!as_row(thplDiameterFromHeightLinear),
                                        "parabolic", !!!as_row(thplDiameterFromHeightParabolic),
                                        "power", !!!as_row(thplDiameterFromHeightPower),
                                        "power AAL", !!!as_row(thplDiameterFromHeightPowerAal),
                                        "power physio", !!!as_row(thplDiameterFromHeightPowerPhysio),
                                        "power RelHt", !!!as_row(thplDiameterFromHeightPowerRelHt),
                                        #bind_rows(c(method = "modified Sharma-Parton", thplDiameterFromHeightSharmaParton$m$getPars())),
                                        "Sibbesen form", !!!as_row(thplDiameterFromHeightSibbesenForm),
                                        "Sibbesen form AAL", !!!as_row(thplDiameterFromHeightSibbesenFormAal),
                                        "Sibbesen form physio", !!!as_row(thplDiameterFromHeightSibbesenFormPhysio),
                                        "Sibbesen form RelHt", !!!as_row(thplDiameterFromHeightSibbesenFormRelHt),
                                        "Chapman-Richards form", !!!as_row(NULL),
                                        "Chapman-Richards form AAL", !!!as_row(NULL),
                                        "Chapman-Richards form BAL", !!!as_row(NULL),
                                        "Chapman-Richards form BAL RelHt", !!!as_row(NULL),
                                        "Chapman-Richards form RelHt", !!!as_row(NULL)) %>%
  mutate(responseVariable = "height", species = "THPL", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
thplDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic)

ggplot(thpl2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = thplDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = thplDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = thplDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = thplDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  #geom_line(aes(x = thplDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -100 * log(1 - pmin(0.015*(TotalHt - 1.37)^1.0, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
  #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.45, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 90, label = "western redcedar, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
thplParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", thplHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", thplHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", thplHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", thplHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", thplHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", thplHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", thplHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", thplHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", thplHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten", thplHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", thplHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", thplHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", thplHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", thplHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", thplHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", thplHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", thplHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", thplHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", thplHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", thplHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", thplHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", thplHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", thplHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", thplHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", thplHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards GNLS", thplHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     bind_rows(c(method = "Chapman-Richards BAL GNLS", thplHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton GNLS", thplHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton BAL GNLS", thplHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang GNLS", thplHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang BAL GNLS", thplHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", thplHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", thplHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", thplDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", thplDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", thplDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", thplDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     bind_rows(c(method = "linear", thplDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "parabolic", thplDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", thplDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", thplDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", thplDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", thplDiameterFromHeightPowerRelHt$m$getPars())),
                                     #bind_rows(c(method = "modified Sharma-Parton", thplDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", thplDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", thplDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", thplDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", thplDiameterFromHeightSibbesenFormRelHt$m$getPars()))) %>%
                                     #bind_rows(c(method = "Chapman-Richards form", thplDiameterFromHeightChapmanForm$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form AAL", thplDiameterFromHeightChapmanFormAal$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form BAL", thplDiameterFromHeightChapmanFormBal$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form BAL RelHt", thplDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form RelHt", thplDiameterFromHeightChapmanFormRelHt$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "THPL",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

