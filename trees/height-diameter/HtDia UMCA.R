# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Oregon myrtle (California bay) height-diameter regression form sweep
# preferred forms: Sibbesen, Korf, Prodan, Ratkowsky
#umcaHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), umca2016, start = list(Ha = 14.4, Hap = -3.6, d = 2.78, dp = -0.99, kU = 0.054, kUp = 0.026), weights = pmin(DBH^-2, 1))
umca2016physio = umca2016 %>% filter(is.na(elevation) == FALSE)
umca2016plantationPhysio = umca2016physio %>% filter(isPlantation)
umcaHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 16.5, b1 = -0.064, b2 = 1.270, b2p = -0.177), weights = pmin(DBH^-2, 1)) # a1p, b1p not significant
umcaHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, umca2016, start = list(a1 = 19.0, a2 = 0.201, a3 = -0.152, b1 = -0.049, b1p = -0.0105, b2 = 1.156), weights = pmin(DBH^-2, 1)) # a2, a2p, a3, a3p, b2p not significant
umcaHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 8.62, a2 = 0.050, a3 = -0.014, a4 = 0.230, a5 = 0.413, a6 = -1.113, a7 = 0, b1 = -0.062, b2 = 1.265, b2p = -0.219), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, a5, a6, b1p not significant
umcaHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = -1.48, a2 = 0.030, a2p = -0.131, a3 = -0.009, a3p = 0.192, a4 = 56.1, a4p = -32.7, b1 = -0.278, b2 = 0.389, b2p = 0.684), weights = pmin(DBH^-2, 1)) # a1p not significant
umcaHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 8.58, a2 = -0.018, a3 = 16.5, a4 = -1.117, a5 = 0.476, a6 = -0.014, b1 = -0.064, b2 = 1.261, b2p = -0.193), weights = pmin(DBH^-2, 1)) # a4, a5, a6, b1p not significant
umcaHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^b1, umca2016, start = list(a1 = 0.909, a1p = 0.221, b1 = 0.219), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), umca2016, start = list(a1 = 21.4, a1p = -4.37, b1 = 43.8, b1p = -19.9, b2 = -1.27), weights = pmin(DBH^-2, 1)) # b2p not significant
umcaHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*DBH^b2), umca2016, start = list(a1 = 49.8, b1 = -5.02, b1p = 0.404, b2 = -0.386), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
umcaHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), umca2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), umca2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
umcaHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016, start = list(a1 = 0.037, a1p = 0.020, a2 = 1.109, a2p = -0.472, a3 = 1.242), weights = pmin(DBH^-2, 1)) # a3p not significant
umcaHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, umca2016, start = list(a1 = 0.821, a1p = 0.206, b1 = 0.810), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(DBH + b2)), umca2016, start = list(a1 = 21.6, b1 = -16.5, b1p = 1.974, b2 = 3.629), weights = pmin(DBH^-2, 1)) # a1p, b2p not significant
umcaHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap * isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016, start = list(Ha = 14.6, Hap = -4.146, d = 2.198, kU = 0.048, kUp = 0.045), weights = pmin(DBH^-2, 1)) # dp not significant
umcaHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, umca2016, start = list(a1 = 5.76, a2 = 0.261, b1 = -0.049, b1p = -0.024, b2 = 0.085, b3 = 1.234), weights = pmin(DBH^-2, 1)) # a1p, a2p, b2p, b3p not significant
umcaHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 9.074, a2 = 0.149, b1 = -0.065, b2 = 0.013, b3 = 1.299, b3p = -0.254), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, b1p, b2p not significant
umcaHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), umca2016physio, start = list(a1 = 9.30, a2 = 0.151, a3 = -0.0007, a4 = 0.0161, a5 = -0.040, a6 = 0.005, b1 = -0.061, b2 = 0.061, b3 = 1.313, b3p = -0.253), weights = pmin(DBH^-2, 1)) # a2, a4, a5, a6, b1p, b2p not significant
umcaHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare))^b2*DBH))^(b3 + b3p * isPlantation), umca2016physio, start = list(a1 = 8.639, a2 = 0.169, a3 = -0.0007, a4 = 0.018, a5 = -0.044, a6 = 0.005, b1 = -0.052, b2 = 0.109, b3 = 1.315, b3p = -0.247), weights = pmin(DBH^-2, 1)) # a2p, a4, a5, a6, b1p, b2, b2p not significant
umcaHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 13.6, a2 = 0.054, b1 = -0.030, b2 = 0.125, b3 = 1.295, b3p = -0.232), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, b1p, b2p not significant
umcaHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, umca2016, start = list(a1 = 8.742, a2 = 0.174, a3 = -0.0001, a3p = -0.007, b1 = -0.039, b2 = 0.070, b2p = 0.090, b3 = 1.258), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, a3, b1p, b3p not significant
umcaHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.346, a1p = 0.153, b1 = 1.783, b2 = -0.150, b2p = -0.033), weights = pmin(DBH^-2, 1)) # b1p not significant
umcaHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = list(a1 = 16.7, a1p = -3.627, b1 = -0.0315, b1p = -0.0256, b2 = 1.174), weights = pmin(DBH^-2, 1)) # b2p not significant
umcaHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = list(a1 = 18.6, a2 = 0.195, a3 = -0.148, b1 = -0.031, b1p = -0.008, b2 = 1.118), weights = pmin(DBH^-2, 1)) # a1p, a2, a2p, a3p, b2p not significant
umcaHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), umca2016, start = list(a1 = -1.490, a2 = 0.029, a2p = -0.127, a3 = -0.009, a3p = 0.195, a4 = 55.7, a4p = -31.3, b1 = -0.628, b1p = 0.366, b2 = 0.913), weights = pmin(DBH^-2, 1)) # a1p, a3, b2p not significant
#confint2(umcaHeightFromDiameterWeibullBalRelHt, level = 0.99)

umcaHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", umcaHeightFromDiameterChapmanRichards, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BAL", umcaHeightFromDiameterChapmanRichardsBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BAL physio", umcaHeightFromDiameterChapmanRichardsBalPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BAL RelHt", umcaHeightFromDiameterChapmanRichardsBalRelHt, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", umcaHeightFromDiameterChapmanRichardsPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterCurtis = get_height_error("Curtis", umcaHeightFromDiameterCurtis, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", umcaHeightFromDiameterHossfeld, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterKorf = get_height_error("Korf", umcaHeightFromDiameterKorf, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterLinear = get_height_error("linear", umcaHeightFromDiameterLinear, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", umcaHeightFromDiameterMichaelisMenten, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterParabolic = get_height_error("parabolic", umcaHeightFromDiameterParabolic, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterProdan = get_height_error("Prodan", umcaHeightFromDiameterProdan, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterPower = get_height_error("power", umcaHeightFromDiameterPower, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", umcaHeightFromDiameterRatkowsky, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterRichards = get_height_error("unified Richards", umcaHeightFromDiameterRichards, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", umcaHeightFromDiameterSharmaParton, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BAL", umcaHeightFromDiameterSharmaPartonBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BAL physio", umcaHeightFromDiameterSharmaPartonBalPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", umcaHeightFromDiameterSharmaPartonPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", umcaHeightFromDiameterSharmaZhang, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BAL", umcaHeightFromDiameterSharmaZhangBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterSibbesen = get_height_error("Sibbesen", umcaHeightFromDiameterSibbesen, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibull = get_height_error("Weibull", umcaHeightFromDiameterWeibull, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBal = get_height_error("Weibull BAL", umcaHeightFromDiameterWeibullBal, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull RelHt", umcaHeightFromDiameterWeibullBalRelHt, umca2016, umca2016natural, umca2016plantation)

umcaHeightFromDiameterResults = bind_rows(as_row(umcaHeightFromDiameterChapmanRichards),
                                          as_row(umcaHeightFromDiameterChapmanRichardsBal),
                                          as_row(umcaHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row(umcaHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(umcaHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(umcaHeightFromDiameterCurtis),
                                          as_row(umcaHeightFromDiameterHossfeld),
                                          as_row(umcaHeightFromDiameterKorf),
                                          as_row(umcaHeightFromDiameterLinear),
                                          as_row(umcaHeightFromDiameterMichaelisMenten),
                                          as_row(umcaHeightFromDiameterParabolic),
                                          as_row(umcaHeightFromDiameterPower),
                                          as_row(umcaHeightFromDiameterProdan),
                                          as_row(umcaHeightFromDiameterRatkowsky),
                                          as_row(umcaHeightFromDiameterRichards),
                                          as_row(umcaHeightFromDiameterSharmaParton),
                                          as_row(umcaHeightFromDiameterSharmaPartonBal),
                                          as_row(umcaHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(umcaHeightFromDiameterSharmaPartonPhysio),
                                          as_row(umcaHeightFromDiameterSharmaZhang),
                                          as_row(umcaHeightFromDiameterSharmaZhangBal),
                                          as_row(umcaHeightFromDiameterSibbesen),
                                          as_row(umcaHeightFromDiameterWeibull),
                                          as_row(umcaHeightFromDiameterWeibullBal),
                                          as_row(umcaHeightFromDiameterWeibullBalRelHt)) %>%
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
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = umca2016$isPlantation)) +
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

#umcaHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", umcaHeightFromDiameterChapmanRichardsGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BAL GNLS", umcaHeightFromDiameterChapmanRichardsBalGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", umcaHeightFromDiameterSharmaPartonGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BAL GNLS", umcaHeightFromDiameterSharmaPartonBalGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", umcaHeightFromDiameterSharmaZhangGnls, umca2016, umca2016natural, umca2016plantation)
#umcaHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BAL GNLS", umcaHeightFromDiameterSharmaZhangBalGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", umcaHeightFromDiameterWeibullGnls, umca2016, umca2016natural, umca2016plantation)
umcaHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BAL GNLS", umcaHeightFromDiameterWeibullBalGnls, umca2016, umca2016natural, umca2016plantation)

umcaHeightFromDiameterResultsGnls = bind_rows(as_row(name = "Chapman-Richards GNLS"),
                                              as_row(umcaHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(name = "Sharma-Parton GNLS"),
                                              as_row(name = "Sharma-Parton BAL GNLS"),
                                              as_row(name = "Sharma-Zhang GNLS"),
                                              as_row(name = "Sharma-Zhang BAL GNLS"),
                                              as_row(umcaHeightFromDiameterWeibullGnls),
                                              as_row(umcaHeightFromDiameterWeibullBalGnls)) %>%
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
#umcaDiameterFromHeightChapmanFormAat = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 15, a2 = 0, b1 = 0.1, b2 = 0.5), weights = pmin(TotalHt^-2, 0.5), control = nls.control(maxiter = 500)) # step factor with nlrob()
#umcaDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, umca2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.01, a2 = -1, b1 = -10, b2 = -0.5, b3 = 0.2),
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 10, b2 = 0.5, b3 = 1.5), modelweights = pmin(TotalHt^-2, 0.5))
#umcaDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 89, a2 = 0.61, a2p = 0.05, b1 = 0.0002, b2 = -0.73, b2p = -0.79, b3 = 0.34, b3p = -0.03), weights = pmin(TotalHt^-2, 0.5)) # singnular gradient with nls()
umcaDiameterFromHeightChapmanForm = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 50000, b1 = 0.00001, b2 = 1.034), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf with nls() at multiple nls_multstart() points, NaN-inf with nlrob()
umcaDiameterFromHeightChapmanFormAat = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 604, a2 = 0.51, b1 = 0.0004, b2 = 1.01), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with nls(), step factor with nlrob()
umcaDiameterFromHeightChapmanFormBal = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 1000, a2 = -16, b1 = 0.001, b2 = 1.02), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with nls(), step factor with nlrob()
umcaDiameterFromHeightChapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 1500, a2 = -400.0, a3 = 520, a4 = -2000, b1 = 0.0006, b2 = 1.02), weights = pmin(TotalHt^-2, 0.5)) # step factor with nls() and nlrob()
umcaDiameterFromHeightChapmanFormRelHt = gsl_nls(DBH ~ (a1 + a2 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 665, a2 = -406, b1 = 0.0034, b2 = 1.12), weights = pmin(TotalHt^-2, 0.5), control = nls.control(maxiter = 50)) # step factor with nls(), step factor with nlrob()
umcaDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 33.9, b1 = -0.047, b2 = 1.43, b2p = -0.15), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 34.4, a2 = 0.004, b1 = -0.047, b2 = 1.41, b2p = -0.13), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), umca2016physio, start = list(a1 = 58.6, a1p = 166.3, a2 = 0.0133, a3 = -52.7, a4 = 5.394, a5 = -0.837, a6 = 0.602, b1 = -0.0637, b1p = 0.0584, b2 = 1.329), weights = pmin(TotalHt^-2, 0.5)) # a2, a3, a4, a5 not significant
umcaDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 39.0, a2 = -2.62, b1 = -0.043, b2 = 1.40, b2p = -0.132), weights = pmin(TotalHt^-2, 0.5)) # a1p convergence questionable, b1p not significant
umcaDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), umca2016, weights = TotalHt^-2)
umcaDiameterFromHeightMichaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), umca2016, start = list(a1 = 100, a2 = 100, b1 = 1), weights = TotalHt^-2) # collapses to linear
umcaDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), umca2016, start = list(a1 = 4.3, a1p = -1.8, a2 = -0.14, a2p = -0.038), weights = TotalHt^-2)
umcaDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), umca2016, weights = TotalHt^-2) # collapses to linear since (TotalHt - 1.37)^2 and isPlantation*(TotalHt - 1.37)^2 not significant
umcaDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332), weights = pmin(TotalHt^-2, 0.5))
umcaDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p not significant
umcaDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016physio, start = list(a1 = 5.12, a1p = -0.67, a2 = 0.0013, a3 = -4.15, a4 = 0.46, a5 = -0.11, a6 = 0.044, b1 = 0.93), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2, b1p not significant
umcaDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 2.37, a1p = -0.887, a2 = -1.508, a2p = 1513, b1 = 1.123), weights = pmin(TotalHt^-2, 0.5)) 
umcaDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 2.48, b1 = 1.14, b1p = -0.57, b2 = -0.019, b2p = 0.092), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
umcaDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), umca2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 177), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with nlrob()
umcaDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, umca2016, start = list(a1 = 104, a2 = 0.66, b1 = 0.0001, b2 = -1.17, b3 = 0.31), weights = pmin(TotalHt^-2, 0.5)) # a2p, b2p, b3p not significant, NaN-inf with nls() at nls_multistart() point, NaN-inf or code error with nlrob()
umcaDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
umcaDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
umcaDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), umca2016physio, start = list(a1 = 3.688, a1p = -0.999, a2 = 0.00089, a3 = -2.834, a4 = -0.327, a5 = -0.091, a6 = 0.028, b1 = 1.278, b2 = -0.0715, b2p = 0.0496), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5 not significant
umcaDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.496, a2 = -0.188, b1 = 2.31, b2 = -0.12), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
umcaDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, umca2016, start = list(a1 = -3800, b1 = 0.0006, b2 = 1.03), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf with nlrob()
#confint2(umcaDiameterFromHeightWeibull, level = 0.99)

umcaDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", umcaDiameterFromHeightChapmanForm, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AAT", umcaDiameterFromHeightChapmanFormAat, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BAL", umcaDiameterFromHeightChapmanFormBal, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BAL RelHt", umcaDiameterFromHeightChapmanFormBalRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", umcaDiameterFromHeightChapmanFormRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", umcaDiameterFromHeightChapmanRichards, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AAT", umcaDiameterFromHeightChapmanRichardsAat, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", umcaDiameterFromHeightChapmanRichardsPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", umcaDiameterFromHeightChapmanRichardsRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightLinear = get_dbh_error("linear", umcaDiameterFromHeightLinear, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", umcaDiameterFromHeightMichaelisMentenForm, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightNaslund = get_dbh_error("Näslund", umcaDiameterFromHeightNaslund, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightParabolic = get_dbh_error("parabolic", umcaDiameterFromHeightParabolic, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPower = get_dbh_error("power", umcaDiameterFromHeightPower, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerAat = get_dbh_error("power AAT", umcaDiameterFromHeightPowerAat, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerPhysio = get_dbh_error("power physio", umcaDiameterFromHeightPowerPhysio, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", umcaDiameterFromHeightPowerRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightRuark = get_dbh_error("Ruark", umcaDiameterFromHeightRuark, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSchnute = get_dbh_error("Schnute", umcaDiameterFromHeightSchnute, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", umcaDiameterFromHeightSharmaParton, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", umcaDiameterFromHeightSibbesenForm, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AAT", umcaDiameterFromHeightSibbesenFormAat, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", umcaDiameterFromHeightSibbesenFormPhysio, umca2016physio, umca2016natural, umca2016plantationPhysio)
umcaDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", umcaDiameterFromHeightSibbesenFormRelHt, umca2016, umca2016natural, umca2016plantation)
umcaDiameterFromHeightWeibull = get_dbh_error("Weibull", umcaDiameterFromHeightWeibull, umca2016, umca2016natural, umca2016plantation)

umcaDiameterFromHeightResults = bind_rows(as_row(umcaDiameterFromHeightChapmanRichards),
                                          as_row(umcaDiameterFromHeightChapmanRichardsAat),
                                          as_row(umcaDiameterFromHeightChapmanRichardsPhysio),
                                          as_row(umcaDiameterFromHeightChapmanRichardsRelHt),
                                          as_row(umcaDiameterFromHeightChapmanForm),
                                          as_row(umcaDiameterFromHeightChapmanFormAat),
                                          as_row(umcaDiameterFromHeightChapmanFormBal),
                                          as_row(umcaDiameterFromHeightChapmanFormBalRelHt),
                                          as_row(umcaDiameterFromHeightChapmanFormRelHt),
                                          as_row(umcaDiameterFromHeightLinear),
                                          as_row(umcaDiameterFromHeightMichaelisMentenForm),
                                          as_row(umcaDiameterFromHeightNaslund),
                                          as_row(umcaDiameterFromHeightParabolic, significant = FALSE),
                                          as_row(umcaDiameterFromHeightPower),
                                          as_row(umcaDiameterFromHeightPowerAat),
                                          as_row(umcaDiameterFromHeightPowerPhysio),
                                          as_row(umcaDiameterFromHeightPowerRelHt),
                                          as_row(umcaDiameterFromHeightRuark),
                                          as_row(umcaDiameterFromHeightSchnute),
                                          as_row(umcaDiameterFromHeightSharmaParton),
                                          as_row(umcaDiameterFromHeightSibbesenForm),
                                          as_row(umcaDiameterFromHeightSibbesenFormAat),
                                          as_row(umcaDiameterFromHeightSibbesenFormPhysio),
                                          as_row(umcaDiameterFromHeightSibbesenFormRelHt),
                                          as_row(umcaDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "height", species = "UMCA", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(umcaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(umca2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman-Richards form approximate BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umcaDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = umcaDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  #geom_line(aes(x = umcaDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAT", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
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


## collect model parameters
umcaParameters = bind_rows(bind_rows(get_coefficients(umcaHeightFromDiameterChapmanRichards),
                                     get_coefficients(umcaHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(umcaHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(umcaHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(umcaHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(umcaHeightFromDiameterCurtis),
                                     get_coefficients(umcaHeightFromDiameterHossfeld),
                                     get_coefficients(umcaHeightFromDiameterKorf),
                                     get_coefficients(umcaHeightFromDiameterLinear),
                                     get_coefficients(umcaHeightFromDiameterMichaelisMenten),
                                     get_coefficients(umcaHeightFromDiameterParabolic),
                                     get_coefficients(umcaHeightFromDiameterPower),
                                     get_coefficients(umcaHeightFromDiameterProdan),
                                     get_coefficients(umcaHeightFromDiameterRatkowsky),
                                     get_coefficients(umcaHeightFromDiameterRichards),
                                     get_coefficients(umcaHeightFromDiameterSharmaParton),
                                     get_coefficients(umcaHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(umcaHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(umcaHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(umcaHeightFromDiameterSharmaZhang),
                                     get_coefficients(umcaHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(umcaHeightFromDiameterSibbesen),
                                     get_coefficients(umcaHeightFromDiameterWeibull),
                                     get_coefficients(umcaHeightFromDiameterWeibullBal),
                                     get_coefficients(umcaHeightFromDiameterWeibullBalRelHt),
                                     #get_coefficients(umcaHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(umcaHeightFromDiameterChapmanRichardsBalGnls),
                                     #get_coefficients(umcaHeightFromDiameterSharmaPartonGnls),
                                     #get_coefficients(umcaHeightFromDiameterSharmaPartonBalGnls),
                                     #get_coefficients(umcaHeightFromDiameterSharmaZhangGnls),
                                     #get_coefficients(umcaHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(umcaHeightFromDiameterWeibullGnls),
                                     get_coefficients(umcaHeightFromDiameterWeibullBalGnls)) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(get_coefficients(umcaDiameterFromHeightChapmanRichards),
                                     get_coefficients(umcaDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(umcaDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(umcaDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(umcaDiameterFromHeightChapmanForm),
                                     get_coefficients(umcaDiameterFromHeightChapmanFormAat),
                                     get_coefficients(umcaDiameterFromHeightChapmanFormBal),
                                     get_coefficients(umcaDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(umcaDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(umcaDiameterFromHeightLinear),
                                     get_coefficients(umcaDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(umcaDiameterFromHeightNaslund),
                                     get_coefficients(umcaDiameterFromHeightParabolic),
                                     get_coefficients(umcaDiameterFromHeightPower),
                                     get_coefficients(umcaDiameterFromHeightPowerAat),
                                     get_coefficients(umcaDiameterFromHeightPowerPhysio),
                                     get_coefficients(umcaDiameterFromHeightPowerRelHt),
                                     get_coefficients(umcaDiameterFromHeightRuark),
                                     get_coefficients(umcaDiameterFromHeightSchnute),
                                     get_coefficients(umcaDiameterFromHeightSharmaParton),
                                     get_coefficients(umcaDiameterFromHeightSibbesenForm),
                                     get_coefficients(umcaDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(umcaDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(umcaDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(umcaDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "UMCA",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
#umcaBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), umca2016, start = list(a1 = 0.3, b1 = 0.0006, b2 = 2.1), weights = pmin(1/basalArea, 1E4)) # step factor with nlrob()
umcaBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), umca2016, start = list(a1 = 1.36, b1 = 0.0002, b2 = 2.06, b2p = -0.27), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor with nlrob()
umcaBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^b1, umca2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.0002, b1 = 2.00), weights = pmin(1/basalArea, 1E4)) # b1p not significant
#confint2(umcaBasalAreaFromHeightKorf, level = 0.99)
#confint_nlrob(umcaBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/umca2016$basalArea, 1E4))

umcaBasalAreaFromHeightKorf$fitted.values = predict(umcaBasalAreaFromHeightKorf, umca2016)
umcaBasalAreaFromHeightKorf$residuals = umcaBasalAreaFromHeightKorf$fitted.values - umca2016$basalArea
umcaBasalAreaFromHeightPower$fitted.values = predict(umcaBasalAreaFromHeightPower, umca2016)
umcaBasalAreaFromHeightPower$residuals = umcaBasalAreaFromHeightPower$fitted.values - umca2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(umcaBasalAreaFromHeightKorf), 100^2 * mean(umcaBasalAreaFromHeightKorf$residuals), mean(abs(umcaBasalAreaFromHeightKorf$residuals)), 1 - sum(umcaBasalAreaFromHeightKorf$residuals^2) / sum((umca2016$basalArea - mean(umca2016$basalArea)^2)),
        "power", AIC(umcaBasalAreaFromHeightPower), 100^2 * mean(umcaBasalAreaFromHeightPower$residuals), mean(abs(umcaBasalAreaFromHeightPower$residuals)), 1 - sum(umcaBasalAreaFromHeightPower$residuals^2) / sum((umca2016$basalArea - mean(umca2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(umca2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = umcaBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = umcaBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "Oregon myrtle height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
