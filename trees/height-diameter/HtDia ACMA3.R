# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## bigleaf maple height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Ratkowsky, Sharma-Parton, Hossfeld, Weibull, Chapman-Richards
#acmaHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + a1*DBH / (a2 + a2p * isPlantation + DBH), acma2016, start = list(a1 = 39.5, a2 = 47.8, a2p = -9.11), weights = pmin(DBH^-2, 1))
#acmaHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), acma2016, start = list(Ha = 23.9, Hap = -1.1, d = 0.67, dp = -0.031, kU = 0.023, kUp = 0.008), weights = pmin(DBH^-2, 1))
acma2016physio = acma2016 %>% filter(is.na(elevation) == FALSE)
acma2016plantationPhysio = acma2016physio %>% filter(isPlantation)
acmaHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^b2, acma2016, start = list(a1 = 26.8, a1p = 4.20, b1 = -0.026, b2 = 0.927), weights = pmin(DBH^-2, 1)) # b1p, b2p not significant
acmaHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016physio, start = list(a1 = 33.2, a2 = 0.156, a3 = -0.009, a4 = -0.130, a5 = 1.872, a6 = 0.276, b1 = -0.024, b2 = 0.988, b2p = -0.129), weights = pmin(DBH^-2, 1)) # a1p, a2p, a3, a4, a6, b1p not significant
acmaHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = -2.4, a1p = 1.51, a2 = -0.023, a2p = 0.159, a3 = 0.045, a4 = 59.4, a4p = -25.9, b1 = -0.033, b2 = 0.018, b2p = 0.406), weights = pmin(DBH^-2, 1)) # a2, a3p not significant, NaN-inf with b1p
acmaHeightFromDiameterChapmanRichardsPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 68.5, a2 = -0.0045, a3 = -8.09, a4 = 0.783, a5 = 0.766, a6 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^-2, 1)) # a1p, a2, a4, b1p not significant
acmaHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^b1, acma2016, start = list(a1 = 1.24, a1p = 0.23, b1 = 0.28), weights = pmin(DBH^-2, 1)) # b1p not significant
acmaHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1 * DBH^b2), acma2016, start = list(a1 = 40.1, a1p = 6.30, b1 = 44.9, b2 = -0.97), weights = pmin(DBH^-2, 1)) # b1p, b2p not significant
acmaHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), acma2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^-2, 1))
acmaHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), acma2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
acmaHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), weights = pmin(DBH^-2, 1)) # a1p, b1p not significant
acmaHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), acma2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
acmaHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1 * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), acma2016, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19), weights = pmin(DBH^-2, 1)) # a1p, a3p not significant
acmaHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, acma2016, start = list(a1 = 1.11, a1p = 0.21, b1 = 0.75), weights = pmin(DBH^-2, 1)) # b1p not significant
acmaHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1/(DBH + b2 + b2p * isPlantation)), acma2016, start = list(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62), weights = pmin(DBH^-2, 1)) # b1p not significant
acmaHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), acma2016, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064), weights = pmin(DBH^-2, 1)) # Hap, dp not significant
acmaHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 5.86, a1p = 39.2, a2 = 0.36, a2p = -0.502, b1 = -0.080, b2 = -0.398, b3 = 0.979, b3p = -0.116), weights = pmin(DBH^-2, 1)) # b1p, b2p not significant
acmaHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 46.0, a2 = -0.12, b1 = -0.054, b2 = -0.375, b3 = 0.960, b3p = -0.102), weights = pmin(DBH^-2, 1)) # a1p, a2p, b1p, b2p not significant
acmaHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), acma2016physio, start = list(a1 = 7.67, a1p = 48.0, a2 = 0.34, a2p = -0.48, a3 = -0.0003, a4 = 0.041, a5 = 0.0067, a6 = -0.0018, b1 = -0.056, b2 = -0.402, b3 = 0.948, b3p = -0.096), weights = pmin(DBH^-2, 1)) # a2, a5, a6, b1p, b2p not significant
acmaHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*(tph/(standBasalAreaPerHectare))^b2*DBH))^(b3 + b3p * isPlantation), acma2016physio, start = list(a1 = 5.21, a1p = 41.5, a2 = 0.416, a2p = -0.540, a3 = -0.0003, a4 = 0.042, a5 = 0.0055, a6 = -0.0019, b1 = -0.083, b2 = -0.419, b3 = 0.968, b3p = -0.122), weights = pmin(DBH^-2, 1)) # a2, a5, a6, b1p, b2p not significant
acmaHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 10.8, a1p = 1.28, a2 = 0.252, b1 = -0.121, b2 = -0.247, b3 = 0.952, b3p = -0.089), weights = pmin(DBH^-2, 1)) # a2p, b1p, b2p not significant
acmaHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 14.1, a2 = 0.189, a3 = 0.003, b1 = -0.137, b2 = -0.283, b3 = 0.966, b3p = -0.137), weights = pmin(DBH^-2, 1)) # a1p, a2p, a3, a3p, b1p, b2p not significant
acmaHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^b2), acma2016, start = list(a1 = 0.752, a1p = 0.120, b1 = 1.180, b2 = -0.087), weights = pmin(DBH^-2, 1)) # b1p, b2p not significant
acmaHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^b2)), acma2016, start = list(a1 = 27.3, a1p = 4.27, b1 = -0.033, b2 = 0.94), weights = pmin(DBH^-2, 1)) # b1p, b2p not significant
acmaHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = 26.6, a2 = 0.095, a3 = 0.073, b1 = -0.0267, b1p = -0.0092, b2 = 0.930), weights = pmin(DBH^-2, 1)) # a1p, a2p, a3p, b2p not significant
acmaHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1  + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * (1 + (b2 + b2p * isPlantation) * pmin(relativeHeight, 1.25)) * DBH^b3)), acma2016, start = list(a1 = -1.78, a1p = 17.5, a2 = 0.043, a3 = -0.020, a3p = 1.11, a4 = 64.0, a4p = -62.3, b1 = -1.926, b1p = 1.921, b2 = -1.169, b2p = 70.9, b3 = 0.253), weights = pmin(DBH^-2, 1)) # a2p, a3 not significant
#confint2(acmaHeightFromDiameterWeibullBalRelHt, level = 0.99)

acmaHeightFromDiameterChapmanRichards = get_height_error(acmaHeightFromDiameterChapmanRichards, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsBal = get_height_error(acmaHeightFromDiameterChapmanRichardsBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(acmaHeightFromDiameterChapmanRichardsBalPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(acmaHeightFromDiameterChapmanRichardsBalRelHt, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterChapmanRichardsPhysio = get_height_error(acmaHeightFromDiameterChapmanRichardsPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterCurtis = get_height_error(acmaHeightFromDiameterCurtis, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterHossfeld = get_height_error(acmaHeightFromDiameterHossfeld, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterKorf = get_height_error(acmaHeightFromDiameterKorf, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterLinear = get_height_error(acmaHeightFromDiameterLinear, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterMichaelisMenten = get_height_error(acmaHeightFromDiameterMichaelisMenten, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterParabolic = get_height_error(acmaHeightFromDiameterParabolic, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterProdan = get_height_error(acmaHeightFromDiameterProdan, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterPower = get_height_error(acmaHeightFromDiameterPower, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterRatkowsky = get_height_error(acmaHeightFromDiameterRatkowsky, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterRichards = get_height_error(acmaHeightFromDiameterRichards, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaParton = get_height_error(acmaHeightFromDiameterSharmaParton, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBal = get_height_error(acmaHeightFromDiameterSharmaPartonBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBalPhysio = get_height_error(acmaHeightFromDiameterSharmaPartonBalPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterSharmaPartonPhysio = get_height_error(acmaHeightFromDiameterSharmaPartonPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaHeightFromDiameterSharmaZhang = get_height_error(acmaHeightFromDiameterSharmaZhang, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaZhangBal = get_height_error(acmaHeightFromDiameterSharmaZhangBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSibbesen = get_height_error(acmaHeightFromDiameterSibbesen, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibull = get_height_error(acmaHeightFromDiameterWeibull, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBal = get_height_error(acmaHeightFromDiameterWeibullBal, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBalRelHt = get_height_error(acmaHeightFromDiameterWeibullBalRelHt, acma2016, acma2016natural, acma2016plantation)

acmaHeightFromDiameterResults = bind_rows(as_row("Chapman-Richards", acmaHeightFromDiameterChapmanRichards),
                                          as_row("Chapman-Richards BAL", acmaHeightFromDiameterChapmanRichardsBal),
                                          as_row("Chapman-Richards BAL physio", acmaHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row("Chapman-Richards BAL RelHt", acmaHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row("Chapman-Richards physio", acmaHeightFromDiameterChapmanRichardsPhysio),
                                          as_row("Curtis", acmaHeightFromDiameterCurtis),
                                          as_row("Hossfeld", acmaHeightFromDiameterHossfeld),
                                          as_row("Korf", acmaHeightFromDiameterKorf),
                                          as_row("linear", acmaHeightFromDiameterLinear),
                                          as_row("generalized Michaelis-Menten", acmaHeightFromDiameterMichaelisMenten),
                                          as_row("parabolic", acmaHeightFromDiameterParabolic),
                                          as_row("power", acmaHeightFromDiameterPower),
                                          as_row("Prodan", acmaHeightFromDiameterProdan),
                                          as_row("Ratkowsky", acmaHeightFromDiameterRatkowsky),
                                          as_row("unified Richards", acmaHeightFromDiameterRichards),
                                          as_row("Sharma-Parton", acmaHeightFromDiameterSharmaParton),
                                          as_row("Sharma-Parton BAL", acmaHeightFromDiameterSharmaPartonBal),
                                          as_row("Sharma-Parton BAL physio", acmaHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row("Sharma-Parton physio", acmaHeightFromDiameterSharmaPartonPhysio),
                                          as_row("Sharma-Zhang", acmaHeightFromDiameterSharmaZhang),
                                          as_row("Sharma-Zhang BAL", acmaHeightFromDiameterSharmaZhangBal),
                                          as_row("Sibbesen", acmaHeightFromDiameterSibbesen),
                                          as_row("Weibull", acmaHeightFromDiameterWeibull),
                                          as_row("Weibull BAL", acmaHeightFromDiameterWeibullBal),
                                          as_row("Weibull BAL RelHt", acmaHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "ACMA3", deltaAic = aic - min(aic)) %>%
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
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = acma2016$isPlantation)) +
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
#acmaHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
##acmaHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), acma2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#acmaHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), acma2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
##save(acmaHeightFromDiameterChapmanRichardsGnls, acmaHeightFromDiameterChapmanRichardsBalGnls, acmaHeightFromDiameterSharmaPartonGnls, acmaHeightFromDiameterSharmaPartonBalGnls, acmaHeightFromDiameterSharmaZhangGnls, acmaHeightFromDiameterSharmaZhangBalGnls, acmaHeightFromDiameterWeibullGnls, acmaHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia ACMA3 GNLS.rdata")
#save(acmaHeightFromDiameterChapmanRichardsGnls, acmaHeightFromDiameterChapmanRichardsBalGnls, acmaHeightFromDiameterSharmaPartonGnls, acmaHeightFromDiameterSharmaPartonBalGnls, acmaHeightFromDiameterSharmaZhangGnls, acmaHeightFromDiameterWeibullGnls, acmaHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia ACMA3 GNLS.rdata")
load("trees/height-diameter/HtDia ACMA3 GNLS.rdata")
acmaHeightFromDiameterWeibullGnls = acmaHeightFromDiameterWykoffGnls # temporary naming error fixup
acmaHeightFromDiameterWeibullBalGnls = acmaHeightFromDiameterWykoffBalGnls

acmaHeightFromDiameterChapmanRichardsGnls = get_height_error(acmaHeightFromDiameterChapmanRichardsGnls, acma2016, acma2016natural, acma2016plantation)
#acmaHeightFromDiameterChapmanRichardsBalGnls = get_height_error(acmaHeightFromDiameterChapmanRichardsBalGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonGnls = get_height_error(acmaHeightFromDiameterSharmaPartonGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaPartonBalGnls = get_height_error(acmaHeightFromDiameterSharmaPartonBalGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterSharmaZhangGnls = get_height_error(acmaHeightFromDiameterSharmaZhangGnls, acma2016, acma2016natural, acma2016plantation)
#acmaHeightFromDiameterSharmaZhangBalGnls = get_height_error(acmaHeightFromDiameterSharmaZhangBalGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullGnls = get_height_error(acmaHeightFromDiameterWeibullGnls, acma2016, acma2016natural, acma2016plantation)
acmaHeightFromDiameterWeibullBalGnls = get_height_error(acmaHeightFromDiameterWeibullBalGnls, acma2016, acma2016natural, acma2016plantation)

acmaHeightFromDiameterResultsGnls = bind_rows(as_row("Chapman-Richards GNLS", acmaHeightFromDiameterChapmanRichardsGnls),
                                              as_row("Chapman-Richards BAL GNLS", NULL),
                                              as_row("Sharma-Parton GNLS", acmaHeightFromDiameterSharmaPartonGnls),
                                              as_row("Sharma-Parton BAL GNLS", acmaHeightFromDiameterSharmaPartonBalGnls),
                                              as_row("Sharma-Zhang GNLS", acmaHeightFromDiameterSharmaZhangGnls),
                                              as_row("Sharma-Zhang BAL GNLS", NULL),
                                              as_row("Weibull GNLS", acmaHeightFromDiameterWeibullGnls),
                                              as_row("Weibull BAL GNLS", acmaHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "ACMA3", deltaAic = aic - min(aic)) %>%
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
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## bigleaf maple diameter-height regressions
#acmaDiameterFromHeightChapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, iter = 10000,
#                                                  start_lower = list(a1 = -15, b1 = -10, b2 = -1), 
#                                                  start_upper = list(a1 = 150, b1 = 100, b2 = 2), modelweights = pmin(TotalHt^-2, 0.5)) # NaN-inf
#acmaDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, acma2016, iter = 100,
#                                                   start_lower = list(a1 = -1, a2 = 0.011, b1 = -1, b2 = -1, b3 = -1), 
#                                                   start_upper = list(a1 = 100, a2 = 3, b1 = 10, b2 = 1, b3 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#acmaDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), acma2016, start = list(a1 = 1, a2 = 1, a2p = 0, b1 = 0.02, b2 = 0.26, b2p = 0, b3 = 0.9, b3p = 0), weights = pmin(TotalHt^-2, 0.5)) # signular gradient
acmaDiameterFromHeightChapmanRichards = nlrob(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), acma2016, start = list(a1 = 23.5, b1 = -0.022, b2 = 2.01, b2p = -0.17), weights = pmin(TotalHt^-2, 0.5)) # b1p not significant
acmaDiameterFromHeightChapmanRichardsAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.999999)), acma2016, start = list(a1 = 24.2, a2 = -0.005, b1 = -0.021, b2 = 2.02, b2p = -0.19), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 500))
acmaDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), acma2016physio, start = list(a1 = 3.30, a1p = 4.985, a2 = 0.00031, a3 = -0.172, a4 = -0.215, a5 = -0.839, a6 = 1.459, b1 = -0.0035, b1p = 0.00193, b2 = 1.798), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5 not significant
acmaDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), acma2016, start = list(a1 = 36.1, a2 = -2.78, b1 = -0.028, b1p = 0.009, b2 = 1.63), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p, b2p not significant
acmaDiameterFromHeightChapmanForm = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 540000, b1 = 0.00001, b2 = 1.11), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls() at multiple nls_multstart() positions, NaN-inf from nlrob()
acmaDiameterFromHeightChapmanFormAal = gsl_nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 46000, a2 = 159, b1 = 0.000013, b2 = 1.11), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls() and nlrob()
acmaDiameterFromHeightChapmanFormBal = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 23600, a2 = -1955, b1 = 0.00001, b2 = 1.07), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # NaN-inf from nls(), step factor from nlrob()
acmaDiameterFromHeightChapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 261000, a2 = -7700, a3 = 4881, a4 = -216000, b1 = 0.00001, b2 = 1.06), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # step factor from nls()
acmaDiameterFromHeightChapmanFormRelHt = gsl_nls(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 241000, a2 = -114000, b1 = 0.000007, b2 = 1.23), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # step factor from nls() and nlrob()
acmaDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), acma2016, weights = TotalHt^-2)
acmaDiameterFromHeightMichaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), acma2016, start = list(a1 = 100, a2 = 100, b1 = 1), weights = TotalHt^-2) # collapses to linear
acmaDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), acma2016, start = list(a1 = 5.1, a1p = -2.0, a2 = -0.12, a2p = -0.023), weights = TotalHt^-2)
acmaDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), acma2016, weights = TotalHt^-2) # (TotalHt - 1.37)^2 not significant (p = 0.058)
acmaDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282), weights = pmin(TotalHt^-2, 0.5))
acmaDiameterFromHeightPowerAal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272), weights = pmin(TotalHt^-2, 0.5)) # a2p not significant
acmaDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016physio, start = list(a1 = 4.38, a1p = -1.88, a2 = 0.00014, a3 = -1.61, a4 = -0.021, a5 = -0.112, a6 = 0.0096, b1 = 0.895, b1p = 0.189), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5, a6 not significant
acmaDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, acma2016, start = list(a1 = 2.37, a1p = -0.887, a2 = -1.508, a2p = 1513, b1 = 1.123), weights = pmin(TotalHt^-2, 0.5)) 
acmaDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
acmaDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 161), weights = pmin(TotalHt^-2, 0.5)) # step factor with nlrob()
acmaDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, acma2016, start = list(a1 = 0.0001, a2 = 3.8, b1 = 0.0174, b2 = 0.112, b3 = -2.4), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf from nls() at nls_multstart() positions, NaN-inf with nlrob()
acmaDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
acmaDiameterFromHeightSibbesenFormAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
acmaDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), acma2016physio, start = list(a1 = 0.792, a1p = -0.272, a2 = 0.00003, a3 = -0.345, a4 = -0.0044, a5 = -0.020, a6 = 0.0021, b1 = 2.433, b2 = -0.164, b2p = 0.0277), weights = pmin(TotalHt^-2, 0.5)) # no physiographic parameters signficant
acmaDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.496, a2 = -0.188, b1 = 2.31, b2 = -0.12), weights = pmin(TotalHt^-2, 0.5)) # a2 not significant
acmaDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, acma2016, start = list(a1 = -185000, b1 = 0.00001, b2 = 1.11), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 50)) # missing value with nlrob()
#confint2(acmaDiameterFromHeightWeibull, level = 0.99)

acmaDiameterFromHeightChapmanRichards = get_dbh_error(acmaDiameterFromHeightChapmanRichards, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanRichardsAal = get_dbh_error(acmaDiameterFromHeightChapmanRichardsAal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(acmaDiameterFromHeightChapmanRichardsPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(acmaDiameterFromHeightChapmanRichardsRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanForm = get_dbh_error(acmaDiameterFromHeightChapmanForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormAal = get_dbh_error(acmaDiameterFromHeightChapmanFormAal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormBal = get_dbh_error(acmaDiameterFromHeightChapmanFormBal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(acmaDiameterFromHeightChapmanFormBalRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightChapmanFormRelHt = get_dbh_error(acmaDiameterFromHeightChapmanFormRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightLinear = get_dbh_error(acmaDiameterFromHeightLinear, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightMichaelisMentenForm = get_dbh_error(acmaDiameterFromHeightMichaelisMentenForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightNaslund = get_dbh_error(acmaDiameterFromHeightNaslund, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightParabolic = get_dbh_error(acmaDiameterFromHeightParabolic, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPower = get_dbh_error(acmaDiameterFromHeightPower, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerAal = get_dbh_error(acmaDiameterFromHeightPowerAal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerPhysio = get_dbh_error(acmaDiameterFromHeightPowerPhysio, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightPowerRelHt = get_dbh_error(acmaDiameterFromHeightPowerRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightRuark = get_dbh_error(acmaDiameterFromHeightRuark, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSchnute = get_dbh_error(acmaDiameterFromHeightSchnute, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSharmaParton = get_dbh_error(acmaDiameterFromHeightSharmaParton, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenForm = get_dbh_error(acmaDiameterFromHeightSibbesenForm, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenFormAal = get_dbh_error(acmaDiameterFromHeightSibbesenFormAal, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightSibbesenFormPhysio = get_dbh_error(acmaDiameterFromHeightSibbesenFormPhysio, acma2016physio, acma2016natural, acma2016plantationPhysio)
acmaDiameterFromHeightSibbesenFormRelHt = get_dbh_error(acmaDiameterFromHeightSibbesenFormRelHt, acma2016, acma2016natural, acma2016plantation)
acmaDiameterFromHeightWeibull = get_dbh_error(acmaDiameterFromHeightWeibull, acma2016, acma2016natural, acma2016plantation)

acmaDiameterFromHeightResults = bind_rows(as_row("Chapman-Richards", acmaDiameterFromHeightChapmanRichards),
                                          as_row("Chapman-Richards AAL", acmaDiameterFromHeightChapmanRichardsAal), # not AIC supported
                                          as_row("Chapman-Richards physio", acmaDiameterFromHeightChapmanRichardsPhysio),
                                          as_row("Chapman-Richards RelHt", acmaDiameterFromHeightChapmanRichardsRelHt),
                                          as_row("Chapman-Richards form", acmaDiameterFromHeightChapmanForm),
                                          as_row("Chapman-Richards form AAL", acmaDiameterFromHeightChapmanFormAal),
                                          as_row("Chapman-Richards form BAL", acmaDiameterFromHeightChapmanFormBal),
                                          as_row("Chapman-Richards form BAL RelHt", acmaDiameterFromHeightChapmanFormBalRelHt),
                                          as_row("Chapman-Richards form RelHt", acmaDiameterFromHeightChapmanFormRelHt),
                                          as_row("linear", acmaDiameterFromHeightLinear),
                                          as_row("generalized Michaelis-Menten form", acmaDiameterFromHeightMichaelisMentenForm),
                                          as_row("Näslund", acmaDiameterFromHeightNaslund),
                                          as_row("parabolic", acmaDiameterFromHeightParabolic),
                                          as_row("power", acmaDiameterFromHeightPower),
                                          as_row("power AAL", acmaDiameterFromHeightPowerAal), # not AIC supported
                                          as_row("power physio", acmaDiameterFromHeightPowerPhysio),
                                          as_row("power RelHt", acmaDiameterFromHeightPowerRelHt),
                                          as_row("Ruark", acmaDiameterFromHeightRuark),
                                          as_row("Schnute", acmaDiameterFromHeightSchnute),
                                          as_row("modified Sharma-Parton", acmaDiameterFromHeightSharmaParton),
                                          as_row("Sibbesen form", acmaDiameterFromHeightSibbesenForm),
                                          as_row("Sibbesen form AAL", acmaDiameterFromHeightSibbesenFormAal),
                                          as_row("Sibbesen form physio", acmaDiameterFromHeightSibbesenFormPhysio),
                                          as_row("Sibbesen form RelHt", acmaDiameterFromHeightSibbesenFormRelHt),
                                          as_row("Weibull", acmaDiameterFromHeightWeibull)) %>%
  mutate(responseVariable = "height", species = "ACMA3", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(acmaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(acma2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = acmaDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman-Richards approximate BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = acmaDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards BAL", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = acmaDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "generalized Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  geom_line(aes(x = acmaDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = acmaDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schute", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  #geom_line(aes(x = acmaDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards form"), na.rm = TRUE) +
  #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 41, label = "bigleaf maple, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
acmaParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", acmaHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", acmaHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", acmaHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", acmaHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", acmaHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", acmaHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", acmaHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", acmaHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", acmaHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten", acmaHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", acmaHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", acmaHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", acmaHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", acmaHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", acmaHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", acmaHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", acmaHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", acmaHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", acmaHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", acmaHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", acmaHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", acmaHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", acmaHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", acmaHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", acmaHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards GNLS", acmaHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     #bind_rows(c(method = "Chapman-Richards BAL GNLS", acmaHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton GNLS", acmaHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton BAL GNLS", acmaHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang GNLS", acmaHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     #bind_rows(c(method = "Sharma-Zhang BAL GNLS", acmaHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", acmaHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", acmaHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", acmaDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", acmaDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", acmaDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", acmaDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form", acmaDiameterFromHeightChapmanForm$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form AAL", acmaDiameterFromHeightChapmanFormAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL", acmaDiameterFromHeightChapmanFormBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL RelHt", acmaDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form RelHt", acmaDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                     bind_rows(c(method = "linear", acmaDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten form", acmaDiameterFromHeightMichaelisMentenForm$m$getPars())),
                                     bind_rows(c(method = "Näslund", acmaDiameterFromHeightNaslund$m$getPars())),
                                     bind_rows(c(method = "parabolic", acmaDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", acmaDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", acmaDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", acmaDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", acmaDiameterFromHeightPowerRelHt$m$getPars())),
                                     bind_rows(c(method = "Ruark", acmaDiameterFromHeightRuark$m$getPars())),
                                     bind_rows(c(method = "Schnute", acmaDiameterFromHeightSchnute$m$getPars())),
                                     bind_rows(c(method = "modified Sharma-Parton", acmaDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", acmaDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", acmaDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", acmaDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", acmaDiameterFromHeightSibbesenFormRelHt$m$getPars())),
                                     bind_rows(c(method = "Weibull", acmaDiameterFromHeightWeibull$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "ACMA3",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


