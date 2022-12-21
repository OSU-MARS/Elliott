# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Oregon myrtle (California bay) height-diameter regression form sweep
#umcaHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), umca2016, start = list(Ha = 14.4, Hap = -3.6, d = 2.78, dp = -0.99, kU = 0.054, kUp = 0.026))
umca2016 = trees2016 %>% filter(Species == "OM", isLiveUnbroken, TotalHt > 0) %>% # live Oregon myrtles measured for height
  mutate(dbhWeight = pmin(DBH^-0.8, 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.6, -1.1), 0.5))
umca2016physio = umca2016 %>% filter(is.na(elevation) == FALSE)
umca2016gamConstraint = c(DBH = -0.8547/0.8790, TotalHt = 1.37, standBasalAreaPerHectare = median(umca2016$standBasalAreaPerHectare), basalAreaLarger = median(umca2016$basalAreaLarger), standBasalAreaApprox = median(umca2016$standBasalAreaApprox), tallerApproxBasalArea = median(umca2016$tallerApproxBasalArea), elevation = median(umca2016physio$elevation), slope = median(umca2016physio$slope), aspect = median(umca2016physio$aspect), topographicShelterIndex = median(umca2016physio$topographicShelterIndex), relativeHeight = median(umca2016$relativeHeight)) # point constraint for mgcv::s()
#umca2016natural = umca2016 %>% filter(isPlantation == FALSE)
#umca2016plantation = umca2016 %>% filter(isPlantation)
#umca2016plantationPhysio = umca2016physio %>% filter(isPlantation)

umcaHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = 16.5, b1 = -0.064, b2 = 1.270, b2p = -0.177))) # a1p, b1p not significant
umcaHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, umca2016, start = list(a1 = 19.0, a2 = 0.201, a3 = -0.152, b1 = -0.049, b1p = -0.0105, b2 = 1.156)) # a2, a2p, a3, a3p, b2p not significant
umcaHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation + a7 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 8.62, a2 = 0.050, a4 = -0.014, a7 = -1.113, b1 = -0.062, b2 = 1.265, b2p = -0.219)) # a1p, a2p, a3, a5, a6, a7, a8, b1p not significant
umcaHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = list(a1 = -1.48, a2 = 0.030, a2p = -0.131, a3 = -0.009, a3p = 0.192, a4 = 56.1, a4p = -32.7, b1 = -0.278, b2 = 0.389, b2p = 0.684)) # a1p not significant
umcaHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a6 * cos(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016physio, start = list(a1 = 8.58, a4 = -0.018, a6 = -1.117, b1 = -0.064, b2 = 1.261, b2p = -0.193)) # a5, a7, a8, b1p not significant
umcaHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^b1, umca2016, start = list(a1 = 0.909, a1p = 0.221, b1 = 0.219)) # b1p not significant
umcaHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 11, pc = umca2016gamConstraint), data = umca2016)
umcaHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15, pc = umca2016gamConstraint), data = umca2016)
umcaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, basalAreaLarger, elevation, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 36, pc = umca2016gamConstraint), data = umca2016physio)
umcaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 35, pc = umca2016gamConstraint), data = umca2016physio)
umcaHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), umca2016, start = list(a1 = 21.4, a1p = -4.37, b1 = 43.8, b1p = -19.9, b2 = -1.27)) # b2p not significant
umcaHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)*DBH^b2), umca2016, start = list(a1 = 49.8, b1 = -5.02, b1p = 0.404, b2 = -0.386)) # a1p, b2p not significant
umcaHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), umca2016)
umcaHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27)) # b1p not significant
umcaHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), umca2016)
umcaHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), umca2016, start = list(a1 = 0.037, a1p = 0.020, a2 = 1.109, a2p = -0.472, a3 = 1.242)) # a3p not significant
umcaHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, umca2016, start = list(a1 = 0.821, a1p = 0.206, b1 = 0.810)) # b1p not significant
umcaHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + a1*exp((b1 + b1p * isPlantation)/(DBH + b2)), umca2016, start = list(a1 = 21.6, b1 = -16.5, b1p = 1.974, b2 = 3.629)) # a1p, b2p not significant
umcaHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap * isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), umca2016, start = list(Ha = 14.6, Hap = -4.146, d = 2.198, kU = 0.048, kUp = 0.045)) # dp not significant
umcaHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, umca2016, start = list(a1 = 5.76, b1 = 0.261, b2 = -0.049, b2p = -0.024, b3 = 0.085, b4 = 1.234)) # a1p, b1p, b3p, b4p not significant
umcaHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 9.074, b1 = 0.149, b2 = -0.065, b3 = 0.013, b4 = 1.299, b4p = -0.254)) # a1p, b1, b1p, b2p, b3p not significant
umcaHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 9.30, a4 = -0.0007, b1 = 0.151, b2 = -0.061, b3 = 0.061, b4 = 1.313, b4p = -0.253)) # a5, a6, a7, a8, b1p, b2p, b3p not significant
umcaHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^(b4 + b4p * isPlantation), umca2016physio, start = list(a1 = 13.8, a4 = -0.01, b1 = 0.08, b2 = -0.05, b3 = 0.011, b4 = 1.2, b4p = -0.23)) # a5, a6, a7, a8, b1p, b2p, b3, b3p not significant
umcaHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), umca2016, start = list(a1 = 13.6, b1 = 0.054, b2 = -0.030, b3 = 0.125, b4 = 1.295, b4p = -0.232)) # a1p, b1, b1p, b2p, b3p not significant
umcaHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, umca2016, start = list(a1 = 8.742, a2 = -0.0001, b1 = 0.174, b2 = -0.039, b3 = 0.070, b3p = 0.090, b4 = 1.258), significant = FALSE) # a1p, a2, a2p, b1, b1p, b2p, b4p not significant
umcaHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), umca2016, start = list(a1 = 0.346, a1p = 0.153, b1 = 1.783, b2 = -0.150, b2p = -0.033)) # b1p not significant
umcaHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = list(a1 = 16.7, a1p = -3.627, b1 = -0.0315, b1p = -0.0256, b2 = 1.174)) # b2p not significant
umcaHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = list(a1 = 18.0, a2 = 0.079, a3 = -0.037, b1 = -0.037, b1p = -0.007, b2 = 1.043)) # a1p, a2, a2p, a3p, b2p not significant
umcaHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare + a9*pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * DBH^b2)), umca2016, start = list(a1 = -2.9, a2 = 0.04, a3 = 0.04, a9 = 53, b1 = -0.55, b1p = 0.35, b2 = 0.39), maxit = 30, control = nls.control(maxiter = 500)) # a1p, a2p, a3, a3p, a4p, b2p not significant
#confint_nlrob(umcaHeightFromDiameter$sharmaZhangBal, level = 0.99, weights = pmin(umca2016$DBH^-0.8, 1))

umcaHeightFromDiameterResults = bind_rows(lapply(umcaHeightFromDiameter, as_row)) %>%
  mutate(responseVariable = "height", species = "UMCA", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
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
#umcaHeightFromDiameterGnls$chapmanRichards = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), umca2016, start = umcaHeightFromDiameter$chapmanRichards$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#umcaHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, umca2016, start = umcaHeightFromDiameter$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE))
#umcaHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, umca2016, start = umcaHeightFromDiameter$sharmaParton$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#umcaHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^(b3 + b3p * isPlantation), umca2016, start = umcaHeightFromDiameter$sharmaPartonBal$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#umcaHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), umca2016, start = umcaHeightFromDiameter$sharmaZhang$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#umcaHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, umca2016, start = umcaHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.005, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation or nlsTol = 0.001
#umcaHeightFromDiameterGnls$Weibull = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = umcaHeightFromDiameter$weibull$m$getPars(), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#umcaHeightFromDiameterGnls$WeibullBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), umca2016, start = umcaHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.40, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#save(umcaHeightFromDiameterGnls, file = "trees/height-diameter/HtDia UMCA GNLS.rdata")

load("trees/height-diameter/HtDia UMCA GNLS.rdata")
umcaHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", umcaHeightFromDiameterGnls$chapmanRichards, umca2016)
umcaHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", umcaHeightFromDiameterGnls$chapmanRichardsBal, umca2016)
umcaHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", umcaHeightFromDiameterGnls$sharmaParton, umca2016)
umcaHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", umcaHeightFromDiameterGnls$sharmaPartonBal, umca2016)
umcaHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", umcaHeightFromDiameterGnls$sharmaZhang, umca2016)
umcaHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", umcaHeightFromDiameterGnls$sharmaZhangBal, umca2016)
umcaHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", umcaHeightFromDiameterGnls$weibull, umca2016)
umcaHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", umcaHeightFromDiameterGnls$weibullBal, umca2016)

umcaHeightFromDiameterResultsGnls = bind_rows(lapply(umcaHeightFromDiameterGnls, as_row)) %>%
  mutate(responseVariable = "height", species = "UMCA", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
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
#umcaDiameterFromHeight$chapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, iter = 100,
#                                                  start_lower = list(a1 = -10, b1 = -2, b2 = -1), 
#                                                  start_upper = list(a1 = 100, b1 = 2.5, b2 = 1), modelweights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.1), 0.5))
#umcaDiameterFromHeight$chapmanFormAat = fit_gsl_nls(DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 15, a2 = 0, b1 = 0.1, b2 = 0.5), control = nls.control(maxiter = 500)) # step factor with nlrob()
#umcaDiameterFromHeight$sharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, umca2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.01, a2 = -1, b1 = -10, b2 = -0.5, b3 = 0.2),
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 10, b2 = 0.5, b3 = 1.5), modelweights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.1), 0.5))
#umcaDiameterFromHeight$sharmaParton = fit_gsl_nls(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), umca2016, start = list(a1 = 89, a2 = 0.61, a2p = 0.05, b1 = 0.0002, b2 = -0.73, b2p = -0.79, b3 = 0.34, b3p = -0.03)) # singnular gradient with nls()
umcaDiameterFromHeight = list(chapmanForm = fit_gsl_nls("Chapman-Richards form", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 50000, b1 = 0.00001, b2 = 1.034), control = gsl_nls_control(maxiter = 50))) # NaN-inf with nls() at multiple nls_multstart() points, NaN-inf with nlrob()
umcaDiameterFromHeight$chapmanFormAat = fit_gsl_nls("Chapman-Richards form AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, umca2016, start = list(a1 = 604, a2 = 0.51, b1 = 0.0004, b2 = 1.01), control = gsl_nls_control(maxiter = 300)) # NaN-inf with nls(), step factor with nlrob()
umcaDiameterFromHeight$chapmanFormBal = fit_gsl_nls("Chapman-Richards form BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 1000, a2 = -16, b1 = 0.001, b2 = 1.02), control = gsl_nls_control(maxiter = 300)) # NaN-inf with nls(), step factor with nlrob()
umcaDiameterFromHeight$chapmanFormBalRelHt = fit_gsl_nls("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 1500, a2 = -400.0, a3 = 520, a9 = -2000, b1 = 0.0006, b2 = 1.02), control = gsl_nls_control(maxiter = 250)) # step factor with nls() and nlrob()
umcaDiameterFromHeight$chapmanFormRelHt = fit_gsl_nls("Chapman-Richards form RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), umca2016, start = list(a1 = 665, a9 = -406, b1 = 0.0034, b2 = 1.12), control = gsl_nls_control(maxiter = 500)) # step factor with nls(), step factor with nlrob()
umcaDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 33.9, b1 = -0.047, b2 = 1.43, b2p = -0.15))
umcaDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 34.4, a2 = 0.004, b1 = -0.047, b2 = 1.41, b2p = -0.13))
umcaDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), umca2016physio, start = list(a1 = 58.6, a1p = 166.3, a4 = 0.0133, a5 = -52.7, a6 = 5.394, a8 = 0.602, b1 = -0.0637, b1p = 0.0584, b2 = 1.329)) # a7 not significant
umcaDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), umca2016, start = list(a1 = 39.0, a9 = -2.62, b1 = -0.043, b2 = 1.40, b2p = -0.132)) # a1p convergence questionable, b1p not significant
umcaDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = umca2016gamConstraint), data = umca2016)
umcaDiameterFromHeight$gamAat = fit_gam("REML GAM AA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 14, pc = umca2016gamConstraint), data = umca2016)
umcaDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM AA+T physio", DBH ~ s(TotalHt, standBasalAreaApprox, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 28, pc = umca2016gamConstraint), data = umca2016physio) # >5 minute Zen 3 3.4 GHz fit time with AA+T, elevation, slope, sin(aspect), and TSI: elevation dropped based on physio AICs, then TSI+aspect
umcaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 60, pc = umca2016gamConstraint), data = umca2016physio) # 72.5% deviance explained (AIC 6089): 75.7% (5995) without elevation, 71.3% (6111) without slope, 72.4% (6043) without aspect, 72.7% (6034) without topographic shelter
umcaDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), umca2016)
umcaDiameterFromHeight$michaelisMentenForm = fit_gsl_nls("Michaelis-Menten form", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), umca2016, start = list(a1 = 100, a2 = 100, b1 = 1), control = gsl_nls_control(maxiter = 75)) # collapses to linear
umcaDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), umca2016, start = list(a1 = 4.3, a1p = -1.8, a2 = -0.14, a2p = -0.038))
umcaDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), umca2016, significant = FALSE) # collapses to linear since (TotalHt - 1.37)^2 and isPlantation*(TotalHt - 1.37)^2 not significant
umcaDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 3.28, a1p = -2.10, b1 = 0.917, b1p = 0.332))
umcaDiameterFromHeight$powerAat = fit_nlrob("power AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), umca2016, start = list(a1 = 2.71, a2 = 0.00054, b1 = 0.975, b1p = -0.0696)) # a1p, a2p not significant
umcaDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, umca2016physio, start = list(a1 = 5.12, a1p = -0.67, a4 = 0.0013, a5 = -4.15, a6 = 0.46, a7 = -0.11, a8 = 0.044, b1 = 0.93)) # a1p, a2, b1p not significant
umcaDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a1p * isPlantation + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, umca2016, start = list(a1 = 2.37, a1p = -0.887, a9 = -1.508, a9p = 1513, b1 = 1.123))
umcaDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 2.48, b1 = 1.14, b1p = -0.57, b2 = -0.019, b2p = 0.092)) # a1p not significant
umcaDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), umca2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 177)) # NaN-inf with nlrob()
umcaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, umca2016, start = list(a1 = 104, b1 = 0.66, b2 = 0.0001, b3 = -1.17, b4 = 0.31)) # a2p, b2p, b3p not significant, NaN-inf with nls() at nls_multistart() point, NaN-inf or code error with nlrob()
umcaDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.43, b1 = 2.45, b2 = -0.15)) # no significant plantation effects
umcaDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.36, a2 = 0.0002, b1 = 2.59, b2 = -0.156)) # a2 not significant
umcaDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016physio, start = list(a1 = 1.53, a1p = -0.38, a4 = 0.001, a5 = -1.15, a6 = 0.14, a8 = 0.02, b1 = 1.95, b2 = -0.15)) # a7 not significant, step factor with b1p
umcaDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), umca2016, start = list(a1 = 0.496, a9 = -0.188, b1 = 2.31, b2 = -0.12)) # a4 not significant
umcaDiameterFromHeight$weibull = fit_gsl_nls("Weibull", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, umca2016, start = list(a1 = -3800, b1 = 0.0006, b2 = 1.03), control = gsl_nls_control(maxiter = 500)) # NaN-inf with nlrob()
#confint_nlrob(umcaDiameterFromHeight$weibull, level = 0.99)

umcaDiameterFromHeightResults = bind_rows(lapply(umcaDiameterFromHeight, as_row)) %>%
  mutate(responseVariable = "DBH", species = "UMCA", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(umcaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

  ggplot(umca2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman-Richards form approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman-Richards form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$schnute), y = TotalHt, color = "Schnute", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$sharmaParton), y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
    #geom_line(aes(x = predict(umcaDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AA+T", group = isPlantation), alpha = 0.5) +
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
}


## collect model parameters
umcaParameters = bind_rows(bind_rows(bind_rows(lapply(umcaHeightFromDiameter, get_coefficients)),
                                     bind_rows(lapply(umcaHeightFromDiameterGnls, get_coefficients))) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(lapply(umcaDiameterFromHeight, get_coefficients)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "UMCA") %>%
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
if (includeInvestigatory)
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