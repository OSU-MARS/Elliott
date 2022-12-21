# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## minority species height-diameter regression form sweep
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 43.0, a1p = 13.5, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875), weights = dbhWeight, control = list(maxiter = 50)) # a1p not significant
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 43.0, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875, b2p = 0), weights = dbhWeight, control = list(maxiter = 500)) # > 500 iterations
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 89.3, a2 = 1.166, a2p = 0, a3 = -0.039, a3p = 0, b1 = -0.00544, b2 = 0.873), weights = dbhWeight) # a2p, a3p not significant
#otherHeightFromDiameter$sharmaZhang = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 58.1, a1p = -47.4, a2 = 0.162, a2p = 0.032, b1 = -0.021, b1p = -0.222, b2 = -0.292, b2p = 0.036, b3 = 0.818, b3p = 0.165), weights = dbhWeight)
#otherHeightFromDiameter$sibbesen = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.714, a1p = 0.088, b1 = 1.172, b2 = -0.074, b2p = -0.0040), weights = dbhWeight) # a1p, b2p not significant
other2016 = trees2016 %>% filter((Species %in% c("DF", "RA", "WH", "BM", "OM", "RC")) == FALSE, isLiveUnbroken, TotalHt > 0) %>%
  mutate(dbhWeight = pmin(DBH^if_else(isPlantation, -1.1, -2.0), 1),
         heightWeight = pmin(TotalHt^-2.4, 0.5))
other2016physio = other2016 %>% filter(is.na(elevation) == FALSE)
other2016gamConstraint = c(DBH = -1.4359/0.8154, TotalHt = 1.37, standBasalAreaPerHectare = median(other2016$standBasalAreaPerHectare), basalAreaLarger = median(other2016$basalAreaLarger), standBasalAreaApprox = median(other2016$standBasalAreaApprox), tallerApproxBasalArea = median(other2016$tallerApproxBasalArea), elevation = median(other2016physio$elevation), slope = median(other2016physio$slope), aspect = median(other2016physio$aspect), topographicShelterIndex = median(other2016physio$topographicShelterIndex), relativeHeight = median(other2016$relativeHeight)) # point constraint for mgcv::s()
#other2016natural = other2016 %>% filter(isPlantation == FALSE)
#other2016plantation = other2016 %>% filter(isPlantation)
#other2016plantationPhysio = other2016physio %>% filter(isPlantation)
#otherConifer2016 = other2016 %>% filter(Species %in% c("XX", "CX", "SS", "PC", "PY", "GF", "LP"))
#otherHardwood2016 = other2016 %>% filter(Species %in% c("CA", "HX", "CH", "PM", "GC", "PD", "TO", "WI", "OA", "WO"))

otherHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016, start = list(a1 = 24.7, b1 = -0.033, b1p = -0.003, b2 = 1.000))) # a1p, b2p not significant
otherHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 81.9, a2 = 1.97, a3 = -1.23, b1 = -0.006, b2 = 0.881), control = list(maxiter = 50)) # a1p, a2p, a2, a3p, a3, b1p, b2p not significant
otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, other2016physio, start = list(a1 = 28.5, a2 = 0.286, a4 = -0.011, a8 = -0.176, b1 = -0.022, b2 = 0.908), maxit = 50) # a1p, a2p, a3, a5, a6, a7, b1p, b2p not significant
otherHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = -1.16, a2 = -0.05, a3 = 0.07, a3p = 0.09, a4 = 16.6, a4p = -6.7, b1 = -0.05, b2 = 0.28, b2p = 0.131)) # a2, a2p, a3, b1p not significant
otherHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016physio, start = list(a1 = 27.0, a4 = -0.015, b1 = -0.037, b1p = -0.006, b2 = 0.999)) # a1p, a4p, a5, a6, a7, a8, b2p not significant
otherHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1*DBH / (1 + DBH)^b1, other2016, start = list(a1 = 1.086, b1 = 0.190)) # a1p, b1p not significant
otherHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 7, pc = other2016gamConstraint), data = other2016)
otherHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15, pc = other2016gamConstraint), data = other2016) # slow with unrestricted k
otherHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 23, pc = other2016gamConstraint), data = other2016physio) # unrestricted k tensor product (te() + te()) impractically slow (>2 h, Zen 3 @ 3.4 GHz)
otherHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 17, pc = other2016gamConstraint), data = other2016physio)
otherHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), other2016, start = list(a1 = 34.6, b1 = 40.2, b2 = -1.03)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), other2016, start = list(a1 = 381, b1 = -6.26, b2 = -0.20), control = list(maxiter = 50)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), other2016)
otherHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 34.6, a2 = 40.2, b1 = 1.03)) # a1p, a2p, b1p not significant
otherHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), other2016)
otherHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), other2016, start = list(a1 = 0.028, a2 = 1.08, a3 = 0.15)) # a1p, a2p, a3p not significant
otherHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, other2016, start = list(a1 = 0.99, b1 = 0.84)) # a1p, b1p not significant
otherHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), other2016, start = list(a1 = 24.9, b1 = -17.8, b2 = 4.65)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), other2016, start = list(Ha = 12.9, d = 2.59, kU = 0.064, kUp = 0.008)) # Hap, dp not significant
otherHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, other2016, maxit = 20, start = list(a1 = 160, b1 = -0.45, b2 = -0.0065, b3 = -0.41, b4 = 0.87)) # a1p, b1p, b2p, b3p, b4p not significant
otherHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, other2016, start = list(a1 = 90.7, b1 = -0.328, b2 = -0.027, b3 = -0.277, b4 = 0.834), maxit = 20) # a1p, b1p, b2p, b3p, b4p not significant
otherHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, other2016physio, start = list(a1 = 108, a4 = -0.05, b1 = -0.40, b2 = -0.08, b3 = -0.40, b4 = 0.83), maxit = 30, control = nls.control(maxiter = 50)) # a4, a5, a6, a7, a8, b1, b2p, b3p, b4p not significant, step factor with fit_nlrob()
otherHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, other2016physio, start = list(a1 = 108, a4 = -0.05, b1 = -0.40, b2 = -0.08, b3 = -0.40, b4 = 0.83), maxit = 20) # a1p, a4, a5, a6, a7, a8, b1p, b2p, b3p, b4p not significant but a4 significant if non-significant a8 is present, step factor with a4p
otherHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), other2016, start = list(a1 = 14.6, b1 = 0.19, b2 = -0.13, b3 = -0.27, b4 = 0.96, b4p = -0.09)) # a1p, b1p, b2p, b3p not significant
otherHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, other2016, start = list(a1 = 62.2, a2 = 0.007, b1 = -0.06, b2 = -0.181, b3 = -0.333, b3p = 0.044, b4 = 0.912), maxit = 50, control = list(maxiter = 50)) # a1p, a2, a2p, b1, b1p, b2p, b4p not significant
otherHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), other2016, start = list(a1 = 0.742, b1 = 1.207, b2 = -0.088)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 24.7, b1 = -0.035, b2 = 0.96, b2p = 0.07)) # a1p, b2p not significant
otherHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = 42.3, a2 = 0.89, a3 = -0.25, b1 = -0.019, b2 = 0.774, b2p = 0.130)) # a1p, a2, a2p, a3, a3p, b1p not significant
otherHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = list(a1 = -1.1, a2 = 0.04, a3 = 0.19, a4 = 25.2, b1 = -0.02, b2 = 0.67, b2p = 0.01), control = list(maxiter = 50)) # a2, a3p, a4p, b2p not significant, a4 evaporates for plantation weight powers smaller than 1.6
#confint_nlrob(otherHeightFromDiameter$sharmaPartonBalPhysio, level = 0.99, weights = pmin(other2016physio$DBH^if_else(other2016physio$isPlantation, -1.0, -1.9), 1))

otherHeightFromDiameterResults = bind_rows(lapply(otherHeightFromDiameter, as_row)) %>%
  mutate(responseVariable = "height", species = "other", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(otherHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = other2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = other2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$curtis), color = "Curtis", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$korf), color = "Korf", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$linear), color = "linear", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$parabolic), color = "parabolic", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$power), color = "power", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$prodan), color = "Prodan", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$richards), color = "unified Richards", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$sibbesen), color = "Sibbesen", group = other2016$isPlantation)) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$weibull), color = "Weibull", group = other2016$isPlantation)) +
    annotate("text", x = 0, y = 70, label = "merged minority species, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 70)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}

  
## other height-diameter GNLS regressions
##otherHeightFromDiameterGnls$chapmanRichards = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3*basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = otherHeightFromDiameter$sharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 1, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 1
#otherHeightFromDiameterGnls$chapmanRichards = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichards$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, logLik() NaN without
#otherHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichardsBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, other2016, start = otherHeightFromDiameter$sharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, other2016, start = otherHeightFromDiameter$sharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, maxIter = 250, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001
#otherHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), other2016, start = otherHeightFromDiameter$sharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001
##otherHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3*basalAreaLarger) * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, other2016, start = otherHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, logLik() NaN without
#otherHeightFromDiameterGnls$weibull = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameter$weibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE))
#otherHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameter$weibullBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot level correlation, with nlsTol = 0.001 without plot level correlation
#save(otherHeightFromDiameterGnls, file = "trees/height-diameter/HtDia other GNLS.rdata")

load("trees/height-diameter/HtDia other GNLS.rdata")
#otherHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", otherHeightFromDiameterGnls$chapmanRichards, other2016)
otherHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", otherHeightFromDiameterGnls$chapmanRichardsBal, other2016)
otherHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", otherHeightFromDiameterGnls$sharmaParton, other2016)
otherHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", otherHeightFromDiameterGnls$sharmaPartonBal, other2016)
otherHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", otherHeightFromDiameterGnls$sharmaZhang, other2016)
#otherHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", otherHeightFromDiameterGnls$sharmaZhangBal, other2016)
otherHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", otherHeightFromDiameterGnls$weibull, other2016)
otherHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", otherHeightFromDiameterGnls$weibullBal, other2016)

otherHeightFromDiameterResultsGnls = bind_rows(bind_rows(lapply(otherHeightFromDiameterGnls, as_row)),
                                               as_row(name = "Chapman-Richards GNLS"),
                                               as_row(name = "Sharma-Zhang BA+L GNLS")) %>%
  mutate(responseVariable = "height", species = "other", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))


if (includeInvestigatory)
{
  otherHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(otherHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(otherHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(otherHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
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
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$weibullBAL), color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = other2016natural$DBH, y = predict(otherHeightFromDiameter$weibullBALnatural), color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = other2016plantation$DBH, y = predict(otherHeightFromDiameter$weibullBALplantation), color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$Base), color = "base")) +
    geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameter$weibull), color = "ElliottWeibull")) +
    annotate("text", x = 0, y = 85, label = "a) minority species, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
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
}


## other species diameter-height regressions
#otherDiameterFromHeight$chapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, iter = 1000,
#                                              start_lower = list(a1 = 0.1, b1 = -1, b2 = -0.5), 
#                                              start_upper = list(a1 = 150, b1 = 1, b2 = 0.5), modelweights = pmin(TotalHt^-3.0, 0.5))
#otherDiameterFromHeight$chapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016, iter = 100,
#                                                             start_lower = list(a1 = 1, a1p = -100, a2 = -0.1, a3 = -100, a4 = -1, a5 = -1, a6 = -1, b1 = -0.1, b2 = -2), 
#                                                             start_upper = list(a1 = 200, a1p = 100, a2 = 0.1, a3 = 100, a4 = 10, a5 = 1, a6 = 1, b1 = 0.1, b2 = 2), modelweights = pmin(TotalHt^-3.0, 0.5))
#otherDiameterFromHeight$powerAat = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.0, a1p = -0.38, a2 = -0.013, b1 = 1.06), weights = pmin(TotalHt^-3.0, 0.5)) # a2p, b2p not significant, fit_nlrob() failure to converge >500 steps
#otherDiameterFromHeight$ruark = gsl_nls(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 1.56, a1p = 0.90, b1 = 1.05, b1p = -0.57, b2 = 0.0057, b2p = 0.053), weights = pmin(TotalHt^-3.0, 0.5)) # step factor with fit_nlrob()
#otherDiameterFromHeight$sharmaParton = fit_nlrob(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, other2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.01, b2 = 0.26, b3 = 0.5), weights = pmin(TotalHt^-3.0, 0.5), control = nls.control(maxiter = 50)) # step factor
#otherDiameterFromHeight$weibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -138, b1 = 0.0116, b2 = 1.00), weights = pmin(TotalHt^-3.0, 0.5)) # fit_nlrob() > 20 steps
otherDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 3.4, b1 = 0.476, b2 = 0.224), maxit = 50, control = nls.control(maxiter = 50))) # a1p, b1p, b2p not significant, NaN-inf with nls(), no convergence from nls_multstart()
otherDiameterFromHeight$chapmanFormAat = fit_nlrob("Chapman-Richards form AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 6.85, a2 = 0.019, b1 = 0.13, b2 = 0.52), maxit = 50) # NaN-inf
otherDiameterFromHeight$chapmanFormBal = fit_nlrob("Chapman-Richards form BA+L", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.490, a2 = -0.022, a3 = 0.018, b1 = 1.72, b2 = 0.27), maxit = 30, control = nls.control(maxiter = 500)) # NaN-inf with nls()
otherDiameterFromHeight$chapmanFormBalRelHt = fit_gsl_nls("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.65, a2 = -0.03, a3 = 0.03, a9 = 0.1, b1 = 1.5, b2 = 0.27), control = nls.control(maxiter = 50)) # step factor with nls(), step factor or singular gradient with fit_nlrob()
otherDiameterFromHeight$chapmanFormRelHt = fit_gsl_nls("Chapman-Richards form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.09, a9 = -0.098, b1 = 1.14, b2 = 0.38), control = nls.control(maxiter = 1500)) # step factor with nls()
otherDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -7.17, b1 = 0.319, b2 = 0.415, b2p = -0.0044), maxit = 20, control = list(maxiter = 50)) # a1p not significant
otherDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -8.91, a2 = -0.0010, b1 = 0.254, b2 = 0.466, b2p = -0.032), maxit = 150, control = nls.control(maxiter = 500)) # a1p, a2p not significant
otherDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = list(a1 = -6.4, a1p = 0.88, a8 = 0.012, b1 = 0.439, b2 = 0.300), maxit = 500, control = list(maxiter = 500), significant = FALSE) # no physiographic effect significant, convergence fails with b1p, unreliable fit_nlrob() convergence
otherDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -21.3, a9 = 1.5, b1 = 0.11, b2 = 0.60, b2p = 0.006), maxit = 20) # a1p, a2p not significant, fit_nlrob() fails to converge from closer positions
otherDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = other2016gamConstraint), data = other2016)
otherDiameterFromHeight$gamAat = fit_gam("REML GAM AA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 15, pc = other2016gamConstraint), data = other2016)
otherDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM AA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 26, pc = other2016gamConstraint), data = other2016physio)
otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 22, pc = other2016gamConstraint), data = other2016physio)
otherDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), other2016)
otherDiameterFromHeight$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = list(a1 = 28.7, a2 = 13.0, b1 = 0.60), maxit = 30, control = nls.control(maxiter = 50)) # a1p, a2p, b1p not significant
otherDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 0.37, a1p = -0.27, a2 = -0.14, a2p = -0.038)) # converges poorly, yields negative values
otherDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2), other2016) # isPlantation*(TotalHt - 1.37)^2 not significant
otherDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 2.47, a1p = 0.09, b1 = 0.73, b1p = -0.14))
otherDiameterFromHeight$powerAat = fit_nlrob("power AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.08, a2 = -0.0004, b1 = 0.81), maxit = 500) # a2p, b2p not significant, fit_nlrob() failure to converge >500 steps with a1p
otherDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = list(a1 = 1.81, a4 = -0.0007, a5 = -1.05, a6 = 0.034, a7 = -0.024, a8 = 0.00049, b1 = 1.21)) # a1p, a6, a7, a8 not significant
otherDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = list(a1 = 1.32, a9 = 0.40, b1 = 1.11, b1p = -0.055), maxit = 50)
otherDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 2.72, b1 = 0.20, b1p = -0.04, b2 = 0.09, b2p = 0.009), maxit = 40) # step factor with a1p
otherDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = list(a1 = 0.00003, a2 = 0.01, b1 = 1.20, Ha = 50), control = gsl_nls_control(maxiter = 100)) # NaN-inf with fit_nlrob()
otherDiameterFromHeight$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1), other2016, start = list(a1 = 32, b1 = -0.6, b2 = 0.1, b3 = -0.07), control = nls.control(maxiter = 1000)) # gsl_nl() parameter evaporation with b4 allowed to vary from 1, step size or singular gradient with nls(), NaN-inf with nlrob()
otherDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.67, a1p = 1.84, b1 = 1.70, b1p = -1.26, b2 = -0.063, b2p = 0.30), maxit = 50)
otherDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form AA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.58, a1p = 1.93, a2 = -0.0005, b1 = 1.93, b1p = -1.49, b2 = -0.085, b2p = 0.32), maxit = 50, significant = FALSE) # a2, a3 not significant
otherDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = list(a1 = 3.47, a1p = -0.43, a8 = -0.001, b1 = 0.34, b2 = 0.28, b2p = 0.019), maxit = 50, significant = FALSE) # no physiographic predictor significant
otherDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 2.84, a9 = -0.007, b1 = 0.407, b1p = -0.121, b2 = 0.260, b2p = 0.122), maxit = 50, significant = FALSE) # a1p, a2p, a9 not significant
otherDiameterFromHeight$weibull = fit_nlrob("Weibull", DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -153, a1p = 37.0, b1 = 0.09, b2 = 0.42), control = nls.control(maxiter = 50)) # NaN-inf with b1p
#confint_nlrob(otherDiameterFromHeight$weibull, level = 0.99, weights = other2016$dbhWeight)

otherDiameterFromHeightResults = bind_rows(bind_rows(lapply(otherDiameterFromHeight, as_row)),
                                           as_row(name = "modified Sharma-Parton")) %>%
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))


if (includeInvestigatory)
{
  print(otherDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(other2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$sharmaParton), y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman-Richards form approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman-Richards form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    geom_line(aes(x = predict(otherDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$schnute), y = TotalHt, color = "Schnute", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
    #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form aBAL", group = isPlantation), alpha = 0.5) +
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
  #  geom_point(aes(x = TotalHt, y = abs(residuals(otherDiameterFromHeight$michaelisMentenForm)), color = "Michaelis-Menten form", group = isPlantation), alpha = 0.1, color = "grey25", shape= 16)
}
 
 
## other diameter-height GNLS regressions
#otherDiameterFromHeightGnls$naslundGnls = gnls(DBH ~ a1* sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 2.36, a2 = -0.13), weights = varPower(1.4, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls = list(chapmanForm = gnls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeight$chapmanForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()))
otherDiameterFromHeightGnls$chapmanFormAat = gnls(DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeight$chapmanFormAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$chapmanFormBal = gnls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanFormBal$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanFormBalRelHt = gnls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanFormBalRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanFormRelHt = gnls(DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanFormRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanRichards = gnls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichards$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanRichardsAat = gnls(DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichardsAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanRichardsPhysio = gnls(DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = otherDiameterFromHeight$chapmanRichardsPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanRichards = gnls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichards$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$chapmanRichardsRelHt = gnls(DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichardsRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$michaelisMentenForm = gnls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = otherDiameterFromHeight$michaelisMentenForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
otherDiameterFromHeightGnls$naslund = gnls(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 2.28, a1p = 0, a2 = -0.13, a2p = 0), weights = varPower(1.4, ~TotalHt | isPlantation)) # step halving with or without isPlantation
otherDiameterFromHeightGnls$power = gnls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeight$power$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$powerAat = gnls(DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, other2016, start = otherDiameterFromHeight$powerAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$powerPhysio = gnls(DBH ~ (a1 + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = otherDiameterFromHeight$powerPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$powerRelHt = gnls(DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeight$powerRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$ruark = gnls(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = otherDiameterFromHeight$ruark$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
#otherDiameterFromHeightGnls$schnute = gnls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = otherDiameterFromHeight$schnute$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation)) # step halving
#otherDiameterFromHeightGnls$sharmaParton = gnls(DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, other2016, start = otherDiameterFromHeight$sharmaParton$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()) # NaN
otherDiameterFromHeightGnls$sibbesenForm = gnls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenForm$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$sibbesenFormAat = gnls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenFormAat$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$sibbesenFormPhysio = gnls(DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = otherDiameterFromHeight$sibbesenFormPhysio$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$sibbesenFormRelHt = gnls(DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenFormRelHt$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
otherDiameterFromHeightGnls$weibull = gnls(DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = otherDiameterFromHeight$weibull$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())

otherDiameterFromHeightGnls$chapmanForm = get_dbh_error("Chapman-Richards form GNLS", otherDiameterFromHeightGnls$chapmanForm, other2016)
otherDiameterFromHeightGnls$chapmanFormAat = get_dbh_error("Chapman-Richards form AA+T GNLS", otherDiameterFromHeightGnls$chapmanFormAat, other2016)
otherDiameterFromHeightGnls$chapmanFormBal = get_dbh_error("Chapman-Richards form BA+L GNLS", otherDiameterFromHeightGnls$chapmanFormBal, other2016)
otherDiameterFromHeightGnls$chapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BA+L RelHt GNLS", otherDiameterFromHeightGnls$chapmanFormBalRelHt, other2016)
otherDiameterFromHeightGnls$chapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt GNLS", otherDiameterFromHeightGnls$chapmanFormRelHt, other2016)
otherDiameterFromHeightGnls$chapmanRichards = get_dbh_error("Chapman-Richards GNLS", otherDiameterFromHeightGnls$chapmanRichards, other2016)
otherDiameterFromHeightGnls$chapmanRichardsAat = get_dbh_error("Chapman-Richards AA+T GNLS", otherDiameterFromHeightGnls$chapmanRichardsAat, other2016)
otherDiameterFromHeightGnls$chapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio GNLS", otherDiameterFromHeightGnls$chapmanRichardsPhysio, other2016physio)
otherDiameterFromHeightGnls$chapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt GNLS", otherDiameterFromHeightGnls$chapmanRichardsRelHt, other2016)
otherDiameterFromHeightGnls$michaelisMentenForm = get_dbh_error("Michaelis-Menten form GNLS", otherDiameterFromHeightGnls$michaelisMentenForm, other2016)
otherDiameterFromHeightGnls$naslund = get_dbh_error("Näslund GNLS", otherDiameterFromHeightGnls$naslund, other2016)
otherDiameterFromHeightGnls$power = get_dbh_error("power GNLS", otherDiameterFromHeightGnls$power, other2016)
otherDiameterFromHeightGnls$powerAat = get_dbh_error("power AA+T GNLS", otherDiameterFromHeightGnls$powerAat, other2016)
otherDiameterFromHeightGnls$powerPhysio = get_dbh_error("power physio GNLS", otherDiameterFromHeightGnls$powerPhysio, other2016physio)
otherDiameterFromHeightGnls$powerRelHt = get_dbh_error("power RelHt GNLS", otherDiameterFromHeightGnls$powerRelHt, other2016)
otherDiameterFromHeightGnls$ruark = get_dbh_error("Ruark GNLS", otherDiameterFromHeightGnls$ruark, other2016)
#otherDiameterFromHeightGnls$schnute = get_dbh_error("Schnute GNLS", otherDiameterFromHeightGnls$schnute, other2016)
#otherDiameterFromHeightGnls$sharmaParton = get_dbh_error("modified Sharma-Parton GNLS", otherDiameterFromHeightGnls$sharmaParton, other2016)
otherDiameterFromHeightGnls$sibbesenForm = get_dbh_error("Sibbesen form GNLS", otherDiameterFromHeightGnls$sibbesenForm, other2016)
otherDiameterFromHeightGnls$sibbesenFormAat = get_dbh_error("Sibbesen form AA+T GNLS", otherDiameterFromHeightGnls$sibbesenFormAat, other2016)
otherDiameterFromHeightGnls$sibbesenFormPhysio = get_dbh_error("Sibbesen form physio GNLS", otherDiameterFromHeightGnls$sibbesenFormPhysio, other2016physio)
otherDiameterFromHeightGnls$sibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt GNLS", otherDiameterFromHeightGnls$sibbesenFormRelHt, other2016)
otherDiameterFromHeightGnls$weibull = get_dbh_error("Weibull GNLS", otherDiameterFromHeightGnls$weibull, other2016)

otherDiameterFromHeightResultsGnls = bind_rows(lapply(otherDiameterFromHeightGnls, as_row)) %>%
                                               #as_row(otherDiameterFromHeight$schnuteGnls),
                                               #as_row(otherDiameterFromHeight$sharmaPartonGnls),
  mutate(responseVariable = "DBH", species = "other", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  ggplot() +
    geom_histogram(aes(x = power), otherDiameterFromHeightResultsGnls, binwidth = 0.1) +
    geom_segment(aes(x = mean(otherDiameterFromHeightResultsGnls$power), xend = mean(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash") +
    geom_segment(aes(x = median(otherDiameterFromHeightResultsGnls$power), xend = median(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash")
    
  otherDiameterFromHeightResultsGnls %>% summarize(power = mean(power), powerPlantation = mean(powerPlantation))
}

  
## collect model parameters
otherParameters = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, get_coefficients)),
                                      bind_rows(lapply(otherHeightFromDiameterGnls, get_coefficients))) %>%
                                      #get_coefficients(otherHeightFromDiameter$chapmanRichardsGnls),
                                      #get_coefficients(otherHeightFromDiameter$sharmaZhangBalGnls),
                              mutate(responseVariable = "height"),
                            bind_rows(bind_rows(lapply(otherDiameterFromHeight, get_coefficients)),
                                      bind_rows(lapply(otherDiameterFromHeightGnls, get_coefficients))) %>%
                                      #get_coefficients(otherDiameterFromHeight$schnuteGnls),
                                      #get_coefficients(otherDiameterFromHeight$sharmaPartonGnls),
                              mutate(responseVariable = "DBH")) %>%
  mutate(species = "other") %>%
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
if (includeInvestigatory)
{
  otherBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), other2016, start = list(a1 = 1.36, b1 = 0.0002, b2 = 2.06, b2p = -0.27), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 50)) # a1p, b1p not significant, step factor with fit_nlrob()
  otherBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p*isPlantation), other2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = -0.00001, b1 = 2.53, b1p = -0.435), maxit = 70, weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 50))
  #confint2(otherBasalAreaFromHeightKorf, level = 0.99)
  #confint_nlrob(otherBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/other2016$basalArea, 1E4))
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(otherBasalAreaFromHeightKorf), 100^2 * mean(residuals(otherBasalAreaFromHeightKorf)), mean(abs(residuals(otherBasalAreaFromHeightKorf))), 1 - sum(residuals(otherBasalAreaFromHeightKorf)^2) / sum((other2016$basalArea - mean(other2016$basalArea)^2)),
          "power", AIC(otherBasalAreaFromHeightPower), 100^2 * mean(residuals(otherBasalAreaFromHeightPower)), mean(abs(residuals(otherBasalAreaFromHeightPower))), 1 - sum(residuals(otherBasalAreaFromHeightPower)^2) / sum((other2016$basalArea - mean(other2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(other2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(otherBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(otherBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "minority species height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}