# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western redcedar height-diameter regression form sweep
#thplHeightFromDiameter$gamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = thpl2016gamConstraint), data = thpl2016physio, select = TRUE, weights = dbhWeight) # bs = "ts" -> 367, gamma = 2 -> 367, k = 169 min vs 367 default, method = "REML" -> 367
#thplHeightFromDiameter$sharmaPartonBalPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), thpl2016physio, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b1p = -0.0143, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = thplHeightFromDiameterWeights)
thpl2016 = trees2016 %>% filter(Species == "RC", isLiveUnbroken, TotalHt > 0) %>% # live western redcedars measured for height
  mutate(dbhWeight = pmin(DBH^-1.2, 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.7, -1.6), 0.5))
thpl2016physio = thpl2016 %>% filter(is.na(elevation) == FALSE)
thpl2016gamConstraint = c(DBH = -1.2264/0.5099, TotalHt = 1.37, standBasalAreaPerHectare = median(thpl2016$standBasalAreaPerHectare), basalAreaLarger = median(thpl2016$basalAreaLarger), standBasalAreaApprox = median(thpl2016$standBasalAreaApprox), tallerApproxBasalArea = median(thpl2016$tallerApproxBasalArea), elevation = median(thpl2016physio$elevation), slope = median(thpl2016physio$slope), aspect = median(thpl2016physio$aspect), topographicShelterIndex = median(thpl2016physio$topographicShelterIndex), relativeHeight = median(thpl2016$relativeHeight)) # point constraint for mgcv::s()
#thpl2016natural = thpl2016 %>% filter(isPlantation == FALSE)
#thpl2016plantation = thpl2016 %>% filter(isPlantation)
#thpl2016plantationPhysio = thpl2016physio %>% filter(isPlantation)

thplHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131))) # a1p, b1p, b2p not significant
thplHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, thpl2016, start = list(a1 = 51.6, a1p = -9.09, a2 = -0.141, a2p = 0.464, a3 = -0.026, b1 = -0.016, b2 = 1.139)) # a3, a3p, b1p, b2p not significant
thplHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, thpl2016physio, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), significant = FALSE) # a2, a4, a5, a6, a7, a8p, b2p not significant
thplHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), thpl2016, start = list(a1 = 8.7, a1p = -0.810, a2 = 21.1, a2p = 0.239, a3 = -0.043, a9 = 64.0, a9p = -42.8, b1 = -0.021, b2 = 0.120, b2p = 0.916)) # a2, a3, a3p, b1p not significant
thplHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, thpl2016physio, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15)) # a4, a5, a6, a7, a8p, b2p not significant
thplHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, thpl2016, start = list(a1 = 0.560, b1 = 0.069)) # a1p, b1p not significant
thplHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = thpl2016gamConstraint), data = thpl2016) # newton() step failure with family = scat, internal code errors with scat(theta = <fixed val>)
thplHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13, pc = thpl2016gamConstraint), data = thpl2016)
thplHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 20, pc = thpl2016gamConstraint), data = thpl2016physio) # slope and elevation not supported, aspect not tested since insufficient data for full model
thplHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 18, pc = thpl2016gamConstraint), data = thpl2016physio) # k reduces from 85 to 18 without aspect
thplHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), thpl2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176)) # b2p not significant
thplHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), thpl2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175)) # a1p, b1p, b2p not significant
thplHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH, thpl2016) # isPlantation*DBH not significant (p = 0.044)
thplHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176)) # b1p not significant
thplHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2), thpl2016) # isPlantation*DBH not quite significant (p = 0.106), isPlantation*DBH^2 not significant
thplHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), thpl2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649)) # a2p, a3p not significant
thplHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, thpl2016, start = list(a1 = 0.542, b1 = 0.939)) # a1p, b1p not significant
thplHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), thpl2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
thplHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), thpl2016, start = list(Ha = 55.3, Hap = -30.1, d = 0.546, dp = 0.298, kU = 0.009, kUp = 0.017))
thplHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, thpl2016, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09)) # a1p, b3p, b4p not significant
thplHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016, start = list(a1 = 50.6, a1p = -15.8, b1 = 0.023, b2 = -0.014, b2p = -0.009, b3 = -0.069, b4 = 1.130)) # b1p, b3p, b4p not significant
thplHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10)) # b1, b1p, a4, a5, a6, a7, b3p, b4p not significant
thplHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, thpl2016physio, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09)) # a4, a5, a5, a6, a7, b1p, b3p, b4p not significant
thplHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), thpl2016, start = list(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097)) # b1, b1p, b2p, b3p not significant
thplHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, thpl2016, start = list(a1 = 53.2, a1p = -8.857, a2 = -0.002, b1 = -0.016, a2p = 0.010, b2 = -0.025, b3 = -0.078, b4 = 1.126)) # a2, a2p, b1, b1p, b3p, b4p not significant
thplHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), thpl2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078)) # a1p, b1p, b2p not significant
thplHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141)) # b2p not significant
thplHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131)) # a1p, a2, a3, b1p, b2p not significant
thplHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), thpl2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a3 = -0.102, a3p = 0.068, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778)) # a1p, a2, a3, b1p, b2p not significant
#confint_nlrob(thplHeightFromDiameter$sharmaPartonPhysio, level = 0.99, weights = pmin(thpl2016$DBH^-1.2, 1))

thplHeightFromDiameterResults = bind_rows(lapply(thplHeightFromDiameter, as_row)) %>%
  mutate(responseVariable = "height", species = "THPL", deltaAicN = aicN - min(aicN)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAicN))

if (includeInvestigatory)
{
  print(thplHeightFromDiameterResults %>% select(-responseVariable, -species, -fixedWeight, -n, -power, -significant, -contains("NaturalRegen"), -contains("Plantation")), n = 30)
  ggplot() +
    geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = thpl2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = thpl2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$curtis), color = "Curtis", group = thpl2016$isPlantation)) +
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$gam), color = "GAM", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$korf), color = "Korf", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$linear), color = "linear", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$parabolic), color = "parabolic", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$power), color = "power", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$prodan), color = "Prodan", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$richards), color = "unified Richards", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$sibbesen), color = "Sibbesen", group = thpl2016$isPlantation)) +
    #geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibull), color = "Weibull", group = thpl2016$isPlantation)) +
    annotate("text", x = 0, y = 65, label = "western redcedar, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 65)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}

# dbhClassSize = 50
# errorByDbhClass = tibble(dbhClass = dbhClassSize*floor(thpl2016$DBH/dbhClassSize) + 0.5*dbhClassSize, fittedValue = predict(thplHeightFromDiameter$gam, thpl2016), height = thpl2016$TotalHt, residual = fittedValue - height) %>%
#  #mutate(residual = residual - if_else(dbhClass == 50, -0.477/376, 0.107/95)) %>%
#  group_by(dbhClass) %>%
#  summarize(n = n(),
#            totalHeight = sum(height),
#            totalFitted = sum(fittedValue),
#            meanBiasPerTree = sum(residual) / n,
#            meanBiasPerTreePct = 100 * sum(residual/height) / n,
#            minError = min(residual),
#            meanError = mean(residual),
#            maxError = max(residual),
#            minPct = 100 * min(residual/height),
#            meanPct = 100 * mean(residual/height),
#            maxPct = 100 * max(residual/height),
#            .groups = "drop") %>%
#  filter(n >= 10)
# errorByDbhClass


## western redcedar height-diameter GNLS regressions
#thplHeightFromDiameter = list(chapmanRichards = gnls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, thpl2016, start = thplHeightFromDiameter$chapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = TRUE))) # step halving with varPower(0.60)
#thplHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, thpl2016, start = thplHeightFromDiameterGnls$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameter$sharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameter$sharmaPartonBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^b2*DBH))^(b3 + b3p * isPlantation), thpl2016, start = thplHeightFromDiameter$sharmaZhang$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, thpl2016, start = thplHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving with plot correlation
#thplHeightFromDiameterGnls$weibull = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), thpl2016, start = thplHeightFromDiameter$weibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE))
#thplHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), thpl2016, start = thplHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.60, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#save(thplHeightFromDiameterGnls, file = "trees/height-diameter/HtDia THPL GNLS.rdata")

load("trees/height-diameter/HtDia THPL GNLS.rdata")
thplHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", thplHeightFromDiameterGnls$chapmanRichards, thpl2016)
thplHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", thplHeightFromDiameterGnls$chapmanRichardsBal, thpl2016)
thplHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", thplHeightFromDiameterGnls$sharmaParton, thpl2016)
thplHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", thplHeightFromDiameterGnls$sharmaPartonBal, thpl2016)
thplHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", thplHeightFromDiameterGnls$sharmaZhang, thpl2016)
thplHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", thplHeightFromDiameterGnls$sharmaZhangBal, thpl2016)
thplHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", thplHeightFromDiameterGnls$weibull, thpl2016)
thplHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", thplHeightFromDiameterGnls$weibullBal, thpl2016)

thplHeightFromDiameterResultsGnls = bind_rows(lapply(thplHeightFromDiameterGnls, as_row)) %>%
  mutate(responseVariable = "height", species = "THPL", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  thplHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(thplHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(thplHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(thplHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
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
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibullBal), color = "Weibull BA+L"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$power), color = "power")) +
    geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameter$weibull), color = "Weibull")) +
    annotate("text", x = 0, y = 85, label = "a) western redcedar, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}
 
 
## western redcedar diameter-height regressions
#thplDiameterFromHeight$chapmanForm = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, iter = 1000,
#                                                  start_lower = list(a1 = -10, b1 = -1, b2 = -1), 
#                                                  start_upper = list(a1 = 50, b1 = 1, b2 = 1), modelweights = heightWeight)
#thplDiameterFromHeight$chapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a4 * pmin(relativeHeight, 1.25))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 115, a2 = 0, a4 = 13.3, b1 = 0.020, b2 = 0.92), control = nls.control(maxiter = 500)) # step size with nls(), step factor with BA+L+nlrob()
#thplDiameterFromHeight$chapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 139, a2 = 10, a3 = -0.066, a4 = 11.6, b1 = 0.016, b2 = 0.944), control = gsl_nls_control(maxiter = 250)) # step factor with BA+L+nlrob()
#thplDiameterFromHeight$chapmanRichards = nls_multstart(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, iter = 100,
#                                                      start_lower = list(a1 = 1, b1 = -0.5, b2 = 0.1), 
#                                                      start_upper = list(a1 = 300, b1 = 0.5, b2 = 2.5), modelweights = heightWeight)
#thplDiameterFromHeight$chapmanRichards = gsl_nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), thpl2016, start = list(a1 = 229.6, b1 = -0.0057, b2 = 1.25, b2p = 0), weights = heightWeight)
#thplDiameterFromHeight$chapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, iter = 100,
#                                                            start_lower = list(a1 = -200, a1p = -100, a2 = -0.1, a3 = -1, a4 = -1, a5 = -1, a6 = -0.1, b1 = -0.1, b1p = -1, b2 = -1), 
#                                                            start_upper = list(a1 = 1, a1p = 100, a2 = 0.1, a3 = 1, a4 = 1, a5 = 1, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 1), modelweights = heightWeight)
#thplDiameterFromHeight$gamPhysio = gam(DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation), k = 85, pc = thpl2016gamConstraint), data = thpl2016physio, method = "REML", select = TRUE, weights = heightWeight)
#thplDiameterFromHeight$sharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, thpl2016, iter = 100,
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = -0.5, b2 = -0.5, b3 = -1), 
#                                                   start_upper = list(a1 = 10, a2 = 1.5, b1 = 0.5, b2 = 0.5, b3 = 1), modelweights = heightWeight)
#thplDiameterFromHeight$schnute = nlrob(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016, start = list(a1 = 0.0003, a2 = 0.1, b1 = 1.5, Ha = 75), weights = heightWeight) # NaN-inf with nlrob()
thplDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 960, b1 = 0.003, b2 = 1.05))) # a1p, b1p, b2p not significant, singular gradient with nls(), no convergence from nls_multstart()
thplDiameterFromHeight$chapmanFormAat = fit_gsl_nls("Chapman-Richards form AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, thpl2016, start = list(a1 = 960, a2 = -100, b1 = 0.001, b2 = 1.07), control = nls.control(maxiter = 500), significant = FALSE) # NaN-inf with nls() and nlrob
thplDiameterFromHeight$chapmanFormBal = fit_gsl_nls("Chapman-Richards form BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 913, a2 = -40, b1 = 0.0003, b2 = 1.04), control = nls.control(maxiter = 300), significant = FALSE) # step size with nls() and nlrob()
thplDiameterFromHeight$chapmanFormBalRelHt = fit_gsl_nls("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.25)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 655, a2 = 0, a3 = 0, a9 = 2.3, b1 = 0.003, b2 = 1.04), control = gsl_nls_control(maxiter = 500), significant = FALSE) # nlrob() step factor with either a2 or a3
thplDiameterFromHeight$chapmanFormRelHt = fit_gsl_nls("Chapman-Richards form RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), thpl2016, start = list(a1 = 655, a9 = 2.3, b1 = 0.003, b2 = 1.04), control = nls.control(maxiter = 500)) # step size with nls(), >50 iterations with nlrob()
thplDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -763, b1 = 0.0027, b2 = 1.05)) # a1p and b2p not significant, poor convergence with b1p
thplDiameterFromHeight$chapmanRichardsAat = fit_gsl_nls("Chapman-Richards AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -1200, a2 = 0, b1 = 0.0016, b2 = 1.06), control = nls.control(maxiter = 500), significant = FALSE) # a1p, b1p not significant, step factor with nlrob()
thplDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), thpl2016physio, start = list(a1 = -480, a1p = 270, a8 = 1.1, b1 = 0.0048, b1p = 0.006, b2 = 1.0), control = list(maxiter = 50)) # no physiographic effects significant
thplDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), thpl2016, start = list(a1 = -350, a9 = -66, b1 = 0.006, b2 = 0.98), control = list(maxiter = 500), significant = FALSE) # step factor with nlrob()
thplDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = thpl2016gamConstraint), data = thpl2016) # newton() step failure with scat()
thplDiameterFromHeight$gamAat = fit_gam("REML GAM AA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = thpl2016gamConstraint), data = thpl2016)
thplDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM AA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 23, pc = thpl2016gamConstraint), data = thpl2016physio)
thplDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 23, pc = thpl2016gamConstraint), data = thpl2016physio)
thplDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), thpl2016) # isPlantation*(TotalHt - 1.37) not significant
thplDiameterFromHeight$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), thpl2016, start = list(a1 = 519, a2 = 237, b1 = 1.00)) # a1p, a2p, b1p not significant
thplDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), thpl2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024))
thplDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), thpl2016) # (TotalHt - 1.37)^2 not significant
thplDiameterFromHeight$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.93, b1 = 1.08)) # no significant plantation effects
thplDiameterFromHeight$powerAat = fit_nlrob("power AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09)) # no significant plantation effects
thplDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, thpl2016physio, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE) # no significant physiographic effects
thplDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, thpl2016, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13)) # a1p and b1p not significant
thplDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.12, b1 = 1.01, b2 = 0.0038, b2p = 0.0025)) # a1p, b1p not significant
thplDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), thpl2016, start = list(a1 = 0.00005, a2 = 0.0074, b1 = 1.12, Ha = 52), control = gsl_nls_control(maxiter = 200)) # singular gradient with nlrob()
thplDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, thpl2016, start = list(a1 = 65, b1 = 0.4, b2 = 0.006, b3 = 0.03, b4 = 0.7), control = gsl_nls_control(maxiter = 200)) # NaN-inf with nls() from nls_multstart() point, NaN-inf, singular gradient, or code syntax error with nlrob()
thplDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.42, b1 = 1.29, b2 = -0.27)) # no significant plantation effects
thplDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form AA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029)) # no significant plantation effects
thplDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), thpl2016physio, start = list(a1 = 2.32, a1p = -0.294, a8 = -0.0056, b1 = 1.084, b2 = -0.0048, b2p = 0.0203), significant = FALSE) # no physiographic effects significant, though reduction makes a8 almost significant
thplDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), thpl2016, start = list(a1 = 1.34, a9 = 0.30, b1 = 1.32, b2 = -0.040))
thplDiameterFromHeight$weibull = fit_nlrob("Weibull", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, thpl2016, start = list(a1 = -330, b1 = 0.0065, b2 = 1.01)) # a1p, b1p, b2p not significant
#confint_nlrob(thplDiameterFromHeight$sibbesenFormPhysio, level = 0.99, weights = pmin(thpl2016$TotalHt^if_else(thpl2016$isPlantation, -1.7, -1.6), 0.5))

thplDiameterFromHeightResults = bind_rows(lapply(thplDiameterFromHeight, as_row)) %>%
                                          #as_row(thplDiameterFromHeight$gamAatPhysio),
  mutate(responseVariable = "DBH", species = "THPL", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))

if (includeInvestigatory)
{
  print(thplDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(thpl2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$sharmaParton), y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman-Richards form approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman-Richards form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$schnute), y = TotalHt, color = "Schnute", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
    #geom_line(aes(x = predict(thplDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -100 * log(1 - pmin(0.015*(TotalHt - 1.37)^1.0, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
    #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.45, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AA+T", group = isPlantation), alpha = 0.5) +
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
}


## collect model parameters
thplParameters = bind_rows(bind_rows(bind_rows(lapply(thplHeightFromDiameter, get_coefficients)),
                                     bind_rows(lapply(thplHeightFromDiameterGnls, get_coefficients))) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(lapply(thplDiameterFromHeight, get_coefficients)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "THPL") %>%
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
if (includeInvestigatory)
{
thplBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), thpl2016, start = list(a1 = 90, b1 = 0.000003, b2 = 2.18), weights = pmin(1/basalArea, 1E4)) # a1p, b1p, b2p not significant, step factor with nlrob()
thplBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), thpl2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, a1p = -0.00002, b1 = 2.14, b1p = 0.34), weights = pmin(1/basalArea, 1E4))
#confint2(thplBasalAreaFromHeightKorf, level = 0.99)
#confint_nlrob(thplBasalAreaFromHeightPower, level = 0.99, weights = pmin(1/thpl2016$basalArea, 1E4))

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(thplBasalAreaFromHeightKorf), 100^2 * mean(residuals(thplBasalAreaFromHeightKorf)), mean(abs(residuals(thplBasalAreaFromHeightKorf))), 1 - sum(residuals(thplBasalAreaFromHeightKorf)^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2)),
        "power", AIC(thplBasalAreaFromHeightPower), 100^2 * mean(residuals(thplBasalAreaFromHeightPower)), mean(abs(residuals(thplBasalAreaFromHeightPower))), 1 - sum(residuals(thplBasalAreaFromHeightPower)^2) / sum((thpl2016$basalArea - mean(thpl2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(thpl2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = predict(thplBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = predict(thplBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "western redcedar height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## exploratory plots
if (includeInvestigatory)
{
  library(GGally)
  ggpairs(thpl2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger, relativeHeight, regeneration), 
          aes(alpha = 0.1, color = regeneration, shape = "16"),
          columnLabels = c("DBH, cm", "height, m", "BA, m² ha⁻¹", "BAL, m² ha⁻¹", "relative height, %", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    #scale_color_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) + # https://github.com/ggobi/ggally/issues/445
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  ggpairs(thpl2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, slope, elevation, topographicShelterIndex, regeneration), 
          aes(alpha = 0.1, color = if_else(thpl2016$isPlantation, "plantation", "natural regen"), shape = "16"), 
          columnLabels = c("DBH, cm", "height, m", "slope, °", "elevation, m", "TSI, °", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  scatterPlotMatrix::scatterPlotMatrix(thpl2016 %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger))
}