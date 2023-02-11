# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western hemlock height-diameter regression form sweep
tshe2016 = trees2016 %>% filter(Species == "WH", isLiveUnbroken, TotalHt > 0) %>% # live western hemlocks measured for height
  mutate(dbhWeight = pmin(DBH^-0.9, 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tshe2016physio = tshe2016 %>% filter(is.na(elevation) == FALSE)
tshe2016gamConstraint = c(DBH = -1.2994/0.6005, TotalHt = 1.37, standBasalAreaPerHectare = median(tshe2016$standBasalAreaPerHectare), basalAreaLarger = median(tshe2016$basalAreaLarger), standBasalAreaApprox = median(tshe2016$standBasalAreaApprox), tallerApproxBasalArea = median(tshe2016$tallerApproxBasalArea), elevation = median(tshe2016physio$elevation), slope = median(tshe2016physio$slope), aspect = median(tshe2016physio$aspect), topographicShelterIndex = median(tshe2016physio$topographicShelterIndex), relativeHeight = median(tshe2016$relativeHeight)) # point constraint for mgcv::s()
#tshe2016natural = tshe2016 %>% filter(isPlantation == FALSE)
#tshe2016plantation = tshe2016 %>% filter(isPlantation)
#tshe2016plantationPhysio = tshe2016physio %>% filter(isPlantation)

tsheHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048)))
tsheHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^b2, tshe2016physio, start = list(a1 = 60.5, a2 = 0.026, a2p = 0.96, a3 = 0, a3p = 0, a4 = -0.012, a5 = -0.19, a6 = 0.44, a7 = 0.80, a8 = 0.15, b1 = -0.022, b2 = 1.31)) # a1p, a6, b1p, b2p not significant
tsheHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a9 = 63.0, a9p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01)) # a3, a3p, b1p not significant
tsheHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58.5, a1p = -3.5, a4 = -0.015, a5 = -8.07, a7 = 0.718, a8 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12)) # a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054)) # b1 not significant
tsheHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = tshe2016gamConstraint), data = tshe2016)
tsheHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 19, pc = tshe2016gamConstraint), data = tshe2016)
tsheHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 90, pc = tshe2016gamConstraint), data = tshe2016physio) # ~2 minutes to fit+evaluate with all predictors at minimum (Zen 3 @ 3.4 GHz, k = 495, edf < 290) -> eliminate cos(aspect) (k = 330, edf < 240) AIC 11023: 10733 without BA, 10887 without BAL, 10729 without elevation, 10714 without slope, 10679 without sin(aspect), 10669 without topographic shelter -> eliminate sin(aspect)
tsheHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio)
tsheHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047))
tsheHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), tshe2016, start = list(a1 = 200, b1 = -7.2, b2 = -0.33)) # a1p, b1p, b2p not significant
tsheHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), tshe2016)
tsheHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), tshe2016)
tsheHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93)) # a3p not significant
tsheHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05))
tsheHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1)) # b2p not significant
tsheHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013)) # dp not significant
tsheHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 30.7, b1 = 0.15, b1p = -0.068, b2 = -0.034, b3 = -0.22, b3p = 0.184, b4 = 1.26)) # a1p, b2p, b3p not significant
tsheHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = list(a1 = 46.7, a1p = -12.6, b1 = 0.054, b2 = -0.021, b2p = -0.013, b3 = -0.054, b4 = 1.26)) # b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26)) # a6, a8, b1p, b3p, b4p not significant
tsheHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = list(a1 = 33.4, a1p = -8.1, b1 = 0.138, b2 = -0.027, b3 = -0.053, b3p = 0.077, b4 = 1.27)) # b1p, b2p, b4p not significant
tsheHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.0005, a2p = 0.0191, b1 = 0.100, b1p = -0.207, b2 = -0.020, b3 = -0.029, b4 = 1.24)) # a3, b1, b2p, b3p, b4p not significant
tsheHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014)) # b1p not significant
tsheHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029))
tsheHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212)) # a1p, b1p, b2p not significant
tsheHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 16.0, a1p = -9.7, a2 = 0.018, a2p = 0.63, a3 = 0.28, a9 = 73.5, b1 = -0.015, b2 = 0.76), control = list(maxiter = 200)) # a3p, a4p, b1p, b2p not significant
#confint_nlrob(tsheHeightFromDiameter$sharmaPartonBalPhysio, level = 0.99)

if (htDiaOptions$includeInvestigatory)
{
  print(tsheHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = tshe2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = tshe2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$curtis), color = "Curtis", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$korf), color = "Korf", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$linear), color = "linear", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$parabolic), color = "parabolic", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$power), color = "power", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$prodan), color = "Prodan", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$richards), color = "unified Richards", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$sibbesen), color = "Sibbesen", group = tshe2016$isPlantation)) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibull), color = "Weibull", group = tshe2016$isPlantation)) +
    annotate("text", x = 0, y = 70, label = "western hemlock, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 70)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}


## western hemlock height-diameter GNLS regressions
if (htDiaOptions$fitSlowGnls)
{
  tsheHeightFromDiameterGnls = list(chapmanRichards = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = tsheHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))) # step halving at nlsTol = 0.001
  tsheHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = tsheHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
  tsheHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
  tsheHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
  tsheHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
  tsheHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^(b1 + b1p * isPlantation) * (1 - exp(b2*tph^b3*DBH))^b4, tshe2016, start = tsheHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
  tsheHeightFromDiameterGnls$weibull = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = tsheHeightFromDiameter$weibull$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-7, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # >250+50 iterations at nlsTol = 0.1, ok at msTol = 1E-3 and tol = 1E-3
  tsheHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = tsheHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
  
  tsheHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", tsheHeightFromDiameterGnls$chapmanRichards, tshe2016)
  tsheHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", tsheHeightFromDiameterGnls$chapmanRichardsBal, tshe2016)
  tsheHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", tsheHeightFromDiameterGnls$sharmaParton, tshe2016)
  tsheHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", tsheHeightFromDiameterGnls$sharmaPartonBal, tshe2016)
  tsheHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", tsheHeightFromDiameterGnls$sharmaZhang, tshe2016)
  tsheHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", tsheHeightFromDiameterGnls$sharmaZhangBal, tshe2016)
  tsheHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", tsheHeightFromDiameterGnls$weibull, tshe2016)
  tsheHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", tsheHeightFromDiameterGnls$weibullBal, tshe2016)

  save(tsheHeightFromDiameterGnls, file = "trees/height-diameter/data/TSHE GNLS.rdata")
} else {
  load("trees/height-diameter/data/TSHE GNLS.rdata")
}

if (htDiaOptions$includeInvestigatory)
{
  tsheHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(tsheHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(tsheHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(tsheHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  
  ggplot() +
    geom_point(aes(x = tshe2016natural$DBH, y = tshe2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = tshe2016natural$DBH, y = tshe2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "western hemlock naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = tshe2016plantation$DBH, y = tshe2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = tshe2016plantation$DBH, y = tshe2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "western hemlock plantation height, m")
  
  ggplot() +
    geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibullBAL), color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016natural$DBH, y = predict(tsheHeightFromDiameter$weibullBALnatural), color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016plantation$DBH, y = predict(tsheHeightFromDiameter$weibullBALplantation), color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterBase), color = "base")) +
    geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameter$weibull), color = "ElliottWeibull")) +
    annotate("text", x = 0, y = 85, label = "a) western hemlock, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## western hemlock diameter-height regressions
tsheDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 136, b1 = 0.0100, b2 = 0.924))) # no significant plantation effects
tsheDiameterFromHeight$chapmanFormAat = fit_nlrob("Chapman-Richards form ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 135, a2 = -0.0045, b1 = 0.010, b2 = 0.924))
tsheDiameterFromHeight$chapmanFormBal = fit_nlrob("Chapman-Richards form BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, b1 = 0.0138, b2 = 0.893)) # a1p not significant, a3 + b2p step size
tsheDiameterFromHeight$chapmanFormBalRelHt = fit_nlrob("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 189, a2 = -3.64, a3 = 2.24, a9 = -46.4, b1 = 0.0137, b2 = 0.847))
tsheDiameterFromHeight$chapmanFormRelHt = fit_nlrob("Chapman-Richards form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 128, a9 = 5.3, b1 = 0.0154, b2 = 0.902))
tsheDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -119, a1p = 26.4, b1 = 0.0196, b1p = 0.0004, b2 = 0.847))
tsheDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -119, a1p = 26.4, a2 = 0, b1 = 0.0196, b1p = 0.0004, b2 = 0.847))
tsheDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope))*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -189, a1p = 71.1, a4 = -0.017, a5 = -0.339, b1 = 0.0095, b1p = 0.00375, b2 = 0.919)) # a6, a7, a8 not significant
tsheDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a9 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912)) # job step factor with nlrob()
tsheDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = tshe2016gamConstraint), data = tshe2016)
tsheDiameterFromHeight$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = tshe2016gamConstraint), data = tshe2016)
tsheDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio, nthreads = 4) # k = 495, ef < 210 with all predictors AIC 13838: 13828 without AAT, 13826 without ABA, 13825 without elevation, 13793 without slope, 13848 without sin(aspect), 13896 without cos(aspect), 13772 without topographic shelter -> eliminate topographic shelter (k = 330, edf < 185) AIC 13772: 13654 without AAT, 13665 without ABA, 13649 without elevation, 13649 without slope, 13635 without sin(aspect), 13656 without cos(aspect) -> eliminate sin(aspect)
tsheDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = tshe2016gamConstraint), data = tshe2016physio)
tsheDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), tshe2016) # isPlantation*(TotalHt - 1.37) not significant(p = 0.036)
tsheDiameterFromHeight$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83)) # a1p, a2p, b1p not significant
tsheDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013))
tsheDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016)
tsheDiameterFromHeight$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04)) # no significant plantation effects
tsheDiameterFromHeight$powerAat = fit_nlrob("power ABA+T", DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047)) # a1p not significant
tsheDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a5 = 0.284, b1 = 1.04)) # a4, a6, a7, a8 not significant
tsheDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a9 = 0.08, b1 = 1.02)) 
tsheDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.31, b1 = 0.818, b2 = 0.010)) # a1p, b1p, b2p not significant
#tsheDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32)) # converges from red alder values but fails to reconverge (singular gradient), NaN-inf or singular gradient with fit_nlrob()
#tsheDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, tshe2016, start = list(a1 = 1.4, b1 = 1.0, b2 = 0.001, b3 = 0.3, b4 = 0.9), control = gsl_nls_control(maxiter = 125, xtol = 1E-6), significant = FALSE) # unreliable convergence and NaN-inf with fit_gsl_nls() due to parameter evaporation and collapse to power form, NaN-inf or singular gradient with fit_nlrob(), singular gradient with nls() even with nls_multstart() parameters
tsheDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057)) # no significant plantation effects
tsheDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578)) # no significant plantation effects
tsheDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), tshe2016physio, start = list(a1 = 2.115, a1p = -0.100, a4 = 0.00017, a5 = 0.486, b1 = 0.736, b2 = 0.0593, b2p = 0.0028)) # a6, a7, a8 not significant
tsheDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a9 = 0.084, b1 = 0.75, b2 = 0.056))
tsheDiameterFromHeight$weibull = fit_nlrob("Weibull", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -225, b1 = 0.011, b2 = 0.88)) # a1p, b1p, b2p not significant
#confint_nlrob(tsheDiameterFromHeight$powerPhysio, level = 0.99)

if (htDiaOptions$includeInvestigatory)
{
  print(tsheDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(tshe2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$sharmaParton), y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman-Richards form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman-Richards form ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen", group = isPlantation)) +
    #geom_line(aes(x = predict(tsheDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    ##geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = -1/0.0005*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute"), alpha = 0.5) +
    geom_line(aes(x = 1.3*(TotalHt - 1.37)^1*exp(0.003*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 81, label = "western hemlock, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## collect model parameters
tsheCoefficients = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_coefficients)),
                                       bind_rows(lapply(tsheHeightFromDiameterGnls, get_model_coefficients)) %>%
                               mutate(responseVariable = "height"),
                             bind_rows(lapply(tsheDiameterFromHeight, get_list_coefficients)) %>%
                                       #get_model_coefficients(tsheDiameterFromHeight$schnute),
                                       #get_model_coefficients(tsheDiameterFromHeight$sharmaParton),
                               mutate(responseVariable = "DBH"))) %>%
  mutate(species = "TSHE") %>% 
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a5, a6, a7, a8, a9, b1, b1p, b2, b2p, b3)

tsheResults = bind_rows(bind_rows(bind_rows(lapply(tsheHeightFromDiameter, get_list_results)),
                                  bind_rows(lapply(tsheHeightFromDiameterGnls, get_model_results))) %>%
                          mutate(responseVariable = "height"),
                        bind_rows(bind_rows(lapply(tsheDiameterFromHeight, get_list_results)),
                                  get_model_results(name = "Schnute", fitting = "gsl_nls"),
                                  get_model_results(name = "modified Sharma-Parton", fitting = "gsl_nls")) %>%
                          mutate(responseVariable = "DBH")) %>%
  mutate(species = "TSHE") %>%
  relocate(responseVariable, species)

save(file = "trees/height-diameter/data/TSHE results.Rdata", tsheCoefficients, tsheResults)
save(file = "trees/height-diameter/data/TSHE models.Rdata", tsheHeightFromDiameter, tsheHeightFromDiameterGnls, tsheDiameterFromHeight)


## preferred forms identified (results.R, Figure 5)
tsheHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1))
tsheHeightFromDiameterPreferred$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 48.3, a1p = -12.9, b1 = 0.090, a4 = -0.011, a5 = 0, a7 = -0.19, b2 = -0.026, b2p = -0.021, b3 = -0.15, b4 = 1.26), folds = 1, repetitions = 1)
tsheHeightFromDiameterPreferred$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a7 * cos(3.14159/180 * aspect))*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, tshe2016physio, start = list(a1 = 43.7, a1p = -11.7, a4 = -0.010, a5 = 0, a7 = -0.14, b1 = 0.114, b2 = -0.028, b2p = -0.022, b3 = -0.14, b4 = 1.26), folds = 1, repetitions = 1)

tsheDiameterFromHeightPreferred = list(chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a9 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912), folds = 1, repetitions = 1))
tsheDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1)
tsheDiameterFromHeightPreferred$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = tshe2016gamConstraint), data = tshe2016, folds = 1, repetitions = 1)
tsheDiameterFromHeightPreferred$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057), folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/TSHE preferred models.Rdata", tsheHeightFromDiameterPreferred, tsheDiameterFromHeightPreferred)


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  #tsheBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ a1*exp(b1*(imputedHeight - 1.37)^b2) - 1, tshe2016, start = list(a1 = 1, b1 = 0.00009, b2 = 2.2), weights = pmin(1/basalArea, 1E4))
  tsheBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 5.5, b1 = 0.00003, b2 = 2.09, b2p = -0.060), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor on b1p without a1p
  tsheBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 11/7 * 0.25 * pi * 0.01^2, a1p = -0.00005, b1 = 2.18, b1p = -0.178), weights = pmin(1/basalArea, 1E4))
  #confint_nlrob(tsheBasalAreaFromHeightKorf, level = 0.99, weights = pmin(1/tshe2016$basalArea, 1E4))
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(tsheBasalAreaFromHeightKorf), 100^2 * mean(residuals(tsheBasalAreaFromHeightKorf)), mean(abs(residuals(tsheBasalAreaFromHeightKorf))), 1 - sum(residuals(tsheBasalAreaFromHeightKorf)^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2)),
          "power", AIC(tsheBasalAreaFromHeightPower), 100^2 * mean(residuals(tsheBasalAreaFromHeightPower)), mean(abs(residuals(tsheBasalAreaFromHeightPower))), 1 - sum(residuals(tsheBasalAreaFromHeightPower)^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(tshe2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(tsheBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(tsheBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "western hemlock height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}
