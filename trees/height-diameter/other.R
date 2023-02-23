# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## minority species height-diameter regression form sweep
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 43.0, a1p = 13.5, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875), weights = dbhWeight, control = gsl_nls_control(maxiter = 50)) # a1p not significant
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = 43.0, a2 = 0.46, a3 = 0.082, b1 = -0.00867, b2 = 0.875, b2p = 0), weights = dbhWeight, control = gsl_nls_control(maxiter = 500)) # > 500 iterations
#otherHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 89.3, a2 = 1.166, a2p = 0, a3 = -0.039, a3p = 0, b1 = -0.00544, b2 = 0.873), weights = dbhWeight) # a2p, a3p not significant
#otherHeightFromDiameter$sharmaZhang = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), other2016, start = list(a1 = 58.1, a1p = -47.4, a2 = 0.162, a2p = 0.032, b1 = -0.021, b1p = -0.222, b2 = -0.292, b2p = 0.036, b3 = 0.818, b3p = 0.165), weights = dbhWeight)
#otherHeightFromDiameter$sibbesen = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.714, a1p = 0.088, b1 = 1.172, b2 = -0.074, b2p = -0.0040), weights = dbhWeight) # a1p, b2p not significant
other2016 = trees2016 %>% filter((Species %in% c("DF", "RA", "WH", "BM", "OM", "RC")) == FALSE, isLiveUnbroken, TotalHt > 0) %>%
  mutate(dbhWeight = pmin(1/(1.13*DBH^0.87), 5),
         heightWeight = pmin(1/(1.31*(TotalHt - 1.37)^1.82), 5))
other2016physio = other2016 %>% filter(is.na(elevation) == FALSE)
other2016gamConstraint = c(DBH = -1.4359/0.8154, TotalHt = 1.37, standBasalAreaPerHectare = median(other2016$standBasalAreaPerHectare), basalAreaLarger = median(other2016$basalAreaLarger), standBasalAreaApprox = median(other2016$standBasalAreaApprox), tallerApproxBasalArea = median(other2016$tallerApproxBasalArea), elevation = median(other2016physio$elevation), slope = median(other2016physio$slope), aspect = median(other2016physio$aspect), topographicShelterIndex = median(other2016physio$topographicShelterIndex), relativeHeight = median(other2016$relativeHeight)) # point constraint for mgcv::s()
#other2016natural = other2016 %>% filter(isPlantation == FALSE)
#other2016plantation = other2016 %>% filter(isPlantation)
#other2016plantationPhysio = other2016physio %>% filter(isPlantation)
#otherConifer2016 = other2016 %>% filter(Species %in% c("XX", "CX", "SS", "PC", "PY", "GF", "LP"))
#otherHardwood2016 = other2016 %>% filter(Species %in% c("CA", "HX", "CH", "PM", "GC", "PD", "TO", "WI", "OA", "WO"))

otherHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1*(1 - exp((b1 + b1p * isPlantation) * DBH))^b2, other2016, start = list(a1 = 24.7, b1 = -0.033, b1p = -0.003, b2 = 1.000))) # a1p, b2p not significant
otherHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = list(a1 = 45, a2 = 1.0, a3 = -0.5, b1 = -0.01, b2 = 0.9)) # a1p, a2p, a2, a3p, a3, b1p, b2p not significant, semi-regular NaN-inf with nlrob()
otherHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * elevation) * (1 - exp(b1*DBH))^b2, other2016physio, start = list(a1 = 80, a2 = 0.3, a4 = 0.005, b1 = -0.005, b2 = 0.8), significant = FALSE) # a1p, a2p, a3, a4, a5, a6, a7, a8, b1p, b2p not significant, occasional step factor with nlrob()
otherHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), other2016, start = list(a1 = -1.16, a2 = -0.05, a3 = 0.07, a3p = 0.09, a4 = 16.6, a4p = -6.7, b1 = -0.05, b2 = 0.28, b2p = 0.131)) # a2, a2p, a3, b1p not significant
otherHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a4 * elevation) * (1 - exp(b1*DBH))^b2, other2016physio, start = list(a1 = 70, a4 = -0.01, b1 = -0.006, b2 = 0.81), control = nls.control(maxiter = 250), significant = FALSE) # a1p, a4, a4p, a5, a6, a7, a8, b1p, b2p not significant
otherHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1*DBH / (1 + DBH)^b1, other2016, start = list(a1 = 1.086, b1 = 0.190)) # a1p, b1p not significant
otherHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 7, pc = other2016gamConstraint), data = other2016)
otherHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15, pc = other2016gamConstraint), data = other2016) # slow with unrestricted k
otherHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 23, pc = other2016gamConstraint), data = other2016physio) # unrestricted k tensor product (te() + te()) impractically slow (>2 h, Zen 3 @ 3.4 GHz)
otherHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 17, pc = other2016gamConstraint), data = other2016physio)
otherHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), other2016, start = list(a1 = 34.6, b1 = 40.2, b2 = -1.03)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), other2016, start = list(a1 = 381, b1 = -6.26, b2 = -0.20)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), other2016)
otherHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 100, a2 = 100, b1 = 0.86)) # a1p, a2p, b1p not significant
otherHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), other2016)
otherHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), other2016, start = list(a1 = 0.028, a2 = 1.08, a3 = 0.15)) # a1p, a2p, a3p not significant
otherHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, other2016, start = list(a1 = 0.99, b1 = 0.84)) # a1p, b1p not significant
otherHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), other2016, start = list(a1 = 40, b1 = -32, b2 = 9)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), other2016, start = list(Ha = 12.9, d = 2.59, kU = 0.064, kUp = 0.008)) # Hap, dp not significant
otherHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^b1*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, other2016, start = list(a1 = 150, b1 = -0.15, b2 = -0.01, b3 = -0.3, b4 = 0.77)) # a1p, b1p, b2p, b3p, b4p not significant, step factor with nlrob()
otherHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, other2016, start = list(a1 = 200, b1 = -0.15, b2 = -0.003, b3 = -0.33, b4 = 0.77), control = gsl_nls_control(maxiter = 500, xtol = 1E-6)) # a1p, b1p, b2p, b3p, b4p not significant, potential a1-b1 evaporation, intermittent NaN-infs with nlrob()
otherHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, other2016physio, start = list(a1 = 300, a4 = 0.0, b1 = -0.2, b2 = -0.003, b3 = -0.34, b4 = 0.76), control = gsl_nls_control(maxiter = 500, xtol = 1E-5), significant = FALSE) # a4, a5, a6, a7, a8, b1, b2p, b3p, b4p not significant, intermittent step factor with nlrob()
otherHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, other2016physio, start = list(a1 = 140, a4 = 0, b1 = -0.14, b2 = -0.01, b3 = -0.3, b4 = 0.78), significant = FALSE) # a1p, a4, a5, a6, a7, a8, b1p, b2p, b3p, b4p not significant but a4 significant if non-significant a8 is present, step factor with a4p, step factor or NaN-inf with nlrob()
otherHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), other2016, start = list(a1 = 30, b1 = 0.23, b2 = -0.023, b3 = -0.25, b4 = 0.9, b4p = -0.08)) # a1p, b1p, b2p, b3p not significant, b4p debatable, intermittent NaN-inf with nlrob()
otherHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, other2016, start = list(a1 = 100, a2 = 0.5, b1 = 0.0, b2 = -0.13, b3 = -0.3, b3p = 0.03, b4 = 0.8), control = gsl_nls_control(maxiter = 500, xtol = 1E-6)) # a1p, a2, a2p, b1, b1p, b2p, b4p not significant, b3p debatable, step factor with nlrob()
otherHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), other2016, start = list(a1 = 0.8, b1 = 1.1, b2 = -0.05)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), other2016, start = list(a1 = 75, b1 = -0.02, b2 = 0.85)) # a1p, b1p, b2p not significant
otherHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp(b1*DBH^b2)), other2016, start = list(a1 = 100, a2 = 0.6, b1 = -0.01, b2 = 0.8), significant = FALSE) # a1p, a2, a2p, a3, a3p, b1p, b2p not significant, step factor with nlrob()
otherHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), other2016, start = list(a1 = -1.1, a3 = 0.3, a4 = 50, b1 = -0.1, b2 = 0.6), control = nls.control(maxiter = 250)) # a2, a3, a3p, a4p, b2p not significant, a4 debatable
#confint_nlrob(otherHeightFromDiameter$sharmaPartonBalPhysio, level = 0.99, weights = pmin(other2016physio$DBH^if_else(other2016physio$isPlantation, -1.0, -1.9), 1))

if (htDiaOptions$includeInvestigatory)
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
if (htDiaOptions$fitGnls)
{
  otherHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + a1*(1 - exp((b1 + b1p*isPlantation)*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.01), folds = 1, repetitions = 1)) # step halving at nlsTol = 1 with corSymm
  otherHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, other2016, start = otherHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # no convergence wtih corSymm
  otherHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, other2016, start = otherHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4), folds = 1, repetitions = 1) # corSymm viable but dropped
  otherHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, other2016, start = otherHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4), folds = 1, repetitions = 1) #  # corSymm viable but dropped, step halving at nlsTol = 0.001
  otherHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), other2016, start = otherHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-4, tolerance = 1E-3, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) #  # corSymm viable but dropped, step halving at nlsTol = 0.005
  #otherHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^b4, other2016, start = otherHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001), folds = 1, repetitions = 1) # step halving with plot level correlation, nlminb() NaN without
  otherHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4), folds = 1, repetitions = 1)  # corSymm viable but dropped
  otherHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH^(b2 + b2p * isPlantation))), other2016, start = otherHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.002), folds = 1, repetitions = 1) # step halving with plot level correlation, with nlsTol = 0.001 without plot level correlation

  save(otherHeightFromDiameterGnls, file = "trees/height-diameter/data/other GNLS.rdata")
}

if (htDiaOptions$includeInvestigatory)
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
otherDiameterFromHeight = list(chapmanReplace = fit_nlrob("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 3.4, b1 = 0.476, b2 = 0.224))) # a1p, b1p, b2p not significant, NaN-inf with nls(), no convergence from nls_multstart()
otherDiameterFromHeight$chapmanReplaceAat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = list(a1 = 6.85, a2 = 0.019, b1 = 0.13, b2 = 0.52)) # occasional step factor with nlrob()
otherDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.0, a3 = 0.003, b1 = 1.25, b2 = 0.33), significant = FALSE) # a1, a2, a3, b2 not significant -> drop BAL on AIC -> no coefficients significant
otherDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 0.65, a2 = -0.03, a3 = 0.03, a9 = 0.1, b1 = 1.5, b2 = 0.27)) # step factor with nls(), step factor or singular gradient with nlrob()
otherDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.0, a9 = 0.0, b1 = 1.3, b2 = 0.35), control = nls.control(maxiter = 150)) # step factor with nls()
otherDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -5, b1 = 0.4, b2 = 0.3, b2p = -0.004), control = nls.control(maxiter = 250)) # a1p not significant, intermittent step factor with nlrob()
otherDiameterFromHeight$chapmanRichardsAat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -7.3, a2 = -0.01, b1 = 0.3, b2 = 0.4, b2p = -0.02), control = gsl_nls_control(maxiter = 500)) # a1p, a2p not significant, occasional step factor with nlrob()
otherDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = list(a1 = -5, a1p = 0.7, a8 = 0.02, b1 = 0.5, b2 = 0.23), maxit = 80, control = nls.control(maxiter = 500), significant = FALSE) # no physiographic effect significant, a1p debatable, convergence fails with b1p, unreliable nlrob() convergence
otherDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = list(a1 = -21.3, a9 = 1.5, b1 = 0.11, b2 = 0.60, b2p = 0.006), control = gsl_nls_control(maxiter = 250), significant = FALSE) # a1p, a2p, a9 not significant, fit_nlrob() fails to converge from closer positions
otherDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = other2016gamConstraint), data = other2016)
otherDiameterFromHeight$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 15, pc = other2016gamConstraint), data = other2016)
otherDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 26, pc = other2016gamConstraint), data = other2016physio)
otherDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 22, pc = other2016gamConstraint), data = other2016physio)
otherDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), other2016)
otherDiameterFromHeight$michaelisMentenReplace = fit_nlrob("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = list(a1 = 17, a2 = 7, b1 = 0.4), control = nls.control(maxiter = 250)) # a1p, a2p, b1p not significant
otherDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ a1*sqrt(TotalHt - 1.37) / (1 + a2*sqrt(TotalHt - 1.37)), other2016, start = list(a1 = 2.3, a2 = -0.13), control = gsl_nls_control(maxiter = 250)) # converges poorly with both a1p and a2p, can yield negative values, a2p dropped on intermittent step factors, prone to >500 iterations with nlrob()
otherDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2), other2016) # isPlantation*(TotalHt - 1.37)^2 not significant
otherDiameterFromHeight$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 3.15, b1 = 0.4)) # a1p, b1p not significant
otherDiameterFromHeight$powerAat = fit_nlrob("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 2.8, a2 = 0.02, b1 = 0.4)) # a2p, b2p not significant, nlrob() failure to converge >500 steps with a1p
otherDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1, other2016physio, start = list(a1 = 2.7, a4 = -0.001, b1 = 0.8), significant = FALSE) # a1p, a4, a5, a6, a7, a8 not significant
otherDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, other2016, start = list(a1 = 3.0, a9 = 2.3, b1 = 0.4), maxit = 80, significant = FALSE) # a1p, a9, b1p not significant
otherDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 2.6, b1 = 0.4, b1p = -0.2, b2 = 0.06, b2p = 0.04)) # step factor with a1p
#otherDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = list(a1 = 0.00003, a2 = 0.01, b1 = 1.20, Ha = 50), control = gsl_nls_control(maxiter = 100)) # NaN-inf with fit_nlrob(), NaNs with gsl_nls()
otherDiameterFromHeight$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1), other2016, start = list(a1 = 32, b1 = -0.6, b2 = 0.1, b3 = -0.07), control = nls.control(maxiter = 1000)) # gsl_nl() parameter evaporation with b4 allowed to vary from 1, step size or singular gradient with nls(), NaN-inf with nlrob()
otherDiameterFromHeight$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = list(a1 = 0.67, a1p = 1.84, b1 = 1.70, b1p = -1.26, b2 = -0.063, b2p = 0.30))
otherDiameterFromHeight$sibbesenReplaceAat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), other2016, start = list(a1 = 2.7, a2 = 0.01, b1 = 0.28, b2 = 0.34), control = nls.control(tol = 1E-3), significant = FALSE) # a1p, a2, a3, b1p, b2p not significant
otherDiameterFromHeight$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = list(a1 = 3.5, a1p = -0.4, a8 = -0.01, b1 = 0.3, b2 = 0.33, b2p = 0.02), significant = FALSE) # b1p, no physiographic predictor significant, b2p debatable
otherDiameterFromHeight$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), other2016, start = list(a1 = 3.0, a9 = -0.7, b1 = 0.3, b2 = 0.348), significant = FALSE) # a1p, a2p, a9, b1p, b2p not significant
otherDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = list(a1 = -133, a1p = 40, b1 = 0.09, b2 = 0.5)) # NaN-inf with b1p, occasional step factor with nlrob()
#confint_nlrob(otherDiameterFromHeight$weibull, level = 0.99, weights = other2016$dbhWeight)

if (htDiaOptions$includeInvestigatory)
{
  print(otherDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(other2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanReplaceAat), y = TotalHt, color = "Chapman-Richards replace approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    geom_line(aes(x = predict(otherDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$schnute), y = TotalHt, color = "Schnute inverse", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen replace", group = isPlantation)) +
    #geom_line(aes(x = predict(otherDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
    #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace aBAL", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 5*(TotalHt - 1.37)^0.5*(exp(0.01*(tph/topHeight)^0.26*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 200*(TotalHt - 1.37)^0.9/(100 - (TotalHt - 1.37)^0.9), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 90, label = "other species, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  #ggplot(other2016) +
  #  geom_point(aes(x = TotalHt, y = abs(residuals(otherDiameterFromHeight$michaelisMentenReplace)), color = "Michaelis-Menten replace", group = isPlantation), alpha = 0.1, color = "grey25", shape= 16)
}
 
 
## other diameter-height GNLS regressions
if (htDiaOptions$fitGnls)
{
  otherDiameterFromHeightGnls = list(chapmanReplace = fit_gnls("Chapman-Richards replace GNLS", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeight$chapmanReplace$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()))
  otherDiameterFromHeightGnls$chapmanReplaceAat = fit_gnls("Chapman-Richards replace ABA+T GNLS", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, other2016, start = otherDiameterFromHeight$chapmanReplaceAat$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$chapmanReplaceBal = fit_gnls("Chapman-Richards replace BA+L GNLS", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanReplaceBal$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  otherDiameterFromHeightGnls$chapmanReplaceBalRelHt = fit_gnls("Chapman-Richards replace BA+L RelHt GNLS", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanReplaceBalRelHt$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  otherDiameterFromHeightGnls$chapmanReplaceRelHt = fit_gnls("Chapman-Richards replace RelHt GNLS", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = otherDiameterFromHeight$chapmanReplaceRelHt$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  otherDiameterFromHeightGnls$chapmanRichards = fit_gnls("Chapman-Richards inverse GNLS", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichards$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  otherDiameterFromHeightGnls$chapmanRichardsAat = fit_gnls("Chapman-Richards inverse ABA+T GNLS", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichardsAat$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  otherDiameterFromHeightGnls$chapmanRichardsPhysio = fit_gnls("Chapman-Richards inverse physio GNLS", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), other2016physio, start = otherDiameterFromHeight$chapmanRichardsPhysio$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl(maxIter = 100))
  otherDiameterFromHeightGnls$chapmanRichardsRelHt = fit_gnls("Chapman-Richards inverse RelHt GNLS", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), other2016, start = otherDiameterFromHeight$chapmanRichardsRelHt$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl(maxIter = 100, msMaxIter = 90))
  otherDiameterFromHeightGnls$michaelisMentenReplace = fit_gnls("Michaelis-Menten replace GNLS", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), other2016, start = otherDiameterFromHeight$michaelisMentenReplace$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl(nlsTol = 0.05)) # step halving at nlsTol = 0.02
  otherDiameterFromHeightGnls$naslund = fit_gnls("Näslund GNLS", DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), other2016, start = otherDiameterFromHeight$naslund$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$power = fit_gnls("power GNLS", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeight$power$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$powerAat = fit_gnls("power ABA+T GNLS", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, other2016, start = otherDiameterFromHeight$powerAat$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$powerPhysio = fit_gnls("power physio GNLS", DBH ~ (a1 + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, other2016physio, start = otherDiameterFromHeight$powerPhysio$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$powerRelHt = fit_gnls("power RelHt GNLS", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), other2016, start = otherDiameterFromHeight$powerRelHt$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$ruark = fit_gnls("Ruark GNLS", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = otherDiameterFromHeight$ruark$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$schnute = fit_gnls("Schnute inverse GNLS", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), other2016, start = otherDiameterFromHeight$schnute$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation)) # step halving
  otherDiameterFromHeightGnls$sharmaParton = fit_gnls("modified Sharma-Parton GNLS", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1), other2016, start = otherDiameterFromHeight$sharmaParton$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl()) # NaN
  otherDiameterFromHeightGnls$sibbesenReplace = fit_gnls("Sibbesen replace GNLS", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenReplace$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$sibbesenReplaceAat = fit_gnls("Sibbesen replace ABA+T GNLS", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenReplaceAat$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$sibbesenReplacePhysio = fit_gnls("Sibbesen replace physio GNLS", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = otherDiameterFromHeight$sibbesenReplacePhysio$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation))
  otherDiameterFromHeightGnls$sibbesenReplaceRelHt = fit_gnls("Sibbesen replace RelHt GNLS", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016, start = otherDiameterFromHeight$sibbesenReplaceRelHt$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation)) # intermittent step halving
  otherDiameterFromHeightGnls$weibull = fit_gnls("Weibull inverse GNLS", DBH ~ ((a1 + a1p*isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, other2016, start = otherDiameterFromHeight$weibull$fit[[1]]$m$getPars(), weights = varPower(1.0, ~TotalHt | isPlantation), control = gnlsControl())
  
  if (exists("otherHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/other GNLS.rdata") }
  save(file = "trees/height-diameter/data/other GNLS.rdata", otherHeightFromDiameterGnls, otherDiameterFromHeightGnls)
}
if (htDiaOptions$includeInvestigatory)
{
  ggplot() +
    geom_histogram(aes(x = power), otherDiameterFromHeightResultsGnls, binwidth = 0.1) +
    geom_segment(aes(x = mean(otherDiameterFromHeightResultsGnls$power), xend = mean(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash") +
    geom_segment(aes(x = median(otherDiameterFromHeightResultsGnls$power), xend = median(otherDiameterFromHeightResultsGnls$power), y = 0, yend = 10), color = "grey70", linetype = "longdash")
    
  otherDiameterFromHeightResultsGnls %>% summarize(power = mean(power), powerPlantation = mean(powerPlantation))
}

  
## collect model parameters
#if (exists("otherHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/other GNLS.rdata") }
otherCoefficients = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, get_list_coefficients))) %>%
                                        #bind_rows(lapply(otherHeightFromDiameterGnls, get_model_coefficients))) %>%
                                        #get_model_coefficients(otherHeightFromDiameter$sharmaZhangBalGnls),
                                mutate(responseVariable = "height"),
                              bind_rows(bind_rows(lapply(otherDiameterFromHeight, get_list_coefficients))) %>%
                                        #bind_rows(lapply(otherDiameterFromHeightGnls, get_model_coefficients))) %>%
                                        #get_model_coefficients(otherDiameterFromHeight$schnuteGnls),
                                mutate(responseVariable = "DBH")) %>%
  mutate(species = "other")

otherResults = bind_rows(bind_rows(bind_rows(lapply(otherHeightFromDiameter, get_list_stats)),
                                   #bind_rows(lapply(otherHeightFromDiameterGnls, get_model_stats)),
                                   get_model_stats(name = "Sharma-Zhang BA+L GNLS", fitting = "gnls")) %>%
                           mutate(responseVariable = "height"),
                         bind_rows(bind_rows(lapply(otherDiameterFromHeight, get_list_stats)),
                                   get_model_stats(name = "Schnute inverse", fitting = "gsl_nls")) %>%
                                   #bind_rows(lapply(otherDiameterFromHeightGnls, get_model_stats)),
                                   #get_model_stats(name = "Schnute GNLS", fitting = "gnls")) %>%
                           mutate(responseVariable = "DBH")) %>%
  mutate(species = "other")

save(file = "trees/height-diameter/data/other results.Rdata", otherCoefficients, otherResults)
if (htDiaOptions$folds * htDiaOptions$repetitions <= htDiaOptions$retainModelThreshold)
{
  save(file = "trees/height-diameter/data/other models and stats.Rdata", otherHeightFromDiameter, otherDiameterFromHeight)
} else {
  save(file = "trees/height-diameter/data/other stats.Rdata", otherHeightFromDiameter, otherDiameterFromHeight)
}


## preferred forms identified (results.R, Figure 8)
# other species  Curtis            REML GAM BA+L                  REML GAM                Sibbesen form physio
#                power             REML GAM BAL+L physio          parabolic               Chapman-Richards form RelHt
#                Korf                                             linear
otherHeightFromDiameterPreferred = list(curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1*DBH / (1 + DBH)^b1, other2016, start = list(a1 = 1.086, b1 = 0.190), folds = 1, repetitions = 1))
otherHeightFromDiameterPreferred$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 15, pc = other2016gamConstraint), data = other2016, folds = 1, repetitions = 1)
otherHeightFromDiameterPreferred$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 23, pc = other2016gamConstraint), data = other2016physio, folds = 1, repetitions = 1)
otherHeightFromDiameterPreferred$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), other2016, start = list(a1 = 381, b1 = -6.26, b2 = -0.20), folds = 1, repetitions = 1)
otherHeightFromDiameterPreferred$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, other2016, start = list(a1 = 0.99, b1 = 0.84), folds = 1, repetitions = 1)

otherDiameterFromHeightPreferred = list(chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), other2016, start = list(a1 = 1.0, a9 = 0.0, b1 = 1.3, b2 = 0.35), control = nls.control(maxiter = 150), folds = 1, repetitions = 1))
otherDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = other2016gamConstraint), data = other2016, folds = 1, repetitions = 1)
otherDiameterFromHeightPreferred$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), other2016, folds = 1, repetitions = 1)
otherDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2), other2016, folds = 1, repetitions = 1)
otherDiameterFromHeightPreferred$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), other2016physio, start = list(a1 = 3.47, a1p = -0.43, a8 = -0.001, b1 = 0.34, b2 = 0.28, b2p = 0.019), significant = FALSE, folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/other preferred models.Rdata", otherHeightFromDiameterPreferred, otherDiameterFromHeightPreferred)


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  otherBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), other2016, start = list(a1 = 1.36, b1 = 0.0002, b2 = 2.06, b2p = -0.27), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor with fit_nlrob()
  otherBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p*isPlantation), other2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = -0.00001, b1 = 2.53, b1p = -0.435), maxit = 70, weights = pmin(1/basalArea, 1E4))
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