# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## bigleaf maple height-diameter regression form sweep
# preferred forms: Sharma-Parton BA+L, Ratkowsky, Sharma-Parton, Hossfeld, Weibull, Chapman-Richards
#acmaHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 31.9, a2 = -0.0044, a3 = -8.22, a4 = 0.198, a5 = 1.387, a6 = 0.014, b1 = -0.031, b2 = 1.02, b2p = -0.108)) # a1p, a2, a4, b1p not significant
#acmaHeightFromDiameter$michaelisMenten = fit_gsl_nls(TotalHt ~ 1.37 + a1*DBH / (a2 + a2p * isPlantation + DBH), acma2016, start = list(a1 = 39.5, a2 = 47.8, a2p = -9.11))
#acmaHeightFromDiameter$richards = fit_gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), acma2016, start = list(Ha = 23.9, Hap = -1.1, d = 0.67, dp = -0.031, kU = 0.023, kUp = 0.008))
acma2016 = trees2016 %>% filter(Species == "BM", isLiveUnbroken, TotalHt > 0) %>% # live bigleaf maples measured for height
  mutate(dbhWeight = pmin(1/(1.19*DBH^0.77), 5),
         heightWeight = pmin(1/(9.70*(TotalHt - 1.37)^1.20), 5))
acma2016physio = acma2016 %>% filter(is.na(elevation) == FALSE)
acma2016gamConstraint = c(DBH = -1.7882/0.7564, TotalHt = 1.37, standBasalAreaPerHectare = median(acma2016$standBasalAreaPerHectare), basalAreaLarger = median(acma2016$basalAreaLarger), standBasalAreaApprox = median(acma2016$standBasalAreaApprox), tallerApproxBasalArea = median(acma2016$tallerApproxBasalArea), elevation = median(acma2016physio$elevation), slope = median(acma2016physio$slope), aspect = median(acma2016physio$aspect), topographicShelterIndex = median(acma2016physio$topographicShelterIndex), relativeHeight = median(acma2016$relativeHeight)) # point constraint for mgcv::s()
#acma2016natural = acma2016 %>% filter(isPlantation == FALSE)
#acma2016plantation = acma2016 %>% filter(isPlantation)
#acma2016plantationPhysio = acma2016physio %>% filter(isPlantation)

acmaHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^b2, acma2016, start = list(a1 = 26.8, a1p = 4.20, b1 = -0.026, b2 = 0.927))) # b1p, b2p not significant
acmaHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29)) # job step factor with nlrob()
acmaHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a5 * slope) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016physio, start = list(a1 = 30, a2 = 0.10, a3 = -0.05, a3p = 0.05, a5 = -0.1, b1 = -0.03, b2 = 1.0, b2p = -0.12)) # a1p, a2, a2p, a3, a4, a7, a8, b1p not significant, a5-a6 not mutually significant, a3p debatable, job step size with nlrob()
acmaHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = -2.4, a1p = 1.51, a2 = -0.023, a2p = 0.159, a3 = 0.045, a4 = 59.4, a4p = -25.9, b1 = -0.033, b2 = 0.018, b2p = 0.406)) # a2, a3p not significant, NaN-inf with b1p, job step factor with nlrob()
acmaHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 32.5, a5 = -8.37, a7 = 1.557, b1 = -0.031, b2 = 1.01, b2p = -0.1-4)) # a1p, a2, a4, a6, a8, b1p not significant
acmaHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^b1, acma2016, start = list(a1 = 1.24, a1p = 0.23, b1 = 0.28)) # b1p not significant
acmaHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 7, pc = acma2016gamConstraint), data = acma2016)
acmaHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 14, pc = acma2016gamConstraint), data = acma2016)
acmaHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, basalAreaLarger, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 57, pc = acma2016gamConstraint), data = acma2016physio) # more data than coefficients -> eliminate topographic shelter AIC 5284: 5288 without cos(aspect), 5273 without sin(aspect), 5261 without slope, 5260 without elevation, 5304 without BAL, 5264 without BA -> eliminate elevation AIC 5260: 5244 without cos(aspect), 5255 without sin(aspect), 5242 without slope, 5255 without BAL, 5228 without BA -> eliminate BA
acmaHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), bs = "ts", by = as.factor(isPlantation), k = 57, pc = acma2016gamConstraint), data = acma2016physio) # AIC 5294 with all physiographic predictors: 5238 without elevation, 5262 without slope, 5271 without sin(aspect), 5256 without cos(aspect), 5248 without topographic shelter -> eliminate shelter AIC 5250: 5253 without elevation, 5262 without slope, 5261 without sin(aspect), 5253 without cos(aspect)
acmaHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + b1 * DBH^b2), acma2016, start = list(a1 = 40.1, a1p = 6.30, b1 = 44.9, b2 = -0.97)) # b1p, b2p not significant
acmaHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), acma2016, start = list(a1 = 70, b1 = -5.0, b2 = -0.35, b2p = -0.02)) # a1p, b1p, b2p not mutually significant
acmaHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), acma2016)
acmaHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986)) # a1p, b1p not significant
acmaHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), acma2016)
acmaHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1 * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), acma2016, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19)) # a1p, a3p not significant
acmaHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1, acma2016, start = list(a1 = 1.11, a1p = 0.21, b1 = 0.75)) # b1p not significant
acmaHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1/(DBH + b2 + b2p * isPlantation)), acma2016, start = list(a1 = 31.8, a1p = 3.00, b1 = -25.7, b2 = 6.70, b2p = 0.62)) # b1p not significant
acmaHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), acma2016, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064)) # Hap, dp not significant
acmaHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 10, a1p = -4, b1 = 0.3, b1p = 0.1, b2 = -0.03, b3 = 0.06, b4 = 1.05, b4p = -0.14)) # b2p, b3p not significant, job step factor with nlrob()
acmaHeightFromDiameter$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 20, b1 = 0.1, b2 = -0.04, b3 = -0.1, b4 = 1.1, b4p = -0.24)) # a1p, b1p, b2p, b3p not significant
acmaHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, acma2016physio, start = list(a1 = 12, a6 = 0.4, b1 = 0.2, b1p = 0.04, b2 = -0.037, b3 = -0.046, b4 = 1.06)) # a1p, a4, a5, a6, a7, a8, b1, b2p, b3p, b4p not significant, job step size with nlrob()
acmaHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a6 * sin(3.14159/180 * aspect))*topHeight^(b1 + b1p * isPlantation) * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, acma2016physio, start = list(a1 = 11, a6 = 0.3, b1 = 0.2, b1p = 0.06, b2 = -0.03, b3 = 0, b4 = 1.0), significant = FALSE) # a4, a5, a6, a7, a8, b1, b2p, b3p, b4p not significant, job step size with nlrob()
acmaHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^b4, acma2016, start = list(a1 = 17, a1p = 2.0, b1 = 0.1, b2 = -0.05, b3 = -0.1, b4 = 1.0)) # b1p, b2p, b3p not significant
acmaHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = list(a1 = 22, a2 = 0.05, b1 = 0.05, b2 = -0.04, b3 = -0.04, b4 = 1.07, b4p = -0.2), significant = FALSE) # a1p, b1p, a2, a2p, b2p, b3p not significant, job step factor with nlrob()
acmaHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^b2), acma2016, start = list(a1 = 0.752, a1p = 0.120, b1 = 1.180, b2 = -0.087)) # b1p, b2p not significant
acmaHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^b2)), acma2016, start = list(a1 = 27, a1p = 3, b1 = -0.033, b2 = 0.94)) # b1p, b2p not significant
acmaHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = 26, a2 = 0.1, b1 = -0.027, b1p = -0.01, b2 = 0.99)) # a1p, a2p, a3p, b2p not significant, a2+a3 not mutually significant
acmaHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = 100, a2 = -20, a4 = -100, b1 = 0.03, b1p = -0.001, b2 = 0.33), control = gsl_nls_control(maxiter = 250, xtol = 0.01)) # a3, a4p, b1p not significant, a1-a4 evaporation: step factor with a1p, a2p, a3p, b2p, b3p
#lapply(acmaHeightFromDiameter$sharmaPartonBalPhysio$fit, confint_nlrob)

if (htDiaOptions$includeInvestigatory)
{
  print(acmaHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = acma2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = acma2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$curtis), color = "Curtis", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$korf), color = "Korf", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$linear), color = "linear", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$parabolic), color = "parabolic", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$power), color = "power", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$prodan), color = "Prodan", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$richards), color = "unified Richards", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$sibbesen), color = "Sibbesen", group = acma2016$isPlantation)) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$weibull), color = "Weibull", group = acma2016$isPlantation)) +
    annotate("text", x = 0, y = 40, label = "bigleaf maple, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 40)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
}


## bigleaf maple height-diameter GNLS regressions
if (htDiaOptions$fitGnls)
{
  acmaHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^b2, acma2016, start = acmaHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1)) # corSymm viable but dropped
  acmaHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), acma2016, start = acmaHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-4, tolerance = 1E-3, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) #  # corSymm viable but dropped, step halving at nlsTol = 0.001, maxiter at default msTol and tolerance
  acmaHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation)*(1 - exp(b2*(tph/standBasalAreaPerHectare)^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = acmaHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.002), folds = 1, repetitions = 1) # foreign NaN-inf with plot correlation, step halving at nlsTol = 0.001
  acmaHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = acmaHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # corSymm viable but dropped, step halving at nlsTol = 0.01
  acmaHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, acma2016, start = acmaHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # corSymm viable but dropped, step halving at nlsTol = 0.005
  acmaHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), acma2016, start = acmaHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50), folds = 1, repetitions = 1) # corSymm viable but dropped, step halving at nlsTol = 0.005
  acmaHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp(b1*DBH^b2)), acma2016, start = acmaHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5), folds = 1, repetitions = 1) # corSymm viable but dropped
  acmaHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + a2*basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = acmaHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5), folds = 1, repetitions = 1) # corSymm viable but dropped

  save(acmaHeightFromDiameterGnls, file = "trees/height-diameter/data/ACMA3 GNLS.rdata")
}

if (htDiaOptions$includeInvestigatory)
{
  acmaHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(acmaHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(acmaHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(acmaHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
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
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$weibullBAL), color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = acma2016natural$DBH, y = predict(acmaHeightFromDiameter$weibullBALnatural), color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = acma2016plantation$DBH, y = predict(acmaHeightFromDiameter$weibullBALplantation), color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterBase), color = "base")) +
    geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameter$weibull), color = "ElliottWeibull")) +
    annotate("text", x = 0, y = 85, label = "a) bigleaf maple, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## bigleaf maple diameter-height regressions
acmaDiameterFromHeight = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 540000, b1 = 0.00001, b2 = 1.11))) # NaN-inf from nls() at multiple nls_multstart() positions, NaN-inf from fit_nlrob()
acmaDiameterFromHeight$chapmanReplaceAat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, acma2016, start = list(a1 = 46000, a2 = 159, b1 = 0.000013, b2 = 1.11)) # NaN-inf from nls() and fit_nlrob()
acmaDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 23600, a2 = -1955, b1 = 0.00001, b2 = 1.07)) # NaN-inf from nls(), step factor from fit_nlrob()
acmaDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 261000, a2 = -7700, a3 = 4881, a9 = -216000, b1 = 0.00001, b2 = 1.06)) # step factor from nls()
acmaDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), acma2016, start = list(a1 = 241000, a9 = -114000, b1 = 0.000007, b2 = 1.23)) # step factor from nls() and fit_nlrob()
acmaDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), acma2016, start = list(a1 = 23.5, b1 = -0.022, b2 = 2.01, b2p = -0.17)) # b1p not significant
acmaDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), acma2016, start = list(a1 = 24.2, a2 = -0.005, b1 = -0.021, b2 = 2.02, b2p = -0.19), control = list(maxiter = 500))
acmaDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), acma2016physio, start = list(a1 = 20, a8 = 0.2, b1 = -0.03, b1p = 0.01, b2 = 1.8), significant = FALSE) # a4, a5, a6, a7, a8, b2p not significant, a1p-b1p not mutually significant 
acmaDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), acma2016, start = list(a1 = 36.1, a9 = -2.78, b1 = -0.028, b1p = 0.009, b2 = 1.63)) # a1p, a2p, b2p not significant
acmaDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = acma2016gamConstraint), data = acma2016)
acmaDiameterFromHeight$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 13, pc = acma2016gamConstraint), data = acma2016)
acmaDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = acma2016gamConstraint), data = acma2016physio) # insufficient data with all predictors -> eliminate elevation AIC 7849: >7900 with any other predictor removed but not always enough data for cross validation without additional removal -> drop cos(aspect) since AIC 7885 is the least increase
acmaDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 57, pc = acma2016gamConstraint), data = acma2016physio) # AIC 7971: 7967 without elevation, 8041 without slope, 7992 without sin(aspect), 7993 without cos(aspect), 8004 without topographic shelter -> eliminate elevation
acmaDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), acma2016)
acmaDiameterFromHeight$michaelisMentenReplace = fit_nlrob("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), acma2016, start = list(a1 = -116, a2 = -128, b1 = 1.50))
acmaDiameterFromHeight$naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), acma2016, start = list(a1 = 5.1, a1p = -2.0, a2 = -0.12, a2p = -0.023))
acmaDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), acma2016) # (TotalHt - 1.37)^2 not significant (p = 0.058)
acmaDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282))
acmaDiameterFromHeight$powerAat = fit_nlrob("power ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.63, a1p = -2.32, a2 = -0.00064, b1 = 0.898, b1p = 0.272)) # a2p not significant
acmaDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, acma2016physio, start = list(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9)) # a4, a6, a7, b1p not significant, a5, a8 debatable
acmaDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1, acma2016, start = list(a1 = 2.5, a9 = -0.6, b1 = 1.0), significant = FALSE) # a1p, a9, a9p not significant
acmaDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037)) # a1p not significant
acmaDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), acma2016, start = list(a1 = 0.000003, a2 = 0.002, b1 = 1.13, Ha = 161)) # singular gradient or step factor with fit_nlrob()
acmaDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(tph/topHeight)^b3*(TotalHt - 1.37)) - 1)^b4, acma2016, start = list(a1 = 3.3, b1 = 0.9, b2 = 0.001, b3 = 1, b4 = 0.5), control = gsl_nls_control(maxiter = 500, xtol = 0.002)) # NaN-inf from nls() at nls_multstart() positions, NaN-inf or code error with fit_nlrob(), a1 -> 0 with gsl_nls() -> eliminate b4
acmaDiameterFromHeight$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 6, b1 = 2, b2 = -0.13)) # no significant plantation effects
acmaDiameterFromHeight$sibbesenReplaceAat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.5, a2 = 0.0002, b1 = 2.2, b2 = -0.1), significant = FALSE) # a2 not significant
acmaDiameterFromHeight$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a5 * sin(3.14159/180 * slope))*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), acma2016physio, start = list(a1 = 0.5, a5 = -0.2, b1 = 2.6, b2 = -0.18, b2p = 0.02), significant = FALSE) # no physiographic parameters significant AIC 7505 -> drop topographic shelter AIC 7494 -> drop sin(aspect) AIC 7498 -> drop cos(aspect) AIC 7508 -> drop elevation AIC 7505, a1, a5 not significant, a1p+b1p not mutually significant
acmaDiameterFromHeight$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), acma2016, start = list(a1 = 0.65, a9 = -0.25, b1 = 2.0, b2 = -0.1))
acmaDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, acma2016, start = list(a1 = -185000, b1 = 0.00001, b2 = 1.05), control = gsl_nls_control(maxiter = 250, xtol = 1E-7)) # a1-b1 evaporation since unfavorable to concave up curvature, NaN-inf with fit_nlrob()
#lapply(acmaDiameterFromHeight$sibbesenReplacePhysio$fit, confint_nlrob)

if (htDiaOptions$includeInvestigatory)
{
  print(acmaDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(acma2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$chapmanReplaceAat), y = TotalHt, color = "Chapman-Richards approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen", group = isPlantation)) +
    #geom_line(aes(x = predict(acmaDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards replace"), na.rm = TRUE) +
    #geom_line(aes(x = 1*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.12*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -116*(TotalHt - 1.37)^1.5/(-128 - (TotalHt - 1.37)^1.5), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation), alpha = 0.5) +#geom_line(aes(x = predict(acmaDiameterFromHeight$schnute), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = -1/0.0003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute inverse"), alpha = 0.5) +
    #geom_line(aes(x = 1*(TotalHt - 1.37)^1.5*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^0.5*(TotalHt - 1.37)))^0.9, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = (-15 * log(1 - pmin(0.01*(TotalHt - 1.37), 0.999)))^-0.5, y = TotalHt, color = "Weibull"), na.rm = TRUE) +
    annotate("text", x = 0, y = 41, label = "bigleaf maple, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


## collect model parameters
#if (exists("acmaHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/ACMA3 GNLS.rdata") }
acmaCoefficients = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_coefficients))) %>%
                                       #bind_rows(lapply(acmaHeightFromDiameterGnls, get_model_coefficients))) %>%
                               mutate(responseVariable = "height"),
                             bind_rows(lapply(acmaDiameterFromHeight, get_list_coefficients)) %>%
                               mutate(responseVariable = "DBH")) %>%
  mutate(species = "ACMA3")

acmaResults = bind_rows(bind_rows(bind_rows(lapply(acmaHeightFromDiameter, get_list_stats))) %>%
                                  #bind_rows(lapply(acmaHeightFromDiameterGnls, get_model_stats))) %>%
                          mutate(responseVariable = "height"),
                        bind_rows(lapply(acmaDiameterFromHeight, get_list_stats)) %>%
                            mutate(responseVariable = "DBH")) %>%
  mutate(species = "ACMA3")

save(file = "trees/height-diameter/data/ACMA3 results.Rdata", acmaCoefficients, acmaResults)
if (htDiaOptions$folds * htDiaOptions$repetitions <= htDiaOptions$retainModelThreshold)
{
  save(file = "trees/height-diameter/data/ACMA3 models and stats.Rdata", acmaHeightFromDiameter, acmaDiameterFromHeight)
} else {
  save(file = "trees/height-diameter/data/ACMA3 stats.Rdata", acmaHeightFromDiameter, acmaDiameterFromHeight)
}


## preferred forms identified (results.R, Figure 6)
acmaHeightFromDiameterPreferred = list(chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a5 * sin(3.14159/180 * slope) + a7 * sin(3.14159/180 * aspect)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), acma2016, start = list(a1 = 32.5, a5 = -8.37, a7 = 1.557, b1 = -0.031, b2 = 1.01, b2p = -0.1-4), folds = 1, repetitions = 1))
acmaHeightFromDiameterPreferred$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), folds = 1, repetitions = 1)
acmaHeightFromDiameterPreferred$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1 * DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), acma2016, start = list(a1 = 0.024, a2 = 1.27, a2p = -0.23, a3 = -0.19), folds = 1, repetitions = 1)
acmaHeightFromDiameterPreferred$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), acma2016, start = list(Ha = 22.8, d = 0.723, kU = 0.025, kUp = 0.0064), folds = 1, repetitions = 1)
acmaHeightFromDiameterPreferred$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), acma2016, start = list(a1 = 26, a2 = 0.1, b1 = -0.027, b1p = -0.01, b2 = 0.99), folds = 1, repetitions = 1)

acmaDiameterFromHeightPreferred = list(chapmanRichardsPhysio = fit_nlrob("Chapman-Richards inverse physio", DBH ~ (a1 + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), acma2016physio, start = list(a1 = 20, a8 = 0.2, b1 = -0.03, b1p = 0.01, b2 = 1.8), significant = FALSE, folds = 1, repetitions = 1))
acmaDiameterFromHeightPreferred$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = acma2016gamConstraint), data = acma2016, folds = 1, repetitions = 1)
acmaDiameterFromHeightPreferred$power = fit_nlrob("power", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 3.57, a1p = -2.30, b1 = 0.894, b1p = 0.282), folds = 1, repetitions = 1)
acmaDiameterFromHeightPreferred$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, acma2016physio, start = list(a1 = 3.7, a1p = -0.9, a5 = -1.2, a8 = 0.02, b1 = 0.9), folds = 1, repetitions = 1)
acmaDiameterFromHeightPreferred$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/ACMA3 preferred models.Rdata", acmaHeightFromDiameterPreferred, acmaDiameterFromHeightPreferred)


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  #acmaBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), acma2016, start = list(a1 = 1.9, b1 = 0.00007, b2 = 2.3), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # a1p, b1p not significant, fit_nlrob() step factor
  acmaBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p*isPlantation)) - 1), acma2016, start = list(a1 = 690, b1 = 0.00007, b2 = 2.2, b2p = -0.15), weights = pmin(1/basalArea, 1E4), control = nls.control(maxiter = 500)) # fit_nlrob() step factor
  acmaBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), acma2016, start = list(a1 = 2/7 * 0.25 * pi * 0.01^2, a1p = -0.0002, b1 = 2.05, b1p = 0.324), weights = pmin(1/basalArea, 1E4))
  #confint2(acmaBasalAreaFromHeightKorf, level = 0.99)
  #confint_nlrob(acmaBasalAreaFromHeightPower, level = 0.99)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(acmaBasalAreaFromHeightKorf), 100^2 * mean(residuals(acmaBasalAreaFromHeightKorf)), mean(abs(residuals(acmaBasalAreaFromHeightKorf))), 1 - sum(residuals(acmaBasalAreaFromHeightKorf)^2) / sum((acma2016$basalArea - mean(acma2016$basalArea)^2)),
          "power", AIC(acmaBasalAreaFromHeightPower), 100^2 * mean(residuals(acmaBasalAreaFromHeightPower)), mean(abs(residuals(acmaBasalAreaFromHeightPower))), 1 - sum(residuals(acmaBasalAreaFromHeightPower)^2) / sum((acma2016$basalArea - mean(acma2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(acma2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(acmaBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(acmaBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "bigleaf maple height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}
