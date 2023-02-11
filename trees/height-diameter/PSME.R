# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Douglas-fir height-diameter regression form sweep
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, TotalHt > 0) %>% # live Douglas-firs measured for height
  mutate(dbhWeight = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016gamConstraint = c(DBH = -1.2240/0.6566, TotalHt = 1.37, standBasalAreaPerHectare = median(psme2016$standBasalAreaPerHectare), basalAreaLarger = median(psme2016$basalAreaLarger), standBasalAreaApprox = median(psme2016$standBasalAreaApprox), tallerApproxBasalArea = median(psme2016$tallerApproxBasalArea), elevation = median(psme2016physio$elevation), slope = median(psme2016physio$slope), aspect = median(psme2016physio$aspect), topographicShelterIndex = median(psme2016physio$topographicShelterIndex), relativeHeight = median(psme2016$relativeHeight)) # point constraint for mgcv::s() where response variable is ignored, zero crossing of height from DBH from fit_lm(TotalHt ~ DBH, data = psme2016 %>% filter(DBH < 6))

psme2016defaultWeight = psme2016 %>%
  mutate(dbhWeight = pmin(DBH^-1, 1))
psme2016defaultWeightPhysio = psme2016defaultWeight %>% filter(is.na(elevation) == FALSE)
#psme2016natural = psme2016 %>% filter(isPlantation == FALSE)
#psme2016plantation = psme2016 %>% filter(isPlantation)
#psme2016plantationPhysio = psme2016physio %>% filter(isPlantation)

psmeOptions = list(fitHeightPrimary = TRUE, # started 2-16 7:50
                   fitHeightGslNls = FALSE, # done
                   fitDbhPrimary = FALSE, # started 2-16 7:58
                   fitDbhGslNls = FALSE, # done
                   fitSlowGams = TRUE)

if (psmeOptions$fitHeightPrimary)
{
  psmeHeightFromDiameter = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31))) # b1p not significant
  psmeHeightFromDiameter$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 72.9, a1p = -11.8, a2 = 0.087, a2p = 0.84, a3 = -0.0021, a3p = -0.073, b1 = -0.016, b2 = 1.26, b2p = -0.054)) # a3 not significant, step factor with b1p
  psmeHeightFromDiameter$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 74.3, a2 = 0.096, a2p = 0.92, a3 = 0, a3p = 0, a4 = -0.015, a5 = -0.101, a6 = 0.793, a7 = 1.695, a8 = 0.183, b1 = -0.018, b1p = 0.005, b2 = 1.30, b2p = -0.154)) # a4 not significant
  psmeHeightFromDiameter$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a4 = 54.0, a4p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37))
  psmeHeightFromDiameter$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31)) # a4p not significant, a5p induces overfitting
  psmeHeightFromDiameter$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156))
  psmeHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = psme2016gamConstraint), data = psme2016)
  psmeHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 26, pc = psme2016gamConstraint), data = psme2016, nthreads = 4)

  if (psmeOptions$fitSlowGams)
  {
    psmeHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 331, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6) # long fit time (k = 455, edf = 405) AIC 140007: 140237 without BAL, 140249 without BA, 140335 without elevation, 140504 without slope, 140129 without sin(aspect), 140126 without cos(aspect), 140979 without topographic shelter -> force drop of cos(aspect) as least disadvantageous option towards viable fitting times -> 140429 without BA, 140496 without BAL, 140448 without elevation, 140589 without slope, 140279 without sin(aspect), 141178 without topographic shelter
    psmeHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6)
    
    psmeHeightFromDiameterGamBalPhysio = psmeHeightFromDiameter$gamBalPhysio
    psmeHeightFromDiameterGamPhysio = psmeHeightFromDiameter$gamPhysio
    save(file = "trees/height-diameter/data/PSME TotalHt primary GAMs.Rdata", psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamPhysio)
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamPhysio)
  } else {
    load("trees/height-diameter/data/PSME TotalHt primary GAMs.Rdata")
    psmeHeightFromDiameter$gamBalPhysio = psmeHeightFromDiameterGamBalPhysio
    psmeHeightFromDiameter$gamPhysio = psmeHeightFromDiameterGamPhysio
    rm(psmeHeightFromDiameterGamBalPhysio, psmeHeightFromDiameterGamPhysio)
  }
  
  psmeHeightFromDiameter$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28))
  psmeHeightFromDiameter$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), control = nls.control(maxiter = 500)) # a1p parameter evaporation, b1p not significant
  psmeHeightFromDiameter$linear = fit_lm("linear", TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016)
  psmeHeightFromDiameter$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30)) # b1p not significant
  psmeHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), psme2016)
  psmeHeightFromDiameter$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6))
  psmeHeightFromDiameter$power = fit_nlrob("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14))
  psmeHeightFromDiameter$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52))
  psmeHeightFromDiameter$richards = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126))
  psmeHeightFromDiameter$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22))
  psmeHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 20, a1p = -9, b1 = 0.3, b1p = 0.07, b2 = -0.023, b2p = -0.02, b3 = 0.0, b4 = 1.53, b4p = -0.15)) # b3p not significant, job step size with nlrob()
  psmeHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51)) # b1p not significant
  psmeHeightFromDiameter$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 20, a4 = -0.0008, a5 = -0.04, a6 = 0.15, a7 = 0.23, a8 = 0.095, b1 = 0.30, b2 = -0.025, b3 = -0.012, b3p = -0.12, b4 = 1.6, b4p = -0.65)) # a1p, b1p, b2p not significant, b3 debatable, nlrob() step factor with a1p
  psmeHeightFromDiameter$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 54, a1p = -33, b1 = 0.05, b1p = 0.2, b2 = -0.03, b2p = -0.05, b3 = -0.04, b3p = -0.16, b4 = 1.56, b4p = -0.48))
  psmeHeightFromDiameter$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 54, a2 = 0.03, a2p = 0.7, b1 = 0.05, b2 = -0.03, b3 = -0.07, b3p = -0.09, b4 = 1.55, b4p = -0.5)) # a1p, a2, b1p, b2p, b3 not significant
  psmeHeightFromDiameter$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050))
  psmeHeightFromDiameter$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64, a1p = -20, b1 = -0.005, b1p = -0.006, b2 = 1.3, b2p = -0.1))
  psmeHeightFromDiameter$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63, a2 = 0.04, a2p = 0.9, a3 = 0.02, a3p = -0.2, b1 = -0.005, b1p = -0.0026, b2 = 1.3, b2p = -0.15)) # a2, a3 debatably significant
  psmeHeightFromDiameter$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 3.0, a2 = 0.22, a2p = 0.95, a3 = 0.12, a4 = 53, b1 = -0.06, b1p = 0.04, b2 = 0.75, b2p = 0.1)) # a3p, a4p not significant
  #lapply(psmeHeightFromDiameter$sharmaPartonBalPhysio$fit, confint_nlrob)

  save(file = "trees/height-diameter/data/PSME TotalHt primary.rdata", psmeHeightFromDiameter)
}

if (psmeOptions$fitHeightGslNls)
{
  if (exists("psmeHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/PSME TotalHt primary.rdata") }
  
  psmeHeightFromDiameterNls = list(chapmanRichardsNls = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars()))
  psmeHeightFromDiameterNls$chapmanRichardsBalNls = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$chapmanRichardsBalPhysioNls = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$chapmanRichardsPhysioNls = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$curtisNls = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = psmeHeightFromDiameter$curtis$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$hossfeldNls = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$hossfeld$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$korfNls = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$korf$fit[[1]]$m$getPars(), control = nls.control(maxiter = 500))
  psmeHeightFromDiameterNls$michaelisMentenNls = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = psmeHeightFromDiameter$michaelisMenten$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$prodanNls = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = psmeHeightFromDiameter$prodan$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$powerNls = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = psmeHeightFromDiameter$power$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$ratkowskyNls = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$ratkowsky$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$richardsNls = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = psmeHeightFromDiameter$richards$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$sharmaPartonNls = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars())
  #psmeHeightFromDiameterNls$sharmaPartonBalNls = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars())
  #psmeHeightFromDiameterNls$sharmaPartonBalPhysioNls = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$sharmaPartonPhysioNls = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$sharmaZhangNls = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$sharmaZhangBalNls = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$sibbesenNls = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$sibbesen$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$weibullNls = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNls$weibullBalNls = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$fit[[1]]$m$getPars())
  
  psmeHeightFromDiameterNlsDwt = list(chapmanRichardsNlsFwt = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars()))
  psmeHeightFromDiameterNlsDwt$chapmanRichardsBalNlsFwt = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$chapmanRichardsBalPhysioNlsFwt = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$chapmanRichardsPhysioNlsFwt = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$curtisNlsFwt = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$curtis$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$hossfeldNlsFwt = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$hossfeld$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$korfNlsFwt = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$korf$fit[[1]]$m$getPars(), control = nls.control(maxiter = 500))
  psmeHeightFromDiameterNlsDwt$michaelisMentenNlsFwt = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016defaultWeight, start = psmeHeightFromDiameter$michaelisMenten$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$prodanNlsFwt = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$prodan$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$powerNlsFwt = fit_gsl_nls("power", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$power$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$ratkowskyNlsFwt = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$ratkowsky$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$richardsNlsFwt = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$richards$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaPartonNlsFwt = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaPartonBalNlsFwt = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaPartonBalPhysioNlsFwt = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaPartonPhysioNlsFwt = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeightPhysio, start = psmeHeightFromDiameter$sharmaPartonPhysio$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaZhangNlsFwt = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sharmaZhangBalNlsFwt = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016defaultWeight, start = psmeHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$sibbesenNlsFwt = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016defaultWeight, start = psmeHeightFromDiameter$sibbesen$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$weibullNlsFwt = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$weibull$fit[[1]]$m$getPars())
  psmeHeightFromDiameterNlsDwt$weibullBalNlsFwt = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016defaultWeight, start = psmeHeightFromDiameter$weibullBal$fit[[1]]$m$getPars())
  
  save(file = "trees/height-diameter/data/PSME TotalHt gsl_nls.rdata", psmeHeightFromDiameterNls, psmeHeightFromDiameterNlsDwt)
}

if (htDiaOptions$includeInvestigatory)
{
  print(psmeHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot() +
    geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$curtis), color = "Curtis", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$korf), color = "Korf", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$linear), color = "linear", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$parabolic), color = "parabolic", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$power), color = "power", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$prodan), color = "Prodan", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$richards), color = "unified Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sibbesen), color = "Sibbesen", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibull), color = "Weibull", group = psme2016$isPlantation)) +
    annotate("text", x = 0, y = 85, label = "Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
  
  ggplot() +
    geom_point(aes(x = psme2016physio$elevation, y = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio)), alpha = 0.15, color = "grey25", shape = 16) +
    geom_smooth(aes(x = psme2016physio$elevation, y = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio)), alpha = 0.15, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
    labs(x = "physiographic", y = "DBH error, cm")
  
  psmeHeightFromDiameter$Efficiency = psmeHeightFromDiameterResults %>% filter(fitting %in% c("nlrob", "gsl_nls"), (fitting == "nlrob") | (is.na(fixedWeight) == FALSE), str_detect(name, "RelHt") == FALSE) %>% 
    select(fitting, name, mae, mape, rmse, rmspe, aicN, nse, pearson) %>% 
    pivot_wider(names_from = fitting, values_from = c(mae, mape, rmse, rmspe, aicN, nse, pearson)) %>%
    mutate(deltaMae = mae_nlrob - mae_nls, deltaMape = mape_nlrob - mape_nls,
           deltaRmse = rmse_nlrob - rmse_nls, deltaRmspe = rmspe_nlrob - rmspe_nls,
           deltaAicN = aicN_nlrob - aicN_nls, deltaNse = nse_nlrob - nse_nls, deltaPearson = pearson_nlrob - pearson_nls)
  
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaMae), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, m", y = "number of regression forms") +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaRmse), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, m", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaAicN), bins = 30) +
    coord_cartesian(xlim = c(-3.75, -3.25), ylim = c(0, 7)) +
    labs(x = "normalized ΔAIC", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaNse), bins = 30) +
    coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(0, 7)) +
    labs(x = "change in model efficiency", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaPearson), bins = 30) +
    coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(0, 7)) +
    labs(x = "change in Pearson's R", y = NULL) +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaMape), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, %", y = "number of regression forms") +
  ggplot(psmeHeightFromDiameter$Efficiency) +
    geom_histogram(aes(x = deltaRmspe), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, %", y = NULL) +
  plot_layout(ncol = 5, nrow = 2)
  
  ggplot(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, error = -residuals(psmeHeightFromDiameter$sharmaParton))) +
    geom_boxplot(aes(x = dbhClass, y = error, group = dbhClass)) +
    coord_cartesian(ylim = c(-40, 40)) +
    labs(x = NULL, y = "fit_nlrob() height residual, m") +
  ggplot(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, error = -residuals(psmeHeightFromDiameter$sharmaPartonNls))) +
    geom_boxplot(aes(x = dbhClass, y = error, group = dbhClass)) +
    coord_cartesian(ylim = c(-40, 40)) +
    labs(x = "DBH, cm", y = "nls() height residual, m") +
  plot_layout(ncol = 1, nrow = 2)
  
  ggplot(bind_rows(tibble(dbhClass = 2.5 * floor(psme2016$DBH / 2.5) + 0.5 * 2.5, 
                          chapmanRichards = -residuals(psmeHeightFromDiameter$chapmanRichards),
                          chapmanRichardsBal = -residuals(psmeHeightFromDiameter$chapmanRichardsBal),
                          curtis = -residuals(psmeHeightFromDiameter$curtis),
                          gam = -residuals(psmeHeightFromDiameter$gam),
                          gamBal = -residuals(psmeHeightFromDiameter$gamBal),
                          hossfeld = -residuals(psmeHeightFromDiameter$hossfeld),
                          korf = -residuals(psmeHeightFromDiameter$korf),
                          linear = -residuals(psmeHeightFromDiameter$linear),
                          michaelisMenten = -residuals(psmeHeightFromDiameter$michaelisMenten),
                          parabolic = -residuals(psmeHeightFromDiameter$parabolic),
                          power = -residuals(psmeHeightFromDiameter$power),
                          prodan = -residuals(psmeHeightFromDiameter$prodan),
                          ratkowsky = -residuals(psmeHeightFromDiameter$ratkowsky),
                          richards = -residuals(psmeHeightFromDiameter$richards),
                          sharmaParton = -residuals(psmeHeightFromDiameter$sharmaParton), 
                          sharmaPartonBal = -residuals(psmeHeightFromDiameter$sharmaPartonBal),
                          sharmaZhang = -residuals(psmeHeightFromDiameter$sharmaZhang),
                          sharmaZhangBal = -residuals(psmeHeightFromDiameter$sharmaZhangBal),
                          sibbesen = -residuals(psmeHeightFromDiameter$sibbesen),
                          weibull = -residuals(psmeHeightFromDiameter$weibull),
                          weibullBal = -residuals(psmeHeightFromDiameter$weibullBal),
                          sharmaPartonNls = -residuals(psmeHeightFromDiameter$sharmaPartonNls)),
                   tibble(dbhClass = 2.5 * floor(psme2016physio$DBH / 2.5) + 0.5 * 2.5,
                          chapmanRichardsBalPhysio = -residuals(psmeHeightFromDiameter$chapmanRichardsBalPhysio),
                          chapmanRichardsPhysio = -residuals(psmeHeightFromDiameter$chapmanRichardsPhysio),
                          gamBalPhysio = -residuals(psmeHeightFromDiameter$gamBalPhysio),
                          gamPhysio = -residuals(psmeHeightFromDiameter$gamPhysio),
                          sharmaPartonBalPhysio = -residuals(psmeHeightFromDiameter$sharmaPartonBalPhysio),
                          sharmaPartonPhysio = -residuals(psmeHeightFromDiameter$sharmaPartonPhysio))) %>%
           pivot_longer(cols = !dbhClass, names_to = "name", values_to = "residual") %>%
           mutate(name = factor(name, levels = c("chapmanRichards", "chapmanRichardsBal", "chapmanRichardsBalPhysio", "chapmanRichardsPhysio", "curtis", "gam", "gamBal", "gamBalPhysio", "gamPhysio", "hossfeld", "korf", "linear", "michaelisMenten", "parabolic", "power", "prodan", "ratkowsky", "richards", "sharmaParton", "sharmaPartonBal", "sharmaPartonBalPhysio", "sharmaPartonPhysio", "sharmaZhang", "sharmaZhangBal", "sibbesen", "weibull", "weibullBal", "sharmaPartonNls"),
                                      labels = c("Chapman-Richards", "Chapman-Richards BA+L", "Chapman-Richards BA+L physio", "Chapman-Richards physio", "Curtis", "GAM", "GAM BAL", "GAM BAL physio", "GAM physio", "Hossfeld", "Korf", "linear", "Michaelis-Menten", "parabolic", "power", "Prodan", "Ratkowsky", "Richards", "Sharma-Parton", "Sharma-Parton BA+L", "Sharma-Parton BA+L physio", "Sharma-Parton physio", "Sharma-Zhang", "Sharma-Zhang BA+L", "Sibbesen", "Weibull", "Weibull BA+L", "Sharma-Parton nls()"))) %>%
           group_by(name, dbhClass) %>%
           summarize(n = n(), bias = mean(residual, na.rm = TRUE), 
                     fitting = if_else(name %in% c("GAM", "GAM BAL", "GAM BAL physio", "GAM"), "GAM", if_else(name %in% c("linear", "parabolic"), "linear", "nonlinear")),
                     .groups = "drop") %>%
           filter(n > 10, is.na(bias) == FALSE)) +
    geom_line(aes(x = dbhClass, y = bias, alpha = fitting, color = fitting, group = name)) +
    coord_cartesian(ylim = c(-5, 5)) +
    labs(x = "DBH, cm", y = "height bias, m", alpha = NULL, color = NULL) +
    scale_alpha_manual(breaks = c("GAM", "nonlinear", "linear"), values = c(1, 0.2, 0.3)) +
    scale_color_discrete(breaks = c("GAM", "nonlinear", "linear")) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))
}


## Douglas-fir height-diameter GNLS regressions
if (htDiaOptions$fitSlowGnls)
{
  if (exists("psmeHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/PSME TotalHt primary.rdata") }

  psmeHeightFromDiameterGnls = list(chapmanRichards = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))) # step halving at nlsTol = 0.05
  psmeHeightFromDiameterGnls$chapmanRichardsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 0.001, tolerance = 0.01, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05, iterations at default msTol and tolerance
  psmeHeightFromDiameterGnls$sharmaParton = gnls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 500, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
  psmeHeightFromDiameterGnls$sharmaPartonBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
  psmeHeightFromDiameterGnls$sharmaZhang = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, msTol = 1E-5, tolerance = 1E-4, returnObject = FALSE)) # step having at nlsTol = 0.1
  psmeHeightFromDiameterGnls$sharmaZhangBal = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, nlsMaxIter = 50, msTol = 1E-7, tolerance = 1E-6, msVerbose = FALSE, returnObject = FALSE)) # step halving factor at nlsTol = 1 with plot correlation, step halving with nlsTol = 0.07
  psmeHeightFromDiameterGnls$weibull = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
  psmeHeightFromDiameterGnls$weibullBal = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 5, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.2
  
  psmeHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", psmeHeightFromDiameterGnls$chapmanRichards, psme2016)
  psmeHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", psmeHeightFromDiameterGnls$chapmanRichardsBal, psme2016)
  psmeHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", psmeHeightFromDiameterGnls$sharmaParton, psme2016)
  psmeHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", psmeHeightFromDiameterGnls$sharmaPartonBal, psme2016)
  psmeHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", psmeHeightFromDiameterGnls$sharmaZhang, psme2016)
  psmeHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", psmeHeightFromDiameterGnls$sharmaZhangBal, psme2016)
  psmeHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", psmeHeightFromDiameterGnls$weibull, psme2016)
  psmeHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", psmeHeightFromDiameterGnls$weibullBal, psme2016)

  save(file = "trees/height-diameter/data/PSME TotalHt GNLS.rdata", psmeHeightFromDiameterGnls)
}

if (htDiaOptions$includeInvestigatory)
{
  psmeHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(psmeHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(psmeHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(psmeHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  
  ggplot() +
    geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsBal), color = "Chapman-Richards BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsPhysio), color = "Chapman-Richards physio", group = psme2016physio$isPlantation), alpha = 0.5) +
    geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameter$chapmanRichardsTopo), color = "Chapman-Richards topo", group = psme2016physio$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaPartonBal), color = "Sharma-Parton BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$sharmaZhangBal), color = "Sharma-Zhang BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibullBal), color = "Weibull BA+L", group = psme2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterGnls$sharmaParton), color = "Sharma-Parton GNLS"), alpha = 0.35) +
    #geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterGnls$weibull), color = "Weibull GNLS"), alpha = 0.35) + 
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = psme2016$isPlantation)) +
    geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameter$weibull), color = "Weibull", group = psme2016$isPlantation)) +
    annotate("text", x = 0, y = 85, label = "a) Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c("ElliottChapmanRichards", "ElliottWeibull", "ElliottChapmanRichardsBal", "ElliottWeibullBalNatural", "ElliottWeibullBalPlantation", "TemesgenWeibull"), labels = c("Chapman-Richards", "Weibull", "Chapman-Richards with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_color_manual(values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot(psme2016) +
    geom_point(aes(x = Elev_Mean, y = -residuals(psmeHeightFromDiameter$weibullBal)), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = Elev_Mean, y = -residuals(psmeHeightFromDiameter$weibullBal)), color = "red", formula = y ~ s(x, k = 10), method = "gam")
  
  ggplot(psme2016) +
    geom_point(aes(x = AspectCos, y = -residuals(psmeHeightFromDiameter$weibullBal)), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = AspectCos, y = -residuals(psmeHeightFromDiameter$weibullBal)), color = "red", formula = y ~ s(x, k = 10), method = "gam")
}


## Douglas-fir diameter-height regressions
if (psmeOptions$fitDbhPrimary)
{
  psmeDiameterFromHeight = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780)))
  psmeDiameterFromHeight$chapmanFormAat = fit_nlrob("Chapman-Richards form ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.4, a1p = -44.6, a2 = 0.0058, a3 = -0.0544, b1 = 0.00166, b1p = 0.017, b2 = 0.788, b2p = -0.056)) # a2, a2p, a3p not significant
  psmeDiameterFromHeight$chapmanFormBal = fit_nlrob("Chapman-Richards form BA+L", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 135, a1p = -37.5, a2 = -1.2, b1 = 0.010, b1p = 0.002, b2 = 0.756, b2p = 0.064), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a3 * BA
  psmeDiameterFromHeight$chapmanFormBalRelHt = fit_nlrob("Chapman-Richards form BA+L RelHt", DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 139, a1p = -38.5, a2 = -1.2, a9 = -5.5, a9p = 0.22, b1 = 0.012, b1p = 0.003, b2 = 0.796, b2p = 0.066), control = nls.control(maxiter = 500)) # a2p not significant, fit_nlrob() step factor with a9 * BA
  psmeDiameterFromHeight$chapmanFormRelHt = fit_nlrob("Chapman-Richards form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 21.1, a1p = -6.56, a9 = 0.278, b1 = 0.170, b2 = 0.580, b2p = 0.0395)) # a4p not significant
  psmeDiameterFromHeight$chapmanRichards = fit_nlrob("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781)) # b2p not significant
  psmeDiameterFromHeight$chapmanRichardsAat = fit_nlrob("Chapman-Richards ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786)) # a2p, a3p, b2p not significant
  psmeDiameterFromHeight$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -13.8, a1p = -0.65, a5 = -3.58, a8 = 0.091, b1 = 0.019, b1p = 0.0062, b2 = 0.30)) # a4, a6, a7 not significant
  psmeDiameterFromHeight$chapmanRichardsRelHt = fit_nlrob("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -1.67, a9 = -9.26, b1 = 0.020, b1p = 0.004, b2 = 0.004, b2p = 0.089)) # a1p not significant
  psmeDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 10, pc = psme2016gamConstraint), data = psme2016)
  psmeDiameterFromHeight$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28, pc = psme2016gamConstraint), data = psme2016, nthreads = 4)

  if (psmeOptions$fitSlowGams)
  {
    psmeDiameterFromHeight$gamAatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 331, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6) # slow since minimum (k = 496, edf = 397, 2x2 cross validation = 18 minutes with Zen 3, 4.6 GHz), AIC 151550: 151552 without AAT, 151700 without ABA, 151783 without elevation, 151533 without slope, 151561 without sin(aspect), 151562 without cos(aspect), 151650 without topographic shelter -> drop slope -> AIC 151800 without AAT, 151920 without ABA, 151926 without elevation, 151742 without sin(aspect), 151693 without cos(aspect), 151902 without topographic shelter
    psmeDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6)
    
    psmeDiameterFromHeightGamAatPhysio = psmeDiameterFromHeight$gamBalPhysio
    psmeDiameterFromHeightGamPhysio = psmeDiameterFromHeight$gamPhysio
    save(file = "trees/height-diameter/data/PSME DBH primary GAMs.Rdata", psmeDiameterFromHeightGamAatPhysio, psmeDiameterFromHeightGamPhysio)
    rm(psmeDiameterFromHeightGamAatPhysio, psmeDiameterFromHeightGamPhysio)
  } else {
    load("trees/height-diameter/data/PSME DBH primary GAMs.Rdata")
    psmeDiameterFromHeight$gamBalPhysio = psmeDiameterFromHeightGamAatPhysio
    psmeDiameterFromHeight$gamPhysio = psmeDiameterFromHeightGamPhysio
    rm(psmeDiameterFromHeightGamAatPhysio, psmeDiameterFromHeightGamPhysio)
  }
  
  psmeDiameterFromHeight$linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016)
  psmeDiameterFromHeight$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08))
  psmeDiameterFromHeight$naslund = fit_nlrob("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018))
  psmeDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016)
  psmeDiameterFromHeight$power = fit_nlrob("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108))
  psmeDiameterFromHeight$powerAat = fit_nlrob("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053)) # a3 not significant
  psmeDiameterFromHeight$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, b1 = 1.03, b1p = -0.102)) # a7, a8 not significant
  psmeDiameterFromHeight$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068)) # a1p, a4p not significant
  psmeDiameterFromHeight$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096))
  #psmeDiameterFromHeight$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 18.6)) # converges from red alder values but fails to reconverge (singular gradient or NaN-inf with nls())
  psmeDiameterFromHeight$sharmaParton = fit_nlrob("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 3.95, b1 = 0.681, b1p = -0.139, b2 = 0.097, b3 = -0.130, b3p = 0.157, b4 = 0.125, b4p = 0.0678), control = list(maxiter = 200)) # a1p NaN-inf (not significant?), b4p significance debatable, singular gradient with all relative height forms attempted
  psmeDiameterFromHeight$sibbesenForm = fit_nlrob("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017))
  psmeDiameterFromHeight$sibbesenFormAat = fit_nlrob("Sibbesen form ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, a3p = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190)) # a2, a2p not significant
  psmeDiameterFromHeight$sibbesenFormPhysio = fit_nlrob("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = list(a1 = 3.812, a1p = -0.925, a4 = 0.00022, a6 = -0.0495, a7 = -0.0151, a8 = -0.00630, b1 = 0.520, b2 = 0.111, b2p = 0.0173)) # a5 not significant
  psmeDiameterFromHeight$sibbesenFormRelHt = fit_nlrob("Sibbesen form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.90, a1p = -0.95, a9 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016)) # a4p not significant
  psmeDiameterFromHeight$weibull = fit_nlrob("Weibull", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81)) # b2p not significant
  #confint_nlrob(psmeDiameterFromHeight$sibbesenFormRelHt, level = 0.99)
  
  save(file = "trees/height-diameter/data/PSME DBH primary.rdata", psmeDiameterFromHeight)
}

if (psmeOptions$fitDbhGslNls)
{
  if (exists("psmeDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/PSME DBH primary.rdata") }
  
  psmeDiameterFromHeightNls = list(chapmanForm = fit_gsl_nls("Chapman-Richards form", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanForm$fit[[1]]$m$getPars()))
  psmeDiameterFromHeightNls$chapmanFormAat = fit_gsl_nls("Chapman-Richards form ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanFormAat$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$chapmanFormRelHt = fit_gsl_nls("Chapman-Richards form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = psmeDiameterFromHeight$chapmanFormRelHt$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$chapmanRichards = fit_gsl_nls("Chapman-Richards", DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichards$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$chapmanRichardsAat = fit_gsl_nls("Chapman-Richards ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsAat$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", DBH ~ (a1 + a1p * isPlantation + a5 * sin(3.14159/180 * slope) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = psmeDiameterFromHeight$chapmanRichardsPhysio$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsRelHt$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$michaelisMentenForm = fit_gsl_nls("Michaelis-Menten form", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = psmeDiameterFromHeight$michaelisMentenForm$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$naslund = fit_gsl_nls("Näslund", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$naslund$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$power = fit_gsl_nls("power", DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$power$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$powerAat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$powerAat$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect))*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = psmeDiameterFromHeight$powerPhysio$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = psmeDiameterFromHeight$powerRelHt$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$ruark$fit[[1]]$m$getPars())
  #psmeDiameterFromHeightNls$schnute = fit_gsl_nls("Schnute", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = psmeDiameterFromHeight$schnute$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 5, b1 = 0.6, b1p = -0.13, b2 = 0.09, b3 = -0.07, b3p = 0.12, b4 = 0.23, b4p = 0.09), control = list(maxiter = 250))
  psmeDiameterFromHeightNls$sibbesenForm = fit_gsl_nls("Sibbesen form", DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenForm$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$sibbesenFormAat = fit_gsl_nls("Sibbesen form ABA+T", DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormAat$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$sibbesenFormPhysio = fit_gsl_nls("Sibbesen form physio", DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = psmeDiameterFromHeight$sibbesenFormPhysio$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$sibbesenFormRelHt = fit_gsl_nls("Sibbesen form RelHt", DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormRelHt$fit[[1]]$m$getPars())
  psmeDiameterFromHeightNls$weibull = fit_gsl_nls("Weibull", DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = psmeDiameterFromHeight$weibull$fit[[1]]$m$getPars())
  
  save(file = "trees/height-diameter/data/PSME DBH gsl_nls.rdata", psmeDiameterFromHeightNls)
}

if (psmeOptions$fitHeightPrimary & psmeOptions$fitHeightGslNls & psmeOptions$fitDbhPrimary & psmeOptions$fitDbhGslNls)
{
  # assemble parameter and results tibbles using incremental load and reduce to fit in memory
  # height primary + GNLS + gsl_nls = 112 GB DDR at 10x10 cross validation
  if (exists("psmeHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/PSME TotalHt primary.rdata") }
  psmeCoefficients = append_model_coefficients(tibble(), psmeHeightFromDiameter, "height")
  psmeResults = append_model_results(tibble(), psmeHeightFromDiameter, "height")
  rm(psmeHeightFromDiameter)
  gc() # force actual release of memory for variables just removed as otherwise R will go over 128 GB (this is slower than loading)
                             
  if (exists("psmeHeightFromDiameterNls") == FALSE) { load("trees/height-diameter/data/PSME TotalHt gsl_nls.rdata") }
  psmeCoefficients = append_model_coefficients(psmeCoefficients, psmeHeightFromDiameterNls, "height", fitSet = "gsl_nls")
  psmeCoefficients = append_model_coefficients(psmeCoefficients, psmeHeightFromDiameterNlsDwt, "height", fitSet = "gsl_nls", fixedWeight = -1)
  psmeResults = append_model_results(psmeResults, psmeHeightFromDiameterNls, "height", fitSet = "gsl_nls")
  psmeResults = append_model_results(psmeResults, psmeHeightFromDiameterNlsDwt, "height", fitSet = "gsl_nls", fixedWeight = -1)
  rm(psmeHeightFromDiameterNls, psmeHeightFromDiameterNlsDwt)
  gc()
    
  if (exists("psmeHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/PSME TotalHt GNLS.rdata") }
  psmeCoefficients = bind_rows(psmeCoefficients,
                             bind_rows(lapply(psmeHeightFromDiameterGnls, get_model_coefficients)) %>%
                              mutate(responseVariable = "height"))
  psmeResults = bind_rows(psmeResults,
                          bind_rows(lapply(psmeHeightFromDiameterGnls, get_model_results)) %>%
                            mutate(responseVariable = "height"))
  rm(psmeHeightFromDiameterGnls) # small, so no need to gc()

  # DBH primary + gsl_nls = 86 GB DDR at 10x10 cross validation
  if (exists("psmeDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/PSME DBH primary.rdata") }
  psmeCoefficients = append_model_coefficients(psmeCoefficients, psmeDiameterFromHeight, "DBH")
  psmeResults = append_model_results(psmeResults, psmeDiameterFromHeight, "DBH")
  rm(psmeDiameterFromHeight)
  gc()
  
  if (exists("psmeDiameterFromHeightNls") == FALSE) { load("trees/height-diameter/data/PSME DBH gsl_nls.rdata") }
  psmeCoefficients = append_model_coefficients(psmeCoefficients, psmeDiameterFromHeightNls, "DBH", fitSet = "gsl_nls") %>%
    mutate(species = "PSME") %>%
    relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)
  psmeResults = append_model_results(psmeResults, psmeDiameterFromHeightNls, "DBH", fitSet = "gsl_nls") %>%
    mutate(species = "PSME") %>%
    relocate(responseVariable, species)
    
  save(file = "trees/height-diameter/data/PSME results.Rdata", psmeCoefficients, psmeResults)
  rm(psmeDiameterFromHeightNls)
  gc()
}

if (htDiaOptions$includeInvestigatory)
{
  print(psmeDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman form", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman form ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanFormBal), y = TotalHt, color = "Chapman form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$linear), y = TotalHt, color = "linear", group = isPlantation)) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$michaelisMentenForm), y = TotalHt, color = "Michaelis-Menten form", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$parabolic), y = TotalHt, color = "parabolic", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$schnute), y = TotalHt, color = "Schnute", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sibbesenForm), y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
    #geom_line(aes(x = predict(psmeDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    geom_line(aes(x = predict(psmeDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -65 * log(1 - pmin((1/85*(TotalHt - 1.37))^0.7, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = -155 * log(1 - pmin((0.00839*(TotalHt - 1.37))^0.970, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
    #geom_line(aes(x = 0.05*topHeight*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman form BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 3.5*sqrt(TotalHt - 1.37) / (1 - 0.1*sqrt(TotalHt - 1.37)), y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = (-204*log(1 - pmin(0.011 * (TotalHt - 1.37), 0.9999)))^0.869, y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
    annotate("text", x = 0, y = 90, label = "Douglas-fir, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
  
  ggplot() +
    #geom_segment(aes(x = 0, y = 0, xend = 80, yend = 50), color = "grey70", linetype = "longdash") +
    #geom_segment(aes(x = 0, y = 0, xend = 80, yend = -50), color = "grey70", linetype = "longdash") +
    geom_point(aes(x = psme2016$tallerTph, y = -residuals(psmeDiameterFromHeight$chapmanRH)), alpha = 0.05, shape = 16) +
    geom_smooth(aes(x = psme2016$tallerTph, y = -residuals(psmeDiameterFromHeight$chapmanRH)), color = "red", fill = "red") +
    labs(x = "height, m", y = "residual DBH error, cm")
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.15, na.rm = TRUE, shape = 16) +
    geom_line(aes(x = DBH, y = 1.37 + 300 * (1 * 0.1 * pmin(relativeHeight, 2)) * exp(-0.5*DBH^-1)), color = "red")
  
  ggplot(psme2016) +
    geom_point(aes(x = DBH, y = heightDiameterRatio, color = Species), alpha = 0.2, na.rm = TRUE, shape = 16) +
    labs(x = "DBH, cm", y = "height-diameter ratio", color = NULL) +
    scale_color_manual(breaks = c("DF", "WH"), values = c("green3", "blue2")) +
    scale_y_continuous(breaks = seq(0, 300, by = 40)) +
    theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.98))
  
  psmeDiameterFromHeight$Efficiency = psmeDiameterFromHeightResults %>% filter(fitting %in% c("nlrob", "nls"), str_detect(name, "BA\\+L") == FALSE) %>% 
    select(fitting, name, mae, mape, rmse, rmspe, aicN, nse, pearson) %>% 
    pivot_wider(names_from = fitting, values_from = c(mae, mape, rmse, rmspe, aicN, nse, pearson)) %>%
    mutate(deltaMae = mae_nlrob - mae_nls, deltaRmse = rmse_nlrob - rmse_nls, deltaAicN = aicN_nlrob - aicN_nls, deltaNse = nse_nlrob - nse_nls, deltaPearson = pearson_nlrob - pearson_nls)
  
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaMae), bins = 30) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 8)) +
    labs(x = "ΔMAE, m", y = "number of regression forms") +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaRmse), bins = 30) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 8)) +
    labs(x = "ΔRMSE, cm", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaAicN), bins = 30) +
    coord_cartesian(xlim = c(-6.5, -6), ylim = c(0, 8)) +
    labs(x = "normalized ΔAIC", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaNse), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 8)) +
    labs(x = "change in model efficiency", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaPearson), bins = 30) +
    coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0, 8)) +
    labs(x = "change in Pearson's R", y = NULL) +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaMape), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔMAE, %", y = "number of regression forms") +
  ggplot(psmeDiameterFromHeight$Efficiency) +
    geom_histogram(aes(x = deltaRmspe), bins = 30) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(0, 7)) +
    labs(x = "ΔRMSE, %", y = NULL) +
  plot_layout(ncol = 5, nrow = 2)
}


## preferred forms identified (results.R, Figure 5)
psmeHeightFromDiameterPreferred = list(gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = psme2016gamConstraint), data = psme2016, folds = 1, repetitions = 1))
psmeHeightFromDiameterPreferred$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), folds = 1, repetitions = 1)
psmeHeightFromDiameterPreferred$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22), folds = 1, repetitions = 1)
psmeHeightFromDiameterPreferred$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 20, a1p = -9, b1 = 0.3, b1p = 0.07, b2 = -0.023, b2p = -0.02, b3 = 0.0, b4 = 1.53, b4p = -0.15), folds = 1, repetitions = 1)
psmeHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51), folds = 1, repetitions = 1)
psmeHeightFromDiameterPreferred$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.0006, a1p = 0.17, b1 = 5.8, b1p = -3.5, b2 = -0.182, b2p = 0.050), folds = 1, repetitions = 1)

psmeDiameterFromHeightPreferred = list(chapmanForm = fit_nlrob("Chapman-Richards form", DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), folds = 1, repetitions = 1))
psmeDiameterFromHeightPreferred$gamAat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 28, pc = psme2016gamConstraint), data = psme2016, nthreads = 4, folds = 1, repetitions = 1)
psmeDiameterFromHeightPreferred$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = psme2016gamConstraint), data = psme2016physio, nthreads = 6, folds = 1, repetitions = 1)
psmeDiameterFromHeightPreferred$michaelisMentenForm = fit_nlrob("Michaelis-Menten form", DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), folds = 1, repetitions = 1)
psmeDiameterFromHeightPreferred$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), folds = 1, repetitions = 1)

save(file = "trees/height-diameter/data/PSME preferred models.Rdata", psmeHeightFromDiameterPreferred, psmeDiameterFromHeightPreferred)


## Q-Q plots
if (htDiaOptions$includeInvestigatory)
{
  ggplot() + # symmetric t fits less well than skewed t, as expected
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'t(df = 7, '*alpha*' = 2.25)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
  ggplot() + # slow! and fits poorly
  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
    annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'pearsonIV(m = 1.72, '*nu*' = 0.656)'", hjust = 0, parse = TRUE, size = 3.5) +
    coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
    labs(x = "theoretical quantile", y = NULL, color = NULL) +
    scale_color_manual(values = dbhColors) +
    theme(legend.justification = c(1, 0), legend.position = "none")
}


## basal area from height
# essentially no difference between fit_gsl_nls() and fit_nlrob() fits
# Chapman-Richards has the wrong curvature
if (htDiaOptions$includeInvestigatory)
{
  #psmeBasalAreaFromHeightKorf = fit_gsl_nls(basalArea ~ a1*exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1, psme2016, start = list(a1 = 1, a1p = 0, b1 = 0.00009, b1p = 0, b2 = 2.2, b2p = 0), weights = pmin(1/basalArea, 1E4)) # a1p not significant
  #psmeBasalAreaFromHeightPower = fit_gsl_nls(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 0.25 * pi * 0.01^2, a1p = 0, b1 = 2.4, b1p = 0), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
  psmeBasalAreaFromHeightKorf = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 0.689, a1p = -0.413, b1 = 0.0003, b1p = 0.0005, b2 = 1.91, b2p = -0.10), weights = pmin(1/basalArea, 1E4)) # a1p not significant
  psmeBasalAreaFromHeightPower = fit_nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = 0.00005, b1 = 2.41, b1p = -0.248), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
  #confint_nlrob(psmeBasalAreaFromHeightKorf, level = 0.99, weights = 1/psme2016$basalArea)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(psmeBasalAreaFromHeightKorf), 100^2 * mean(-residuals(psmeBasalAreaFromHeightKorf)), mean(abs(residuals(psmeBasalAreaFromHeightKorf))), 1 - sum(residuals(psmeBasalAreaFromHeightKorf)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2)),
          "power", AIC(psmeBasalAreaFromHeightPower), 100^2 * mean(-residuals(psmeBasalAreaFromHeightPower)), mean(abs(residuals(psmeBasalAreaFromHeightPower))), 1 - sum(residuals(psmeBasalAreaFromHeightPower)^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(psme2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(psmeBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "measured or imputed Douglas-fir height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
  
  ggplot(psme2016) +
    geom_point(aes(x = basalArea, y = -residuals(psmeBasalAreaFromHeightKorf), color = "Korf"), alpha = 0.1, shape = 16) +
    labs(x = "basal area, m²", y = "residual, m", color = NULL) +
    theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
  ggplot(psme2016) +
    geom_point(aes(x = basalArea, y = -residuals(psmeBasalAreaFromHeightPower), color = "power"), alpha = 0.1, shape = 16) +
    labs(x = "basal area, m²", y = "residual, m", color = NULL) +
    theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
  plot_layout(nrow = 1, ncol = 2)
    
  ggplot(psme2016) +
    geom_point(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.1, color = "grey25", shape = 16) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.10, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.95)) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.25, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.80)) +
    stat_summary_bin(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), alpha = 0.35, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.50)) +
    geom_quantile(aes(x = imputedHeight, y = -residuals(psmeBasalAreaFromHeightKorf)), color = "darkviolet", formula = y ~ x) +
    labs(x = "Douglas-fir height, m", y = "basal area residual, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## exploratory plots
if (htDiaOptions$includeInvestigatory)
{
  ggplot() +
    geom_point(aes(x = psme2016natural$DBH, y = psme2016natural$TotalHt), alpha = 0.10, color = "green4", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = psme2016natural$DBH, y = psme2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "Douglas-fir naturally regenerated height, m") +
  ggplot() +
    geom_point(aes(x = psme2016plantation$DBH, y = psme2016plantation$TotalHt), alpha = 0.10, color = "grey20", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = psme2016plantation$DBH, y = psme2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "Douglas-fir plantation height, m")
  #ggsave("Presentation/Douglas-fir height-diameter natural-plantation.png", width = 12.5, height = 11, units = "cm")
}