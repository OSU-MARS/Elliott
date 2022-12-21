# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Douglas-fir height-diameter regression form sweep
# preferred forms: Sharma-Parton BA+L, Sharma-Parton, Sharma-Zhang BA+L, Chapman-Richards BA+L
# fits without isPlantation as a factor
#psmeHeightFromDiameter$chapmanRichards = gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 74, b1 = -0.014, b2 = 1.14), weights = dbhWeight)
#psmeHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = 0.49, a3 = -0.053, b1 = -0.011, b2 = 1.09), weights = dbhWeight)
#psmeHeightFromDiameter$chapmanRichardsBalRelHt = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a9 * relativeHeight) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 35, a2 = 0.70, a3 = 0.055, a9 = 38.32, b1 = -0.010, b2 = 0.92), weights = dbhWeight)
#psmeHeightFromDiameter$chapmanRichardsPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -9.56, a2 = 0.0003, a2p = -0.011, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = dbhWeight) # b1p not significant
#psmeHeightFromDiameter$chapmanRichardsPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation + a3 * sin(pi/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 66.7, a1p = -9.66, a2 = 0.0003, a2p = -0.011, a3 = -2.06, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = dbhWeight) # a3p not significant, little sensitivity to sin() or cos() of slope but sin() slightly more accurate
#psmeHeightFromDiameter$chapmanRichardsTopo = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * topographicShelterIndex) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, a2 = 0, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = dbhWeight) # accuracy increases through cos(), sin(), linear but differences are limited
#psmeHeightFromDiameter$curtis = gsl_nls(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b2, psme2016, start = list(a1 = 1.6, b2 = 0.24), weights = dbhWeight)
#psmeHeightFromDiameter$hossfeld = gsl_nls(TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), psme2016, start = list(a1 = 99, b1 = 187, b2 = -1.2), weights = dbhWeight)
#psmeHeightFromDiameter$korf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, a1p = 376, b1 = -7.83, b2 = -0.323, b2p = 0.084), weights = dbhWeight, control = nls.control(maxiter = 500)) # b1p not significant
#psmeHeightFromDiameter$michaelisMenten = gsl_nls(TotalHt ~ 1.37 + a1*DBH / (a2 + DBH), psme2016, start = list(a1 = 138, a2 = 156), weights = dbhWeight)
#psmeHeightFromDiameter$prodan = gsl_nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), psme2016, start = list(a1 = 0.008, a2 = 1.0, a3 = 2.7), weights = dbhWeight)
#psmeHeightFromDiameter$power = gsl_nls(TotalHt ~ 1.37 + b0*DBH^b1, psme2016, start = list(b0 = 1.54, b1 = 0.77), weights = dbhWeight)
#psmeHeightFromDiameter$ratkowsky = gsl_nls(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), psme2016, start = list(a1 = 93, b1 = -65, b2 = 15), weights = dbhWeight)
#psmeHeightFromDiameter$richards = gsl_nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 58.3, d = 0.609, kU = 0.0134), weights = dbhWeight)
#psmeHeightFromDiameter$richards = gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 62.6, Hap = -26.5, d = 0.609, kU = 0.0128, kUp = 0.0114), weights = dbhWeight)
#psmeHeightFromDiameter$sharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 55.4, a2 = 0.082, b1 = -0.022, b2 = -0.17, b3 = 1.04), weights = dbhWeight)
#psmeHeightFromDiameter$sharmaPartonBal = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, psme2016, start = list(a1 = 39, a2 = 0.15, b1 = -0.018, b2 = -0.16, b3 = 1.01), weights = dbhWeight)
#psmeHeightFromDiameter$sharmaZhang = gsl_nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 55, a2 = 0.07, b1 = -0.012, b2 = 0.03, b3 = 1.1), weights = dbhWeight)
#psmeHeightFromDiameter$sharmaZhangBal = gsl_nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 66, a2 = 0.056, a3 = 0.006, b1 = -0.022, b2 = -0.13, b3 = 1.04), weights = dbhWeight)
#psmeHeightFromDiameter$sibbesen = gsl_nls(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), psme2016, start = list(a1 = 0.034, b1 = 3.1, b2 = -0.15), weights = dbhWeight)
#psmeHeightFromDiameter$korf = gsl_nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), psme2016, start = list(a1 = 172, b1 = -9.682, b2 = -0.4574), weights = dbhWeight)
#psmeHeightFromDiameter$weibull = gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 69, b1 = -0.0075, b2 = 1.15), weights = dbhWeight)
#psmeHeightFromDiameter$weibullBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 73.7, a2 = 0.38, a3 = -0.007, b1 = -0.008, b2 = 1.09), weights = dbhWeight)
#psmeHeightFromDiameter$weibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation + b2 * pmin(relativeHeight, 1.25))*DBH^(b3 + b3p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a3p = 0.217, a4 = 41.9, a4p = 76.9, b1 = -0.013, b1p = 0.007, b2 = 0, b3 = 0.94, b3p = -0.096), weights = dbhWeight, control = nls.control(maxiter = 500)) # step factor with nlrob() 
# fits with isPlantation
#psmeHeightFromDiameter$chapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 68.6, a1p = -3.2, a2 = 0.024, a2p = 0.96, a3 = 0.014, a3p = -0.15, b1 = -0.018, b1p = 0.0028, b2 = 1.29, b2p = -0.10), weights = dbhWeight)
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, TotalHt > 0) %>% # live Douglas-firs measured for height
  mutate(dbhWeight = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1),
         heightWeight = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016gamConstraint = c(DBH = -1.2240/0.6566, TotalHt = 1.37, standBasalAreaPerHectare = median(psme2016$standBasalAreaPerHectare), basalAreaLarger = median(psme2016$basalAreaLarger), standBasalAreaApprox = median(psme2016$standBasalAreaApprox), tallerApproxBasalArea = median(psme2016$tallerApproxBasalArea), elevation = median(psme2016physio$elevation), slope = median(psme2016physio$slope), aspect = median(psme2016physio$aspect), topographicShelterIndex = median(psme2016physio$topographicShelterIndex), relativeHeight = median(psme2016$relativeHeight)) # point constraint for mgcv::s() where response variable is ignored, zero crossing of height from DBH from lm(TotalHt ~ DBH, data = psme2016 %>% filter(DBH < 6))
#psme2016natural = psme2016 %>% filter(isPlantation == FALSE)
#psme2016plantation = psme2016 %>% filter(isPlantation)
#psme2016plantationPhysio = psme2016physio %>% filter(isPlantation)

psmeHeightFromDiameter = list(chapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = dbhWeight)) # b1p not significant
psmeHeightFromDiameter$chapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 72.9, a1p = -11.8, a2 = 0.087, a2p = 0.84, a3 = -0.0021, a3p = -0.073, b1 = -0.016, b2 = 1.26, b2p = -0.054), weights = dbhWeight) # a3 not significant, step factor with b1p
psmeHeightFromDiameter$chapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 74.3, a2 = 0.096, a2p = 0.92, a3 = 0, a3p = 0, a4 = -0.015, a5 = -0.101, a6 = 0.793, a7 = 1.695, a8 = 0.183, b1 = -0.018, b1p = 0.005, b2 = 1.30, b2p = -0.154), weights = dbhWeight) # a4 not significant
psmeHeightFromDiameter$chapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a4 = 54.0, a4p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37), weights = dbhWeight)
psmeHeightFromDiameter$chapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -13.4, a4 = -0.0045, a5 = -8.09, a6 = 0.783, a7 = 0.766, a8 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = dbhWeight) # a4p not significant, a5p induces overfitting
psmeHeightFromDiameter$curtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156), weights = dbhWeight)
psmeHeightFromDiameter$gam = gam(TotalHt ~ s(DBH, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016, family = scat, select = TRUE, weights = dbhWeight)
psmeHeightFromDiameter$gamBal = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016, family = scat, select = TRUE, weights = dbhWeight)
#psmeHeightFromDiameter$gamBalPhysio = gam(TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016physio, select = TRUE, , weights = dbhWeight)
psmeHeightFromDiameter$gamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016physio, select = TRUE, weights = dbhWeight)
psmeHeightFromDiameter$hossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28), weights = dbhWeight)
psmeHeightFromDiameter$korf = nlrob(TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), weights = dbhWeight, control = nls.control(maxiter = 500)) # a1p parameter evaporation, b1p not significant
psmeHeightFromDiameter$linear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016, offset = breastHeight, weights = dbhWeight)
psmeHeightFromDiameter$michaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30), weights = dbhWeight) # b1p not significant
psmeHeightFromDiameter$parabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), psme2016, offset = breastHeight, weights = dbhWeight)
psmeHeightFromDiameter$prodan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), weights = dbhWeight)
psmeHeightFromDiameter$power = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14), weights = dbhWeight)
psmeHeightFromDiameter$ratkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52), weights = dbhWeight)
psmeHeightFromDiameter$richards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126), weights = dbhWeight)
psmeHeightFromDiameter$sharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 37.66, b1 = 0.19, b1p = -0.123, b2 = -0.017, b2p = -0.026, b3 = 0.061, b3p = -0.259, b4 = 1.33, b4p = -0.22), weights = dbhWeight)
psmeHeightFromDiameter$sharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 18.5, a1p = 11.3, b1 = 0.30, b1p = -0.14, b2 = -0.019, b2p = -0.011, b3 = 0.089, b4 = 1.49, b4p = -0.44), weights = dbhWeight) # b3p not significant
psmeHeightFromDiameter$sharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 52.6, a1p = -0.10, a4 = 0.00004, a5 = 0, a6 = 0.0090, a7 = 0.0032, a8 = 0.0040, b1 = 0.53, b2 = -0.025, b2p = -0.0090, b3 = 0.036, b3p = -0.19, b4 = 1.57, b4p = -0.51), weights = dbhWeight) # a4, a7, b1p not significant
psmeHeightFromDiameter$sharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = list(a1 = 51.9, a4 = 0.0004, a5 = 0, a6 = 0.010, a7 = 0.003, a8 = 0.004, b1 = 0.056, b2 = -0.024, b3 = -0.012, b3p = -0.23, b4 = 1.56, b4p = -0.64), weights = dbhWeight) # a7, b1p, b2p not significant, nlrob() step factor with a1p
psmeHeightFromDiameter$sharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 56.1, a1p = -23.1, b1 = 0.042, b1p = 0.117, b2 = -0.0247, b2p = -0.0131, b3 = -0.0217, b3p = -0.112, b4 = 1.476, b4p = -0.456), weights = dbhWeight) # b1 not significant
psmeHeightFromDiameter$sharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 53.5, a1p = -10.1, a2 = -0.00003, a2p = 0.012, b1 = 0.070, b2 = -0.027, b2p = -0.017, b3 = -0.064, b3p = -0.096, b4 = 1.28, b4p = -0.16), weights = dbhWeight) # a2, b1p, b3 not significant
psmeHeightFromDiameter$sibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.001, a1p = 0.168, b1 = 5.78, b1p = -3.54, b2 = -0.181, b2p = 0.051), weights = dbhWeight, control = list(maxiter = 50))
psmeHeightFromDiameter$weibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = dbhWeight)
psmeHeightFromDiameter$weibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64.1, a2 = -0.007, a2p = 1.05, a3 = 0.032, a3p = -0.193, b1 = -0.006, b1p = -0.001, b2 = 1.231, b2p = -0.070), weights = dbhWeight) # a2 not significant
psmeHeightFromDiameter$weibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a4 = 41.9, b1 = -0.013, b1p = 0.007, b2 = 0.94, b2p = -0.096), weights = dbhWeight) # a3p, a4p not significant
#confint_nlrob(psmeHeightFromDiameter$sharmaPartonBalPhysio, level = 0.99, weights = pmin(psme2016physio$DBH^if_else(psme2016physio$isPlantation, -1.2, -0.4), 1))

psmeHeightFromDiameter$chapmanRichards = get_height_error("Chapman-Richards", psmeHeightFromDiameter$chapmanRichards, psme2016)
psmeHeightFromDiameter$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L", psmeHeightFromDiameter$chapmanRichardsBal, psme2016)
psmeHeightFromDiameter$chapmanRichardsBalPhysio = get_height_error("Chapman-Richards BA+L physio", psmeHeightFromDiameter$chapmanRichardsBalPhysio, psme2016physio)
psmeHeightFromDiameter$chapmanRichardsBalRelHt = get_height_error("Chapman-Richards BA+L RelHt", psmeHeightFromDiameter$chapmanRichardsBalRelHt, psme2016)
psmeHeightFromDiameter$chapmanRichardsPhysio = get_height_error("Chapman-Richards physio", psmeHeightFromDiameter$chapmanRichardsPhysio, psme2016physio)
psmeHeightFromDiameter$curtis = get_height_error("Curtis", psmeHeightFromDiameter$curtis, psme2016)
psmeHeightFromDiameter$gam = get_height_error("GCV GAM", psmeHeightFromDiameter$gam, psme2016)
psmeHeightFromDiameter$gamBal = get_height_error("GCV GAM BA+L", psmeHeightFromDiameter$gamBal, psme2016)
#psmeHeightFromDiameter$gamBalPhysio = get_height_error("GCV GAM BA+L physio", psmeHeightFromDiameter$gamBalPhysio, psme2016physio) # slow
#psmeHeightFromDiameter$gamPhysio = get_height_error("GCV GAM physio", psmeHeightFromDiameter$gamPhysio, psme2016physio)
#save(psmeHeightFromDiameter$gamBalPhysio, psmeHeightFromDiameter$gamPhysio, file = "trees/height-diameter/HtDia PSME spline height.rdata")
load("trees/height-diameter/HtDia PSME spline height.rdata")
psmeHeightFromDiameter$hossfeld = get_height_error("Hossfeld IV", psmeHeightFromDiameter$hossfeld, psme2016)
psmeHeightFromDiameter$korf = get_height_error("Korf", psmeHeightFromDiameter$korf, psme2016)
psmeHeightFromDiameter$linear = get_height_error("linear", psmeHeightFromDiameter$linear, psme2016)
psmeHeightFromDiameter$michaelisMenten = get_height_error("Michaelis-Menten", psmeHeightFromDiameter$michaelisMenten, psme2016)
psmeHeightFromDiameter$parabolic = get_height_error("parabolic", psmeHeightFromDiameter$parabolic, psme2016)
psmeHeightFromDiameter$power = get_height_error("power", psmeHeightFromDiameter$power, psme2016)
psmeHeightFromDiameter$prodan = get_height_error("Prodan", psmeHeightFromDiameter$prodan, psme2016)
psmeHeightFromDiameter$ratkowsky = get_height_error("Ratkowsky", psmeHeightFromDiameter$ratkowsky, psme2016)
psmeHeightFromDiameter$richards = get_height_error("unified Richards", psmeHeightFromDiameter$richards, psme2016)
psmeHeightFromDiameter$sharmaParton = get_height_error("Sharma-Parton", psmeHeightFromDiameter$sharmaParton, psme2016)
psmeHeightFromDiameter$sharmaPartonBal = get_height_error("Sharma-Parton BA+L", psmeHeightFromDiameter$sharmaPartonBal, psme2016)
psmeHeightFromDiameter$sharmaPartonBalPhysio = get_height_error("Sharma-Parton BA+L physio", psmeHeightFromDiameter$sharmaPartonBalPhysio, psme2016physio)
psmeHeightFromDiameter$sharmaPartonPhysio = get_height_error("Sharma-Parton physio", psmeHeightFromDiameter$sharmaPartonPhysio, psme2016physio)
psmeHeightFromDiameter$sharmaZhang = get_height_error("Sharma-Zhang", psmeHeightFromDiameter$sharmaZhang, psme2016)
psmeHeightFromDiameter$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L", psmeHeightFromDiameter$sharmaZhangBal, psme2016)
psmeHeightFromDiameter$sibbesen = get_height_error("Sibbesen", psmeHeightFromDiameter$sibbesen, psme2016)
psmeHeightFromDiameter$weibull = get_height_error("Weibull", psmeHeightFromDiameter$weibull, psme2016)
psmeHeightFromDiameter$weibullBal = get_height_error("Weibull BA+L", psmeHeightFromDiameter$weibullBal, psme2016)
psmeHeightFromDiameter$weibullBalRelHt = get_height_error("Weibull BA+L RelHt", psmeHeightFromDiameter$weibullBalRelHt, psme2016)

psmeHeightFromDiameterNls$chapmanRichardsNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$chapmanRichardsBalNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$chapmanRichardsBalPhysioNls = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$chapmanRichardsPhysioNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$curtisNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = psmeHeightFromDiameter$curtis$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$hossfeldNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$hossfeld$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$korfNls = gsl_nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$korf$m$getPars(), weights = dbhWeight, control = nls.control(maxiter = 500))
psmeHeightFromDiameterNls$michaelisMentenNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = psmeHeightFromDiameter$michaelisMenten$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$prodanNls = gsl_nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = psmeHeightFromDiameter$prodan$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$powerNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = psmeHeightFromDiameter$power$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$ratkowskyNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$ratkowsky$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$richardsNls = gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = psmeHeightFromDiameter$richards$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaPartonNls = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaPartonBalNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaPartonBalPhysioNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaPartonPhysioNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonPhysio$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaZhangNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sharmaZhangBalNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$sibbesenNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$sibbesen$m$getPars(), weights = dbhWeight, control = list(maxiter = 50))
psmeHeightFromDiameterNls$weibullNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$m$getPars(), weights = dbhWeight)
psmeHeightFromDiameterNls$weibullBalNls = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$m$getPars(), weights = dbhWeight)

psmeHeightFromDiameterNls$chapmanRichardsNls = get_height_error("Chapman-Richards", psmeHeightFromDiameterNls$chapmanRichards, psme2016)
psmeHeightFromDiameterNls$chapmanRichardsBalNls = get_height_error("Chapman-Richards BA+L", psmeHeightFromDiameterNls$chapmanRichardsBal, psme2016)
psmeHeightFromDiameterNls$chapmanRichardsBalPhysioNls = get_height_error("Chapman-Richards BA+L physio", psmeHeightFromDiameterNls$chapmanRichardsBalPhysio, psme2016physio)
psmeHeightFromDiameterNls$chapmanRichardsPhysioNls = get_height_error("Chapman-Richards physio", psmeHeightFromDiameterNls$chapmanRichardsPhysio, psme2016physio)
psmeHeightFromDiameterNls$curtisNls = get_height_error("Curtis", psmeHeightFromDiameterNls$curtis, psme2016)
psmeHeightFromDiameterNls$hossfeldNls = get_height_error("Hossfeld IV", psmeHeightFromDiameterNls$hossfeld, psme2016)
psmeHeightFromDiameterNls$korfNls = get_height_error("Korf", psmeHeightFromDiameterNls$korf, psme2016)
psmeHeightFromDiameterNls$michaelisMentenNls = get_height_error("Michaelis-Menten", psmeHeightFromDiameterNls$michaelisMenten, psme2016)
psmeHeightFromDiameterNls$powerNls = get_height_error("power", psmeHeightFromDiameterNls$power, psme2016)
psmeHeightFromDiameterNls$prodanNls = get_height_error("Prodan", psmeHeightFromDiameterNls$prodan, psme2016)
psmeHeightFromDiameterNls$ratkowskyNls = get_height_error("Ratkowsky", psmeHeightFromDiameterNls$ratkowsky, psme2016)
psmeHeightFromDiameterNls$richardsNls = get_height_error("unified Richards", psmeHeightFromDiameterNls$richards, psme2016)
psmeHeightFromDiameterNls$sharmaPartonNls = get_height_error("Sharma-Parton", psmeHeightFromDiameterNls$sharmaParton, psme2016)
psmeHeightFromDiameterNls$sharmaPartonBalNls = get_height_error("Sharma-Parton BA+L", psmeHeightFromDiameterNls$sharmaPartonBal, psme2016)
psmeHeightFromDiameterNls$sharmaPartonBalPhysioNls = get_height_error("Sharma-Parton BA+L physio", psmeHeightFromDiameterNls$sharmaPartonBalPhysio, psme2016physio)
psmeHeightFromDiameterNls$sharmaPartonPhysioNls = get_height_error("Sharma-Parton physio", psmeHeightFromDiameterNls$sharmaPartonPhysio, psme2016physio)
psmeHeightFromDiameterNls$sharmaZhangNls = get_height_error("Sharma-Zhang", psmeHeightFromDiameterNls$sharmaZhang, psme2016)
psmeHeightFromDiameterNls$sharmaZhangBalNls = get_height_error("Sharma-Zhang BA+L", psmeHeightFromDiameterNls$sharmaZhangBal, psme2016)
psmeHeightFromDiameterNls$sibbesenNls = get_height_error("Sibbesen", psmeHeightFromDiameterNls$sibbesen, psme2016)
psmeHeightFromDiameterNls$weibullNls = get_height_error("Weibull", psmeHeightFromDiameterNls$weibull, psme2016)
psmeHeightFromDiameterNls$weibullBalNls = get_height_error("Weibull BA+L", psmeHeightFromDiameterNls$weibullBal, psme2016)

psmeHeightFromDiameterDefaultWeight = pmin(DBH^-1, 1)
psmeHeightFromDiameterNlsFwt = list(chapmanRichardsNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$m$getPars(), weights = pmin(DBH^-1, 1)))
psmeHeightFromDiameterNlsFwt$chapmanRichardsBalNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$chapmanRichardsBalPhysioNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsBalPhysio$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$chapmanRichardsPhysioNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$chapmanRichardsPhysio$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$curtisNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = psmeHeightFromDiameter$curtis$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$hossfeldNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$hossfeld$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$korfNlsFwt = gsl_nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$korf$m$getPars(), weights = pmin(DBH^-1, 1), control = nls.control(maxiter = 500))
psmeHeightFromDiameterNlsFwt$michaelisMentenNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = psmeHeightFromDiameter$michaelisMenten$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$prodanNlsFwt = gsl_nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = psmeHeightFromDiameter$prodan$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$powerNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = psmeHeightFromDiameter$power$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$ratkowskyNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$ratkowsky$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$richardsNlsFwt = gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = psmeHeightFromDiameter$richards$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaPartonNlsFwt = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaPartonBalNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(b1 + b1p * isPlantation) * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaPartonBalPhysioNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonBalPhysio$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaPartonPhysioNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a4 * elevation + a5 * slope + a6 * sin(3.14159/180 * aspect) + a7 * cos(3.14159/180 * aspect) + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp(b2*(tph/standBasalAreaPerHectare)^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016physio, start = psmeHeightFromDiameter$sharmaPartonPhysio$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaZhangNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sharmaZhangBalNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp((b2 + b2p * isPlantation)*tph^(b3 + b3p * isPlantation)*DBH))^(b4 + b4p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$sibbesenNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = psmeHeightFromDiameter$sibbesen$m$getPars(), weights = pmin(DBH^-1, 1), control = nls.control(maxiter = 100))
psmeHeightFromDiameterNlsFwt$weibullNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$m$getPars(), weights = pmin(DBH^-1, 1))
psmeHeightFromDiameterNlsFwt$weibullBalNlsFwt = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$m$getPars(), weights = pmin(DBH^-1, 1))

psmeHeightFromDiameterNlsFwt$chapmanRichardsNlsFwt = get_height_error("Chapman-Richards", psmeHeightFromDiameterNlsFwt$chapmanRichards, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$chapmanRichardsBalNlsFwt = get_height_error("Chapman-Richards BA+L", psmeHeightFromDiameterNlsFwt$chapmanRichardsBal, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$chapmanRichardsBalPhysioNlsFwt = get_height_error("Chapman-Richards BA+L physio", psmeHeightFromDiameterNlsFwt$chapmanRichardsBalPhysio, psme2016physio, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$chapmanRichardsPhysioNlsFwt = get_height_error("Chapman-Richards physio", psmeHeightFromDiameterNlsFwt$chapmanRichardsPhysio, psme2016physio, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$curtisNlsFwt = get_height_error("Curtis", psmeHeightFromDiameterNlsFwt$curtis, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$hossfeldNlsFwt = get_height_error("Hossfeld IV", psmeHeightFromDiameterNlsFwt$hossfeld, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$korfNlsFwt = get_height_error("Korf", psmeHeightFromDiameterNlsFwt$korf, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$michaelisMentenNlsFwt = get_height_error("Michaelis-Menten", psmeHeightFromDiameterNlsFwt$michaelisMenten, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$powerNlsFwt = get_height_error("power", psmeHeightFromDiameterNlsFwt$power, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$prodanNlsFwt = get_height_error("Prodan", psmeHeightFromDiameterNlsFwt$prodan, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$ratkowskyNlsFwt = get_height_error("Ratkowsky", psmeHeightFromDiameterNlsFwt$ratkowsky, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$richardsNlsFwt = get_height_error("unified Richards", psmeHeightFromDiameterNlsFwt$richards, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaPartonNlsFwt = get_height_error("Sharma-Parton", psmeHeightFromDiameterNlsFwt$sharmaParton, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaPartonBalNlsFwt = get_height_error("Sharma-Parton BA+L", psmeHeightFromDiameterNlsFwt$sharmaPartonBal, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaPartonBalPhysioNlsFwt = get_height_error("Sharma-Parton BA+L physio", psmeHeightFromDiameterNlsFwt$sharmaPartonBalPhysio, psme2016physio, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaPartonPhysioNlsFwt = get_height_error("Sharma-Parton physio", psmeHeightFromDiameterNlsFwt$sharmaPartonPhysio, psme2016physio, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaZhangNlsFwt = get_height_error("Sharma-Zhang", psmeHeightFromDiameterNlsFwt$sharmaZhang, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sharmaZhangBalNlsFwt = get_height_error("Sharma-Zhang BA+L", psmeHeightFromDiameterNlsFwt$sharmaZhangBal, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$sibbesenNlsFwt = get_height_error("Sibbesen", psmeHeightFromDiameterNlsFwt$sibbesen, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$weibullNlsFwt = get_height_error("Weibull", psmeHeightFromDiameterNlsFwt$weibull, psme2016, weights = psme2016$dbhWeightDefault)
psmeHeightFromDiameterNlsFwt$weibullBalNlsFwt = get_height_error("Weibull BA+L", psmeHeightFromDiameterNlsFwt$weibullBal, psme2016, weights = psme2016$dbhWeightDefault)

psmeHeightFromDiameterResults = bind_rows(as_row(psmeHeightFromDiameter$chapmanRichards),
                                          as_row(psmeHeightFromDiameter$chapmanRichardsBal),
                                          as_row(psmeHeightFromDiameter$chapmanRichardsBalPhysio),
                                          as_row(psmeHeightFromDiameter$chapmanRichardsBalRelHt),
                                          as_row(psmeHeightFromDiameter$chapmanRichardsPhysio),
                                          as_row(psmeHeightFromDiameter$curtis),
                                          as_row(psmeHeightFromDiameter$gam),
                                          as_row(psmeHeightFromDiameter$gamBal),
                                          as_row(psmeHeightFromDiameter$gamBalPhysio),
                                          as_row(psmeHeightFromDiameter$gamPhysio),
                                          as_row(psmeHeightFromDiameter$hossfeld),
                                          as_row(psmeHeightFromDiameter$korf),
                                          as_row(psmeHeightFromDiameter$linear),
                                          as_row(psmeHeightFromDiameter$michaelisMenten),
                                          as_row(psmeHeightFromDiameter$parabolic),
                                          as_row(psmeHeightFromDiameter$power),
                                          as_row(psmeHeightFromDiameter$prodan),
                                          as_row(psmeHeightFromDiameter$ratkowsky),
                                          as_row(psmeHeightFromDiameter$richards),
                                          as_row(psmeHeightFromDiameter$sharmaParton),
                                          as_row(psmeHeightFromDiameter$sharmaPartonBal),
                                          as_row(psmeHeightFromDiameter$sharmaPartonBalPhysio),
                                          as_row(psmeHeightFromDiameter$sharmaPartonPhysio),
                                          as_row(psmeHeightFromDiameter$sharmaZhang),
                                          as_row(psmeHeightFromDiameter$sharmaZhangBal),
                                          as_row(psmeHeightFromDiameter$sibbesen),
                                          as_row(psmeHeightFromDiameter$weibull),
                                          as_row(psmeHeightFromDiameter$weibullBal),
                                          as_row(psmeHeightFromDiameter$weibullBalRelHt),
                                          as_row(psmeHeightFromDiameterNls$chapmanRichardsNls),
                                          as_row(psmeHeightFromDiameterNls$chapmanRichardsBalNls),
                                          as_row(psmeHeightFromDiameterNls$chapmanRichardsBalPhysioNls),
                                          as_row(psmeHeightFromDiameterNls$chapmanRichardsPhysioNls),
                                          as_row(psmeHeightFromDiameterNls$curtis),
                                          as_row(psmeHeightFromDiameterNls$hossfeld),
                                          as_row(psmeHeightFromDiameterNls$korf),
                                          as_row(psmeHeightFromDiameterNls$michaelisMenten),
                                          as_row(psmeHeightFromDiameterNls$power),
                                          as_row(psmeHeightFromDiameterNls$prodan),
                                          as_row(psmeHeightFromDiameterNls$ratkowsky),
                                          as_row(psmeHeightFromDiameterNls$richards),
                                          as_row(psmeHeightFromDiameterNls$sharmaParton),
                                          as_row(psmeHeightFromDiameterNls$sharmaPartonBal),
                                          as_row(psmeHeightFromDiameterNls$sharmaPartonBalPhysio),
                                          as_row(psmeHeightFromDiameterNls$sharmaPartonPhysio),
                                          as_row(psmeHeightFromDiameterNls$sharmaZhang),
                                          as_row(psmeHeightFromDiameterNls$sharmaZhangBal),
                                          as_row(psmeHeightFromDiameterNls$sibbesen),
                                          as_row(psmeHeightFromDiameterNls$weibull),
                                          as_row(psmeHeightFromDiameterNls$weibullBal),
                                          as_row(psmeHeightFromDiameterNlsFwt$chapmanRichards, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$chapmanRichardsBal, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$chapmanRichardsBalPhysio, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$chapmanRichardsPhysio, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$curtis, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$hossfeld, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$korf, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$michaelisMenten, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$power, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$prodan, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$ratkowsky, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$richards, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaParton, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaPartonBal, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaPartonBalPhysio, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaPartonPhysio, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaZhang, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sharmaZhangBal, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$sibbesen, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$weibull, fixedWeight = -1),
                                          as_row(psmeHeightFromDiameterNlsFwt$weibullBal, fixedWeight = -1)) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(psmeHeightFromDiameter$results %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

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

psmeHeightFromDiameter$Efficiency = psmeHeightFromDiameter$results %>% filter(fitting %in% c("nlrob", "gsl_nls"), (fitting == "nlrob") | (is.na(fixedWeight) == FALSE), str_detect(name, "RelHt") == FALSE) %>% 
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
  labs(x = NULL, y = "nlrob() height residual, m") +
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
  

#psmeHeightFromDiameter$chapmanRichardsBalPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * Elev_Mean + a5 * SlopeMean + a6 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 87, a2 = 0.49, a3 = -0.032, a4 = -0.0025, a5 = -0.16, a6 = 3.63, b1 = -0.011, b2 = 1.09), weights = dbhWeight)
#psmeHeightFromDiameter$chapmanRichardsNatural = gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016natural, start = list(a1 = 66, b1 = -0.02, b2 = 1.5), weights = pmin(1/psme2016natural$DBH, 1))
#psmeHeightFromDiameter$chapmanRichardsPhysiographic = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * Elev_Mean + a3 * SlopeMean + a4 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = -0.002, a3 = -0.11, a4 = 1.78, b1 = -0.014, b2 = 1.13), weights = dbhWeight)
#psmeHeightFromDiameter$chapmanRichardsPlantation = gsl_nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016plantation, start = list(a1 = 51, b1 = -0.02, b2 = 1.2), weights = pmin(1/psme2016plantation$DBH, 1))
#psmeHeightFromDiameter$chapmanRichardsTopHeight = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 28.4, a2 = 0.21, b1 = -0.017, b2 = 1.08), weights = dbhWeight)
#psmeHeightFromDiameter$sharmaPartonBal = gsl_nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * basalAreaLarger) *(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 43, a2 = 0.13, a3 = 0.002, b1 = -0.017, b2 = -0.11, b3 = 1.00), weights = dbhWeight)
#psmeHeightFromDiameter$weibullBalNatural = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016natural, start = list(a1 = 62.2, a2 = 0.038, a3 = 0.025, b1 = -0.0045, b2 = 1.32), weights = pmin(1/psme2016natural$DBH, 1))
#psmeHeightFromDiameter$weibullBalPlantation = gsl_nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016plantation, start = list(a1 = 73.2, a2 = 0.93, a3 = -0.19, b1 = -0.008, b2 = 1.07), weights = pmin(1/psme2016plantation$DBH, 1))

#psmeHeightFromDiameter$chapmanRichardsBalPhysio = get_height_error(psmeHeightFromDiameter$chapmanRichardsBalPhysio, psme2016)
#psmeHeightFromDiameter$chapmanRichardsPhysiographic = get_height_error(psmeHeightFromDiameter$chapmanRichardsPhysiographic, psme2016)
#psmeHeightFromDiameter$chapmanRichardsNatural = get_height_error(psmeHeightFromDiameter$chapmanRichardsNatural, psme2016)
#psmeHeightFromDiameter$chapmanRichardsPlantation = get_height_error(psmeHeightFromDiameter$chapmanRichardsPlantation, psme2016)
#psmeHeightFromDiameter$chapmanRichardsTopHeight = get_height_error(psmeHeightFromDiameter$chapmanRichardsTopHeight, psme2016)
#psmeHeightFromDiameter$weibullBalNatural = get_height_error(psmeHeightFromDiameter$weibullBalNatural, psme2016natural)
#psmeHeightFromDiameter$weibullBalPlantation = get_height_error(psmeHeightFromDiameter$weibullBalPlantation, psme2016plantation)

## Douglas-fir height-diameter GNLS regressions
#psmeHeightFromDiameter$chapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameter$chapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameter$chapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 0.001, tolerance = 0.01, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05, iterations at default msTol and tolerance
#psmeHeightFromDiameter$sharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 500, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameter$sharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#psmeHeightFromDiameter$sharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, msTol = 1E-5, tolerance = 1E-4, returnObject = FALSE)) # step having at nlsTol = 0.1
#psmeHeightFromDiameter$sharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameter$sharmaZhangBal$m$getPars(), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, nlsMaxIter = 50, msTol = 1E-7, tolerance = 1E-6, msVerbose = FALSE, returnObject = FALSE)) # step halving factor at nlsTol = 1 with plot correlation, step halving with nlsTol = 0.07
#psmeHeightFromDiameter$weibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameter$weibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameter$weibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 5, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.2
#save(psmeHeightFromDiameter$chapmanRichardsGnls, psmeHeightFromDiameter$chapmanRichardsBalGnls, psmeHeightFromDiameter$sharmaPartonGnls, psmeHeightFromDiameter$sharmaPartonBalGnls, psmeHeightFromDiameter$sharmaZhangGnls, psmeHeightFromDiameter$sharmaZhangBalGnls, psmeHeightFromDiameter$weibullGnls, psmeHeightFromDiameter$weibullBalGnls, file = "trees/height-diameter/HtDia PSME GNLS.rdata")

load("trees/height-diameter/HtDia PSME GNLS.rdata")
psmeHeightFromDiameterGnls$chapmanRichards = get_height_error("Chapman-Richards GNLS", psmeHeightFromDiameterGnls$chapmanRichards, psme2016)
psmeHeightFromDiameterGnls$chapmanRichardsBal = get_height_error("Chapman-Richards BA+L GNLS", psmeHeightFromDiameterGnls$chapmanRichardsBal, psme2016)
psmeHeightFromDiameterGnls$sharmaParton = get_height_error("Sharma-Parton GNLS", psmeHeightFromDiameterGnls$sharmaParton, psme2016)
psmeHeightFromDiameterGnls$sharmaPartonBal = get_height_error("Sharma-Parton BA+L GNLS", psmeHeightFromDiameterGnls$sharmaPartonBal, psme2016)
psmeHeightFromDiameterGnls$sharmaZhang = get_height_error("Sharma-Zhang GNLS", psmeHeightFromDiameterGnls$sharmaZhang, psme2016)
psmeHeightFromDiameterGnls$sharmaZhangBal = get_height_error("Sharma-Zhang BA+L GNLS", psmeHeightFromDiameterGnls$sharmaZhangBal, psme2016)
psmeHeightFromDiameterGnls$weibull = get_height_error("Weibull GNLS", psmeHeightFromDiameterGnls$weibull, psme2016)
psmeHeightFromDiameterGnls$weibullBal = get_height_error("Weibull BA+L GNLS", psmeHeightFromDiameterGnls$weibullBal, psme2016)

psmeHeightFromDiameterResultsGnls = bind_rows(as_row(psmeHeightFromDiameterGnls$chapmanRichards),
                                              as_row(psmeHeightFromDiameterGnls$chapmanRichardsBal),
                                              as_row(psmeHeightFromDiameterGnls$sharmaParton),
                                              as_row(psmeHeightFromDiameterGnls$sharmaPartonBal),
                                              as_row(psmeHeightFromDiameterGnls$sharmaZhang),
                                              as_row(psmeHeightFromDiameterGnls$sharmaZhangBal),
                                              as_row(psmeHeightFromDiameterGnls$weibull),
                                              as_row(psmeHeightFromDiameterGnls$weibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
psmeHeightFromDiameter$resultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

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

## Douglas-fir diameter-height regressions
# Chapman-Richards, Gomperz, Prodan, Ratkowsky, and Korf run opposite the desired curvature
# Huang has DBH as ln of height (https://doi.org/10.1093/forestry/cpz002, Eq. 1)
# Hossfeld is sensitive to going negative, extrapolating poorly
# Prodan sometimes produces negative DBH
# Sharma-Parton and Sharma-Zhang have wrong curvature, simply modified forms fail to start or fail to step
# Sibbesen and Korf are near identical
#psmeDiameterFromHeight$chapmanRichards = gsl_nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999999)), psme2016, start = list(a1 = -65, b1 = 1/85, b2 = 0.75), weights = heightWeight)
#psmeDiameterFromHeight$chapmanFormBalRelHt = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 109, a1p = -40.0, a2 = -1.42, a3 = 0.490, a3p = 0.120, a4 = -45.1, a4p = 36.5, b1 = 0.041, b2 = 0.763, b2p = -0.0092), weights = heightWeight) # a2p not significant
#psmeDiameterFromHeight$chapmanRichardsPhysio = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016physio, start = list(a1 = -20.8, a1p = 5.59, a2 = 0.00218, a3 = -3.39, a4 = 0.209, a5 = 0.109, a6 = 0.101, b1 = 0.0176, b1p = 0.00776, b2 = 0.329), weights = heightWeight) # a5 not significant
#psmeDiameterFromHeight$michaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), psme2016, start = list(a1 = 108, a2 = 46.3, b1 = 0.77), weights = TotalHt^-2)
#psmeDiameterFromHeight$naslund = gsl_nls(DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 3.7, a2 = -0.09), weights = TotalHt^-2) # diverges if a2 start value is less than ~-0.1
#psmeDiameterFromHeight$power = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^b1, psme2016, start = list(a1 = 1.26, b1 = 1.08), weights = heightWeight)
#psmeDiameterFromHeight$prodan = gsl_nls(DBH ~ (TotalHt - 1.37)^2 / (a0 + a1 * (TotalHt - 1.37) + a2 * (TotalHt - 1.37)^2), psme2016, start = list(a0 = -0.7, a1 = 0.8, a2 = -0.004), weights = heightWeight) # AIC 
#psmeDiameterFromHeight$prodan = get_dbh_error(psmeDiameterFromHeight$prodan, psme2016)
#psmeDiameterFromHeight$sharmaParton = gsl_nls(DBH ~ a1*topHeight^a2*exp(b1*(tph/standBasalAreaPerHectare)^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 1, a2 = 1, b1 = -0.01, b2 = 1, b3 = 1), weights = heightWeight)
#psmeDiameterFromHeight$sharmaParton = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 57.1, a1p = -25.5, a2 = -0.241, a2p = 0.060, b1 = 0.0261, b2 = -0.0535, b2p = 0.0906, b3 = 0.689, b3p = 0.0829), weights = heightWeight)
#psmeDiameterFromHeight$sharmaParton = gsl_nls(DBH ~ a1*topHeight^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, psme2016, start = list(a1 = 39.3, a2 = 0.127, b1 = 0.0180, b2 = -0.010, b3 = 0.824), weights = heightWeight)
#psmeDiameterFromHeight$sharmaZhang = gsl_nls(DBH ~ a1*standBasalAreaPerHectare^a2*exp(b1*tph^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.0005, b2 = 0.5, b3 = 1), weights = heightWeight)
#psmeDiameterFromHeight$wykoff = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.081, a1p = -0.390, b1 = 1.576, b2 = 0.263, b2p = 0.0237), weights = heightWeight)
#psmeDiameterFromHeight$wykoff = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 15, b1 = 0.1, b2 = 0.35), weights = heightWeight) # NaN or inf 
#psmeDiameterFromHeight$wykoff = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 36.6, b1 = 0.059, b2 = 0.77), weights = heightWeight)
#psmeDiameterFromHeight$wykoff = gsl_nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 63.0, b1 = 0.018, b2 = 0.87), weights = heightWeight)
#psmeDiameterFromHeight$wykoff = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 17.8, a1p = -6.45, b1 = 0.216, b2 = 0.541, b2p = 0.048), weights = heightWeight)
#psmeDiameterFromHeight$wykoffBal = gsl_nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 121.7, a2 = -2.17, a3 = 1.04, b1 = 0.017, b2 = 0.85), weights = heightWeight, control = list(maxiter = 500))
#psmeDiameterFromHeight$wykoffBalRelHt = gsl_nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.340, a1p = -0.175, a2 = -0.0046, a2p = 0.00162, a3 = 0.00142, a3p = -0.00058, a4 = -0.132, a4p = 0.119, b1 = 3.040, b2 = 0.173, b2p = 0.0129), weights = heightWeight) # >500 iterations with b1p
#psmeDiameterFromHeight$wykoffRelHt = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.065, a1p = -0.386, a2 = 0.0071, b1 = 1.59, b2 = 0.262, b2p = 0.0234), weights = heightWeight)
#psmeDiameterFromHeight$weibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -204, b1 = 0.011, b2 = 0.87), weights = heightWeight)
#psmeDiameterFromHeight$weibull = gsl_nls(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -374, a1p = 163, b1 = 0.009, b1p = 0.003, b2 = 0.82), weights = heightWeight) # b2p not significant
#cor(cbind(dbh = psme2016$DBH, bal = psme2016$basalAreaLarger, height = psme2016$TotalHt, relHt = psme2016$relativeHeight, topHt = psme2016$topHeight, deltaHt = psme2016$TotalHt - psme2016$topHeight, treesPerPlot = psme2016$treesPerPlot, tph = psme2016$tph, tallerTph = psme2016$tallerTph * psme2016$topHeight, tallerQuasiBA = psme2016$tallerQuasiBA))
psmeDiameterFromHeight$chapmanForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), weights = heightWeight)
psmeDiameterFromHeight$chapmanFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.4, a1p = -44.6, a2 = 0.0058, a3 = -0.0544, b1 = 0.00166, b1p = 0.017, b2 = 0.788, b2p = -0.056), weights = heightWeight) # a2, a2p, a3p not significant
psmeDiameterFromHeight$chapmanFormBal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 135, a1p = -37.5, a2 = -1.2, b1 = 0.010, b1p = 0.002, b2 = 0.756, b2p = 0.064), weights = heightWeight, control = nls.control(maxiter = 500)) # a2p not significant, nlrob() step factor with a3 * BA
psmeDiameterFromHeight$chapmanFormBalRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 139, a1p = -38.5, a2 = -1.2, a9 = -5.5, a9p = 0.22, b1 = 0.012, b1p = 0.003, b2 = 0.796, b2p = 0.066), weights = heightWeight, control = nls.control(maxiter = 500)) # a2p not significant, nlrob() step factor with a9 * BA
psmeDiameterFromHeight$chapmanFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 21.1, a1p = -6.56, a9 = 0.278, b1 = 0.170, b2 = 0.580, b2p = 0.0395), weights = heightWeight) # a4p not significant
psmeDiameterFromHeight$chapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781), weights = heightWeight) # b2p not significant
psmeDiameterFromHeight$chapmanRichardsAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786), weights = heightWeight) # a2p, a3p, b2p not significant
psmeDiameterFromHeight$chapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -17.9, a1p = -0133, a4 = 0.0024, a5 = -7.189, a6 = 0.359, a7 = 0.170, a8 = 0.180, b1 = 0.0176, b1p = 0.0058, b2 = 0.404), maxit = 20, weights = heightWeight, control = nls.control(maxiter = 50)) # a7 not significant
psmeDiameterFromHeight$chapmanRichardsRelHt = nlrob(DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -1.22, a9 = -3.11, b1 = 0.018, b1p = 0.006, b2 = 0.003, b2p = 0.052), weights = heightWeight, maxit = 20) # a1p not significant
psmeDiameterFromHeight$gam = gam(DBH ~ s(TotalHt, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016, family = scat, select = TRUE, weights = heightWeights)
psmeDiameterFromHeight$gamAat = gam(DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016, select = TRUE, weights = heightWeights)
#psmeDiameterFromHeight$gamAatPhysio = gam(DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016physio, select = TRUE, weights = heightWeights) # warnings
#psmeDiameterFromHeight$gamPhysio = gam(DBH ~ s(TotalHt, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, by = as.factor(isPlantation), pc = psme2016gamConstraint), data = psme2016physio, select = TRUE, weights = heightWeights)
psmeDiameterFromHeight$linear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016, select = TRUE, weights = heightWeight)
psmeDiameterFromHeight$michaelisMentenForm = nlrob(DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), weights = heightWeight)
psmeDiameterFromHeight$naslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018), weights = heightWeight)
psmeDiameterFromHeight$parabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016, weights = heightWeight)
psmeDiameterFromHeight$power = nlrob(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108), weights = heightWeight)
psmeDiameterFromHeight$powerAat = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053), weights = heightWeight) # a3 not significant
psmeDiameterFromHeight$powerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.630, a1p = 0.284, a4 = 0.00001, a5 = -0.082, a6 = -0.019, a7 = 0.0028, a8 = -0.0018, b1 = 1.03, b1p = -0.102), weights = heightWeight) # a7 not significant
psmeDiameterFromHeight$powerRelHt = nlrob(DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a9 = 0.361, b1 = 0.943, b1p = -0.068), weights = heightWeight) # a1p, a4p not significant
psmeDiameterFromHeight$ruark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), weights = heightWeight)
#psmeDiameterFromHeight$schnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 18.6), weights = heightWeight) # converges from red alder values but fails to reconverge (singular gradient or NaN-inf with nls())
psmeDiameterFromHeight$sharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = list(a1 = 3.95, b1 = 0.681, b1p = -0.139, b2 = 0.097, b3 = -0.130, b3p = 0.157, b4 = 0.125, b4p = 0.0678), weights = heightWeight, control = list(maxiter = 200)) # a1p NaN-inf (not significant?), singular gradient with all relative height forms attempted
psmeDiameterFromHeight$sibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017), weights = heightWeight)
psmeDiameterFromHeight$sibbesenFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, a3p = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190), weights = heightWeight) # a2, a2p not significant
psmeDiameterFromHeight$sibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = list(a1 = 3.812, a1p = -0.925, a4 = 0.00022, a5 = 0.0663, a6 = -0.0495, a7 = -0.0151, a8 = -0.00630, b1 = 0.520, b2 = 0.111, b2p = 0.0173), weights = heightWeight) # a5, a7 not significant
psmeDiameterFromHeight$sibbesenFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.90, a1p = -0.95, a9 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016), weights = heightWeight) # a4p not significant
psmeDiameterFromHeight$weibull = nlrob(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81), weights = heightWeight) # b2p not significant
#confint_nlrob(psmeDiameterFromHeight$sibbesenFormRelHt, level = 0.99, weights = pmin(psme2016$TotalHt^if_else(psme2016$isPlantation, -1.6, -1.7), 0.5))

psmeDiameterFromHeight$chapmanForm = get_dbh_error("Chapman-Richards form", psmeDiameterFromHeight$chapmanForm, psme2016)
psmeDiameterFromHeight$chapmanFormAat = get_dbh_error("Chapman-Richards form AA+T", psmeDiameterFromHeight$chapmanFormAat, psme2016)
psmeDiameterFromHeight$chapmanFormBal = get_dbh_error("Chapman-Richards form BA+L", psmeDiameterFromHeight$chapmanFormBal, psme2016)
psmeDiameterFromHeight$chapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BA+L RelHt", psmeDiameterFromHeight$chapmanFormBalRelHt, psme2016)
psmeDiameterFromHeight$chapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", psmeDiameterFromHeight$chapmanFormRelHt, psme2016)
psmeDiameterFromHeight$chapmanRichards = get_dbh_error("Chapman-Richards", psmeDiameterFromHeight$chapmanRichards, psme2016)
psmeDiameterFromHeight$chapmanRichardsAat = get_dbh_error("Chapman-Richards AA+T", psmeDiameterFromHeight$chapmanRichardsAat, psme2016)
psmeDiameterFromHeight$chapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", psmeDiameterFromHeight$chapmanRichardsPhysio, psme2016physio)
psmeDiameterFromHeight$chapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", psmeDiameterFromHeight$chapmanRichardsRelHt, psme2016)
psmeDiameterFromHeight$gam = get_dbh_error("GCV GAM", psmeDiameterFromHeight$gam, psme2016)
psmeDiameterFromHeight$gamAat = get_dbh_error("GCV GAM AA+T", psmeDiameterFromHeight$gamAat, psme2016)
#psmeDiameterFromHeight$gamAatPhysio = get_dbh_error("GCV GAM AA+T physio", psmeDiameterFromHeight$gamAatPhysio, psme2016)
#psmeDiameterFromHeight$gamPhysio = get_dbh_error("GCV GAM physio", psmeDiameterFromHeight$gamPhysio, psme2016)
#save(psmeDiameterFromHeight$gamAatPhysio, psmeDiameterFromHeight$gamPhysio, file = "trees/height-diameter/HtDia PSME spline DBH.rdata")
load("trees/height-diameter/HtDia PSME spline DBH.rdata")
psmeDiameterFromHeight$linear = get_dbh_error("linear", psmeDiameterFromHeight$linear, psme2016)
psmeDiameterFromHeight$michaelisMentenForm = get_dbh_error("Michaelis-Menten form", psmeDiameterFromHeight$michaelisMentenForm, psme2016)
psmeDiameterFromHeight$naslund = get_dbh_error("Näslund", psmeDiameterFromHeight$naslund, psme2016)
psmeDiameterFromHeight$parabolic = get_dbh_error("parabolic", psmeDiameterFromHeight$parabolic, psme2016)
psmeDiameterFromHeight$power = get_dbh_error("power", psmeDiameterFromHeight$power, psme2016)
psmeDiameterFromHeight$powerAat = get_dbh_error("power AA+T", psmeDiameterFromHeight$powerAat, psme2016)
psmeDiameterFromHeight$powerPhysio = get_dbh_error("power physio", psmeDiameterFromHeight$powerPhysio, psme2016physio)
psmeDiameterFromHeight$powerRelHt = get_dbh_error("power RelHt", psmeDiameterFromHeight$powerRelHt, psme2016)
psmeDiameterFromHeight$ruark = get_dbh_error("Ruark", psmeDiameterFromHeight$ruark, psme2016)
#psmeDiameterFromHeight$schnute = get_dbh_error("Schnute", psmeDiameterFromHeight$schnute, psme2016)
psmeDiameterFromHeight$sharmaParton = get_dbh_error("modified Sharma-Parton", psmeDiameterFromHeight$sharmaParton, psme2016)
psmeDiameterFromHeight$sibbesenForm = get_dbh_error("Sibbesen form", psmeDiameterFromHeight$sibbesenForm, psme2016)
psmeDiameterFromHeight$sibbesenFormAat = get_dbh_error("Sibbesen form AA+T", psmeDiameterFromHeight$sibbesenFormAat, psme2016)
psmeDiameterFromHeight$sibbesenFormPhysio = get_dbh_error("Sibbesen form physio", psmeDiameterFromHeight$sibbesenFormPhysio, psme2016physio)
psmeDiameterFromHeight$sibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", psmeDiameterFromHeight$sibbesenFormRelHt, psme2016)
psmeDiameterFromHeight$weibull = get_dbh_error("Weibull", psmeDiameterFromHeight$weibull, psme2016)

psmeDiameterFromHeight$chapmanFormNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanForm$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$chapmanFormAatNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = psmeDiameterFromHeight$chapmanFormAat$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$chapmanFormRelHtNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = psmeDiameterFromHeight$chapmanFormRelHt$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$chapmanRichardsNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichards$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$chapmanRichardsAatNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + a3 * standBasalAreaApprox)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsAat$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$chapmanRichardsPhysioNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = psmeDiameterFromHeight$chapmanRichardsPhysio$m$getPars(), weights = heightWeight, control = nls.control(maxiter = 50))
psmeDiameterFromHeight$chapmanRichardsRelHtNls = gsl_nls(DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = psmeDiameterFromHeight$chapmanRichardsRelHt$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$michaelisMentenFormNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = psmeDiameterFromHeight$michaelisMentenForm$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$naslundNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$naslund$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$powerNls = gsl_nls(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$power$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$powerAatNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerApproxBasalArea + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = psmeDiameterFromHeight$powerAat$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$powerPhysioNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = psmeDiameterFromHeight$powerPhysio$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$powerRelHtNls = gsl_nls(DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = psmeDiameterFromHeight$powerRelHt$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$ruarkNls = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = psmeDiameterFromHeight$ruark$m$getPars(), weights = heightWeight)
#psmeDiameterFromHeight$schnuteNls = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = psmeDiameterFromHeight$schnute$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$sharmaPartonNls = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation)*(exp(b2*(tph/topHeight)^(b3 + b3p * isPlantation)*(TotalHt - 1.37)) - 1)^(b4 + b4p * isPlantation), psme2016, start = psmeDiameterFromHeight$sharmaParton$m$getPars(), weights = heightWeight, control = list(maxiter = 200))
psmeDiameterFromHeight$sibbesenFormNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenForm$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$sibbesenFormAatNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerApproxBasalArea + (a3 + a3p * isPlantation) * standBasalAreaApprox)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormAat$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$sibbesenFormPhysioNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a4 * elevation + a5 * sin(3.14159/180 * slope) + a6 * cos(3.14159/180 * aspect) + a7 * sin(3.14159/180 * aspect) + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = psmeDiameterFromHeight$sibbesenFormPhysio$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$sibbesenFormRelHtNls = gsl_nls(DBH ~ (a1 + a1p * isPlantation + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = psmeDiameterFromHeight$sibbesenFormRelHt$m$getPars(), weights = heightWeight)
psmeDiameterFromHeight$weibullNls = gsl_nls(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = psmeDiameterFromHeight$weibull$m$getPars(), weights = heightWeight)

psmeDiameterFromHeight$chapmanFormNls = get_dbh_error("Chapman-Richards form", psmeDiameterFromHeight$chapmanFormNls, psme2016)
psmeDiameterFromHeight$chapmanFormAatNls = get_dbh_error("Chapman-Richards form AA+T", psmeDiameterFromHeight$chapmanFormAatNls, psme2016)
psmeDiameterFromHeight$chapmanFormRelHtNls = get_dbh_error("Chapman-Richards form RelHt", psmeDiameterFromHeight$chapmanFormRelHtNls, psme2016)
psmeDiameterFromHeight$chapmanRichardsNls = get_dbh_error("Chapman-Richards", psmeDiameterFromHeight$chapmanRichardsNls, psme2016)
psmeDiameterFromHeight$chapmanRichardsAatNls = get_dbh_error("Chapman-Richards AA+T", psmeDiameterFromHeight$chapmanRichardsAatNls, psme2016)
psmeDiameterFromHeight$chapmanRichardsPhysioNls = get_dbh_error("Chapman-Richards physio", psmeDiameterFromHeight$chapmanRichardsPhysioNls, psme2016physio)
psmeDiameterFromHeight$chapmanRichardsRelHtNls = get_dbh_error("Chapman-Richards RelHt", psmeDiameterFromHeight$chapmanRichardsRelHtNls, psme2016)
psmeDiameterFromHeight$michaelisMentenFormNls = get_dbh_error("Michaelis-Menten form", psmeDiameterFromHeight$michaelisMentenFormNls, psme2016)
psmeDiameterFromHeight$naslundNls = get_dbh_error("Näslund", psmeDiameterFromHeight$naslundNls, psme2016)
psmeDiameterFromHeight$powerNls = get_dbh_error("power", psmeDiameterFromHeight$powerNls, psme2016)
psmeDiameterFromHeight$powerAatNls = get_dbh_error("power AA+T", psmeDiameterFromHeight$powerAatNls, psme2016)
psmeDiameterFromHeight$powerPhysioNls = get_dbh_error("power physio", psmeDiameterFromHeight$powerPhysioNls, psme2016physio)
psmeDiameterFromHeight$powerRelHtNls = get_dbh_error("power RelHt", psmeDiameterFromHeight$powerRelHtNls, psme2016)
psmeDiameterFromHeight$ruarkNls = get_dbh_error("Ruark", psmeDiameterFromHeight$ruarkNls, psme2016)
#psmeDiameterFromHeight$schnuteNls = get_dbh_error("Schnute", psmeDiameterFromHeight$schnuteNls, psme2016)
psmeDiameterFromHeight$sharmaPartonNls = get_dbh_error("modified Sharma-Parton", psmeDiameterFromHeight$sharmaPartonNls, psme2016)
psmeDiameterFromHeight$sibbesenFormNls = get_dbh_error("Sibbesen form", psmeDiameterFromHeight$sibbesenFormNls, psme2016)
psmeDiameterFromHeight$sibbesenFormAatNls = get_dbh_error("Sibbesen form AA+T", psmeDiameterFromHeight$sibbesenFormAatNls, psme2016)
psmeDiameterFromHeight$sibbesenFormPhysioNls = get_dbh_error("Sibbesen form physio", psmeDiameterFromHeight$sibbesenFormPhysioNls, psme2016physio)
psmeDiameterFromHeight$sibbesenFormRelHtNls = get_dbh_error("Sibbesen form RelHt", psmeDiameterFromHeight$sibbesenFormRelHtNls, psme2016)
psmeDiameterFromHeight$weibullNls = get_dbh_error("Weibull", psmeDiameterFromHeight$weibullNls, psme2016)

psmeDiameterFromHeightResults = bind_rows(as_row(psmeDiameterFromHeight$chapmanRichards),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsAat),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsPhysio),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsRelHt),
                                          as_row(psmeDiameterFromHeight$chapmanForm),
                                          as_row(psmeDiameterFromHeight$chapmanFormAat),
                                          as_row(psmeDiameterFromHeight$chapmanFormBal),
                                          as_row(psmeDiameterFromHeight$chapmanFormBalRelHt),
                                          as_row(psmeDiameterFromHeight$chapmanFormRelHt),
                                          as_row(psmeDiameterFromHeight$gam),
                                          as_row(psmeDiameterFromHeight$gamAat),
                                          as_row(psmeDiameterFromHeight$gamAatPhysio),
                                          as_row(psmeDiameterFromHeight$gamPhysio),
                                          as_row(psmeDiameterFromHeight$linear),
                                          as_row(psmeDiameterFromHeight$michaelisMentenForm),
                                          as_row(psmeDiameterFromHeight$naslund),
                                          as_row(psmeDiameterFromHeight$parabolic),
                                          as_row(psmeDiameterFromHeight$power),
                                          as_row(psmeDiameterFromHeight$powerAat),
                                          as_row(psmeDiameterFromHeight$powerPhysio),
                                          as_row(psmeDiameterFromHeight$powerRelHt),
                                          as_row(psmeDiameterFromHeight$ruark),
                                          as_row(name = "Schnute"),
                                          as_row(psmeDiameterFromHeight$sharmaParton),
                                          as_row(psmeDiameterFromHeight$sibbesenForm),
                                          as_row(psmeDiameterFromHeight$sibbesenFormAat),
                                          as_row(psmeDiameterFromHeight$sibbesenFormPhysio),
                                          as_row(psmeDiameterFromHeight$sibbesenFormRelHt),
                                          as_row(psmeDiameterFromHeight$weibull),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsNls),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsAatNls),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsPhysioNls),
                                          as_row(psmeDiameterFromHeight$chapmanRichardsRelHtNls),
                                          as_row(psmeDiameterFromHeight$chapmanFormNls),
                                          as_row(psmeDiameterFromHeight$chapmanFormAatNls),
                                          as_row(psmeDiameterFromHeight$chapmanFormRelHtNls),
                                          as_row(psmeDiameterFromHeight$michaelisMentenFormNls),
                                          as_row(psmeDiameterFromHeight$naslundNls),
                                          as_row(psmeDiameterFromHeight$powerNls),
                                          as_row(psmeDiameterFromHeight$powerAatNls),
                                          as_row(psmeDiameterFromHeight$powerPhysioNls),
                                          as_row(psmeDiameterFromHeight$powerRelHtNls),
                                          as_row(psmeDiameterFromHeight$ruarkNls),
                                          as_row(psmeDiameterFromHeight$sharmaPartonNls),
                                          as_row(psmeDiameterFromHeight$sibbesenFormNls),
                                          as_row(psmeDiameterFromHeight$sibbesenFormAatNls),
                                          as_row(psmeDiameterFromHeight$sibbesenFormPhysioNls),
                                          as_row(psmeDiameterFromHeight$sibbesenFormRelHtNls),
                                          as_row(psmeDiameterFromHeight$weibullNls)) %>% 
  mutate(responseVariable = "DBH", species = "PSME", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(psmeDiameterFromHeight$results %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(psme2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanForm), y = TotalHt, color = "Chapman form", group = isPlantation)) +
  #geom_line(aes(x = predict(psmeDiameterFromHeight$chapmanFormAat), y = TotalHt, color = "Chapman form AA+T", group = isPlantation), alpha = 0.5) +
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

psmeDiameterFromHeight$Efficiency = psmeDiameterFromHeight$results %>% filter(fitting %in% c("nlrob", "nls"), str_detect(name, "BA\\+L") == FALSE) %>% 
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

## Q-Q plots
#ggplot() + # symmetric t fits less well than skewed t, as expected
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'t(df = 7, '*alpha*' = 2.25)'", hjust = 0, parse = TRUE, size = 3.5) +
#  coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
#  labs(x = "theoretical quantile", y = NULL, color = NULL) +
#  scale_color_manual(values = dbhColors) +
#  theme(legend.justification = c(1, 0), legend.position = "none")
#ggplot() + # slow! and fits poorly
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanRichards), color = "Chapman-Richards"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$chapmanForm), color = "Chapman-Richards form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$ruark), color = "Ruark"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = -residuals(psmeDiameterFromHeight$sibbesenForm), color = "Sibbesen form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'pearsonIV(m = 1.72, '*nu*' = 0.656)'", hjust = 0, parse = TRUE, size = 3.5) +
#  coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
#  labs(x = "theoretical quantile", y = NULL, color = NULL) +
#  scale_color_manual(values = dbhColors) +
#  theme(legend.justification = c(1, 0), legend.position = "none")


## collect model parameters
psmeParameters = bind_rows(bind_rows(get_coefficients(psmeHeightFromDiameter$chapmanRichards),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsBal),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsBalPhysio),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsBalRelHt),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsPhysio),
                                     get_coefficients(psmeHeightFromDiameter$curtis),
                                     get_coefficients(psmeHeightFromDiameter$gam),
                                     get_coefficients(psmeHeightFromDiameter$gamBal),
                                     get_coefficients(psmeHeightFromDiameter$gamBalPhysio),
                                     get_coefficients(psmeHeightFromDiameter$gamPhysio),
                                     get_coefficients(psmeHeightFromDiameter$hossfeld),
                                     get_coefficients(psmeHeightFromDiameter$korf),
                                     get_coefficients(psmeHeightFromDiameter$linear),
                                     get_coefficients(psmeHeightFromDiameter$michaelisMenten),
                                     get_coefficients(psmeHeightFromDiameter$parabolic),
                                     get_coefficients(psmeHeightFromDiameter$power),
                                     get_coefficients(psmeHeightFromDiameter$prodan),
                                     get_coefficients(psmeHeightFromDiameter$ratkowsky),
                                     get_coefficients(psmeHeightFromDiameter$richards),
                                     get_coefficients(psmeHeightFromDiameter$sharmaParton),
                                     get_coefficients(psmeHeightFromDiameter$sharmaPartonBal),
                                     get_coefficients(psmeHeightFromDiameter$sharmaPartonBalPhysio),
                                     get_coefficients(psmeHeightFromDiameter$sharmaPartonPhysio),
                                     get_coefficients(psmeHeightFromDiameter$sharmaZhang),
                                     get_coefficients(psmeHeightFromDiameter$sharmaZhangBal),
                                     get_coefficients(psmeHeightFromDiameter$sibbesen),
                                     get_coefficients(psmeHeightFromDiameter$weibull),
                                     get_coefficients(psmeHeightFromDiameter$weibullBal),
                                     get_coefficients(psmeHeightFromDiameter$weibullBalRelHt),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsGnls),
                                     get_coefficients(psmeHeightFromDiameter$chapmanRichardsBalGnls),
                                     get_coefficients(psmeHeightFromDiameter$sharmaPartonGnls),
                                     get_coefficients(psmeHeightFromDiameter$sharmaPartonBalGnls),
                                     get_coefficients(psmeHeightFromDiameter$sharmaZhangGnls),
                                     get_coefficients(psmeHeightFromDiameter$sharmaZhangBalGnls),
                                     get_coefficients(psmeHeightFromDiameter$weibullGnls),
                                     get_coefficients(psmeHeightFromDiameter$weibullBalGnls)) %>%
                              mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(psmeDiameterFromHeight$chapmanRichards),
                                     get_coefficients(psmeDiameterFromHeight$chapmanRichardsAat),
                                     get_coefficients(psmeDiameterFromHeight$chapmanRichardsPhysio),
                                     get_coefficients(psmeDiameterFromHeight$chapmanRichardsRelHt),
                                     get_coefficients(psmeDiameterFromHeight$chapmanForm),
                                     get_coefficients(psmeDiameterFromHeight$chapmanFormAat),
                                     get_coefficients(psmeDiameterFromHeight$chapmanFormBal),
                                     get_coefficients(psmeDiameterFromHeight$chapmanFormBalRelHt),
                                     get_coefficients(psmeDiameterFromHeight$chapmanFormRelHt),
                                     get_coefficients(psmeDiameterFromHeight$gam),
                                     get_coefficients(psmeDiameterFromHeight$gamAat),
                                     get_coefficients(psmeDiameterFromHeight$gamAatPhysio),
                                     get_coefficients(psmeDiameterFromHeight$gamPhysio),
                                     get_coefficients(psmeDiameterFromHeight$linear),
                                     get_coefficients(psmeDiameterFromHeight$michaelisMentenForm),
                                     get_coefficients(psmeDiameterFromHeight$naslund),
                                     get_coefficients(psmeDiameterFromHeight$parabolic),
                                     get_coefficients(psmeDiameterFromHeight$power),
                                     get_coefficients(psmeDiameterFromHeight$powerAat),
                                     get_coefficients(psmeDiameterFromHeight$powerPhysio),
                                     get_coefficients(psmeDiameterFromHeight$powerRelHt),
                                     get_coefficients(psmeDiameterFromHeight$ruark),
                                     #get_coefficients(psmeDiameterFromHeight$schnute),
                                     get_coefficients(psmeDiameterFromHeight$sharmaParton),
                                     get_coefficients(psmeDiameterFromHeight$sibbesenForm),
                                     get_coefficients(psmeDiameterFromHeight$sibbesenFormAat),
                                     get_coefficients(psmeDiameterFromHeight$sibbesenFormPhysio),
                                     get_coefficients(psmeDiameterFromHeight$sibbesenFormRelHt),
                                     get_coefficients(psmeDiameterFromHeight$weibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "PSME", 
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, fitting, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
# essentially no difference between gsl_nls() and nlrob() fits
# Chapman-Richards has the wrong curvature
#psmeBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1, psme2016, start = list(a1 = 1, a1p = 0, b1 = 0.00009, b1p = 0, b2 = 2.2, b2p = 0), weights = pmin(1/basalArea, 1E4)) # a1p not significant
#psmeBasalAreaFromHeightPower = gsl_nls(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 0.25 * pi * 0.01^2, a1p = 0, b1 = 2.4, b1p = 0), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
psmeBasalAreaFromHeightKorf = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 0.689, a1p = -0.413, b1 = 0.0003, b1p = 0.0005, b2 = 1.91, b2p = -0.10), weights = pmin(1/basalArea, 1E4)) # a1p not significant
psmeBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = 0.00005, b1 = 2.41, b1p = -0.248), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
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


## exploratory plots
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
ggsave("Presentation/Douglas-fir height-diameter natural-plantation.png", width = 12.5, height = 11, units = "cm")
