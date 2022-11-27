# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## Douglas-fir height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Parton, Sharma-Zhang BAL, Chapman-Richards BAL
# fits without isPlantation as a factor
#psmeHeightFromDiameterChapmanRichards = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 74, b1 = -0.014, b2 = 1.14), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = 0.49, a3 = -0.053, b1 = -0.011, b2 = 1.09), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterChapmanRichardsBalRelHt = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 35, a2 = 0.70, a3 = 0.055, a4 = 38.32, b1 = -0.010, b2 = 0.92), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterChapmanRichardsPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -9.56, a2 = 0.0003, a2p = -0.011, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^-2, 1)) # b1p not significant
#psmeHeightFromDiameterChapmanRichardsPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation)*elevation + a3 * sin(pi/180 * slope)) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 66.7, a1p = -9.66, a2 = 0.0003, a2p = -0.011, a3 = -2.06, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = pmin(DBH^-2, 1)) # a3p not significant, little sensitivity to sin() or cos() of slope but sin() slightly more accurate
#psmeHeightFromDiameterChapmanRichardsTopo = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * topographicShelterIndex) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, a2 = 0, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = pmin(DBH^-2, 1)) # accuracy increases through cos(), sin(), linear but differences are limited
#psmeHeightFromDiameterCurtis = nls(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b2, psme2016, start = list(a1 = 1.6, b2 = 0.24), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterHossfeld = nls(TotalHt ~ 1.37 + a1 / (1 + b1*DBH^b2), psme2016, start = list(a1 = 99, b1 = 187, b2 = -1.2), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, a1p = 376, b1 = -7.83, b2 = -0.323, b2p = 0.084), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1), control = nls.control(maxiter = 500)) # b1p not significant
#psmeHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + a1*DBH / (a2 + DBH), psme2016, start = list(a1 = 138, a2 = 156), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterProdan = nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), psme2016, start = list(a1 = 0.008, a2 = 1.0, a3 = 2.7), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterPower = nls(TotalHt ~ 1.37 + b0*DBH^b1, psme2016, start = list(b0 = 1.54, b1 = 0.77), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRatkowsky = nls(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), psme2016, start = list(a1 = 93, b1 = -65, b2 = 15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 58.3, d = 0.609, kU = 0.0134), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 62.6, Hap = -26.5, d = 0.609, kU = 0.0128, kUp = 0.0114), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 55.4, a2 = 0.082, b1 = -0.022, b2 = -0.17, b3 = 1.04), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, psme2016, start = list(a1 = 39, a2 = 0.15, b1 = -0.018, b2 = -0.16, b3 = 1.01), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 55, a2 = 0.07, b1 = -0.012, b2 = 0.03, b3 = 1.1), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaZhangBal = nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 66, a2 = 0.056, a3 = 0.006, b1 = -0.022, b2 = -0.13, b3 = 1.04), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), psme2016, start = list(a1 = 0.034, b1 = 3.1, b2 = -0.15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterKorf = nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), psme2016, start = list(a1 = 172, b1 = -9.682, b2 = -0.4574), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWeibull = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 69, b1 = -0.0075, b2 = 1.15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWeibullBal = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 73.7, a2 = 0.38, a3 = -0.007, b1 = -0.008, b2 = 1.09), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation + b2 * pmin(relativeHeight, 1.25))*DBH^(b3 + b3p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a3p = 0.217, a4 = 41.9, a4p = 76.9, b1 = -0.013, b1p = 0.007, b2 = 0, b3 = 0.94, b3p = -0.096), weights = pmin(DBH^-2, 1), control = nls.control(maxiter = 500)) # step factor with nlrob() 
# fits with isPlantation
#psmeHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 68.6, a1p = -3.2, a2 = 0.024, a2p = 0.96, a3 = 0.014, a3p = -0.15, b1 = -0.018, b1p = 0.0028, b2 = 1.29, b2p = -0.10), weights = pmin(DBH^-2, 1))
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, TotalHt > 0) # live Douglas-firs measured for height
psme2016natural = psme2016 %>% filter(isPlantation == FALSE)
psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016plantation = psme2016 %>% filter(isPlantation)
psme2016plantationPhysio = psme2016physio %>% filter(isPlantation)

psmeHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # b1p not significant
psmeHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 72.9, a1p = -11.8, a2 = 0.087, a2p = 0.84, a3 = -0.0021, a3p = -0.073, b1 = -0.016, b2 = 1.26, b2p = -0.054), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a3 not significant, step factor with b1p
psmeHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 74.3, a2 = 0.096, a2p = 0.92, a3 = -0.015, a4 = -0.101, a5 = 0.793, a6 = 1.695, a7 = 0.183, b1 = -0.018, b1p = 0.005, b2 = 1.30, b2p = -0.154), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a4 not significant
psmeHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a4 = 54.0, a4p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016physio, start = list(a1 = 68.5, a1p = -13.4, a2 = -0.0045, a3 = -8.09, a4 = 0.783, a5 = 0.766, a6 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a4p not significant, a5p induces overfitting
psmeHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + a1*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, b1 = -7.83, b2 = -0.323, b2p = 0.084), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1), control = nls.control(maxiter = 500)) # a1p parameter evaporation, b1p not significant
psmeHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # b1p not significant
psmeHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), psme2016, offset = breastHeight, weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.15, a1p = -0.422, b1 = 0.85, b1p = 0.14), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 37.66, a2 = 0.19, a2p = -0.123, b1 = -0.017, b1p = -0.026, b2 = 0.061, b2p = -0.259, b3 = 1.33, b3p = -0.22), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016physio, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016physio, start = list(a1 = 38.4, a2 = 0.21, a2p = -0.122, a3 = -0.0002, a4 = 0.012, a5 = 0.007, a6 = 0.0006, b1 = -0.018, b1p = -0.024, b2 = 0.033, b2p = -0.225, b3 = 1.32, b3p = -0.21), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # nlrob() step factor with a1p
psmeHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a2 not significant
psmeHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 53.5, a1p = -10.1, a2 = 0.070, a3 = -0.00003, a3p = 0.012, b1 = -0.027, b1p = -0.017, b2 = -0.064, b2p = -0.096, b3 = 1.28, b3p = -0.16), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a2p, a3, b2 not significant
psmeHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.001, a1p = 0.168, b1 = 5.78, b1p = -3.54, b2 = -0.181, b2p = 0.051), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1), control = list(maxiter = 50))
psmeHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
psmeHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 64.1, a2 = -0.007, a2p = 1.05, a3 = 0.032, a3p = -0.193, b1 = -0.006, b1p = -0.001, b2 = 1.231, b2p = -0.070), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a2 not significant
psmeHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 35.2, a2 = 0.10, a2p = 2.05, a3 = -0.035, a3p = 0.217, a4 = 41.9, a4p = 76.9, b1 = -0.013, b1p = 0.007, b2 = 0.94, b2p = -0.096), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1)) # a2 not significant
#confint2(psmeHeightFromDiameterWeibullBalRelHt, level = 0.99)

psmeHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", psmeHeightFromDiameterChapmanRichards, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BAL", psmeHeightFromDiameterChapmanRichardsBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BAL physio", psmeHeightFromDiameterChapmanRichardsBalPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BAL RelHt", psmeHeightFromDiameterChapmanRichardsBalRelHt, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", psmeHeightFromDiameterChapmanRichardsPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterCurtis = get_height_error("Curtis", psmeHeightFromDiameterCurtis, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", psmeHeightFromDiameterHossfeld, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterKorf = get_height_error("Korf", psmeHeightFromDiameterKorf, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterLinear = get_height_error("linear", psmeHeightFromDiameterLinear, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", psmeHeightFromDiameterMichaelisMenten, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterParabolic = get_height_error("parabolic", psmeHeightFromDiameterParabolic, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterPower = get_height_error("power", psmeHeightFromDiameterPower, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterProdan = get_height_error("Prodan", psmeHeightFromDiameterProdan, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", psmeHeightFromDiameterRatkowsky, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterRichards = get_height_error("unified Richards", psmeHeightFromDiameterRichards, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", psmeHeightFromDiameterSharmaParton, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BAL", psmeHeightFromDiameterSharmaPartonBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BAL physio", psmeHeightFromDiameterSharmaPartonBalPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", psmeHeightFromDiameterSharmaPartonPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", psmeHeightFromDiameterSharmaZhang, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BAL", psmeHeightFromDiameterSharmaZhangBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSibbesen = get_height_error("Sibbesen", psmeHeightFromDiameterSibbesen, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibull = get_height_error("Weibull", psmeHeightFromDiameterWeibull, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBal = get_height_error("Weibull BAL", psmeHeightFromDiameterWeibullBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BAL RelHt", psmeHeightFromDiameterWeibullBalRelHt, psme2016, psme2016natural, psme2016plantation)

psmeHeightFromDiameterResults = bind_rows(as_row(psmeHeightFromDiameterChapmanRichards),
                                          as_row(psmeHeightFromDiameterChapmanRichardsBal),
                                          as_row(psmeHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row(psmeHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(psmeHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(psmeHeightFromDiameterCurtis),
                                          as_row(psmeHeightFromDiameterHossfeld),
                                          as_row(psmeHeightFromDiameterKorf),
                                          as_row(psmeHeightFromDiameterLinear),
                                          as_row(psmeHeightFromDiameterMichaelisMenten),
                                          as_row(psmeHeightFromDiameterParabolic),
                                          as_row(psmeHeightFromDiameterPower),
                                          as_row(psmeHeightFromDiameterProdan),
                                          as_row(psmeHeightFromDiameterRatkowsky),
                                          as_row(psmeHeightFromDiameterRichards),
                                          as_row(psmeHeightFromDiameterSharmaParton),
                                          as_row(psmeHeightFromDiameterSharmaPartonBal),
                                          as_row(psmeHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(psmeHeightFromDiameterSharmaPartonPhysio),
                                          as_row(psmeHeightFromDiameterSharmaZhang),
                                          as_row(psmeHeightFromDiameterSharmaZhangBal),
                                          as_row(psmeHeightFromDiameterSibbesen),
                                          as_row(psmeHeightFromDiameterWeibull),
                                          as_row(psmeHeightFromDiameterWeibullBal),
                                          as_row(psmeHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(psmeHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterKorf$fitted.values, color = "Korf", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterLinear$fitted.values, color = "linear", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterPower$fitted.values, color = "power", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterProdan$fitted.values, color = "Prodan", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = psme2016$isPlantation)) +
  annotate("text", x = 0, y = 85, label = "Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

ggplot() +
  geom_point(aes(x = psme2016physio$elevation, y = psmeHeightFromDiameterChapmanRichardsBalPhysio$residuals), alpha = 0.15, color = "grey25", shape = 16) +
  geom_smooth(aes(x = psme2016physio$elevation, y = psmeHeightFromDiameterChapmanRichardsBalPhysio$residuals), alpha = 0.15, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  labs(x = "physiographic", y = "DBH error, cm")

#psmeHeightFromDiameterChapmanRichardsBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * Elev_Mean + a5 * SlopeMean + a6 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 87, a2 = 0.49, a3 = -0.032, a4 = -0.0025, a5 = -0.16, a6 = 3.63, b1 = -0.011, b2 = 1.09), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
#psmeHeightFromDiameterChapmanRichardsNatural = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016natural, start = list(a1 = 66, b1 = -0.02, b2 = 1.5), weights = pmin(1/psme2016natural$DBH, 1))
#psmeHeightFromDiameterChapmanRichardsPhysiographic = nls(TotalHt ~ 1.37 + (a1 + a2 * Elev_Mean + a3 * SlopeMean + a4 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = -0.002, a3 = -0.11, a4 = 1.78, b1 = -0.014, b2 = 1.13), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
#psmeHeightFromDiameterChapmanRichardsPlantation = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016plantation, start = list(a1 = 51, b1 = -0.02, b2 = 1.2), weights = pmin(1/psme2016plantation$DBH, 1))
#psmeHeightFromDiameterChapmanRichardsTopHeight = nls(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 28.4, a2 = 0.21, b1 = -0.017, b2 = 1.08), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
#psmeHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * basalAreaLarger) *(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 43, a2 = 0.13, a3 = 0.002, b1 = -0.017, b2 = -0.11, b3 = 1.00), weights = pmin(DBH^if_else(isPlantation, -1.2, -0.4), 1))
#psmeHeightFromDiameterWeibullBalNatural = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016natural, start = list(a1 = 62.2, a2 = 0.038, a3 = 0.025, b1 = -0.0045, b2 = 1.32), weights = pmin(1/psme2016natural$DBH, 1))
#psmeHeightFromDiameterWeibullBalPlantation = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016plantation, start = list(a1 = 73.2, a2 = 0.93, a3 = -0.19, b1 = -0.008, b2 = 1.07), weights = pmin(1/psme2016plantation$DBH, 1))

#psmeHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(psmeHeightFromDiameterChapmanRichardsBalPhysio, psme2016, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterChapmanRichardsPhysiographic = get_height_error(psmeHeightFromDiameterChapmanRichardsPhysiographic, psme2016, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterChapmanRichardsNatural = get_height_error(psmeHeightFromDiameterChapmanRichardsNatural, psme2016, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterChapmanRichardsPlantation = get_height_error(psmeHeightFromDiameterChapmanRichardsPlantation, psme2016, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterChapmanRichardsTopHeight = get_height_error(psmeHeightFromDiameterChapmanRichardsTopHeight, psme2016, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterWeibullBalNatural = get_height_error(psmeHeightFromDiameterWeibullBalNatural, psme2016natural, psme2016natural, psme2016plantation)
#psmeHeightFromDiameterWeibullBalPlantation = get_height_error(psmeHeightFromDiameterWeibullBalPlantation, psme2016plantation, psme2016natural, psme2016plantation)

## Douglas-fir height-diameter GNLS regressions
#psmeHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-5, tolerance = 1E-4, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = psmeHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 0.001, tolerance = 0.01, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05, iterations at default msTol and tolerance
#psmeHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameterSharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 500, nlsMaxIter = 50, msTol = 1E-6, tolerance = 1E-5, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameterSharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#psmeHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameterSharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.2, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, msTol = 1E-5, tolerance = 1E-4, returnObject = FALSE)) # step having at nlsTol = 0.1
#psmeHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = psmeHeightFromDiameterSharmaZhangBal$m$getPars(), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.08, maxIter = 250, nlsMaxIter = 50, msTol = 1E-7, tolerance = 1E-6, msVerbose = FALSE, returnObject = FALSE)) # step halving factor at nlsTol = 1 with plot correlation, step halving with nlsTol = 0.07
#psmeHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.05
#psmeHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = psmeHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 5, maxIter = 250, nlsMaxIter = 50, msTol = 1E-5, tolerance = 1E-4, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.2
#save(psmeHeightFromDiameterChapmanRichardsGnls, psmeHeightFromDiameterChapmanRichardsBalGnls, psmeHeightFromDiameterSharmaPartonGnls, psmeHeightFromDiameterSharmaPartonBalGnls, psmeHeightFromDiameterSharmaZhangGnls, psmeHeightFromDiameterSharmaZhangBalGnls, psmeHeightFromDiameterWeibullGnls, psmeHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia PSME GNLS.rdata")

load("trees/height-diameter/HtDia PSME GNLS.rdata")
psmeHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", psmeHeightFromDiameterChapmanRichardsGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BAL GNLS", psmeHeightFromDiameterChapmanRichardsBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", psmeHeightFromDiameterSharmaPartonGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BAL GNLS", psmeHeightFromDiameterSharmaPartonBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", psmeHeightFromDiameterSharmaZhangGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BAL GNLS", psmeHeightFromDiameterSharmaZhangBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", psmeHeightFromDiameterWeibullGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BAL GNLS", psmeHeightFromDiameterWeibullBalGnls, psme2016, psme2016natural, psme2016plantation)

psmeHeightFromDiameterResultsGnls = bind_rows(as_row(psmeHeightFromDiameterChapmanRichardsGnls),
                                              as_row(psmeHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(psmeHeightFromDiameterSharmaPartonGnls),
                                              as_row(psmeHeightFromDiameterSharmaPartonBalGnls),
                                              as_row(psmeHeightFromDiameterSharmaZhangGnls),
                                              as_row(psmeHeightFromDiameterSharmaZhangBalGnls),
                                              as_row(psmeHeightFromDiameterWeibullGnls),
                                              as_row(psmeHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
psmeHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(psmeHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(psmeHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(psmeHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
#  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
#  select(-bal, -balN, -balP)

ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterChapmanRichardsBal$fitted.values, color = "Chapman-Richards BAL", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016physio$DBH, y = psmeHeightFromDiameterChapmanRichardsPhysio$fitted.values, color = "Chapman-Richards physio", group = psme2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = psme2016physio$DBH, y = psmeHeightFromDiameterChapmanRichardsTopo$fitted.values, color = "Chapman-Richards topo", group = psme2016physio$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaPartonBal$fitted.values, color = "Sharma-Parton BAL", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaZhangBal$fitted.values, color = "Sharma-Zhang BAL", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterWeibullBal$fitted.values, color = "Weibull BAL", group = psme2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSharmaPartonGnls$fitted.values, color = "Sharma-Parton GNLS"), alpha = 0.35) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterWeibullGnls$fitted.values, color = "Weibull GNLS"), alpha = 0.35) + 
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = psme2016$isPlantation)) +
  annotate("text", x = 0, y = 85, label = "a) Douglas-fir, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c("ElliottChapmanRichards", "ElliottWeibull", "ElliottChapmanRichardsBal", "ElliottWeibullBalNatural", "ElliottWeibullBalPlantation", "TemesgenWeibull"), labels = c("Chapman-Richards", "Weibull", "Chapman-Richards with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_color_manual(values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

ggplot(psme2016) +
  geom_point(aes(x = Elev_Mean, y = psmeHeightFromDiameterWeibullBal$residuals), alpha = 0.15, na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = Elev_Mean, y = psmeHeightFromDiameterWeibullBal$residuals), color = "red", formula = y ~ s(x, k = 10), method = "gam")

ggplot(psme2016) +
  geom_point(aes(x = AspectCos, y = psmeHeightFromDiameterWeibullBal$residuals), alpha = 0.15, na.rm = TRUE, shape = 16) +
  geom_smooth(aes(x = AspectCos, y = psmeHeightFromDiameterWeibullBal$residuals), color = "red", formula = y ~ s(x, k = 10), method = "gam")

## Douglas-fir diameter-height regressions
# Chapman-Richards, Gomperz, Prodan, Ratkowsky, and Korf run opposite the desired curvature
# Huang has DBH as ln of height (https://doi.org/10.1093/forestry/cpz002, Eq. 1)
# Hossfeld is sensitive to going negative, extrapolating poorly
# Prodan sometimes produces negative DBH
# Sharma-Parton and Sharma-Zhang have wrong curvature, simply modified forms fail to start or fail to step
# Sibbesen and Korf are near identical
#psmeDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999999)), psme2016, start = list(a1 = -65, b1 = 1/85, b2 = 0.75), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 109, a1p = -40.0, a2 = -1.42, a3 = 0.490, a3p = 0.120, a4 = -45.1, a4p = 36.5, b1 = 0.041, b2 = 0.763, b2p = -0.0092), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a2p not significant
#psmeDiameterFromHeightChapmanRichardsPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016physio, start = list(a1 = -20.8, a1p = 5.59, a2 = 0.00218, a3 = -3.39, a4 = 0.209, a5 = 0.109, a6 = 0.101, b1 = 0.0176, b1p = 0.00776, b2 = 0.329), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a5 not significant
#psmeDiameterFromHeightMichaelisMentenForm = gsl_nls(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), psme2016, start = list(a1 = 108, a2 = 46.3, b1 = 0.77), weights = TotalHt^-2)
#psmeDiameterFromHeightNaslund = gsl_nls(DBH ~ a1 * sqrt(TotalHt - 1.37) / (1 + a2 * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 3.7, a2 = -0.09), weights = TotalHt^-2) # diverges if a2 start value is less than ~-0.1
#psmeDiameterFromHeightPower = nls(DBH ~ a1*(TotalHt - 1.37)^b1, psme2016, start = list(a1 = 1.26, b1 = 1.08), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightProdan = nls(DBH ~ (TotalHt - 1.37)^2 / (a0 + a1 * (TotalHt - 1.37) + a2 * (TotalHt - 1.37)^2), psme2016, start = list(a0 = -0.7, a1 = 0.8, a2 = -0.004), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # AIC 
#psmeDiameterFromHeightProdan = get_dbh_error(psmeDiameterFromHeightProdan, psme2016, psme2016natural, psme2016plantation)
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ a1*topHeight^a2*exp(b1*(tph/standBasalAreaPerHectare)^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 1, a2 = 1, b1 = -0.01, b2 = 1, b3 = 1), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 57.1, a1p = -25.5, a2 = -0.241, a2p = 0.060, b1 = 0.0261, b2 = -0.0535, b2p = 0.0906, b3 = 0.689, b3p = 0.0829), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ a1*topHeight^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, psme2016, start = list(a1 = 39.3, a2 = 0.127, b1 = 0.0180, b2 = -0.010, b3 = 0.824), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightSharmaZhang = nls(DBH ~ a1*standBasalAreaPerHectare^a2*exp(b1*tph^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.0005, b2 = 0.5, b3 = 1), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ (a1 + a1p * isPlantation)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.081, a1p = -0.390, b1 = 1.576, b2 = 0.263, b2p = 0.0237), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 15, b1 = 0.1, b2 = 0.35), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # NaN or inf 
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 36.6, b1 = 0.059, b2 = 0.77), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 63.0, b1 = 0.018, b2 = 0.87), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ (a1 + a1p * isPlantation)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 17.8, a1p = -6.45, b1 = 0.216, b2 = 0.541, b2p = 0.048), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWykoffBal = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 121.7, a2 = -2.17, a3 = 1.04, b1 = 0.017, b2 = 0.85), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), control = list(maxiter = 500))
#psmeDiameterFromHeightWykoffBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.340, a1p = -0.175, a2 = -0.0046, a2p = 0.00162, a3 = 0.00142, a3p = -0.00058, a4 = -0.132, a4p = 0.119, b1 = 3.040, b2 = 0.173, b2p = 0.0129), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # >500 iterations with b1p
#psmeDiameterFromHeightWykoffRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.065, a1p = -0.386, a2 = 0.0071, b1 = 1.59, b2 = 0.262, b2p = 0.0234), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWeibull = gsl_nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -204, b1 = 0.011, b2 = 0.87), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightWeibull = gsl_nls(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -374, a1p = 163, b1 = 0.009, b1p = 0.003, b2 = 0.82), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # b2p not significant
#cor(cbind(dbh = psme2016$DBH, bal = psme2016$basalAreaLarger, height = psme2016$TotalHt, relHt = psme2016$relativeHeight, topHt = psme2016$topHeight, deltaHt = psme2016$TotalHt - psme2016$topHeight, treesPerPlot = psme2016$treesPerPlot, tph = psme2016$tph, tallerTph = psme2016$tallerTph * psme2016$topHeight, tallerQuasiBA = psme2016$tallerQuasiBA))
psmeDiameterFromHeightChapmanForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightChapmanFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.4, a1p = -44.6, a2 = 0.0058, a3 = -0.0544, b1 = 0.00166, b1p = 0.017, b2 = 0.788, b2p = -0.056), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a2, a2p, a3p not significant
psmeDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 135, a1p = -37.5, a2 = -1.2, b1 = 0.010, b1p = 0.002, b2 = 0.756, b2p = 0.064), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), control = nls.control(maxiter = 500)) # a2p not significant, nlrob() step factor with a3 * BA
psmeDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * relativeHeight) * (exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 139, a1p = -38.5, a2 = -1.2, a3 = -5.5, a3p = 0.22, b1 = 0.012, b1p = 0.003, b2 = 0.796, b2p = 0.066), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), control = nls.control(maxiter = 500)) # a2p not significant, nlrob() step factor with a3 * BA
psmeDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 21.1, a1p = -6.56, a2 = 0.278, b1 = 0.170, b2 = 0.580, b2p = 0.0395), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a2p not significant
psmeDiameterFromHeightChapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # b2p not significant
psmeDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a2p, a3p, b2p not significant
psmeDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.9999)), psme2016physio, start = list(a1 = -17.9, a1p = -0133, a2 = 0.0024, a3 = -7.189, a4 = 0.359, a5 = 0.170, a6 = 0.180, b1 = 0.0176, b1p = 0.0058, b2 = 0.404), maxit = 20, weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), control = nls.control(maxiter = 50)) # a5 not significant
psmeDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), psme2016, start = list(a1 = -1.22, a2 = -3.11, b1 = 0.018, b1p = 0.006, b2 = 0.003, b2p = 0.052), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), maxit = 20) # a1p not significant
psmeDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016, weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation) / (a2 + a2p * isPlantation - (TotalHt - 1.37)^(b1 + b1p * isPlantation)), psme2016, start = list(a1 = 190, a1p = -118, a2 = 67.3, a2p = -38.3, b1 = 0.78, b1p = -0.08), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), psme2016, start = list(a1 = 5.0, a1p = -1.6, a2 = -0.085, a2p = -0.018), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016, weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a3 not significant
psmeDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016physio, start = list(a1 = 1.630, a1p = 0.284, a2 = 0.00001, a3 = -0.082, a4 = -0.019, a5 = 0.0028, a6 = -0.0018, b1 = 1.03, b1p = -0.102), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a5 not significant
psmeDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a2 = 0.361, b1 = 0.943, b1p = -0.068), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a1p, a2p not significant
psmeDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
#psmeDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), psme2016, start = list(a1 = 0.002, a2 = 0.055, b1 = 1.05, Ha = 18.6), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # converges from red alder values but fails to reconverge (singular gradient or NaN-inf with nls())
psmeDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 3.95, a2 = 0.681, a2p = -0.139, b1 = 0.097, b2 = -0.130, b2p = 0.157, b3 = 0.125, b3p = 0.0678), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5), control = list(maxiter = 200)) # a1p NaN-inf (not significant?), singular gradient with all relative height forms attempted
psmeDiameterFromHeightSibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5))
psmeDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + (a3 + a3p * isPlantation) * standQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, a3p = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a2, a2p not significant
psmeDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016physio, start = list(a1 = 3.812, a1p = -0.925, a2 = 0.00022, a3 = 0.0663, a4 = -0.0495, a5 = -0.0151, a6 = -0.00630, b1 = 0.520, b2 = 0.111, b2p = 0.0173), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # a3, a5 not significant
psmeDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.90, a1p = -0.95, a2 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # neither a2 nor a2p significant
psmeDiameterFromHeightWeibull = nlrob(DBH ~ ((a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37), 0.9999)))^b2, psme2016, start = list(a1 = -347, a1p = 128, b1 = 0.010, b1p = 0.0027, b2 = 0.81), weights = pmin(TotalHt^if_else(isPlantation, -1.6, -1.7), 0.5)) # b2p not significant
#confint2(psmeDiameterFromHeightWeibull, level = 0.99)

psmeDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", psmeDiameterFromHeightChapmanForm, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AAT", psmeDiameterFromHeightChapmanFormAat, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BAL", psmeDiameterFromHeightChapmanFormBal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BAL RelHt", psmeDiameterFromHeightChapmanFormBalRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", psmeDiameterFromHeightChapmanFormRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", psmeDiameterFromHeightChapmanRichards, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AAT", psmeDiameterFromHeightChapmanRichardsAat, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", psmeDiameterFromHeightChapmanRichardsPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", psmeDiameterFromHeightChapmanRichardsRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightLinear = get_dbh_error("linear", psmeDiameterFromHeightLinear, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", psmeDiameterFromHeightMichaelisMentenForm, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightNaslund = get_dbh_error("Näslund", psmeDiameterFromHeightNaslund, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightParabolic = get_dbh_error("parabolic", psmeDiameterFromHeightParabolic, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPower = get_dbh_error("power", psmeDiameterFromHeightPower, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPowerAat = get_dbh_error("power AAT", psmeDiameterFromHeightPowerAat, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPowerPhysio = get_dbh_error("power physio", psmeDiameterFromHeightPowerPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", psmeDiameterFromHeightPowerRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightRuark = get_dbh_error("Ruark", psmeDiameterFromHeightRuark, psme2016, psme2016natural, psme2016plantation)
#psmeDiameterFromHeightSchnute = get_dbh_error("Schnute", psmeDiameterFromHeightSchnute, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", psmeDiameterFromHeightSharmaParton, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", psmeDiameterFromHeightSibbesenForm, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AAT", psmeDiameterFromHeightSibbesenFormAat, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", psmeDiameterFromHeightSibbesenFormPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", psmeDiameterFromHeightSibbesenFormRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightWeibull = get_dbh_error("Weibull", psmeDiameterFromHeightWeibull, psme2016, psme2016natural, psme2016plantation)

psmeDiameterFromHeightResults = bind_rows(as_row(psmeDiameterFromHeightChapmanRichards),
                                          as_row(psmeDiameterFromHeightChapmanRichardsAat),
                                          as_row(psmeDiameterFromHeightChapmanRichardsPhysio),
                                          as_row(psmeDiameterFromHeightChapmanRichardsRelHt),
                                          as_row(psmeDiameterFromHeightChapmanForm),
                                          as_row(psmeDiameterFromHeightChapmanFormAat),
                                          as_row(psmeDiameterFromHeightChapmanFormBal),
                                          as_row(psmeDiameterFromHeightChapmanFormBalRelHt),
                                          as_row(psmeDiameterFromHeightChapmanFormRelHt),
                                          as_row(psmeDiameterFromHeightLinear),
                                          as_row(psmeDiameterFromHeightMichaelisMentenForm),
                                          as_row(psmeDiameterFromHeightNaslund),
                                          as_row(psmeDiameterFromHeightParabolic),
                                          as_row(psmeDiameterFromHeightPower),
                                          as_row(psmeDiameterFromHeightPowerAat),
                                          as_row(psmeDiameterFromHeightPowerPhysio),
                                          as_row(psmeDiameterFromHeightPowerRelHt),
                                          as_row(psmeDiameterFromHeightRuark),
                                          as_row(name = "Schnute"),
                                          as_row(psmeDiameterFromHeightSharmaParton),
                                          as_row(psmeDiameterFromHeightSibbesenForm),
                                          as_row(psmeDiameterFromHeightSibbesenFormAat),
                                          as_row(psmeDiameterFromHeightSibbesenFormPhysio),
                                          as_row(psmeDiameterFromHeightSibbesenFormRelHt),
                                          as_row(psmeDiameterFromHeightWeibull)) %>% 
  mutate(responseVariable = "DBH", species = "PSME", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(psmeDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(psme2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman form", group = isPlantation)) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman form AAT", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = psmeDiameterFromHeightLinear$fitted.values, y = TotalHt, color = "linear", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = psmeDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightParabolic$fitted.values, y = TotalHt, color = "parabolic", group = isPlantation)) +
  #geom_line(aes(x = psmeDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightSchnute$fitted.values, y = TotalHt, color = "Schnute", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  #geom_line(aes(x = psmeDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = psmeDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = -65 * log(1 - pmin((1/85*(TotalHt - 1.37))^0.7, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = -155 * log(1 - pmin((0.00839*(TotalHt - 1.37))^0.970, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = 0.05*topHeight*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman form BAL", group = isPlantation), alpha = 0.5) +
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
  geom_point(aes(x = psme2016$tallerTph, y = psmeDiameterFromHeightChapmanRH$residuals), alpha = 0.05, shape = 16) +
  geom_smooth(aes(x = psme2016$tallerTph, y = psmeDiameterFromHeightChapmanRH$residuals), color = "red", fill = "red") +
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


## Q-Q plots
#ggplot() + # symmetric t fits less well than skewed t, as expected
#  geom_qq_line(aes(sample = psmeDiameterFromHeightChapmanRichards$residuals, color = "Chapman-Richards"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightChapmanForm$residuals, color = "Chapman-Richards form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightRuark$residuals, color = "Ruark"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightSibbesenForm$residuals, color = "Sibbesen form"), alpha = 0.4, distribution = qt, dparams = list(df = 7)) +
#  geom_qq(aes(sample = psmeDiameterFromHeightChapmanRichards$residuals, color = "Chapman-Richards"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightChapmanForm$residuals, color = "Chapman-Richards form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightRuark$residuals, color = "Ruark"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightSibbesenForm$residuals, color = "Sibbesen form"), alpha = 0.8, distribution = qt, dparams = list(df = 7), geom = "line") +
#  annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'t(df = 7, '*alpha*' = 2.25)'", hjust = 0, parse = TRUE, size = 3.5) +
#  coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
#  labs(x = "theoretical quantile", y = NULL, color = NULL) +
#  scale_color_manual(values = dbhColors) +
#  theme(legend.justification = c(1, 0), legend.position = "none")
#ggplot() + # slow! and fits poorly
#  geom_qq_line(aes(sample = psmeDiameterFromHeightChapmanRichards$residuals, color = "Chapman-Richards"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightChapmanForm$residuals, color = "Chapman-Richards form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightRuark$residuals, color = "Ruark"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq_line(aes(sample = psmeDiameterFromHeightSibbesenForm$residuals, color = "Sibbesen form"), alpha = 0.4, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34)) +
#  geom_qq(aes(sample = psmeDiameterFromHeightChapmanRichards$residuals, color = "Chapman-Richards"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightChapmanForm$residuals, color = "Chapman-Richards form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightRuark$residuals, color = "Ruark"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  geom_qq(aes(sample = psmeDiameterFromHeightSibbesenForm$residuals, color = "Sibbesen form"), alpha = 0.8, distribution = PearsonDS::qpearsonIV, dparams = list(m = 1.72, nu = 0.656, scale = 1, location = -0.34), geom = "line") +
#  annotate("text", x = -10, y = 160, label = "'d) Douglas-fir DBH, '*epsilon~'~'~'pearsonIV(m = 1.72, '*nu*' = 0.656)'", hjust = 0, parse = TRUE, size = 3.5) +
#  coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
#  labs(x = "theoretical quantile", y = NULL, color = NULL) +
#  scale_color_manual(values = dbhColors) +
#  theme(legend.justification = c(1, 0), legend.position = "none")


## collect model parameters
psmeParameters = bind_rows(bind_rows(get_coefficients(psmeHeightFromDiameterChapmanRichards),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(psmeHeightFromDiameterCurtis),
                                     get_coefficients(psmeHeightFromDiameterHossfeld),
                                     get_coefficients(psmeHeightFromDiameterKorf),
                                     get_coefficients(psmeHeightFromDiameterLinear),
                                     get_coefficients(psmeHeightFromDiameterMichaelisMenten),
                                     get_coefficients(psmeHeightFromDiameterParabolic),
                                     get_coefficients(psmeHeightFromDiameterPower),
                                     get_coefficients(psmeHeightFromDiameterProdan),
                                     get_coefficients(psmeHeightFromDiameterRatkowsky),
                                     get_coefficients(psmeHeightFromDiameterRichards),
                                     get_coefficients(psmeHeightFromDiameterSharmaParton),
                                     get_coefficients(psmeHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(psmeHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(psmeHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(psmeHeightFromDiameterSharmaZhang),
                                     get_coefficients(psmeHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(psmeHeightFromDiameterSibbesen),
                                     get_coefficients(psmeHeightFromDiameterWeibull),
                                     get_coefficients(psmeHeightFromDiameterWeibullBal),
                                     get_coefficients(psmeHeightFromDiameterWeibullBalRelHt),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(psmeHeightFromDiameterChapmanRichardsBalGnls),
                                     get_coefficients(psmeHeightFromDiameterSharmaPartonGnls),
                                     get_coefficients(psmeHeightFromDiameterSharmaPartonBalGnls),
                                     get_coefficients(psmeHeightFromDiameterSharmaZhangGnls),
                                     get_coefficients(psmeHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(psmeHeightFromDiameterWeibullGnls),
                                     get_coefficients(psmeHeightFromDiameterWeibullBalGnls)) %>%
                              mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(psmeDiameterFromHeightChapmanRichards),
                                     get_coefficients(psmeDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(psmeDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(psmeDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(psmeDiameterFromHeightChapmanForm),
                                     get_coefficients(psmeDiameterFromHeightChapmanFormAat),
                                     get_coefficients(psmeDiameterFromHeightChapmanFormBal),
                                     get_coefficients(psmeDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(psmeDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(psmeDiameterFromHeightLinear),
                                     get_coefficients(psmeDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(psmeDiameterFromHeightNaslund),
                                     get_coefficients(psmeDiameterFromHeightParabolic),
                                     get_coefficients(psmeDiameterFromHeightPower),
                                     get_coefficients(psmeDiameterFromHeightPowerAat),
                                     get_coefficients(psmeDiameterFromHeightPowerPhysio),
                                     get_coefficients(psmeDiameterFromHeightPowerRelHt),
                                     get_coefficients(psmeDiameterFromHeightRuark),
                                     #get_coefficients(psmeDiameterFromHeightSchnute),
                                     get_coefficients(psmeDiameterFromHeightSharmaParton),
                                     get_coefficients(psmeDiameterFromHeightSibbesenForm),
                                     get_coefficients(psmeDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(psmeDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(psmeDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(psmeDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "PSME", 
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
# essentially no difference between gsl_nls() and nlrob() fits
# Chapman-Richards has the wrong curvature
#psmeBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1, psme2016, start = list(a1 = 1, a1p = 0, b1 = 0.00009, b1p = 0, b2 = 2.2, b2p = 0), weights = pmin(1/basalArea, 1E4)) # a1p not significant
#psmeBasalAreaFromHeightPower = gsl_nls(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 0.25 * pi * 0.01^2, a1p = 0, b1 = 2.4, b1p = 0), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
psmeBasalAreaFromHeightKorf = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(exp((b1 + b1p * isPlantation)*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 0.689, a1p = -0.413, b1 = 0.0003, b1p = 0.0005, b2 = 1.91, b2p = -0.10), weights = pmin(1/basalArea, 1E4)) # a1p not significant
psmeBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 4/7 * 0.25 * pi * 0.01^2, a1p = 0.00005, b1 = 2.41, b1p = -0.248), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²
#confint_nlrob(psmeBasalAreaFromHeightKorf, level = 0.99, weights = 1/psme2016$basalArea)

psmeBasalAreaFromHeightKorf$fitted.values = predict(psmeBasalAreaFromHeightKorf, psme2016)
psmeBasalAreaFromHeightKorf$residuals = psmeBasalAreaFromHeightKorf$fitted.values - psme2016$basalArea
psmeBasalAreaFromHeightPower$fitted.values = predict(psmeBasalAreaFromHeightPower, psme2016)
psmeBasalAreaFromHeightPower$residuals = psmeBasalAreaFromHeightPower$fitted.values - psme2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(psmeBasalAreaFromHeightKorf), 100^2 * mean(psmeBasalAreaFromHeightKorf$residuals), mean(abs(psmeBasalAreaFromHeightKorf$residuals)), 1 - sum(psmeBasalAreaFromHeightKorf$residuals^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2)),
        "power", AIC(psmeBasalAreaFromHeightPower), 100^2 * mean(psmeBasalAreaFromHeightPower$residuals), mean(abs(psmeBasalAreaFromHeightPower$residuals)), 1 - sum(psmeBasalAreaFromHeightPower$residuals^2) / sum((psme2016$basalArea - mean(psme2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(psme2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = psmeBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "measured or imputed Douglas-fir height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))

ggplot(psme2016) +
  geom_point(aes(x = basalArea, y = psmeBasalAreaFromHeightKorf$residuals, color = "Korf"), alpha = 0.1, shape = 16) +
  labs(x = "basal area, m²", y = "residual, m", color = NULL) +
  theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
ggplot(psme2016) +
  geom_point(aes(x = basalArea, y = psmeBasalAreaFromHeightPower$residuals, color = "power"), alpha = 0.1, shape = 16) +
  labs(x = "basal area, m²", y = "residual, m", color = NULL) +
  theme(legend.justification = c(0, 0), legend.position = c(0.02, 0.02)) +
plot_layout(nrow = 1, ncol = 2)
  
ggplot(psme2016) +
  geom_point(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.10, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.95)) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.25, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.80)) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.35, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.50)) +
  geom_quantile(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), color = "darkviolet", formula = y ~ x) +
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
