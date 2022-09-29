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
#psmeHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + a1*DBH / (a2 + DBH), psme2016, start = list(a1 = 138, a2 = 156), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterProdan = nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + a2*DBH + a3), psme2016, start = list(a1 = 0.008, a2 = 1.0, a3 = 2.7), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterPower = nls(TotalHt ~ 1.37 + b0*DBH^b1, psme2016, start = list(b0 = 1.54, b1 = 0.77), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRatkowsky = nls(TotalHt ~ 1.37 + a1*exp(b1/(DBH + b2)), psme2016, start = list(a1 = 93, b1 = -65, b2 = 15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 58.3, d = 0.609, kU = 0.0134), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), psme2016, start = list(Ha = 62.6, Hap = -26.5, d = 0.609, kU = 0.0128, kUp = 0.0114), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaParton = nls(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 31, a2 = 0.19, b1 = -0.020, b2 = -0.13, b3 = 0.99), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 - exp(b1*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, psme2016, start = list(a1 = 39, a2 = 0.15, b1 = -0.018, b2 = -0.16, b3 = 1.01), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 55, a2 = 0.07, b1 = -0.012, b2 = 0.03, b3 = 1.1), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaZhangBal = nls(TotalHt ~ 1.37 + a1*standBasalAreaPerHectare^a2 * (1 + a3 * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, psme2016, start = list(a1 = 66, a2 = 0.056, a3 = 0.006, b1 = -0.022, b2 = -0.13, b3 = 1.04), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), psme2016, start = list(a1 = 0.034, b1 = 3.1, b2 = -0.15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterKorf = nls(TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), psme2016, start = list(a1 = 172, b1 = -9.682, b2 = -0.4574), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWykoff = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 69, b1 = -0.0075, b2 = 1.15), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWykoffBal = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), psme2016, start = list(a1 = 73.7, a2 = 0.38, a3 = -0.007, b1 = -0.008, b2 = 1.09), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterWykoffBalRelHt = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1 * (1 + b2 * pmin(relativeHeight, 1.25)) * DBH^b3)), psme2016, start = list(a1 = 220, a2 = 3.38, a3 = -0.119, a4 = -15.3, b1 = -0.0024, b2 = 2.89, b3 = 0.737), weights = pmin(DBH^-2, 1))
# fits with isPlantation
psmeHeightFromDiameterChapmanRichards = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), weights = pmin(DBH^-2, 1)) # b1p not significant
psmeHeightFromDiameterChapmanRichardsBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterChapmanRichardsBalPhysio = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(pi/180 * aspect) + a6 * cos(pi/180 * aspect)) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 69.4, a2 = 0.077, a2p = 0.72, a3 = -0.008, a4 = -0.031, a5 = 0.548, a6 = 0.980, b1 = -0.021, b1p = 0.0075, b2 = 1.49, b2p = -0.370), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterChapmanRichardsBalRelHt = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 8.8, a1p = 11.0, a2 = 0.18, a2p = 0.42, a3 = -0.0083, a3p = 0.070, a4 = 54.0, a4p = -28.3, b1 = -0.021, b2 = 0.65, b2p = 0.37), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterChapmanRichardsPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 68.5, a1p = -13.4, a2 = -0.0045, a3 = -8.09, a4 = 0.783, a5 = 0.766, a6 = 0.213, b1 = -0.022, b2 = 1.50, b2p = -0.31), weights = pmin(DBH^-2, 1)) # a4p not significant, a5p induces overfitting
psmeHeightFromDiameterCurtis = nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation) * DBH / (1 + DBH)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.409, a1p = -0.685, b1 = 0.200, b1p = -0.156), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterHossfeld = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) * DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 75.4, a1p = -11.4, b1 = 462, b1p = -322, b2 = -1.54, b2p = 0.28), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterKorf = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 320, a1p = 376, b1 = -7.83, b2 = -0.323, b2p = 0.084), weights = pmin(DBH^-2, 1)) # b1p not significant
psmeHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), psme2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 87.8, a1p = -26.9, a2 = 236, a2p = -92.0, b1 = 1.30), weights = pmin(DBH^-2, 1)) # b1p not significant
psmeHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), psme2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterProdan = nls(TotalHt ~ 1.37 + DBH^2 / (a1*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3 + a3p* isPlantation), psme2016, start = list(a1 = 0.012, a2 = 0.41, a2p = 0.47, a3 = 17.9, a3p = -14.6), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterPower = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 2.9, a1p = -1.6, b1 = 0.63, b1p = 0.18), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterRatkowsky = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), psme2016, start = list(a1 = 90.0, a1p = -25.8, b1 = -55.2, b1p = 14.5, b2 = 10.0, b2p = -1.52), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - (d + dp*isPlantation)))))^(1/(1 - (d + dp*isPlantation))), psme2016, start = list(Ha = 65.3, Hap = -29.3, d = 0.574, dp = 0.151, kU = 0.0118, kUp = 0.0126), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaParton = nls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaPartonBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(pi/180 * aspect) + a5 * cos(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaPartonPhysio = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(pi/180 * aspect) + a5 * cos(pi/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaZhang = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSharmaZhangBal = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterSibbesen = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.122, a1p = 0.157, b1 = 2.38, b1p = -0.552, b2 = -0.127, b2p = 0.022), weights = pmin(DBH^-2, 1), control = list(maxiter = 200))
psmeHeightFromDiameterWeibull = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterWeibullBal = nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), weights = pmin(DBH^-2, 1))
psmeHeightFromDiameterWeibullBalRelHt = nls(TotalHt ~ 1.37 + (a1  + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation) * (1 + (b2 + b2p * isPlantation) * pmin(relativeHeight, 1.25)) * DBH^b3)), psme2016, start = list(a1 = -3.82, a1p = 87.0, a2 = 0.117, a2p = 1.86, a3 = -0.00061, a3p = -0.0657, a4 = 68.1, a4p = -54.0, b1 = -0.146, b1p = 0.140, b2 = -0.610, b2p = 2.25, b3 = 0.800), weights = pmin(DBH^-2, 1))
#confint2(psmeHeightFromDiameterKorf, level = 0.99)

psme2016physio = psme2016 %>% filter(is.na(elevation) == FALSE)
psme2016plantationPhysio = psme2016physio %>% filter(isPlantation)
psmeHeightFromDiameterChapmanRichards = get_height_error(psmeHeightFromDiameterChapmanRichards, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBal = get_height_error(psmeHeightFromDiameterChapmanRichardsBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(psmeHeightFromDiameterChapmanRichardsBalPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(psmeHeightFromDiameterChapmanRichardsBalRelHt, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsPhysio = get_height_error(psmeHeightFromDiameterChapmanRichardsPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterCurtis = get_height_error(psmeHeightFromDiameterCurtis, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterHossfeld = get_height_error(psmeHeightFromDiameterHossfeld, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterKorf = get_height_error(psmeHeightFromDiameterKorf, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterLinear = get_height_error(psmeHeightFromDiameterLinear, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterMichaelisMenten = get_height_error(psmeHeightFromDiameterMichaelisMenten, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterParabolic = get_height_error(psmeHeightFromDiameterParabolic, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterProdan = get_height_error(psmeHeightFromDiameterProdan, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterPower = get_height_error(psmeHeightFromDiameterPower, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterRatkowsky = get_height_error(psmeHeightFromDiameterRatkowsky, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterRichards = get_height_error(psmeHeightFromDiameterRichards, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaParton = get_height_error(psmeHeightFromDiameterSharmaParton, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBal = get_height_error(psmeHeightFromDiameterSharmaPartonBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBalPhysio = get_height_error(psmeHeightFromDiameterSharmaPartonBalPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterSharmaPartonPhysio = get_height_error(psmeHeightFromDiameterSharmaPartonPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeHeightFromDiameterSharmaZhang = get_height_error(psmeHeightFromDiameterSharmaZhang, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangBal = get_height_error(psmeHeightFromDiameterSharmaZhangBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSibbesen = get_height_error(psmeHeightFromDiameterSibbesen, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibull = get_height_error(psmeHeightFromDiameterWeibull, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBal = get_height_error(psmeHeightFromDiameterWeibullBal, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBalRelHt = get_height_error(psmeHeightFromDiameterWeibullBalRelHt, psme2016, psme2016natural, psme2016plantation)

psmeHeightFromDiameterResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(psmeHeightFromDiameterChapmanRichards),
                                        "Chapman-Richards BAL", !!!as_row(psmeHeightFromDiameterChapmanRichardsBal),
                                        "Chapman-Richards BAL physio", !!!as_row(psmeHeightFromDiameterChapmanRichardsBalPhysio),
                                        "Chapman-Richards BAL RelHt", !!!as_row(psmeHeightFromDiameterChapmanRichardsBalRelHt),
                                        "Chapman-Richards physio", !!!as_row(psmeHeightFromDiameterChapmanRichardsPhysio),
                                        "Curtis", !!!as_row(psmeHeightFromDiameterCurtis),
                                        "Hossfeld", !!!as_row(psmeHeightFromDiameterHossfeld),
                                        "Korf", !!!as_row(psmeHeightFromDiameterKorf),
                                        "linear", !!!as_row(psmeHeightFromDiameterLinear),
                                        "Michaelis-Menten", !!!as_row(psmeHeightFromDiameterMichaelisMenten),
                                        "parabolic", !!!as_row(psmeHeightFromDiameterParabolic),
                                        "power", !!!as_row(psmeHeightFromDiameterPower),
                                        "Prodan", !!!as_row(psmeHeightFromDiameterProdan),
                                        "Ratkowsky", !!!as_row(psmeHeightFromDiameterRatkowsky),
                                        "unified Richards", !!!as_row(psmeHeightFromDiameterRichards),
                                        "Sharma-Parton", !!!as_row(psmeHeightFromDiameterSharmaParton),
                                        "Sharma-Parton BAL", !!!as_row(psmeHeightFromDiameterSharmaPartonBal),
                                        "Sharma-Parton BAL physio", !!!as_row(psmeHeightFromDiameterSharmaPartonBalPhysio),
                                        "Sharma-Parton physio", !!!as_row(psmeHeightFromDiameterSharmaPartonPhysio),
                                        "Sharma-Zhang", !!!as_row(psmeHeightFromDiameterSharmaZhang),
                                        "Sharma-Zhang BAL", !!!as_row(psmeHeightFromDiameterSharmaZhangBal),
                                        "Sibbesen", !!!as_row(psmeHeightFromDiameterSibbesen),
                                        "Weibull", !!!as_row(psmeHeightFromDiameterWeibull),
                                        "Weibull BAL", !!!as_row(psmeHeightFromDiameterWeibullBal),
                                        "Weibull BAL RelHt", !!!as_row(psmeHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "PSME", deltaAic = aic - min(aic)) %>%
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
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = psme2016$isPlantation)) +
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

#psmeHeightFromDiameterChapmanRichardsBalPhysio = nls(TotalHt ~ 1.37 + (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * Elev_Mean + a5 * SlopeMean + a6 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 87, a2 = 0.49, a3 = -0.032, a4 = -0.0025, a5 = -0.16, a6 = 3.63, b1 = -0.011, b2 = 1.09), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterChapmanRichardsNatural = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016natural, start = list(a1 = 66, b1 = -0.02, b2 = 1.5), weights = pmin(1/psme2016natural$DBH, 1))
#psmeHeightFromDiameterChapmanRichardsPhysiographic = nls(TotalHt ~ 1.37 + (a1 + a2 * Elev_Mean + a3 * SlopeMean + a4 * AspectCos) * (1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 80, a2 = -0.002, a3 = -0.11, a4 = 1.78, b1 = -0.014, b2 = 1.13), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterChapmanRichardsPlantation = nls(TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, psme2016plantation, start = list(a1 = 51, b1 = -0.02, b2 = 1.2), weights = pmin(1/psme2016plantation$DBH, 1))
#psmeHeightFromDiameterChapmanRichardsTopHeight = nls(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*DBH))^b2, psme2016, start = list(a1 = 28.4, a2 = 0.21, b1 = -0.017, b2 = 1.08), weights = pmin(DBH^-2, 1))
#psmeHeightFromDiameterSharmaPartonBal = nls(TotalHt ~ 1.37 + a1*topHeight^a2 * (1 + a3 * basalAreaLarger) *(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 43, a2 = 0.13, a3 = 0.002, b1 = -0.017, b2 = -0.11, b3 = 1.00), weights = pmin(DBH^-2, 1))
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
#psmeHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterChapmanRichardsGnlsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterSharmaPartonGnlsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterSharmaZhangGnlsBal = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#psmeHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), psme2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(psmeHeightFromDiameterChapmanRichardsGnls, psmeHeightFromDiameterChapmanRichardsBalGnls, psmeHeightFromDiameterSharmaPartonGnls, psmeHeightFromDiameterSharmaPartonBalGnls, psmeHeightFromDiameterSharmaZhangGnls, psmeHeightFromDiameterSharmaZhangBalGnls, psmeHeightFromDiameterWeibullGnls, psmeHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia PSME GNLS.rdata")
load("trees/height-diameter/HtDia PSME GNLS.rdata")
psmeHeightFromDiameterWeibullGnls = psmeHeightFromDiameterWykoffGnls # temporary naming error fixup
psmeHeightFromDiameterWeibullBalGnls = psmeHeightFromDiameterWykoffBalGnls

psmeHeightFromDiameterChapmanRichardsGnls = get_height_error(psmeHeightFromDiameterChapmanRichardsGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterChapmanRichardsBalGnls = get_height_error(psmeHeightFromDiameterChapmanRichardsBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonGnls = get_height_error(psmeHeightFromDiameterSharmaPartonGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaPartonBalGnls = get_height_error(psmeHeightFromDiameterSharmaPartonBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangGnls = get_height_error(psmeHeightFromDiameterSharmaZhangGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterSharmaZhangBalGnls = get_height_error(psmeHeightFromDiameterSharmaZhangBalGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullGnls = get_height_error(psmeHeightFromDiameterWeibullGnls, psme2016, psme2016natural, psme2016plantation)
psmeHeightFromDiameterWeibullBalGnls = get_height_error(psmeHeightFromDiameterWeibullBalGnls, psme2016, psme2016natural, psme2016plantation)

psmeHeightFromDiameterResultsGnls = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                            "Chapman-Richards GNLS", !!!as_row(psmeHeightFromDiameterChapmanRichardsGnls),
                                            "Chapman-Richards BAL GNLS", !!!as_row(psmeHeightFromDiameterChapmanRichardsBalGnls),
                                            "Sharma-Parton GNLS", !!!as_row(psmeHeightFromDiameterSharmaPartonGnls),
                                            "Sharma-Parton BAL GNLS", !!!as_row(psmeHeightFromDiameterSharmaPartonBalGnls),
                                            "Sharma-Zhang GNLS", !!!as_row(psmeHeightFromDiameterSharmaZhangGnls),
                                            "Sharma-Zhang BAL GNLS", !!!as_row(psmeHeightFromDiameterSharmaZhangBalGnls),
                                            "Weibull GNLS", !!!as_row(psmeHeightFromDiameterWeibullGnls),
                                            "Weibull BAL GNLS", !!!as_row(psmeHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "PSME", deltaAic = aic - min(aic)) %>%
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
#psmeDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999999)), psme2016, start = list(a1 = -65, b1 = 1/85, b2 = 0.75), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightPower = nls(DBH ~ a1*(TotalHt - 1.37)^b1, psme2016, start = list(a1 = 1.26, b1 = 1.08), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightProdan = nls(DBH ~ (TotalHt - 1.37)^2 / (a0 + a1 * (TotalHt - 1.37) + a2 * (TotalHt - 1.37)^2), psme2016, start = list(a0 = -0.7, a1 = 0.8, a2 = -0.004), weights = pmin(TotalHt^-2, 0.5)) # AIC 
#psmeDiameterFromHeightProdan = get_dbh_error(psmeDiameterFromHeightProdan, psme2016, psme2016natural, psme2016plantation)
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ a1*topHeight^a2*exp(b1*(tph/standBasalAreaPerHectare)^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 1, a2 = 1, b1 = -0.01, b2 = 1, b3 = 1), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 57.1, a1p = -25.5, a2 = -0.241, a2p = 0.060, b1 = 0.0261, b2 = -0.0535, b2p = 0.0906, b3 = 0.689, b3p = 0.0829), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightSharmaParton = nls(DBH ~ a1*topHeight^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, psme2016, start = list(a1 = 39.3, a2 = 0.127, b1 = 0.0180, b2 = -0.010, b3 = 0.824), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightSharmaZhang = nls(DBH ~ a1*standBasalAreaPerHectare^a2*exp(b1*tph^b2*(TotalHt - 1.37))^b3, psme2016, start = list(a1 = 5, a2 = 0.5, b1 = 0.0005, b2 = 0.5, b3 = 1), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ (a1 + a1p * isPlantation)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.081, a1p = -0.390, b1 = 1.576, b2 = 0.263, b2p = 0.0237), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 15, b1 = 0.1, b2 = 0.35), weights = pmin(TotalHt^-2, 0.5)) # NaN or inf 
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 36.6, b1 = 0.059, b2 = 0.77), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, psme2016, start = list(a1 = 63.0, b1 = 0.018, b2 = 0.87), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightWykoff = nls(DBH ~ (a1 + a1p * isPlantation)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 17.8, a1p = -6.45, b1 = 0.216, b2 = 0.541, b2p = 0.048), weights = pmin(TotalHt^-2, 0.5))
#psmeDiameterFromHeightWykoffBal = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), psme2016, start = list(a1 = 121.7, a2 = -2.17, a3 = 1.04, b1 = 0.017, b2 = 0.85), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 500))
#psmeDiameterFromHeightWykoffBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 0.340, a1p = -0.175, a2 = -0.0046, a2p = 0.00162, a3 = 0.00142, a3p = -0.00058, a4 = -0.132, a4p = 0.119, b1 = 3.040, b2 = 0.173, b2p = 0.0129), weights = pmin(TotalHt^-2, 0.5)) # >500 iterations with b1p
#psmeDiameterFromHeightWykoffRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 1.065, a1p = -0.386, a2 = 0.0071, b1 = 1.59, b2 = 0.262, b2p = 0.0234), weights = pmin(TotalHt^-2, 0.5))
#cor(cbind(dbh = psme2016$DBH, bal = psme2016$basalAreaLarger, height = psme2016$TotalHt, relHt = psme2016$relativeHeight, topHt = psme2016$topHeight, deltaHt = psme2016$TotalHt - psme2016$topHeight, treesPerPlot = psme2016$treesPerPlot, tph = psme2016$tph, tallerTph = psme2016$tallerTph * psme2016$topHeight, tallerQuasiBA = psme2016$tallerQuasiBA))
psmeDiameterFromHeightChapmanForm = nls(DBH ~ (a1 + a1p * isPlantation)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.6, a1p = -47.4, b1 = 0.016, b1p = 0.020, b2 = 0.792, b2p = -0.0780), weights = pmin(TotalHt^-2, 0.5))
psmeDiameterFromHeightChapmanFormAal = nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*(exp((b1 + b1p * isPlantation)*(TotalHt - 1.37)) - 1)^(b2 + b2p * isPlantation), psme2016, start = list(a1 = 75.4, a1p = -44.6, a2 = 0.0058, a3 = -0.0544, b1 = 0.00166, b1p = 0.017, b2 = 0.788, b2p = -0.056), weights = pmin(TotalHt^-2, 0.5)) # a2, a2p, a3p not significant
psmeDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 149.5, a1p = -74.0, a2 = -2.68, a2p = -0.91, a3 = 1.13, a3p = -0.243, b1 = 0.0363, b2 = 0.636, b2p = 0.0888), weights = pmin(TotalHt^-2, 0.5)) # step factor with b1p
psmeDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 109, a1p = -40.0, a2 = -1.42, a3 = 0.490, a3p = 0.120, a4 = -45.1, a4p = 36.5, b1 = 0.041, b2 = 0.763, b2p = -0.0092), weights = pmin(TotalHt^-2, 0.5)) # a2p not significant
psmeDiameterFromHeightChapmanFormRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), psme2016, start = list(a1 = 21.1, a1p = -6.56, a2 = 0.278, b1 = 0.170, b2 = 0.580, b2p = 0.0395), weights = pmin(TotalHt^-2, 0.5)) # a2p not significant
psmeDiameterFromHeightChapmanRichards = nls(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016, start = list(a1 = -121, a1p = 48.7, b1 = 0.00866, b1p = 0.00364, b2 = 0.781), weights = pmin(TotalHt^-2, 0.5)) # b2p not significant
psmeDiameterFromHeightChapmanRichardsAal = nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016, start = list(a1 = -136, a1p = 59.2, a2 = 0.109, a3 = 0.0684, b1 = 0.00811, b1p = 0.00403, b2 = 0.786), weights = pmin(TotalHt^-2, 0.5)) # a2p, a3p, b2p not significant
psmeDiameterFromHeightChapmanRichardsPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^b2, 0.999)), psme2016, start = list(a1 = -20.8, a1p = 5.59, a2 = 0.00218, a3 = -3.39, a4 = 0.209, a5 = 0.109, a6 = 0.101, b1 = 0.0176, b1p = 0.00776, b2 = 0.329), weights = pmin(TotalHt^-2, 0.5)) # a5 not significant
psmeDiameterFromHeightChapmanRichardsRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.999)), psme2016, start = list(a1 = -5.45, a2 = -10.4, b1 = 0.0214, b1p = 0.00209, b2 = 0.0560, b2p = 0.150), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
psmeDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), psme2016, weights = TotalHt^-2)
psmeDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), psme2016, weights = TotalHt^-2)
psmeDiameterFromHeightPower = nls(DBH ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 1.57, a1p = 0.327, b1 = 1.04, b1p = -0.108), weights = pmin(TotalHt^-2, 0.5))
psmeDiameterFromHeightPowerAal = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), psme2016, start = list(a1 = 2.14, a1p = -0.051, a2 = -0.0065, a2p = -0.0038, a3 = 0.00085, b1 = 0.963, b1p = -0.053), weights = pmin(TotalHt^-2, 0.5)) # a3 not significant
psmeDiameterFromHeightPowerPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.630, a1p = 0.284, a2 = 0.00001, a3 = -0.082, a4 = -0.019, a5 = 0.0028, a6 = -0.0018, b1 = 1.03, b1p = -0.102), weights = pmin(TotalHt^-2, 0.5)) # a5 not significant
psmeDiameterFromHeightPowerRelHt = nls(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), psme2016, start = list(a1 = 1.95, a2 = 0.361, b1 = 0.943, b1p = -0.068), weights = pmin(TotalHt^-2, 0.5)) # a1p, a2p not significant
psmeDiameterFromHeightSharmaParton = nls(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 3.95, a2 = 0.681, a2p = -0.139, b1 = 0.097, b2 = -0.130, b2p = 0.157, b3 = 0.125, b3p = 0.0678), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 200)) # a1p NaN-inf (not significant?), singular gradient with all relative height forms attempted
psmeDiameterFromHeightSibbesenForm = nls(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.89, a1p = -0.922, b1 = 0.519, b2 = 0.111, b2p = 0.017), weights = pmin(TotalHt^-2, 0.5))
psmeDiameterFromHeightSibbesenFormAal = nls(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea + (a3 + a3p * isPlantation) * standQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.898, a1p = -0.879, a2 = 0.00198, a3 = -0.00386, a3p = -0.00386, b1 = 0.527, b2 = 0.111, b2p = 0.0190), weights = pmin(TotalHt^-2, 0.5)) # a2, a2p not significant
psmeDiameterFromHeightSibbesenFormPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(pi/180 * slope) + a4 * cos(pi/180 * aspect) + a5 * sin(pi/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.812, a1p = -0.925, a2 = 0.00022, a3 = 0.0663, a4 = -0.0495, a5 = -0.0151, a6 = -0.00630, b1 = 0.520, b2 = 0.111, b2p = 0.0173), weights = pmin(TotalHt^-2, 0.5)) # a3, a5 not significant
psmeDiameterFromHeightSibbesenFormRelHt = nls(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), psme2016, start = list(a1 = 3.90, a1p = -0.95, a2 = 0.085, b1 = 0.520, b2 = 0.109, b2p = 0.016), weights = pmin(TotalHt^-2, 0.5)) # neither a2 nor a2p significant
#confint2(psmeDiameterFromHeightPower, level = 0.99)

psmeDiameterFromHeightChapmanForm = get_dbh_error(psmeDiameterFromHeightChapmanForm, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormAal = get_dbh_error(psmeDiameterFromHeightChapmanFormAal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormBal = get_dbh_error(psmeDiameterFromHeightChapmanFormBal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(psmeDiameterFromHeightChapmanFormBalRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanFormRelHt = get_dbh_error(psmeDiameterFromHeightChapmanFormRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichards = get_dbh_error(psmeDiameterFromHeightChapmanRichards, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichardsAal = get_dbh_error(psmeDiameterFromHeightChapmanRichardsAal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(psmeDiameterFromHeightChapmanRichardsPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(psmeDiameterFromHeightChapmanRichardsRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightLinear = get_dbh_error(psmeDiameterFromHeightLinear, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightParabolic = get_dbh_error(psmeDiameterFromHeightParabolic, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPower = get_dbh_error(psmeDiameterFromHeightPower, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPowerAal = get_dbh_error(psmeDiameterFromHeightPowerAal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightPowerPhysio = get_dbh_error(psmeDiameterFromHeightPowerPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightPowerRelHt = get_dbh_error(psmeDiameterFromHeightPowerRelHt, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSharmaParton = get_dbh_error(psmeDiameterFromHeightSharmaParton, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenForm = get_dbh_error(psmeDiameterFromHeightSibbesenForm, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenFormAal = get_dbh_error(psmeDiameterFromHeightSibbesenFormAal, psme2016, psme2016natural, psme2016plantation)
psmeDiameterFromHeightSibbesenFormPhysio = get_dbh_error(psmeDiameterFromHeightSibbesenFormPhysio, psme2016physio, psme2016natural, psme2016plantationPhysio)
psmeDiameterFromHeightSibbesenFormRelHt = get_dbh_error(psmeDiameterFromHeightSibbesenFormRelHt, psme2016, psme2016natural, psme2016plantation)

psmeDiameterFromHeightResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        #"Prodan", !!!as_row(psmeDiameterFromHeightProdan),
                                        "Chapman-Richards", !!!as_row(psmeDiameterFromHeightChapmanRichards),
                                        "Chapman-Richards AAL", !!!as_row(psmeDiameterFromHeightChapmanRichardsAal), # not AIC supported
                                        "Chapman-Richards physio", !!!as_row(psmeDiameterFromHeightChapmanRichardsPhysio),
                                        "Chapman-Richards RelHt", !!!as_row(psmeDiameterFromHeightChapmanRichardsRelHt), # not AIC supported
                                        "Chapman-Richards form", !!!as_row(psmeDiameterFromHeightChapmanForm),
                                        "Chapman-Richards form AAL", !!!as_row(psmeDiameterFromHeightChapmanFormAal), # not AIC supported
                                        "Chapman-Richards form BAL", !!!as_row(psmeDiameterFromHeightChapmanFormBal),
                                        "Chapman-Richards form BAL RelHt", !!!as_row(psmeDiameterFromHeightChapmanFormBalRelHt),
                                        "Chapman-Richards form RelHt", !!!as_row(psmeDiameterFromHeightChapmanFormRelHt), # not AIC supported                              "linear", !!!as_row(psmeDiameterFromHeightLinear),
                                        "linear", !!!as_row(psmeDiameterFromHeightLinear),
                                        "parabolic", !!!as_row(psmeDiameterFromHeightParabolic),
                                        "power", !!!as_row(psmeDiameterFromHeightPower),
                                        "power AAL", !!!as_row(psmeDiameterFromHeightPowerAal),
                                        "power physio", !!!as_row(psmeDiameterFromHeightPowerPhysio),
                                        "power RelHt", !!!as_row(psmeDiameterFromHeightPowerRelHt),
                                        "modified Sharma-Parton", !!!as_row(psmeDiameterFromHeightSharmaParton),
                                        "Sibbesen form", !!!as_row(psmeDiameterFromHeightSibbesenForm),
                                        "Sibbesen form AAL", !!!as_row(psmeDiameterFromHeightSibbesenFormAal), # not AIC supported
                                        "Sibbesen form physio", !!!as_row(psmeDiameterFromHeightSibbesenFormPhysio),
                                        "Sibbesen form RelHt", !!!as_row(psmeDiameterFromHeightSibbesenFormRelHt)) %>%  # not AIC supported
  mutate(responseVariable = "height", species = "PSME", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
psmeDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic)


ggplot(psme2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  geom_line(aes(x = psmeDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman form", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightLinear$fitted.values, y = TotalHt, color = "linear", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightParabolic$fitted.values, y = TotalHt, color = "parabolic", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen form", group = isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman form AAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = psmeDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -65 * log(1 - pmin((1/85*(TotalHt - 1.37))^0.7, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = -155 * log(1 - pmin((0.00839*(TotalHt - 1.37))^0.970, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = 0.05*topHeight*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman top height", group = isPlantation), alpha = 0.5) +
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

## basal area from height
# Chapman-Richards has the wrong curvature
psmeBasalAreaFromHeightKorf = nls(basalArea ~ a1*exp(b1*(imputedHeight - 1.37)^b2) - 1, psme2016, start = list(a1 = 1, b1 = 0.00009, b2 = 2.2), weights = pmin(1/basalArea, 1E4)) # 
psmeBasalAreaFromHeightPower = nls(basalArea ~ a1*(imputedHeight - 1.37)^b1, psme2016, start = list(a1 = 0.25 * pi * 0.01^2, b1 = 2.4), weights = pmin(1/basalArea, 1E4)) # 0.25 * pi * (1/height-diameter ratio)²

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
  geom_path(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$fitted.values, color = "Korf")) +
  geom_path(aes(x = imputedHeight, y = psmeBasalAreaFromHeightPower$fitted.values, color = "power")) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "measured or imputed height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))

ggplot(psme2016) +
  geom_point(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.10, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.95)) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.25, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.80)) +
  stat_summary_bin(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), alpha = 0.35, binwidth = 2.5, fill = "violet", fun.data = median_hilow, geom = "ribbon", fun.args = list(conf.int = 0.50)) +
  geom_quantile(aes(x = imputedHeight, y = psmeBasalAreaFromHeightKorf$residuals), color = "darkviolet", formula = y ~ x) +
  labs(x = "measured or imputed height, m", y = "basal area residual, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))

## collect model parameters
psmeParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", psmeHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", psmeHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", psmeHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", psmeHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", psmeHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", psmeHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", psmeHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", psmeHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", psmeHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "Michaelis-Menten", psmeHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", psmeHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", psmeHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", psmeHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", psmeHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", psmeHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Richards", psmeHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Richards", psmeHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", psmeHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", psmeHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", psmeHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", psmeHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", psmeHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", psmeHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", psmeHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", psmeHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", psmeHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", psmeHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards GNLS", psmeHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     bind_rows(c(method = "Chapman-Richards BAL GNLS", psmeHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton GNLS", psmeHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton BAL GNLS", psmeHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang GNLS", psmeHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang BAL GNLS", psmeHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", psmeHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", psmeHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                              mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", psmeDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", psmeDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", psmeDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", psmeDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form", psmeDiameterFromHeightChapmanForm$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form AAL", psmeDiameterFromHeightChapmanFormAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL", psmeDiameterFromHeightChapmanFormBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL RelHt", psmeDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form RelHt", psmeDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                     bind_rows(c(method = "modified Sharma-Parton", psmeDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "linear", psmeDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "parabolic", psmeDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", psmeDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", psmeDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", psmeDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", psmeDiameterFromHeightPowerRelHt$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", psmeDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", psmeDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", psmeDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", psmeDiameterFromHeightSibbesenFormRelHt$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "PSME", 
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


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
