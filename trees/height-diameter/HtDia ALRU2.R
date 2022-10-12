# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## red alder height-diameter regression form sweep
# R options for nonlinear least squares
#   gslnls::gsl_nls() - fast but fixed weights only, defaults to Levenberg-Marquadt
#   minpack.lm::nlsLM() - weights = wfcs() fails with index out of range, Levenberg-Marquadt only
#   nlme::gnls() - varPower() based weighting fragile andconvergence usually fails
#   robustbase::nlrob() - defaults to iterative reweighted least squares with stats::nls()
#   stats::nls() - fixed weights only, defaults to Gauss-Newton, NL2SOL with algorithm = "port"
# Difference between stats::nls() with fixed weighting and robustbase::nlrob() is relatively small but nlrob() 
# does noticeably decrease asymptoticity of most regression forms (slight decreases occur in some cases).
#alruHeightFromDiameterMichaelisMenten = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH / (a2 + a2p * isPlantation + DBH), alru2016, start = list(a1 = 45.5, a1p = 18.7, a2 = 49.1, a2p = 20.5), weights = pmin(DBH^-2, 1))
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, d = 1.17, kU = 0.0366), weights = pmin(DBH^-2, 1))
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, kU = 0.0366, kUp = 0), weights = pmin(DBH^-2, 1)) # confint2() step factor
#alruHeightFromDiameterRichards = nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - (d + dp*isPlantation)) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp*isPlantation)^((d + dp*isPlantation)/(1 - d))))^(1/(1 - (d + dp*isPlantation))), alru2016, start = list(Ha = 22.8, Hap = 0, d = 1.172, dp = 0, kU = 0.0366, kUp = 0), weights = pmin(DBH^-2, 1)) # confint2() NaN-inf
alru2016physio = alru2016 %>% filter(is.na(elevation) == FALSE)
alru2016plantationPhysio = alru2016physio %>% filter(isPlantation)
alruHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 26.4, a1p = 2.74, b1 = -0.041, b2 = 1.11, b2p = 0.027), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 28.0, a1p = 5.93, a2 = 0.187, a2p = 0.061, a3 = -0.165, a3p = 0.118, b1 = -0.0438, b1p = 0.0164, b2 = 1.174, b2p = -0.119), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect)) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 35.1, a1p = 6.02, a2 = 0.040, a2p = 0.131, a3 = -0.0050, a4 = -0.241, a5 = 0.364, a6 = 0.102, b1 = -0.040, b1p = -0.040, b2 = 1.09), weights = pmin(DBH^-2, 1)) # a5, a6 not significant
alruHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = -1.39, a1p = 12.5, a2 = 0.020, a2p = 0.368, a3 = -0.00222, a4 = 56.2, a4p = -27.6, b1 = -0.022, b2 = 0.026, b2p = 0.750), weights = pmin(DBH^-2, 1)) # a3p not significant
alruHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, alru2016physio, start = list(a1 = 33.7, a1p = 4.73, a2 = -0.00025, a3 = -0.195, a4 = 0.195, a5 = -0.113, a6 = 0.200, b1 = -0.044, b1p = 0.0054, b2 = 1.14), weights = pmin(DBH^-2, 1)) # b2p not significant
alruHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, alru2016, start = list(a1 = 1.6, b1 = 0.24), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 30.7, a1p = 13.3, b1 = 49.1, b1p = 7.95, b2 = -1.25, b2p = 0.12), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), alru2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
alruHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 30.4, a1p = 11.2, a2 = 54.1, b1 = 1.29, b1p = -0.152), weights = pmin(DBH^-2, 1)) # a2p not significant
alruHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I((isPlantation*DBH)^2), alru2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
alruHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + a2*DBH + a3 + a3p* isPlantation), alru2016, start = list(a1 = 0.0249, a1p = -0.0054, a2 = 0.938, a3 = 0.761, a3p = -0.177), weights = pmin(DBH^-2, 1)) # a2p not significant
alruHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + a1*DBH^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 1.13, b1 = 0.786, b1p = 0.047), weights = pmin(DBH^-2, 1)) # a1p not significant
alruHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), alru2016, start = list(a1 = 32.9, a1p = -7.85, b1 = -21.2, b1p = -5.41, b2 = 4.96, b2p = 1.58), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), alru2016, start = list(Ha = 22.3, Hap = 0.972, d = 1.196, kU = 0.0361), weights = pmin(DBH^-2, 1))
# below coefficients not synced
alruHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 16.5, a1p = 12.4, a2 = 0.33, a2p = -0.162, a3 = -0.00008, a4 = 0.0090, a5 = 0.0045, a6 = 0.00256, b1 = -0.020, b1p = -0.0091, b2 = 0.062, b2p = -0.235, b3 = 1.50, b3p = -0.45), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016physio, start = list(a1 = 17.6, a1p = 5.06, a2 = 0.31, a2p = -0.106, a3 = -0.00008, a4 = 0.0092, a5 = 0.0046, a6 = 0.00257, b1 = -0.023, b1p = -0.012, b2 = 0.0010, b2p = -0.159, b3 = 1.52, b3p = -0.46), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.0019, a1p = 0.163, b1 = 4.98, b1p = -2.72, b2 = -0.175, b2p = 0.0427), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 102, a1p = 92, b1 = -17.7, b1p = 10.3, b2 = -0.725, b2p = 0.365), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), weights = pmin(DBH^-2, 1))
alruHeightFromDiameterWeibullBalRelHt = gsl_nls(TotalHt ~ 1.37 + (a1  + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp(b1 * (1 + (b2 + b2p * isPlantation) * pmin(relativeHeight, 1.25)) * DBH^b3)), alru2016, start = list(a1 = 5.30, a1p = -4.53, a2 = 0.194, a2p = -0.1068, a3 = -0.155, a3p = 0.256, a4 = 52.7, a4p = -21.2, b1 = -0.136, b2 = -0.740, b2p = 0.105, b3 = 0.878), weights = pmin(DBH^-2, 1)) # nlrob() failed to converge
#confint2(alruHeightFromDiameterChapmanRichardsBalPhysio, level = 0.99)

alruHeightFromDiameterChapmanRichards = get_height_error(alruHeightFromDiameterChapmanRichards, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBal = get_height_error(alruHeightFromDiameterChapmanRichardsBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(alruHeightFromDiameterChapmanRichardsBalPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(alruHeightFromDiameterChapmanRichardsBalRelHt, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsPhysio = get_height_error(alruHeightFromDiameterChapmanRichardsPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterCurtis = get_height_error(alruHeightFromDiameterCurtis, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterHossfeld = get_height_error(alruHeightFromDiameterHossfeld, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterKorf = get_height_error(alruHeightFromDiameterKorf, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterLinear = get_height_error(alruHeightFromDiameterLinear, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterMichaelisMenten = get_height_error(alruHeightFromDiameterMichaelisMenten, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterParabolic = get_height_error(alruHeightFromDiameterParabolic, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterProdan = get_height_error(alruHeightFromDiameterProdan, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterPower = get_height_error(alruHeightFromDiameterPower, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterRatkowsky = get_height_error(alruHeightFromDiameterRatkowsky, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterRichards = get_height_error(alruHeightFromDiameterRichards, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaParton = get_height_error(alruHeightFromDiameterSharmaParton, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBal = get_height_error(alruHeightFromDiameterSharmaPartonBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBalPhysio = get_height_error(alruHeightFromDiameterSharmaPartonBalPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterSharmaPartonPhysio = get_height_error(alruHeightFromDiameterSharmaPartonPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruHeightFromDiameterSharmaZhang = get_height_error(alruHeightFromDiameterSharmaZhang, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangBal = get_height_error(alruHeightFromDiameterSharmaZhangBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSibbesen = get_height_error(alruHeightFromDiameterSibbesen, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibull = get_height_error(alruHeightFromDiameterWeibull, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBal = get_height_error(alruHeightFromDiameterWeibullBal, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBalRelHt = get_height_error(alruHeightFromDiameterWeibullBalRelHt, alru2016, alru2016natural, alru2016plantation)

alruHeightFromDiameterResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(alruHeightFromDiameterChapmanRichards),
                                        "Chapman-Richards BAL", !!!as_row(alruHeightFromDiameterChapmanRichardsBal),
                                        "Chapman-Richards BAL physio", !!!as_row(alruHeightFromDiameterChapmanRichardsBalPhysio),
                                        "Chapman-Richards BAL RelHt", !!!as_row(alruHeightFromDiameterChapmanRichardsBalRelHt),
                                        "Chapman-Richards physio", !!!as_row(alruHeightFromDiameterChapmanRichardsPhysio),
                                        "Curtis", !!!as_row(alruHeightFromDiameterCurtis),
                                        "Hossfeld", !!!as_row(alruHeightFromDiameterHossfeld),
                                        "Korf", !!!as_row(alruHeightFromDiameterKorf),
                                        "linear", !!!as_row(alruHeightFromDiameterLinear),
                                        "generalized Michaelis-Menten", !!!as_row(alruHeightFromDiameterMichaelisMenten),
                                        "parabolic", !!!as_row(alruHeightFromDiameterParabolic),
                                        "power", !!!as_row(alruHeightFromDiameterPower),
                                        "Prodan", !!!as_row(alruHeightFromDiameterProdan),
                                        "Ratkowsky", !!!as_row(alruHeightFromDiameterRatkowsky),
                                        "unified Richards", !!!as_row(alruHeightFromDiameterRichards),
                                        "Sharma-Parton", !!!as_row(alruHeightFromDiameterSharmaParton),
                                        "Sharma-Parton BAL", !!!as_row(alruHeightFromDiameterSharmaPartonBal),
                                        "Sharma-Parton BAL physio", !!!as_row(alruHeightFromDiameterSharmaPartonBalPhysio),
                                        "Sharma-Parton physio", !!!as_row(alruHeightFromDiameterSharmaPartonPhysio),
                                        "Sharma-Zhang", !!!as_row(alruHeightFromDiameterSharmaZhang),
                                        "Sharma-Zhang BAL", !!!as_row(alruHeightFromDiameterSharmaZhangBal),
                                        "Sibbesen", !!!as_row(alruHeightFromDiameterSibbesen),
                                        "Weibull", !!!as_row(alruHeightFromDiameterWeibull),
                                        "Weibull BAL", !!!as_row(alruHeightFromDiameterWeibullBal),
                                        "Weibull BAL RelHt", !!!as_row(alruHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(alruHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterKorf$fitted.values, color = "Korf", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterLinear$fitted.values, color = "linear", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterPower$fitted.values, color = "power", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterProdan$fitted.values, color = "Prodan", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = alru2016$isPlantation)) +
  annotate("text", x = 0, y = 50, label = "red alder, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## red alder height-diameter GNLS regressions
#alruHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), alru2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), alru2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(alruHeightFromDiameterChapmanRichardsGnls, alruHeightFromDiameterChapmanRichardsBalGnls, alruHeightFromDiameterSharmaPartonGnls, alruHeightFromDiameterSharmaPartonBalGnls, alruHeightFromDiameterSharmaZhangGnls, alruHeightFromDiameterSharmaZhangBalGnls, alruHeightFromDiameterWeibullGnls, alruHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia ALRU2 GNLS.rdata")
load("trees/height-diameter/HtDia ALRU2 GNLS.rdata")
alruHeightFromDiameterWeibullGnls = alruHeightFromDiameterWykoffGnls # temporary naming error fixup
alruHeightFromDiameterWeibullBalGnls = alruHeightFromDiameterWykoffBalGnls

alruHeightFromDiameterChapmanRichardsGnls = get_height_error(alruHeightFromDiameterChapmanRichardsGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterChapmanRichardsBalGnls = get_height_error(alruHeightFromDiameterChapmanRichardsBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonGnls = get_height_error(alruHeightFromDiameterSharmaPartonGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaPartonBalGnls = get_height_error(alruHeightFromDiameterSharmaPartonBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangGnls = get_height_error(alruHeightFromDiameterSharmaZhangGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterSharmaZhangBalGnls = get_height_error(alruHeightFromDiameterSharmaZhangBalGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullGnls = get_height_error(alruHeightFromDiameterWeibullGnls, alru2016, alru2016natural, alru2016plantation)
alruHeightFromDiameterWeibullBalGnls = get_height_error(alruHeightFromDiameterWeibullBalGnls, alru2016, alru2016natural, alru2016plantation)

alruHeightFromDiameterResultsGnls = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                            "Chapman-Richards GNLS", !!!as_row(alruHeightFromDiameterChapmanRichardsGnls),
                                            "Chapman-Richards BAL GNLS", !!!as_row(alruHeightFromDiameterChapmanRichardsBalGnls),
                                            "Sharma-Parton GNLS", !!!as_row(alruHeightFromDiameterSharmaPartonGnls),
                                            "Sharma-Parton BAL GNLS", !!!as_row(alruHeightFromDiameterSharmaPartonBalGnls),
                                            "Sharma-Zhang GNLS", !!!as_row(alruHeightFromDiameterSharmaZhangGnls),
                                            "Sharma-Zhang BAL GNLS", !!!as_row(alruHeightFromDiameterSharmaZhangBalGnls),
                                            "Weibull GNLS", !!!as_row(alruHeightFromDiameterWeibullGnls),
                                            "Weibull BAL GNLS", !!!as_row(alruHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "ALRU2", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
alruHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichardsGnls$fitted.values, color = "Chapman-Richards GNLS", group = alru2016$isPlantation)) +
  annotate("text", x = 0, y = 50, label = "a) red alder, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 115), ylim = c(0, 50)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

## red alder diameter-height regressions
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, start = list(a1 = -18.3, b1 = 0.0305, b2 = 0.427), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50))
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.9999)), alru2016, start = list(a1 = -15.8, b1 = 0.0304, b2 = 0.362), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 50))
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ a1*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, algorithm = "port",
#                                            lower = list(a1 = -150, b1 = 0.02, b2 = 0.1), 
#                                            start = list(a1 = -50, b1 = 0.03, b2 = 0.7), 
#                                            upper = list(a1 = -1, b1 = 0.05, b2 = 2), weights = pmin(TotalHt^-2, 0.5), control = list(warnOnly = TRUE))
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ a1 * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100,
#                                                      start_lower = list(a1 = -20, b1 = -0.03, b2 = -1), 
#                                                      start_upper = list(a1 = -10, b1 = 0.03, b2 = 1), modelweights = pmin(TotalHt^-2, 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 2), modelweights = pmin(TotalHt^-2, 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^-2, 0.5)) # b2p not significant
#alruDiameterFromHeightChapmanRichards = nls_multstart(DBH ~ (a1 + a1p * isPlantation) * log(1 - pmin(((b1 + b1p * isPlantation)*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.99)), alru2016, iter = 100,
#                                                      lower = c(a1 = -25, a1p = -10, b1 = 0.015, b1p = -0.02, b2 = 0.01, b2p = -0.01),
#                                                      start_lower = list(a1 = -20, a1p = -1, b1 = 0.02, b1p = 0.00, b2 = 0.3, b2p = -0.1), 
#                                                      start_upper = list(a1 = -10, a1p = 1, b1 = 0.03, b1p = 0.01, b2 = 1.5, b2p = 0), 
#                                                      upper = c(a1 = -5, a1p = 10, b1 = 0.035, b1p = 0.02, b2 = 3, b2p = 0.01), modelweights = pmin(TotalHt^-2, 0.5)) # b1p not significant
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-1, 0.7))
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2, 0.7))
#alruDiameterFromHeightChapmanRichardsMod = nls(DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.999)), alru2016, algorithm = "port",
#                                               lower = list(a1 = -200, b1 = 0.02, b2 = 0.1),
#                                               start = list(a1 = -30, b1 = 0.03, b2 = 1.0), weights = pmin(TotalHt^-2.6, 0.7))
#alruDiameterFromHeightChapmanRichardsPhysio = nls_multstart(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.1415926/180 * slope) + a4 * cos(3.1415926/180 * aspect) + a5 * sin(3.1415926/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999)), alru2016, iter = 100,
#                                                            start_lower = list(a1 = 10, a1p = 5, a2 = -0.01, a3 = -2, a4 = -2, a5 = -0.5, a6 = -0.1, b1 = -0.1, b1p = -0.1, b2 = -2), # 3.14159, if used, is seen as a mutable parameter
#                                                            start_upper = list(a1 = 20, a1p = 15, a2 = 0.01, a3 = 1, a4 = 2, a5 = 0.5, a6 = 0.1, b1 = 0.1, b1p = 0.1, b2 = 2), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightChapmanRichardsRelHt = nls_multstart(DBH ~ (a1 + a2 * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^b2, 0.999)), alru2016, iter = 100, 
#                                                           lower = c(a1 = -100, a2 = -10, b1 = 1/60, b2 = 0),
#                                                           start_lower = list(a1 = -35, a2 = -1, b1 = 0.025, b2 = 0.4),
#                                                           start_upper = list(a1 = -25, a2 = 1, b1 = 0.035, b2 = 0.8), 
#                                                           upper = c(a1 = 0.1, a2 = 100, b1 = 1/25, b2 = 2), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightKorf = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), lower = list(a1 = 0.5, b1 = 1, b2 = 0.3), weights = pmin(TotalHt^-2, 0.5), algorithm = "port")
#alruDiameterFromHeightKorf = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = pmin(TotalHt^-1.3, 0.7), control = list(maxiter = 200, tol = 0.001, warnOnly = TRUE))
#alruDiameterFromHeightKorf = gnls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 1, b1 = 1, b2 = 0.4), weights = varPower(0.65, ~TotalHt), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#alruDiameterFromHeightLinear = lm(DBH ~ exp(1*(TotalHt - 1.37)^0.4) + 0, alru2016, weights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightLinear = lm(DBH ~ I(TotalHt - 1.37) + 0, alru2016, weights = pmin(TotalHt^-2, 0.5))
#summary(alruDiameterFromHeightLinear)
#alruDiameterFromHeightPower = nls(DBH ~ a1*(TotalHt - 1.37)^b1, alru2016, algorithm = "port",
#                                  lower = list(a1 = 0, b1 = 1.2), # b1 < 1 is asymptotic in diameter rather than height
#                                  start = list(a1 = 0.4, b1 = 1.5), 
#                                  upper = list(a1 = 10, b2 = 3), weights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightRuark = nls(DBH ~ a1*(TotalHt - 1.37)^b1*exp(b2*(TotalHt - 1.37)), alru2016, start = list(a1 = 1, b1 = 0.04, b2 = 0.04), weights = pmin(TotalHt^-2, 0.5)) # diverges to concave up
#alruDiameterFromHeightSchnute = nls(DBH ~ a1*log(1 - b1*(TotalHt^b2 - 1.37^b2) / (40^b2 - 1.37^b2)), alru2016, start = list(a1 = -150, b1 = 0.5, b2 = 1.5), weights = pmin(TotalHt^-2, 0.5)) # diverges to concave up
#alruDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, alru2016, iter = 100, 
#                                                   start_lower = list(a1 = 0.1, a2 = 0.1, b1 = 0.001, b2 = 0.1, b3 = 0.1), 
#                                                   start_upper = list(a1 = 100, a2 = 2, b1 = 0.1, b2 = 1, b3 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightSibbesen = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), alru2016, iter = 100, 
#                                               lower = c(a1 = 0.01, b1 = 0.2, b2 = 0.2),
#                                               start_lower = list(a1 = 2, b1 = 0.25, b2 = 0.25),
#                                               start_upper = list(a1 = 4, b1 = 0.30, b2 = 0.30), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightWeibull = nls(DBH ~ a1*exp(b1*(TotalHt - 1.37)^b2), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53), weights = pmin(TotalHt^-2, 0.5), trace = TRUE) # diverges to concave up
#alruDiameterFromHeight = nls(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = 5, b1 = 0.4, b2 = 0.53), weights = pmin(TotalHt^-2, 0.5)) # step factor
#alruDiameterFromHeightWykoff = nls_multstart(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                             lower = c(a1 = 1, b1 = 1/55, b2 = 0.2),
#                                             start_lower = list(a1 = 3, b1 = 0.2, b2 = 0.5), 
#                                             start_upper = list(a1 = 7, b1 = 0.4, b2 = 0.7), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightWykoffAal = nls_multstart(DBH ~ (a1 + a2 * tallerQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a2 = -1, b1 = -0.05, b2 = -1), 
#                                                start_upper = list(a1 = 1, a2 = 1, b1 = 0.05, b2 = 1), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightWykoffAal = nls_multstart(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea + (a3 + a3p * isPlantation) * standQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, iter = 100,
#                                                start_lower = list(a1 = -100, a1p = -100, a2 = -1, a2p = -1, a3 = -1, a3p = -1, b1 = -1, b2 = -1, b2p = -1), 
#                                                start_upper = list(a1 = 1, a1p = 1, a2 = 1, a2p = 1, a3 = 1, a3p = 1, b1 = 1, b2 = 2, b2p = 1), modelweights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightYang = nls(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, alru2016, algorithm = "port",
#                                 lower = list(a1 = -1000, b1 = 0.001, b2 = 0.4), # b1 constraint prevents NaN-inf
#                                 start = list(a1 = -150, b1 = 0.02, b2 = 0.75), weights = pmin(TotalHt^-2, 0.5)) # collapses to linear fit
# coefficients for physiologically plausible curvature but poor accuracy
#alruDiameterFromHeightChapmanRichards = nls(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -16.9, a1p = 2.676, b1 = 0.0307, b2 = 0.0304, b2p = 0.0723), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 210)) # b1p not significant, >200 iterations to converge to start point
#alruDiameterFromHeightChapmanRichardsAal = nls(DBH ~ (a1 + a2 * tallerQuasiBasalArea + a3 * standQuasiBasalArea)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -15.4, a2 = 0.234, a3 = 0.233, b1 = 0.0342, b2 = 0.263, b2p = 0.164), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 500)) # a1p, a2, a2p, a3, a3p not significant
#alruDiameterFromHeightChapmanRichardsPhysio = nls(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -10.5, a1p = 2.93, a2 = 0.0030, a3 = -18.3, a4 = -0.095, a5 = -0.42, a6 = 0.137, b1 = 0.0307, b2 = 0.342, b2p = 0.0683), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 200)) # a2, a4, a5, a6, b1p not significant
#alruDiameterFromHeightChapmanRichardsRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*log(1 - pmin((b1*(TotalHt - 1.37))^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = -36.7, a1p = 30.7, a2 = 33.2, a2p = -40.8, b1 = 0.0307, b2 = 0.562, b2p = -0.360), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 200))
# coefficients for physically implausible accuracy
alruDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -57.5, b1 = -0.024, b2 = 1.38, b2p = -0.22), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
alruDiameterFromHeightChapmanFormAal = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerQuasiBasalArea) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = 51.1, a2 = -0.057, a2p = 0, b1 = -0.024, b2 = 1.31, b2p = 0), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant, all NaN-inf with a3 * standQuasiBasalArea
alruDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -32.1, a1p = 4.8, a2 = 1.72, a2p = -0.45, a3 = -1.46, a3p = 0.45, b1 = -0.040, b2 = 1.31), weights = pmin(TotalHt^-2, 0.5)) # b2p not significant
alruDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), alru2016, start = list(a1 = -45.2, a1p = 3.3, a2 = 2.21, a2p = -0.70, a3 = -1.85, a3p = 0.743, a4 = 29.8, a4p = -14.8, b1 = -0.0302, b2 = 1.30), weights = pmin(TotalHt^-2, 0.5)) # NaN-inf
alruDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), alru2016, start = list(a1 = -54.2, a2 = 0.49, b1 = -0.020, b2 = 1.49, b2p = -0.19), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
alruDiameterFromHeightChapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 39.1, b1 = -0.028, b2 = 2.68, b2p = -1.50), weights = pmin(TotalHt^-2, 0.7)) # b1p not significant
alruDiameterFromHeightChapmanRichardsAal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016, start = list(a1 = 9.7, a1p = 42.7, a2 = -0.028, b1 = -0.026, b2 = 2.81, b2p = -1.64), weights = pmin(TotalHt^-2, 0.7), control = list(maxiter = 50)) # a2p not significant
alruDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation), 0.9999)), alru2016physio, start = list(a1 = 6.3, a1p = 30.9, a2 = -0.0004, a3 = 7.37, a4 = 0.41, a5 = 0.52, a6 = 0.037, b1 = -0.031, b2 = 2.47, b2p = -1.26), weights = pmin(TotalHt^-2, 0.7)) # a2, a4, a6 not significant
alruDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), alru2016, start = list(a1 = 33.0, a1p = -7.77, a2 = -0.95, b1 = -0.043, b2 = 1.39), weights = pmin(TotalHt^-2, 0.7)) # a2p, b1p not significant
alruDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)), alru2016, weights = TotalHt^-2)
alruDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), alru2016, weights = TotalHt^-2)
alruDiameterFromHeightPower = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.98, a1p = -2.03, b1 = 0.78, b1p = 0.15), weights = pmin(TotalHt^-2, 0.5))
alruDiameterFromHeightPowerAal = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 4.15, a1p = -2.22, a2 = -0.0012, a2p = 0.00032, b1 = 0.79, b1p = 0.147), weights = pmin(TotalHt^-2, 0.5))
alruDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation  + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, alru2016physio, start = list(a1 = 2.38, a1p = -0.70, a2 = -0.00002, a3 = 1.26, a4 = 0.078, a5 = 0.0075, a6 = 0.00017, b1 = 0.88), weights = pmin(TotalHt^-2, 0.5)) # a2, a6, b1p not significant
alruDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), alru2016, start = list(a1 = 3.55, a1p = -1.73, a2 = -0.98, a2p = 0.65, b1 = 0.86, b1p = 0.14), weights = pmin(TotalHt^-2, 0.5))
#alruDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, alru2016, start = list(a1 = 0.0005, a2 = 3.35, b1 = 0.019, b2 = 0.0758, b3 = -2.12), weights = pmin(TotalHt^-2, 0.5)) # NaN-infinity even with parameters from nls_multstart()
alruDiameterFromHeightSibbesenForm = nlrob(DBH ~ (a1 + a1p * isPlantation)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.263, a1p = 0.522, b1 = 3.349, b1p = -1.629, b2 = -0.226, b2p = 0.119), weights = pmin(TotalHt^-2, 0.5))
alruDiameterFromHeightSibbesenFormAal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^((b1 + b1p * isPlantation)*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.25, a1p = 0.56, a2 = -0.0001, b1 = 3.43, b1p = -1.74, b2 = -0.23, b2p = 0.12), weights = pmin(TotalHt^-2, 0.5))  # a2 not significant
alruDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016physio, start = list(a1 = 0.716, a1p = -0.336, a2 = -0.00002, a3 = 0.231, a4 = 0.191, a5 = 0.018, a6 = 0.0007, b1 = 2.13, b2 = -0.163, b2p = 0.0197), weights = pmin(TotalHt^-2, 0.5)) # a2, a6 not significant
alruDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.712, a1p = -0.341, a2 = -0.0437, b1 = 2.379, b2 = -0.182, b2p = 0.0348), weights = pmin(TotalHt^-2, 0.5)) # neither a2 not significant
#confint2(alruDiameterFromHeightChapmanFormBalRelHt, level = 0.99)

alruDiameterFromHeightChapmanForm = get_dbh_error(alruDiameterFromHeightChapmanForm, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormAal = get_dbh_error(alruDiameterFromHeightChapmanFormAal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormBal = get_dbh_error(alruDiameterFromHeightChapmanFormBal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(alruDiameterFromHeightChapmanFormBalRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanFormRelHt = get_dbh_error(alruDiameterFromHeightChapmanFormRelHt, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichards = get_dbh_error(alruDiameterFromHeightChapmanRichards, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichardsAal = get_dbh_error(alruDiameterFromHeightChapmanRichardsAal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(alruDiameterFromHeightChapmanRichardsPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(alruDiameterFromHeightChapmanRichardsRelHt, alru2016, alru2016natural, alru2016plantation)
#alruDiameterFromHeightKorf = get_dbh_error(alruDiameterFromHeightKorf, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightLinear = get_dbh_error(alruDiameterFromHeightLinear, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightParabolic = get_dbh_error(alruDiameterFromHeightParabolic, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPower = get_dbh_error(alruDiameterFromHeightPower, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPowerAal = get_dbh_error(alruDiameterFromHeightPowerAal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightPowerPhysio = get_dbh_error(alruDiameterFromHeightPowerPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightPowerRelHt = get_dbh_error(alruDiameterFromHeightPowerRelHt, alru2016, alru2016natural, alru2016plantation)
#alruDiameterFromHeightSharmaParton = get_dbh_error(alruDiameterFromHeightSharmaParton, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenForm = get_dbh_error(alruDiameterFromHeightSibbesenForm, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenFormAal = get_dbh_error(alruDiameterFromHeightSibbesenFormAal, alru2016, alru2016natural, alru2016plantation)
alruDiameterFromHeightSibbesenFormPhysio = get_dbh_error(alruDiameterFromHeightSibbesenFormPhysio, alru2016physio, alru2016natural, alru2016plantationPhysio)
alruDiameterFromHeightSibbesenFormRelHt = get_dbh_error(alruDiameterFromHeightSibbesenFormRelHt, alru2016, alru2016natural, alru2016plantation)

alruDiameterFromHeightResults = tribble(~method, ~pae, ~paeNR, ~paePl, ~bias, ~biasNR, ~biasPl, ~mae, ~maeNR, ~maePl, ~rmse, ~rmseNR, ~rmsePl, ~nse, ~nseNR, ~nsePl, ~pearson, ~pearsonNR, ~pearsonPl, ~aic, ~bic, ~power,
                                        "Chapman-Richards", !!!as_row(alruDiameterFromHeightChapmanRichards),
                                        "Chapman-Richards AAL", !!!as_row(alruDiameterFromHeightChapmanRichardsAal), # not accuracy supported
                                        "Chapman-Richards physio", !!!as_row(alruDiameterFromHeightChapmanRichardsPhysio),
                                        "Chapman-Richards RelHt", !!!as_row(alruDiameterFromHeightChapmanRichardsRelHt), # not AIC supported
                                        #"Korf", !!!as_row(alruDiameterFromHeightKorf),
                                        "linear", !!!as_row(alruDiameterFromHeightLinear),
                                        "parabolic", !!!as_row(alruDiameterFromHeightParabolic),
                                        "power", !!!as_row(alruDiameterFromHeightPower),
                                        "power AAL", !!!as_row(alruDiameterFromHeightPowerAal),
                                        "power physio", !!!as_row(alruDiameterFromHeightPowerPhysio),
                                        "power RelHt", !!!as_row(alruDiameterFromHeightPowerRelHt),
                                        #"Ruark", !!!as_row(alruDiameterFromHeightRuark),
                                        "modified Sharma-Parton", !!!as_row(NULL),
                                        "Sibbesen form", !!!as_row(alruDiameterFromHeightSibbesenForm),
                                        "Sibbesen form AAL", !!!as_row(alruDiameterFromHeightSibbesenFormAal), # not accuracy supported
                                        "Sibbesen form physio", !!!as_row(alruDiameterFromHeightSibbesenFormPhysio),
                                        "Sibbesen form RelHt", !!!as_row(alruDiameterFromHeightSibbesenFormRelHt), # not AIC supported
                                        "Chapman-Richards form", !!!as_row(alruDiameterFromHeightChapmanForm),
                                        "Chapman-Richards form AAL", !!!as_row(alruDiameterFromHeightChapmanFormAal), # not accuracy supported
                                        "Chapman-Richards form BAL", !!!as_row(alruDiameterFromHeightChapmanFormBal),
                                        "Chapman-Richards form BAL RelHt", !!!as_row(alruDiameterFromHeightChapmanFormBalRelHt),
                                        "Chapman-Richards form RelHt", !!!as_row(alruDiameterFromHeightChapmanFormRelHt)) %>%# not AIC supported
  mutate(responseVariable = "height", species = "ALRU2", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
alruDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic)

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_smooth(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  geom_line(aes(x = alruDiameterFromHeightChapmanRichards$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsAal$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards AAL", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightChapmanRichardsRelHt$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightKorf$fitted.values, y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightLinear$fitted.values, y = alru2016$TotalHt, color = "linear", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightParabolic$fitted.values, y = alru2016$TotalHt, color = "parabolic", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightPower$fitted.values, y = alru2016$TotalHt, color = "power", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alruDiameterFromHeightRuark$fitted.values, y = alru2016$TotalHt, color = "Ruark form", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightSibbesenForm$fitted.values, y = alru2016$TotalHt, color = "Sibbesen form", group = alru2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightChapmanForm$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards form", group = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = -30*log(1 - pmin((0.03*(alru2016$TotalHt - 1.37))^0.85, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards")) +
  #geom_line(aes(x = (-16.95 + 2.676 * alru2016$isPlantation)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^(0.304 + 0.0723 * alru2016$isPlantation), 0.9999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99999", group = alru2016$isPlantation)) +
  #geom_line(aes(x = (-15.886)*log(1 - pmin((0.0307*(alru2016$TotalHt - 1.37))^0.304, 0.999)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.999")) +
  #geom_line(aes(x = (-18.329)*log(1 - pmin((0.0305*(alru2016$TotalHt - 1.37))^0.427, 0.99)), y = alru2016$TotalHt, color = "Chapman-Richards fit 0.99")) +
  #geom_line(aes(x = -30*log(1 - pmin(0.03*(alru2016$TotalHt - 1.37)^1.0, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
  #geom_line(aes(x = -150*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.79, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 1")) +
  #geom_line(aes(x = -119*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.86, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 2")) +
  #geom_line(aes(x = -108*log(1 - pmin(0.02*(alru2016$TotalHt - 1.37)^0.89, 0.999)), y = alru2016$TotalHt, color = "modified Chapman-Richards 3")) +
  #geom_line(aes(x = -33 + 3*alru2016$relativeHeight)*log(1 - pmin((0.03*(alru2016$TotalHt - 1.37))^0.85, 0.999), y = alru2016$TotalHt, color = "Chapman-Richards RelHt", group = alru2016$isPlantation)) +
  #geom_line(aes(x = (0.016*(alru2016$TotalHt - 1.37))^0.027 / (1 - (0.016*(alru2016$TotalHt - 1.37))^0.027), y = alru2016$TotalHt, color = "Curtis")) +
  #geom_line(aes(x = 1*exp(0.7*(alru2016$TotalHt - 1.37)^0.5), y = alru2016$TotalHt, color = "Korf")) +
  #geom_line(aes(x = -1/0.08*log((40 - alru2016$TotalHt)/alru2016$TotalHt) + 30, y = alru2016$TotalHt, color = "logistic")) +
  #geom_line(aes(x = 10 / (1 + 1*exp(-0.03*(alru2016$TotalHt - 1.37))), y = alru2016$TotalHt, color = "Pearl-Reed")) +
  #geom_line(aes(x = 0.3*(alru2016$TotalHt - 1.37)^1.5, y = alru2016$TotalHt, color = "power")) +
  #geom_line(aes(x = 1*(alru2016$TotalHt - 1.37)^1.2, y = alru2016$TotalHt, color = "power")) +
  #geom_line(aes(x = 30*(alru2016$TotalHt - 1.37) / (60 - (alru2016$TotalHt - 1.37)), y = alru2016$TotalHt, color = "Ratkowsky")) +
  #geom_line(aes(x = 2*(alru2016$TotalHt - 1.37)^0.5 / (1 - 0.12*(alru2016$TotalHt - 1.37)^0.5), y = alru2016$TotalHt, color = "Ratkowsky")) +
  #geom_line(aes(x = 100*(alru20 16$TotalHt - 1.37)^0.01*(exp(0.01*(alru2016$TotalHt - 1.37)) - 1), y = alru2016$TotalHt, color = "Ruark form")) +
  #geom_line(aes(x = -150*log(1 - 0.4*(alru2016$TotalHt^1.5 - 1.3^1.5) / (40^1.5 - 1.3^1.5)), y = alru2016$TotalHt, color = "Schnute")) +
  #geom_line(aes(x = 3*(alru2016$TotalHt - 1.37)^(0.30*(alru2016$TotalHt - 1.37)^0.30), y = alru2016$TotalHt, color = "Sibbesen form")) +
  #geom_line(aes(x = 5.3*(alru2016$TotalHt - 1.37)^(0.35*(alru2016$TotalHt - 1.37)^0.20), y = alru2016$TotalHt, color = "Sibbesen form")) +
  #geom_line(aes(x = 5*(exp(0.33*(alru2016$TotalHt - 1.37)^0.60) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
  #geom_line(aes(x = 39*(exp(0.1*(alru2016$TotalHt - 1.37)^0.66) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
  #geom_line(aes(x = 125.8*(exp(0.0182*(alru2016$TotalHt - 1.37)^0.880) - 1), y = alru2016$TotalHt, color = "Chapman-Richards form")) +
  #geom_line(aes(x = (-150*log(1 - pmin(0.013*(alru2016$TotalHt - 1.37), 0.9999)))^0.75, y = alru2016$TotalHt, color = "Yang")) +
  geom_line(aes(x = (-1000*log(1 - pmin(0.002*(alru2016$TotalHt - 1.37), 0.9999)))^0.95, y = alru2016$TotalHt, color = "Yang")) +
  annotate("text", x = 0, y = 50, label = "red alder, diameter from height", hjust = 0, size = 3.5) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

ggplot(alru2016) +
  geom_point(aes(x = TotalHt, y = DBH), alpha = 0.10, color = "grey25", shape = 16) +
  geom_line(aes(x = TotalHt, y = -50*log(1 - pmin(0.004*(TotalHt - 1.37)^1.5, 0.9)), color = "Chapman-Richards")) +
  annotate("text", x = 0, y = 90, label = "red alder, diameter from height", hjust = 0, size = 3.5) +
  labs(x = "height, m", y = "DBH, cm", color = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))

ggplot(alru2016) +
  geom_bin_2d(aes(x = DBH, y = TotalHt), binwidth = c(2.5, 1)) +
  #geom_line(aes(x = 1.246*exp(1*(TotalHt - 1.37)^0.4), y = TotalHt, color = "linear")) +
  #geom_line(aes(x = alruDiameterFromHeightKorf$fitted.values, y = TotalHt, color = "Korf")) +
  geom_line(aes(x = 1.246*exp(1*(TotalHt - 1.37)^0.4), y = TotalHt, color = "Korf manual")) +
  geom_line(aes(x = 0.155*exp(2.708*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf unweighted")) +
  geom_line(aes(x = 0.146*exp(2.732*(TotalHt - 1.37)^0.22), y = TotalHt, color = "Korf weighted")) +
  geom_line(aes(x = 0.019*exp(4.691*(TotalHt - 1.37)^0.15), y = TotalHt, color = "Korf weighted to step factor")) +
  geom_line(aes(x = 0.051*exp(3.724*(TotalHt - 1.37)^0.18), y = TotalHt, color = "Korf GNLS")) +
  labs(x = "DBH, cm", y = "height, m", color = "regression\nform", fill = "stems\nmeasured") +
  scale_fill_viridis_c(trans = "log10")

ggplot() +
  geom_point(aes(x = alru2016$TotalHt, y = alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$TotalHt, y = alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  #geom_point(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "grey25", shape = 16) +
  #geom_smooth(aes(x = alru2016$TotalHt, y = 1/alru2016$TotalHt * alruDiameterFromHeightKorf$residuals), alpha = 0.1, color = "red", formula = y ~ s(x, k = 20), method = "gam") +
  labs(x = "height, m", y = "DBH error residual, cm")

ggplot(alru2016) + 
  geom_histogram(aes(x = DBH, y = 100 * ..count../sum(..count..)), binwidth = 2.5) +
  labs(x = "DBH, cm", y = "percentage of stems measured") +
ggplot(alru2016) + 
  geom_histogram(aes(x = 100 * ..count../sum(..count..), y = TotalHt), binwidth = 1) +
  labs(x = "percentage of stems measured", y = "height, m")

ggplot(alru2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = DBH, y = 40 / (1 + exp(-0.08*(DBH - 30))), color = "logistic")) +
  geom_line(aes(x = DBH, y = 1.37 + (40 - 1.37) / (1 + exp(-0.08*(DBH - 30)))^1.1, color = "generalized logistic")) +
  #geom_line(aes(x = seq(0, 100, by = 0.1), y = 1.37 + 40 * (1 - exp(-0.01*seq(0, 100, by = 0.1))^1.5), color = "Korf")) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## collect model parameters
alruParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", alruHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", alruHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", alruHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", alruHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", alruHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", alruHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", alruHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", alruHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", alruHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten", alruHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", alruHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", alruHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", alruHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", alruHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", alruHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", alruHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", alruHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", alruHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", alruHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", alruHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", alruHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", alruHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", alruHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", alruHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", alruHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards GNLS", alruHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     bind_rows(c(method = "Chapman-Richards BAL GNLS", alruHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton GNLS", alruHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton BAL GNLS", alruHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang GNLS", alruHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang BAL GNLS", alruHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", alruHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", alruHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", alruDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", alruDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", alruDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", alruDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form", alruDiameterFromHeightChapmanForm$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form AAL", alruDiameterFromHeightChapmanFormAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL", alruDiameterFromHeightChapmanFormBal$m$getPars())),
                                     #bind_rows(c(method = "Chapman-Richards form BAL RelHt", alruDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form RelHt", alruDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                     bind_rows(c(method = "linear", alruDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "parabolic", alruDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", alruDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", alruDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", alruDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", alruDiameterFromHeightPowerRelHt$m$getPars())),
                                     #bind_rows(c(method = "modified Sharma-Parton", alruDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", alruDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", alruDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", alruDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", alruDiameterFromHeightSibbesenFormRelHt$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "ALRU2",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

