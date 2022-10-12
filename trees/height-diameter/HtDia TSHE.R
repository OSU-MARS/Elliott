# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western hemlock height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Zhang, Sharma-Parton, Chapman-Richards BAL
#tsheHeightFromDiameterRichards = gsl_nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 34.0, d = 1.38, kU = 0.0272), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterChapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterRichards = gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d - dp * isPlantation) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp * isPlantation)^((d + dp * isPlantation)/(1 - d - dp * isPlantation))))^(1/(1 - d - dp * isPlantation)), tshe2016, start = list(Ha = 47.7, Hap = -16.0, d = 0.993, dp = 0, kU = 0.0165, kUp = 0.0124), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterSharmaPartonBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterSharmaZhang = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 42.7, a1p = -18.0, a2 = 0.068, a2p = 0.076, b1 = -0.0188, b1p = -0.0110, b2 = -0.0134, b2p = -0.0032, b3 = 1.295, b3p = -0.023), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterSharmaZhangBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), weights = pmin(DBH^-2, 1))
#tsheHeightFromDiameterWeibullBal = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 47.2, a2 = -0.0585, a2p = 0.715, a3 = 0.110, a3p = -0.214, b1 = -0.00758, b1p = -0.00758, b2 = 1.238, b2p = -0.027), weights = pmin(DBH^-2, 1))
tshe2016physio = tshe2016 %>% filter(is.na(elevation) == FALSE)
tshe2016plantationPhysio = tshe2016physio %>% filter(isPlantation)
tsheHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048), weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24), weights = pmin(DBH^-2, 1)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect)) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 62.5, a2 = 0.092, a2p = 0.51, a3 = -0.021, a4 = -0.136, a5 = 0.057, a6 = 0.957, b1 = -0.022, b1p = 0.003, b2 = 1.31, b2p = -0.090), weights = pmin(DBH^-2, 1)) # a1p, a5, a6 not significant
tsheHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a4 = 63.0, a4p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01), weights = pmin(DBH^-2, 1)) # a3, a3p, b1p not significant
tsheHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58.5, a1p = -3.5, a2 = -0.015, a3 = -8.07, a4 = 0.761, a5 = 0.718, a6 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12), weights = pmin(DBH^-2, 1)) # a4, a5, a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054), weights = pmin(DBH^-2, 1)) # b1 not significant
tsheHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047), weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 815, a1p = 590, b1 = -8.23, b2 = -0.234, b2p = 0.040), weights = pmin(DBH^-2, 1)) # b1p not significant
tsheHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), tshe2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264), weights = pmin(DBH^-2, 1)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), tshe2016, offset = breastHeight, weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93), weights = pmin(DBH^-2, 1)) # a3p not significant
tsheHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05), weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1), weights = pmin(DBH^-2, 1)) # b2p not significant
tsheHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013), weights = pmin(DBH^-2, 1)) # dp not significant
tsheHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = list(a1 = 30.7, a2 = 0.15, a2p = -0.068, b1 = -0.034, b2 = -0.22, b2p = 0.184, b3 = 1.26), weights = pmin(DBH^-2, 1)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, tshe2016, start = list(a1 = 46.7, a1p = -12.6, a2 = 0.054, b1 = -0.021, b1p = -0.013, b2 = -0.054, b3 = 1.26), weights = pmin(DBH^-2, 1)) # a2p, b2p, b3p not significant
tsheHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, tshe2016physio, start = list(a1 = 48.7, a1p = -11.5, a2 = 0.075, a3 = -0.0003, a4 = 0.0089, a5 = 0.0042, a6 = -0.0031, b1 = -0.022, b1p = -0.013, b2 = -0.056, b3 = 1.26), weights = pmin(DBH^-2, 1)) # a2p, a4, a5, b2p, b3p not significant
tsheHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^b2*DBH))^b3, tshe2016physio, start = list(a1 = 42.8, a1p = -9.7, a2 = 0.101, a3 = -0.0003, a4 = 0.0098, a5 = 0.0051, a6 = -0.0032, b1 = -0.022, b1p = -0.012, b2 = -0.038, b3 = 1.27), weights = pmin(DBH^-2, 1)) # a2p, a4, a5, b2p, b3p not significant
tsheHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = list(a1 = 33.4, a1p = -8.1, a2 = 0.138, b1 = -0.027, b2 = -0.053, b2p = 0.077, b3 = 1.27), weights = pmin(DBH^-2, 1)) # a2p, b1p, b3p not significant
tsheHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, tshe2016, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.100, a2p = -0.207, a3 = 0.0005, a3p = 0.0191, b1 = -0.020, b2 = -0.029, b3 = 1.24), weights = pmin(DBH^-2, 1)) # a2, a3, b1p, b2p, b3p not significant
tsheHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014), weights = pmin(DBH^-2, 1)) # b1p not significant
tsheHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029), weights = pmin(DBH^-2, 1))
tsheHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212), weights = pmin(DBH^-2, 1)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * pmin(relativeHeight, 1.25)) * (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), tshe2016, start = list(a1 = -1.59, a1p = 14.5, a2 = 0.005, a2p = 0.351, a3 = 0.002, a4 = 62.8, a4p = -32.6, b1 = -1.169, b1p = 1.14, b2 = 1.02), weights = pmin(DBH^-2, 1), control = list(maxiter = 200)) # a3p, b2p not significant
#confint2(tsheHeightFromDiameterWeibullBalRelHt, level = 0.99)

tsheHeightFromDiameterChapmanRichards = get_height_error(tsheHeightFromDiameterChapmanRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBal = get_height_error(tsheHeightFromDiameterChapmanRichardsBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBalPhysio = get_height_error(tsheHeightFromDiameterChapmanRichardsBalPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterChapmanRichardsBalRelHt = get_height_error(tsheHeightFromDiameterChapmanRichardsBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsPhysio = get_height_error(tsheHeightFromDiameterChapmanRichardsPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterCurtis = get_height_error(tsheHeightFromDiameterCurtis, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterHossfeld = get_height_error(tsheHeightFromDiameterHossfeld, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterKorf = get_height_error(tsheHeightFromDiameterKorf, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterLinear = get_height_error(tsheHeightFromDiameterLinear, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterMichaelisMenten = get_height_error(tsheHeightFromDiameterMichaelisMenten, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterParabolic = get_height_error(tsheHeightFromDiameterParabolic, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterProdan = get_height_error(tsheHeightFromDiameterProdan, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterPower = get_height_error(tsheHeightFromDiameterPower, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterRatkowsky = get_height_error(tsheHeightFromDiameterRatkowsky, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterRichards = get_height_error(tsheHeightFromDiameterRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaParton = get_height_error(tsheHeightFromDiameterSharmaParton, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBal = get_height_error(tsheHeightFromDiameterSharmaPartonBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBalPhysio = get_height_error(tsheHeightFromDiameterSharmaPartonBalPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterSharmaPartonPhysio = get_height_error(tsheHeightFromDiameterSharmaPartonPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterSharmaZhang = get_height_error(tsheHeightFromDiameterSharmaZhang, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangBal = get_height_error(tsheHeightFromDiameterSharmaZhangBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSibbesen = get_height_error(tsheHeightFromDiameterSibbesen, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibull = get_height_error(tsheHeightFromDiameterWeibull, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBal = get_height_error(tsheHeightFromDiameterWeibullBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBalRelHt = get_height_error(tsheHeightFromDiameterWeibullBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)

tsheHeightFromDiameterResults = bind_rows(as_row("Chapman-Richards", tsheHeightFromDiameterChapmanRichards),
                                          as_row("Chapman-Richards BAL", tsheHeightFromDiameterChapmanRichardsBal),
                                          as_row("Chapman-Richards BAL physio", tsheHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row("Chapman-Richards BAL RelHt", tsheHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row("Chapman-Richards physio", tsheHeightFromDiameterChapmanRichardsPhysio),
                                          as_row("Curtis", tsheHeightFromDiameterCurtis),
                                          as_row("Hossfeld", tsheHeightFromDiameterHossfeld),
                                          as_row("Korf", tsheHeightFromDiameterKorf),
                                          as_row("linear", tsheHeightFromDiameterLinear),
                                          as_row("generalized Michaelis-Menten", tsheHeightFromDiameterMichaelisMenten),
                                          as_row("parabolic", tsheHeightFromDiameterParabolic),
                                          as_row("power", tsheHeightFromDiameterPower),
                                          as_row("Prodan", tsheHeightFromDiameterProdan),
                                          as_row("Ratkowsky", tsheHeightFromDiameterRatkowsky),
                                          as_row("unified Richards", tsheHeightFromDiameterRichards),
                                          as_row("Sharma-Parton", tsheHeightFromDiameterSharmaParton),
                                          as_row("Sharma-Parton BAL", tsheHeightFromDiameterSharmaPartonBal),
                                          as_row("Sharma-Parton BAL physio", tsheHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row("Sharma-Parton physio", tsheHeightFromDiameterSharmaPartonPhysio),
                                          as_row("Sharma-Zhang", tsheHeightFromDiameterSharmaZhang),
                                          as_row("Sharma-Zhang BAL", tsheHeightFromDiameterSharmaZhangBal),
                                          as_row("Sibbesen", tsheHeightFromDiameterSibbesen),
                                          as_row("Weibull", tsheHeightFromDiameterWeibull),
                                          as_row("Weibull BAL", tsheHeightFromDiameterWeibullBal),
                                          as_row("Weibull BAL RelHt", tsheHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "DBH", species = "TSHE", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
print(tsheHeightFromDiameterResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterSharmaZhang$fitted.values, color = "Sharma-Zhang", group = tshe2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterSharmaParton$fitted.values, color = "Sharma-Parton", group = tshe2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterCurtis$fitted.values, color = "Curtis", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterKorf$fitted.values, color = "Korf", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterLinear$fitted.values, color = "linear", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterMichaelisMenten$fitted.values, color = "generalized Michaelis-Menten", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterParabolic$fitted.values, color = "parabolic", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterPower$fitted.values, color = "power", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterProdan$fitted.values, color = "Prodan", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterRatkowsky$fitted.values, color = "Ratkowsky", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterRichards$fitted.values, color = "unified Richards", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = tshe2016$isPlantation)) +
  annotate("text", x = 0, y = 70, label = "western hemlock, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(ylim = c(0, 70)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))

## western hemlock height-diameter GNLS regressions
#tsheHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 65.3, a1p = -13.1, b1 = -0.022, b2 = 1.51, b2p = -0.31), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 22.6, a2 = 0.26, a2p = -0.050, b1 = -0.021, b1p = -0.014, b2 = 0.025, b2p = -0.187, b3 = 1.51, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 56.1, a1p = -23.1, a2 = 0.042, a2p = 0.117, b1 = -0.0247, b1p = -0.0131, b2 = -0.0217, b2p = -0.112, b3 = 1.476, b3p = -0.456), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 63.6, a1p = -12.7, b1 = -0.00516, b1p = -0.00652, b2 = 1.29, b2p = -0.16), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#tsheHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 63.6, a2 = 0.035, a2p = 0.832, a3 = 0.0120, a3p = -0.184, b1 = -0.0052, b1p = -0.0024, b2 = 1.281, b2p = -0.133), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH), control = gnlsControl(nlsTol = 0.1, tolerance = 0.001, msTol = 0.001, msVerbose = FALSE, returnObject = TRUE))
#save(tsheHeightFromDiameterChapmanRichardsGnls, tsheHeightFromDiameterChapmanRichardsBalGnls, tsheHeightFromDiameterSharmaPartonGnls, tsheHeightFromDiameterSharmaPartonBalGnls, tsheHeightFromDiameterSharmaZhangGnls, tsheHeightFromDiameterSharmaZhangBalGnls, tsheHeightFromDiameterWeibullGnls, tsheHeightFromDiameterWeibullBalGnls, file = "Timber Inventory/HtDia TSHE GNLS.rdata")
load("trees/height-diameter/HtDia TSHE GNLS.rdata")
tsheHeightFromDiameterWeibullGnls = tsheHeightFromDiameterWykoffGnls # temporary naming error fixup
tsheHeightFromDiameterWeibullBalGnls = tsheHeightFromDiameterWykoffBalGnls

tsheHeightFromDiameterChapmanRichardsGnls = get_height_error(tsheHeightFromDiameterChapmanRichardsGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBalGnls = get_height_error(tsheHeightFromDiameterChapmanRichardsBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonGnls = get_height_error(tsheHeightFromDiameterSharmaPartonGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBalGnls = get_height_error(tsheHeightFromDiameterSharmaPartonBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangGnls = get_height_error(tsheHeightFromDiameterSharmaZhangGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangBalGnls = get_height_error(tsheHeightFromDiameterSharmaZhangBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullGnls = get_height_error(tsheHeightFromDiameterWeibullGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBalGnls = get_height_error(tsheHeightFromDiameterWeibullBalGnls, tshe2016, tshe2016natural, tshe2016plantation)

tsheHeightFromDiameterResultsGnls = bind_rows(as_row("Chapman-Richards GNLS", tsheHeightFromDiameterChapmanRichardsGnls),
                                              as_row("Chapman-Richards BAL GNLS", tsheHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row("Sharma-Parton GNLS", tsheHeightFromDiameterSharmaPartonGnls),
                                              as_row("Sharma-Parton BAL GNLS", tsheHeightFromDiameterSharmaPartonBalGnls),
                                              as_row("Sharma-Zhang GNLS", tsheHeightFromDiameterSharmaZhangGnls),
                                              as_row("Sharma-Zhang BAL GNLS", tsheHeightFromDiameterSharmaZhangBalGnls),
                                              as_row("Weibull GNLS", tsheHeightFromDiameterWeibullGnls),
                                              as_row("Weibull BAL GNLS", tsheHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "DBH", species = "TSHE", deltaAic = aic - min(aic)) %>%
  relocate(responseVariable, species) %>%
  arrange(desc(deltaAic))
tsheHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)


#bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(tsheHeightFromDiameterWeibullBAL, level = 0.99), balN = confint2(tsheHeightFromDiameterWeibullBalNatural, level = 0.99), balP = confint2(tsheHeightFromDiameterWeibullBalPlantation, level = 0.99)) %>%
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
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterWeibullBAL$fitted.values, color = "ElliottBAL"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = tshe2016natural$DBH, y = tsheHeightFromDiameterWeibullBALnatural$fitted.values, color = "ElliottBALn"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = tshe2016plantation$DBH, y = tsheHeightFromDiameterWeibullBALplantation$fitted.values, color = "ElliottBALp"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterBase$fitted.values, color = "base")) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterWeibull$fitted.values, color = "ElliottWeibull")) +
  annotate("text", x = 0, y = 85, label = "a) western hemlock, height from diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BAL", "Weibull with BAL, natural regeneration", "Weibull with BAL, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## western hemlock diameter-height regressions
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, b2p = 0, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 200)) # step factor
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, a3 = 0, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 400)) # step factor
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 130, a1p = -55.1, a2 = -1.63, a2p = -1.06, a3 = 0.61, a3p = -0.4, b1 = 0.10, b2 = 0.50, b2p = 0.16), weights = pmin(TotalHt^-2, 0.5), control = list(maxiter = 400), trace = TRUE) # step factor
#tsheDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 985, a1p = 262, a2 = -10.6, a2p = -16.6, a3 = 4.3, a3p = 3.0, a4 = -384, a4p = 79, b1 = 0.022, b2 = 1.02, b2p = -0.100), weights = pmin(TotalHt^-2, 0.5)) # plantation effects not significant
#tsheDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, tshe2016, iter = 100,
#                                                   start_lower = list(a1 = 0.01, a2 = 0.01, b1 = -10, b2 = -1, b3 = -1), 
#                                                   start_upper = list(a1 = 10, a2 = 3, b1 = 10, b2 = 1, b3 = 2), modelweights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), tshe2016, start = list(a1 = -119, a1p = 26.4, b1 = 0.0196, b1p = 0.0004, b2 = 0.847), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanRichardsAal = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), tshe2016, start = list(a1 = -119, a1p = 26.4, a2 = 0, b1 = 0.0196, b1p = 0.0004, b2 = 0.847), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), tshe2016physio, start = list(a1 = -189, a1p = 71.1, a2 = -0.017, a3 = -0.339, a4 = -0.129, a5 = 0.427, a6 = -0.008, b1 = 0.0095, b1p = 0.00375, b2 = 0.919), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5, a6 not significant
tsheDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.999999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a2 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 136, b1 = 0.0100, b2 = 0.924), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
tsheDiameterFromHeightChapmanFormAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 135, a2 = -0.0045, b1 = 0.010, b2 = 0.924), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant, a3 + b2p step size
tsheDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 189, a2 = -3.64, a3 = 2.24, a4 = -46.4, b1 = 0.0137, b2 = 0.847), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 128, a2 = 5.3, b1 = 0.0154, b2 = 0.902), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37), tshe2016, weights = TotalHt^-2) # isPlantation*(TotalHt - 1.37) not significant(p = 0.036)
tsheDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83), weights = TotalHt^-2) # a1p, a2p, b1p not significant
tsheDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013), weights = TotalHt^-2)
tsheDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016, weights = TotalHt^-2)
tsheDiameterFromHeightPower = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
tsheDiameterFromHeightPowerAal = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047), weights = pmin(TotalHt^-2, 0.5)) # a1p not significant
tsheDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a2 = 0.00006, a3 = 0.284, a4 = 0.00044, a5 = -0.00056, a6 = 0.00045, b1 = 1.04), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5 a6 not significant
tsheDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a2 = 0.08, b1 = 1.02), weights = pmin(TotalHt^-2, 0.5)) 
tsheDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.31, b1 = 0.818, b2 = 0.010), weights = pmin(TotalHt^-2, 0.5)) # a1p, b1p, b2p not significant
#tsheDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32), weights = pmin(TotalHt^-2, 0.5)) # converges from red alder values but fails to reconverge (singular gradient)
#tsheDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, tshe2016, start = list(a1 = 350, a2 = 1.3, b1 = 0.0025, b2 = -0.015, b3 = 0.702), weights = pmin(TotalHt^-2, 0.5), control = gsl_nls_control(maxiter = 250)) # unlreiable convergence and NaN-inf with gsl_nls(), NaN-inf with nlrob(), singular gradient with nls() even with nls_multstart() parameters
#tsheDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 1.901, a2 = 0.934, a2p = 0, b1 = 19.35, b2 = -0.395, b2p = 0, b3 = 0.00052, b3p = 0), weights = pmin(TotalHt^-2, 0.5)) # singular gradient
tsheDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
tsheDiameterFromHeightSibbesenFormAal = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578), weights = pmin(TotalHt^-2, 0.5)) # no significant plantation effects
tsheDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), tshe2016physio, start = list(a1 = 2.115, a1p = -0.100, a2 = 0.00017, a3 = 0.486, a4 = 0.00326, a5 = -0.0051, a6 = -0.0001, b1 = 0.736, b2 = 0.0593, b2p = 0.0028), weights = pmin(TotalHt^-2, 0.5)) # a2, a4, a5, a6 not significant
tsheDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a2 = 0.084, b1 = 0.75, b2 = 0.056), weights = pmin(TotalHt^-2, 0.5))
tsheDiameterFromHeightWeibull = nlrob(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -225, b1 = 0.011, b2 = 0.88), weights = pmin(TotalHt^-2, 0.5)) # a1p, b1p, b2p not significant
#confint2(tsheDiameterFromHeightPowerPhysio, level = 0.99)

tsheDiameterFromHeightChapmanRichards = get_dbh_error(tsheDiameterFromHeightChapmanRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanRichardsAal = get_dbh_error(tsheDiameterFromHeightChapmanRichardsAal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanRichardsPhysio = get_dbh_error(tsheDiameterFromHeightChapmanRichardsPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightChapmanRichardsRelHt = get_dbh_error(tsheDiameterFromHeightChapmanRichardsRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanForm = get_dbh_error(tsheDiameterFromHeightChapmanForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormAal = get_dbh_error(tsheDiameterFromHeightChapmanFormAal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormBal = get_dbh_error(tsheDiameterFromHeightChapmanFormBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormBalRelHt = get_dbh_error(tsheDiameterFromHeightChapmanFormBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormRelHt = get_dbh_error(tsheDiameterFromHeightChapmanFormRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightLinear = get_dbh_error(tsheDiameterFromHeightLinear, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightMichaelisMentenForm = get_dbh_error(tsheDiameterFromHeightMichaelisMentenForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightNaslund = get_dbh_error(tsheDiameterFromHeightNaslund, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightParabolic = get_dbh_error(tsheDiameterFromHeightParabolic, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPower = get_dbh_error(tsheDiameterFromHeightPower, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPowerAal = get_dbh_error(tsheDiameterFromHeightPowerAal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPowerPhysio = get_dbh_error(tsheDiameterFromHeightPowerPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightPowerRelHt = get_dbh_error(tsheDiameterFromHeightPowerRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightRuark = get_dbh_error(tsheDiameterFromHeightRuark, tshe2016, tshe2016natural, tshe2016plantation)
#tsheDiameterFromHeightSchnute = get_dbh_error(tsheDiameterFromHeightSchnute, tshe2016, tshe2016natural, tshe2016plantation)
#tsheDiameterFromHeightSharmaParton = get_dbh_error(tsheDiameterFromHeightSharmaParton, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenForm = get_dbh_error(tsheDiameterFromHeightSibbesenForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenFormAal = get_dbh_error(tsheDiameterFromHeightSibbesenFormAal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenFormPhysio = get_dbh_error(tsheDiameterFromHeightSibbesenFormPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightSibbesenFormRelHt = get_dbh_error(tsheDiameterFromHeightSibbesenFormRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightWeibull = get_dbh_error(tsheDiameterFromHeightWeibull, tshe2016, tshe2016natural, tshe2016plantation)

tsheDiameterFromHeightResults = bind_rows(as_row("Chapman-Richards", tsheDiameterFromHeightChapmanRichards),
                                          as_row("Chapman-Richards AAL", tsheDiameterFromHeightChapmanRichardsAal),
                                          as_row("Chapman-Richards physio", tsheDiameterFromHeightChapmanRichardsPhysio),
                                          as_row("Chapman-Richards RelHt", tsheDiameterFromHeightChapmanRichardsRelHt),
                                          as_row("Chapman-Richards form", tsheDiameterFromHeightChapmanForm),
                                          as_row("Chapman-Richards form AAL", tsheDiameterFromHeightChapmanFormAal), # not AIC supported
                                          as_row("Chapman-Richards form BAL", tsheDiameterFromHeightChapmanFormBal),
                                          as_row("Chapman-Richards form BAL RelHt", tsheDiameterFromHeightChapmanFormBalRelHt),
                                          as_row("Chapman-Richards form RelHt", tsheDiameterFromHeightChapmanFormRelHt), # not AIC supported
                                          as_row("linear", tsheDiameterFromHeightLinear),
                                          as_row("generalized Michaelis-Menten form", tsheDiameterFromHeightMichaelisMentenForm),
                                          as_row("Näslund", tsheDiameterFromHeightNaslund),
                                          as_row("parabolic", tsheDiameterFromHeightParabolic),
                                          as_row("power", tsheDiameterFromHeightPower),
                                          as_row("power AAL", tsheDiameterFromHeightPowerAal),
                                          as_row("power physio", tsheDiameterFromHeightPowerPhysio),
                                          as_row("power RelHt", tsheDiameterFromHeightPowerRelHt),
                                          as_row("modified Sharma-Parton", NULL),
                                          as_row("Ruark", tsheDiameterFromHeightRuark),
                                          as_row("Schnute", NULL),
                                          as_row("Sibbesen form", tsheDiameterFromHeightSibbesenForm),
                                          as_row("Sibbesen form AAL", tsheDiameterFromHeightSibbesenFormAal), # not AIC supported
                                          as_row("Sibbesen form physio", tsheDiameterFromHeightSibbesenFormPhysio),
                                          as_row("Sibbesen form RelHt", tsheDiameterFromHeightSibbesenFormRelHt), # not AIC supported
                                          as_row("Weibull", tsheDiameterFromHeightWeibull)) %>% 
  mutate(responseVariable = "height", species = "TSHE", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(tsheDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(tshe2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = tsheDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanFormAal$fitted.values, y = TotalHt, color = "Chapman-Richards form AAL", group = isPlantation), alpha = 0.5) +
  geom_line(aes(x = tsheDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "generalized Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards form top height", group = isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 81, label = "western hemlock, diameter from height", hjust = 0, size = 3.5) +
  #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
  #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))


## collect model parameters
tsheParameters = bind_rows(bind_rows(bind_rows(c(method = "Chapman-Richards", tsheHeightFromDiameterChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL", tsheHeightFromDiameterChapmanRichardsBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL physio", tsheHeightFromDiameterChapmanRichardsBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards BAL RelHt", tsheHeightFromDiameterChapmanRichardsBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", tsheHeightFromDiameterChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Curtis", tsheHeightFromDiameterCurtis$m$getPars())),
                                     bind_rows(c(method = "Hossfeld", tsheHeightFromDiameterHossfeld$m$getPars())),
                                     bind_rows(c(method = "Korf", tsheHeightFromDiameterKorf$m$getPars())),
                                     bind_rows(c(method = "linear", tsheHeightFromDiameterLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten", tsheHeightFromDiameterMichaelisMenten$m$getPars())),
                                     bind_rows(c(method = "parabolic", tsheHeightFromDiameterParabolic$coefficients)),
                                     bind_rows(c(method = "power", tsheHeightFromDiameterPower$m$getPars())),
                                     bind_rows(c(method = "Prodan", tsheHeightFromDiameterProdan$m$getPars())),
                                     bind_rows(c(method = "Ratkowsky", tsheHeightFromDiameterRatkowsky$m$getPars())),
                                     bind_rows(c(method = "Richards", tsheHeightFromDiameterRichards$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton", tsheHeightFromDiameterSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL", tsheHeightFromDiameterSharmaPartonBal$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton BAL physio", tsheHeightFromDiameterSharmaPartonBalPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Parton physio", tsheHeightFromDiameterSharmaPartonPhysio$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang", tsheHeightFromDiameterSharmaZhang$m$getPars())),
                                     bind_rows(c(method = "Sharma-Zhang BAL", tsheHeightFromDiameterSharmaZhangBal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen", tsheHeightFromDiameterSibbesen$m$getPars())),
                                     bind_rows(c(method = "Weibull", tsheHeightFromDiameterWeibull$m$getPars())),
                                     bind_rows(c(method = "Weibull BAL", tsheHeightFromDiameterWeibullBal$m$getPars())),
                                     bind_rows(c(method = "Weibull RelHt", tsheHeightFromDiameterWeibullBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards GNLS", tsheHeightFromDiameterChapmanRichardsGnls$coefficients)),
                                     bind_rows(c(method = "Chapman-Richards BAL GNLS", tsheHeightFromDiameterChapmanRichardsBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton GNLS", tsheHeightFromDiameterSharmaPartonGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Parton BAL GNLS", tsheHeightFromDiameterSharmaPartonBalGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang GNLS", tsheHeightFromDiameterSharmaZhangGnls$coefficients)),
                                     bind_rows(c(method = "Sharma-Zhang BAL GNLS", tsheHeightFromDiameterSharmaZhangBalGnls$coefficients)),
                                     bind_rows(c(method = "Weibull GNLS", tsheHeightFromDiameterWeibullGnls$coefficients)),
                                     bind_rows(c(method = "Weibull BAL GNLS", tsheHeightFromDiameterWeibullBalGnls$coefficients))) %>%
                             mutate(responseVariable = "DBH"),
                           bind_rows(bind_rows(c(method = "Chapman-Richards", tsheDiameterFromHeightChapmanRichards$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards AAL", tsheDiameterFromHeightChapmanRichardsAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards physio", tsheDiameterFromHeightChapmanRichardsPhysio$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards RelHt", tsheDiameterFromHeightChapmanRichardsRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form", tsheDiameterFromHeightChapmanForm$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form AAL", tsheDiameterFromHeightChapmanFormAal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL", tsheDiameterFromHeightChapmanFormBal$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form BAL RelHt", tsheDiameterFromHeightChapmanFormBalRelHt$m$getPars())),
                                     bind_rows(c(method = "Chapman-Richards form RelHt", tsheDiameterFromHeightChapmanFormRelHt$m$getPars())),
                                     bind_rows(c(method = "linear", tsheDiameterFromHeightLinear$coefficients)),
                                     bind_rows(c(method = "generalized Michaelis-Menten form", tsheDiameterFromHeightMichaelisMentenForm$m$getPars())),
                                     bind_rows(c(method = "Näslund", tsheDiameterFromHeightNaslund$m$getPars())),
                                     bind_rows(c(method = "parabolic", tsheDiameterFromHeightParabolic$coefficients)),
                                     bind_rows(c(method = "power", tsheDiameterFromHeightPower$m$getPars())),
                                     bind_rows(c(method = "power AAL", tsheDiameterFromHeightPowerAal$m$getPars())),
                                     bind_rows(c(method = "power physio", tsheDiameterFromHeightPowerPhysio$m$getPars())),
                                     bind_rows(c(method = "power RelHt", tsheDiameterFromHeightPowerRelHt$m$getPars())),
                                     bind_rows(c(method = "Ruark", tsheDiameterFromHeightRuark$m$getPars())),
                                     #bind_rows(c(method = "Schnute", tsheDiameterFromHeightSchnute$m$getPars())),
                                     #bind_rows(c(method = "modified Sharma-Parton", tsheDiameterFromHeightSharmaParton$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form", tsheDiameterFromHeightSibbesenForm$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form AAL", tsheDiameterFromHeightSibbesenFormAal$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form physio", tsheDiameterFromHeightSibbesenFormPhysio$m$getPars())),
                                     bind_rows(c(method = "Sibbesen form RelHt", tsheDiameterFromHeightSibbesenFormRelHt$m$getPars())),
                                     bind_rows(c(method = "Weibull", tsheDiameterFromHeightWeibull$m$getPars()))) %>%
                             mutate(responseVariable = "height")) %>%
  mutate(species = "TSHE",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3), b3p = as.numeric(b3p)) %>%
  relocate(responseVariable, species, method, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)

