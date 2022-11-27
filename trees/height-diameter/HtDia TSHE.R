# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R

## western hemlock height-diameter regression form sweep
# preferred forms: Sharma-Parton BAL, Sharma-Zhang, Sharma-Parton, Chapman-Richards BAL
#tsheHeightFromDiameterRichards = gsl_nls(TotalHt ~ 1.37 + Ha * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-kU * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 34.0, d = 1.38, kU = 0.0272), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterChapmanRichardsBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 64.9, a1p = 3.8, a2 = 0.023, a2p = 0.92, a3 = 0.022, a3p = -0.22, b1 = -0.021, b1p = 0.0066, b2 = 1.47, b2p = -0.29), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterRichards = gsl_nls(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d - dp * isPlantation) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/(d + dp * isPlantation)^((d + dp * isPlantation)/(1 - d - dp * isPlantation))))^(1/(1 - d - dp * isPlantation)), tshe2016, start = list(Ha = 47.7, Hap = -16.0, d = 0.993, dp = 0, kU = 0.0165, kUp = 0.0124), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterSharmaPartonBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 18.5, a1p = 11.3, a2 = 0.30, a2p = -0.14, b1 = -0.019, b1p = -0.011, b2 = 0.089, b2p = -0.266, b3 = 1.49, b3p = -0.44), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterSharmaZhang = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 42.7, a1p = -18.0, a2 = 0.068, a2p = 0.076, b1 = -0.0188, b1p = -0.0110, b2 = -0.0134, b2p = -0.0032, b3 = 1.295, b3p = -0.023), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterSharmaZhangBal = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp((b1 + b1p * isPlantation)*tph^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 56.3, a1p = 14.7, a2 = 0.0412, a2p = -0.0535, a3 = 0.0146, a3p = 0.0146, b1 = -0.0249, b1p = -0.00024, b2 = -0.0240, b2p = -0.0969, b3 = 1.48, b3p = -0.370), weights = pmin(DBH^-0.9, 1))
#tsheHeightFromDiameterWeibullBal = gsl_nls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 47.2, a2 = -0.0585, a2p = 0.715, a3 = 0.110, a3p = -0.214, b1 = -0.00758, b1p = -0.00758, b2 = 1.238, b2p = -0.027), weights = pmin(DBH^-0.9, 1))
tshe2016 = trees2016 %>% filter(Species == "WH", isLiveUnbroken, TotalHt > 0) # live western hemlocks measured for height
tshe2016natural = tshe2016 %>% filter(isPlantation == FALSE)
tshe2016physio = tshe2016 %>% filter(is.na(elevation) == FALSE)
tshe2016plantation = tshe2016 %>% filter(isPlantation)
tshe2016plantationPhysio = tshe2016physio %>% filter(isPlantation)

tsheHeightFromDiameterChapmanRichards = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * (1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = 53.7, a1p = -10.7, b1 = -0.021, b1p = -0.006, b2 = 1.30, b2p = -0.048), weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterChapmanRichardsBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = list(a1 = 55.4, a2 = -0.036, a2p = 0.851, a3 = 0.086, a3p = -0.240, b1 = -0.018, b2 = 1.24), weights = pmin(DBH^-0.9, 1)) # a1p, a2, a3, b1p, b2p not significant
tsheHeightFromDiameterChapmanRichardsBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * elevation + a4 * slope + a5 * sin(3.14159/180 * aspect) + a6 * cos(3.14159/180 * aspect) + a7 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 62.5, a2 = 0.092, a2p = 0.51, a3 = -0.021, a4 = -0.136, a5 = 0.057, a6 = 0.957, a7 = 0, b1 = -0.022, b1p = 0.003, b2 = 1.31, b2p = -0.090), weights = pmin(DBH^-0.9, 1)) # a1p, a5, a6 not significant
tsheHeightFromDiameterChapmanRichardsBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016, start = list(a1 = -1.62, a1p = 14.7, a2 = 0.006, a2p = 0.353, a3 = -0.0011, a4 = 63.0, a4p = -32.7, b1 = -0.027, b2 = 0.007, b2p = 1.01), weights = pmin(DBH^-0.9, 1)) # a3, a3p, b1p not significant
tsheHeightFromDiameterChapmanRichardsPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), tshe2016physio, start = list(a1 = 58.5, a1p = -3.5, a2 = -0.015, a3 = -8.07, a4 = 0.761, a5 = 0.718, a6 = -0.058, b1 = -0.026, b2 = 1.36, b2p = -0.12), weights = pmin(DBH^-0.9, 1)) # a4, a5, a6 not significant, b1p+b2p not both significant
tsheHeightFromDiameterCurtis = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH / (1 + DBH)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.55, a1p = 0.16, b1 = -0.021, b1p = 0.054), weights = pmin(DBH^-0.9, 1)) # b1 not significant
tsheHeightFromDiameterHossfeld = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 69.3, a1p = -11.6, b1 = 196, b1p = -73., b2 = -1.30, b2p = 0.047), weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterKorf = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 815, a1p = 590, b1 = -8.23, b2 = -0.234, b2p = 0.040), weights = pmin(DBH^-0.9, 1)) # b1p not significant
tsheHeightFromDiameterLinear = lm(TotalHt ~ 0 + DBH + I(isPlantation*DBH), tshe2016, offset = breastHeight, weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterMichaelisMenten = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264), weights = pmin(DBH^-0.9, 1)) # {a1p, a2p}+b1p not mutually significant
tsheHeightFromDiameterParabolic = lm(TotalHt ~ 0 + DBH + I(DBH^2) + I(isPlantation*DBH) + I(isPlantation*DBH^2), tshe2016, offset = breastHeight, weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterProdan = nlrob(TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation)*DBH^2 + (a2 + a2p * isPlantation)*DBH + a3), tshe2016, start = list(a1 = 0.007, a1p = 0.005, a2 = 1.237, a2p = -0.247, a3 = 1.93), weights = pmin(DBH^-0.9, 1)) # a3p not significant
tsheHeightFromDiameterPower = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 0.59, a1p = -0.15, b1 = 1.02, b1p = -0.05), weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterRatkowsky = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2)), tshe2016, start = list(a1 = 64.1, a1p = -6.0, b1 = -42.1, b1p = 3.8, b2 = 8.1), weights = pmin(DBH^-0.9, 1)) # b2p not significant
tsheHeightFromDiameterRichards = nlrob(TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/Ha)^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), tshe2016, start = list(Ha = 46.7, Hap = -16.7, d = 1.03, kU = 0.017, kUp = 0.013), weights = pmin(DBH^-0.9, 1)) # dp not significant
tsheHeightFromDiameterSharmaParton = nlrob(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = list(a1 = 30.7, a2 = 0.15, a2p = -0.068, b1 = -0.034, b2 = -0.22, b2p = 0.184, b3 = 1.26), weights = pmin(DBH^-0.9, 1)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterSharmaPartonBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, tshe2016, start = list(a1 = 46.7, a1p = -12.6, a2 = 0.054, b1 = -0.021, b1p = -0.013, b2 = -0.054, b3 = 1.26), weights = pmin(DBH^-0.9, 1)) # a2p, b2p, b3p not significant
tsheHeightFromDiameterSharmaPartonBalPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, tshe2016physio, start = list(a1 = 48.7, a1p = -11.5, a2 = 0.075, a3 = -0.0003, a4 = 0.0089, a5 = 0.0042, a6 = -0.0031, b1 = -0.022, b1p = -0.013, b2 = -0.056, b3 = 1.26), weights = pmin(DBH^-0.9, 1)) # a2p, a4, a5, b2p, b3p not significant
tsheHeightFromDiameterSharmaPartonPhysio = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare))^b2*DBH))^b3, tshe2016physio, start = list(a1 = 42.8, a1p = -9.7, a2 = 0.101, a3 = -0.0003, a4 = 0.0098, a5 = 0.0051, a6 = -0.0032, b1 = -0.022, b1p = -0.012, b2 = -0.038, b3 = 1.27), weights = pmin(DBH^-0.9, 1)) # a2p, a4, a5, b2p, b3p not significant
tsheHeightFromDiameterSharmaZhang = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2*(1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = list(a1 = 33.4, a1p = -8.1, a2 = 0.138, b1 = -0.027, b2 = -0.053, b2p = 0.077, b3 = 1.27), weights = pmin(DBH^-0.9, 1)) # a2p, b1p, b3p not significant
tsheHeightFromDiameterSharmaZhangBal = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, tshe2016, start = list(a1 = 42.0, a1p = 32.9, a2 = 0.100, a2p = -0.207, a3 = 0.0005, a3p = 0.0191, b1 = -0.020, b2 = -0.029, b3 = 1.24), weights = pmin(DBH^-0.9, 1)) # a2, a3, b1p, b2p, b3p not significant
tsheHeightFromDiameterSibbesen = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1*DBH^(b2 + b2p * isPlantation)), tshe2016, start = list(a1 = 0.271, a1p = 0.071, b1 = 1.68, b2 = -0.085, b2p = -0.014), weights = pmin(DBH^-0.9, 1)) # b1p not significant
tsheHeightFromDiameterWeibull = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = list(a1 = 49.7, a1p = -10.3, b1 = -0.0076, b1p = -0.0046, b2 = 1.25, b2p = -0.029), weights = pmin(DBH^-0.9, 1))
tsheHeightFromDiameterWeibullBal = nlrob(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 49.4, a2 = -0.055, a2p = 0.781, a3 = 0.091, a3p = -0.231, b1 = -0.008, b2 = 1.212), weights = pmin(DBH^-0.9, 1)) # a1p, b1p, b2p not significant
tsheHeightFromDiameterWeibullBalRelHt = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * pmin(relativeHeight, 1.25)) * (1 - exp(b1*DBH^b2)), tshe2016, start = list(a1 = 16.0, a1p = -9.7, a2 = 0.018, a2p = 0.63, a3 = 0.28, a4 = 73.5, b1 = -0.015, b2 = 0.76), weights = pmin(DBH^-0.9, 1), control = list(maxiter = 200)) # a3p, a4p, b1p, b2p not significant
#confint_nlrob(tsheHeightFromDiameterWeibullBalRelHt, level = 0.99, weights = tshe2016$DBH^-2)

tsheHeightFromDiameterChapmanRichards = get_height_error("Chapman-Richards", tsheHeightFromDiameterChapmanRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBal = get_height_error("Chapman-Richards BAL", tsheHeightFromDiameterChapmanRichardsBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBalPhysio = get_height_error("Chapman-Richards BAL physio", tsheHeightFromDiameterChapmanRichardsBalPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterChapmanRichardsBalRelHt = get_height_error("Chapman-Richards BAL RelHt", tsheHeightFromDiameterChapmanRichardsBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsPhysio = get_height_error("Chapman-Richards physio", tsheHeightFromDiameterChapmanRichardsPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterCurtis = get_height_error("Curtis", tsheHeightFromDiameterCurtis, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterHossfeld = get_height_error("Hossfeld IV", tsheHeightFromDiameterHossfeld, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterKorf = get_height_error("Korf", tsheHeightFromDiameterKorf, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterLinear = get_height_error("linear", tsheHeightFromDiameterLinear, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterMichaelisMenten = get_height_error("Michaelis-Menten", tsheHeightFromDiameterMichaelisMenten, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterParabolic = get_height_error("parabolic", tsheHeightFromDiameterParabolic, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterProdan = get_height_error("Prodan", tsheHeightFromDiameterProdan, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterPower = get_height_error("power", tsheHeightFromDiameterPower, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterRatkowsky = get_height_error("Ratkowsky", tsheHeightFromDiameterRatkowsky, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterRichards = get_height_error("unified Richards", tsheHeightFromDiameterRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaParton = get_height_error("Sharma-Parton", tsheHeightFromDiameterSharmaParton, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBal = get_height_error("Sharma-Parton BAL", tsheHeightFromDiameterSharmaPartonBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBalPhysio = get_height_error("Sharma-Parton BAL physio", tsheHeightFromDiameterSharmaPartonBalPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterSharmaPartonPhysio = get_height_error("Sharma-Parton physio", tsheHeightFromDiameterSharmaPartonPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheHeightFromDiameterSharmaZhang = get_height_error("Sharma-Zhang", tsheHeightFromDiameterSharmaZhang, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangBal = get_height_error("Sharma-Zhang BAL", tsheHeightFromDiameterSharmaZhangBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSibbesen = get_height_error("Sibbesen", tsheHeightFromDiameterSibbesen, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibull = get_height_error("Weibull", tsheHeightFromDiameterWeibull, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBal = get_height_error("Weibull BAL", tsheHeightFromDiameterWeibullBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBalRelHt = get_height_error("Weibull BAL RelHt", tsheHeightFromDiameterWeibullBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)

tsheHeightFromDiameterResults = bind_rows(as_row(tsheHeightFromDiameterChapmanRichards),
                                          as_row(tsheHeightFromDiameterChapmanRichardsBal),
                                          as_row(tsheHeightFromDiameterChapmanRichardsBalPhysio),
                                          as_row(tsheHeightFromDiameterChapmanRichardsBalRelHt),
                                          as_row(tsheHeightFromDiameterChapmanRichardsPhysio),
                                          as_row(tsheHeightFromDiameterCurtis),
                                          as_row(tsheHeightFromDiameterHossfeld),
                                          as_row(tsheHeightFromDiameterKorf),
                                          as_row(tsheHeightFromDiameterLinear),
                                          as_row(tsheHeightFromDiameterMichaelisMenten),
                                          as_row(tsheHeightFromDiameterParabolic),
                                          as_row(tsheHeightFromDiameterPower),
                                          as_row(tsheHeightFromDiameterProdan),
                                          as_row(tsheHeightFromDiameterRatkowsky),
                                          as_row(tsheHeightFromDiameterRichards),
                                          as_row(tsheHeightFromDiameterSharmaParton),
                                          as_row(tsheHeightFromDiameterSharmaPartonBal),
                                          as_row(tsheHeightFromDiameterSharmaPartonBalPhysio),
                                          as_row(tsheHeightFromDiameterSharmaPartonPhysio),
                                          as_row(tsheHeightFromDiameterSharmaZhang),
                                          as_row(tsheHeightFromDiameterSharmaZhangBal),
                                          as_row(tsheHeightFromDiameterSibbesen),
                                          as_row(tsheHeightFromDiameterWeibull),
                                          as_row(tsheHeightFromDiameterWeibullBal),
                                          as_row(tsheHeightFromDiameterWeibullBalRelHt)) %>%
  mutate(responseVariable = "height", species = "TSHE", deltaAic = aic - min(aic)) %>%
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
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterMichaelisMenten$fitted.values, color = "Michaelis-Menten", group = tshe2016$isPlantation)) +
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
#tsheHeightFromDiameterChapmanRichardsGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) *(1 - exp((b1 + b1p*isPlantation)*DBH))^(b2 + b2p * isPlantation), tshe2016, start = tsheHeightFromDiameterChapmanRichards$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.002, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.001
#tsheHeightFromDiameterChapmanRichardsBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, tshe2016, start = tsheHeightFromDiameterChapmanRichardsBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.05, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.02
#tsheHeightFromDiameterSharmaPartonGnls = gnls(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = tsheHeightFromDiameterSharmaParton$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE))
#tsheHeightFromDiameterSharmaPartonBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^a2 * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b2*DBH))^b3, tshe2016, start = tsheHeightFromDiameterSharmaPartonBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msVerbose = FALSE, returnObject = FALSE))
#tsheHeightFromDiameterSharmaZhangGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^a2 * (1 - exp(b1*tph^(b2 + b2p * isPlantation)*DBH))^b3, tshe2016, start = tsheHeightFromDiameterSharmaZhang$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#tsheHeightFromDiameterSharmaZhangBalGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^(a2 + a2p * isPlantation) * (1 + (a3 + a3p * isPlantation) * basalAreaLarger) * (1 - exp(b1*tph^b2*DBH))^b3, tshe2016, start = tsheHeightFromDiameterSharmaZhangBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.02, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, msVerbose = FALSE, returnObject = FALSE)) # step halving at nlsTol = 0.01
#tsheHeightFromDiameterWeibullGnls = gnls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation))), tshe2016, start = tsheHeightFromDiameterWeibull$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.1, msTol = 1E-7, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE)) # >250+50 iterations at nlsTol = 0.1, ok at msTol = 1E-3 and tol = 1E-3
#tsheHeightFromDiameterWeibullBalGnls = gnls(TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), tshe2016, start = tsheHeightFromDiameterWeibullBal$m$getPars(), correlation = corSymm(value = numeric(0.1), form = ~ 0 | PlotID), weights = varPower(0.50, ~DBH | isPlantation), control = gnlsControl(nlsTol = 0.001, msTol = 1E-6, tolerance = 1E-5, maxIter = 250, nlsMaxIter = 50, msVerbose = FALSE, returnObject = FALSE))
#save(tsheHeightFromDiameterChapmanRichardsGnls, tsheHeightFromDiameterChapmanRichardsBalGnls, tsheHeightFromDiameterSharmaPartonGnls, tsheHeightFromDiameterSharmaPartonBalGnls, tsheHeightFromDiameterSharmaZhangGnls, tsheHeightFromDiameterSharmaZhangBalGnls, tsheHeightFromDiameterWeibullGnls, tsheHeightFromDiameterWeibullBalGnls, file = "trees/height-diameter/HtDia TSHE GNLS.rdata")

load("trees/height-diameter/HtDia TSHE GNLS.rdata")
tsheHeightFromDiameterChapmanRichardsGnls = get_height_error("Chapman-Richards GNLS", tsheHeightFromDiameterChapmanRichardsGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterChapmanRichardsBalGnls = get_height_error("Chapman-Richards BAL GNLS", tsheHeightFromDiameterChapmanRichardsBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonGnls = get_height_error("Sharma-Parton GNLS", tsheHeightFromDiameterSharmaPartonGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaPartonBalGnls = get_height_error("Sharma-Parton BAL GNLS", tsheHeightFromDiameterSharmaPartonBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangGnls = get_height_error("Sharma-Zhang GNLS", tsheHeightFromDiameterSharmaZhangGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterSharmaZhangBalGnls = get_height_error("Sharma-Zhang BAL GNLS", tsheHeightFromDiameterSharmaZhangBalGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullGnls = get_height_error("Weibull GNLS", tsheHeightFromDiameterWeibullGnls, tshe2016, tshe2016natural, tshe2016plantation)
tsheHeightFromDiameterWeibullBalGnls = get_height_error("Weibull BAL GNLS", tsheHeightFromDiameterWeibullBalGnls, tshe2016, tshe2016natural, tshe2016plantation)

tsheHeightFromDiameterResultsGnls = bind_rows(as_row(tsheHeightFromDiameterChapmanRichardsGnls),
                                              as_row(tsheHeightFromDiameterChapmanRichardsBalGnls),
                                              as_row(tsheHeightFromDiameterSharmaPartonGnls),
                                              as_row(tsheHeightFromDiameterSharmaPartonBalGnls),
                                              as_row(tsheHeightFromDiameterSharmaZhangGnls),
                                              as_row(tsheHeightFromDiameterSharmaZhangBalGnls),
                                              as_row(tsheHeightFromDiameterWeibullGnls),
                                              as_row(tsheHeightFromDiameterWeibullBalGnls)) %>%
  mutate(responseVariable = "height", species = "TSHE", deltaAic = aic - min(aic)) %>%
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
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, b2p = 0, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5), control = list(maxiter = 200)) # step factor
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, a3 = 0, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5), control = list(maxiter = 400)) # step factor
#tsheDiameterFromHeightChapmanFormBal = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 130, a1p = -55.1, a2 = -1.63, a2p = -1.06, a3 = 0.61, a3p = -0.4, b1 = 0.10, b2 = 0.50, b2p = 0.16), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5), control = list(maxiter = 400), trace = TRUE) # step factor
#tsheDiameterFromHeightChapmanFormBalRelHt = nls(DBH ~ (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare + (a4 + a4p * isPlantation) * relativeHeight) * (exp(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 985, a1p = 262, a2 = -10.6, a2p = -16.6, a3 = 4.3, a3p = 3.0, a4 = -384, a4p = 79, b1 = 0.022, b2 = 1.02, b2p = -0.100), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # plantation effects not significant
#tsheDiameterFromHeightSharmaParton = nls_multstart(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, tshe2016, iter = 100,
#                                                   start_lower = list(a1 = 0.01, a2 = 0.01, b1 = -10, b2 = -1, b3 = -1), 
#                                                   start_upper = list(a1 = 10, a2 = 3, b1 = 10, b2 = 1, b3 = 2), modelweights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
#tsheDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, tshe2016, start = list(a1 = 0.7, a2 = 1, b1 = 0.003, b2 = 0.26, b3 = 0.9), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # step factor
#tsheDiameterFromHeightSharmaParton = nlrob(DBH ~ a1*(TotalHt - 1.37)^(a2 + a2p * isPlantation)*(exp(b1*(tph/topHeight)^(b2 + b2p * isPlantation)*(TotalHt - 1.37)) - 1)^(b3 + b3p * isPlantation), tshe2016, start = list(a1 = 0.7, a2 = 1, a2p = 0, b1 = 0.003, b2 = 0.26, b2p = 0, b3 = 0.9, b3p = 0), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanRichards = nlrob(DBH ~ (a1 + a1p * isPlantation)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -119, a1p = 26.4, b1 = 0.0196, b1p = 0.0004, b2 = 0.847), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanRichardsAat = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * tallerQuasiBasalArea)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -119, a1p = 26.4, a2 = 0, b1 = 0.0196, b1p = 0.0004, b2 = 0.847), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanRichardsPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016physio, start = list(a1 = -189, a1p = 71.1, a2 = -0.017, a3 = -0.339, a4 = -0.129, a5 = 0.427, a6 = -0.008, b1 = 0.0095, b1p = 0.00375, b2 = 0.919), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a2, a4, a5, a6 not significant
tsheDiameterFromHeightChapmanRichardsRelHt = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * relativeHeight)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), tshe2016, start = list(a1 = -322, a1p = 17.7, a2 = -58.4, b1 = 0.0062, b1p = -0.0001, b2 = 0.912), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanForm = nlrob(DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 136, b1 = 0.0100, b2 = 0.924), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # no significant plantation effects
tsheDiameterFromHeightChapmanFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, tshe2016, start = list(a1 = 135, a2 = -0.0045, b1 = 0.010, b2 = 0.924), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanFormBal = nlrob(DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 172, a2 = -1.144, b1 = 0.0138, b2 = 0.893), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a1p not significant, a3 + b2p step size
tsheDiameterFromHeightChapmanFormBalRelHt = nlrob(DBH ~ (a1 + a2 * basalAreaLarger + a3 * standBasalAreaPerHectare + a4 * relativeHeight) * (exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 189, a2 = -3.64, a3 = 2.24, a4 = -46.4, b1 = 0.0137, b2 = 0.847), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightChapmanFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(exp(b1*(TotalHt - 1.37)^b2) - 1), tshe2016, start = list(a1 = 128, a2 = 5.3, b1 = 0.0154, b2 = 0.902), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightLinear = lm(DBH ~ 0 + I(TotalHt - 1.37), tshe2016, weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # isPlantation*(TotalHt - 1.37) not significant(p = 0.036)
tsheDiameterFromHeightMichaelisMentenForm = nlrob(DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), tshe2016, start = list(a1 = 153, a2 = 68, b1 = 0.83), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a1p, a2p, b1p not significant
tsheDiameterFromHeightNaslund = nlrob(DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), tshe2016, start = list(a1 = 3.6, a1p = -0.47, a2 = -0.10, a2p = -0.013), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightParabolic = lm(DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I((TotalHt - 1.37)^2) + I(isPlantation*(TotalHt - 1.37)^2), tshe2016, weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightPower = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, b1 = 1.04), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # no significant plantation effects
tsheDiameterFromHeightPowerAat = nlrob(DBH ~ (a1 + (a2 + a2p * isPlantation) * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 1.70, a2 = -0.00038, a2p = -0.0037, b1 = 1.02, b1p = -0.0047), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a1p not significant
tsheDiameterFromHeightPowerPhysio = nlrob(DBH ~ (a1 + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^b1, tshe2016physio, start = list(a1 = 1.33, a2 = 0.00006, a3 = 0.284, a4 = 0.00044, a5 = -0.00056, a6 = 0.00045, b1 = 1.04), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a2, a4, a5 a6 not significant
tsheDiameterFromHeightPowerRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^b1, tshe2016, start = list(a1 = 1.52, a2 = 0.08, b1 = 1.02), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) 
tsheDiameterFromHeightRuark = nlrob(DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.31, b1 = 0.818, b2 = 0.010), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a1p, b1p, b2p not significant
#tsheDiameterFromHeightSchnute = gsl_nls(DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), tshe2016, start = list(a1 = 0.00108, a2 = 0.058, b1 = 0.96, Ha = 32), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # converges from red alder values but fails to reconverge (singular gradient), NaN-inf or singular gradient with nlrob()
tsheDiameterFromHeightSharmaParton = gsl_nls(DBH ~ a1*(TotalHt - 1.37)^a2*(exp(b1*(tph/topHeight)^b2*(TotalHt - 1.37)) - 1)^b3, tshe2016, start = list(a1 = 0.84, a2 = 1.04, b1 = 0.001, b2 = 0.34, b3 = 0.9), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5), control = gsl_nls_control(maxiter = 125)) # unreliable convergence and NaN-inf with gsl_nls() due to parameter evaporation and collapse to power form, NaN-inf or singular gradient with nlrob(), singular gradient with nls() even with nls_multstart() parameters
tsheDiameterFromHeightSibbesenForm = nlrob(DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, b1 = 0.750, b2 = 0.057), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # no significant plantation effects
tsheDiameterFromHeightSibbesenFormAat = nlrob(DBH ~ (a1 + a2 * tallerQuasiBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.33, a2 = -0.00001, b1 = 0.748, b2 = 0.0578), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # no significant plantation effects
tsheDiameterFromHeightSibbesenFormPhysio = nlrob(DBH ~ (a1 + a1p * isPlantation + a2 * elevation + a3 * sin(3.14159/180 * slope) + a4 * cos(3.14159/180 * aspect) + a5 * sin(3.14159/180 * aspect) + a6 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^(b2 + b2p * isPlantation)), tshe2016physio, start = list(a1 = 2.115, a1p = -0.100, a2 = 0.00017, a3 = 0.486, a4 = 0.00326, a5 = -0.0051, a6 = -0.0001, b1 = 0.736, b2 = 0.0593, b2p = 0.0028), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a2, a4, a5, a6 not significant
tsheDiameterFromHeightSibbesenFormRelHt = nlrob(DBH ~ (a1 + a2 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), tshe2016, start = list(a1 = 2.32, a2 = 0.084, b1 = 0.75, b2 = 0.056), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5))
tsheDiameterFromHeightWeibull = nlrob(DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, tshe2016, start = list(a1 = -225, b1 = 0.011, b2 = 0.88), weights = pmin(TotalHt^if_else(isPlantation, -1.9, -1.8), 0.5)) # a1p, b1p, b2p not significant
#confint2(tsheDiameterFromHeightPowerPhysio, level = 0.99)

tsheDiameterFromHeightChapmanRichards = get_dbh_error("Chapman-Richards", tsheDiameterFromHeightChapmanRichards, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanRichardsAat = get_dbh_error("Chapman-Richards AAT", tsheDiameterFromHeightChapmanRichardsAat, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanRichardsPhysio = get_dbh_error("Chapman-Richards physio", tsheDiameterFromHeightChapmanRichardsPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightChapmanRichardsRelHt = get_dbh_error("Chapman-Richards RelHt", tsheDiameterFromHeightChapmanRichardsRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanForm = get_dbh_error("Chapman-Richards form", tsheDiameterFromHeightChapmanForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormAat = get_dbh_error("Chapman-Richards form AAT", tsheDiameterFromHeightChapmanFormAat, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormBal = get_dbh_error("Chapman-Richards form BAL", tsheDiameterFromHeightChapmanFormBal, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormBalRelHt = get_dbh_error("Chapman-Richards form BAL RelHt", tsheDiameterFromHeightChapmanFormBalRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightChapmanFormRelHt = get_dbh_error("Chapman-Richards form RelHt", tsheDiameterFromHeightChapmanFormRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightLinear = get_dbh_error("linear", tsheDiameterFromHeightLinear, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightMichaelisMentenForm = get_dbh_error("Michaelis-Menten form", tsheDiameterFromHeightMichaelisMentenForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightNaslund = get_dbh_error("Näslund", tsheDiameterFromHeightNaslund, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightParabolic = get_dbh_error("parabolic", tsheDiameterFromHeightParabolic, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPower = get_dbh_error("power", tsheDiameterFromHeightPower, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPowerAat = get_dbh_error("power AAT", tsheDiameterFromHeightPowerAat, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightPowerPhysio = get_dbh_error("power physio", tsheDiameterFromHeightPowerPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightPowerRelHt = get_dbh_error("power RelHt", tsheDiameterFromHeightPowerRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightRuark = get_dbh_error("Ruark", tsheDiameterFromHeightRuark, tshe2016, tshe2016natural, tshe2016plantation)
#tsheDiameterFromHeightSchnute = get_dbh_error("Schnute", tsheDiameterFromHeightSchnute, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSharmaParton = get_dbh_error("modified Sharma-Parton", tsheDiameterFromHeightSharmaParton, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenForm = get_dbh_error("Sibbesen form", tsheDiameterFromHeightSibbesenForm, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenFormAat = get_dbh_error("Sibbesen form AAT", tsheDiameterFromHeightSibbesenFormAat, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightSibbesenFormPhysio = get_dbh_error("Sibbesen form physio", tsheDiameterFromHeightSibbesenFormPhysio, tshe2016physio, tshe2016natural, tshe2016plantationPhysio)
tsheDiameterFromHeightSibbesenFormRelHt = get_dbh_error("Sibbesen form RelHt", tsheDiameterFromHeightSibbesenFormRelHt, tshe2016, tshe2016natural, tshe2016plantation)
tsheDiameterFromHeightWeibull = get_dbh_error("Weibull", tsheDiameterFromHeightWeibull, tshe2016, tshe2016natural, tshe2016plantation)

tsheDiameterFromHeightResults = bind_rows(as_row(tsheDiameterFromHeightChapmanRichards),
                                          as_row(tsheDiameterFromHeightChapmanRichardsAat),
                                          as_row(tsheDiameterFromHeightChapmanRichardsPhysio),
                                          as_row(tsheDiameterFromHeightChapmanRichardsRelHt),
                                          as_row(tsheDiameterFromHeightChapmanForm),
                                          as_row(tsheDiameterFromHeightChapmanFormAat),
                                          as_row(tsheDiameterFromHeightChapmanFormBal),
                                          as_row(tsheDiameterFromHeightChapmanFormBalRelHt),
                                          as_row(tsheDiameterFromHeightChapmanFormRelHt),
                                          as_row(tsheDiameterFromHeightLinear),
                                          as_row(tsheDiameterFromHeightMichaelisMentenForm),
                                          as_row(tsheDiameterFromHeightNaslund),
                                          as_row(tsheDiameterFromHeightParabolic),
                                          as_row(tsheDiameterFromHeightPower),
                                          as_row(tsheDiameterFromHeightPowerAat),
                                          as_row(tsheDiameterFromHeightPowerPhysio),
                                          as_row(tsheDiameterFromHeightPowerRelHt),
                                          as_row(tsheDiameterFromHeightSharmaParton, significant = FALSE),
                                          as_row(tsheDiameterFromHeightRuark),
                                          as_row(name = "Schnute"),
                                          as_row(tsheDiameterFromHeightSibbesenForm),
                                          as_row(tsheDiameterFromHeightSibbesenFormAat),
                                          as_row(tsheDiameterFromHeightSibbesenFormPhysio),
                                          as_row(tsheDiameterFromHeightSibbesenFormRelHt),
                                          as_row(tsheDiameterFromHeightWeibull)) %>% 
  mutate(responseVariable = "DBH", species = "TSHE", deltaAic = aic - min(aic, na.rm = TRUE)) %>%
  arrange(desc(deltaAic))
print(tsheDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)

ggplot(tshe2016) +
  geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
  #geom_line(aes(x = tsheDiameterFromHeightSharmaParton$fitted.values, y = TotalHt, color = "adapted Sharma-Parton", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanFormBal$fitted.values, y = TotalHt, color = "Chapman-Richards form BAL", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanFormAat$fitted.values, y = TotalHt, color = "Chapman-Richards form AAT", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanRichards$fitted.values, y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightChapmanForm$fitted.values, y = TotalHt, color = "Chapman-Richards form", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightMichaelisMentenForm$fitted.values, y = TotalHt, color = "Michaelis-Menten form", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightNaslund$fitted.values, y = TotalHt, color = "Näslund", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightPower$fitted.values, y = TotalHt, color = "power", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightRuark$fitted.values, y = TotalHt, color = "Ruark", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightSibbesenForm$fitted.values, y = TotalHt, color = "Sibbesen", group = isPlantation)) +
  #geom_line(aes(x = tsheDiameterFromHeightWeibull$fitted.values, y = TotalHt, color = "Weibull", group = isPlantation)) +
  #geom_line(aes(x = 1*topHeight^1*(1 - exp(-0.01 * (tph/standBasalAreaPerHectare)^1*(TotalHt - 1.37)))^1, y = TotalHt, color = "Sharma-Parton"), alpha = 0.5) +
  #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
  #geom_line(aes(x = -70 * log(1 - pmin(0.01*(TotalHt - 1.37)^1.1, 0.999)), y = TotalHt, color = "Chapman-Richards"), na.rm = TRUE) +
  ##geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards form", group = isPlantation), alpha = 0.5) +
  #geom_line(aes(x = (1.75 + 0.000001 * tallerQuasiBasalArea + -0.000001 * standQuasiBasalArea) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards form AAT", group = isPlantation), alpha = 0.5) +
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


## collect model parameters
tsheParameters = bind_rows(bind_rows(get_coefficients(tsheHeightFromDiameterChapmanRichards),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsBal),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsBalPhysio),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsBalRelHt),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsPhysio),
                                     get_coefficients(tsheHeightFromDiameterCurtis),
                                     get_coefficients(tsheHeightFromDiameterHossfeld),
                                     get_coefficients(tsheHeightFromDiameterKorf),
                                     get_coefficients(tsheHeightFromDiameterLinear),
                                     get_coefficients(tsheHeightFromDiameterMichaelisMenten),
                                     get_coefficients(tsheHeightFromDiameterParabolic),
                                     get_coefficients(tsheHeightFromDiameterPower),
                                     get_coefficients(tsheHeightFromDiameterProdan),
                                     get_coefficients(tsheHeightFromDiameterRatkowsky),
                                     get_coefficients(tsheHeightFromDiameterRichards),
                                     get_coefficients(tsheHeightFromDiameterSharmaParton),
                                     get_coefficients(tsheHeightFromDiameterSharmaPartonBal),
                                     get_coefficients(tsheHeightFromDiameterSharmaPartonBalPhysio),
                                     get_coefficients(tsheHeightFromDiameterSharmaPartonPhysio),
                                     get_coefficients(tsheHeightFromDiameterSharmaZhang),
                                     get_coefficients(tsheHeightFromDiameterSharmaZhangBal),
                                     get_coefficients(tsheHeightFromDiameterSibbesen),
                                     get_coefficients(tsheHeightFromDiameterWeibull),
                                     get_coefficients(tsheHeightFromDiameterWeibullBal),
                                     get_coefficients(tsheHeightFromDiameterWeibullBalRelHt),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsGnls),
                                     get_coefficients(tsheHeightFromDiameterChapmanRichardsBalGnls),
                                     get_coefficients(tsheHeightFromDiameterSharmaPartonGnls),
                                     get_coefficients(tsheHeightFromDiameterSharmaPartonBalGnls),
                                     get_coefficients(tsheHeightFromDiameterSharmaZhangGnls),
                                     get_coefficients(tsheHeightFromDiameterSharmaZhangBalGnls),
                                     get_coefficients(tsheHeightFromDiameterWeibullGnls),
                                     get_coefficients(tsheHeightFromDiameterWeibullBalGnls)) %>%
                             mutate(responseVariable = "height"),
                           bind_rows(get_coefficients(tsheDiameterFromHeightChapmanRichards),
                                     get_coefficients(tsheDiameterFromHeightChapmanRichardsAat),
                                     get_coefficients(tsheDiameterFromHeightChapmanRichardsPhysio),
                                     get_coefficients(tsheDiameterFromHeightChapmanRichardsRelHt),
                                     get_coefficients(tsheDiameterFromHeightChapmanForm),
                                     get_coefficients(tsheDiameterFromHeightChapmanFormAat),
                                     get_coefficients(tsheDiameterFromHeightChapmanFormBal),
                                     get_coefficients(tsheDiameterFromHeightChapmanFormBalRelHt),
                                     get_coefficients(tsheDiameterFromHeightChapmanFormRelHt),
                                     get_coefficients(tsheDiameterFromHeightLinear),
                                     get_coefficients(tsheDiameterFromHeightMichaelisMentenForm),
                                     get_coefficients(tsheDiameterFromHeightNaslund),
                                     get_coefficients(tsheDiameterFromHeightParabolic),
                                     get_coefficients(tsheDiameterFromHeightPower),
                                     get_coefficients(tsheDiameterFromHeightPowerAat),
                                     get_coefficients(tsheDiameterFromHeightPowerPhysio),
                                     get_coefficients(tsheDiameterFromHeightPowerRelHt),
                                     get_coefficients(tsheDiameterFromHeightRuark),
                                     #get_coefficients(tsheDiameterFromHeightSchnute),
                                     get_coefficients(tsheDiameterFromHeightSharmaParton),
                                     get_coefficients(tsheDiameterFromHeightSibbesenForm),
                                     get_coefficients(tsheDiameterFromHeightSibbesenFormAat),
                                     get_coefficients(tsheDiameterFromHeightSibbesenFormPhysio),
                                     get_coefficients(tsheDiameterFromHeightSibbesenFormRelHt),
                                     get_coefficients(tsheDiameterFromHeightWeibull)) %>%
                             mutate(responseVariable = "DBH")) %>%
  mutate(species = "TSHE",
         a1 = as.numeric(a1), a1p = as.numeric(a1p), a2 = as.numeric(a2), a2p = as.numeric(a2p), a3 = as.numeric(a3), a3p = as.numeric(a3p),
         a4 = as.numeric(a4), a4p = as.numeric(a4p), a5 = as.numeric(a5), a6 = as.numeric(a6), 
         b1 = as.numeric(b1), b1p = as.numeric(b1p), b2 = as.numeric(b2), b2p = as.numeric(b2p), b3 = as.numeric(b3))
if ("b3p" %in% names(tsheParameters))
{
  tsheParameters$b3p = as.numeric(tsheParameters$b3p)
} else
{
  tsheParameters$b3p = NA_real_
}
tsheParameters %<>% relocate(responseVariable, species, name, a1, a1p, a2, a2p, a3, a3p, a4, a4p, a5, a6, b1, b1p, b2, b2p, b3, b3p)


## basal area from height
#tsheBasalAreaFromHeightKorf = nlrob(basalArea ~ a1*exp(b1*(imputedHeight - 1.37)^b2) - 1, psme2016, start = list(a1 = 1, b1 = 0.00009, b2 = 2.2), weights = pmin(1/basalArea, 1E4))
tsheBasalAreaFromHeightKorf = nlrob(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^(b2 + b2p * isPlantation)) - 1), tshe2016, start = list(a1 = 5.5, b1 = 0.00003, b2 = 2.09, b2p = -0.060), weights = pmin(1/basalArea, 1E4)) # a1p, b1p not significant, step factor on b1p without a1p
tsheBasalAreaFromHeightPower = nlrob(basalArea ~ (a1 + a1p*isPlantation)*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), tshe2016, start = list(a1 = 11/7 * 0.25 * pi * 0.01^2, a1p = -0.00005, b1 = 2.18, b1p = -0.178), weights = pmin(1/basalArea, 1E4))
#confint_nlrob(tsheBasalAreaFromHeightKorf, level = 0.99, weights = pmin(1/tshe2016$basalArea, 1E4))

tsheBasalAreaFromHeightKorf$fitted.values = predict(tsheBasalAreaFromHeightKorf, tshe2016)
tsheBasalAreaFromHeightKorf$residuals = tsheBasalAreaFromHeightKorf$fitted.values - tshe2016$basalArea
tsheBasalAreaFromHeightPower$fitted.values = predict(tsheBasalAreaFromHeightPower, tshe2016)
tsheBasalAreaFromHeightPower$residuals = tsheBasalAreaFromHeightPower$fitted.values - tshe2016$basalArea

tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
        "Korf", AIC(tsheBasalAreaFromHeightKorf), 100^2 * mean(tsheBasalAreaFromHeightKorf$residuals), mean(abs(tsheBasalAreaFromHeightKorf$residuals)), 1 - sum(tsheBasalAreaFromHeightKorf$residuals^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2)),
        "power", AIC(tsheBasalAreaFromHeightPower), 100^2 * mean(tsheBasalAreaFromHeightPower$residuals), mean(abs(tsheBasalAreaFromHeightPower$residuals)), 1 - sum(tsheBasalAreaFromHeightPower$residuals^2) / sum((tshe2016$basalArea - mean(tshe2016$basalArea)^2))) %>%
  mutate(deltaAIC = aic - min(aic)) %>%
  arrange(desc(deltaAIC))

ggplot(tshe2016) +
  geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = imputedHeight, y = tsheBasalAreaFromHeightKorf$fitted.values, color = "Korf", group = isPlantation)) +
  geom_line(aes(x = imputedHeight, y = tsheBasalAreaFromHeightPower$fitted.values, color = "power", group = isPlantation)) +
  #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
  labs(x = "western hemlock height, m", y = "basal area, m²", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
