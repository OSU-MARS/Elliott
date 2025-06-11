# load libraries, functions, and trees2016 from Elliott Stand Data Feb2022.R


## cascara buckthorn height-diameter regression form sweep

#rhpu HeightFromDiameter$gamPhysio = gam(TotalHt ~ s(DBH, elevation, slope, sin(3.14159/180 * aspect), cos(3.14159/180 * aspect), topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 85, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, select = TRUE, weights = dbhWeight) 

# bs= "ts" -> 367, gamma = 2 -> 367, k = 169 min vs 367 default, method = "REML" -> 367

#rhpu HeightFromDiameter$sharmaPartonBalPhysio = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^(a2 + a2p * isPlantation) * (1 + a3 * elevation + a4 * sin(3.14159/180 * aspect) + a5 * cos(3.14159/180 * aspect) + a6 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), rhpu2016, start = list(a1 = 39.8, a1p = -12.3, a2 = 0.52, a2p = 0.0027, a3 = 0.00001, a4 = 0.0131, a5 = 0.0046, a6 = 0.0060, b1 = -0.0098, b1p = -0.0143, b2 = 0.125, b2p = -0.186, b3 = 1.12, b3p = 0.0086), weights = rhpuHeightFromDiameterWeights)

rhpu2016 = trees2016 %>%
  filter(Species == "RHPU", isLiveUnbroken, is.na(TotalHt) == FALSE,is.na(elevation)==FALSE) %>% # live cascara buckthorns measured for height
  mutate(dbhWeight = pmin(TreeCount/(0.14*DBH^1.20), 5*TreeCount),
         heightWeight = pmin(TreeCount/(2.29*(TotalHt - 1.37)^1.45), 5*TreeCount))

# no trees without physiographic variables
rhpu2016gamConstraint = c(DBH = -1.2264/0.5099, TotalHt = 1.37, standBasalAreaPerHectare = median(rhpu2016$standBasalAreaPerHectare), basalAreaLarger = median(rhpu2016$basalAreaLarger), standBasalAreaApprox = median(rhpu2016$standBasalAreaApprox), tallerApproxBasalArea = median(rhpu2016$tallerApproxBasalArea), elevation = median(rhpu2016$elevation), slope = median(rhpu2016$slope), aspect = median(rhpu2016$aspect), topographicShelterIndex = median(rhpu2016$topographicShelterIndex), relativeHeight = median(rhpu2016$relativeHeight), relativeDiameter = median(rhpu2016$relativeDiameter)) # point constraint for mgcv::s()

rhpu2016defaultWeight = rhpu2016 %>% mutate(dbhWeight = pmin(TreeCount/DBH, 5*TreeCount),
                                            heightWeight = pmin(TreeCount/TotalHt, 5*TreeCount))
rhpu2016defaultWeightPhysio = rhpu2016defaultWeight %>% filter(is.na(elevation) == FALSE)

# rhpuOptions = tibble(fitHeight = TRUE, 
#                      fitHeightNlrob = FALSE,
#                      fitHeightGnls = FALSE,
#                      fitHeightMixed = FALSE,
#                      fitDbh = TRUE,
#                      fitDbhNlrob = FALSE,
#                      fitDbhMixed = FALSE)

#make a tibble to store the parameters and later call them into the code
rhpuOptions = tibble(fitHeight = TRUE, #non-linear least square, height as response
                     fitHeightNlrob = FALSE, #robust non-linear least square, height as response
                     fitHeightGnls = TRUE, #generalized least square, height as response
                     fitHeightMixed = FALSE, #non-linear mixed effects, height as response
                     fitDbh = TRUE, #non-linear least square, dbh as response
                     fitDbhNlrob = FALSE, #robust non-linear least square, dbh as response
                     fitDbhMixed = FALSE, #non-linear mixed effects, dbh as response
                     includeInvestigatory = TRUE #added investigatory plots, and figures of the results
)

if (rhpuOptions$fitHeight) { #if the value in the column fitHeight of rhpuOptions table is TRUE execute the expression within the curly braces.
  rhpuHeightFromDiameter = list(linear = fit_lm("linear", TotalHt ~ 0 + DBH, rhpu2016)) # isPlantation*DBH not significant (p = 0.044) #creates output from the model fitting and validation 10*10=100 rows of all combination of folds and repetition and stores it as a list and adds all the following models in a similar fashion to the existing list (notice the $ sign in the code from the second line withing the curly braces)
  rhpuHeightFromDiameter$parabolic = fit_lm("parabolic", TotalHt ~ 0 + DBH + I(DBH^2), rhpu2016) # isPlantation*DBH not quite significant (p = 0.106), isPlantation*DBH^2 not significant
  
  rhpuHeightFromDiameter$chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131)) # a1p, b1p, b2p not significant
  rhpuHeightFromDiameter$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1)) # a3, a3p, b1p, b2p not significant
  rhpuHeightFromDiameter$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), significant = FALSE) # a2, a3, a4, a5, a6, a7, a8p, b2p not significant
  #rhpuHeightFromDiameter$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, start = list(a1 = 58, a1p = -16, a2 = 0, a2p = 0.4, a8 = 0.3, a10 = -1.3, b1 = -0.012, b1p = -0.003, b2 = 1.13), significant = FALSE) # a2, a10, a10p not significant
  rhpuHeightFromDiameter$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 61, a1p = -9, a2 = -0.1, a2p = 0.6, a10 = -1.3, b1 = -0.012, b2 = 1.1), significant = FALSE) # a2, a10, a10p not significant
  #rhpuHeightFromDiameter$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), rhpu2016, start = list(a1 = 7, a1p = 5, a2 = 0.2, a2p = 0.24, a9 = 47, a9p = -27, b1 = -0.021, b2 = 0.8, b2p = 0.2)) # a2, a3, a3p, b1p not significant, job step factor with nlrob()
  rhpuHeightFromDiameter$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15)) # a4, a5, a6, a7, a8p, b2p not significant
  rhpuHeightFromDiameter$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 72, a10 = -3.2, b1 = -0.012, b2 = 1.09)) # a10p not significant
  #rhpuHeightFromDiameter$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), significant = FALSE) # a10, a10p not significant
  rhpuHeightFromDiameter$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, rhpu2016, start = list(a1 = 0.560, b1 = 0.069)) # a1p, b1p not significant
  rhpuHeightFromDiameter$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176)) # b2p not significant
  rhpuHeightFromDiameter$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), rhpu2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175)) # a1p, b1p, b2p not significant
  rhpuHeightFromDiameter$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176)) # b1p not significant
  rhpuHeightFromDiameter$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), rhpu2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649)) # a2p, a3p not significant
  rhpuHeightFromDiameter$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^b1, rhpu2016, start = list(a1 = 0.542, b1 = 0.939)) # a1p, b1p not significant
  rhpuHeightFromDiameter$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), rhpu2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
  #rhpuHeightFromDiameter$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), rhpu2016, start = list(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008)) # dp not significant, susceptible to NaN-inf
  #rhpuHeightFromDiameter$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09)) # a1p, b3p, b4p not significant
  rhpuHeightFromDiameter$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 38, b1 = 0.1, b2 = -0.013, b3 = -0.1, b4 = 1.03)) # a1p, b1p, b2p, b3p, b4p not significant
  #rhpuHeightFromDiameter$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10)) # b1, b1p, a4, a5, a6, a7, b3p, b4p not significant
  #rhpuHeightFromDiameter$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 25, a1p = -6, a8 = 0.12, a10 = -0.7, b1 = 0.21, b2 = -0.008, b2p = -0.011, b3 = -0.01, b4 = 1.12), significant = FALSE) # a10, a10p not significant
  rhpuHeightFromDiameter$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 39, a10 = -1.7, b1 = 0.12, b2 = -0.01, b3 = 0, b4 = 1.07), significant = FALSE) # a10, a10p not significant
  #rhpuHeightFromDiameter$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09)) # a4, a5, a5, a6, a7, b1p, b3p, b4p not significant
  #rhpuHeightFromDiameter$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, start = list(a1 = 21, a10 = 0, b1 = 0.25, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = 0, b4 = 1.12), significant = FALSE) # a10, a10p not significant
  #rhpuHeightFromDiameter$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, start = list(a1 = 36, a8 = 0.18, a10 = -2, b1 = 0.13, b2 = -0.01, b3 = -0.03, b4 = 1.09), significant = FALSE) # a1p, a10, a10p, b2p not significant
  #rhpuHeightFromDiameter$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), rhpu2016, start = list(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097)) # b1, b1p, b2p, b3p not significant
  rhpuHeightFromDiameter$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, rhpu2016, start = list(a1 = 45, a1p = -7, a2 = -0.1, a2p = 0.4, b1 = -0.05, b2 = -0.02, b3 = -0.078, b4 = 1.08)) # a2, b1, b1p, b3, b3p, b4p not significant
  rhpuHeightFromDiameter$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), rhpu2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078)) # a1p, b1p, b2p not significant
  rhpuHeightFromDiameter$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), rhpu2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141)) # b2p not significant
  rhpuHeightFromDiameter$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), rhpu2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131)) # a1p, a2, a3, b1p, b2p not significant
  rhpuHeightFromDiameter$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), rhpu2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778)) # a1p, a2, a3, a3p, b1p, b2p not significant
  
  if (rhpuOptions$fitHeightNlrob)
  {
    rhpuHeightFromDiameterNlrob = list(chapmanRichards = fit_nlrob("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 5.2, b1 = -0.015, b2 = 1.131)))
    #rhpuHeightFromDiameterNlrob$chapmanRichardsBal = fit_nlrob("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1))
    #rhpuHeightFromDiameterNlrob$chapmanRichardsBalPhysio = fit_nlrob("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # job step factor
    #rhpuHeightFromDiameterNlrob$chapmanRichardsBalPhysioRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, start = list(a1 = 58, a1p = -14, a2 = -0.11, a2p = 0.5, a8 = 0.3, a10 = -1.8, b1 = -0.012, b1p = -0.003, b2 = 1.14), significant = FALSE)
    #rhpuHeightFromDiameterNlrob$chapmanRichardsBalRelDbh = fit_nlrob("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 61, a1p = -8, a2 = -0.13, a2p = 0.6, a10 = -1.4, b1 = -0.012, b2 = 1.12), significant = FALSE)
    rhpuHeightFromDiameterNlrob$chapmanRichardsBalRelHt = fit_nlrob("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), rhpu2016, start = list(a1 = 0, a1p = 17, a2 = 0, a2p = 0.25, a3 = 0.02, a9 = 38, a9p = -28, b1 = -0.023, b2 = 0.4, b2p = 0.9), control = nls.control(tol = 0.01)) # job step factor
    #rhpuHeightFromDiameterNlrob$chapmanRichardsPhysio = fit_nlrob("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15))
    rhpuHeightFromDiameterNlrob$chapmanRichardsRelDbh = fit_nlrob("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 66, a10 = -3.2, b1 = -0.010, b2 = 1.09))
    #rhpuHeightFromDiameterNlrob$chapmanRichardsRelDbhPhysio = fit_nlrob("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), control = nls.control(tol = 1E-4), significant = FALSE)
    rhpuHeightFromDiameterNlrob$curtis = fit_nlrob("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, rhpu2016, start = list(a1 = 0.560, b1 = 0.069))
    rhpuHeightFromDiameterNlrob$hossfeld = fit_nlrob("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176))
    rhpuHeightFromDiameterNlrob$korf = fit_nlrob("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), rhpu2016, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175))
    #rhpuHeightFromDiameterNlrob$michaelisMenten = fit_nlrob("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176))
    rhpuHeightFromDiameterNlrob$prodan = fit_nlrob("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), rhpu2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649))
    rhpuHeightFromDiameterNlrob$power = fit_nlrob("power", TotalHt ~ 1.37 + a1*DBH^b1, rhpu2016, start = list(a1 = 0.542, b1 = 0.939))
    rhpuHeightFromDiameterNlrob$ratkowsky = fit_nlrob("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), rhpu2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
    #rhpuHeightFromDiameterNlrob$richardsW = fit_nlrob("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), rhpu2016, start = list(Ha = 43, Hap = -10, d = 0.9, kU = 0.012, kUp = 0.004), control = nls.control(tol = 0.001)) # job step factor
    #rhpuHeightFromDiameterNlrob$sharmaParton = fit_nlrob("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09), control = nls.control(tol = 0.001)) # job step factor
    rhpuHeightFromDiameterNlrob$sharmaPartonBal = fit_nlrob("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 44, b1 = 0.07, b2 = -0.013, b3 = -0.10, b4 = 1.03), control = nls.control(maxiter = 100, tol = 0.001)) # job step factor
    #rhpuHeightFromDiameterNlrob$sharmaPartonBalPhysio = fit_nlrob("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10), control = nls.control(tol = 1E-4)) # job step factor
    #rhpuHeightFromDiameterNlrob$sharmaPartonBalPhysioRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 25, a1p = -8, a8 = 0.13, a10 = -0.9, b1 = 0.18, b2 = -0.011, b2p = -0.010, b3 = 0, b4 = 1.13), control = nls.control(tol = 0.001), significant = FALSE) # job step factor
    #rhpuHeightFromDiameterNlrob$sharmaPartonBalRelDbh = fit_nlrob("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 50, a10 = -3, b1 = 0.12, b2 = -0.01, b3 = 0, b4 = 1.07), control = nls.control(maxiter = 100, tol = 0.001), significant = FALSE) # step factor
    #rhpuHeightFromDiameterNlrob$sharmaPartonPhysio = fit_nlrob("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, start = list(a1 = 28, a1p = -10, a8 = 0.13, b1 = 0.16, b2 = -0.011, b2p = -0.01, b3 = 0, b4 = 1.1), control = nls.control(tol = 0.01)) # b3 not significant, job step factor
    #rhpuHeightFromDiameterNlrob$sharmaPartonRelDbh = fit_nlrob("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, start = list(a1 = 27, a10 = -0.7, b1 = 0.22, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = -0.003, b4 = 1.12), control = nls.control(tol = 0.001), significant = FALSE)
    #rhpuHeightFromDiameterNlrob$sharmaPartonRelDbhPhysio = fit_nlrob("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, start = list(a1 = 41, a8 = 0.2, a10 = -2, b1 = 0.13, b2 = -0.01, b3 = 0, b4 = 1.09), significant = FALSE)
    #rhpuHeightFromDiameterNlrob$sharmaZhang = fit_nlrob("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), rhpu2016, start = list(a1 = 36, a1p = -3.0, b1 = 0.1, b2 = -0.02, b3 = 0, b4 = 1.2, b4p = -0.2)) # b3 not significant
    #rhpuHeightFromDiameterNlrob$sharmaZhangBal = fit_nlrob("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, rhpu2016, start = list(a1 = 44, a1p = -7, a2 = -0.12, a2p = 0.45, b1 = 0.05, b2 = -0.017, b3 = -0.02, b4 = 1.1), control = nls.control(maxiter = 100, tol = 0.001)) # b3 not significant, job step factor
    rhpuHeightFromDiameterNlrob$sibbesen = fit_nlrob("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), rhpu2016, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078))
    rhpuHeightFromDiameterNlrob$weibull = fit_nlrob("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)* (1 - exp((b1 + b1p * isPlantation)*DBH^b2)), rhpu2016, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141), control = nls.control(maxiter = 100, tol = 1E-4)) # job step factor
    #rhpuHeightFromDiameterNlrob$weibullBal = fit_nlrob("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), rhpu2016, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131))
    rhpuHeightFromDiameterNlrob$weibullBalRelHt = fit_nlrob("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), rhpu2016, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778))
    lapply(rhpuHeightFromDiameterNlrob$sharmaPartonPhysio$fit, confint_nlrob, level = 0.99)
  } else {
    rhpuHeightFromDiameterNlrob = list()
  }
  #fitting models with defaultweight 'rhpu2016defaultWeight' which were fitted earlier without any weights.
  #commented out models without any note at the end did not converge, for others the note is given mentioning the reason of error.
  rhpuHeightFromDiameterGslNlsDefault = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, rhpu2016defaultWeight, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131)))
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsBal = fit_gsl_nls("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, rhpu2016defaultWeight, start = list(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1))
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysio = fit_gsl_nls("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016defaultWeightPhysio, start = list(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12), significant = FALSE)
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelHt = fit_gsl_nls("Chapman-Richards BA+L RelHt", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3 * standBasalAreaPerHectare + (a9 + a9p * isPlantation) * relativeHeight) * (1 - exp(b1*DBH))^(b2 + b2p * isPlantation), rhpu2016defaultWeight, start = list(a1 = 7, a1p = 5, a2 = 0.2, a2p = 0.24, a3 = -0.03, a9 = 47, a9p = -27, b1 = -0.021, b2 = 0.8, b2p = 0.2))
  #rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsBalPhysioRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, start = list(a1 = 58, a1p = -14, a2 = -0.07, a2p = 0.48, a8 = 0.27, a10 = -1.3, b1 = -0.012, b1p = -0.003, b2 = 1.13), significant = FALSE)
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsBalRelDbh = fit_gsl_nls("Chapman-Richards BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 61, a1p = -9, a2 = -0.1, a2p = 0.6, a10 = -1.4, b1 = -0.012, b2 = 1.11), significant = FALSE)
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016defaultWeightPhysio, start = list(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15))
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbh = fit_gsl_nls("Chapman-Richards RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter) * (1 - exp(b1*DBH))^b2, rhpu2016defaultWeight, start = list(a1 = 74, a10 = -3.2, b1 = -0.011, b2 = 1.09))
  rhpuHeightFromDiameterGslNlsDefault$chapmanRichardsRelDbhPhysio = fit_gsl_nls("Chapman-Richards RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016defaultWeightPhysio, start = list(a1 = 63, a1p = -17, a8 = 0.3, a10 = -2.1, b1 = -0.011, b1p = -0.006, b2 = 1.15), significant = FALSE) #produced NA or infinity
  rhpuHeightFromDiameterGslNlsDefault$curtis = fit_gsl_nls("Curtis", TotalHt ~ 1.37 + a1 * DBH / (1 + DBH)^b1, rhpu2016defaultWeight, start = list(a1 = 0.560, b1 = 0.069))
  rhpuHeightFromDiameterGslNlsDefault$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), rhpu2016defaultWeight, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176))
  rhpuHeightFromDiameterGslNlsDefault$korf = fit_gsl_nls("Korf", TotalHt ~ 1.37 + a1*exp(b1*DBH^b2), rhpu2016defaultWeight, start = list(a1 = 1825, b1 = -8.726, b2 = -0.175))
  rhpuHeightFromDiameterGslNlsDefault$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), rhpu2016defaultWeight, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176))
  rhpuHeightFromDiameterGslNlsDefault$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), rhpu2016defaultWeight, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649))
  rhpuHeightFromDiameterGslNlsDefault$power = fit_gsl_nls("power", TotalHt ~ 1.37 + a1*DBH^b1, rhpu2016defaultWeight, start = list(a1 = 0.542, b1 = 0.939))
  rhpuHeightFromDiameterGslNlsDefault$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), rhpu2016defaultWeight, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151))
  #rhpuHeightFromDiameterGslNlsDefault$richardsW = fit_gsl_nls("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation) * (1 + ((1.37/(Ha + Hap*isPlantation))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), rhpu2016defaultWeight, start = list(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008)) #produced NA or infinity
  #rhpuHeightFromDiameterGslNlsDefault$sharmaParton = fit_gsl_nls("Sharma-Parton", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016defaultWeight, start = list(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09))#produced NA or infinity
  rhpuHeightFromDiameterGslNlsDefault$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + a1*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016defaultWeight, start = list(a1 = 38, b1 = 0.12, b2 = -0.013, b3 = -0.1, b4 = 1.02))
  #rhpuHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016defaultWeightPhysio, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10))
  rhpuHeightFromDiameterGslNlsDefault$sharmaPartonBalPhysioRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016defaultWeightPhysio, start = list(a1 = 23, a1p = -6, a8 = 0.12, a10 = -0.7, b1 = 0.21, b2 = -0.01, b2p = -0.010, b3 = -0.012, b4 = 1.14), significant = FALSE)#produced NA or infinity
  rhpuHeightFromDiameterGslNlsDefault$sharmaPartonBalRelDbh = fit_gsl_nls("Sharma-Parton BA+L RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016defaultWeight, start = list(a1 = 35, a10 = -1.4, b1 = 0.15, b2 = -0.01, b3 = 0, b4 = 1.07), significant = FALSE)#produced NA or infinity
  #rhpuHeightFromDiameterGslNlsDefault$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016defaultWeightPhysio, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09))#produced NA or infinity
  #rhpuHeightFromDiameterGslNlsDefault$sharmaPartonRelDbh = fit_gsl_nls("Sharma-Parton RelDbh", TotalHt ~ 1.37 + (a1 + a10 * relativeDiameter)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016defaultWeight, start = list(a1 = 19, a10 = -0.3, b1 = 0.29, b1p = -0.09, b2 = -0.013, b2p = -0.011, b3 = -0.03, b4 = 1.13), significant = FALSE)#produced NA or infinity
  rhpuHeightFromDiameterGslNlsDefault$sharmaPartonRelDbhPhysio = fit_gsl_nls("Sharma-Parton RelDbh physio", TotalHt ~ 1.37 + (a1 + a8 * topographicShelterIndex + a10 * relativeDiameter)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016defaultWeightPhysio, start = list(a1 = 36, a8 = 0.18, a10 = 0, b1 = 0.2, b2 = -0.01, b3 = 0.03, b4 = 1.09), significant = FALSE)
  #rhpuHeightFromDiameterGslNlsDefault$sharmaZhang = fit_gsl_nls("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), rhpu2016defaultWeight, start = list(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097))
  rhpuHeightFromDiameterGslNlsDefault$sharmaZhangBal = fit_gsl_nls("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, rhpu2016defaultWeight, start = list(a1 = 53.2, a1p = -8.857, a2 = -0.002, a2p = 0.10, b1 = -0.016, b2 = -0.025, b3 = -0.078, b4 = 1.126))
  rhpuHeightFromDiameterGslNlsDefault$sibbesen = fit_gsl_nls("Sibbesen", TotalHt ~ 1.37 + a1*DBH^(b1*DBH^b2), rhpu2016defaultWeight, start = list(a1 = 0.302, b1 = 1.495, b2 = -0.078))
  rhpuHeightFromDiameterGslNlsDefault$weibull = fit_gsl_nls("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), rhpu2016defaultWeight, start = list(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141))
  rhpuHeightFromDiameterGslNlsDefault$weibullBal = fit_gsl_nls("Weibull BA+L", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), rhpu2016defaultWeight, start = list(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131))
  rhpuHeightFromDiameterGslNlsDefault$weibullBalRelHt = fit_gsl_nls("Weibull BA+L RelHt", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a9 + a9p * isPlantation) * pmin(relativeHeight, 1.5)) * (1 - exp(b1*DBH^b2)), rhpu2016defaultWeight, start = list(a1 = 18.9, a2 = 0.171, a2p = 0.166, a9 = 46.6, a9p = -9.98, b1 = -0.019, b2 = 0.778))

    rhpuHeightFromDiameter$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # newton() step failure with family = scat, internal code errors with scat(theta = <fixed val>), see https://stats.stackexchange.com/questions/410515/how-different-are-restricted-cubic-splines-and-penalized-splines for discusson of thin plate versus other spline types
    rhpuHeightFromDiameter$gamBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
    rhpuHeightFromDiameter$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 20, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # slope and elevation not supported, aspect not tested since insufficient data for full model
    rhpuHeightFromDiameter$gamBalPhysioRelDbh = fit_gam("REML GAM BA+L RelDbh physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 57, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
    rhpuHeightFromDiameter$gamBalRelDbh = fit_gam("REML GAM BA+L RelDbh", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 22, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
    rhpuHeightFromDiameter$gamPhysio = fit_gam("REML GAM physio", TotalHt ~ s(DBH, elevation, slope, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 18, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # k reduces from 85 to 18 without aspect
    rhpuHeightFromDiameter$gamRelDbh = fit_gam("REML GAM RelDbh", TotalHt ~ s(DBH, relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
    rhpuHeightFromDiameter$gamRelDbhPhysio = fit_gam("REML GAM RelDbh physio", TotalHt ~ s(DBH, elevation, slope, topographicShelterIndex, relativeDiameter, bs = "ts", k = 57, by = as.factor(isPlantation), pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)

  save(file = "data/rhpu TotalHt.Rdata", rhpuHeightFromDiameter, rhpuHeightFromDiameterNlrob, rhpuHeightFromDiameterGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory)
{
  #print(rhpuHeightFromDiameterResults %>% select(-responseVariable, -species, -fixedWeight, -n, -power, -significant, -contains("NaturalRegen"), -contains("Plantation")), n = 30)
  ggplot() +
    geom_point(aes(x = rhpu2016$DBH, y = rhpu2016$TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$sharmaZhang), color = "Sharma-Zhang", group = rhpu2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$sharmaParton), color = "Sharma-Parton", group = rhpu2016$isPlantation), alpha = 0.5) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$chapmanRichards), color = "Chapman-Richards", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$curtis), color = "Curtis", group = rhpu2016$isPlantation)) +
    geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$gam), color = "GAM", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$korf), color = "Korf", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$linear), color = "linear", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$michaelisMenten), color = "Michaelis-Menten", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$parabolic), color = "parabolic", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$power), color = "power", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$prodan), color = "Prodan", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$ratkowsky), color = "Ratkowsky", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$richardsW), color = "unified Richards", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$sibbesen), color = "Sibbesen", group = rhpu2016$isPlantation)) +
    #geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$weibull), color = "Weibull", group = rhpu2016$isPlantation)) +
    annotate("text", x = 0, y = 65, label = "cascara buckthorn, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(ylim = c(0, 65)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
  
  # dbhClassSize = 50
  # errorByDbhClass = tibble(dbhClass = dbhClassSize*floor(rhpu2016$DBH/dbhClassSize) + 0.5*dbhClassSize, fittedValue = predict(rhpuHeightFromDiameter$gam, rhpu2016), height = rhpu2016$TotalHt, residual = fittedValue - height) %>%
  #  #mutate(residual = residual - if_else(dbhClass == 50, -0.477/376, 0.107/95)) %>%
  #  group_by(dbhClass) %>%
  #  summarize(n = n(),
  #            totalHeight = sum(height),
  #            totalFitted = sum(fittedValue),
  #            meanBiasPerTree = sum(residual) / n,
  #            meanBiasPerTreePct = 100 * sum(residual/height) / n,
  #            minError = min(residual),
  #            meanError = mean(residual),
  #            maxError = max(residual),
  #            minPct = 100 * min(residual/height),
  #            meanPct = 100 * mean(residual/height),
  #            maxPct = 100 * max(residual/height),
  #            .groups = "drop") %>%
  #  filter(n >= 10)
  # errorByDbhClass
}


## cascara buckthorn height-diameter GNLS regressions
if (rhpuOptions$fitHeightGnls)
{
  rhpuHeightFromDiameterGnls = list(chapmanRichards = fit_gnls("Chapman-Richards GNLS", TotalHt ~ 1.37 + a1*(1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131), control = gnlsControl(nlsTol = 0.001))) # step halving at nlsTol = 1 with corSymm
  #rhpuHeightFromDiameterGnls$chapmanRichardsBal = fit_gnls("Chapman-Richards BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger + a3*standBasalAreaPerHectare) * (1 - exp(b1*DBH))^b2, rhpu2016, start = rhpuHeightFromDiameter$chapmanRichardsBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving at nlsTol = 0.2 with corSymm
  #rhpuHeightFromDiameterGnls$sharmaParton = fit_gnls("Sharma-Parton GNLS", TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, start = rhpuHeightFromDiameter$sharmaParton$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving at nlsTol = 0.2 with corSymm
  #rhpuHeightFromDiameterGnls$sharmaPartonBal = fit_gnls("Sharma-Parton BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = rhpuHeightFromDiameter$sharmaPartonBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving with plot correlation
  #rhpuHeightFromDiameterGnls$sharmaZhang = fit_gnls("Sharma-Zhang GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), rhpu2016, start = rhpuHeightFromDiameter$sharmaZhang$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001)) # step halving with plot correlation
  #rhpuHeightFromDiameterGnls$sharmaZhangBal = fit_gnls("Sharma-Zhang BA+L GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, rhpu2016, start = rhpuHeightFromDiameter$sharmaZhangBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # step halving with plot correlation
  #rhpuHeightFromDiameterGnls$weibull = fit_gnls("Weibull GNLS", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), rhpu2016, start = rhpuHeightFromDiameter$weibull$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001, maxIter = 250, nlsMaxIter = 50)) # corSymm() viable but dropped
  #rhpuHeightFromDiameterGnls$weibullBal = fit_gnls("Weibull BA+L GNLS", TotalHt ~ 1.37 + (a1 + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), rhpu2016, start = rhpuHeightFromDiameter$weibullBal$fit[[1]]$m$getPars(), control = gnlsControl(nlsTol = 0.001)) # step halving at nlsTol = 1 with corSymm
  
  save(file = "data/rhpu TotalHt gnls.Rdata", rhpuHeightFromDiameterGnls)
}
if (htDiaOptions$includeInvestigatory)
{
  rhpuHeightFromDiameterResultsGnls %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic) %>% arrange(method)
  
  #bind_cols(parameter = c("a1", "a2", "a3", "b1", "b2"), bal = confint2(rhpuHeightFromDiameter$weibullBAL, level = 0.99), balN = confint2(rhpuHeightFromDiameter$weibullBalNatural, level = 0.99), balP = confint2(rhpuHeightFromDiameter$weibullBalPlantation, level = 0.99)) %>%
  #  mutate(bal005 = bal[, 1], bal995 = bal[, 2], balN005 = balN[, 1], balN995 = balN[, 2], balP005 = balP[, 1], balP995 = balP[, 2]) %>%
  #  select(-bal, -balN, -balP)
  ggplot() +
    geom_point(aes(x = rhpu2016natural$DBH, y = rhpu2016natural$TotalHt), alpha = 0.15, color = "navyblue", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = rhpu2016natural$DBH, y = rhpu2016natural$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "natural regeneration DBH, cm", y = "cascara buckthorn naturally regenerated height, m") +
    ggplot() +
    geom_point(aes(x = rhpu2016plantation$DBH, y = rhpu2016plantation$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_smooth(aes(x = rhpu2016plantation$DBH, y = rhpu2016plantation$TotalHt), alpha = 0.20, color = "red", formula = y ~ s(x, k = 20), method = "gam", size = 0.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "plantation DBH, cm", y = "cascara buckthorn plantation height, m")
  
  ggplot() +
    geom_point(aes(x = rhpu2016$DBH, y = rhpu2016$TotalHt), alpha = 0.15, color = "black", na.rm = TRUE, shape = 16) +
    geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$weibullBal), color = "Weibull BA+L"), alpha = 0.5) + # Temesgen et al. 2007, Eq. 5
    geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$power), color = "power")) +
    geom_line(aes(x = rhpu2016$DBH, y = predict(rhpuHeightFromDiameter$weibull), color = "Weibull")) +
    annotate("text", x = 0, y = 85, label = "a) cascara buckthorn, height from diameter", hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    scale_color_manual(breaks = c("base", "ElliottWeibull", "ElliottBAL", "ElliottBALn", "ElliottBALp", "TemesgenWeibull"), labels = c(bquote("1.37 + b"[0]*"DBH"^{b[1]}), "Weibull", "Weibull with BA+L", "Weibull with BA+L, natural regeneration", "Weibull with BA+L, plantation", "Weibull, Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}


if (rhpuOptions$fitHeightMixed){ #fitting height diameter using mixed effect models
  rhpuHeightFromDiameterMixed = list(chapmanRichards = fit_nlme("Chapman-Richards", TotalHt ~ 1.37 + (a1 + a1r)*(1 - exp(b1*DBH))^b2, rhpu2016, 
                                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                start = list(fixed = c(a1 = 48.2, b1 = -0.015, b2 = 1.131)), control = nlmeControl(maxIter = 250)))
  rhpuHeightFromDiameterMixed$chapmanRichardsBal = fit_nlme("Chapman-Richards BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger) * (1 - exp(b1*DBH))^b2, rhpu2016, 
                                                            fixedFormula = a1 + a1p + a2 + a2p + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                            start = list(fixed = c(a1 = 55, a1p = -10, a2 = -0.1, a2p = 0.6, b1 = -0.012, b2 = 1.1)), control = nlmeControl(maxIter = 500))
  rhpuHeightFromDiameterMixed$chapmanRichardsBalPhysio = fit_nlme("Chapman-Richards BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation)*DBH))^b2, rhpu2016, 
                                                                  fixedFormula = a1 + a1p + a2 + a2p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                                  start = list(fixed = c(a1 = 50.8, a1p = -14.4, a2 = -0.09, a2p = 0.47, a8 = 0.23, b1 = -0.013, b1p = -0.003, b2 = 1.12)), control = nlmeControl(maxIter = 250), significant = FALSE)
  rhpuHeightFromDiameterMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex) * (1 - exp((b1 + b1p * isPlantation) * DBH))^b2, rhpu2016, 
                                                               fixedFormula = a1 + a1p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1, start = list(fixed = c(a1 = 52.0, a1p = -19.0, a8 = 0.20, b1 = -0.013, b1p = -0.009, b2 = 1.15)))
  # rhpuHeightFromDiameterMixed$curtis = fit_nlme("Curtis", TotalHt ~ 1.37 + (a1 + a1r) * DBH / (1 + DBH)^b1, rhpu2016, 
  #                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                               start = list(fixed = c(a1 = 0.560, b1 = 0.069)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # max iterations in job
  rhpuHeightFromDiameterMixed$hossfeld = fit_nlme("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) / (1 + (b1 + b1p * isPlantation) *DBH^b2), rhpu2016, 
                                                  fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                  start = list(fixed = c(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176)), control = nlmeControl(maxIter = 250))
  rhpuHeightFromDiameterMixed$korf = fit_nlme("Korf", TotalHt ~ 1.37 + (a1 + a1r)*exp(b1*DBH^b2), rhpu2016, 
                                              fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                              start = list(fixed = c(a1 = 1825, b1 = -8.726, b2 = -0.175)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  rhpuHeightFromDiameterMixed$michaelisMenten = fit_nlme("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), rhpu2016, 
                                                         fixedFormula = a1 + a1p + a2 + a2p + b1 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176)), control = nlmeControl(maxIter = 250)) # job >100 iterations
  rhpuHeightFromDiameterMixed$prodan = fit_nlme("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3 + a3r), rhpu2016, 
                                                fixedFormula = a1 + a1p + a2 + a3 ~ 1, randomFormula = a3r ~ 1,
                                                start = list(fixed = c(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649)))
  rhpuHeightFromDiameterMixed$power = fit_nlme("power", TotalHt ~ 1.37 + (a1 + a1r)*DBH^b1, rhpu2016, 
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1,
                                               start = list(fixed = c(a1 = 0.542, b1 = 0.939)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5)) # job >500 iterations without relaxed tolerances
  #rhpuHeightFromDiameterMixed$ratkowsky = fit_nlme("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), rhpu2016, 
  # fixedFormula = a1 + a1p + b1 + b1p + b2 + b2p ~ 1, randomFormula = a1r ~ 1,
  #start = list(fixed = c(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151)))
  #rhpuHeightFromDiameterMixed$richardsW = fit_nlme("unified Richards", TotalHt ~ 1.37 + (Ha + Hap*isPlantation + Har) * (1 + ((1.37/(Ha + Hap*isPlantation + Har))^(1 - d) - 1) * exp((-(kU + kUp * isPlantation) * DBH)/d^(d/(1 - d))))^(1/(1 - d)), rhpu2016, 
  #fixedFormula = Ha + Hap + d + kU + kUp ~ 1, randomFormula = Har ~ 1,
  #start = list(fixed = c(Ha = 52, Hap = -20, d = 0.5, kU = 0.008, kUp = 0.008)))
  rhpuHeightFromDiameterMixed$sharmaParton = fit_nlme("Sharma-Parton", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, rhpu2016, 
                                                      fixedFormula = a1 + b1 + b1p + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                      start = list(fixed = c(a1 = 38.0, b1 = 0.131, b1p = -0.135, b2 = -0.015, b2p = -0.011, b3 = -0.114, b4 = 1.09)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  rhpuHeightFromDiameterMixed$sharmaPartonBal = fit_nlme("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1r)*topHeight^b1 * (1 - exp(b2*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, 
                                                         fixedFormula = a1 + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                         start = list(fixed = c(a1 = 38, b1 = 0.1, b2 = -0.013, b3 = -0.1, b4 = 1.03)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve, step halving
  rhpuHeightFromDiameterMixed$sharmaPartonBalPhysio = fit_nlme("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, 
                                                               fixedFormula = a1 + a1p + a8 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                               start = list(fixed = c(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  rhpuHeightFromDiameterMixed$sharmaPartonPhysio = fit_nlme("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, 
                                                            fixedFormula = a1 + a1p + a8 + b1 + b2 + b2p + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                            start = list(fixed = c(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 1E-3)) # singular precision matrix, step halving
  rhpuHeightFromDiameterMixed$sharmaZhang = fit_nlme("Sharma-Zhang", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*standBasalAreaPerHectare^b1*(1 - exp(b2*tph^b3*DBH))^(b4 + b4p * isPlantation), rhpu2016, 
                                                     fixedFormula = a1 + a1p + b1 + b2 + b3 + b4 + b4p ~ 1, randomFormula = a1r ~ 1,
                                                     start = list(fixed = c(a1 = 40.1, a1p = -4.259, b1 = 0.040, b2 = -0.042, b3 = -0.148, b4 = 1.190, b4p = -0.097)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001))  # singularity in backsolve
  rhpuHeightFromDiameterMixed$sharmaZhangBal = fit_nlme("Sharma-Zhang BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r + (a2 + a2p * isPlantation) * basalAreaLarger)*standBasalAreaPerHectare^b1 * (1 - exp(b2*tph^b3*DBH))^b4, rhpu2016, 
                                                        fixedFormula = a1 + a1p + a2 + a2p + b1 + b2 + b3 + b4 ~ 1, randomFormula = a1r ~ 1,
                                                        start = list(fixed = c(a1 = 45, a1p = -7, a2 = -0.1, a2p = 0.4, b1 = -0.05, b2 = -0.02, b3 = -0.078, b4 = 1.08)), control = nlmeControl(tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  #rhpuHeightFromDiameterMixed$sibbesen = fit_nlme("Sibbesen", TotalHt ~ 1.37 + a1*DBH^((b1 + b1r)*DBH^b2), rhpu2016, 
  #                                                fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = b1r ~ 1,
  #                                                start = list(fixed = c(a1 = 0.302, b1 = 1.495, b2 = -0.078))) # a1r: step halving, singular precision
  rhpuHeightFromDiameterMixed$weibull = fit_nlme("Weibull", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a1r)*(1 - exp((b1 + b1p * isPlantation)*DBH^b2)), rhpu2016, 
                                                 fixedFormula = a1 + a1p + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1,
                                                 start = list(fixed = c(a1 = 49.3, a1p = -13.8, b1 = -0.007, b1p = -0.004, b2 = 1.141)))
  rhpuHeightFromDiameterMixed$weibullBal = fit_nlme("Weibull BA+L", TotalHt ~ 1.37 + (a1 + a1r + (a2 + a2p * isPlantation) * basalAreaLarger + (a3 + a3p * isPlantation) * standBasalAreaPerHectare) * (1 - exp(b1*DBH^b2)), rhpu2016,
                                                    fixedFormula = a1 + a2 + a2p + a3 + a3p + b1 + b2 ~ 1, randomFormula = a1r ~ 1,
                                                    start = list(fixed = c(a1 = 45.4, a2 = -0.178, a2p = 0.581, a3 = 0.096, a3p = -0.258, b1 = -0.008, b2 = 1.131)))
  
  rhpuHeightFromDiameterMixed$gamm = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8) + s(StandID, bs = "re"), data = rhpu2016, mixed = TRUE)
  rhpuHeightFromDiameterMixed$gammBal = fit_gam("REML GAM BA+L", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 13) + s(StandID, bs = "re"), data = rhpu2016, mixed = TRUE)
  
  save(file = "data/rhpu TotalHt mixed.Rdata", rhpuHeightFromDiameterMixed)
}


## Cascara buckthorn diameter-height regressions
if (rhpuOptions$fitDbh) {
  rhpuDiameterFromHeight = list(linear = fit_lm("linear", DBH ~ 0 + I(TotalHt - 1.37), rhpu2016)) # isPlantation*(TotalHt - 1.37) not significant
  rhpuDiameterFromHeight$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), rhpu2016) # (TotalHt - 1.37)^2 not significant
  
  rhpuDiameterFromHeight$chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016, start = list(a1 = 200, b1 = 0.01, b2 = 0.95), control = gsl_nls_control(maxiter = 500, xtol = 1E-5)) # a1p, b1p, b2p not significant, a1-b1 parameter evaporation: singular gradient with nls(), no convergence from nls_multstart(), NaN-inf with nlrob()
  #rhpuDiameterFromHeight$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016, start = list(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE) # NaN-inf with nls() and nlrob
   rhpuDiameterFromHeight$chapmanReplaceBal = fit_gsl_nls("Chapman-Richards replace BA+L", DBH ~ (a1 + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, start = list(a1 = 200, a2 = -10, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 300), significant = FALSE) # step size with nls() and nlrob()
   rhpuDiameterFromHeight$chapmanReplaceBalRelHt = fit_gsl_nls("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a2 * basalAreaLarger + a9 * pmin(relativeHeight, 1.5)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, start = list(a1 = 10, a2 = 0, a9 = 2.3, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250, xtol = 0.001), significant = FALSE) # a2, a3 not significant, a1-b1 parameter evaporation: nlrob() step factor with either a2 or a3
   rhpuDiameterFromHeight$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, start = list(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8), control = gsl_nls_control(maxiter = 500)) # step size with nls(), >500 iterations with nlrob()
  rhpuDiameterFromHeight$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, start = list(a1 = -200, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250)) # a1p and b2p not significant, poor convergence with b1p, step factor with nlrob()
  rhpuDiameterFromHeight$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, start = list(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE) # a1p, b1p not significant, step factor with nlrob()
  rhpuDiameterFromHeight$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, start = list(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 5E-5)) # no physiographic effects significant, a1-b1 parameter evaporation: step factor with nlrob()
  rhpuDiameterFromHeight$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, start = list(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9), control = gsl_nls_control(maxiter = 500), significant = FALSE) # step factor with nlrob()
  rhpuDiameterFromHeight$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), rhpu2016, start = list(a1 = 519, a2 = 237, b1 = 1.00)) # a1p, a2p, b1p not significant, singular gradient with nlrob()
  #rhpuDiameterFromHeight$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), rhpu2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024))
  rhpuDiameterFromHeight$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.93, b1 = 1.08)) # no significant plantation effects
  rhpuDiameterFromHeight$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09)) # no significant plantation effects
  rhpuDiameterFromHeight$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE) # no significant physiographic effects
  rhpuDiameterFromHeight$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13)) # a1p and b1p not significant
  rhpuDiameterFromHeight$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01)) # a1p, b1p, b2p not significant
  #rhpuDiameterFromHeight$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), significant = FALSE) # a2, a2p, a3, a3p, b1p, b2p not significant
  rhpuDiameterFromHeight$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.9, a2 = -0.005, a4 = -0.001, b1 = 0.93, b2 = 0.006), significant = FALSE) # a2, a3 not significant, no AIC discrimination
  rhpuDiameterFromHeight$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = -1, b1 = 0.9, b2 = 0), significant = FALSE) # a2, a3, a4, a9, b2 not significant, drop ABA on AIC
  #rhpuDiameterFromHeight$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.7, a3 = 0, a9 = 0, b1 = 0.95, b2 = 0.005), significant = FALSE) # a9, a9p, b2 not significant
  rhpuDiameterFromHeight$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE) # a1p, a5, a6, a7, a8, b1p, b2p not significant
  rhpuDiameterFromHeight$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE) # a9, a9p, b1p, b2, b2p not significant
  rhpuDiameterFromHeight$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 3.2, a4 = 0, a9 = -1, b1 = 0.9, b2 = 0.01), significant = FALSE) # a4, a9 not significant
  #rhpuDiameterFromHeight$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), rhpu2016, start = list(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30), control = gsl_nls_control(maxiter = 200)) # singular gradient with nlrob() and gsl_nls()
  #rhpuDiameterFromHeight$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, rhpu2016, start = list(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1), control = gsl_nls_control(maxiter = 250, xtol = 0.025)) # a1-b2 evaporation, b1, b3 not significant, NaN-inf with nls() from nls_multstart() point, NaN-inf, singular gradient, or code syntax error with nlrob()
  rhpuDiameterFromHeight$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12)) # no significant plantation effects
  rhpuDiameterFromHeight$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029), significant = FALSE) # no significant plantation effects
  rhpuDiameterFromHeight$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1), significant = FALSE) # a2, a3, a8 not significant, drop ABA on AIC
  rhpuDiameterFromHeight$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.3, a2 = 0, a8 = -0.017, a9 = 1.0, b1 = 0.7, b2 = 0), significant = FALSE) # a2, a3, a8, a9, b2 not significant, no a2-a3 AIC discrimination
  rhpuDiameterFromHeight$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.5, a2 = 0, a9 = 0, a9p = 0, b1 = 0.6, b2 = 0.12), significant = FALSE) # a2, a9, a9p, b2 not significant
  rhpuDiameterFromHeight$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE) # a1p, no physiographic effects significant
  rhpuDiameterFromHeight$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  rhpuDiameterFromHeight$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.6, a8 = -0.01, a9 = 0.7, b1 = 0.7, b2 = 0.1), significant = FALSE) # a9 not significant
  #rhpuDiameterFromHeight$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, rhpu2016, start = list(a1 = -300, b1 = 0.04, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 1E-4)) # a1p, b1p, b2p not significant, a1-b1 parameter evaporation: NaN-inf with nlrob()
  lapply(rhpuDiameterFromHeight$chapmanReplaceAbat$fit, confint2, level = 0.99)
  lapply(rhpuDiameterFromHeight$chapmanReplaceAbat$fit, get_model_coefficients)
  
  if (rhpuOptions$fitDbhNlrob)
  {
    rhpuDiameterFromHeightNlrob = list(naslund = fit_nlrob("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), rhpu2016, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024)))
    #rhpuDiameterFromHeightNlrob$power = fit_nlrob("power", DBH ~ a1*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.93, b1 = 1.08))
    #rhpuDiameterFromHeightNlrob$powerAbat = fit_nlrob("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09))
    #rhpuDiameterFromHeightNlrob$powerPhysio = fit_nlrob("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE)
    #rhpuDiameterFromHeightNlrob$powerRelHt = fit_nlrob("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13))
    rhpuDiameterFromHeightNlrob$ruark = fit_nlrob("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01))
    rhpuDiameterFromHeightNlrob$ruarkAbat = fit_nlrob("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # occasional job step factor
    rhpuDiameterFromHeightNlrob$ruarkAbatPhysio = fit_nlrob("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.9, a2 = 0, a4 = -0.001, b1 = 0.93, b2 = 0.006), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # step factor
    rhpuDiameterFromHeightNlrob$ruarkAbatPhysioRelHt = fit_nlrob("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = 2.3, b1 = 0.9, b2 = 0), control = nls.control(maxiter = 100, tol = 1E-4), significant = FALSE) # job step factor
    rhpuDiameterFromHeightNlrob$ruarkAbatRelHt = fit_nlrob("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.4, a3 = 0, a9 = 3, b1 = 0.8, b2 = 0.01), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
    rhpuDiameterFromHeightNlrob$ruarkPhysio = fit_nlrob("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE)
    rhpuDiameterFromHeightNlrob$ruarkRelHt = fit_nlrob("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE)
    rhpuDiameterFromHeightNlrob$ruarkRelHtPhysio = fit_nlrob("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, start = list(a1 = 3.2, a4 = -0.001, a9 = 1, b1 = 0.8, b2 = 0), significant = FALSE) # a4, a9, b2 not significant
    rhpuDiameterFromHeightNlrob$sibbesenReplace = fit_nlrob("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12))
    rhpuDiameterFromHeightNlrob$sibbesenReplaceAbat = fit_nlrob("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.1, a2 = -0.004, b1 = 0.7, b2 = 0.1), control = nls.control(tol = 1E-4)) # job step factor
    rhpuDiameterFromHeightNlrob$sibbesenReplaceAbatPhysio = fit_nlrob("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1), significant = FALSE)
    rhpuDiameterFromHeightNlrob$sibbesenReplaceAbatPhysioRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.3, a2 = 0, a8 = -0.01, a9 = 0.5, b1 = 0.7, b2 = 0.1), control = nls.control(tol = 1E-4), significant = FALSE) # job step factor
    rhpuDiameterFromHeightNlrob$sibbesenReplaceAbatRelHt = fit_nlrob("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.1, a2 = 0, a9 = 0, a9p = 0.7, b1 = 0.6, b2 = 0.12), significant = FALSE)
    rhpuDiameterFromHeightNlrob$sibbesenReplacePhysio = fit_nlrob("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
    rhpuDiameterFromHeightNlrob$sibbesenReplaceRelHt = fit_nlrob("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
    rhpuDiameterFromHeightNlrob$weibull = fit_nlrob("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, rhpu2016, start = list(a1 = -250, b1 = 0.043, b2 = 0.58), control = nls.control(maxiter = 500))
    confint_nlrob(rhpuDiameterFromHeight$sibbesenReplacePhysio, level = 0.99, weights = pmin(rhpu2016$TotalHt^if_else(rhpu2016$isPlantation, -1.7, -1.6), 0.5))
  } else {
    rhpuDiameterFromHeightNlrob = list()
  }
  rhpuDiameterFromHeightGslNlsDefault = list(chapmanReplace = fit_gsl_nls("Chapman-Richards replace", DBH ~ a1*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016defaultWeight, start = list(a1 = 200, b1 = 0.01, b2 = 0.95), control = gsl_nls_control(maxiter = 250, xtol = 1E-5)))
  rhpuDiameterFromHeightGslNlsDefault$chapmanReplaceAbat = fit_gsl_nls("Chapman-Richards replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016defaultWeight, start = list(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016defaultWeight, start = list(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8), control = gsl_nls_control(maxiter = 500))
  rhpuDiameterFromHeightGslNlsDefault$chapmanRichards = fit_gsl_nls("Chapman-Richards inverse", DBH ~ a1*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016defaultWeight, start = list(a1 = -200, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 250))
  rhpuDiameterFromHeightGslNlsDefault$chapmanRichardsAbat = fit_gsl_nls("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016defaultWeight, start = list(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$chapmanRichardsPhysio = fit_gsl_nls("Chapman-Richards inverse physio", DBH ~ (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016defaultWeightPhysio, start = list(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 5E-5))
  rhpuDiameterFromHeightGslNlsDefault$chapmanRichardsRelHt = fit_gsl_nls("Chapman-Richards inverse RelHt", DBH ~ (a1 + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016defaultWeight, start = list(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9), control = gsl_nls_control(maxiter = 500), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$michaelisMentenReplace = fit_gsl_nls("Michaelis-Menten replace", DBH ~ a1 * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), rhpu2016defaultWeight, start = list(a1 = 519, a2 = 237, b1 = 1.00))
  rhpuDiameterFromHeightGslNlsDefault$naslund = fit_gsl_nls("Näslund inverse", DBH ~ (a1 + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), rhpu2016defaultWeight, start = list(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024))
  rhpuDiameterFromHeightGslNlsDefault$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, rhpu2016defaultWeight, start = list(a1 = 1.93, b1 = 1.08))
  rhpuDiameterFromHeightGslNlsDefault$powerAbat = fit_gsl_nls("power ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, rhpu2016defaultWeight, start = list(a1 = 1.94, a2 = -0.00051, b1 = 1.09))
  rhpuDiameterFromHeightGslNlsDefault$powerPhysio = fit_gsl_nls("power physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, rhpu2016defaultWeightPhysio, start = list(a1 = 2.26, a8 = -0.0060, b1 = 1.08), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$powerRelHt = fit_gsl_nls("power RelHt", DBH ~ (a1 + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, rhpu2016defaultWeight, start = list(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13))
  rhpuDiameterFromHeightGslNlsDefault$ruark = fit_gsl_nls("Ruark", DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeight, start = list(a1 = 2.8, b1 = 0.9, b2 = 0.01))
  rhpuDiameterFromHeightGslNlsDefault$ruarkAbat = fit_gsl_nls("Ruark ABA+T", DBH ~ (a1 + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeight, start = list(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$ruarkAbatPhysio = fit_gsl_nls("Ruark ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeightPhysio, start = list(a1 = 1.6, a2 = -0.01, a4 = -0.0006, b1 = 1.2, b2 = -0.009), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$ruarkAbatPhysioRelHt = fit_gsl_nls("Ruark ABA+T RelHt physio", DBH ~ (a1 + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeightPhysio, start = list(a1 = 1.6, a3 = -0.003, a4 = -0.0006, a9 = 0.4, b1 = 1.27, b2 = -0.01), significant = FALSE) 
  rhpuDiameterFromHeightGslNlsDefault$ruarkAbatRelHt = fit_gsl_nls("Ruark ABA+T RelHt", DBH ~ (a1 + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeight, start = list(a1 = 1.3, a3 = -0.003, a9 = 0.25, b1 = 1.3, b2 = -0.008), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$ruarkPhysio = fit_gsl_nls("Ruark physio", DBH ~ (a1 + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeightPhysio, start = list(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$ruarkRelHt = fit_gsl_nls("Ruark RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeight, start = list(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$ruarkRelHtPhysio = fit_gsl_nls("Ruark RelHt physio", DBH ~ (a1 + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016defaultWeightPhysio, start = list(a1 = 1.6, a4 = -0.0005, a9 = -0.4, b1 = 1.2, b2 = -0.01), significant = FALSE) # a4, a9 not significant
  #rhpuDiameterFromHeightGslNlsDefault$schnute = fit_gsl_nls("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/(Ha^b1 - 1.3^b1)), rhpu2016defaultWeight, start = list(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30), control = gsl_nls_control(maxiter = 200))
  #rhpuDiameterFromHeightGslNlsDefault$sharmaParton = fit_gsl_nls("modified Sharma-Parton", DBH ~ a1*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, rhpu2016defaultWeight, start = list(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1), control = gsl_nls_control(maxiter = 250, xtol = 0.025))
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplace = fit_gsl_nls("Sibbesen replace", DBH ~ a1*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeight, start = list(a1 = 3.4, b1 = 0.8, b2 = 0.12))
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceAbat = fit_gsl_nls("Sibbesen replace ABA+T", DBH ~ (a1 + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeight, start = list(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029))
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysio = fit_gsl_nls("Sibbesen replace ABA+T physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeightPhysio, start = list(a1 = 1.5, a2 = -0.009, a8 = -0.005, b1 = 1.2, b2 = -0.04), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatPhysioRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeightPhysio, start = list(a1 = 1.41, a2 = -0.009, a8 = -0.005, a9 = 0, b1 = 1.4, b2 = -0.05), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceAbatRelHt = fit_gsl_nls("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeight, start = list(a1 = 1.5, a2 = -0.008, a9 = 0, a9p = 0, b1 = 1.4, b2 = 0), significant = FALSE)
  #rhpuDiameterFromHeightGslNlsDefault$sibbesenReplacePhysio = fit_gsl_nls("Sibbesen replace physio", DBH ~ (a1 + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeightPhysio, start = list(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHt = fit_gsl_nls("Sibbesen replace RelHt", DBH ~ (a1 + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeight, start = list(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1), significant = FALSE)
  #rhpuDiameterFromHeightGslNlsDefault$sibbesenReplaceRelHtPhysio = fit_gsl_nls("Sibbesen replace RelHt physio", DBH ~ (a1 + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016defaultWeightPhysio, start = list(a1 = 1.4, a8 = 0, a9 = 0.3, b1 = 1.3, b2 = -0.035), significant = FALSE) # a8, a9 not significant
  #rhpuDiameterFromHeightGslNlsDefault$weibull = fit_gsl_nls("Weibull inverse", DBH ~ (a1*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, rhpu2016defaultWeight, start = list(a1 = -300, b1 = 0.04, b2 = 0.55), control = gsl_nls_control(maxiter = 250, xtol = 1E-4))
  
  # individual term selection: TotalHt by = isPlantation only, AAT retained by AIC but not significant (p = 0.38), #did not run this part of the code because the variable 'pc=gamConstraint' is not defined, or not any information on how it may be defined
  rhpuDiameterFromHeight$gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # newton() step failure with scat()
  rhpuDiameterFromHeight$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
  rhpuDiameterFromHeight$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
  rhpuDiameterFromHeight$gamAbatPhysioRelHt = fit_gam("REML GAM ABA+T RelHt physio", DBH ~ s(TotalHt, standBasalAreaApprox, topographicShelterIndex, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 22, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # drop ABA and elevation on AIC
  rhpuDiameterFromHeight$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # drop elevation and topographic shelter on AIC
  rhpuDiameterFromHeight$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint)
  rhpuDiameterFromHeight$gamRelHtPhysio = fit_gam("REML GAM RelHt physio", DBH ~ s(TotalHt, slope, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 57, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint) # drop elevation and aspect on AIC
  
  save(file = "data/rhpu DBH.Rdata", rhpuDiameterFromHeight, rhpuDiameterFromHeightNlrob, rhpuDiameterFromHeightGslNlsDefault)
}
if (htDiaOptions$includeInvestigatory) {
  print(rhpuDiameterFromHeightResults %>% select(-responseVariable, -species, -biasNR, -biasPl, -rmse, -rmseNR, -rmsePl, -pearsonNR, -pearsonPl, -aic, -bic), n = 25)
  
  ggplot(rhpu2016) +
    geom_point(aes(x = DBH, y = TotalHt), alpha = 0.10, color = "grey25", shape = 16) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$sharmaParton), y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$chapmanReplace), y = TotalHt, color = "Chapman-Richards replace", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$chapmanReplaceAbat), y = TotalHt, color = "Chapman-Richards replace approximate BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$chapmanReplaceBal), y = TotalHt, color = "Chapman-Richards replace BA+L", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$chapmanRichards), y = TotalHt, color = "Chapman-Richards", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$michaelisMentenReplace), y = TotalHt, color = "Michaelis-Menten replace", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$naslund), y = TotalHt, color = "Näslund", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$power), y = TotalHt, color = "power", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$ruark), y = TotalHt, color = "Ruark", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$schnute), y = TotalHt, color = "Schnute inverse", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$sibbesenReplace), y = TotalHt, color = "Sibbesen replace", group = isPlantation)) +
    #geom_line(aes(x = predict(rhpuDiameterFromHeight$weibull), y = TotalHt, color = "Weibull", group = isPlantation)) +
    #geom_line(aes(x = 5*standBasalAreaPerHectare^0.5 * exp(0.0005*tph^0.5*(TotalHt - 1.37))^1, y = TotalHt, color = "Sharma-Zhang"), alpha = 0.5) +
    #geom_line(aes(x = -100 * log(1 - pmin(0.015*(TotalHt - 1.37)^1.0, 0.999)), y = TotalHt, color = "Chapman-Richards inversion"), na.rm = TRUE) +
    #geom_line(aes(x = 0.5*(TotalHt - 1.37)^1*exp(0.02*(tph/topHeight)^0.26*(TotalHt - 1.37))^0.9, y = TotalHt, color = "modified Sharma-Parton", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.45, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 15 * (exp(0.1*(TotalHt - 1.37)) - 1)^0.35, y = TotalHt, color = "Chapman-Richards replace", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = (1.75 + 0.000001 * tallerApproxBasalArea + -0.000001 * standBasalAreaApprox) * exp(1.46*(TotalHt - 1.37)^0.280), y = TotalHt, color = "Chapman-Richards replace ABA+T", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = 0.03*topHeight*exp(1.6*(TotalHt - 1.37)^0.26), y = TotalHt, color = "Chapman-Richards replace top height", group = isPlantation), alpha = 0.5) +
    #geom_line(aes(x = -1/0.0003*log(1 - (1 - exp(-0.1))*(TotalHt^1.5 - 1.37^1.5)/(75^1.5 - 1.37^1.5)), y = TotalHt, color = "Schnute inverse"), alpha = 0.5) +
    geom_line(aes(x = 30*topHeight^0.5*(exp(0.01 * (tph/standBasalAreaPerHectare)^0.25*(TotalHt - 1.37)) - 1)^0.5, y = TotalHt, color = "modified Sharma-Parton"), alpha = 0.5) +
    annotate("text", x = 0, y = 62, label = "cascara buckthorn, diameter from height", hjust = 0, size = 3.5) +
    #coord_cartesian(xlim = c(0, 250), ylim = c(0, 90)) +
    labs(x = "DBH, cm", y = "height, m", color = NULL) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("grey25", "transparent", "red")) +
    #scale_color_manual(breaks = c(FALSE, TRUE, "Chapman-Richards"), values = c("transparent", "grey25", "red")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.03))
}

if (rhpuOptions$fitDbhMixed) {
  rhpuDiameterFromHeightMixed = list(chapmanReplace = fit_nlme("Chapman-Richards replace", DBH ~ (a1 + a1r)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016, 
                                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                               start = list(fixed = c(a1 = 200, b1 = 0.01, b2 = 0.95)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001))) # singularity in backsolve, max iterations
  rhpuDiameterFromHeightMixed$chapmanReplaceAbat = fit_nlme("Chapman-Richards replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(exp(b1*(TotalHt - 1.37)) - 1)^b2, rhpu2016, 
                                                            fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                            start = list(fixed = c(a1 = 200, a2 = 0, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singularity in backsolve
  #rhpuDiameterFromHeightMixed$chapmanReplaceBal = fit_nlme("Chapman-Richards replace BA+L", DBH ~ (a1 + a1r + a2 * basalAreaLarger) * (exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, 
  #                                                         fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                         start = list(fixed = c(a1 = 200, a2 = -10, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 300, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving
  #rhpuDiameterFromHeightMixed$chapmanReplaceBalRelHt = fit_nlme("Chapman-Richards replace BA+L RelHt", DBH ~ (a1 + a1r + a2 * basalAreaLarger + a9 * pmin(relativeHeight, 1.5)) * (exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, 
  #                                                              fixedFormula = a1 + a2 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                              start = list(fixed = c(a1 = 10, a2 = 0, a9 = 2.3, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 250, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # singularity in backsolve
  #rhpuDiameterFromHeightMixed$chapmanReplaceRelHt = fit_nlme("Chapman-Richards replace RelHt", DBH ~ (a1 + a1r + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, 
  #                                                           fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                           start = list(fixed = c(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  rhpuDiameterFromHeightMixed$chapmanRichards = fit_nlme("Chapman-Richards inverse", DBH ~ (a1 + a1r)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = -200, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations
  rhpuDiameterFromHeightMixed$chapmanRichardsAbat = fit_nlme("Chapman-Richards inverse ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                             start = list(fixed = c(a1 = -200, a2 = 0, b1 = 0.01, b2 = 1.0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations, step halving
  #rhpuDiameterFromHeightMixed$chapmanRichardsPhysio = fit_nlme("Chapman-Richards inverse physio", DBH ~ (a1 + a1r + a1p * isPlantation + a8 * topographicShelterIndex)*log(1 - pmin((b1 + b1p * isPlantation)*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, 
  #                                                             fixedFormula = a1 + a1p + a8 + b1 + b1p + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                             start = list(fixed = c(a1 = -70, a1p = 40, a8 = 0.3, b1 = 0.01, b1p = 0.03, b2 = 0.55)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations, step halving
  rhpuDiameterFromHeightMixed$chapmanRichardsRelHt = fit_nlme("Chapman-Richards inverse RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*log(1 - pmin(b1*(TotalHt - 1.37)^b2, 0.9999)), rhpu2016, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                              start = list(fixed = c(a1 = -200, a9 = -70, b1 = 0.01, b2 = 0.9)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # step halving, singularity in backsolve
  rhpuDiameterFromHeightMixed$michaelisMentenReplace = fit_nlme("Michaelis-Menten replace", DBH ~ (a1 + a1r) * (TotalHt - 1.37)^b1 / (a2 - (TotalHt - 1.37)^b1), rhpu2016, 
                                                                fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1, 
                                                                start = list(fixed = c(a1 = 519, a2 = 237, b1 = 1.00)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # max iterations, step halving
  rhpuDiameterFromHeightMixed$naslund = fit_nlme("Näslund inverse", DBH ~ (a1 + a1r + a1p * isPlantation) * sqrt(TotalHt - 1.37) / (1 + (a2 + a2p * isPlantation) * sqrt(TotalHt - 1.37)), rhpu2016, 
                                                 fixedFormula = a1 + a1p + a2 + a2p ~ 1, randomFormula = a1r ~ 1, 
                                                 start = list(fixed = c(a1 = 5.1, a1p = -1.6, a2 = -0.11, a2p = -0.024)))
  rhpuDiameterFromHeightMixed$power = fit_nlme("power", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1, rhpu2016, 
                                               fixedFormula = a1 + b1 ~ 1, randomFormula = a1r ~ 1, 
                                               start = list(fixed = c(a1 = 1.93, b1 = 1.08)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4)) # job max iterations
  #rhpuDiameterFromHeightMixed$powerAbat = fit_nlme("power ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^b1, rhpu2016, 
  #                                                 fixedFormula = a1 + a2 + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                 start = list(fixed = c(a1 = 1.94, a2 = -0.00051, b1 = 1.09)))
  #rhpuDiameterFromHeightMixed$powerPhysio = fit_nlme("power physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex)*(TotalHt - 1.37)^b1, rhpu2016, 
  #                                                   fixedFormula = a1 + a8 + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                   start = list(fixed = c(a1 = 2.26, a8 = -0.0060, b1 = 1.08)), significant = FALSE)
  #rhpuDiameterFromHeightMixed$powerRelHt = fit_nlme("power RelHt", DBH ~ (a1 + a1r + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^b1, rhpu2016, 
  #                                                  fixedFormula = a1 + a9 + a9p + b1 ~ 1, randomFormula = a1r ~ 1, 
  #                                                  start = list(fixed = c(a1 = 1.68, a9 = -0.11, a9p = 0.23, b1 = 1.13)))
  rhpuDiameterFromHeightMixed$ruark = fit_nlme("Ruark", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                               fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                               start = list(fixed = c(a1 = 2.8, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations
  rhpuDiameterFromHeightMixed$ruarkAbat = fit_nlme("Ruark ABA+T", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                                   fixedFormula = a1 + a3 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                   start = list(fixed = c(a1 = 2.7, a3 = -0.003, b1 = 0.95, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # job max iterations
  rhpuDiameterFromHeightMixed$ruarkAbatPhysio = fit_nlme("Ruark ABA+T physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                                         fixedFormula = a1 + a2 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = 2.9, a2 = -0.005, a4 = -0.001, b1 = 0.93, b2 = 0.006)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # job max iterations
  #rhpuDiameterFromHeightMixed$ruarkAbatPhysioRelHt = fit_nlme("Ruark ABA+T RelHt physio", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
  #                                                            fixedFormula = a1 + a3 + a4 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                            start = list(fixed = c(a1 = 3.2, a3 = 0, a4 = -0.002, a9 = -1, b1 = 0.9, b2 = 0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations
  rhpuDiameterFromHeightMixed$ruarkAbatRelHt = fit_nlme("Ruark ABA+T RelHt", DBH ~ (a1 + a1r + a3 * standBasalAreaApprox + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                                        fixedFormula = a1 + a3 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                        start = list(fixed = c(a1 = 2.7, a3 = 0, a9 = 0, b1 = 0.95, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  rhpuDiameterFromHeightMixed$ruarkPhysio = fit_nlme("Ruark physio", DBH ~ (a1 + a1r + a4 * elevation)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                                     fixedFormula = a1 + a4 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                     start = list(fixed = c(a1 = 2.9, a4 = -0.001, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  #rhpuDiameterFromHeightMixed$ruarkRelHt = fit_nlme("Ruark RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
  #                                                  fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                  start = list(fixed = c(a1 = 2.8, a9 = 0.5, b1 = 0.9, b2 = 0.005)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # job max iterations
  rhpuDiameterFromHeightMixed$ruarkRelHtPhysio = fit_nlme("Ruark RelHt physio", DBH ~ (a1 + a1r + a4 * elevation + a9 * relativeHeight)*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), rhpu2016, 
                                                          fixedFormula = a1 + a4 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                          start = list(fixed = c(a1 = 3.2, a4 = 0, a9 = -1, b1 = 0.9, b2 = 0.01)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations, false convergence
  #rhpuDiameterFromHeightMixed$schnute = fit_nlme("Schnute inverse", DBH ~ -1/a1 * log(1 - (1 - exp(-a2))*(TotalHt^b1 - 1.37^b1)/((Ha + Har)^b1 - 1.3^b1)), rhpu2016, 
  #                                               fixedFormula = a1 + a2 + b1 + Ha ~ 1, randomFormula = Har ~ 1, 
  #                                               start = list(fixed = c(a1 = 0.00005, a2 = 0.001, b1 = 1.05, Ha = 30)), control = nlmeControl(maxIter = 100, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # step halving
  #rhpuDiameterFromHeightMixed$sharmaParton = fit_nlme("modified Sharma-Parton", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^b1*(exp(b2*(TotalHt - 1.37)) - 1)^b4, rhpu2016, 
  #                                                    fixedFormula = a1 + b1 + b2 + b4 ~ 1, randomFormula = a1r ~ 1, 
  #                                                    start = list(fixed = c(a1 = 100, b1 = -0.15, b2 = 0.01, b4 = 1.1)), control = nlmeControl(maxIter = 250, tolerance = 0.1, pnlsTol = 1, msTol = 0.01)) # singularity in backsolve
  rhpuDiameterFromHeightMixed$sibbesenReplace = fit_nlme("Sibbesen replace", DBH ~ (a1 + a1r)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                         fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                         start = list(fixed = c(a1 = 3.4, b1 = 0.8, b2 = 0.12)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # job max iterations
  rhpuDiameterFromHeightMixed$sibbesenReplaceAbat = fit_nlme("Sibbesen replace ABA+T", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                             fixedFormula = a1 + a2 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                             start = list(fixed = c(a1 = 1.39, a2 = -0.00036, b1 = 1.31, b2 = -0.029)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5))
  rhpuDiameterFromHeightMixed$sibbesenReplaceAbatPhysio = fit_nlme("Sibbesen replace ABA+T physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                                   fixedFormula = a1 + a2 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                   start = list(fixed = c(a1 = 3.6, a2 = 0, a8 = -0.01, b1 = 0.7, b2 = 0.1)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  #rhpuDiameterFromHeightMixed$sibbesenReplaceAbatPhysioRelHt = fit_nlme("Sibbesen replace ABA+T RelHt physio", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
  #                                                                      fixedFormula = a1 + a2 + a8 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
  #                                                                      start = list(fixed = c(a1 = 3.3, a2 = 0, a8 = -0.017, a9 = 1.0, b1 = 0.7, b2 = 0)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001), significant = FALSE) # max iterations
  rhpuDiameterFromHeightMixed$sibbesenReplaceAbatRelHt = fit_nlme("Sibbesen replace ABA+T RelHt", DBH ~ (a1 + a1r + a2 * tallerApproxBasalArea + (a9 + a9p * isPlantation) * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                                  fixedFormula = a1 + a2 + a9 + a9p + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                  start = list(fixed = c(a1 = 3.5, a2 = 0, a9 = 0, a9p = 0, b1 = 0.6, b2 = 0.12)), control = nlmeControl(maxIter = 500, tolerance = 1E-4, pnlsTol = 0.01, msTol = 1E-5), significant = FALSE) # singular precision matrix
  rhpuDiameterFromHeightMixed$sibbesenReplacePhysio = fit_nlme("Sibbesen replace physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                               fixedFormula = a1 + a8 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                               start = list(fixed = c(a1 = 3.6, a8 = 0, b1 = 0.6, b2 = 0.1)), control = nlmeControl(maxIter = 500), significant = FALSE)
  rhpuDiameterFromHeightMixed$sibbesenReplaceRelHt = fit_nlme("Sibbesen replace RelHt", DBH ~ (a1 + a1r + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                              fixedFormula = a1 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                              start = list(fixed = c(a1 = 3.3, a9 = 0, b1 = 0.6, b2 = 0.1)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  rhpuDiameterFromHeightMixed$sibbesenReplaceRelHtPhysio = fit_nlme("Sibbesen replace RelHt physio", DBH ~ (a1 + a1r + a8 * topographicShelterIndex + a9 * relativeHeight)*(TotalHt - 1.37)^(b1*(TotalHt - 1.37)^b2), rhpu2016, 
                                                                    fixedFormula = a1 + a8 + a9 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                                    start = list(fixed = c(a1 = 3.0, a8 = -0.01, a9 = 0, b1 = 0.73, b2 = 0.07)), control = nlmeControl(maxIter = 500, tolerance = 0.001, pnlsTol = 0.1, msTol = 1E-4), significant = FALSE) # max iterations
  rhpuDiameterFromHeightMixed$weibull = fit_nlme("Weibull inverse", DBH ~ ((a1 + a1r)*log(1 - pmin(b1*(TotalHt - 1.37), 0.9999)))^b2, rhpu2016, 
                                                 fixedFormula = a1 + b1 + b2 ~ 1, randomFormula = a1r ~ 1, 
                                                 start = list(fixed = c(a1 = -300, b1 = 0.04, b2 = 0.55)), control = nlmeControl(maxIter = 500, tolerance = 0.01, pnlsTol = 1, msTol = 0.001)) # singularity in backsolve
  
  rhpuDiameterFromHeightMixed$gamm = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9) + s(StandID, bs = "re"), data = rhpu2016, mixed = TRUE)
  rhpuDiameterFromHeightMixed$gammAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16) + s(StandID, bs = "re"), data = rhpu2016, mixed = TRUE)
  rhpuDiameterFromHeightMixed$gammRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9) + s(StandID, bs = "re"), data = rhpu2016, mixed = TRUE)
  
  save(file = "data/rhpu DBH mixed.Rdata", rhpuDiameterFromHeightMixed)
}


# ## collect model results and parameters
# if (rhpuOptions$fitHeight & rhpuOptions$fitHeightMixed & rhpuOptions$fitDbh & rhpuOptions$fitDbhMixed) {
#   if (exists("rhpuHeightFromDiameter") == FALSE) { load("trees/height-diameter/data/rhpu TotalHt.Rdata") }
#   #if (exists("rhpuHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/rhpu TotalHt gnls.Rdata") }
#   if (exists("rhpuHeightFromDiameterMixed") == FALSE) { load("trees/height-diameter/data/rhpu TotalHt mixed.Rdata") }
#   if (exists("rhpuDiameterFromHeight") == FALSE) { load("trees/height-diameter/data/rhpu DBH.Rdata") }
#   if (exists("rhpuDiameterFromHeightMixed") == FALSE) { load("trees/height-diameter/data/rhpu DBH mixed.Rdata") }

## collect model results and parameters
if (rhpuOptions$fitHeight & rhpuOptions$fitHeightMixed & rhpuOptions$fitDbh & rhpuOptions$fitDbhMixed) {
  if (exists("rhpuHeightFromDiameter") == FALSE) { load("data/rhpu TotalHt.Rdata") }
  #if (exists("rhpuHeightFromDiameterGnls") == FALSE) { load("trees/height-diameter/data/rhpu TotalHt gnls.Rdata") }
  if (exists("rhpuHeightFromDiameterMixed") == FALSE) { load("data/rhpu TotalHt mixed.Rdata") }
  if (exists("rhpuDiameterFromHeight") == FALSE) { load("data/rhpu DBH.Rdata") }
  if (exists("rhpuDiameterFromHeightMixed") == FALSE) { load("data/rhpu DBH mixed.Rdata") }
  rhpuCoefficients = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_coefficients)),
                                         #bind_rows(lapply(rhpuHeightFromDiameterGnls, get_model_coefficients)),
                                         bind_rows(lapply(rhpuHeightFromDiameterGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(rhpuHeightFromDiameterMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(rhpuHeightFromDiameterNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_coefficients)),
                                         bind_rows(lapply(rhpuDiameterFromHeightGslNlsDefault, get_list_coefficients, fitSet = "gsl_nls", fixedWeight = -1)),
                                         bind_rows(lapply(rhpuDiameterFromHeightMixed, get_list_coefficients, fitSet = "mixed")),
                                         bind_rows(lapply(rhpuDiameterFromHeightNlrob, get_list_coefficients, fitSet = "nlrob"))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  rhpuResults = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_stats)),
                                    #bind_rows(lapply(rhpuHeightFromDiameterGnls, get_stats)),
                                    bind_rows(lapply(rhpuHeightFromDiameterGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(rhpuHeightFromDiameterMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(rhpuHeightFromDiameterNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitSet = "primary", fittingMethod = "gsl_nls"),
                                    bind_rows(lapply(rhpuDiameterFromHeightGslNlsDefault, get_list_stats, fitSet = "gsl_nls", fixedWeight = -1)),
                                    bind_rows(lapply(rhpuDiameterFromHeightMixed, get_list_stats, fitSet = "mixed")),
                                    bind_rows(lapply(rhpuDiameterFromHeightNlrob, get_list_stats, fitSet = "nlrob"))) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  
  check_plot_results(rhpuResults)
  save(file = "data/rhpu results.Rdata", rhpuCoefficients, rhpuResults)
} else if (rhpuOptions$fitHeight & rhpuOptions$fitHeightMixed & rhpuOptions$fitDbh & rhpuOptions$fitDbhMixed)
{
  if (exists("rhpuHeightFromDiameter") == FALSE) { load("data/rhpu TotalHt.Rdata") }
  if (exists("rhpuDiameterFromHeight") == FALSE) { load("data/rhpu DBH.Rdata") }
  
  rhpuCoefficients = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  rhpuResults = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  
  check_plot_results(rhpuResults)
  save(file = "data/rhpu results.Rdata", rhpuCoefficients, rhpuResults)
}else(rhpuOptions$fitHeight &rhpuOptions$fitDbh) #added for height and diameter fit only.
{
  if (exists("rhpuHeightFromDiameter") == FALSE) { load("data/rhpu TotalHt.Rdata") }
  if (exists("rhpuDiameterFromHeight") == FALSE) { load("data/rhpu DBH.Rdata") }
  
  rhpuCoefficients = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_coefficients))) %>%
                                 mutate(responseVariable = "height"),
                               bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_coefficients))) %>%
                                 mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  rhpuResults = bind_rows(bind_rows(bind_rows(lapply(rhpuHeightFromDiameter, get_list_stats))) %>%
                            mutate(responseVariable = "height"),
                          bind_rows(bind_rows(lapply(rhpuDiameterFromHeight, get_list_stats)),
                                    create_model_stats(name = "Schnute inverse", fitting = "gsl_nls", fitSet = "primary")) %>%
                            mutate(responseVariable = "DBH")) %>%
    mutate(species = "rhpu")
  
  check_plot_results(rhpuResults)
  save(file = "data/rhpu results.Rdata", rhpuCoefficients, rhpuResults)
}

## preferred forms identified (results.R, Figure 8)
if (rhpuOptions$fitHeight & rhpuOptions$fitDbh)
{
  rhpuHeightFromDiameterPreferred = list(chapmanRichards = fit_gsl_nls("Chapman-Richards", TotalHt ~ 1.37 + a1 * (1 - exp(b1*DBH))^b2, rhpu2016, start = list(a1 = 48.2, b1 = -0.015, b2 = 1.131), folds = 1, repetitions = 1))
  rhpuHeightFromDiameterPreferred$gam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  #rhpuHeightFromDiameterPreferred$gamBalPhysio = fit_gam("REML GAM BA+L physio", TotalHt ~ s(DBH, standBasalAreaPerHectare, basalAreaLarger, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 20, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  rhpuHeightFromDiameterPreferred$hossfeld = fit_gsl_nls("Hossfeld IV", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) / (1 + (b1 + b1p * isPlantation) *DBH^b2), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, b1 = 200, b1p = -68.2, b2 = -1.176), folds = 1, repetitions = 1)
  rhpuHeightFromDiameterPreferred$michaelisMenten = fit_gsl_nls("Michaelis-Menten", TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), rhpu2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176), folds = 1, repetitions = 1)
  rhpuHeightFromDiameterPreferred$prodan = fit_gsl_nls("Prodan", TotalHt ~ 1.37 + DBH^2 / ((a1 + a1p * isPlantation) * DBH^2 + a2*DBH + a3), rhpu2016, start = list(a1 = 0.011, a1p = 0.002, a2 = 1.600, a3 = 1.649), folds = 1, repetitions = 1)
  #rhpuHeightFromDiameterPreferred$sharmaPartonBal = fit_gsl_nls("Sharma-Parton BA+L", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 50.6, a1p = -15.8, b1 = 0.023, b2 = -0.014, b2p = -0.009, b3 = -0.069, b4 = 1.130), folds = 1, repetitions = 1)
  #rhpuHeightFromDiameterPreferred$sharmaPartonBalPhysio = fit_gsl_nls("Sharma-Parton BA+L physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare + basalAreaLarger))^b3*DBH))^b4, rhpu2016, start = list(a1 = 37.0, a1p = -13.4, a8 = 0.13, b1 = 0.11, b2 = -0.013, b2p = -0.012, b3 = -0.10, b4 = 1.10), folds = 1, repetitions = 1)
  #rhpuHeightFromDiameterPreferred$sharmaPartonPhysio = fit_gsl_nls("Sharma-Parton physio", TotalHt ~ 1.37 + (a1 + a1p * isPlantation + a8 * topographicShelterIndex)*topHeight^b1 * (1 - exp((b2 + b2p * isPlantation)*(tph/(standBasalAreaPerHectare))^b3*DBH))^b4, rhpu2016, start = list(a1 = 32.7, a1p = -11.6, a8 = 0.11, b1 = 0.13, b2 = -0.014, b2p = -0.014, b3 = -0.11, b4 = 1.09), folds = 1, repetitions = 1)
  rhpuHeightFromDiameterPreferred$ratkowsky = fit_gsl_nls("Ratkowsky", TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*exp((b1 + b1p * isPlantation)/(DBH + b2 + b2p * isPlantation)), rhpu2016, start = list(a1 = 62.9, a1p = -19.3, b1 = -61.8, b1p = 23.1, b2 = 13.3, b2p = -5.151), folds = 1, repetitions = 1)
  AIC(rhpuHeightFromDiameterPreferred$hossfeld, rhpuHeightFromDiameterPreferred$michaelisMenten, rhpuHeightFromDiameterPreferred$prodan, rhpuHeightFromDiameterPreferred$ratkowsky)
  
  rhpuDiameterFromHeightPreferred = list(gam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1))
  #rhpuDiameterFromHeightPreferred$chapmanReplaceRelHt = fit_gsl_nls("Chapman-Richards replace RelHt", DBH ~ (a1 + a9 * pmin(relativeHeight, 1.5))*(exp(b1*(TotalHt - 1.37)^b2) - 1), rhpu2016, start = list(a1 = 100, a9 = 2.3, b1 = 0.01, b2 = 0.8), control = gsl_nls_control(maxiter = 500), folds = 1, repetitions = 1)
  rhpuDiameterFromHeightPreferred$parabolic = fit_lm("parabolic", DBH ~ 0 + I(TotalHt - 1.37) + I(isPlantation*(TotalHt - 1.37)) + I(isPlantation*(TotalHt - 1.37)^2), rhpu2016, folds = 1, repetitions = 1)
  rhpuDiameterFromHeightPreferred$power = fit_gsl_nls("power", DBH ~ a1*(TotalHt - 1.37)^b1, rhpu2016, start = list(a1 = 1.93, b1 = 1.08), folds = 1, repetitions = 1)
  #rhpuDiameterFromHeightPreferred$gamAbat = fit_gam("REML GAM ABA+T", DBH ~ s(TotalHt, tallerApproxBasalArea, standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  #rhpuDiameterFromHeightPreferred$gamAbatPhysio = fit_gam("REML GAM ABA+T physio", DBH ~ s(TotalHt, tallerApproxBasalArea, slope, bs = "ts", by = as.factor(isPlantation), k = 16, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  #rhpuDiameterFromHeightPreferred$gamPhysio = fit_gam("REML GAM physio", DBH ~ s(TotalHt, topographicShelterIndex, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  rhpuDiameterFromHeightPreferred$gamRelHt = fit_gam("REML GAM RelHt", DBH ~ s(TotalHt, relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 9, pc = gamConstraint), data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  
  save(file = "data/rhpu preferred models.Rdata", rhpuHeightFromDiameterPreferred, rhpuDiameterFromHeightPreferred)
}


## basal area from height
if (htDiaOptions$includeInvestigatory)
{
  rhpuBasalAreaFromHeightKorf = gsl_nls(basalArea ~ a1*(exp(b1*(imputedHeight - 1.37)^b2) - 1), rhpu2016, start = list(a1 = 90, b1 = 0.000003, b2 = 2.18), weights = heightWeight^2) # a1p, b1p, b2p not significant
  rhpuBasalAreaFromHeightPower = gsl_nls(basalArea ~ a1*(imputedHeight - 1.37)^(b1 + b1p * isPlantation), rhpu2016, start = list(a1 = 3/7 * 0.25 * pi * 0.01^2, b1 = 2.14, b1p = 0.34), weights = heightWeight^2) # a1p not significant
  #confint2(rhpuBasalAreaFromHeightPower, level = 0.99)
  
  tribble(~method, ~aic, ~biasCm2, ~maeM2, ~nse,
          "Korf", AIC(rhpuBasalAreaFromHeightKorf), 100^2 * mean(residuals(rhpuBasalAreaFromHeightKorf)), mean(abs(residuals(rhpuBasalAreaFromHeightKorf))), 1 - sum(residuals(rhpuBasalAreaFromHeightKorf)^2) / sum((rhpu2016$basalArea - mean(rhpu2016$basalArea)^2)),
          "power", AIC(rhpuBasalAreaFromHeightPower), 100^2 * mean(residuals(rhpuBasalAreaFromHeightPower)), mean(abs(residuals(rhpuBasalAreaFromHeightPower))), 1 - sum(residuals(rhpuBasalAreaFromHeightPower)^2) / sum((rhpu2016$basalArea - mean(rhpu2016$basalArea)^2))) %>%
    mutate(deltaAIC = aic - min(aic)) %>%
    arrange(desc(deltaAIC))
  
  ggplot(rhpu2016) +
    geom_point(aes(x = imputedHeight, y = 0.25*pi*(0.01*DBH)^2), alpha = 0.1, color = "grey25", shape = 16) +
    geom_line(aes(x = imputedHeight, y = predict(rhpuBasalAreaFromHeightKorf), color = "Korf", group = isPlantation)) +
    geom_line(aes(x = imputedHeight, y = predict(rhpuBasalAreaFromHeightPower), color = "power", group = isPlantation)) +
    #geom_path(aes(x = imputedHeight, y = 10*(1 - exp(-0.1*(imputedHeight - 1.37)))^1.2, color = "Chapman-Richards")) +
    labs(x = "cascara buckthorn height, m", y = "basal area, m²", color = NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
}


## exploratory plots
if (htDiaOptions$includeInvestigatory)
{
  library(GGally)
  ggpairs(rhpu2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger, relativeHeight, regeneration), 
          aes(alpha = 0.1, color = regeneration, shape = "16"),
          columnLabels = c("DBH, cm", "height, m", "BA, m² ha⁻¹", "BAL, m² ha⁻¹", "relative height, %", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    #scale_color_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) + # https://github.com/ggobi/ggally/issues/445
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  ggpairs(rhpu2016 %>% mutate(regeneration = if_else(isPlantation, "plantation", "natural regen")) %>% select(TotalHt, DBH, slope, elevation, topographicShelterIndex, regeneration), 
          aes(alpha = 0.1, color = if_else(rhpu2016$isPlantation, "plantation", "natural regen"), shape = "16"), 
          columnLabels = c("DBH, cm", "height, m", "slope, °", "elevation, m", "TSI, °", "stand type"),
          upper = list(continuous = wrap("cor", size = 3)),
          lower = list(combo = wrap("facethist", bins = 30))) +
    scale_color_discrete(type = c("forestgreen", "darkviolet")) +
    scale_fill_manual(breaks = c("natural regen", "plantation"), values = c("forestgreen", "darkviolet")) +
    theme(strip.background = element_blank())
  scatterPlotMatrix::scatterPlotMatrix(rhpu2016 %>% select(TotalHt, DBH, standBasalAreaPerHectare, basalAreaLarger))
}


## GAM smooth effects
if (htDiaOptions$includeInvestigatory)
{
  #rhpuInteraction = lm(TotalHt ~ DBH*standBasalAreaPerHectare + DBH:basalAreaLarger + standBasalAreaPerHectare:basalAreaLarger, rhpu2016)
  #summary(rhpuInteraction)
  #ggplot() +
  #  geom_point(aes(x = DBH, y = basalAreaLarger, color = rhpuInteraction$residuals), rhpu2016, shape = 16) +
  #  labs(x = "DBH, cm", y = bquote("BAL, m"^2*" ha"^-1), color = "height\nresidual, m") +
  #  scale_color_scico(palette = "bam", limits = c(-20, 20))
  rhpuHeightGam = fit_gam("REML GAM", TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint) + 
                            #s(standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint) + # not significant
                            #s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint) + # not significant
                            #s(elevation, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                            #s(slope, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                            #s(aspect, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                            s(topographicShelterIndex, bs = "ts", k = 5, pc = gamConstraint) + 
                            s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint), 
                          data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  #rhpuHeightGam = fit_gam("REML GAM", TotalHt ~ s(DBH, standBasalAreaPerHectare, bs = "ts", by = as.factor(isPlantation), k = 11, pc = gamConstraint) + 
  #                          #s(basalAreaLarger, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint) + # not significant
  #                          #s(elevation, bs = "ts", k = 3, pc = gamConstraint) + # not significant
  #                          #s(slope, bs = "ts", k = 3, pc = gamConstraint) + # not significant
  #                          #s(aspect, bs = "ts", k = 3, pc = gamConstraint) + # not significant
  #                          s(topographicShelterIndex, bs = "ts", k = 5, pc = gamConstraint),
  #                          #s(relativeDiameter, bs = "ts", by = as.factor(isPlantation), k = 4, pc = gamConstraint), # not significant
  #                        data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  k.check(rhpuHeightGam)
  summary(rhpuHeightGam)
  par(mfrow = c(2, 3), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(rhpuHeightGam, scale = 0, scheme = 2)
  
  rhpuDbhGam = fit_gam("REML GAM", DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 8, pc = gamConstraint),
                       #s(standBasalAreaApprox, bs = "ts", by = as.factor(isPlantation), k = 3, pc = gamConstraint) + # not significant
                       #s(tallerApproxBasalArea, bs = "ts", by = as.factor(isPlantation), k = 3, pc = gamConstraint) + # not significant
                       #s(elevation, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                       #s(slope, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                       #s(aspect, bs = "ts", k = 3, pc = gamConstraint) + # not significant
                       #s(topographicShelterIndex, bs = "ts", k = 3, pc = gamConstraint), # not significant
                       #s(relativeHeight, bs = "ts", by = as.factor(isPlantation), k = 3, pc = gamConstraint), # not significant
                       data = rhpu2016, constraint = rhpu2016gamConstraint, folds = 1, repetitions = 1)
  k.check(rhpuDbhGam)
  summary(rhpuDbhGam)
  par(mfrow = c(1, 4), mar = c(2.2, 2.2, 0.5, 0) + 0.1, mgp = c(1.5, 0.4, 0))
  plot.gam(rhpuDbhGam, scale = 0)
}


## random forest regression
if (htDiaOptions$includeInvestigatory)
{
  library(caret)
  library(ranger)
  #rhpuForest = ranger(TotalHt ~ DBH + standBasalAreaPerHectare + basalAreaLarger, rhpu2016, classification = TRUE, num.threads = 12)
  repeatedCrossValidation = trainControl(method = "repeatedcv", number = htDiaOptions$folds, repeats = htDiaOptions$repetitions, verboseIter = FALSE)
  rhpuHeightForest = train(TotalHt ~ DBH + standBasalAreaPerHectare + basalAreaLarger + elevation + slope + aspect + topographicShelterIndex + relativeDiameter, data = rhpu2016, method = "ranger", trControl = repeatedCrossValidation, 
                           importance = "impurity_corrected",
                           tuneGrid = expand.grid(mtry = c(6, 8),
                                                  splitrule = "variance",
                                                  min.node.size = c(1, 2)))
  rhpuHeightForest
  varImp(rhpuHeightForest)
  
  rhpuDbhForest = train(DBH ~ TotalHt + standBasalAreaApprox + tallerApproxBasalArea + elevation + slope + aspect + topographicShelterIndex + relativeHeight, data = rhpu2016, method = "ranger", trControl = repeatedCrossValidation, 
                        importance = "impurity_corrected",
                        tuneGrid = expand.grid(mtry = c(7, 8),
                                               splitrule = "variance",
                                               min.node.size = c(2, 3, 4)))
  rhpuDbhForest
  varImp(rhpuDbhForest)
}

