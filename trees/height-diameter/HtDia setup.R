library(dplyr)
library(ggplot2)
library(gslnls)
library(leaps)
library(magrittr)
library(nlme)
library(nls.multstart)
library(nlstools)
library(patchwork)
library(readxl)
library(robustbase)
library(stringr)
library(tibble)
library(tidyr)
library(writexl)

theme_set(theme_bw() + theme(axis.line = element_line(size = 0.3),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             legend.key.height = unit(0.85, "line"),
                             legend.spacing.y = unit(0, "line"),
                             panel.border = element_blank()))

as_row = function(regression)
{
  if (is.null(regression))
  {
    return(rep(NA_real_, 21))
  }
  
  #return(c(pae = regression$pae, paeNR = regression$paeNaturalRegen, paePl = regression$paePlantation,
  #         bias = regression$bias, biasNR = regression$biasNaturalRegen, biasPl = regression$biasPlantation,
  #         mae = regression$mae, maeNR = regression$maeNaturalRegen, maePl = regression$maePlantation,
  #         rmse = regression$rmse, rmseNR = regression$rmseNaturalRegen, rmsePl = regression$rmsePlantation,
  #         nse = regression$nse, nseNR = regression$nseNaturalRegen, nsePl = regression$nsePlantation,
  #         pearson = regression$pearson, pearsonNR = regression$pearsonNaturalRegen, pearsonPL = regression$pearsonPlantation,
  #         aic = regression$aic, bic = regression$bic))
  # omit names for compatibility with !!! operator
  power = NA_real_
  if (is.null(regression$modelStruct$varStruct) == FALSE)
  {
    power = regression$modelStruct$varStruct
  }
  
  row = c(regression$pae, regression$paeNaturalRegen, regression$paePlantation,
          regression$bias, regression$biasNaturalRegen, regression$biasPlantation,
          regression$mae, regression$maeNaturalRegen, regression$maePlantation,
          regression$rmse, regression$rmseNaturalRegen, regression$rmsePlantation,
          regression$nse, regression$nseNaturalRegen, regression$nsePlantation,
          regression$pearson, regression$pearsonNaturalRegen, regression$pearsonPlantation,
          regression$aic, regression$bic, power)
  if (length(row) != 21)
  {
    stop("Did not find all 21 expected summary fields attached to regression.")
  }
  return(row)
}

get_dbh_error = function(regression, data, dataNaturalRegen, dataPlantation)
{
  regression$fitted.values = predict(regression, data)
  regression$residuals = data$DBH - regression$fitted.values

  regression$aic = AIC(regression)
  regression$bias = mean(regression$residuals)
  regression$bic = BIC(regression)
  regression$mae = mean(abs(regression$residuals))
  regression$nse = 1 - sum(regression$residuals^2) / sum((data$DBH - mean(data$DBH))^2)
  regression$pae = 100 * mean(abs(regression$residuals / data$DBH))
  regression$pearson = cor(regression$fitted.values, data$DBH)
  regression$rmse = sqrt(mean(regression$residuals^2))
  
  dbhNaturalRegen = predict(regression, dataNaturalRegen)
  residualsNaturalRegen = dbhNaturalRegen - dataNaturalRegen$DBH
  regression$biasNaturalRegen = mean(residualsNaturalRegen)
  regression$maeNaturalRegen = mean(abs(residualsNaturalRegen))
  regression$nseNaturalRegen = 1 - sum(residualsNaturalRegen^2) / sum((dataNaturalRegen$DBH - mean(dataNaturalRegen$DBH))^2)
  regression$paeNaturalRegen = 100 * mean(abs(residualsNaturalRegen / dataNaturalRegen$DBH))
  regression$pearsonNaturalRegen = cor(residualsNaturalRegen, dataNaturalRegen$DBH)
  regression$rmseNaturalRegen = sqrt(mean(residualsNaturalRegen^2))
  
  dbhPlantation = predict(regression, dataPlantation)
  residualsPlantation = dbhPlantation - dataPlantation$DBH
  regression$biasPlantation = mean(residualsPlantation)
  regression$maePlantation = mean(abs(residualsPlantation))
  regression$nsePlantation = 1 - sum(residualsPlantation^2) / sum((dataPlantation$DBH - mean(dataPlantation$DBH))^2)
  regression$paePlantation = 100 * mean(abs(residualsPlantation / dataPlantation$DBH))
  regression$pearsonPlantation = cor(residualsPlantation, dataPlantation$DBH)
  regression$rmsePlantation = sqrt(mean(residualsPlantation^2))
  
  return(regression)
}

get_height_error = function(regression, data, dataNaturalRegen, dataPlantation)
{
  regression$fitted.values = predict(regression, data)
  regression$residuals = data$TotalHt - regression$fitted.values
  
  regression$aic = AIC(regression)
  regression$bias = mean(regression$residuals)
  regression$bic = BIC(regression)
  regression$mae = mean(abs(regression$residuals))
  regression$nse = 1 - sum(regression$residuals^2) / sum((data$TotalHt - mean(data$TotalHt))^2)
  regression$pae = 100 * mean(abs(regression$residuals / data$TotalHt))
  regression$pearson = cor(regression$fitted.values, data$TotalHt)
  regression$rmse = sqrt(mean(regression$residuals^2))
  
  heightNaturalRegen = predict(regression, dataNaturalRegen)
  residualsNaturalRegen = heightNaturalRegen - dataNaturalRegen$TotalHt
  regression$biasNaturalRegen = mean(residualsNaturalRegen)
  regression$maeNaturalRegen = mean(abs(residualsNaturalRegen))
  regression$nseNaturalRegen = 1 - sum(residualsNaturalRegen^2) / sum((dataNaturalRegen$TotalHt - mean(dataNaturalRegen$TotalHt))^2)
  regression$paeNaturalRegen = 100 * mean(abs(residualsNaturalRegen / dataNaturalRegen$TotalHt))
  regression$pearsonNaturalRegen = cor(residualsNaturalRegen, dataNaturalRegen$TotalHt)
  regression$rmseNaturalRegen = sqrt(mean(residualsNaturalRegen^2))
  
  heightPlantation = predict(regression, dataPlantation)
  residualsPlantation = heightPlantation - dataPlantation$TotalHt
  regression$biasPlantation = mean(residualsPlantation)
  regression$maePlantation = mean(abs(residualsPlantation))
  regression$nsePlantation = 1 - sum(residualsPlantation^2) / sum((dataPlantation$TotalHt - mean(dataPlantation$TotalHt))^2)
  regression$paePlantation = 100 * mean(abs(residualsPlantation / dataPlantation$TotalHt))
  regression$pearsonPlantation = cor(residualsPlantation, dataPlantation$TotalHt)
  regression$rmsePlantation = sqrt(mean(residualsPlantation^2))
  
  return(regression)
}

impute_height = function(Species, DBH, isPlantation)
{
  # for simplicity and for now, central Chapman-Richards fits for all species
  # Other forms have slightly lower error, see HtDia.xlsx.
  # switch fails with multi-argument returns not permitted
  return(recode(Species,
                DF = 1.37 + (65.30943 - 13.13382 * isPlantation) * (1 - exp(-0.02209 * DBH))^(1.50887 - 0.31001 * isPlantation),
                RA = 1.37 + (24.70557 + 2.93605 * isPlantation) * (1 - exp(-0.04842 * DBH))^(1.17287 - 0.07762 * isPlantation),
                WH = 1.37 + (50.20361 - 6.13208 * isPlantation) * (1 - exp(-0.0260 * DBH))^(1.43803 - 0.23284 * isPlantation),
                BM = 1.37 + (26.01242 + 1.31902 * isPlantation) * (1 - exp(-0.03156 * DBH))^(1.03369 - 0.09901 * isPlantation),
                OM = 1.37 + (17.65233 - 1.52591 * isPlantation) * (1 - exp(-0.05721 * DBH))^(1.23109 - 0.25308 * isPlantation),
                RC = 1.37 + (53.09458 - 8.88290 * isPlantation) * (1 - exp(-0.01434 * DBH))^(1.22216 - 0.17121 * isPlantation),
                .default = 1.37 + (57.89105 - 9.52781 * isPlantation) * (1 - exp(-0.01026  * DBH))^(0.93663 - 0.12157 * isPlantation)))
}

## load data
stands2022 = read_xlsx("GIS/Planning/Elliott Stand Data Feb2022.xlsx") %>% 
  mutate(Cruised_Si = na_if(Cruised_Si, 0),
         ODSL_Site_ = na_if(ODSL_Site_, 0),
         siteSpecies = if_else(startsWith(ODSL_VEG_L, "1W") | startsWith(ODSL_VEG_L, "WX"), "hemlock", 
                               if_else(startsWith(ODSL_VEG_L, "1H") | startsWith(ODSL_VEG_L, "HX"), "hardwood",
                                       if_else(startsWith(ODSL_VEG_L, "OT"), "other",
                                               "Douglas-fir"))))

plots2016 = read_xlsx("GIS/Trees/2015-16 cruise/CruisePlots_All_20151211.xlsx")

trees2016 = left_join(left_join(read_xlsx("trees/Elliott final cruise records 2015-16.xlsx", sheet = "CRUISERECS"),
                                          stands2022 %>% select(StandID, GrossAc, Age_2020, Cruised_Si, Elev_Mean, SlopeMean, AspectSin, AspectCos),
                                          by = c("StandID")),
                                plots2016 %>% select(STAND, PltInteger, elevation, slope, aspect, topographicShelterIndex) %>% rename(PlotID = PltInteger),
                                by = c("PlotID")) %>%
  rename(standAge2020 = Age_2020) %>%
  mutate(isLiveUnbroken = (CompCode %in% c("BT", "D.", "SN")) == FALSE,
         BHAge = na_if(BHAge, 0), # years
         DBH = na_if(2.54 * DBH, 0), # inches to cm
         Dia1 = na_if(2.54 * Dia1, 0),
         Ht1 = na_if(0.3048 * Ht1, 0), # feet to m
         Ht2 = na_if(0.3048 * Ht2, 0),
         TotalHt = na_if(0.3048 * TotalHt, 0),
         basalArea = 0.25 * pi * (0.01*DBH)^2, # m² 
         breastHeight = 1.37,
         isPlantation = standAge2020 < 75,
         SampleFactor = 2.47105 * SampleFactor, # trees per acre to trees per hectare
         standArea = 0.404686 * GrossAc,  # ac to ha
         standSampleFactor = mean(SampleFactor),
         treeBasalAreaPerHectare = SampleFactor * if_else(SamplingMethod == "BAF", 0.092903, basalArea), # m²/ha, conversion factor is either 2.47105 * 0.092903 = 0.229568 m²/ha / ft²/ac or 2.47105 ac/ha
         heightDiameterRatio = TotalHt / (0.01 * DBH), # (DBH conversion from cm to m)
         imputedHeight = if_else(is.na(TotalHt) == FALSE, TotalHt, if_else(is.na(DBH) == FALSE, impute_height(Species, DBH, isPlantation), NA_real_)), # where possible, perform basic height imputation
         quasiBasalArea = SampleFactor * if_else(SamplingMethod == "BAF", 0.092903, 2.617E-05*imputedHeight^2.543)) %>% # stack basal area regression on height regression when possible; BAF has to be used with prism trees
  select(-GrossAc) %>%
  group_by(PlotID) %>%
  mutate(plotTrees = sum(SampleFactor),
         plotTreesWithDbh = sum(if_else(is.na(DBH), 0, SampleFactor)),
         plotContributionToStandBasalArea = 0,
         plotContributionToStandBasalArea = replace(plotContributionToStandBasalArea, 1, sum(treeBasalAreaPerHectare)),
         plotContributionToStandQuasiBasalArea = 0,
         plotContributionToStandQuasiBasalArea = replace(plotContributionToStandQuasiBasalArea, 1, sum(quasiBasalArea))) %>%
  group_by(StandID) %>%
  arrange(desc(isLiveUnbroken), desc(DBH), .by_group = TRUE) %>%
  mutate(plotsInStand = length(unique(PlotID)),
         standBasalAreaPerHectare = sum(plotContributionToStandBasalArea) / plotsInStand, # m²/ha
         standQuasiBasalArea = sum(plotContributionToStandQuasiBasalArea) / plotsInStand,
         basalAreaLarger = (cumsum(isLiveUnbroken * treeBasalAreaPerHectare) - treeBasalAreaPerHectare[1]) / plotsInStand, # m²/ha
         tph = sum(SampleFactor) / plotsInStand) %>%
  arrange(desc(isLiveUnbroken), desc(imputedHeight), .by_group = TRUE) %>% # put tallest live trees without broken tops first in each stand (numbers sort before NA)
  mutate(topHeight = mean(if_else(row_number() < 100 * standArea / standSampleFactor, TotalHt, NA_real_), na.rm = TRUE), # m, tallest 100 trees per hectare
         topHeight = if_else(is.na(topHeight), mean(if_else(row_number() < 100 * standArea / standSampleFactor, imputedHeight, NA_real_), na.rm = TRUE), topHeight), # fall back to imputed heights if no trees in stand were measured for height
         relativeHeight = TotalHt / topHeight, # individual trees' heights as a fraction of top height, may be greater than 1, especially for retention trees
         tallerQuasiBasalArea = (cumsum(isLiveUnbroken * quasiBasalArea) - quasiBasalArea[1]) / plotsInStand,
         tallerTph = cumsum(isLiveUnbroken * SampleFactor) / plotsInStand) %>% 
  ungroup()

#print(trees2016 %>% filter(is.na(elevation)) %>% group_by(PlotID) %>% summarize(trees = n(), .groups = "drop"), n = 51)
#ggplot(trees2016) + geom_histogram(aes(x = standQuasiBasalArea))

## species unpacking
psme2016 = trees2016 %>% filter(Species == "DF", isLiveUnbroken, TotalHt > 0) # live Douglas-firs measured for height
psme2016natural = psme2016 %>% filter(isPlantation == FALSE)
psme2016plantation = psme2016 %>% filter(isPlantation)

alru2016 = trees2016 %>% filter(Species == "RA", isLiveUnbroken, TotalHt > 0) # live red alders measured for height
alru2016natural = alru2016 %>% filter(isPlantation == FALSE)
alru2016plantation = alru2016 %>% filter(isPlantation)

tshe2016 = trees2016 %>% filter(Species == "WH", isLiveUnbroken, TotalHt > 0) # live western hemlocks measured for height
tshe2016natural = tshe2016 %>% filter(isPlantation == FALSE)
tshe2016plantation = tshe2016 %>% filter(isPlantation)

acma2016 = trees2016 %>% filter(Species == "BM", isLiveUnbroken, TotalHt > 0) # live bigleaf maples measured for height
acma2016natural = acma2016 %>% filter(isPlantation == FALSE)
acma2016plantation = acma2016 %>% filter(isPlantation)

umca2016 = trees2016 %>% filter(Species == "OM", isLiveUnbroken, TotalHt > 0) # live Oregon myrtles measured for height
umca2016natural = umca2016 %>% filter(isPlantation == FALSE)
umca2016plantation = umca2016 %>% filter(isPlantation)

thpl2016 = trees2016 %>% filter(Species == "RC", isLiveUnbroken, TotalHt > 0) # live western redcedars measured for height
thpl2016natural = thpl2016 %>% filter(isPlantation == FALSE)
thpl2016plantation = thpl2016 %>% filter(isPlantation)

other2016 = trees2016 %>% filter((Species %in% c("DF", "RA", "WH", "BM", "OM", "RC")) == FALSE, isLiveUnbroken, TotalHt > 0) # live western redcedars measured for height
other2016natural = other2016 %>% filter(isPlantation == FALSE)
other2016plantation = other2016 %>% filter(isPlantation)

#otherConifer2016 = other2016 %>% filter(Species %in% c("XX", "CX", "SS", "PC", "PY", "GF", "LP"))
#otherHardwood2016 = other2016 %>% filter(Species %in% c("CA", "HX", "CH", "PM", "GC", "PD", "TO", "WI", "OA", "WO"))


## data tabulation and basic plotting
trees2016summary = trees2016 %>% 
  mutate(isLive = CompCode %in% c("D.", "SN")) %>%
  #mutate(speciesClassification = if_else(Species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), Species, "other")) %>%
  #group_by(speciesClassification) %>%
  group_by(Species) %>% 
  summarize(percentage = 100 * n() / nrow(trees2016), 
            trees = n(),
            live = sum(isLive), plantation = sum(isPlantation), retention = sum(CompCode == "RT"), snag = sum(isLive == FALSE), 
            dbh = sum(DBH > 0, na.rm = TRUE), 
            height = sum(TotalHt > 0, na.rm = TRUE), 
            age = sum(BHAge > 0, na.rm = TRUE), 
            crownRatio = sum(CrownRatio > 0, na.rm = TRUE), 
            dia1 = sum(Dia1 > 0, na.rm = TRUE), height1 = sum(Ht1 > 0, na.rm = TRUE), 
            height2 = sum(Ht2 > 0, na.rm = TRUE), 
            .groups = "drop") %>%
  arrange(desc(trees))
print(trees2016summary, n = 25)
trees2016 %>% filter((CompCode %in% c("D.", "SN")) == FALSE, DBH > 0, Ht2 > 0) %>% summarize(n = n())
trees2016 %>% filter((CompCode %in% c("D.", "SN")) == FALSE, TotalHt > 0 | Ht2 > 0) %>% summarize(n = n())
trees2016 %>% filter(CompCode == OS) %>% group_by(Stand) summarize(n = n())
print(trees2016 %>% filter(CompCode == "RT", isPlantation == FALSE) %>% select(StandID, Species, DBH, standAge2020), n = 35)


trees2016 %>% group_by(StandID) %>% summarize(standArea = standArea[1], tph = tph[1]) %>%
  summarize(trees = sum(standArea * tph))

ggplot(trees2016 %>% filter(is.na(Ht1) == FALSE)) +
  geom_histogram(aes(y = 100 * Ht1 / TotalHt, x = 100 * ..count.. / sum(..count..)), binwidth = 1) +
  labs(x = "percentage of trees, %", y = "taper measurement's relative height, %")

ggplot(trees2016 %>% filter(CompCode == "OS", BHAge > 0) %>% group_by(StandID) %>% summarize(siteTrees = n())) +
  geom_histogram(aes(x = siteTrees, y = 100 * ..count.. / sum(..count..)), binwidth = 1) +
  labs(x = "site trees", y = "percentage of stands") +
  scale_x_continuous(breaks = seq(1, 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))

## Douglas-fir site index regression: not enough data for other species
# site species  number of stands
# PSME          412
# hardwood      26
# hemlock       6
# other         5
psmeStands2022 = stands2022 %>% filter(is.na(Cruised_Si) == FALSE, siteSpecies == "Douglas-fir") %>%
  mutate(Elev_MeanSquared = Elev_Mean^2, SlopeMeanPercent = 100 * tan(pi/180 * SlopeMean), SlopeMeanSquared = SlopeMean^2, SlopeMeanPercentSquared = SlopeMeanPercent^2, AWS100squared = AWS100^2, planted = Age_2020 < 100)

psmeSiteIndexPredictorSubsets = regsubsets(x = as.matrix(psmeStands2022 %>% select(Elev_Mean, Elev_MeanSquared, SlopeMeanPercent, SlopeMeanPercentSquared, AspectMean, AspectSin, AspectCos, 
                                                                                   PrecipNorm, AWS025, AWS050, AWS100, AWS100squared, AWS150, planted, X, Y,
                                                                                   Age_2015, TPA_Total, BA_Total, QMD_Total, BFperAcre_, BA_DF, QMD_DF, LeafRetnDF, BA_WH, QMD_WH, Shape_Area)),
                                           y = psmeStands2022$Cruised_Si)
plot(psmeSiteIndexPredictorSubsets)

#psmeSiteIndexModelLinear = lm(Cruised_Si ~ Elev_Mean + Elev_MeanSquared + SlopeMeanPercent + SlopeMeanPercentSquared + AspectSin + AspectCos + PrecipNorm + AWS100 + TPA_Total + BA_Total + QMD_Total + QMD_DF + LeafRetnDF + planted, psmeStands2022)
#psmeSiteIndexModelLinear = lm(Cruised_Si ~ Elev_MeanSquared + SlopeMeanPercent + SlopeMeanPercentSquared + AspectSin + AspectCos + AWS100 + planted, psmeStands2022)
psmeSiteIndexModelLinear = lm(Cruised_Si ~ Elev_MeanSquared + SlopeMeanPercentSquared + QMD_DF + planted, psmeStands2022)
summary(psmeSiteIndexModelLinear)
psmeSiteIndexModelNonlinear = nls(Cruised_Si ~ b0 + b1*Elev_Mean^b2 + b3*SlopeMeanPercent^b4 + b5*planted, psmeStands2022, start = list(b0 = 120, b1 = -1E-6, b2 = 2, b3 = -2E-3, b4 = 5, b5 = 20), control = list(maxiter = 250))
c(linear = AIC(psmeSiteIndexModelLinear), nonlinear = AIC(psmeSiteIndexModelNonlinear))

ggplot(psmeStands2022) + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "longdash") + 
  geom_point(aes(x = Cruised_Si, y = psmeSiteIndexModelLinear$fitted.values, color = planted), alpha = 0.3) +
  labs(x = "measured 50-year site index, feet", y = "linear model prediction, feet", color = NULL) +
  theme(legend.position = "none") +
ggplot(psmeStands2022) + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "longdash") + 
  geom_point(aes(x = Cruised_Si, y = predict(psmeSiteIndexModelNonlinear, psmeStands2022), color = planted), alpha = 0.3) +
  labs(x = "measured 50-year site index, feet", y = "nonlinear model prediction, feet", color = "stand age") +
  scale_color_discrete(breaks = c(FALSE, TRUE), labels = c("≥100 years", "<100 years")) +
  theme(legend.justification = c(1, 0), legend.position = c(0.98, 0.02))

ggplot(psmeStands2022) +
  geom_point(aes(x = Cruised_Si, y = -psmeSiteIndexModelLinear$residuals, color = planted), alpha = 0.3, shape = 16) +
  labs(x = "measured 50-year site index, feet", y = "linear model error, feet", color = NULL) +
  theme(legend.position = "none") +
ggplot(psmeStands2022) +
  geom_point(aes(x = Cruised_Si, y = predict(psmeSiteIndexModelNonlinear, psmeStands2022) - Cruised_Si, color = planted), alpha = 0.3, shape = 16) +
  labs(x = "measured 50-year site index, feet", y = "nonlinear model error, feet", color = "stand age") +
  scale_color_discrete(breaks = c(FALSE, TRUE), labels = c("≥100 years", "<100 years")) +
  theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.98))

ggplot(psmeStands2022) + geom_point(aes(x = Elev_Mean, y = Cruised_Si), alpha = 0.3, shape = 16) +
ggplot(psmeStands2022) + geom_point(aes(x = SlopeMeanPercent, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = AspectSin, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = AspectCos, y = Cruised_Si), alpha = 0.3, shape = 16) +  labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = TPA_Total, y = Cruised_Si), alpha = 0.3, shape = 16) + 
ggplot(psmeStands2022) + geom_point(aes(x = BA_Total, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = QMD_Total, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) + 
ggplot(psmeStands2022) + geom_point(aes(x = BA_DF / BA_Total, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = PrecipNorm, y = Cruised_Si), alpha = 0.3, shape = 16) + 
ggplot(psmeStands2022) + geom_point(aes(x = AWS100, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) +
ggplot(psmeStands2022) + geom_point(aes(x = QMD_DF, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL) + 
ggplot(psmeStands2022) + geom_point(aes(x = QMD_WH, y = Cruised_Si), alpha = 0.3, shape = 16) + labs(y = NULL)
  

## aggregate tree distribution plots
liveTrees2016 = trees2016 %>% filter(live) %>% mutate(speciesGroup = factor(if_else(Species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), Species, "other"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other")))
ggplot(liveTrees2016) +
  geom_histogram(aes(x = DBH, y = 100 * ..count../sum(..count..), fill = speciesGroup, alpha = isPlantation), binwidth = 2.5, na.rm = TRUE) +
  coord_cartesian(ylim = c(0, 4.4)) +
  labs(x = "DBH, cm", y = "percentage of live stems measured", alpha = NULL, fill = NULL) +
  scale_alpha_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c(1, 0.7)) +
  scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), values = c("green3", "red2", "blue2", "cyan2", "darkorchid3", "firebrick", "grey35")) +
  theme(legend.position = "none") +
ggplot(liveTrees2016) +
  geom_histogram(aes(x = 100 * ..count../sum(..count..), y = TotalHt, fill = speciesGroup, alpha = isPlantation), binwidth = 1, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 4.4)) +
  labs(x = "percentage of live stems measured", y = "height, m", alpha = NULL, fill = NULL) +
  scale_alpha_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c(1, 0.7)) +
  scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other"), values = c("green3", "red2", "blue2", "cyan2", "darkorchid3", "firebrick", "grey35")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.spacing.y = unit(0.3, "line"))
  

## site index plots
ggplot(stands2022 %>% filter(Cruised_Si > 0)) +
  geom_point(aes(x = Age_2020, y = Cruised_Si, color = siteSpecies), alpha = 0.6, shape = 16) +
  labs(x = "stand age in 2020, years", y = "50-year site index measured in 2015-2016, feet", color = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 0.8))) +
  scale_color_manual(breaks = c("Douglas-fir", "hemlock", "hardwood", "other"), values = c("green4", "blue2", "gold1", "purple1")) +
  theme(legend.justification = c(1, 0), legend.position = c(0.98, 0.02))

ggplot() +
  geom_hline(yintercept = 75, color = "grey70", linetype = "longdash") +
  geom_hline(yintercept = 95, color = "grey70", linetype = "longdash") +
  geom_hline(yintercept = 115, color = "grey70", linetype = "longdash") +
  geom_hline(yintercept = 135, color = "grey70", linetype = "longdash") +
  geom_histogram(aes(y = Cruised_Si, weight = GrossAc), stands2022 %>% filter(Age_2020 < 100, Cruised_Si > 0, siteSpecies == "Douglas-fir"), binwidth = 1, color = "white", fill = "green4") +
  annotate("text", x = 1150, y = 145, label = "I", color = "grey70", size = 3) +
  annotate("text", x = 1150, y = 125, label = "II", color = "grey70", size = 3) +
  annotate("text", x = 1150, y = 105, label = "III", color = "grey70", size = 3) +
  annotate("text", x = 1150, y = 85, label = "IV", color = "grey70", size = 3) +
  annotate("text", x = 1150, y = 65, label = "V", color = "grey70", size = 3) +
  labs(x = "acres of Douglas-fir majority stands cruised winter 2015-2016", y = "Douglas-fir 50-year site index, feet", color = NULL) +
  scale_x_continuous(breaks = seq(0, 2000, by = 250)) +
  scale_y_continuous(breaks = seq(0, 200, by = 10))

ggplot() +
  geom_histogram(aes(x = 100 * ..count../sum(..count..), y = siteClass, weight = GrossAc), 
                 stands2022 %>% filter(Age_2020 < 100, Cruised_Si > 0, siteSpecies == "Douglas-fir") %>% mutate(siteClass = factor(if_else(Cruised_Si >= 135, 1, if_else(Cruised_Si >= 115, 2, if_else(Cruised_Si >= 95, 3, if_else(Cruised_Si >= 75, 4, 5)))), levels = seq(5, 1, by = -1), labels = c("V", "IV", "III", "II", "I"))),
                 fill = "green4", stat = "count") +
  labs(x = "percentage of Douglas-fir majority stand area cruised in 2015-2016", y = "Douglas-fir 50-year site class", color = NULL) +
  scale_x_continuous(breaks = seq(0, 100, by = 10))

ggplot(stands2022 %>% filter(siteSpecies == "Douglas-fir")) +
  geom_histogram(aes(x = Age_2020, y = ..density..), binwidth = 5, fill = "green4") +
  labs(x = "stand age in 2020, years", y = "probability", color = NULL)
