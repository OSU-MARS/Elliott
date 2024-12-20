# load data from setup.R, get regressions and summaries from <species>.R
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 7)

figureDpi = 500
speciesGroupColors = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")

#rm(psmeResults, alruResults, tsheResults, acmaResults, umcaResults, thplResults, otherResults, psmeCoefficients, alruCoefficients, tsheCoefficients, acmaCoefficients, umcaCoefficients, thplCoefficients, otherCoefficients)
if (exists("psmeResults") == FALSE) { load("trees/height-diameter/data/PSME results.Rdata") }
if (exists("alruResults") == FALSE) { load("trees/height-diameter/data/ALRU2 results.Rdata") }
if (exists("tsheResults") == FALSE) { load("trees/height-diameter/data/TSHE results.Rdata") }
if (exists("acmaResults") == FALSE) { load("trees/height-diameter/data/ACMA3 results.Rdata") }
if (exists("thplResults") == FALSE) { load("trees/height-diameter/data/THPL results.Rdata") }
if (exists("umcaResults") == FALSE) { load("trees/height-diameter/data/UMCA results.Rdata") }
if (exists("otherResults") == FALSE) { load("trees/height-diameter/data/other results.Rdata") }


## assemble results tibbles from individual species data
heightDiameterResults = bind_rows(psmeResults, alruResults, tsheResults, acmaResults,
                                  umcaResults, thplResults, otherResults) %>%
  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name)),
         species = factor(species, labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"), levels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other")),
         speciesFraction = recode(species, "Douglas-fir" = 0.750, "red alder" = 0.101, "western hemlock" = 0.056, "bigleaf maple" = 0.029, "Oregon myrtle" = 0.025, "western redcedar" = 0.013, "other species" = 0.017),
         isBaseForm = (str_detect(name, "Sharma-") == FALSE) & (str_detect(name, "ABA\\+T") == FALSE) & (str_detect(name, "BA\\+L") == FALSE) & (str_detect(name, "physio") == FALSE) & (str_detect(name, "RelDbh") == FALSE) & (str_detect(name, "RelHt") == FALSE),
         hasPhysio = str_detect(name, "physio"),
         hasStand = str_detect(name, "ABA\\+T") | str_detect(name, "BA\\+L"),
         hasRelative = str_detect(name, "RelDbh") | str_detect(name, "RelHt"),
         significant = as.logical(significant), # since R lacks NA_logical_ significant can end up being either of type double (0/1/NA_real_) or logical (TRUE/FALSE), standardize back to logical (TRUE/FALSE/NA)
         weighting = if_else(fitting %in% c("gnls", "nlrob"), "reweighted", "fixed weights"),
         sizeShapeAlpha = as.factor(if_else(significant == TRUE, weighting, "not significant"))) %>%
  group_by(fitSet, fixedWeight, responseVariable, species) %>%
  mutate(nFits = n(),
         deltaAicN = aic/nValidation - min(aic/nValidation, na.rm = TRUE)) %>% # ΔAIC within response variable and species, needed for AUCs and figures
  ungroup()
# report duplicate naming and fit failures
heightDiameterResults %>% group_by(fitSet, responseVariable, species, name) %>% summarize(n = n(), .groups = "drop") %>% filter(n != htDiaOptions$folds * htDiaOptions$repetitions)

heightDiameterCoefficients = left_join(bind_rows(psmeCoefficients, alruCoefficients, tsheCoefficients, acmaCoefficients, # ~28 s
                                                 umcaCoefficients, thplCoefficients, otherCoefficients) %>%
                                         mutate(species = factor(species, labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"), levels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other"))),
                                       heightDiameterResults %>% select(-fitting, -fixedWeight, -significant), # no need to join duplicate columns
                                       by = join_by(fitSet, responseVariable, species, name, repetition, fold)) %>%
  mutate(isConverged = as.logical(isConverged)) %>%
  select(-weighting, -sizeShapeAlpha, -nFits, -nTaperImplausible, -speciesFraction) %>%
  relocate(responseVariable, species, fitSet, fixedWeight, name, significant, isBaseForm, hasRelative, hasStand, hasPhysio, fitting, repetition, fold, nObservations, nValidation, fitTimeInS, isConverged, effectiveDegreesOfFreedom, nNonPhysical, mab, mapb, mae, mape, rmse, rmspe, aic, deltaAicN, nse, meanAbsolutePlantationEffect,	meanAbsolutePercentPlantationEffect, a0, a1, a1p, a2, a2p, a3, a3p, a4, a5, a6, a7, a8, a9, a9p, a10, a10p, b1, b1p, b2, b2p, b3, b3p, b4, b4p)
#write_xlsx(heightDiameterCoefficients %>% 
#             filter(fitSet == "primary", is.na(fixedWeight)) %>%
#             select(-baseName, -fitSet, -fixedWeight, -aict, -bic, -bict, -bias, -ends_with("NaturalRegen"), -ends_with("Plantation"), -adaptiveWeightFraction) %>% # drop diagnostic columns
#             select(-a0, -starts_with("s(")), # drop columns for GAM intercept and smooth coefficients since 1) they're not meaningful outside of mgcv and 2) doing so reduces the primary fit set's .xlsx from 30+ MB to 15.8 MB
#             # select(-starts_with("a1r"), -starts_with("a3r")), -starts_with("Har"), -X, -starts_with("Xr")), %>% # drop random effects if mixed models are included
#           "trees/height-diameter/data/Elliott height-diameter model coefficients.xlsx") # slow on full results (74+ MB .xlsx)


# extract subset of results for analysis
primaryResults = heightDiameterResults %>% 
  filter(fitSet == "primary", is.na(fixedWeight), fitting != "gnls",  # exclude fits from gnls()
         (responseVariable != "height") | (str_detect(name, "RelHt") == FALSE), # exclude height control forms using relative height
         (responseVariable != "DBH") | (str_detect(name, "BA\\+L") == FALSE)) # exclude diameter control forms using basal area
#print(primaryResults %>% group_by(fitSet, responseVariable) %>% reframe(n = n(), names = unique(name)), n = 60)
#primaryResults %>% group_by(fitSet, species) %>% summarize(deltaAicN = sum(is.na(deltaAicN)), mab = sum(is.na(mab)), mae = sum(is.na(mae)), nse = sum(is.na(nse)), rmse = sum(is.na(rmse)))

# rank model forms by estimated prediction ability (using AUC) for form selection
#          runtime, seconds
# workers  group_modify()   group_by() %>% group_split() %>% future_map(), including 5-6s startup latency
#  1       157              104
#  4                         30.8
#  6                         23.4
#  7                         20.6
#  8                         20.9
# 16                         20.4
with_progress({
  crossValidatedModelCount = primaryResults %>% group_by(responseVariable, species) %>% summarize(n = n_distinct(name), .groups = "drop")
  progressBar = progressor(steps = sum(crossValidatedModelCount$n))
  
  heightDiameterModelAucs = primaryResults %>%
    group_by(responseVariable, species, name) %>%
    group_split() %>%
    future_map_dfr(function(fitResults)
    {
      if ((nrow(fitResults) == 1) | all(is.na(fitResults$nse)))
      {
        # no distribution to compare to since this model has only a no fit result or wasn't cross validated
        progressBar(str_pad(paste(fitResults$responseVariable[1], fitResults$species[1], fitResults$name[1]), 60, "right"))
        return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1],
                      otherModelName = NA_character_, fitting = fitResults$fitting[1], isBaseForm = fitResults$isBaseForm[1], hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                      aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucNse = NA_real_, aucRmse = NA_real_,
                      speciesFraction = fitResults$speciesFraction[1]))
      }
      
      # get all other cross-validation results for this response variable and species
      # Assumes no names are shared across fittings in the results set.
      matchingFitResults = primaryResults %>% filter(responseVariable == fitResults$responseVariable[1], species == fitResults$species[1])
      matchingModelNames = unique(matchingFitResults$name)
      pairwiseAucs = bind_rows(lapply(matchingModelNames, function(otherModelName) 
      {
        otherFitResults = matchingFitResults %>% filter(name == otherModelName)
        
        if ((nrow(otherFitResults) == 1) | all(is.na(otherFitResults$nse)))
        {
          # also no distribution to compare to if the other model has only a no fit result or wasn't cross validated
          return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1],
                        otherModelName = otherModelName, fitting = fitResults$fitting[1], isBaseForm = fitResults$isBaseForm[1], hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                        aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucNse = NA_real_, aucRmse = NA_real_,
                        speciesFraction = fitResults$speciesFraction[1]))
        }
        
        # find AUCs: for species with few samples (e.g. THPL) mean absolute bias may be all NA, which case WeightedROC() errors
        # AUC is the probability a sample from one empirical distribution is greater than a sample from an another 
        # distribution (see AUC.R). If the labels are reversed, the probability changes to a sample from the first distribution
        # being less than a sample from the other distribution.
        # The labels here are thus set up to measure the probabilities goodness of fit metrics from fitResults are preferable
        # to those from otherFitResults.
        lowerIsBetter = factor(c(rep(0, nrow(fitResults)), rep(1, nrow(otherFitResults))), levels = c(0, 1))
        higherIsBetter = factor(c(rep(1, nrow(fitResults)), rep(0, nrow(otherFitResults))), levels = c(0, 1))
        #if ((length(lowerIsBetter) < 2) | (length(higherIsBetter) < 2) | (n_distinct(lowerIsBetter) < 2) | (n_distinct(higherIsBetter) < 2))
        #{
        #  stop(paste0("ROC label formation error with name = ", otherModelName, " for ", fitResults$species[1], " ", fitResults$responseVariable[1], ".  nrow(fitResults) = ", nrow(fitResults), ", nrow(otherFitResults) = ", nrow(otherFitResults), "."))
        #}
        
        # estimate AUC for bias based on what data is available, which is potentially none for species with few measurements
        availableMabData = tibble(guess = c(fitResults$mab, otherFitResults$mab), label = lowerIsBetter) %>% drop_na()
        aucMab = NA_real_
        if ((nrow(availableMabData) > 1) & (n_distinct(availableMabData$label) > 1)) # unlikely but possible that availableMabData ends up with a single row, also possible one set of fits has MAB values but the other does not
        {
          #if ((nrow(availableMabData) < 2) | (n_distinct(availableMabData$label) < 2))
          #{
          #  stop(paste0("MAB ROC label formation error with name = ", otherModelName, " for ", fitResults$species[1], " ", fitResults$responseVariable[1], ".  nrow(fitResults) = ", nrow(fitResults), ", nrow(otherFitResults) = ", nrow(otherFitResults), ", nrow(availableMabData) = ", nrow(availableMabData), "."))
          #}
          aucMab = WeightedAUC(WeightedROC(guess = availableMabData$guess, label = availableMabData$label))
        }
        
        if ((sum(is.na(fitResults$mae)) + sum(is.na(otherFitResults$mae)) + sum(is.na(fitResults$nse)) + sum(is.na(otherFitResults$nse)) + 
             sum(is.na(fitResults$rmse)) + sum(is.na(otherFitResults$rmse))) > 0)
        {
          # WeightedROC() errors on NAs in its arguments but can't do so informatively
          # Fail informatively instead so that investigation is possible.
          stop(paste0(fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " (", nrow(fitResults), ") x ", otherModelName, " (", nrow(otherFitResults), "):",
                      " mab ", sum(is.na(fitResults$mab)), " ", sum(is.na(otherFitResults$mab)),
                      " mae ", sum(is.na(fitResults$mae)), " ", sum(is.na(otherFitResults$mae)),
                      " nse ", sum(is.na(fitResults$nse)), " ", sum(is.na(otherFitResults$nse)),
                      " rmse ", sum(is.na(fitResults$rmse)), " ", sum(is.na(otherFitResults$rmse))), quote = FALSE)
        }
        
        deltaAicNguess = c(fitResults$deltaAicN, otherFitResults$deltaAicN)
        maeGuess = c(fitResults$mae, otherFitResults$mae)
        nseGuess = c(fitResults$nse, otherFitResults$nse)
        rmseGuess = c(fitResults$rmse, otherFitResults$rmse)
        expectedLength = nrow(fitResults) + nrow(otherFitResults)
        if ((length(deltaAicNguess) != expectedLength) | (length(maeGuess) != expectedLength) | (length(nseGuess) != expectedLength) | 
            (length(rmseGuess) != expectedLength))
        {
          stop(paste0(fitResults$species[1], " ", fitResults$responseVariable[1], " ", fitResults$name[1], " (", nrow(fitResults), ") x ", otherModelName, " (", nrow(otherFitResults), ") NULLs in data:", 
                      " AICn ", length(deltaAicNguess), 
                      " MAE ", length(maeGuess),
                      " NSE ", length(nseGuess),
                      " RMSE ", length(rmseGuess)))
        }
        return(tibble(responseVariable = fitResults$responseVariable[1], species = fitResults$species[1], name = fitResults$name[1], fitting = fitResults$fitting[1], significant = fitResults$significant[1], isBaseForm = fitResults$isBaseForm[1],
                      otherModelName = otherModelName, otherModelSignificant = otherFitResults$significant[1],
                      hasPhysio = fitResults$hasPhysio[1], hasStand = fitResults$hasStand[1], hasRelative = fitResults$hasRelative[1],
                      aucDeltaAicN = WeightedAUC(WeightedROC(guess = deltaAicNguess, label = lowerIsBetter)),
                      aucMab = aucMab,
                      aucMabN = nrow(availableMabData),
                      aucMae = WeightedAUC(WeightedROC(guess = maeGuess, label = lowerIsBetter)),
                      aucNse = WeightedAUC(WeightedROC(guess = nseGuess, label = higherIsBetter)),
                      aucRmse = WeightedAUC(WeightedROC(guess = rmseGuess, label = lowerIsBetter)),
                      speciesFraction = fitResults$speciesFraction[1]))
      }))
      
      progressBar(str_pad(paste(fitResults$responseVariable[1], fitResults$species[1], fitResults$name[1]), 60, "right"))
      return(pairwiseAucs)
    })
})
#heightDiameterModelAucs %>% group_by(responseVariable, species) %>%
#  summarize(matrixElements = n(), mabN = sum(is.na(aucMab) == FALSE), mabInputN = sum(aucMabN, na.rm = TRUE), maeN = sum(is.na(aucMae) == FALSE), rmseN = sum(is.na(aucRmse) == FALSE), nseN = sum(is.na(aucNse) == FALSE), .groups = "drop")

# take median AUCs over otherModelName, excluding self and unsuccessful fits
# Exclusion of unsuccessful fits is debatable. If a model could not be fit then its AUC could reasonably be taken to be
# zero rather than NA, implying all models which could be fit have an AUC of 1 in comparison.
heightDiameterModelRanking = heightDiameterModelAucs %>% 
  group_by(responseVariable, species, name) %>%
  summarize(fitting = fitting[1], isBaseForm = isBaseForm[1], hasPhysio = hasPhysio[1], hasStand = hasStand[1], hasRelative = hasRelative[1],
            aucDeltaAicN = median(if_else(name != otherModelName, aucDeltaAicN, NA_real_), na.rm = TRUE),
            aucMab = median(if_else(name != otherModelName, aucMab, NA_real_), na.rm = TRUE),
            aucMae = median(if_else(name != otherModelName, aucMae, NA_real_), na.rm = TRUE),
            aucNse = median(if_else(name != otherModelName, aucNse, NA_real_), na.rm = TRUE),
            aucRmse = median(if_else(name != otherModelName, aucNse, NA_real_), na.rm = TRUE),
            significant = significant[1],
            speciesFraction = speciesFraction[1],
            .groups = "drop_last") %>%
  group_by(responseVariable, species, isBaseForm) %>%
  mutate(aucBlended = 0.1 * aucMab + 0.3 * if_else(responseVariable == "height", aucMae, aucRmse) + 0.1 * if_else(responseVariable == "DBH", aucMae, aucRmse) + 0.4 * aucDeltaAicN + 0.1 * aucNse,
         aucDeltaAicNRank = min_rank(desc(round(aucDeltaAicN * significant, 3))), 
         aucMabRank = min_rank(desc(round(aucMab * significant, 3))), 
         aucMaeRank = min_rank(desc(round(aucMae * significant, 3))), 
         aucNseRank = min_rank(desc(round(aucNse * significant, 3))), 
         aucRmseRank = min_rank(desc(round(aucRmse * significant, 3)))) %>%
  ungroup()
heightDiameterModelDisplaySort = heightDiameterModelRanking %>%
  group_by(responseVariable, name) %>%
  summarize(penalizedBlendedAuc = sum(speciesFraction * if_else(is.na(aucBlended), 0, aucBlended)), .groups = "drop") %>% # can also weight non-significant fits differently
  arrange(responseVariable, penalizedBlendedAuc)
#heightDiameterModelRanking %>% group_by(responseVariable, species) %>%
#  summarize(n = n(), mabN = sum(is.na(aucMab) == FALSE), maeN = sum(is.na(aucMae) == FALSE), rmseN = sum(is.na(aucRmse) == FALSE), nseN = sum(is.na(aucNse) == FALSE), .groups = "drop")

# find preferred model forms
nPreferredModelForms = 4
preferredForms = full_join(full_join(full_join(heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMab, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucMab)) %>% mutate(mabName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, mabName, aucMab),
                                               heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMae, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucMae)) %>% mutate(maeName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, maeName, aucMae),
                                               by = c("responseVariable", "species", "isBaseForm", "rank")), # full join since western redcedar mean absolute bias is all NA
                                     full_join(heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucRmse, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucRmse)) %>% mutate(rmseName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, rmseName, aucRmse),
                                               heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucDeltaAicN, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucDeltaAicN)) %>% mutate(aicName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, aicName, aucDeltaAicN),
                                               by = c("responseVariable", "species", "isBaseForm", "rank")),
                                     by = c("responseVariable", "species", "isBaseForm", "rank")),
                           full_join(heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucNse, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucNse)) %>% mutate(nseName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, nseName, aucNse),
                                     heightDiameterModelRanking %>% filter(significant) %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucBlended, n = nPreferredModelForms, na_rm = TRUE) %>% arrange(desc(aucBlended)) %>% mutate(blendedName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, blendedName, aucBlended),
                                     by = c("responseVariable", "species", "isBaseForm", "rank")),
                           by = c("responseVariable", "species", "isBaseForm", "rank")) %>%
  arrange(desc(responseVariable), species, desc(isBaseForm), rank) %>%
  relocate(responseVariable, species, isBaseForm, rank) %>%
  ungroup()
#preferredForms %>% group_by(responseVariable) %>% mutate(speciesPresent = n_distinct(species)) %>% group_by(responseVariable, species) %>% summarize(speciesPresent = speciesPresent[1], nPreferredForms = n(), .groups = "drop")

# summarize predictor variable selection
predictorVariableResults = heightDiameterCoefficients %>% 
  filter(fitSet == "primary", is.na(fixedWeight), fitting %in% c("gsl_nls", "nlrob")) %>% # exclude GAMs and linear controls
  mutate(hasBasalArea = (is.na(a2) == FALSE) | (is.na(a3) == FALSE) | ((responseVariable == "height") & str_detect(name, "(Sharma-Parton|Sharma-Zhang)")),
         hasPhysio = (is.na(a4) == FALSE) | (is.na(a5) == FALSE) | (is.na(a6) == FALSE) | (is.na(a7) == FALSE) | (is.na(a8) == FALSE),
         hasPlantation = (is.na(a1p) == FALSE) | (is.na(a2p) == FALSE) | (is.na(a3p) == FALSE) | (is.na(a9p) == FALSE) | (is.na(a10p) == FALSE) | (is.na(b1p) == FALSE) | (is.na(b2p) == FALSE) | (is.na(b3p) == FALSE) | (is.na(b4p) == FALSE),
         hasRelDbh = is.na(a10) == FALSE,
         hasRelHt = is.na(a9) == FALSE,
         hasSignificantBasalArea = hasBasalArea * significant,
         hasSignificantPhysio = hasPhysio * significant,
         hasSignificantRelative = (hasRelDbh | hasRelHt) * significant,
         nBasalAreaSignificant = ((is.na(a2) == FALSE) + (is.na(a3) == FALSE)) * significant, # Sharma-Parton and Sharma-Zhang height forms are not counted here as there is no test for whether BA or BAL are significant predictors (modified Sharma-Parton for DBH does not use basal area)
         nPhysioSignificant = ((is.na(a4) == FALSE) + (is.na(a5) == FALSE) + (is.na(a6) == FALSE) + (is.na(a7) == FALSE) + (is.na(a8) == FALSE)) * significant,
         significantBalOrAat = (is.na(a2) == FALSE) * significant,
         significantBAorAba = (is.na(a3) == FALSE) * significant,
         significantElevation = (is.na(a4) == FALSE) * significant,
         significantSlope = (is.na(a5) == FALSE) * significant,
         significantSinAspect = (is.na(a6) == FALSE) * significant,
         significantCosAspect = (is.na(a7) == FALSE) * significant,
         significantTopographicShelter = (is.na(a8) == FALSE) * significant,
         significantRelHt = (is.na(a9) == FALSE) * significant,
         significantRelDbh = (is.na(a10) == FALSE) * significant,
         a1 = is.na(a1) == FALSE, a1p = is.na(a1p) == FALSE,
         a2 = is.na(a2) == FALSE, a2p = is.na(a2p) == FALSE,
         a3 = is.na(a3) == FALSE, a3p = is.na(a3p) == FALSE,
         a4 = is.na(a4) == FALSE, a5 = is.na(a5) == FALSE, # a5-a8: plantation effects not tested on physiographic predictors
         a6 = is.na(a6) == FALSE, a7 = is.na(a7) == FALSE, a8 = is.na(a8) == FALSE, 
         a9 = is.na(a9) == FALSE, a9p = is.na(a9p) == FALSE,
         a10 = is.na(a10) == FALSE, a10p = is.na(a10p) == FALSE, 
         b1 = is.na(b1) == FALSE, b1p = is.na(b1p) == FALSE,
         b2 = is.na(b2) == FALSE, b2p = is.na(b2p) == FALSE,
         b3 = is.na(b3) == FALSE, b3p = is.na(b3p) == FALSE,
         b4 = is.na(b4) == FALSE, b4p = is.na(b4p) == FALSE)
predictorVariableStats = predictorVariableResults %>% 
  group_by(responseVariable, species) %>%
  summarize(n = n(),
            hasBasalArea = sum(hasBasalArea),
            hasPhysio = sum(hasPhysio),
            hasPlantation = sum(hasPlantation),
            hasRelDbh = sum(hasRelDbh),
            hasRelHt = sum(hasRelHt),
            hasSignificantBasalArea = sum(hasSignificantBasalArea),
            hasSignificantPhysio = sum(hasSignificantPhysio),
            hasSignificantRelative = sum(hasSignificantRelative),
            nBasalAreaSignificant = sum(nBasalAreaSignificant),
            nPhysioSignificant = sum(nPhysioSignificant),
            significantBAorAba = sum(significantBAorAba),
            significantBalOrAat = sum(significantBalOrAat),
            significantElevation = sum(significantElevation),
            significantSlope = sum(significantSlope),
            significantSinAspect = sum(significantSinAspect),
            significantCosAspect = sum(significantCosAspect),
            significantTopographicShelter = sum(significantTopographicShelter),
            significantRelHt = sum(significantRelHt),
            significantRelDbh = sum(significantRelDbh),
            significant = sum(significant),
            a1 = sum(a1), a1p = sum(a1p),
            a2 = sum(a2), a2p = sum(a2p),
            a3 = sum(a3), a3p = sum(a3p),
            a4 = sum(a4), a5 = sum(a5), # a5-a8: plantation effects not tested on physiographic predictors
            a6 = sum(a6), a7 = sum(a7), a8 = sum(a8), 
            a9 = sum(a9), a9p = sum(a9p),
            a10 = sum(a10), a10p = sum(a10p), 
            b1 = sum(b1), b1p = sum(b1p),
            b2 = sum(b2), b2p = sum(b2p),
            b3 = sum(b3), b3p = sum(b3p),
            b4 = sum(b4), b4p = sum(b4p),
            totalC = a1 + a2 + a3 + a4 + a9 + a10 + b1 + b2 + b3 + b4,
            totalCp = a1p + a2p + a3p + a9p + a10p + b1p + b2p + b3p + b4p,
            plantationPct = 100 * hasPlantation / n, 
            significantPct = 100 * significant / n,
            significantRelDbhPct = 100 * significantRelDbh / hasRelDbh,
            significantRelHtPct = 100 * significantRelHt / hasRelHt,
            .groups = "drop")


## summary for Abstract (plus AUC counting for Section 2.3)
# form counts, changes in model efficiency with inclusion of generalizing predictors
# 39 height forms + 38 DBH forms = 77 + 3 control forms => 539 * k * r fits
#
# number of AUCs per k x r cross validation combination (Section 2.3)
#   number of unique pairs in set = 1/2 * n * (n - 1) => 741 height pairs + 703 DBH pairs
#   number of AUCs = pairs * five goodness of fit stats * seven species = (741 + 703) * 5 * 7 = 50540 max, 49905 actual
#   but actual number is reduced
(modelCounts = heightDiameterResults %>% 
  filter(fitSet == "primary", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% 
  group_by(responseVariable, species) %>%
  summarize(nForms = length(unique(name)), nFits = n(), nAucs = 0.5 * nForms * (nForms - 1) * 5))
modelCounts %>% summarize(nAucsIfAllFitsSucceeded = sum(nAucs)) %>% mutate(totalAucsIfAllFitsSucceeded = sum(nAucsIfAllFitsSucceeded))

(actualAucs = heightDiameterModelAucs %>% filter(name != otherModelName) %>%
    group_by(responseVariable, species) %>% 
    mutate(otherBaseName = if_else(word(otherModelName) %in% c("REML", "modified", "unified"), paste(word(otherModelName, 1), word(otherModelName, 2)), word(otherModelName)),
           otherIsBaseForm = (str_detect(otherBaseName, "Sharma-") == FALSE) & (str_detect(otherBaseName, "ABA\\+T") == FALSE) & (str_detect(otherBaseName, "BA\\+L") == FALSE) & (str_detect(otherBaseName, "physio") == FALSE) & (str_detect(otherBaseName, "RelDbh") == FALSE) & (str_detect(otherBaseName, "RelHt") == FALSE),
           isBaseToBaseOrSharmaToSharma = (isBaseForm | (name %in% c("Sharma-Parton", "Sharma-Zhang", "modified Sharma-Parton"))) & (otherIsBaseForm | (name %in% c("Sharma-Parton", "Sharma-Zhang", "modified Sharma-Parton")))) %>%
    summarize(nPairs = n(), nBaseSharma = sum(isBaseToBaseOrSharmaToSharma, na.rm = TRUE)) %>% 
    mutate(nAucs = ceiling(0.5 * nPairs) * 5,
           nBaseSharmaAucs = ceiling(0.5 * nBaseSharma) * 5,
           nTotalAucs = sum(nAucs),
           nTotalBaseSharmaAucs = sum(nBaseSharmaAucs)))

#redAlder = heightDiameterResults %>% filter(responseVariable == "DBH", species == "red alder", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% reframe(name = unique(name))
#hemlock = heightDiameterResults %>% filter(responseVariable == "DBH", species == "western hemlock", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% reframe(name = unique(name))
#maple = heightDiameterResults %>% filter(responseVariable == "DBH", species == "bigleaf maple", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% reframe(name = unique(name))
#myrtle = heightDiameterResults %>% filter(responseVariable == "DBH", species == "Oregon myrtle", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% reframe(name = unique(name))
#redcedar = heightDiameterResults %>% filter(responseVariable == "DBH", species == "western redcedar", if_else(responseVariable == "height", str_detect(name, "RelHt") == FALSE, str_detect(name, "BA\\+L") == FALSE)) %>% reframe(name = unique(name))
#setdiff(redAlder, hemlock) # Schnute inverse
#setdiff(redAlder, maple) # Chapman-Richards replace ABA+T
#setdiff(redAlder, myrtle) # modified Sharma-Parton
#setdiff(redAlder, redcedar) # Schnute inverse

heightDiameterResults %>% 
  filter(fitting %in% c("gsl_nls", "nlrob"), significant) %>% 
  mutate(deltaAicN = aic/nValidation - min(aic/nValidation, na.rm = TRUE)) %>% # change to unrestricted ΔAICn
  filter(baseName %in% c("Chapman-Richards", "Ruark", "Sharma-Parton", "Sharma-Zhang", "Sibbesen", "Weibull")) %>% # remove base forms which weren't generalized
  group_by(responseVariable, species, baseName) %>%
  mutate(nseBase = median(if_else(isBaseForm | (baseName == "Sharma-Parton") | (baseName == "Sharma-Zhang"), nse, NA_real_), na.rm = TRUE),
         nseStandRelDelta = if_else((hasStand | hasRelative) & (hasPhysio == FALSE), nse, NA_real_) - nseBase,
         nsePhysioDelta = if_else(hasPhysio & (hasStand == FALSE) & (hasRelative == FALSE), nse, NA_real_) - nseBase,
         nseCombinedDelta = if_else(hasPhysio & (hasStand | hasRelative), nse, NA_real_) - nseBase) %>%
  #group_by(fitSet, fixedWeight, responseVariable, species) %>%
  #slice_max(nse, n = 3 * 10 * 10) %>% # nmodels * nfolds * nrepetitions
  group_by(fitSet, fixedWeight, responseVariable, baseName) %>%
  summarize(nForms = n_distinct(name),
            nFits = n(),
            mae = mean(mae),
            mape = mean(mape),
            rmse = mean(rmse),
            #deltaAicN = mean(deltaAicN),
            nseMean = mean(nse),
            nseBase = nseBase[1],
            nseStandRelDelta = median(nseStandRelDelta, na.rm = TRUE),
            nsePhysioDelta = median(nsePhysioDelta, na.rm = TRUE),
            nseCombinedDelta = median(nseCombinedDelta, na.rm = TRUE),
            .groups = "drop") %>%
  filter(baseName %in% c("Chapman-Richards", "Sharma-Parton", "Ruark", "Sibbesen"), (responseVariable != "height") | (baseName != "Sibbesen")) %>%
  arrange(desc(fitSet), desc(responseVariable), is.na(fixedWeight) == FALSE) %>%
  select(-fixedWeight)

# highest ranked base forms
(topBaseForms = left_join(left_join(left_join(left_join(preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucMab, n = 1) %>% select(responseVariable, species, mabName),
                                                        preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucMae, n = 1) %>% select(responseVariable, species, maeName),
                                                        by = join_by(responseVariable, species)),
                                              preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucRmse, n = 1) %>% select(responseVariable, species, rmseName),
                                              by = join_by(responseVariable, species)),
                                    preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucDeltaAicN, n = 1) %>% select(responseVariable, species, aicName),
                                    by = join_by(responseVariable, species)),
                          preferredForms %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% slice_max(aucNse, n = 1) %>% select(responseVariable, species, nseName),
                          by = join_by(responseVariable, species)) %>%
    mutate(isGam = (mabName == "REML GAM") + (maeName == "REML GAM") + (rmseName == "REML GAM") + (aicName == "REML GAM") + (nseName == "REML GAM")) %>%
    arrange(desc(responseVariable)))
(pctGamPreference = 100 * sum(topBaseForms$isGam) / (5 * nrow(topBaseForms))) # GAM preference percentage for abstract


## summaries for Results
# fraction of nonlinear base forms which are more accurate than a base GAM
improvementThresholdProbability = 0.5
heightDiameterModelAucs %>% filter(name != "REML GAM", otherModelName == "REML GAM", isBaseForm) %>%
  group_by(responseVariable, species) %>% # restriction to primary fit set and base forms is above
  summarize(deltaAic = 100 * sum(aucDeltaAicN > improvementThresholdProbability) / n(),
            mab = 100 * sum(aucMab > improvementThresholdProbability) / n(),
            mae = 100 * sum(aucMae > improvementThresholdProbability) / n(),
            nse = 100 * sum(aucNse > improvementThresholdProbability) / n(),
            rmse = 100 * sum(aucRmse > improvementThresholdProbability) / n(),
            aicMostPreferred = name[which.max(aucDeltaAicN)],
            .groups = "drop") %>%
  arrange(desc(responseVariable))

primaryResults %>% group_by(responseVariable, species) %>% # plantation effect quantiles
  reframe(quantiles = c(0.05, 0.50, 0.95), plantFx = quantile(meanAbsolutePlantationEffect, probs = quantiles, na.rm = TRUE), plantFxPct = quantile(meanAbsolutePercentPlantationEffect, probs = quantiles, na.rm = TRUE)) %>%
  pivot_wider(id_cols = c("responseVariable", "species"), names_from = "quantiles", values_from = c("plantFx", "plantFxPct")) %>%
  arrange(desc(responseVariable))

# AIC selection of generalizing predictors based on ΔAICn AUC
improvementThresholdProbability = 0.5
heightDiameterModelAucs %>% 
  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name))) %>%
  filter(isBaseForm == FALSE, fitting %in% c("gsl_nls", "nlrob"), significant, baseName %in% c("Chapman-Richards", "Sharma-Parton", "Ruark", "Sibbesen"), (responseVariable == "height") | (baseName != "Sibbesen")) %>%
  group_by(responseVariable) %>%
  summarize(n = sum(is.na(aucDeltaAicN) == FALSE), 
            baseBetter = sum(aucDeltaAicN <= improvementThresholdProbability), # base form is ΔAICn preferable
            generalizedBetter = sum(aucDeltaAicN > improvementThresholdProbability), # a generalized form is preferable
            physioAicMedianAuc = median(if_else(hasPhysio & (hasStand == FALSE) & (hasRelative == FALSE), aucDeltaAicN, NA_real_), na.rm = TRUE), 
            standRelAicMedianAuc = median(if_else((hasPhysio == FALSE) & (hasStand | hasRelative), NA_real_, aucDeltaAicN), na.rm = TRUE),
            bothAicMedianAuc = median(if_else(hasPhysio & (hasStand | hasRelative), NA_real_, aucDeltaAicN), na.rm = TRUE)) %>%
  mutate(baseBetterPct = 100 * baseBetter / n, 
         genBetterPct = 100 * generalizedBetter / n) %>%
  arrange(desc(responseVariable))

# coefficient counts
#heightDiameterCoefficientCounts = heightDiameterCoefficients %>% rowwise() %>% mutate(coefficients = ncol(heightDiameterCoefficients) - sum(is.na(c_across(is.numeric)))) %>%
#  select(responseVariable, species, name, coefficients, coefficients) # slow
#print(heightDiameterCoefficientCounts %>% filter(fitting == "gam"), n = 50)

# frequency of generalizing predictor selection
predictorVariableStats %>%
  mutate(balAat = significantBalOrAat / hasSignificantBasalArea, baAba = significantBAorAba / hasSignificantBasalArea,
         elev = significantElevation / hasSignificantPhysio, slope = significantSlope / hasSignificantPhysio, sinAspect = significantSinAspect / hasSignificantPhysio, cosAspect = significantCosAspect / hasSignificantPhysio, tsi = significantTopographicShelter / hasSignificantPhysio) %>%
  select(responseVariable, species, balAat, baAba, elev, slope, sinAspect, cosAspect, tsi) %>%
  arrange(desc(responseVariable)) # since relative DBH and height are singleton groups they are always selected with probability 1 if when significant

# model efficiency of generalized predictors
#generalizedPredictors = primaryResults %>%
#  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name))) %>%
#  filter(baseName %in% c("Chapman-Richards", "Ruark", "Sharma-Parton", "Sibbesen"), (responseVariable != "height")  | (baseName != "Sibbesen"), significant) # remove base forms which weren't generalized and non-significant forms
#generalizedPredictorMapbModel = lm(mapb ~ responseVariable + species + baseName + hasPhysio + responseVariable:hasRelative + hasStand, generalizedPredictors)
#summary(generalizedPredictorMapbModel)
#generalizedPredictorMapeModel = lm(mape ~ responseVariable + species + baseName + hasPhysio + responseVariable:hasRelative + hasStand, generalizedPredictors)
#summary(generalizedPredictorMapeModel)
#generalizedPredictorRmspeModel = lm(rmspe ~ responseVariable + species + baseName + hasPhysio + responseVariable:hasRelative + hasStand, generalizedPredictors)
#summary(generalizedPredictorRmspeModel)
#generalizedPredictorAicNModel = lm(deltaAicN ~ responseVariable + species + baseName + hasPhysio + responseVariable:hasRelative + hasStand, generalizedPredictors)
#summary(generalizedPredictorAicNModel)
#generalizedPredictorNseModel = lm(nse ~ responseVariable + species + baseName + hasPhysio + responseVariable:hasRelative + hasStand, generalizedPredictors)
#summary(generalizedPredictorNseModel)

# AUCs of generalized predictors
#generalizedAucs = heightDiameterModelAucs %>%
#  filter(fitting %in% c("gsl_nls", "nlrob"), (responseVariable != "height") | (str_detect(name, "RelHt") == FALSE), isBaseForm == FALSE, significant, otherModelSignificant) %>% 
#  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name))) %>%
#  filter(baseName %in% c("Chapman-Richards", "Ruark", "Sharma-Parton", "Sibbesen"), (responseVariable != "height")  | (baseName != "Sibbesen"), # remove base forms which weren't generalized
#         str_starts(otherModelName, baseName))
#generalizedAucModel = lm(aucNse ~ hasPhysio + hasRelative + hasStand, generalizedAucs) # response variable, species, and base form not significant: normalized out
#summary(generalizedAucModel)
#confint(generalizedAucModel, level = 0.99)

# highest ranked model forms by goodness of fit metric
print(left_join(left_join(left_join(left_join(preferredForms %>% group_by(responseVariable, species) %>% slice_max(aucMab, n = 1) %>% select(mabName),
                                              preferredForms %>% group_by(responseVariable, species) %>% slice_max(aucMae, n = 1) %>% select(maeName),
                                              by = join_by(responseVariable, species)),
                                    preferredForms %>% group_by(responseVariable, species) %>% slice_max(aucRmse, n = 1) %>% select(rmseName),
                                    by = join_by(responseVariable, species)),
                          preferredForms %>% group_by(responseVariable, species) %>% slice_max(aucDeltaAicN, n = 1) %>% select(aicName),
                          by = join_by(responseVariable, species)),
                preferredForms %>% group_by(responseVariable, species) %>% slice_max(aucNse, n = 1) %>% select(nseName),
                by = join_by(responseVariable, species)) %>%
        pivot_longer(cols = ends_with("Name"), names_to = "statistic", values_to = "name"), # %>%
      #group_by(responseVariable, species) %>%
      #summarize(modelsSelected = length(unique(name)))
      n = 70)

# summarize fitting success
# Show results for all fit sets here, whether or not primary.
heightDiameterResults %>% group_by(fitSet, responseVariable, species, name) %>%
  # reduce cross validated fits to a single record
  summarize(nlrob = fitting[1] == "nlrob", 
            gslNls = fitting[1] == "gsl_nls", 
            gnls = (fitting[1] == "gnls") | (is.na(fitting[1]) & str_detect(name[1], "GNLS")),
            lm = fitting[1] == "lm", 
            gam = fitting[1] == "gam",
            fail = is.na(nse[1]),
            .groups = "drop") %>%
  # summarize across fitting attempts at equal weight, regardless of cross validation settings
  group_by(fitSet) %>%
  summarize(n = n(),
            nlrob = sum(nlrob, na.rm = TRUE), 
            gslNls = sum(gslNls, na.rm = TRUE), 
            gnls = sum(gnls, na.rm = TRUE), 
            lm = sum(lm, na.rm = TRUE), 
            gam = sum(gam, na.rm = TRUE),
            fail = sum(fail),
            uniqueNonlinear = n() - gnls - lm - gam) %>%
  mutate(nlrobPct = 100 * nlrob / uniqueNonlinear, 
         gslNlsPct = nlrobPct + 100 * gslNls / uniqueNonlinear, 
         totalFailPct = 100 * fail / uniqueNonlinear) %>%
  arrange(desc(fitSet))
# fittings which failed
heightDiameterResults %>% filter(is.na(nse)) %>% select(fitSet, responseVariable, species, name) %>% arrange(desc(fitSet))

# generalization's effects on MAB versus plantation effect sizes
deltaMapb = heightDiameterResults %>% filter(fitSet == "primary", significant) %>%
  group_by(responseVariable, species, baseName) %>%
  mutate(deltaMapb = mapb - median(na.omit(if_else(isBaseForm | (name == "Sharma-Parton"), mapb, NA_real_)))) %>% # find MAB relative to median MAB (%) of all base fits
  filter(isBaseForm == FALSE, baseName %in% c("REML GAM", "Chapman-Richards", "Ruark", "Sharma-Parton", "Sibbesen"))
print(deltaMapb %>% group_by(responseVariable, species, baseName) %>%
  reframe(quantiles = c(0, 0.25, 0.5, 0.75, 1), deltaMapb = quantile(deltaMapb, probs = quantiles)) %>%
  pivot_wider(id_cols = c("responseVariable", "species", "baseName"), names_from = "quantiles", names_prefix = "generalizationDeltaMapbQ", values_from = "deltaMapb"), n = 65)

ggplot() +
  geom_violin(aes(x = deltaMapb, y = species, color = responseVariable), deltaMapb %>% filter(responseVariable == "height"), draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "generalization ΔMAB, %", y = NULL, color = NULL, title = "a) height") +
ggplot() +
  geom_violin(aes(x = deltaMapb, y = species, color = responseVariable), deltaMapb %>% filter(responseVariable == "DBH"), draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "generalization ΔMAB, %", y = NULL, color = NULL, title = "b) DBH") +
  theme(axis.text.y = element_blank()) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  coord_cartesian(xlim = c(-20, 100)) &
  scale_color_manual(breaks = c("height", "DBH"), values = c("green4", "burlywood4")) &
  scale_y_discrete(limits = rev)

ggplot() +
  geom_point(aes(x = deltaMapb, y = species, color = baseName), deltaMapb %>% filter(responseVariable == "height") %>% group_by(species, baseName) %>% summarize(deltaMapb = mean(abs(deltaMapb)), .groups = "drop"), alpha = 0.6, size = 2) +
  labs(x = "mean absolute ΔMAB, %", y = NULL, color = "height", title = "a) height") +
ggplot() +
  geom_point(aes(x = deltaMapb, y = species, color = baseName), deltaMapb %>% filter(responseVariable == "DBH") %>% group_by(species, baseName) %>% summarize(deltaMapb = mean(abs(deltaMapb)), .groups = "drop"), alpha = 0.6, size = 2) + # Oregon myrtle Chapman-Richards max @ 11.4%
  labs(x = "mean absolute ΔMAB, %", y = NULL, color = "DBH", title = "b) DBH") +
  theme(axis.text.y = element_blank()) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  coord_cartesian(xlim = c(0, 12)) &
  guides(color = guide_legend(override.aes = list(alpha = 0.8))) &
  scale_x_continuous(breaks = seq(0, 12, by = 3)) &
  scale_y_discrete(limits = rev)


## Figure 1: overall dataset summary
plot_exploratory(trees2016 %>% filter(isLiveUnbroken, isConifer), speciesLabel = "conifer", maxTreesMeasured = 170, titleLetters = c(plotLetters[1], "", ""), omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(trees2016 %>% filter(isLiveUnbroken, isConifer == FALSE), speciesLabel = "broadleaf", maxTreesMeasured = 170, titleLetters = c(plotLetters[2], "", "")) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.png", height = 13, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.tif", height = 13, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.pdf", height = 13, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 2: height-diameter AUCs
heightFromDiameterModelComparison = heightDiameterModelRanking %>% filter(responseVariable == "height") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "height"))$name)))
plot_auc_bank(heightFromDiameterModelComparison, legendHjustification = 1.05, plotRightMargin = 10)
#ggsave("trees/height-diameter/figures/Figure 02 height accuracy AUC median.png", height = 17, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 02 height accuracy AUC median.tif", height = 20, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 02 height accuracy AUC median.pdf", height = 20, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 3: diameter AUCs
diameterFromHeightModelComparison = heightDiameterModelRanking %>% filter(responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "DBH"))$name)))
plot_auc_bank(diameterFromHeightModelComparison, legendHjustification = 1.02)
#ggsave("trees/height-diameter/figures/Figure 03 diameter accuracy AUC median.png", height = 17.5, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 03 diameter accuracy AUC median.tif", height = 20.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 03 diameter accuracy AUC median.pdf", height = 20.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 4: model efficiency
# Long tail of negative model efficiencies is suppressed as histogram plotting becomes very slow even if bins
# are not within the display window set by coord_cartesian().
smallTreeEfficiency = left_join(trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE) %>% mutate(species = factor(speciesGroup, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other"))) %>%
                                  group_by(species) %>%
                                  summarize(stems = sum(TreeCount), smallStems = sum(TreeCount * (DBH <= 20.0)), pctSmall = 100 * smallStems / stems),
                                primaryResults,
                                by = c("species"))

ggplot(primaryResults %>% filter(responseVariable == "height", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = NULL, y = "fraction of fits, %", fill = NULL, title = bquote(.(plotLetters[1])~"base form height prediction")) +
  #labs(x = NULL, y = "fraction of fits, %", fill = NULL, title = bquote(bold(.(plotLetters[1]))~"base form height prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(responseVariable == "DBH", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[2])~"base form DBH prediction")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[3]))~"base form DBH prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#ggplot(smallTreeEfficiency) + 
#  #geom_violin(aes(x = pctSmall, y = nse, color = species, group = species), draw_quantiles = c(0.25, 0.50, 0.75), fill = alpha("white", 0.5), position = position_identity(), width = 2) + # too stretched to be useful due to tail of negative model efficiencies
#  geom_boxplot(aes(x = pctSmall, y = nse, color = species, group = species), fill = alpha("white", 0.5), position = position_identity(), na.rm = TRUE, outlier.alpha = 0.1, outlier.size = 0.1, width = 4) +
#  coord_equal(ratio = 60, xlim = c(0, 100), ylim = c(0, 21.5/20)) +
#  labs(x = "stems less than 20 cm DBH, %", y = "model efficiency", color = NULL, title = "c) lack of small stem effect") +
#  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
#  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
#  theme(plot.title = element_text(vjust = -27)) +
ggplot(primaryResults %>% filter(responseVariable == "height", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = "model efficiency", y = "fraction of fits, %", fill = NULL, title = bquote(.(plotLetters[3])~"generalized height prediction")) +
  #labs(x = "model efficiency", y = "fraction of fits, %", fill = NULL, title = bquote(bold(.(plotLetters[3]))~"generalized height prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(responseVariable == "DBH", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 24)) +
  labs(x = "model efficiency", y = NULL, fill = NULL, title = bquote(.(plotLetters[4])~"generalized DBH prediction")) +
  #labs(x = "model efficiency", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[4]))~"generalized DBH prediction")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
#plot_layout(nrow = 3, ncol = 3, design = "123\n453", guides = "collect") &
plot_layout(nrow = 2, ncol = 2, guides = "collect") &
  guides(color = "none", fill = guide_legend(ncol = 7)) &
  scale_color_manual(breaks = levels(primaryResults$species), limits = levels(primaryResults$species), values = speciesGroupColors) &
  scale_fill_manual(breaks = levels(primaryResults$species), limits = levels(primaryResults$species), values = speciesGroupColors) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.6, "line"), legend.position = "bottom", legend.text = element_text(margin = margin(l = 1, r = 6.5)))
#ggsave("trees/height-diameter/figures/Figure 04 model efficiency.png", height = 10, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 04 model efficiency.tif", height = 10, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 04 model efficiency.pdf", height = 10, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 5: fitting success, model significance, and predictor variable use
modelFitStats = primaryResults %>% group_by(responseVariable, species, name) %>%
  summarize(modelFit = is.na(nse[1]) == FALSE, .groups = "drop_last") %>% 
  summarize(fitPct = 100 * mean(modelFit), .groups = "drop") %>%
  mutate(species = factor(species, levels = levels(predictorVariableResults$species)))
  
ggplot(modelFitStats %>% filter(responseVariable == "height")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[1])~"height fits")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[1]))~"height fits")) +
  scale_y_discrete(limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableStats %>% filter(responseVariable == "height")) +
  geom_col(aes(x = significantPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[2])~"height forms")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[2]))~"height forms")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(primaryResults %>% filter(responseVariable == "height", significant)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 125)) +
  labs(x = NULL, y = NULL, color = NULL, title = bquote(.(plotLetters[3])~"height")) +
  #labs(x = NULL, y = NULL, color = NULL, title = bquote(bold(.(plotLetters[3]))~"height")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[4])~"BA+L")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[4]))~"BA+L")) +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasPhysio)) +
  geom_count(aes(x = nPhysioSignificant, y = species, color = species), shape = 18) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[5])~"height physiographic")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[5]))~"height physiographic")) +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "height", hasRelDbh)) +
  geom_count(aes(x = significantRelDbh, y = species, color = species)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = bquote(.(plotLetters[6])~"RelDbh")) +
  #labs(x = NULL, y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[6]))~"RelDbh")) +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = seq(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(modelFitStats %>% filter(responseVariable == "DBH")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "model forms\nfit, %", y = NULL, fill = NULL, title = bquote(.(plotLetters[7])~"DBH fits")) +
  #labs(x = "model forms\nfit, %", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[7]))~"DBH fits")) +
  scale_y_discrete(limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableStats %>% filter(responseVariable == "DBH")) +
  geom_col(aes(x = significantPct, y = species, fill = species), width = 0.6) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "significant\nforms, %", y = NULL, fill = NULL, title = bquote(.(plotLetters[8])~"DBH forms")) +
  #labs(x = "significant\nforms, %", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[8]))~"DBH forms")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(primaryResults %>% filter(responseVariable == "DBH", significant)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 125)) +
  labs(x = "mean absolute\nplantation effect, %", y = NULL, color = NULL, title = bquote(.(plotLetters[9])~"DBH")) +
  #labs(x = "mean absolute\nplantation effect, %", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[9]))~"DBH")) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[10])~"ABA+T")) +
  #labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[10]))~"ABA+T")) +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasPhysio)) +
  geom_count(aes(x = nPhysioSignificant, y = species, color = species), shape = 18) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[11])~"DBH physiographic")) +
  #labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[11]))~"DBH physiographic")) +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
ggplot(predictorVariableResults %>% filter(responseVariable == "DBH", hasRelHt)) +
  geom_count(aes(x = significantRelHt, y = species, color = species)) + 
  labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(.(plotLetters[12])~"RelHt")) +
  #labs(x = "significant\nvariables", y = NULL, fill = NULL, title = bquote(bold(.(plotLetters[12]))~"RelHt")) +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = c(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(levels(predictorVariableResults$species))) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 6, widths = c(1.7, 1.7, 2.2, 1.7, 2.7, 1.2), guides = "collect") &
  scale_color_manual(breaks = levels(predictorVariableStats$species), limits = levels(predictorVariableResults$species), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) &
  scale_fill_manual(breaks = levels(predictorVariableStats$species), limits = levels(predictorVariableResults), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) &
  scale_size_area(max_size = 4.5) & # 4.5 for Cairo .pdf, 5.0 for .tif
  theme(legend.position = "none")
#ggsave("trees/height-diameter/figures/Figure 05 predictor variables.png", height = 10, width = 20, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/Figure 05 predictor variables.tif", height = 10, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 05 predictor variables.pdf", height = 10, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 6: Douglas-fir and red alder preferred models
# preferred forms from 10x10 cross validation
#print(preferredForms %>% filter(species %in% c("Douglas-fir", "red alder")) %>% select(-mabName, -aucMab, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other"), levels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 33)
print(heightDiameterModelRanking %>% filter(significant, isBaseForm) %>% select(-aucBlended, -starts_with("has"), -speciesFraction) %>%
        pivot_longer(starts_with("auc"), names_to = "statistic", values_to = "auc") %>% mutate(statistic = factor(statistic, levels = c("aucMab", "aucMae", "aucRmse", "aucDeltaAicN", "aucNse"))) %>%
        group_by(species, responseVariable, statistic) %>%
        slice_max(auc, n = 1) %>% arrange(species, desc(responseVariable), statistic), n = 70)

if (exists("psmeHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/PSME preferred models.Rdata") }
if (exists("alruHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/ALRU2 preferred models.Rdata") }

# Temesgen et al. 2007 height = 1.3 + exp(b1 - b2 * DBH^b3) => b1 - b2 * DBH^b3 = ln(height - 1.3) => DBH^b3 = 1/b2 * (b1 - ln(height - 1.3))
#                      DBH = (1/b2 * (b1 - ln(height - 1.3)))^(1/b3)
psmeReference = bind_rows(bind_rows(psmeHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    psmeHeightFromDiameterPreferred$ratkowsky$stats %>% mutate(model = "base form 2"),
                                    psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Temesgen et al. 2007", "height", psme2016, 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), 4, psme2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(psmeDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    psmeDiameterFromHeightPreferred$sibbesenReplace$stats %>% mutate(model = "base form 2"),
                                    psmeDiameterFromHeightPreferred$gamAbatPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Temesgen et al. 2007", "DBH", psme2016, (1/6.7792 * (5.7567 - log(psme2016$TotalHt - 1.3)))^(1/-0.2795), 4, psme2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "PSME", adaptiveWeightFraction = 0)
psmeReferenceDbh = seq(0, 241) # Temesgen et al. 2007 dataset limit of 190 cm, extended of necessity

# Hanus et al. height = h0 + a1 * exp(b1 - b2 * (b3 * DBH)^b4) => b1 - b2 * (b3 * DBH)^b4 = ln(1/a1 * (height - h0)) => (b3 * DBH)^b4 = 1/b2 * (b1 - ln(1/a1 * (height - h0)))
#              DBH = 1/b3 * (1/b2 * (b1 - ln(1/a1 * (height - h0)))^(1/b4))
alruReference = bind_rows(bind_rows(alruHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    alruHeightFromDiameterPreferred$hossfeld$stats %>% mutate(model = "base form 2"),
                                    alruHeightFromDiameterPreferred$sharmaPartonPhysio$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Hanus et al. 1999", "height", alru2016, 1.37 + 0.3048 * exp(4.41820972 - 12.00274935 * (0.393701 * alru2016$DBH)^-2.13835482), 4, alru2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(alruDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    alruDiameterFromHeightPreferred$chapmanRichards$stats %>% mutate(model = "base form 2"),
                                    alruDiameterFromHeightPreferred$ruarkAbatPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Hanus et al. 1999", "DBH", alru2016, 1/0.393701 * (1/12.00274935 * (4.41820972 - log(1/0.3048 * (alru2016$TotalHt - 1.37))))^(1/-2.13835482), 4, alru2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "ALRU2", adaptiveWeightFraction = 0)
alruPreviousDbhFailures = sum(alru2016$TreeCount * is.na(1/0.393701 * (1/12.00274935 * (4.41820972 - log(1/0.3048 * (alru2016$TotalHt - 1.37))))^(1/-2.13835482))) # 484 of 3307 trees = 14.6% numerical failure due to negative power of negative number
alruReferenceDbh = seq(0, 135) # Hanus et al. 1999 dataset limit of 50 cm, extended of necessity plus for plot visiblity

ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh), color = "Sharma-Parton BA+L RelDbh physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$gam), color = "REML GAM (MAB, MAE, RMSE, ME)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$ratkowsky), color = "Ratkowsky (ΔAICn)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) + # or Prodan
  geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(.(plotLetters[1])~"Douglas-fir height")) +
  #labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[1]))~"Douglas-fir height")) +
  scale_color_manual(breaks = c("REML GAM (MAB, MAE, RMSE, ME)", "Ratkowsky (ΔAICn)", "Sharma-Parton BA+L RelDbh physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey60")) +
  theme(legend.key = element_rect(fill = alpha("white", 0.5)), legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gamAbatPhysio), y = psme2016physio$TotalHt, color = "REML GAM ABA+T physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gam), y = psme2016$TotalHt, color = "REML GAM (MAB, MAE, RMSE, ME)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$sibbesenReplace), y = psme2016$TotalHt, color = "Sibbesen replace (ΔAICn)", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(.(plotLetters[2])~"Douglas-fir DBH")) +
  #labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[2]))~"Douglas-fir DBH")) +
  scale_color_manual(breaks = c("REML GAM (MAB, MAE, RMSE, ME)", "Sibbesen replace (ΔAICn)", "REML GAM ABA+T physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = alru2016physio$DBH, y = predict(alruHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$gam), color = "REML GAM (MAE, RMSE, ME)", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$hossfeld), color = "Weibull (ΔAICn)", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  #geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = alruReferenceDbh, y = 1.37 + 0.3048 * exp(4.41820972 - 12.00274935 * (0.393701 * alruReferenceDbh)^-2.13835482), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL, title = bquote(.(plotLetters[3])~"red alder height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[3]))~"red alder height")) +
  scale_color_manual(breaks = c("REML GAM (MAE, RMSE, ME)", "Weibull (ΔAICn)", "Sharma-Parton physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  #geom_line(aes(x = predict(alruDiameterFromHeightPreferred$gamAbatPhysioRelHt), y = alru2016physio$TotalHt, color = "REML GAM ABA+T RelHt physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$ruarkAbatPhysio), y = alru2016physio$TotalHt, color = "Ruark ABA+T physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$chapmanReplace), y = alru2016$TotalHt, color = "Chapman-Richards replace (RMSE, ME)", group = alru2016$isPlantation, linetype = alru2016$isPlantation), orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$gam), y = alru2016$TotalHt, color = "REML GAM (MAB, MAE, ΔAICn)", group = alru2016$isPlantation, linetype = alru2016$isPlantation), orientation = "y") +
  #geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = alruReferenceDbh, y = 1.37 + 0.3048 * exp(4.41820972 - 12.00274935 * (0.393701 * alruReferenceDbh)^-2.13835482), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL, title = bquote(.(plotLetters[4])~"red alder DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[4]))~"red alder DBH")) +
  scale_color_manual(breaks = c("REML GAM (MAB, MAE, ΔAICn)", "Chapman-Richards replace (RMSE, ME)", "Ruark ABA+T physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 06 PSME-ALRU2 curves.png", height = 12, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 06 PSME-ALRU2 curves.tif", height = 12, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 06 PSME-ALRU2 curves 1000.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 7: western hemlock and bigleaf maple preferred models
#print(preferredForms %>% filter(species %in% c("western hemlock", "bigleaf maple")) %>% select(-mabName, -aucMab, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other"), levels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 32)
if (exists("tsheHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/TSHE preferred models.Rdata") }
if (exists("acmaHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/ACMA3 preferred models.Rdata") }

tsheReference = bind_rows(bind_rows(tsheHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    tsheHeightFromDiameterPreferred$chapmanRichards$stats %>% mutate(model = "base form 2"),
                                    tsheHeightFromDiameterPreferred$chapmanRichardsBalPhysio$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Hanus et al. 1999", "height", tshe2016, 1.37 + 0.3048 * exp(6.555344622 - 5.137174162 * (0.393701 * tshe2016$DBH)^-0.3645508), 4, tshe2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(tsheDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    tsheDiameterFromHeightPreferred$linear$stats %>% mutate(model = "base form 2"),
                                    tsheDiameterFromHeightPreferred$gamRelHt$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Hanus et al. 1999", "DBH", tshe2016, 1/0.393701 * (1/5.137174162 * (6.555344622 - log(1/0.3048 * (tshe2016$TotalHt - 1.37))))^(1/-0.3645508), 4, tshe2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "TSHE", adaptiveWeightFraction = 0)
tsheReferenceDbh = seq(0, 185) # Hanus et al. 1999 dataset limit of 78 cm, extended of necessity plus for plot visibility

acmaReference = bind_rows(bind_rows(acmaHeightFromDiameterPreferred$ratkowsky$stats %>% mutate(model = "base form 1"),
                                    acmaHeightFromDiameterPreferred$michaelisMenten$stats %>% mutate(model = "base form 2"),
                                    acmaHeightFromDiameterPreferred$weibullBal$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Wang and Hann 1988", "height", acma2016, 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acma2016$DBH)^-0.354756), 4, acma2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(acmaDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    acmaDiameterFromHeightPreferred$ruark$stats %>% mutate(model = "base form 2"),
                                    acmaDiameterFromHeightPreferred$ruarkAbatPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Wang and Hann 1988", "DBH", acma2016, 1/0.393701 * (1/2.70252 * (5.21462 - log(1/1.37 * (acma2016$TotalHt - 1.37))))^(1/-0.354756), 4, acma2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "ACMA3", adaptiveWeightFraction = 0)
acmaReferenceDbh = seq(0, 182) # Wang and Hann 1988 dataset limit of 130 cm, extended of necessity plus for plot visibility

ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = tshe2016physio$DBH, y = predict(tsheHeightFromDiameterPreferred$chapmanRichardsBalPhysio), color = "Chapman-Richards BA+L physio", group = tshe2016physio$isPlantation, linetype = tshe2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterPreferred$gam), color = "REML GAM (all five statistics)", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterPreferred$chapmanRichards), color = "Hossfeld IV", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes (x = tsheReferenceDbh, y = 1.37 + 0.3048 * exp(6.555344622 - 5.137174162 * (0.393701 * tsheReferenceDbh)^-0.3645508), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(.(plotLetters[1])~"western hemlock height")) +
  #labs(x = NULL, y = "height, m", color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[1]))~"western hemlock height")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "Hossfeld IV", "Chapman-Richards BA+L physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$gamRelHt), y = tshe2016$TotalHt, color = "REML GAM RelHt", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$parabolic), y = tshe2016$TotalHt, color = "parabolic", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$gam), y = tshe2016$TotalHt, color = "REML GAM (all five statistics)", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = tsheReferenceDbh, y = 1.37 + 0.3048 * exp(6.555344622 - 5.137174162 * (0.393701 * tsheReferenceDbh)^-0.3645508), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(.(plotLetters[2])~"western hemlock DBH")) +
  #labs(x = NULL, y = NULL, color = NULL, linetype = NULL, title = bquote(bold(.(plotLetters[2]))~"western hemlock DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "parabolic", "REML GAM RelHt", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$weibullBal), color = "Weibull BA+L", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$michaelisMenten), color = "Michaelis-Menten (ΔAICn)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$ratkowsky), color = "Ratkowsky (RMSE, ME)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = acmaReferenceDbh, y = 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acmaReferenceDbh)^-0.354756), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[3])~"bigleaf maple height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[3]))~"bigleaf maple height")) +
  scale_color_manual(breaks = c("Ratkowsky (RMSE, ME)", "Michaelis-Menten (ΔAICn)", "Weibull BA+L", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  #geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$gamAbatPhysioRelHt), y = acma2016$TotalHt, color = "REML GAM ABA+T RelHt physio", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$ruarkAbatPhysio), y = acma2016$TotalHt, color = "Ruark ABA+T physio", group = acma2016$isPlantation, linetype = acma2016$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$gam), y = acma2016$TotalHt, color = "REML GAM (all five statistics)", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$ruark), y = acma2016$TotalHt, color = "Ruark", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = acmaReferenceDbh, y = 1.37 + 0.3048 * exp(5.21462 - 2.70252 * (0.393701 * acmaReferenceDbh)^-0.354756), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[4])~"bigleaf maple DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[4]))~"bigleaf maple DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "Ruark", "Ruark ABA+T physio", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 07 TSHE-ACMA3 curves.png", height = 12, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 07 TSHE-ACMA3 curves.tif", height = 12, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 07 TSHE-ACMA3 curves.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 8: Oregon myrtle and western redcedar preferred models
#print(preferredForms %>% filter(species %in% c("Oregon myrtle", "western redcedar")) %>% select(-mabName, -aucMab, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other"), levels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 32)
if (exists("umcaHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/UMCA preferred models.Rdata") }
if (exists("thplHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/THPL preferred models.Rdata") }

# no prior height model for Oregon myrtle, so no UMCA reference or reference DBH
umcaReference = bind_rows(bind_rows(umcaHeightFromDiameterPreferred$michaelisMenten$stats %>% mutate(model = "base form 1"),
                                    umcaHeightFromDiameterPreferred$linear$stats %>% mutate(model = "base form 2"),
                                    umcaHeightFromDiameterPreferred$sharmaPartonPhysio$stats %>% mutate(model = "generalized height")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(umcaDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    umcaDiameterFromHeightPreferred$sibbesenReplace$stats %>% mutate(model = "base form 2"),
                                    umcaDiameterFromHeightPreferred$gamRelHtPhysio$stats %>% mutate(model = "generalized DBH")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "UMCA", adaptiveWeightFraction = 0)

#print(bind_rows(thplHeightFromDiameterPreferred$gam$stats,
#                get_prediction_stats("height test", "height", thpl2016, predict(thplHeightFromDiameterPreferred$gam), thplHeightFromDiameterPreferred$gam$stats$effectiveDegreesOfFreedom, thpl2016)) %>%
#  mutate(), width = Inf)
#print(bind_rows(thplDiameterFromHeightPreferred$chapmanReplaceRelHt$stats,
#                get_prediction_stats("DBH test", "DBH", thpl2016, predict(thplDiameterFromHeightPreferred$chapmanReplaceRelHt), thplDiameterFromHeightPreferred$chapmanReplaceRelHt$stats$effectiveDegreesOfFreedom, thpl2016)) %>%
#  mutate(), width = Inf)

thplReference = bind_rows(bind_rows(thplHeightFromDiameterPreferred$hossfeld$stats %>% mutate(model = "base form 1"),
                                    thplHeightFromDiameterPreferred$chapmanRichards$stats %>% mutate(model = "base form 2"),
                                    thplHeightFromDiameterPreferred$sharmaPartonPhysio$stats %>% mutate(model = "generalized height"),
                                    get_prediction_stats("Hanus et al. 1999", "height", thpl2016, 1.37 + 0.3048 * exp(7.232880669 - 5.746899904 * (0.393701 * thpl2016$DBH)^-0.271564741), 4, thpl2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "height"),
                          bind_rows(thplDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                    thplDiameterFromHeightPreferred$parabolic$stats %>% mutate(model = "base form 2"),
                                    thplDiameterFromHeightPreferred$gamPhysio$stats %>% mutate(model = "generalized DBH"),
                                    get_prediction_stats("Hanus et al. 1999", "DBH", thpl2016, 1/0.393701 * (1/5.746899904 * (7.232880669 - log(1/0.3048 * (thpl2016$TotalHt - 1.37))))^(1/-0.271564741), 4, thpl2016) %>% mutate(model = "previous model")) %>% 
                            mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "THPL", adaptiveWeightFraction = 0)
thplReferenceDbh = seq(0, 175) # Hanus et al. 1999 dataset limit of 50 cm, extended of necessity plus for plot visibility


ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = umca2016$isPlantation, linetype = umca2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameterPreferred$michaelisMenten), color = "Michaelis-Menten (MAE, RMSE, ME)", group = umca2016$isPlantation, linetype = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = if_else(predict(umcaHeightFromDiameterPreferred$linear) <= 59, predict(umcaHeightFromDiameterPreferred$linear), NA_real_), color = "linear (ΔAICn)", group = umca2016$isPlantation, linetype = umca2016$isPlantation), na.rm = TRUE) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, title = bquote(.(plotLetters[1])~"Oregon myrtle height")) +
  #labs(x = NULL, y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[1]))~"Oregon myrtle height")) +
  scale_color_manual(breaks = c("Michaelis-Menten (MAE, RMSE, ME)", "linear (ΔAICn)", "Sharma-Parton physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$gamRelHt), y = umca2016$TotalHt, color = "REML GAM RelHt", group = umca2016$isPlantation, linetype = umca2016$isPlantation), alpha = 0.5, orientation = "y") +
  #geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$gamRelHtPhysio), y = umca2016$TotalHt, color = "REML GAM RelHt physio", group = umca2016$isPlantation, linetype = umca2016$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$gam), y = umca2016$TotalHt, color = "REML GAM (all five statistics)", group = umca2016$isPlantation, linetype = umca2016$isPlantation)) +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$linear), y = umca2016$TotalHt, color = "linear", group = umca2016$isPlantation, linetype = umca2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, title = bquote(.(plotLetters[2])~"Oregon myrtle DBH")) +
  #labs(x = NULL, y = NULL, color = NULL, title = bquote(bold(.(plotLetters[2]))~"Oregon myrtle DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "linear", "REML GAM RelHt", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position.inside = c(1, 1)) +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameterPreferred$hossfeld), color = "Hossfeld IV (MAE, RMSE, ME)", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation), alpha = 0.8) + # Hossfeld and Michaelis-Menten fits are visually indistinguishable
  geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameterPreferred$chapmanRichards), color = "Chapman-Richards (ΔAICn)", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation), alpha = 0.8) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = thplReferenceDbh, y = 1.37 + 0.3048 * exp(7.232880669 - 5.746899904 * (0.393701 * thplReferenceDbh)^-0.271564741), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[3])~"western redcedar height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[3]))~"western redcedar height")) +
  scale_color_manual(breaks = c("Hossfeld IV (MAE, RMSE, ME)", "Chapman-Richards (ΔAICn)", "Sharma-Parton physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$gamPhysio), y = thpl2016$TotalHt, color = "REML GAM physio", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation)) +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$gam), y = thpl2016$TotalHt, color = "REML GAM (all but MAB)", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation)) +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$parabolic), y = thpl2016$TotalHt, color = "parabolic (MAB)", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  geom_line(aes(x = thplReferenceDbh, y = 1.37 + 0.3048 * exp(7.232880669 - 5.746899904 * (0.393701 * thplReferenceDbh)^-0.271564741), linetype = "previous model"), color = "grey70") + # Hanus et al. 1999, Eq. 1
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[4])~"western redcedar DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[4]))~"western redcedar DBH")) +
  scale_color_manual(breaks = c("REML GAM (all but MAB)", "parabolic (MAB)", "REML GAM physio", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 08 UMCA-THPL curves.png", height = 12, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 08 UMCA-THPL curves.tif", height = 12, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 08 UMCA-THPL curves.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 9: preferred models for other species group
# preferred forms from 10x10 cross validation
#                height                                           DBH
#                base              generalized                    base                    generalized
# other species  Curtis            REML GAM BA+L                  REML GAM                Sibbesen form physio
#                power             REML GAM BAL+L physio          parabolic               Chapman-Richards form RelHt
#                Korf                                             linear
#print(preferredForms %>% filter(species == "other species") %>% select(-mabName, -aucMab, -nseName, -aucNse) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "ALRU2", "TSHE", "ALRU2", "ACMA3", "UMCA", "THPL", "other"), levels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")), maeName = str_trunc(maeName, 28, ellipsis = ""), rmseName = str_trunc(rmseName, 28, ellipsis = ""), aicName = str_trunc(aicName, 28, ellipsis = "")), n = 16)
if (exists("otherHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/other preferred models.Rdata") }

# no prior height model for other species, so no reference or reference DBH
otherReference = bind_rows(bind_rows(otherHeightFromDiameterPreferred$gam$stats %>% mutate(model = "base form 1"),
                                     otherHeightFromDiameterPreferred$linear$stats %>% mutate(model = "base form 2"),
                                     otherHeightFromDiameterPreferred$gamRelDbhPhysio$stats %>% mutate(model = "generalized height")) %>% 
                             mutate(responseVariable = "height"),
                           bind_rows(otherDiameterFromHeightPreferred$gam$stats %>% mutate(model = "base form 1"),
                                     otherDiameterFromHeightPreferred$sibbesenReplace$stats %>% mutate(model = "base form 2"),
                                     otherDiameterFromHeightPreferred$gamRelHt$stats %>% mutate(model = "generalized DBH")) %>% 
                             mutate(responseVariable = "DBH")) %>%
  mutate(fitSet = "primary", species = "other", adaptiveWeightFraction = 0)

ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = other2016physio$DBH, y = predict(otherHeightFromDiameterPreferred$gamRelDbhPhysio), color = "REML GAM RelDbh physio", group = other2016physio$isPlantation, linetype = other2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$gam), color = "REML GAM (MAE, RMSE, ME)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$linear), color = "linear (ΔAICn)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(.(plotLetters[1])~"other species height")) +
  #labs(x = "DBH, cm", y = "height, m", color = NULL, title = bquote(bold(.(plotLetters[1]))~"other species height")) +
  scale_color_manual(breaks = c("REML GAM (MAE, RMSE, ME)", "linear (ΔAICn)", "REML GAM RelDbh physio", "previous model"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$gamRelHt), y = other2016$TotalHt, color = "REML GAM RelHt", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$gam), y = other2016$TotalHt, color = "REML GAM (all five statistics)", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$sibbesenReplace), y = other2016$TotalHt, color = "Sibbesen replace", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = predict(otherDiameterFromHeightPreferred$parabolic), y = other2016$TotalHt, color = "parabolic", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  #geom_line(aes(x = psmeReferenceDbh, y = 1.3 + exp(5.7567 - 6.7792*psmeReferenceDbh^-0.2795), linetype = "previous model"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(color = guide_legend(position = "inside", override.aes = list(alpha = 0.8)), linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(.(plotLetters[2])~"other species DBH")) +
  #labs(x = "DBH, cm", y = NULL, color = NULL, title = bquote(bold(.(plotLetters[2]))~"other species DBH")) +
  scale_color_manual(breaks = c("REML GAM (all five statistics)", "Sibbesen replace", "REML GAM RelHt", "previous model"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position.inside = c(1, 0.01)) +
get_preferred_model_linetype_legend(noPreviousModel = TRUE) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n33", heights = c(1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "previous model"), labels = c("natural regeneration", "plantation", "previous model"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 09 other species curves.png", height = 17/3 + 0.5, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 09 other species curves.tif", height = 17/3 + 0.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 09 other species curves.pdf", height = 17/3 + 0.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure 10: comparison of goodness of fit statistics between reference and revised models
selectedModels = bind_rows(psmeReference, alruReference, tsheReference, acmaReference, umcaReference, thplReference, otherReference) %>% 
  select(-coefficients) %>%
  mutate(deltaAicN = aic/nValidation - min(aic/nValidation, na.rm = TRUE),
         model = as.factor(model),
         species = factor(species, labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"), levels = c("PSME", "ALRU2", "TSHE", "ACMA3", "UMCA", "THPL", "other")))
#write_xlsx(list(preferredModelStats = selectedModels %>% select(species, responseVariable, name, mab, mapb, mae, mape, rmse, rmpse, deltaAicN, nse), "trees/height-diameter/figures/selected model stats.xlsx")
# selectedModels %>% filter(species == "bigleaf maple") %>% select(responseVariable, name, mab, mae, rmse, deltaAicN, nse)

ggplot() +
  geom_point(aes(x = mab, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAB, m", y = NULL, color = NULL, shape = NULL, size = NULL, title = paste0(plotLetters[1], " height prediction")) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(limits = rev) +
ggplot() +
  geom_point(aes(x = mae, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAE, m", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = rmse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "RMSE, m", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = deltaAicN, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "ΔAICn", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 6)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = nse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "height"), alpha = 0.6) +
  coord_cartesian(xlim = c(0.1, 1)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "model efficiency", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1.0), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = mab, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAB, cm", y = NULL, color = NULL, shape = NULL, size = NULL, title = paste0(plotLetters[2], " DBH prediction")) +
  scale_y_discrete(limits = rev) +
ggplot() +
  geom_point(aes(x = mae, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "MAE, cm", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = rmse, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(0, NA)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "RMSE, cm", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = deltaAicN, y = species, color = model, shape = model, size = model), selectedModels %>% filter(responseVariable == "DBH"), alpha = 0.6, na.rm = TRUE) + # NaN in red alder
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  guides(color = "none", shape = "none", size = "none") +
  labs(x = "ΔAICn", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 6)) +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot() +
  geom_point(aes(x = nse, y = species, color = model, shape = model, size = model), selectedModels %>% filter((responseVariable == "DBH") | ((species == "bigleaf maple") & (model == "generalized height"))) %>% mutate(nse = if_else(model != "generalized height", nse, -1)), alpha = 0.6) + # flow fake height data point off plot to get ggplot2 3.5.1 to draw the generalized height legend key, for some reason inserting a fake row with add_row() causes ggplot to draw additional red alder points
  annotate("text", x = 0.06, y = "bigleaf maple", label = "🡄 −1.47", alpha = 0.8, color = "grey60", hjust = 0, size = 2.0, vjust = 0.4) +
  coord_cartesian(xlim = c(0.1, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 0.8))) +
  labs(x = "model efficiency", y = NULL, color = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1.0), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, guides = "collect") &
  scale_color_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c("dodgerblue2", "red2", "green2", "cyan", "grey60"), drop = FALSE) &
  scale_shape_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c(15, 18, 16, 16, 17), drop = FALSE) &
  scale_size_manual(breaks = c("base form 1", "base form 2", "generalized height", "generalized DBH", "previous model"), values = c(1.9, 2.5, 2.2, 2.2, 1.7), drop = FALSE) &
  theme(plot.title.position = "plot")
#ggsave("trees/height-diameter/figures/Figure 10 model accuracy.png", height = 9, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure 10 model accuracy.tif", height = 9, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure 10 model accuracy.pdf", height = 9, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure S1: comparison of accuracy metrics
accuracyCorrelation = bind_rows(as.data.frame(cor(primaryResults %>% filter(responseVariable == "height") %>% 
                                                    select(mapb, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "height"),
                                as.data.frame(cor(primaryResults %>% filter(responseVariable == "DBH") %>% 
                                                    select(mapb, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "DBH")) %>%
  mutate(metricX = factor(metricX, levels = c("mapb", "mape", "rmse", "deltaAicN", "nse"), labels = c("MAB", "MAE", "RMSE", "ΔAICn", "model\nefficiency")),
         metricY = factor(metricY, levels = c("mapb", "mape", "rmse", "deltaAicN", "nse"), labels = c("MAB", "MAE", "RMSE", "ΔAICn", "model efficiency")))

ggplot(accuracyCorrelation %>% filter(responseVariable == "height")) + 
  coord_equal() +
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  labs(x = NULL, y = NULL, fill = "correlation", title = bquote(.(plotLetters[1])~"height prediction")) +
  #labs(x = NULL, y = NULL, fill = "correlation", title = bquote(bold(.(plotLetters[1]))~"height prediction")) +
  scale_y_discrete(limits = rev) +
ggplot(accuracyCorrelation %>% filter(responseVariable == "DBH")) + 
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "correlation", title = bquote(.(plotLetters[2])~"DBH prediction")) +
  #labs(x = NULL, y = NULL, fill = "correlation", title = bquote(bold(.(plotLetters[2]))~"DBH prediction")) +
  scale_y_discrete(labels = NULL, limits = rev) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
  plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  scale_fill_scico(palette = "vik", limits = c(-1, 1)) &
  theme(legend.spacing.y = unit(0.5, "line"))
#ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.png", height = 8.5, width = 20, units = "cm")
#ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.tif", height = 8.5, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure S01 goodness of fit correlation.pdf", height = 8.5, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


# Figure S2: height prediction accuracy
# If axis limits, vertical dodges, or aesthetics are changed Figure S6 should be updated.
heightFromDiameterAccuracyLevels = primaryResults %>% 
  filter(responseVariable == "height") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmspe = mean(if_else(is.na(rmspe), 100, rmse)),
            .groups = "drop_last") %>% # change to grouping only by name
  summarize(weightedMae = sum(speciesFraction * meanPenalizedMae),
            .groups = "drop") %>%
  arrange(weightedMae)
heightFromDiameterResults = primaryResults %>% 
  filter(responseVariable == "height", str_detect(name, "GNLS") == FALSE, str_detect(name, "RelHt") == FALSE) %>%
  mutate(name = factor(name, levels = heightFromDiameterAccuracyLevels$name)) %>%
  group_by(species, name) %>%
  reframe(quantiles = c(0.025, 0.5, 0.975),
          mapb = quantile(mapb, probs = quantiles, na.rm = TRUE), 
          mape = quantile(mape, probs = quantiles, na.rm = TRUE), 
          rmspe = quantile(rmspe, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE), 
          sizeShapeAlpha = sizeShapeAlpha[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("mapb", "mape", "rmspe", "deltaAicN", "nse")) %>%
  mutate(verticalDodge = recode(species, "Douglas-fir" = 0.3, "red alder" = 0.2, "western hemlock" = 0.1, "bigleaf maple" = 0.0, "Oregon myrtle" = -0.1, "western redcedar" = -0.2, "other species" = -0.3))

ggplot(heightFromDiameterResults) +
  geom_errorbarh(aes(xmin = mapb_q0.025, xmax = mapb_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = mapb_q0.5, y = name, color = species, alpha = sizeShapeAlpha, size = sizeShapeAlpha, shape = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "MAB, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbarh(aes(xmin = mape_q0.025, xmax = mape_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = mape_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbarh(aes(xmin = rmspe_q0.025, xmax = rmspe_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = rmspe_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "RMSE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  #scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  #geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM", yend = "REML GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM BA+L", yend = "REML GAM BA+L", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM physio", yend = "REML GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = deltaAicN_q0.025, xmax = deltaAicN_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = deltaAicN_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  labs(x = "ΔAICn", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  coord_cartesian(xlim = c(0, 10)) + # exclude high AIC of linear, parabolic, and Douglas-fir power+Curtis fits to avoid squashing of nonlinear differences
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbarh(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0.1, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 1), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(0.75, 0.75, 0.3)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(18, 16, 3)) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(1.5, 1.9, 1.4)) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure S02 height accuracy.png", height = 16, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure S02 height accuracy.tif", height = 22, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure S02 height accuracy.pdf", height = 22, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


# Figure S3: DBH prediction accuracy
# If axis limits, vertical dodges, or aesthetics are changed Figure S4 should be updated.
diameterFromHeightAccuracyLevels = primaryResults %>% 
  filter(responseVariable == "DBH") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmspe = mean(if_else(is.na(rmspe), 100, rmspe)),
            .groups = "drop_last") %>% # change to grouping only by name
  summarize(weightedRmse = sum(speciesFraction * meanPenalizedRmspe), 
            .groups = "drop") %>%
  arrange(weightedRmse)
diameterFromHeightResults = primaryResults %>% 
  filter(responseVariable == "DBH", str_detect(name, "BA\\+L") == FALSE, str_detect(name, "GNLS") == FALSE, name != "REML GAM ABA+T physio") %>%
  mutate(name = factor(name, levels = diameterFromHeightAccuracyLevels$name)) %>%
  group_by(species, name) %>%
  reframe(quantiles = c(0.025, 0.5, 0.975),
          mapb = quantile(mapb, probs = quantiles, na.rm = TRUE), 
          mape = quantile(mape, probs = quantiles, na.rm = TRUE), 
          rmspe = quantile(rmspe, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE), 
          sizeShapeAlpha = sizeShapeAlpha[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("mapb", "mape", "rmspe", "deltaAicN", "nse")) %>%
  mutate(verticalDodge = recode(species, "Douglas-fir" = 0.3, "red alder" = 0.2, "western hemlock" = 0.1, "bigleaf maple" = 0.0, "Oregon myrtle" = -0.1, "western redcedar" = -0.2, "other species" = -0.3))

ggplot(diameterFromHeightResults) +
  geom_errorbarh(aes(xmin = mapb_q0.025, xmax = mapb_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = mapb_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "MAB, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbarh(aes(xmin = mape_q0.025, xmax = mape_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = mape_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbarh(aes(xmin = rmspe_q0.025, xmax = rmspe_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = rmspe_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "RMSE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  #geom_segment(x = 0.255, xend = 0.275, y = "Schnute inverse", yend = "Schnute inverse", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM", yend = "REML GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM ABA+T", yend = "REML GAM ABA+T", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM physio", yend = "REML GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = deltaAicN_q0.025, xmax = deltaAicN_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = deltaAicN_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "ΔAICn", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), trans = scales::pseudo_log_trans()) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbarh(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0.1, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.4, 0.7, 1), minor_breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 1), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(0.75, 0.75, 0.3)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(16, 18, 3)) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), labels = c("NLME/nlrob", "form significant", "not significant"), values = c(1.5, 1.9, 1.4)) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.png", height = 16, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.tif", height = 22, width = 20, units = "cm", dpi = figureDpi, compression = "lzw+p")
ggsave("trees/height-diameter/figures/Figure S03 DBH accuracy.pdf", height = 22, width = 20, units = "cm", dpi = figureDpi, device = cairo_pdf, fallback_resolution = figureDpi)


## Figure S4 in residuals.R


## supplementary spreadsheet: model forms by species and use 
preferredModelForms = heightDiameterModelRanking %>% filter(significant) %>% 
  mutate(species = fct_relevel(species, c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"))) %>%
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-fitting, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, isBaseForm, significant, aucMab, aucMae, aucRmse, aucDeltaAicN, aucNse)
preferredModelFormsHeightNonPhysio = heightDiameterModelRanking %>% filter(significant, hasPhysio == FALSE) %>% 
  mutate(species = fct_relevel(species, c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"))) %>%
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-fitting, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, isBaseForm, significant, aucMab, aucMae, aucRmse, aucDeltaAicN, aucNse)
preferredModelFormsInitialDbh = heightDiameterModelRanking %>% filter(significant, hasStand == FALSE) %>% 
  mutate(species = fct_relevel(species, c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"))) %>%
  group_by(responseVariable, species) %>% 
  arrange(desc(aucBlended)) %>% 
  mutate(rank = row_number()) %>% 
  arrange(desc(responseVariable), species, isBaseForm, rank) %>%
  select(-fitting, -starts_with("has"), -speciesFraction) %>%
  relocate(responseVariable, species, rank, name, isBaseForm, significant, aucMab, aucMae, aucRmse, aucDeltaAicN, aucNse)
#write_xlsx(list(preferred = preferredModelForms, heightNonPhysio = preferredModelFormsHeightNonPhysio, initialDbh = preferredModelFormsInitialDbh), "trees/height-diameter/figures/Elliott height-diameter model AUCs.xlsx")


## additional info
# nonphysical predictions
heightDiameterCoefficients %>% summarize(pctFitsWithNonPhysical = 100 * sum(nNonPhysical > 0) / n())
nonPhysical = heightDiameterCoefficients %>% filter(fitSet == "primary", significant) %>% 
  group_by(responseVariable, species, name) %>%
  summarize(baseName = baseName[1], isBaseForm = isBaseForm[1], pctWithNonPhysical = 100 * sum(nNonPhysical > 0) / n(), pctNonPhysical = mean(100 * nNonPhysical / nValidation), .groups = "drop")

nonPhysical %>% group_by(responseVariable, species) %>% 
  summarize(pBaseNonPhysical = sum(isBaseForm * pctNonPhysical > 0) / sum(isBaseForm),
            pGeneralizedNonPhysical = sum((isBaseForm == FALSE) * pctNonPhysical > 0) / sum(isBaseForm == FALSE))

print(nonPhysical %>% group_by(responseVariable) %>% # nonzero for 157 of 400 significant combinations
  slice_max(pctNonPhysical, n = 20) %>%
  arrange(desc(responseVariable), desc(pctNonPhysical)), n = 40)
print(nonPhysical %>% filter(pctNonPhysical > 0) %>% group_by(responseVariable) %>%
        slice_min(pctNonPhysical, n = 20) %>%
        arrange(desc(responseVariable), desc(pctNonPhysical)), n = 40)
print(nonPhysical %>% filter(isBaseForm, pctNonPhysical == 0), n = 250)
print(nonPhysical %>% filter(baseName == "Ruark"), n = 30)

ggplot() +
  geom_histogram(aes(x = pctWithNonPhysical, fill = responseVariable, group = responseVariable), nonPhysical, binwidth = 1) +
  coord_trans(y = scales::pseudo_log_trans()) +
  labs(x = "fraction of fits with non-physical trees, %", y = "model fits", fill = NULL) +
  scale_y_continuous(breaks = c(0, 1, 10, 100, 1000))

ggplot() +
  stat_ecdf(aes(x = pctNonPhysical, color = responseVariable, group = responseVariable), nonPhysical) +
  labs(x = "non-physical fraction of trees, %", y = "fraction of fits", color = NULL)

# AUC principal components
aucComponents = prcomp(~ aucMab + aucMae + aucRmse + aucDeltaAicN + aucNse, heightDiameterModelAucs)
summary(aucComponents) # 82% in first component: 0.5 MAB + 0.5 MAE + 0.49 RMSE + 0.13 AIC + 0.49 model efficiency

# GAM effective degrees of freedom
ggplot() +
  geom_boxplot(aes(x = effectiveDegreesOfFreedom, y = name, color = species, group = paste(name, species)), heightDiameterResults %>% filter(fitSet == "primary", fitting == "gam"), alpha = 0.2) +
  coord_trans(x = scales::transform_pseudo_log()) +
  guides(color = guide_legend(override.aes = list(alpha = 0.9))) +
  labs(x = "effective degrees of freedom", y = NULL, color= NULL) +
  scale_color_manual(breaks = levels(heightDiameterResults$species), limits = levels(heightDiameterResults$species), values = speciesGroupColors) +
  scale_x_continuous(breaks = c(0, 1, 2, 10, 20, 100, 200, 1000), minor_breaks = c(3, 4, 5, 6, 7, 8, 9, 30, 40, 50, 60, 70, 80, 300, 400, 500, 600, 700, 800, 900))