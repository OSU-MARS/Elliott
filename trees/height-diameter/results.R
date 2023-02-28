# load data from setup.R, get regressions and summaries from <species>.R
handlers(global = TRUE)
handlers("progress")

if (exists("psmeResults") == FALSE) { load("trees/height-diameter/data/PSME results.Rdata") }
if (exists("alruResults") == FALSE) { load("trees/height-diameter/data/ALRU2 results.Rdata") }
if (exists("tsheResults") == FALSE) { load("trees/height-diameter/data/TSHE results.Rdata") }
if (exists("acmaResults") == FALSE) { load("trees/height-diameter/data/ACMA3 results.Rdata") }
if (exists("thplResults") == FALSE) { load("trees/height-diameter/data/THPL results.Rdata") }
if (exists("umcaResults") == FALSE) { load("trees/height-diameter/data/UMCA results.Rdata") }
if (exists("otherResults") == FALSE) { load("trees/height-diameter/data/other results.Rdata") }

# form regression summary
heightDiameterCoefficients = bind_rows(psmeCoefficients, alruCoefficients, tsheCoefficients, acmaCoefficients,
                                       umcaCoefficients, thplCoefficients, otherCoefficients) %>% 
  relocate(responseVariable, species, fitSet, fixedWeight, name, fitting, repetition, fold, a0, a1, a1p, a2, a2p, a3, a3p, a4, a5, a6, a7, a8, a9, a9p, b1, b1p, b2, b2p, b3, b3p, b4, b4p)
#write_xlsx(heightDiameterCoefficients, "trees/height-diameter/HtDia coefficients.xlsx")

heightDiameterResults = bind_rows(psmeResults, alruResults, tsheResults, acmaResults,
                                  umcaResults, thplResults, otherResults) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"), levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         speciesFraction = recode(species, "Douglas-fir" = 0.750, "red alder" = 0.101, "western hemlock" = 0.056, "bigleaf maple" = 0.029, "Oregon myrtle" = 0.025, "western redcedar" = 0.013, "other species" = 0.017),
         isBaseForm = (str_detect(name, "Sharma-") == FALSE) & (str_detect(name, "ABA\\+T") == FALSE) & (str_detect(name, "BA\\+L") == FALSE) & (str_detect(name, "physio") == FALSE) & (str_detect(name, "RelHt") == FALSE),
         hasPhysio = str_detect(name, "physio"),
         hasStand = str_detect(name, "ABA\\+T") | str_detect(name, "BA\\+L") | str_detect(name, "RelHt"),
         weighting = if_else(fitting %in% c("gnls", "nlrob"), "reweighted", "fixed weights"),
         sizeShapeAlpha = as.factor(if_else(significant == TRUE, weighting, "not significant"))) %>%
  group_by(fitSet, fixedWeight, responseVariable, species) %>%
  mutate(nFits = n()) %>%
  ungroup()

speciesGroupColors = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")

# extract subset of results for analysis
primaryResults = heightDiameterResults %>% 
  filter(fitSet == "primary", is.na(fixedWeight), 
         str_detect(name, "GNLS") == FALSE,  # exclude fits from gnls(), for now can't filter on $fitting as it was left NA for failed fits and results files haven't been rebuilt with the fix
         (responseVariable != "height") | (str_detect(name, "RelHt") == FALSE), # exclude height control forms using relative height
         (responseVariable != "DBH") | (str_detect(name, "BA\\+L") == FALSE)) %>% # exclude diameter control forms using basal area
  group_by(responseVariable, species) %>%
  mutate(deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) %>% # ΔAIC within response variable and species, needed for AUCs and figures
  ungroup()
#print(primaryResults %>% group_by(fitSet, responseVariable) %>% summarize(n = n(), names = unique(name)), n = 60)
#primaryResults %>% group_by(fitSet, species) %>% summarize(deltaAicN = sum(is.na(deltaAicN)), mab = sum(is.na(mab)), mae = sum(is.na(mae)), pearson = sum(is.na(pearson)), nse = sum(is.na(nse)), rmse = sum(is.na(rmse)))

# rank model forms by estimated prediction ability (using AUC) for form selection
# This is somewhat computationally intensive, hence the progress bar, but still only takes a minute or so single threaded.
with_progress({
  crossValidatedModelCount = primaryResults %>% group_by(responseVariable, species) %>% summarize(n = length(unique(name)), .groups = "drop")
  progressBar = progressor(steps = sum(crossValidatedModelCount$n))
  
  heightDiameterModelAucs = primaryResults %>% 
    group_by(responseVariable, species, name) %>% # for now, look across fit set and weighting as GAMs are only in the primary set with estimated weights
    group_modify(~{
      if ((nrow(.x) == 1) | all(is.na(.x$nse)))
      {
        # no distribution to compare to since this model has only a no fit result or wasn't cross validated
        progressBar(str_pad(paste(.y$responseVariable, .y$species, .y$name), 60, "right"))
        return(tibble(otherModelName = NA_character_, fitting = .x$fitting[1], isBaseForm = .x$isBaseForm[1], hasPhysio = .x$hasPhysio[1], hasStand = .x$hasStand[1],
                      aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucPearson = NA_real_, aucNse = NA_real_, aucRmse = NA_real_))
      }
      
      # get all other cross-validation results for this response variable and species
      # Assumes no names are shared across fittings in the results set.
      matchingResults = primaryResults %>% filter(responseVariable == .y$responseVariable, species == .y$species)
      matchingModelNames = unique(matchingResults$name)
      pairwiseAucs = bind_rows(lapply(matchingModelNames, function(otherModelName) 
        {
          otherModelResults = matchingResults %>% filter(responseVariable == .y$responseVariable, species == .y$species, name == otherModelName)
          
          if ((nrow(otherModelResults) == 1) | all(is.na(otherModelResults$nse)))
          {
            # also no distribution to compare to if the other model has only a no fit result or wasn't cross validated
            return(tibble(otherModelName = otherModelName, fitting = .x$fitting[1], isBaseForm = .x$isBaseForm[1], hasPhysio = .x$hasPhysio[1], hasStand = .x$hasStand[1],
                          aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucPearson = NA_real_, aucNse = NA_real_, aucRmse = NA_real_))
          }
          
          # find AUCs: for species with few samples (e.g. THPL) mean absolute bias may be all NA, which case WeightedROC() errors
          # AUC is the probability a sample from one empirical distribution is greater than a sample from an another 
          # distribution (see AUC.R). If the labels are reversed, the probability changes to the sample being less than.
          # The AUCs found in the group_modify() below are thus the probability a base form improves on a base GAM.
          # Base GAMs are retained as controls and, since the GAM's goodness of fit metrics are compared with themselves,
          # should always have AUCs of 0.5.
          lowerIsBetter = factor(c(rep(0, nrow(.x)), rep(1, nrow(otherModelResults))), levels = c(0, 1))
          higherIsBetter = factor(c(rep(1, nrow(.x)), rep(0, nrow(otherModelResults))), levels = c(0, 1))
          
          # estimate AUC for bias based on what data is available, which is potentially none for species with few measurements
          availableMabData = tibble(guess = c(.x$mab, otherModelResults$mab), label = lowerIsBetter) %>% drop_na()
          aucMab = NA_real_
          if (nrow(availableMabData) > 0)
          {
            aucMab = WeightedAUC(WeightedROC(guess = availableMabData$guess, label = availableMabData$label))
          }

          if ((sum(is.na(.x$mae)) + sum(is.na(otherModelResults$mae)) + sum(is.na(.x$pearson)) + sum(is.na(otherModelResults$pearson)) + 
               sum(is.na(.x$nse)) + sum(is.na(otherModelResults$nse)) + sum(is.na(.x$rmse)) + sum(is.na(otherModelResults$rmse))) > 0)
          {
            # WeightedROC() errors on NAs in its arguments but can't do so informatively
            # Fail informatively instead so that investigation is possible.
            stop(paste0(.y$species[1], " ", .y$responseVariable[1], " ", .y$name[1], " (", nrow(.x), ") x ", otherModelName, " (", nrow(otherModelResults), "):",
                        " mab ", sum(is.na(.x$mab)), " ", sum(is.na(otherModelResults$mab)),
                        " mae ", sum(is.na(.x$mae)), " ", sum(is.na(otherModelResults$mae)),
                        " pearson ", sum(is.na(.x$pearson)), " ", sum(is.na(otherModelResults$pearson)),
                        " nse ", sum(is.na(.x$nse)), " ", sum(is.na(otherModelResults$nse)),
                        " rmse ", sum(is.na(.x$rmse)), " ", sum(is.na(otherModelResults$rmse))), quote = FALSE)
          }
          return(tibble(otherModelName = otherModelName, fitting = .x$fitting[1], isBaseForm = .x$isBaseForm[1],
                        hasPhysio = .x$hasPhysio[1], hasStand = .x$hasStand[1],
                        aucDeltaAicN = WeightedAUC(WeightedROC(guess = c(.x$deltaAicN, otherModelResults$deltaAicN), label = lowerIsBetter)),
                        aucMab = aucMab,
                        aucMabN = nrow(availableMabData),
                        aucMae = WeightedAUC(WeightedROC(guess = c(.x$mae, otherModelResults$mae), label = lowerIsBetter)),
                        aucPearson = WeightedAUC(WeightedROC(guess = c(.x$pearson, otherModelResults$pearson), label = higherIsBetter)),
                        aucNse = WeightedAUC(WeightedROC(guess = c(.x$nse, otherModelResults$nse), label = higherIsBetter)),
                        aucRmse = WeightedAUC(WeightedROC(guess = c(.x$rmse, otherModelResults$rmse), label = lowerIsBetter))))
        }))
      
      progressBar(str_pad(paste(.y$responseVariable, .y$species, .y$name), 60, "right"))
      return(pairwiseAucs)
    }) %>%
    ungroup()
})
#heightDiameterModelAucs %>% group_by(responseVariable, species) %>%
#  summarize(matrixElements = n(), mabN = sum(is.na(aucMab) == FALSE), mabInputN = sum(aucMabN, na.rm = TRUE), maeN = sum(is.na(aucMae) == FALSE), rmseN = sum(is.na(aucRmse) == FALSE), nseN = sum(is.na(aucNse) == FALSE), pearsonN = sum(is.na(aucNse) == FALSE), .groups = "drop")

# take mean AUCs over otherModelName, excluding self and unsuccessful fits
heightDiameterModelRanking = heightDiameterModelAucs %>% 
  group_by(responseVariable, species, name) %>%
  summarize(fitting = fitting[1], isBaseForm = isBaseForm[1], hasPhysio = hasPhysio[1], hasStand = hasStand[1],
            aucDeltaAicN = median(if_else(name != otherModelName, aucDeltaAicN, NA_real_), na.rm = TRUE),
            aucMab = median(if_else(name != otherModelName, aucMab, NA_real_), na.rm = TRUE),
            aucMae = median(if_else(name != otherModelName, aucMae, NA_real_), na.rm = TRUE),
            aucPearson = median(if_else(name != otherModelName, aucPearson, NA_real_), na.rm = TRUE),
            aucNse = median(if_else(name != otherModelName, aucNse, NA_real_), na.rm = TRUE),
            aucRmse = median(if_else(name != otherModelName, aucNse, NA_real_), na.rm = TRUE),
            .groups = "drop_last") %>%
  mutate(aucBlended = 0.5 * if_else(responseVariable == "height", aucMae, aucRmse) + 0.2 * if_else(responseVariable == "DBH", aucMae, aucRmse) + 0.3 * aucDeltaAicN) %>%
  ungroup()
heightDiameterModelRanking %>% group_by(responseVariable, species) %>%
  summarize(n = n(), mabN = sum(is.na(aucMab) == FALSE), maeN = sum(is.na(aucMae) == FALSE), rmseN = sum(is.na(aucRmse) == FALSE), nseN = sum(is.na(aucNse) == FALSE), pearsonN = sum(is.na(aucNse) == FALSE), .groups = "drop")

# find preferred model forms
nPreferredModelForms = 4
preferredForms = full_join(full_join(full_join(heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMab, n = nPreferredModelForms) %>% arrange(desc(aucMab)) %>% mutate(mabName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, mabName, aucMab),
                                               heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucMae, n = nPreferredModelForms) %>% arrange(desc(aucMae)) %>% mutate(maeName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, maeName, aucMae),
                                               by = c("responseVariable", "species", "isBaseForm", "rank")), # full join since western redcedar mean absolute bias is all NA
                                     full_join(heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucRmse, n = nPreferredModelForms) %>% arrange(desc(aucRmse)) %>% mutate(rmseName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, rmseName, aucRmse),
                                               heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucDeltaAicN, n = nPreferredModelForms) %>% arrange(desc(aucDeltaAicN)) %>% mutate(aicName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, aicName, aucDeltaAicN),
                                               by = c("responseVariable", "species", "isBaseForm", "rank")),
                                     by = c("responseVariable", "species", "isBaseForm", "rank")),
                           full_join(full_join(heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucNse, n = nPreferredModelForms) %>% arrange(desc(aucNse)) %>% mutate(nseName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, nseName, aucNse),
                                               heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucPearson, n = nPreferredModelForms) %>% arrange(desc(aucPearson)) %>% mutate(pearsonName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, pearsonName, aucPearson),
                                               by = c("responseVariable", "species", "isBaseForm", "rank")),
                                     heightDiameterModelRanking %>% group_by(responseVariable, species, isBaseForm) %>% slice_max(aucBlended, n = nPreferredModelForms) %>% arrange(desc(aucBlended)) %>% mutate(blendedName = name, rank = row_number()) %>% select(responseVariable, species, isBaseForm, rank, blendedName, aucBlended),
                                     by = c("responseVariable", "species", "isBaseForm", "rank")),
                           by = c("responseVariable", "species", "isBaseForm", "rank")) %>%
  arrange(desc(responseVariable), species, desc(isBaseForm), rank) %>%
  relocate(responseVariable, species, isBaseForm, rank) %>%
  ungroup()
#preferredForms %>% group_by(responseVariable) %>% mutate(speciesPresent = length(unique(species))) %>% group_by(responseVariable, species) %>% summarize(speciesPresent = speciesPresent[1], nPreferredForms = n(), .groups = "drop")


## summary for Abstract
# form counts, changes in model efficiency with inclusion of generalizing predictors
heightDiameterResults %>% filter(fitting %in% c("nlrob", "gsl_nls"), responseVariable != "height" | str_detect(name, "RelHt") == FALSE) %>% 
  mutate(baseName = if_else(word(name) %in% c("REML", "modified", "unified"), paste(word(name, 1), word(name, 2)), word(name)),
         deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) %>% # unrestricted ΔAICn
  group_by(responseVariable, species, baseName) %>%
  mutate(nseBase = first(na.omit(if_else(isBaseForm, nse, NA_real_))),
         nseStandDelta = if_else(hasStand & (hasPhysio == FALSE), nse, NA_real_) - nseBase,
         nsePhysioDelta = if_else(hasPhysio & (hasStand == FALSE), nse, NA_real_) - nseBase,
         nseCombinedDelta = if_else(hasPhysio & hasStand, nse, NA_real_) - nseBase) %>%
  group_by(fitSet, fixedWeight, responseVariable, species) %>%
  slice_max(nse, n = 3 * 10 * 10) %>% # nmodels * nfolds * nrepetitions
  group_by(fitSet, fixedWeight, responseVariable) %>%
  summarize(nForms = length(unique(name)),
            nFits = n(),
            mae = mean(mae),
            mape = mean(mape),
            rmse = mean(rmse),
            deltaAicN = mean(deltaAicN),
            nse = mean(nse),
            pearson = mean(pearson),
            nseBase = median(if_else((is.na(nseStandDelta) == FALSE) | (is.na(nsePhysioDelta) == FALSE), nseBase, NA_real_), na.rm = TRUE),
            nseStandDelta = median(nseStandDelta, na.rm = TRUE),
            nsePhysioDelta = median(nsePhysioDelta, na.rm = TRUE),
            nseCombinedDelta = median(nseCombinedDelta, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(fitSet), desc(responseVariable), is.na(fixedWeight) == FALSE)


## summaries for Results
# fraction of nonlinear base forms which are more accurate than a base GAM
improvementThresholdProbability = 0.5
heightDiameterModelAucs %>% filter(name != "REML GAM", otherModelName == "REML GAM") %>%
  group_by(responseVariable, species) %>% # restriction to primary fit set and base forms is above
  summarize(deltaAic = 100 * sum(aucDeltaAicN > improvementThresholdProbability) / n(),
            mae = 100 * sum(aucMae > improvementThresholdProbability) / n(),
            pearson = 100 * sum(aucPearson > improvementThresholdProbability) / n(),
            nse = 100 * sum(aucNse > improvementThresholdProbability) / n(),
            rmse = 100 * sum(aucRmse > improvementThresholdProbability) / n(),
            .groups = "drop") %>%
  arrange(desc(responseVariable))

heightDiameterCoefficients %>% filter(fitting %in% c("lm", "nlrob", "gsl_nls")) %>% # exclude GAMs
  #group_by(responseVariable) %>%
  summarise(n = n(),
            anyP = sum((is.na(a1p) == FALSE) | (is.na(a2p) == FALSE) | (is.na(a3p) == FALSE) | (is.na(a4p) == FALSE) | (is.na(b1p) == FALSE) | (is.na(b2p) == FALSE) | (is.na(b3p) == FALSE) | (is.na(b4p) == FALSE)),
            a1 = sum(is.na(a1) == FALSE), a1p = sum(is.na(a1p) == FALSE),
            a2 = sum(is.na(a2) == FALSE), a2p = sum(is.na(a2p) == FALSE),
            a3 = sum(is.na(a3) == FALSE), a3p = sum(is.na(a3p) == FALSE),
            a4 = sum(is.na(a4) == FALSE), a4p = sum(is.na(a4p) == FALSE), # plantation effects (a5-a8) not tested
            a9 = sum(is.na(a9) == FALSE), a9p = sum(is.na(a9p) == FALSE), 
            b1 = sum(is.na(b1) == FALSE), b1p = sum(is.na(b1p) == FALSE),
            b2 = sum(is.na(b2) == FALSE), b2p = sum(is.na(b2p) == FALSE),
            b3 = sum(is.na(b3) == FALSE), b3p = sum(is.na(b3p) == FALSE),
            b4 = sum(is.na(b3) == FALSE), b4p = sum(is.na(b4p) == FALSE)) %>%
  mutate(total = a1 + a2 + a3 + a4 + a9 + b1 + b2 + b3,
         totalP = a1p + a2p + a3p + a4p + a9 + b1p + b2p + b3p + b4p,
         hasPpct = 100 * anyP/a1)

# AIC selection of generalizing predictors based on ΔAICn AUC
improvementThresholdProbability = 0.5
heightDiameterModelAucs %>% filter(isBaseForm == FALSE, fitting %in% c("nlrob", "gsl_nls")) %>%
  group_by(responseVariable) %>%
  summarize(n = n(), 
            baseBetter = sum(aucDeltaAicN <= improvementThresholdProbability), # base form is ΔAICn preferable
            generalizedBetter = sum(aucDeltaAicN > improvementThresholdProbability), # a generalized form is preferable
            physioDeltaAicNmed = median(if_else(hasPhysio & (hasStand == FALSE), aucDeltaAicN, NA_real_), na.rm = TRUE), 
            standDeltaAicNmed = median(if_else((hasPhysio == FALSE) & hasStand, NA_real_, aucDeltaAicN), na.rm = TRUE),
            bothDeltaAicNmed = median(if_else(hasPhysio & hasStand, NA_real_, aucDeltaAicN), na.rm = TRUE)) %>%
  mutate(baseBetterPct = 100 * baseBetter / n, 
         generalizedBetterPct = 100 * generalizedBetter / n) %>%
  arrange(desc(responseVariable))

# coefficient counts
#heightDiameterCoefficientCounts = heightDiameterCoefficients %>% rowwise() %>% mutate(coefficients = ncol(heightDiameterCoefficients) - sum(is.na(c_across(is.numeric)))) %>%
#  select(responseVariable, species, name, coefficients, coefficients) # slow
#print(heightDiameterCoefficientCounts %>% filter(fitting == "gam"), n = 50)

# summarize fitting success for Section S1
# Show results for all fit sets here, whether or not primary.
heightDiameterResults %>% group_by(fitSet, responseVariable, species, name) %>%
  # reduce cross validated fits to a single record
  summarize(nlrob = fitting[1] == "nlrob", 
            gslNls = fitting[1] == "gsl_nls", 
            gnls = (fitting[1] == "gnls") | (is.na(fitting[1]) & str_detect(name[1], "GNLS")),
            lm = fitting[1] == "lm", 
            gam = fitting[1] == "gam",
            fail = is.na(fitting[1]),
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
  mutate(nlrobPct = 100 * nlrob / uniqueNonlinear, gslNlsPct = nlrobPct + 100 * gslNls / uniqueNonlinear, totalFailPct = 100 * fail / uniqueNonlinear) %>%
  arrange(desc(fitSet))
# fittings which failed
heightDiameterResults %>% filter(is.na(nse)) %>% select(fitSet, responseVariable, species, name) %>% arrange(desc(fitSet))


## Figure 1: overall dataset summary
liveUnbrokenTrees2016 = trees2016 %>% filter(isLiveUnbroken)
plot_exploratory(liveUnbrokenTrees2016 %>% filter(isConifer), speciesLabel = "conifer", maxTreesMeasured = 170, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(isConifer == FALSE), speciesLabel = "broadleaf", maxTreesMeasured = 170, plotLetters = c("d)", "e)", "f)")) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.png", height = 13, width = 22, units = "cm", dpi = 150)
#plot_exploratory(liveUnbrokenTrees2016) + plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure 01 height-diameter distribution.png", height = 10, width = 22, units = "cm", dpi = 150)


## Figure 2: height-diameter error summary
# If axis limits, vertical dodges, or aesthetics are changed Figure 3 should be updated.
heightFromDiameterAccuracyLevels = primaryResults %>% 
  filter(responseVariable == "height") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmse = mean(if_else(is.na(rmse), 100, rmse)),
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
          rmse = quantile(rmse, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE), 
          sizeShapeAlpha = sizeShapeAlpha[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("mapb", "mape", "rmse", "deltaAicN", "nse")) %>%
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
  geom_errorbarh(aes(xmin = rmse_q0.025, xmax = rmse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = rmse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 11)) +
  labs(x = "RMSE, m", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  #geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM", yend = "REML GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM BA+L", yend = "REML GAM BA+L", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  #geom_segment(x = 0.255, xend = 0.275, y = "REML GAM physio", yend = "REML GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = deltaAicN_q0.025, xmax = deltaAicN_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = deltaAicN_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  labs(x = bquote("normalized "*Delta*"AIC"), y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  coord_cartesian(xlim = c(0, 0.3)) + # exclude high AIC of linear, parabolic, and Douglas-fir power+Curtis fits to avoid squashing of nonlinear differences
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_errorbarh(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
#ggplot(heightFromDiameterResults) +
#  geom_point(aes(x = nse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = heightFromDiameterResults$verticalDodge)) +
#  coord_cartesian(xlim = c(0, 1)) +
#  labs(x = "Pearson's R", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
#  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 2), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure 02 height accuracy.png", height = 12, width = 20, units = "cm", dpi = 150)

heightFromDiameterModelComparison = heightDiameterModelRanking %>% filter(responseVariable == "height") %>%
  mutate(name = factor(name, levels = heightFromDiameterAccuracyLevels$name))
ggplot(heightFromDiameterModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucMab)) +
  labs(title = "a) MAB", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(limits = rev) +
ggplot(heightFromDiameterModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucMae)) +
  labs(title = "b) MAE", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(heightFromDiameterModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucRmse)) +
  labs(title = "c) RMSE", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(heightFromDiameterModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucDeltaAicN)) +
  labs(title = "d) AIC", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(heightFromDiameterModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucNse)) +
  labs(title = "e) model efficiency", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin =  margin())) +
plot_layout(nrow = 1, guides = "collect") &
  scale_fill_scico(palette = "bam", limits = c(0, 1), na.value = rgb(0.9642, 0.9444, 0.9435)) &
  scale_x_discrete(labels = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), limits = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")) &
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.spacing.y = unit(0.3, "line"), title = element_text(size = 8))
#ggsave("trees/height-diameter/figures/Figure 02 height accuracy AUC median.png", height = 12, width = 20, units = "cm", dpi = 150)


## Figure 3: diameter-height error summary
# If axis limits, vertical dodges, or aesthetics are changed Figure 2 should be updated.
diameterFromHeightAccuracyLevels = primaryResults %>% 
  filter(responseVariable == "DBH") %>%
  group_by(name, species) %>%
  summarize(speciesFraction = speciesFraction[1],
            meanPenalizedMae = mean(if_else(is.na(mae), 100, mae)),
            meanPenalizedMape = mean(if_else(is.na(mape), 100, mape)),
            meanPenalizedRmse = mean(if_else(is.na(rmse), 100, rmse)),
            .groups = "drop_last") %>% # change to grouping only by name
  summarize(weightedRmse = sum(speciesFraction * meanPenalizedRmse), 
                               .groups = "drop") %>%
  arrange(weightedRmse)
diameterFromHeightResults = primaryResults %>% 
  filter(responseVariable == "DBH", str_detect(name, "BA\\+L") == FALSE, str_detect(name, "GNLS") == FALSE, name != "REML GAM ABA+T physio") %>%
  mutate(name = factor(name, levels = diameterFromHeightAccuracyLevels$name)) %>%
  group_by(species, name) %>%
  reframe(quantiles = c(0.025, 0.5, 0.975),
          mapb = quantile(mapb, probs = quantiles, na.rm = TRUE), 
          mape = quantile(mape, probs = quantiles, na.rm = TRUE), 
          rmse = quantile(rmse, probs = quantiles, na.rm = TRUE), 
          deltaAicN = quantile(deltaAicN, probs = quantiles, na.rm = TRUE), 
          nse = quantile(nse, probs = quantiles, na.rm = TRUE), 
          sizeShapeAlpha = sizeShapeAlpha[1]) %>%
  pivot_wider(names_from = quantiles, names_sep = "_q", values_from = c("mapb", "mape", "rmse", "deltaAicN", "nse")) %>%
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
  geom_errorbarh(aes(xmin = rmse_q0.025, xmax = rmse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = rmse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "RMSE, cm", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
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
  coord_cartesian(xlim = c(0, 0.3)) +
  labs(x = bquote("normalized "*Delta*"AIC"), y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  #scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), trans = scales::pseudo_log_trans()) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_errorbarh(aes(xmin = nse_q0.025, xmax = nse_q0.975, y = name, color = species, alpha = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge), height = 0.1, linewidth = 0.5) +
  geom_point(aes(x = nse_q0.5, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
#ggplot(diameterFromHeightResults) +
#  geom_point(aes(x = pearson, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE, position = position_nudge(y = diameterFromHeightResults$verticalDodge)) +
#  coord_cartesian(xlim = c(0, 1)) +
#  labs(x = "Pearson's R", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
#  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 2), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = speciesGroupColors) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure 03 diameter accuracy.png", height = 12, width = 20, units = "cm", dpi = 150)

diameterFromHeightModelComparison = heightDiameterModelRanking %>% filter(responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = diameterFromHeightAccuracyLevels$name))
ggplot(diameterFromHeightModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucMab)) +
  labs(title = "a) MAB", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(limits = rev) +
ggplot(diameterFromHeightModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucMae)) +
  labs(title = "b) MAE", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(diameterFromHeightModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucRmse)) +
  labs(title = "c) RMSE", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(diameterFromHeightModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucDeltaAicN)) +
  labs(title = "d) AIC", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
ggplot(diameterFromHeightModelComparison) +
  geom_raster(aes(x = species, y = name, fill = aucNse)) +
  labs(title = "e) model efficiency", x = NULL, y = NULL, fill = "median\nAUC") +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin =  margin())) +
plot_layout(nrow = 1, guides = "collect") &
  scale_fill_scico(palette = "bam", limits = c(0, 1), na.value = rgb(0.9642, 0.9444, 0.9435)) &
  scale_x_discrete(labels = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), limits = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")) &
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.spacing.y = unit(0.3, "line"), title = element_text(size = 8))
#ggsave("trees/height-diameter/figures/Figure 03 diameter accuracy AUC median.png", height = 12, width = 20, units = "cm", dpi = 150)


## Figure 4: comparison of accuracy metrics
accuracyCorrelation = bind_rows(as.data.frame(cor(primaryResults %>% filter(responseVariable == "height") %>% 
                                                    select(mapb, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "height"),
                                as.data.frame(cor(primaryResults %>% filter(responseVariable == "DBH") %>% 
                                                    select(mapb, mape, rmse, deltaAicN, nse), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "DBH")) %>%
  mutate(metricX = factor(metricX, levels = c("mapb", "mape", "rmse", "deltaAicN", "nse", "pearson"), labels = c("MAB,\n%", "MAE,\n%", "RMSE", "normalized\nΔAIC", "model\nefficiency", "Pearson's\nR")),
         metricY = factor(metricY, levels = c("mapb", "mape", "rmse", "deltaAicN", "nse", "pearson"), labels = c("MAB, %", "MAE, %", "RMSE", "normalize ΔAIC", "model efficiency", "Pearson's R")))

ggplot(accuracyCorrelation %>% filter(responseVariable == "height")) + 
  coord_equal() +
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  labs(x = NULL, y = NULL, fill = "correlation", title = "a) height prediction") +
ggplot(accuracyCorrelation %>% filter(responseVariable == "DBH")) + 
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "correlation", title = "b) DBH prediction") +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  scale_fill_scico(palette = "vik", limits = c(-1, 1)) &
  theme(legend.spacing.y = unit(0.5, "line"), title = element_text(size = 8))
ggsave("trees/height-diameter/figures/Figure 04 accuracy metric correlation.png", height = 9, width = 20, units = "cm")


## Figure 5: Douglas-fir and red alder model comparison
# preferred forms from 10x10 cross validation
#             height                                           DBH
#             base              generalized                    base                    generalized
# Douglas-fir REML GAM          Sharma-Parton BA+L physio      Michaelis-Menten form   REML GAM ABA+T
#             Sibbesen          Sharma-Parton physio           Ruark                   REML GAM physio
#             Prodan            Sharma-Parton                  Chapman-Richards form   
# red alder   Weibull           Sharma-Parton physio           REML GAM                Chapman-Richards physio
#             Chapman-Richards  Sharma-Parton BA+L physio      parabolic               Sibbesen form physio
#             Michaelis-Menten  Chapman-Richards BA+L physio
#print(preferredForms %>% filter(species %in% c("Douglas-fir", "red alder")) %>% select(-mabName, -aucMab, -nseName, -aucNse, -pearsonName, -aucPearson) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other"), levels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species")), maeName = str_trunc(maeName, 20, ellipsis = ""), rmseName = str_trunc(rmseName, 20, ellipsis = ""), aicName = str_trunc(aicName, 20, ellipsis = "")), n = 32)
if (exists("psmeHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/PSME preferred models.Rdata") }
if (exists("alruHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/ALRU2 preferred models.Rdata") }

ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysio), color = "Sharma-Parton BA+L physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$gam), color = "REML GAM", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$sibbesen), color = "Sibbesen", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) Douglas-fir height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Sibbesen", "Sharma-Parton BA+L physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey60")) +
  theme(legend.key = element_rect(fill = alpha("white", 0.5)), legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gamAat), y = psme2016$TotalHt, color = "REML GAM ABA+T", group = psme2016$isPlantation, linetype = psme2016$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$michaelisMentenForm), y = psme2016$TotalHt, color = "Michaelis-Menten form", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$ruark), y = psme2016$TotalHt, color = "Ruark", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) Douglas-fir DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("Michaelis-Menten form", "Ruark", "REML GAM ABA+T", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alru2016physio$DBH, y = predict(alruHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$chapmanRichards), color = "Chapman-Richards", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$weibull), color = "Weibull", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) red alder height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Weibull", "Sharma-Parton physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$chapmanRichardsPhysio), y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$gam), y = alru2016$TotalHt, color = "REML GAM", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$parabolic), y = alru2016$TotalHt, color = "parabolic", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) red alder DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "parabolic", "Chapman-Richards physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 05 PSME-ALRU2 curves.png", height = 12, width = 20, units = "cm")


## Figure 6: western hemlock and bigleaf maple model comparison
# preferred forms from 10x10 cross validation
#                  height                                           DBH
#                  base              generalized                    base                    generalized
# western hemlock  REML GAM          Sharma-Parton physio           REML GAM                REML GAM ABA+T
#                  Hossfeld IV       Sharma-Parton BA+L physio      Sibbesen form           Chapman-Richards RelHt
#                  Prodan
# bigleaf maple    unified Richards  Weibull BA+L                   REML GAM                power physio
#                  Michaelis-Menten  Chapman-Richards physio        Ruark                   Chapman-Richards physio
#                  Prodan                                           power
#print(preferredForms %>% filter(species %in% c("western hemlock", "bigleaf maple")) %>% select(-mabName, -aucMab, -nseName, -aucNse, -pearsonName, -aucPearson) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other"), levels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species")), maeName = str_trunc(maeName, 20, ellipsis = ""), rmseName = str_trunc(rmseName, 20, ellipsis = ""), aicName = str_trunc(aicName, 20, ellipsis = "")), n = 32)
if (exists("tsheHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/TSHE preferred models.Rdata") }
if (exists("acmaHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/ACMA3 preferred models.Rdata") }

ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = tshe2016physio$DBH, y = predict(tsheHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = tshe2016physio$isPlantation, linetype = tshe2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterPreferred$gam), color = "REML GAM", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = predict(tsheHeightFromDiameterPreferred$hossfeld), color = "Hossfeld IV", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) western hemlock height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Hossfeld IV", "Sharma-Parton physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$gamAat), y = tshe2016$TotalHt, color = "REML GAM ABA+T", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$gam), y = tshe2016$TotalHt, color = "REML GAM", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = predict(tsheDiameterFromHeightPreferred$sibbesenForm), y = tshe2016$TotalHt, color = "Sibbesen form", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) western hemlock DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Sibbesen form", "REML GAM ABA+T", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = acma2016physio$DBH, y = predict(acmaHeightFromDiameterPreferred$weibullBal), color = "Weibull BA+L", group = acma2016physio$isPlantation, linetype = acma2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$michaelisMenten), color = "Michaelis-Menten", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = predict(acmaHeightFromDiameterPreferred$richards), color = "unified Richards", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) bigleaf maple height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Michaelis-Menten", "unified Richards", "Weibull BA+L", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.90)) +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$powerPhysio), y = acma2016physio$TotalHt, color = "power physio", group = acma2016physio$isPlantation, linetype = acma2016physio$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$gam), y = acma2016$TotalHt, color = "REML GAM", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = predict(acmaDiameterFromHeightPreferred$ruark), y = acma2016$TotalHt, color = "Ruark", group = acma2016$isPlantation, linetype = acma2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) bigleaf maple diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Ruark", "power physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.90)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 06 TSHE-ACMA3 curves.png", height = 12, width = 20, units = "cm")


## Figure 7: Oregon myrtle and western redcedar model comparison
# preferred forms from 10x10 cross validation
#                   height                                           DBH
#                   base              generalized                    base                    generalized
# Oregon myrtle     Prodan            Chapman-Richards physio        REML GAM                Chapman-Richards physio
#                   Sibbesen          Chapman-Richards BA+L physio   power                   Sibbesen form Physio
#                   REML GAM                                         parabolic
# western redcedar  Hossfeld IV       Sharma-Parton physio           parabolic               power ABA+T
#                   Michaelis-Menten  Sharma-Parton BA+L physio      Ruark                   Chapman-Richards RelHt
#                   Prodan            Sharma-Parton BA+L             REML GAM
#print(preferredForms %>% filter(species %in% c("Oregon myrtle", "western redcedar")) %>% select(-mabName, -aucMab, -nseName, -aucNse, -pearsonName, -aucPearson) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other"), levels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species")), maeName = str_trunc(maeName, 22, ellipsis = ""), rmseName = str_trunc(rmseName, 22, ellipsis = ""), aicName = str_trunc(aicName, 22, ellipsis = "")), n = 32)
if (exists("umcaHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/UMCA preferred models.Rdata") }
if (exists("thplHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/THPL preferred models.Rdata") }

ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = umca2016physio$DBH, y = predict(umcaHeightFromDiameterPreferred$chapmanRichardsPhysio), color = "Chapman-Richards physio", group = umca2016physio$isPlantation, linetype = umca2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameterPreferred$prodan), color = "Prodan", group = umca2016$isPlantation, linetype = umca2016physio$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = predict(umcaHeightFromDiameterPreferred$sibbesen), color = "Sibbesen", group = umca2016$isPlantation, linetype = umca2016physio$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) Oregon myrtle height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Prodan", "Sibbesen", "Chapman-Richards physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.90)) +
ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$chapmanRichardsPhysio), y = umca2016physio$TotalHt, color = "Chapman-Richards physio", group = umca2016physio$isPlantation, linetype = umca2016physio$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$gam), y = umca2016$TotalHt, color = "REML GAM", group = umca2016$isPlantation, linetype = umca2016physio$isPlantation)) +
  geom_line(aes(x = predict(umcaDiameterFromHeightPreferred$power), y = umca2016$TotalHt, color = "power", group = umca2016$isPlantation, linetype = umca2016physio$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) Oregon myrtle diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("REML GAM", "power", "Chapman-Richards physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.90)) +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = thpl2016physio$DBH, y = predict(thplHeightFromDiameterPreferred$sharmaPartonPhysio), color = "Sharma-Parton physio", group = thpl2016physio$isPlantation, linetype = thpl2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameterPreferred$hossfeld), color = "Hossfeld IV", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation), alpha = 0.8) + # Hossfeld and Michaelis-Menten fits are visually indistinguishable
  geom_line(aes(x = thpl2016$DBH, y = predict(thplHeightFromDiameterPreferred$michaelisMenten), color = "Michaelis-Menten", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation), alpha = 0.8) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) western redcedar height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Hossfeld IV", "Michaelis-Menten", "Sharma-Parton physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$powerAat), y = thpl2016physio$TotalHt, color = "power ABA+T", group = thpl2016physio$isPlantation, linetype = thpl2016physio$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$parabolic), y = thpl2016$TotalHt, color = "parabolic", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation)) +
  geom_line(aes(x = predict(thplDiameterFromHeightPreferred$ruark), y = thpl2016$TotalHt, color = "Ruark", group = thpl2016$isPlantation, linetype = thpl2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) western redcedar diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  scale_color_manual(breaks = c("parabolic", "Ruark", "power ABA+T", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 07 UMCA-THPL curves.png", height = 12, width = 20, units = "cm")


## Figure 8: minority species model comparison
# preferred forms from 10x10 cross validation
#                height                                           DBH
#                base              generalized                    base                    generalized
# other species  Curtis            REML GAM BA+L                  REML GAM                Sibbesen form physio
#                power             REML GAM BAL+L physio          parabolic               Chapman-Richards form RelHt
#                Korf                                             linear
#print(preferredForms %>% filter(species == "other species") %>% select(-mabName, -aucMab, -nseName, -aucNse, -pearsonName, -aucPearson) %>% rename(respVar = responseVariable, base = isBaseForm, aucAic = aucDeltaAicN) %>% mutate(species = factor(species, labels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other"), levels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species")), maeName = str_trunc(maeName, 22, ellipsis = ""), rmseName = str_trunc(rmseName, 22, ellipsis = ""), aicName = str_trunc(aicName, 22, ellipsis = "")), n = 16)
if (exists("otherHeightFromDiameterPreferred") == FALSE) { load("trees/height-diameter/data/other preferred models.Rdata") }

ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$gamBal), color = "REML GAM BA+L", group = other2016$isPlantation, linetype = other2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$curtis), color = "Curtis", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = predict(otherHeightFromDiameterPreferred$power), color = "power", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) other species height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Curtis", "power", "REML GAM BA+L", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$sibbesenFormPhysio), y = other2016physio$TotalHt, color = "Sibbesen form physio", group = other2016physio$isPlantation, linetype = other2016physio$isPlantation), alpha = 0.5, orientation = "y") +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$gam), y = other2016$TotalHt, color = "REML GAM", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = predict(otherDiameterFromHeightPreferred$parabolic), y = other2016$TotalHt, color = "parabolic", group = other2016$isPlantation, linetype = other2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) other species diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  scale_color_manual(breaks = c("REML GAM", "parabolic", "Sibbesen form physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n33", heights = c(1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/Figure 08 other species curves.png", height = 17/3 + 0.5, width = 20, units = "cm")


## Figure 9: model efficiency
# Long tail of negative model efficiencies is suppressed as histogram plotting becomes very slow even if bins
# are not within the display window set by coord_cartesian().
# liveUnbrokenTrees2016 is from Figure 1
smallTreeEfficiency = left_join(liveUnbrokenTrees2016 %>% filter(is.na(TotalHt) == FALSE) %>% mutate(species = factor(speciesGroup, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other"))) %>%
                                  group_by(species) %>%
                                  summarize(stems = n(), smallStems = sum(DBH <= 20.0), pctSmall = 100 * smallStems / stems),
                                primaryResults,
                                by = c("species"))

ggplot(primaryResults %>% filter(responseVariable == "height")) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  annotate("text", x = 0, y = 18.5, label = "a) height prediction", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 19)) +
  labs(x = "model efficiency", y = "fraction of models, %", fill = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(responseVariable == "DBH")) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  annotate("text", x = 0, y = 18.5, label = "b) DBH prediction", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 19)) +
  labs(x = "model efficiency", y = NULL, fill = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(smallTreeEfficiency) + 
  #geom_violin(aes(x = pctSmall, y = nse, color = species, group = species), draw_quantiles = c(0.25, 0.50, 0.75), fill = alpha("white", 0.5), position = position_identity(), width = 2) +
  geom_boxplot(aes(x = pctSmall, y = nse, color = species, group = species), fill = alpha("white", 0.5), position = position_identity(), na.rm = TRUE, width = 4) +
  annotate("text", x = 0, y = 21.5/20, label = "c) lack of small stem effect", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 22/20)) +
  labs(x = "stems less than 20 cm DBH, %", y = "model efficiency", color = NULL) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 3, guides = "collect") &
  guides(color = "none", fill = guide_legend(ncol = 7)) &
  scale_color_manual(breaks = levels(primaryResults$species), limits = levels(primaryResults$species), values = speciesGroupColors) &
  scale_fill_manual(breaks = levels(primaryResults$species), limits = levels(primaryResults$species), values = speciesGroupColors) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.position = "bottom")
#ggsave("trees/height-diameter/figures/Figure 09 model efficiency.png", height = 7, width = 20, units = "cm")


## Figures S1-3: species level exploratory plots
# liveUnbrokenTrees2016 is from Figure 1
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "DF"), speciesLabel = "Douglas-fir", maxTreesMeasured = 150, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "RA"), speciesLabel = "red alder", maxTreesMeasured = 150, plotLetters = c("d)", "e)", "f)"), omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "WH"), speciesLabel = "western hemlock", maxTreesMeasured = 150, plotLetters = c("g)", "h)", "i)"), omitLegends = TRUE) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure S01 PSME-ALRU2-TSHE.png", height = 18, width = 22, units = "cm", dpi = 150)

plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "BM"), speciesLabel = "bigleaf maple", maxTreesMeasured = 150, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "OM"), speciesLabel = "Oregon myrtle", maxTreesMeasured = 150, distributionLegendPositionY = 0.92, plotLetters = c("d)", "e)", "f)"), omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "RC"), speciesLabel = "western redcedar", maxTreesMeasured = 150, plotLetters = c("g)", "h)", "i)"), omitLegends = TRUE) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure S02 ACMA3-UMCA-THPL.png", height = 18, width = 22, units = "cm", dpi = 150)

plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "other"), speciesLabel = "other species ", distributionLegendPositionY = 0.92) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/figures/Figure S03 other species.png", height = 1/3*(18 - 1) + 1, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 1 Douglas-fir.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 2 red alder.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 3 western hemlock.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 4 bigleaf maple.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 5 Oregon myrtle.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 6 western redcedar.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/figures/exploratory 7 minority species.png", height = 10, width = 22, units = "cm", dpi = 150)


## Figures S4-6: estimates of residual variance
# heightIqr, diameterIqr, and residualPower are calculated in residuals.R
ggplot() +
  geom_line(aes(x = seq(0, 250), y = seq(0, 250)^0.5), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = DBH, y = iqr, color = species, group = paste(species, name, isPlantation)), heightIqr %>% filter(name == "Michaelis-Menten", isPlantation == FALSE)) +
  annotate("text", x = 0, y = 24, label = "a) natural regeneration", hjust = 0, size = 3.5) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "DBH, cm", y = "Michaelis-Menten interquartile height range, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = speciesGroupColors) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_line(aes(x = seq(0, 250), y = seq(0, 250)^0.5), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = DBH, y = iqr, color = species, group = paste(species, name, isPlantation)), heightIqr %>% filter(name == "Michaelis-Menten", isPlantation == TRUE)) +
  annotate("text", x = 0, y = 24, label = "b) plantations", hjust = 0, size = 3.5) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = speciesGroupColors) +
  theme(legend.justification = c(1, 0), legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
ggsave("trees/height-diameter/figures/Figure S04 height interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

ggplot() +
  geom_line(aes(x = seq(0, 100), y = 80/80*seq(0, 100)^0.9), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = TotalHt, y = iqr, color = species, group = paste(species, name, isPlantation)), diameterIqr %>% filter(isPlantation == FALSE)) +
  annotate("text", x = 0, y = 62, label = "a) natural regeneration", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 82), ylim = c(0, 61)) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "height, m", y = "Ruark interquartile DBH range, cm", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = speciesGroupColors) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_line(aes(x = TotalHt, y = iqr, color = species, group = paste(species, name, isPlantation)), diameterIqr %>% filter(isPlantation == TRUE)) +
  geom_line(aes(x = seq(0, 100), y = 30/80*seq(0, 100)^1.2), color = "grey70", linetype = "longdash") +
  annotate("text", x = 0, y = 62, label = "b) plantations", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 82), ylim = c(0, 61)) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "height, m", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = speciesGroupColors) +
  theme(legend.justification = c(1, 0), legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
ggsave("trees/height-diameter/figures/Figure S05 diameter interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

ggplot(residualPower) +
  geom_errorbar(aes(xmin = ht_b1_min, xmax = ht_b1_max, y = htName, color = species), width = 0.25) +
  geom_errorbar(aes(xmin = ht_b1p_min, xmax = ht_b1p_max, y = htName, color = species), width = 0.25, na.rm = TRUE) +
  geom_point(aes(x = ht_b1, y = htName, color = species, shape = "b1"), size = 1.7) +
  geom_point(aes(x = ht_b1 + ht_b1p, y = htName, color = species, shape = "b1p"), size = 1.7, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 2.1)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  labs(x = "power fitted to height residuals", y = NULL, color = NULL, shape = NULL) +
  scale_y_discrete(limits = rev) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
ggplot(residualPower) +
  geom_segment(x = 2.01, xend = 2.11, y = "Ruark", yend = "Ruark", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", data = tibble(species = factor("other species", levels = levels(residualPower$species))), linewidth = 0.4) +
  geom_errorbar(aes(xmin = dia_b1_min, xmax = dia_b1_max, y = diaName, color = species), width = 0.25) +
  geom_errorbar(aes(xmin = dia_b1p_min, xmax = dia_b1p_max, y = diaName, color = species), width = 0.25) +
  geom_point(aes(x = dia_b1, y = diaName, color = species, shape = "b1"), size = 1.7) +
  geom_point(aes(x = dia_b1 + dia_b1p, y = diaName, color = species, shape = "b1p"), size = 1.7) +
  coord_cartesian(xlim = c(0, 2.1)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  labs(x = "power fitted to DBH residuals", y = NULL, color = NULL, shape = NULL) +
  scale_y_discrete(limits = rev) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_blank()) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, guides = "collect") &
  guides(color = "none") &
  scale_color_manual(breaks = levels(residualPower$species), limits = levels(residualPower$species), values = speciesGroupColors) &
  scale_shape_manual(breaks = c("b1", "b1p"), labels = c("base power", "plantations"), values = c(15, 17)) &
  theme(legend.position = "bottom")
ggsave("trees/height-diameter/figures/Figure S06 residual power estimates.png", height = 14, width = 20, units = "cm", dpi = 150)


## Figures S7-12: Q-Q plots
plot_qq(psmeHeightFromDiameterPreferred$chapmanRichards, psmeHeightFromDiameterPreferred$michaelisMenten, psmeHeightFromDiameterPreferred$sharmaParton, psmeHeightFromDiameterPreferred$sharmaZhang,
        psmeDiameterFromHeightPreferred$ChapmanRichards, psmeDiameterFromHeightPreferred$chapmanForm, psmeDiameterFromHeightPreferred$ruark, psmeDiameterFromHeightPreferred$sibbesenForm,
        "Douglas-fir")
ggsave("trees/height-diameter/figures/Figure S07 PSME Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(alruHeightFromDiameterPreferred$chapmanRichards, alruHeightFromDiameterPreferred$michaelisMenten, alruHeightFromDiameterPreferred$sharmaParton, alruHeightFromDiameterPreferred$sharmaZhang,
        alruDiameterFromHeightPreferred$chapmanRichards, alruDiameterFromHeightPreferred$chapmanForm, alruDiameterFromHeightPreferred$ruark, alruDiameterFromHeightPreferred$sibbesenForm,
        "red alder", tDegreesOfFreedom = 8, tSkew = 2.1)
ggsave("trees/height-diameter/figures/Figure S08 ALRU2 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(tsheHeightFromDiameterPreferred$chapmanRichards, tsheHeightFromDiameterPreferred$michaelisMenten, tsheHeightFromDiameterPreferred$sharmaParton, tsheHeightFromDiameterPreferred$sharmaZhang,
        tsheDiameterFromHeightPreferred$chapmanRichards, tsheDiameterFromHeightPreferred$chapmanForm, tsheDiameterFromHeightPreferred$ruark, tsheDiameterFromHeightPreferred$sibbesenForm,
        "western hemlock", tDegreesOfFreedom = 7, tSkew = 2)
ggsave("trees/height-diameter/figures/Figure S08 TSHE Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(acmaHeightFromDiameterPreferred$chapmanRichards, acmaHeightFromDiameterPreferred$michaelisMenten, acmaHeightFromDiameterPreferred$sharmaParton, acmaHeightFromDiameterPreferred$sharmaZhang,
        acmaDiameterFromHeightPreferred$chapmanRichards, acmaDiameterFromHeightPreferred$chapmanForm, acmaDiameterFromHeightPreferred$ruark, acmaDiameterFromHeightPreferred$sibbesenForm,
        "bigleaf maple", tDegreesOfFreedom = 10, tSkew = 4)
ggsave("trees/height-diameter/figures/Figure S09 ACMA3 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(umcaHeightFromDiameterPreferred$chapmanRichards, umcaHeightFromDiameterPreferred$michaelisMenten, umcaHeightFromDiameterPreferred$sharmaParton, umcaHeightFromDiameterPreferred$sharmaZhang,
        umcaDiameterFromHeightPreferred$chapmanRichards, umcaDiameterFromHeightPreferred$chapmanForm, umcaDiameterFromHeightPreferred$ruark, umcaDiameterFromHeightPreferred$sibbesenForm,
        "Oregon myrtle", tDegreesOfFreedom = 8, tSkew = 4)
ggsave("trees/height-diameter/figures/Figure S10 UMCA Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(thplHeightFromDiameterPreferred$chapmanRichards, thplHeightFromDiameterPreferred$michaelisMenten, thplHeightFromDiameterPreferred$sharmaParton, thplHeightFromDiameterPreferred$sharmaZhang,
        thplDiameterFromHeightPreferred$chapmanRichards, thplDiameterFromHeightPreferred$chapmanForm, thplDiameterFromHeightPreferred$ruark, thplDiameterFromHeightPreferred$sibbesenForm,
        "western redcedar", tDegreesOfFreedom = 8, tSkew = 3)
ggsave("trees/height-diameter/figures/Figure S11 THPL Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(otherHeightFromDiameterPreferred$chapmanRichards, otherHeightFromDiameterPreferred$michaelisMenten, otherHeightFromDiameterPreferred$sharmaParton, otherHeightFromDiameterPreferred$sharmaZhang,
        otherDiameterFromHeightPreferred$chapmanRichards, otherDiameterFromHeightPreferred$chapmanForm, otherDiameterFromHeightPreferred$ruark, otherDiameterFromHeightPreferred$sibbesenForm,
        "other species", tDegreesOfFreedom = 3, tSkew = 1)
ggsave("trees/height-diameter/figures/Figure S12 other Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)


## companion variance structure summary from GNLS
# Use all available results here as gnls() is excluded from primary results set by default.
ggplot(heightDiameterResults %>% filter(is.na(power) == FALSE) %>% mutate(name = str_replace(name, " GNLS", ""))) +
  geom_point(aes(x = power, y = name, color = species)) +
  coord_cartesian(xlim = c(0, 2)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  guides(color = guide_legend(ncol = 7)) +
  labs(x = "power fitted to height residuals", y = NULL, color = NULL) +
  scale_color_manual(breaks = levels(heightDiameterResults$species), limits = levels(heightDiameterResults$species), values = speciesGroupColors) +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")

ggplot(heightDiameterResults %>% filter(responseVariable == "height", fitting != "gnls")) +
  geom_boxplot(aes(x = bias, y = species, color = species), width = 0.6) +
  coord_cartesian(xlim = c(-2, 2)) +
  guides(color = "none") +
  labs(x = "nlrob() or gsl_nls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightDiameterResults$species), limits = levels(heightDiameterResults$species), values = speciesGroupColors) +
  theme(legend.position = "none") +
ggplot(heightDiameterResults %>% filter(responseVariable == "height", fitting == "gnls")) +
  geom_boxplot(aes(x = bias, y = species, color = species), width = 0.6) +
  coord_cartesian(xlim = c(-2, 2)) +
  guides(color = "none") +
  labs(x = "gnls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightDiameterResults$species), limits = levels(heightDiameterResults$species), values = speciesGroupColors) +
  scale_y_discrete(labels = NULL) +
  theme(legend.position = "none")
  