# get results from results.R
handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 7)


## find AUCs across fit sets and weightings
# fit set      primary    gsl_nls    mixed      nlrob
# weightings   variable   default    variable   variable
with_progress({
  crossFitSetModelCount = primaryResults %>% group_by(responseVariable, species) %>% summarize(n = n_distinct(name), .groups = "drop")
  progressBar = progressor(steps = sum(crossFitSetModelCount$n))
  
  crossFitSetAucs = primaryResults %>%
    group_by(responseVariable, species, name) %>%
    group_split() %>%
    future_map_dfr(function(primaryFitSetResults)
    {
      if ((nrow(primaryFitSetResults) == 1) | all(is.na(primaryFitSetResults$nse)))
      {
        # no distribution to compare to since this model has only a no fit result or wasn't cross validated
        progressBar(str_pad(paste(primaryFitSetResults$responseVariable[1], primaryFitSetResults$species[1], primaryFitSetResults$name[1]), 60, "right"))
        return(tibble(responseVariable = primaryFitSetResults$responseVariable[1], species = primaryFitSetResults$species[1], name = primaryFitSetResults$name[1],
                      otherModelName = NA_character_, fitting = primaryFitSetResults$fitting[1], isBaseForm = primaryFitSetResults$isBaseForm[1], hasPhysio = primaryFitSetResults$hasPhysio[1], hasStand = primaryFitSetResults$hasStand[1],
                      aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucPearson = NA_real_, aucNse = NA_real_, aucRmse = NA_real_,
                      speciesFraction = primaryFitSetResults$speciesFraction[1]))
      }
      
      # get all other cross-validation results for this response variable, species, and model
      # Assumes no names are shared across fittings in the results set.
      matchingFitResults = heightDiameterResults %>% filter(fitSet != "primary", responseVariable == primaryFitSetResults$responseVariable[1], species == primaryFitSetResults$species[1]) %>%
        mutate(deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) %>%
        filter(name == primaryFitSetResults$name[1])
      matchingFitSets = unique(matchingFitResults$fitSet)
      pairwiseAucs = bind_rows(lapply(matchingFitSets, function(otherFitSet) 
      {
        otherFitSetResults = matchingFitResults %>% filter(fitSet == otherFitSet)
        
        if ((nrow(otherFitSetResults) == 1) | all(is.na(otherFitSetResults$nse)))
        {
          # also no distribution to compare to if the other model has only a no fit result or wasn't cross validated
          return(tibble(responseVariable = primaryFitSetResults$responseVariable[1], species = primaryFitSetResults$species[1], name = primaryFitSetResults$name[1], fixedWeight = primaryFitSetResults$fixedWeight[1],
                        otherFitSet = otherFitSet, otherFitting = otherFitSetResults$fitting[1], otherFixedWeight = otherFitSetResults$fixedWeight[1], hasPhysio = primaryFitSetResults$hasPhysio[1], hasStand = primaryFitSetResults$hasStand[1],
                        aucDeltaAicN = NA_real_, aucMab = NA_real_, aucMae = NA_real_, aucPearson = NA_real_, aucNse = NA_real_, aucRmse = NA_real_,
                        speciesFraction = primaryFitSetResults$speciesFraction[1], isBaseForm = primaryFitSetResults$isBaseForm[1]))
        }
        
        # find AUCs: Labels are set so that AUCs > 0.5 indicate models not from from the primary fit set perform
        # better than models in the primary fit set. Thus, labeling is reversed from this code's equivalent in results.R.
        lowerIsBetter = factor(c(rep(1, nrow(primaryFitSetResults)), rep(0, nrow(otherFitSetResults))), levels = c(0, 1))
        higherIsBetter = factor(c(rep(0, nrow(primaryFitSetResults)), rep(1, nrow(otherFitSetResults))), levels = c(0, 1))
        
        # estimate AUC for bias based on what data is available, which is potentially none for species with few measurements
        availableMabData = tibble(guess = c(primaryFitSetResults$mab, otherFitSetResults$mab), label = lowerIsBetter) %>% drop_na()
        aucMab = NA_real_
        if (nrow(availableMabData) > 0)
        {
          aucMab = WeightedAUC(WeightedROC(guess = availableMabData$guess, label = availableMabData$label))
        }
        
        if ((sum(is.na(primaryFitSetResults$mae)) + sum(is.na(otherFitSetResults$mae)) + sum(is.na(primaryFitSetResults$pearson)) + sum(is.na(otherFitSetResults$pearson)) + 
             sum(is.na(primaryFitSetResults$nse)) + sum(is.na(otherFitSetResults$nse)) + sum(is.na(primaryFitSetResults$rmse)) + sum(is.na(otherFitSetResults$rmse))) > 0)
        {
          # WeightedROC() errors on NAs in its arguments but can't do so informatively
          # Fail informatively instead so that investigation is possible.
          print(otherFitSetResults %>% filter(is.na(mae)))
          stop(paste0(primaryFitSetResults$species[1], " ", primaryFitSetResults$responseVariable[1], " ", primaryFitSetResults$name[1], " (", nrow(primaryFitSetResults), ") x ", otherFitSet, " (", nrow(otherFitSetResults), "):",
                      " mab ", sum(is.na(primaryFitSetResults$mab)), " ", sum(is.na(otherFitSetResults$mab)),
                      " mae ", sum(is.na(primaryFitSetResults$mae)), " ", sum(is.na(otherFitSetResults$mae)),
                      " pearson ", sum(is.na(primaryFitSetResults$pearson)), " ", sum(is.na(otherFitSetResults$pearson)),
                      " nse ", sum(is.na(primaryFitSetResults$nse)), " ", sum(is.na(otherFitSetResults$nse)),
                      " rmse ", sum(is.na(primaryFitSetResults$rmse)), " ", sum(is.na(otherFitSetResults$rmse))), quote = FALSE)
        }
        
        deltaAicNguess = c(primaryFitSetResults$deltaAicN, otherFitSetResults$deltaAicN)
        maeGuess = c(primaryFitSetResults$mae, otherFitSetResults$mae)
        pearsonGuess = c(primaryFitSetResults$pearson, otherFitSetResults$pearson)
        nseGuess = c(primaryFitSetResults$nse, otherFitSetResults$nse)
        rmseGuess = c(primaryFitSetResults$rmse, otherFitSetResults$rmse)
        expectedLength = nrow(primaryFitSetResults) + nrow(otherFitSetResults)
        if ((length(deltaAicNguess) != expectedLength) | (length(maeGuess) != expectedLength) | (length(pearsonGuess) != expectedLength) |
            (length(nseGuess) != expectedLength) | (length(rmseGuess) != expectedLength))
        {
          stop(paste0(primaryFitSetResults$species[1], " ", primaryFitSetResults$responseVariable[1], " ", primaryFitSetResults$name[1], " (", nrow(primaryFitSetResults), ") x ", otherFitSet, " (", nrow(otherFitSetResults), ") NULLs in data:", 
                      " AICn ", length(deltaAicNguess), 
                      " MAE ", length(maeGuess),
                      " Pearson ", length(pearsonGuess),
                      " NSE ", length(nseGuess),
                      " RMSE ", length(rmseGuess)))
        }
        return(tibble(responseVariable = primaryFitSetResults$responseVariable[1], species = primaryFitSetResults$species[1], 
                      name = primaryFitSetResults$name[1], fitting = primaryFitSetResults$fitting[1], fixedWeight = primaryFitSetResults$fixedWeight[1],
                      otherFitSet = otherFitSet, otherFitting = otherFitSetResults$fitting[1], otherFixedWeight = otherFitSetResults$fixedWeight[1],
                      hasPhysio = primaryFitSetResults$hasPhysio[1], hasStand = primaryFitSetResults$hasStand[1],
                      aucDeltaAicN = WeightedAUC(WeightedROC(guess = deltaAicNguess, label = lowerIsBetter)),
                      aucMab = aucMab,
                      aucMabN = nrow(availableMabData),
                      aucMae = WeightedAUC(WeightedROC(guess = maeGuess, label = lowerIsBetter)),
                      aucPearson = WeightedAUC(WeightedROC(guess = pearsonGuess, label = higherIsBetter)),
                      aucNse = WeightedAUC(WeightedROC(guess = nseGuess, label = higherIsBetter)),
                      aucRmse = WeightedAUC(WeightedROC(guess = rmseGuess, label = lowerIsBetter)),
                      speciesFraction = primaryFitSetResults$speciesFraction[1], isBaseForm = primaryFitSetResults$isBaseForm[1]))
      }))
      
      progressBar(str_pad(paste(primaryFitSetResults$responseVariable[1], primaryFitSetResults$species[1], primaryFitSetResults$name[1]), 60, "right"))
      return(pairwiseAucs)
    })
})
#heightDiameterResults %>% filter(species == "Douglas-fir", name == "Weibull BA+L") %>%
#  group_by(fitSet) %>%
#  reframe(quantile = c(0.1, 0.5, 0.9), 
#          mapb = quantile(mapb, probs = c(0.1, 0.5, 0.9)), 
#          mape = quantile(mape, probs = c(0.1, 0.5, 0.9)), 
#          rmse = quantile(rmse, probs = c(0.1, 0.5, 0.9)), 
#          aic = quantile(aic, probs = c(0.1, 0.5, 0.9)), 
#          nse = quantile(nse, probs = c(0.1, 0.5, 0.9)))


## primary fits to gsl_nls(weights = default) comparison
# heightDiameterModelDisplaySort from results.R
gslNlsHeightAucsDefaultWeight = crossFitSetAucs %>% filter(otherFitSet == "gsl_nls", responseVariable == "height") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "height"))$name)))
plot_auc_bank(gslNlsHeightAucsDefaultWeight, fillLabel = "pairwise\nAUC")

gslNlsDbhAucsDefaultWeight = crossFitSetAucs %>% filter(otherFitSet == "gsl_nls", responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "DBH"))$name)))
plot_auc_bank(gslNlsDbhAucsDefaultWeight, fillLabel = "pairwise\nAUC")


## primary fits to nlrob() comparison
nlrobHeightAucs = crossFitSetAucs %>% filter(otherFitSet == "nlrob", responseVariable == "height") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "height"))$name)))
plot_auc_bank(nlrobHeightAucs, fillLabel = "pairwise\nAUC")

nlrobDbhAucs = crossFitSetAucs %>% filter(otherFitSet == "nlrob", responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "DBH"))$name)))
plot_auc_bank(nlrobDbhAucs, fillLabel = "pairwise\nAUC")


## primary fits to mixed effects comparison
mixedHeightAucs = crossFitSetAucs %>% filter(otherFitSet == "mixed", responseVariable == "height") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "height"))$name)))
plot_auc_bank(mixedHeightAucs, fillLabel = "pairwise\nAUC")

mixedDbhAucs = crossFitSetAucs %>% filter(otherFitSet == "mixed", responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "DBH"))$name)))
plot_auc_bank(mixedDbhAucs, fillLabel = "pairwise\nAUC")
