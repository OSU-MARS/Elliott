# assumes library()s, functions, and dataset from treetops.R setup
library(rsample) # blocking by merge cluster not needed

get_radius_accuracy_dsm_quadratic = function(data, coefficients)
{
  a0 = coefficients[1]
  a1 = coefficients[2]
  a2 = coefficients[3]
  data %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height + a2 * height^2), levels = c(FALSE, TRUE))) 
  return(caret::confusionMatrix(data$isTreetopRadiusParameterSearch, data$isTreetopRadiusDsm)$overall["Accuracy"])
}

get_radius_accuracy_dsm_power = function(data, coefficients)
{
  a0 = coefficients[1]
  a1 = coefficients[2]
  b1 = coefficients[3]
  data %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height^b1), levels = c(FALSE, TRUE))) 
  return(caret::confusionMatrix(data$isTreetopRadiusParameterSearch, data$isTreetopRadiusDsm)$overall["Accuracy"])
}

rename_vfold_cv_ids = function(splitsAndFits)
{
  if ("id2" %in% names(splitsAndFits))
  {
    splitsAndFits %<>% rename(repetition = id, fold = id2) %>% mutate(repetition = as.integer(str_remove(repetition, "Repeat")),
                                                                      fold = as.integer(str_remove(fold, "Fold")))
  } else { # 
    splitsAndFits %<>% rename(fold = id) %>% mutate(repetition = 1,
                                                    fold = as.integer(str_remove(fold, "Fold"))) %>%
      relocate(repetition, fold)
  }
  return(splitsAndFits)
}

# linear decision boundary
# Optimums found
#             DSM                                                         accuracy
# linear      0.419059 + 0.0381804 h                                      0.929002 Nelder-Mead
# quadratic   0.450256 + 0.0381384 h - 2.835942e-05 h²                    0.929273 Nelder-Mead (BFGS diverges to positive h²)
#             0.420000 + 0.0374000 h - 4.702588e-05 h²                    0.929389 conjugate gradient
# power       0.42184646 + 0.03624721 h^0.99839497                        0.929327 BFGS -> 2x25 cross validation + Kolmogorov-Smirnov selected
# cubic       0.437904 + 0.0468697 h - 5.203931e-04 h² + 4.891224e-06 h³  0.928593 Nelder-Mead but of questionable utility
# TODO: fit CHM once maxima are available?
dsmFitStart = Sys.time()
dsmFitPower = optim(par = c(0.420, 0.038, 1), method = "BFGS", control = list(fnscale = -1, trace = 1), fn = function(coefficients)
  {
    return(get_radius_accuracy_dsm_power(treetopData, coefficients))
  })
dsmFitPower
Sys.time() - dsmFitStart

dsmFitQuadratic = optim(par = c(0.420, 0.0374, -0.00005), get_radius_accuracy_dsm_quadratic, method = "CG", control = list(fnscale = -1, trace = 1)) # negative fnscale maximizes
dsmFitQuadratic

if (treetopOptions$includeSetup)
{
  handlers(global = TRUE)
  handlers("cli")
  plan(multisession, workers = 8)
  
  fit_radius_dsm_power = function(maximaData, startingParameters, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
  {
    progressBar = progressor(steps = folds * repetitions)
    
    fitFunction = function(dataFold)
    {
      #dataFold = splitsAndFits$splits[[1]]
      foldTrainingData = analysis(dataFold)
      get_radius_accuracy_dsm_power_cv = function(coefficients)
      {
        return(get_radius_accuracy_dsm_power(foldTrainingData, coefficients))
      }
      
      fit = optim(par = startingParameters, get_radius_accuracy_dsm_power_cv, method = "BFGS", control = list(fnscale = -1))
      a0 = fit$par[1]
      a1 = fit$par[2]
      b1 = fit$par[3]
      
      validationData = assessment(dataFold)
      treetopPrediction = factor(validationData$radius > (a0 + a1 * validationData$height^b1), levels = c(FALSE, TRUE))

      confusion = caret::confusionMatrix(validationData$isTreetopRadiusDsm, treetopPrediction)
      overallAccuracyByHeightClass = validationData %>% select(isTreetopRadiusDsm, height) %>% 
        mutate(heightClass = round(height),
               treetopPrediction = treetopPrediction) %>%
        group_by(heightClass) %>%
        summarize(n = n(), overallAccuracy = sum(isTreetopRadiusDsm == treetopPrediction) / n(), .groups = "drop")
      
      progressBar()
      return(tibble(fit = list(fit),
                    overallAccuracy = confusion$overall["Accuracy"],
                    confusionMatrix = list(confusion),
                    overallAccuracyByHeight = overallAccuracyByHeightClass))
    }
    
    splitsAndFits = vfold_cv(maximaData, v = folds, repeats = repetitions) %>%
      mutate(fit = future_map(splits, fitFunction)) %>% 
      select(-splits) %>% 
      unnest(fit)
    return(rename_vfold_cv_ids(splitsAndFits))
  }
  
  fit_radius_dsm_quadratic = function(maximaData, startingParameters, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
  {
    progressBar = progressor(steps = folds * repetitions)
    
    fitFunction = function(dataFold)
    {
      #dataFold = splitsAndFits$splits[[1]]
      foldTrainingData = analysis(dataFold)
      get_radius_accuracy_dsm_quadratic_cv = function(coefficients)
      {
        return(get_radius_accuracy_dsm_quadratic(foldTrainingData, coefficients))
      }
      
      fit = optim(par = startingParameters, get_radius_accuracy_dsm_quadratic_cv, method = "BFGS", control = list(fnscale = -1))
      a0 = fit$par[1]
      a1 = fit$par[2]
      a2 = fit$par[3]
      
      validationData = assessment(dataFold)
      treetopPrediction = factor(validationData$radius > (a0 + a1 * validationData$height + a2 * validationData$height^2), levels = c(FALSE, TRUE))
      
      confusion = caret::confusionMatrix(validationData$isTreetopRadiusDsm, treetopPrediction)
      overallAccuracyByHeightClass = validationData %>% select(isTreetopRadiusDsm, height) %>% 
        mutate(heightClass = round(height),
               treetopPrediction = treetopPrediction) %>%
        group_by(heightClass) %>%
        summarize(n = n(), overallAccuracy = sum(isTreetopRadiusDsm == treetopPrediction) / n(), .groups = "drop")
      
      progressBar()
      return(tibble(fit = list(fit),
                    overallAccuracy = confusion$overall["Accuracy"],
                    confusionMatrix = list(confusion),
                    overallAccuracyByHeight = overallAccuracyByHeightClass))
    }
    
    splitsAndFits = vfold_cv(maximaData, v = folds, repeats = repetitions)
      mutate(fit = future_map(splits, fitFunction)) %>% 
      select(-splits) %>% 
      unnest(fit)
    return(rename_vfold_cv_ids(splitsAndFits))
  }

  powerStart = Sys.time() # ~63s, 9900X
  radiusDsmAccuracyPower = fit_radius_dsm_power(treetopData, c(0.420, 0.038, 1))
  Sys.time() - powerStart
  saveRDS(radiusDsmAccuracyPower, "trees/segmentation/treetops/radius DSM power s4268 458k 2x25.Rds")
  
  quadraticStart = Sys.time()
  radiusDsmAccuracyQuadratic = fit_radius_dsm_quadratic(treetopData, c(0.420, 0.0374, -0.00005))
  Sys.time() - quadraticStart
  
  tibble(power = tibble::num(mean(radiusDsmAccuracyPower$overallAccuracy), digits = 5), 
         quad = tibble::num(mean(radiusDsmAccuracyQuadratic$overallAccuracy), digits = 5), 
         kolmoP = ks.test(radiusDsmAccuracyPower$overallAccuracy, radiusDsmAccuracyQuadratic$overallAccuracy)$p.value)
}

if (treetopOptions$includeInvestigatory)
{
  get_radius_accuracy_dsm_linear = function(coefficients)
  {
    a0 = coefficients[1]
    a1 = coefficients[2]
    treetopData %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height), levels = c(FALSE, TRUE))) 
    confusion = caret::confusionMatrix(treetopData$isTreetopRadiusParameterSearch, treetopData$isTreetopRadiusDsm)
    return(confusion$overall["Accuracy"])
  }
  
  get_radius_accuracy_dsm_cubic = function(coefficients)
  {
    a0 = coefficients[1]
    a1 = coefficients[2]
    a2 = coefficients[3]
    a3 = coefficients[4]
    treetopData %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height + a2 * height^2 + a3 * height^3), levels = c(FALSE, TRUE))) 
    confusion = caret::confusionMatrix(treetopData$isTreetopRadiusParameterSearch, treetopData$isTreetopRadiusDsm)
    return(confusion$overall["Accuracy"])
  }
  
  dsmFitLinear = optim(par = c(0.420, 0.362), get_radius_accuracy_dsm_linear, control = list(fnscale = -1))
  dsmFitLinear
  
  dsmFitCubic = optim(par = c(0.420, 0.0374, -0.00054, 0), get_radius_accuracy_dsm_cubic, control = list(fnscale = -1, trace = 1))
  dsmFitCubic
  
  # grid search for optim() verification
  # searchAccuracy = crossing(a0 = seq(0.41, 0.43, length.out = 6), 
  #                           a1 = seq(0.037, 0.039, length.out = 7), 
  #                           a2 = seq(-0.00008, -0.00004, length.out = 7)) %>%
  # searchAccuracy = crossing(a0 = tibble::num(seq(0.41, 0.43, length.out = 11), digits = 4), 
  #                           a1 = tibble::num(seq(0.036, 0.038, length.out = 11), digits = 4),
  #                           a2 = tibble::num(seq(-0.00006, -0.00004, length.out = 21), digits = 7),
  #                           a3 = tibble::num(seq(-0.000005, -0.000005, length.out = 21), digits = 9)) %>%
  #     mutate(overallAccuracy = NA_real_, confusionMatrix = vector(mode = "list", length = n()))
  # 
  # searchStart = Sys.time()
  # for (index in 1:nrow(searchAccuracy)) # ~0.1s per parameterization
  # {
  #   a0 = searchAccuracy$a0[index]
  #   a1 = searchAccuracy$a1[index]
  #   a2 = searchAccuracy$a2[index]
  #   a3 = searchAccuracy$a2[index]
  #   treetopData %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height + a2 * height^2), levels = c(FALSE, TRUE))) 
  #   #confusion = caret::confusionMatrix(treetopData$isTreetopRadiusParameterSearch, treetopData$isTreetopRadiusChm)
  #   confusion = caret::confusionMatrix(treetopData$isTreetopRadiusParameterSearch, treetopData$isTreetopRadiusDsm)
  #   searchAccuracy$overallAccuracy[index] = confusion$overall["Accuracy"]
  #   searchAccuracy$confusionMatrix[index] = list(confusion)
  # }
  # Sys.time() - searchStart
  # 
  # ggplot() +
  #   geom_tile(aes(x = a0, y = a1, fill = overallAccuracy), searchAccuracy, linewidth = 1) +
  #   geom_tile(aes(x = a0, y = a1, color = overallAccuracy == max(overallAccuracy)), searchAccuracy, fill = "transparent", linewidth = 0.5) +
  #   labs(x = "a0", y = "a1", color = NULL, fill = "overall\naccuracy") +
  # ggplot() +
  #   geom_tile(aes(x = a1, y = a2, fill = overallAccuracy), searchAccuracy, linewidth = 1) +
  #   geom_tile(aes(x = a1, y = a2, color = overallAccuracy == max(overallAccuracy)), searchAccuracy, fill = "transparent", linewidth = 0.5) +
  #   labs(x = "a1", y = "a2", color = NULL, fill = "overall\naccuracy") +
  # plot_annotation(theme = theme(plot.margin = margin())) +
  # plot_layout(guides = "collect") &
  #   scale_color_manual(breaks = c(TRUE, FALSE), labels = c("max", ""), values = c("blue", "transparent")) &
  #   scale_fill_viridis_c()
  # 
  # print(searchAccuracy %>% slice_max(overallAccuracy, n = 4))
  # (searchAccuracy %>% slice_max(overallAccuracy, n = 4))$overallAccuracy
}