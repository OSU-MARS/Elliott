#.libPaths(.libPaths()[2])
#install.packages(c("dplyr", "forcats", "furrr", "ggplot2", "gslnls", "magrittr", "mgcv", "nlstools", "patchwork", "progressr", "ragg", "readxl", "remotes", "robustbase", "rsample", "scales", "scico", "sn", "stringr", "tibble", "tidyr", "WeightedROC", "writexl"))
#remotes::install_github("bbolker/nlme")
library(dplyr)
library(forcats)
library(furrr)
library(ggplot2)
library(gslnls)
library(magrittr)
library(mgcv)
library(nlme)
library(nlstools)
library(patchwork)
library(progressr)
library(readxl)
library(robustbase)
library(rsample)
library(scico)
library(sn)
library(stringr)
library(tibble)
library(tidyr)
library(WeightedROC)
library(writexl)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             legend.key.height = unit(0.85, "line"),
                             legend.spacing.y = unit(0, "line"),
                             legend.title = element_text(size = 10),
                             panel.border = element_blank(), 
                             plot.title = element_text(size = 10)))
htDiaOptions = tibble(folds = 10,
                      repetitions = 10,
                      includeInvestigatory = FALSE, # default to excluding plotting and other add ons in species scripts
                      retainModelThreshold = 10) # cross validation retains model objects if folds * repetitions is less than or equal to this threshold, e.g. 25 = retaining models up to and including 5x5 cross validation but sufficient DDR for loading all results may be an issue (5x5 easily exceeds 90 GB)

append_model_results = function(loadedResults, modelList, responseVariable, fitSet = "primary", fixedWeight = NA_real_)
{
  return(bind_rows(loadedResults,
                   bind_rows(lapply(modelList, get_list_stats, fitSet = fitSet, fixedWeight = fixedWeight)) %>%
                     mutate(responseVariable = responseVariable)))
}

check_plot_results = function(results)
{
  ggplot(results) +
    geom_histogram(aes(x = bias, fill = responseVariable)) +
    labs(x = "bias, m or cm") +
  ggplot(results) +
    geom_histogram(aes(x = mab, fill = responseVariable)) +
    labs(x = "MAB, m or cm") +
  ggplot(results) +
    geom_histogram(aes(x = mapb, fill = responseVariable)) +
    labs(x = "MAB, %") +
  ggplot(results) +
    geom_histogram(aes(x = mae, fill = responseVariable)) +
    labs(x = "MAE, m or cm") +
  ggplot(results) +
    geom_histogram(aes(x = mape, fill = responseVariable)) +
    labs(x = "MAE, %") +
  ggplot(results) +
    geom_histogram(aes(x = nse, fill = responseVariable)) +
    labs(x = "model efficiency") +
    scale_x_continuous(trans = scales::pseudo_log_trans()) +
  ggplot(results) +
    geom_histogram(aes(x = rmse, fill = responseVariable)) +
    labs(x = "RMSE, m or cm") +
  ggplot(results) +
    geom_histogram(aes(x = rmspe, fill = responseVariable)) +
    labs(x = "RMSE, %") +
  ggplot(results) +
    geom_histogram(aes(x = aic, fill = responseVariable)) +
    labs(x = "AIC") +
  ggplot(results) +
    geom_histogram(aes(x = aict, fill = responseVariable)) +
    labs(x = "AIC, ~t") +
  ggplot(results) +
    geom_histogram(aes(x = bic, fill = responseVariable)) +
    labs(x = "BIC") +
  ggplot(results) +
    geom_histogram(aes(x = bict, fill = responseVariable)) +
    labs(x = "BIC, ~t") +
  ggplot(results) +
    geom_histogram(aes(x = meanAbsolutePlantationEffect, fill = responseVariable)) +
    labs(x = "MAPE, m or cm") +
  ggplot(results) +
    geom_histogram(aes(x = meanAbsolutePercentPlantationEffect, fill = responseVariable)) +
    labs(x = "MAPE, %") +
  ggplot(results) +
    geom_histogram(aes(x = fitTimeInS, fill = responseVariable)) +
    labs(x = "fit time, s") +
    scale_x_log10() +
  guide_area() +
    theme(legend.justification = c(1, 0), legend.position = c(0.9, 0.1)) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    labs(y = "model fits", fill = NULL) &
    scale_fill_discrete(breaks = c("height", "DBH"))
}

confint_nlrob = function(regression, level = 0.99, df = df.residual(regression), trainingWeights = regression$weights)
{
  if (is.null(regression$rweights))
  {
    stop("regression$rweights is not set. Is the regression of type nlrob?")
  }
  if (is.null(trainingWeights))
  {
    stop("Either regression$weights is not set or trainingWeights was not specified.")
  }
  squaredDeviation = sum(trainingWeights * regression$rweights * residuals(regression)^2) / df.residual(regression)
  gradient = regression$m$gradient()
  levels = c((1 - level)/2, 1 - (1 - level)/2)
  parameterValues = regression$m$getPars()
  confidenceInterval = parameterValues + sqrt(diag(squaredDeviation * solve(t(gradient) %*% gradient))) %o% qt(p = levels, df = df)
  colnames(confidenceInterval) = sprintf("%g%%", 100*levels)
  rownames(confidenceInterval) = names(parameterValues)
  return(confidenceInterval)
}

create_model_stats = function(name, fittingMethod, fitSet = NA_character_, fixedWeight = NA_real_)
{
  if (is.null(name) | is.na(fittingMethod))
  {
    stop("Name and fittingMethod must be specified if modelStats is NULL.")
  }
  return(tibble(fitSet = fitSet, fitting = fittingMethod, fixedWeight = fixedWeight, name = name, 
                nNonPhysical = NA_real_, nObservations = NA_real_, nTaperImplausible = NA_real_, nValidation = NA_real_,
                isConverged = NA_real_, significant = NA_real_,
                aic = NA_real_, aict = NA_real_,
                bias = NA_real_, biasNaturalRegen = NA_real_, biasPlantation = NA_real_,
                bic = NA_real_, bict = NA_real_,
                mab = NA_real_, mapb = NA_real_,
                mape = NA_real_, mapeNaturalRegen = NA_real_, mapePlantation = NA_real_,
                meanAbsolutePlantationEffect = NA_real_, meanAbsolutePercentPlantationEffect = NA_real_,
                nse = NA_real_, nseNaturalRegen = NA_real_, nsePlantation = NA_real_,
                rmse = NA_real_, rmseNaturalRegen = NA_real_, rmsePlantation = NA_real_,
                rmspe = NA_real_, rmspeNaturalRegen = NA_real_, rmspePlantation = NA_real_,
                power = NA_real_, powerPlantation = NA_real_, adaptiveWeightFraction = NA_real_))
}

# wrap calls to fitting functions for consistency of arguments and use of get_*_error()
# This is a little fragile from a code maintenance perspective as the weights need to be column in data for R to flow them
# correctly but the risk appears low and worth the simplification elsewhere.
# future_map() sometimes works but doesn't reliably pass smooth parameters like constraints
fit_gam = function(name, formula, data, constraint = c(), family = gaussian(), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, bam = FALSE, mixed = FALSE, nthreads = 1, significant = TRUE, tDegreesOfFreedom = 8)
{
  if (bam & mixed)
  {
    stop("bam() does not support fixed effects. One of fit_gam()'s bam or mixed arguments can be true but not both.")
  }
  if (mixed & (nthreads > 1))
  {
    stop("gamm() does not support nthreads.")
  }
  
  responseVariable = formula[2] # displays as TotalHt or DBH but compares at TotalHt() or DBH()
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using ", if_else(mixed, "gamm", if_else(bam, "bam", "gam")), "()..."))
  progressBar = progressor(steps = folds * repetitions)
  
  # work around https://github.com/HenrikBengtsson/globals/issues/87 to enable GAM fitting using future_map()
  localFormula = local({ gamConstraint = constraint
                         formula(paste(deparse(formula), collapse = " ")) })
  if (responseVariable == "TotalHt()")
  {
    if (bam)
    {
      if ((folds == 1) & (repetitions == 1))
      {
        startFit = Sys.time()
        allFit = bam(formula = localFormula, data = data, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
        allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        allFitStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(allFit, allFitStats, returnModel))
      }
      
      fitFunction = function(dataFold)
      {
        startFit = Sys.time()
        trainingData = analysis(dataFold)
        model = bam(formula = localFormula, data = trainingData, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
        modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        modelStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(model, modelStats, returnModel))
      }
    }
    else
    {
      if ((folds == 1) & (repetitions == 1))
      {
        startFit = Sys.time()
        if (mixed)
        {
          allFit = gamm(formula = localFormula, data = data, family = family, method = "REML", weights = dbhWeight, verbosePQL = FALSE)
        } else {
          allFit = gam(formula = localFormula, data = data, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
        }
        allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        allFitStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(allFit, allFitStats, returnModel))
      }
      
      fitFunction = function(dataFold)
      {
        startFit = Sys.time()
        trainingData = analysis(dataFold)
        if (mixed)
        {
          model = gamm(formula = localFormula, data = trainingData, family = family, method = "REML", weights = dbhWeight, verbosePQL = FALSE)
        } else {
          model = gam(formula = localFormula, data = trainingData, family = family, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
        }
        modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        modelStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(model, modelStats, returnModel))
      }
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    if (bam)
    {
      if ((folds == 1) & (repetitions == 1))
      {
        startFit = Sys.time()
        allFit = bam(formula = localFormula, data = data, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
        allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        allFitStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(allFit, allFitStats, returnModel))
      }
      
      fitFunction = function(dataFold)
      {
        startFit = Sys.time()
        trainingData = analysis(dataFold)
        model = bam(formula = localFormula, data = trainingData, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
        modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        modelStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(model, modelStats, returnModel))
      }
    }
    else
    {
      if ((folds == 1) & (repetitions == 1))
      {
        startFit = Sys.time()
        if (mixed)
        {
          allFit = gamm(formula = localFormula, data = data, family = family, method = "REML", weights = heightWeight, verbosePQL = FALSE)
        } else {
          allFit = gam(formula = localFormula, data = data, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
        }
        allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        allFitStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(allFit, allFitStats, returnModel))
      }
      
      fitFunction = function(dataFold)
      {
        startFit = Sys.time()
        trainingData = analysis(dataFold)
        if (mixed)
        {
          model = gamm(formula = localFormula, data = trainingData, family = family, method = "REML", weights = heightWeight, verbosePQL = FALSE)
        } else {
          model = gam(formula = localFormula, data = trainingData, family = family, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
        }
        modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
        modelStats$fitTimeInS = get_elapsed_time(startFit)
        progressBar()
        return(get_fit_return_value(model, modelStats, returnModel))
      }
    }
  }
  
  # use map() instead of future_map() since GAM constraints fail to flow with future_map()
  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction, .options = furrr_options(seed = TRUE)))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

# nlme::gnls() has the same call capturing issue as nlme::lme() and nlme()
fit_gnls = function(name, modelFormula, data, start, control = gnlsControl(maxIter = 100), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = modelFormula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using gnls()..."))
  progressBar = progressor(steps = folds * repetitions)

  if (responseVariable == "TotalHt()")
  {
    startFit = Sys.time()
    allFit = do.call(gnls, list(model = modelFormula, data = data, start = start, weights = varPower(0.50, ~DBH | isPlantation), control = control))
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = varWeights(allFit$modelStruct$varStruct)
      allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$coefficients
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = do.call(gnls, list(model = modelFormula, data = trainingData, start = allFitParameters, weights = varPower(0.50, ~DBH | isPlantation), control = control))
      model$weights = varWeights(model$modelStruct$varStruct)
      modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    startFit = Sys.time()
    allFit = do.call(gnls, list(model = modelFormula, data = data, start = start, weights = varPower(0.50, ~TotalHt | isPlantation), control = control))
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = varWeights(allFit$modelStruct$varStruct)
      allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$coefficients
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = do.call(gnls, list(model = modelFormula, data = trainingData, start = allFitParameters, weights = varPower(0.50, ~TotalHt | isPlantation), control = control))
      model$weights = varWeights(model$modelStruct$varStruct)
      modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  
  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_gsl_nls = function(name, formula, data, start, control = gsl_nls_control(maxiter = 100), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using gsl_nls()..."))
  progressBar = progressor(steps = folds * repetitions)
  
  if (responseVariable == "TotalHt()")
  {
    startFit = Sys.time()
    allFit = gsl_nls(fn = formula, data = data, start = start, weights = dbhWeight, control = control)
    if ((folds == 1) & (repetitions == 1))
    {
      allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$m$getPars()
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = gsl_nls(fn = formula, data = trainingData, start = allFitParameters, weights = dbhWeight, control = control)
      modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    startFit = Sys.time()
    allFit = gsl_nls(fn = formula, data = data, start = start, weights = heightWeight, control = control)
    if ((folds == 1) & (repetitions == 1))
    {
      allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$m$getPars()
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = gsl_nls(fn = formula, data = trainingData, start = allFitParameters, weights = heightWeight, control = control)
      modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  
  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_lm = function(name, formula, data, folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using lm()..."))
  progressBar = progressor(steps = folds * repetitions)
  
  if (responseVariable == "TotalHt()")
  {
    if ((folds == 1) & (repetitions == 1))
    {
      startFit = Sys.time()
      allFit = lm(formula = formula, data = data, weights = dbhWeight)
      allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = lm(formula = formula, data = trainingData, offset = breastHeight, weights = dbhWeight)
      modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    if ((folds == 1) & (repetitions == 1))
    {
      startFit = Sys.time()
      allFit = lm(formula = formula, data = data, weights = heightWeight)
      allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = lm(formula = formula, data = trainingData, weights = heightWeight)
      modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  
  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

# weighted nlme() always raises warnings in CRAN version due to longstanding bug
# conLin$Xy * varWeights(object): longer object length is not a multiple of shorter object length
# https://stackoverflow.com/questions/74235304/error-while-using-the-weights-option-in-nlme-in-r
# https://github.com/bbolker/nlme
fit_nlme = function(name, modelFormula, data, fixedFormula, randomFormula, start, control = nlmeControl(maxIter = 100), folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = modelFormula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using nlme()..."))
  progressBar = progressor(steps = folds * repetitions)
  
  if (responseVariable == "TotalHt()")
  {
    startFit = Sys.time()
    # https://stackoverflow.com/questions/11778773/using-predict-in-a-function-call-with-nlme-objects-and-a-formula
    allFit = do.call(nlme, list(model = modelFormula, data = data, fixed = fixedFormula, random = randomFormula, groups = ~StandID, start = start, weights = varFixed(~1/dbhWeight), control = control))
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = varWeights(allFit$modelStruct$varStruct)
      allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$coefficients$fixed
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = do.call(nlme, list(model = modelFormula, data = trainingData, fixed = fixedFormula, random = randomFormula, groups = ~StandID, start = allFitParameters, weights = varFixed(~1/dbhWeight), control = control))
      model$weights = varWeights(model$modelStruct$varStruct)
      modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }
    
    startFit = Sys.time()
    allFit = do.call(nlme, list(model = modelFormula, data = data, fixed = fixedFormula, random = randomFormula, groups = ~StandID, start = start, weights = varFixed(~1/heightWeight), control = control))
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = varWeights(allFit$modelStruct$varStruct)
      allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$coefficients$fixed
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = do.call(nlme, list(model = modelFormula, data = trainingData, fixed = fixedFormula, random = randomFormula, groups = ~StandID, start = allFitParameters, weights = varFixed(~1/heightWeight), control = control))
      model$weights = varWeights(model$modelStruct$varStruct)
      modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  
  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

fit_nlrob = function(name, formula, data, start, control = nls.control(maxiter = 100), maxit = 50, folds = htDiaOptions$folds, repetitions = htDiaOptions$repetitions, returnModel = folds * repetitions <= htDiaOptions$retainModelThreshold, significant = TRUE, tDegreesOfFreedom = 8)
{
  responseVariable = formula[2]
  message(paste0("Fitting ", name, " for ", folds, "x", repetitions, " ", responseVariable, " using nlrob()..."))
  progressBar = progressor(steps = folds * repetitions)

  if (responseVariable == "TotalHt()")
  {
    startFit = Sys.time()
    allFit = nlrob(formula = formula, data = data, maxit = maxit, start = start, weights = dbhWeight, control = control)
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = data$dbhWeight
      allFitStats = get_height_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$m$getPars()
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = nlrob(formula = formula, data = trainingData, maxit = maxit, start = allFitParameters, weights = dbhWeight, control = control)
      model$weights = trainingData$dbhWeight
      modelStats = get_height_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }
  else
  {
    if (responseVariable != "DBH()")
    {
      stop("Expected response variable to be DBH.")
    }

    startFit = Sys.time()
    allFit = nlrob(formula = formula, data = data, maxit = maxit, start = start, weights = heightWeight, control = control)
    if ((folds == 1) & (repetitions == 1))
    {
      allFit$weights = data$heightWeight
      allFitStats = get_dbh_stats(name = name, model = allFit, trainingData = data, validationData = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      allFitStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(allFit, allFitStats, returnModel))
    }
    
    allFitParameters = allFit$m$getPars()
    fitFunction = function(dataFold)
    {
      startFit = Sys.time()
      trainingData = analysis(dataFold)
      model = nlrob(formula = formula, data = trainingData, maxit = maxit, start = allFitParameters, weights = heightWeight, control = control)
      model$weights = trainingData$heightWeight
      modelStats = get_dbh_stats(name = name, model = model, trainingData = trainingData, validationData = assessment(dataFold), significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
      modelStats$fitTimeInS = get_elapsed_time(startFit)
      progressBar()
      return(get_fit_return_value(model, modelStats, returnModel))
    }
  }

  splitsAndFits = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, fitFunction))
  return(get_cross_validation_return_value(splitsAndFits, returnModel))
}

get_adaptive_weighting = function(model)
{
  switch(class(model)[1],
         "gnls" = 
           { 
             # for now, assume a varStruct was used
             if (is.null(model$modelStruct) | is.null(model$modelStruct$varStruct) | (length(model$modelStruct$varStruct) < 1))
             {
               stop("gnls() model does not have a varStruct.")
             }
             adaptiveWeightFraction = 1 
           },
         "nlrob" = 
           { 
             adaptiveWeightFraction = sum(model$rweights != 1) / nobs(model) 
           },
         { 
           # assume fixed weights with gam(), gsl_nls, and nlme()
           adaptiveWeightFraction = 0 
         })
  return(adaptiveWeightFraction)
}

get_cross_validation_return_value = function(splitsAndFits, returnModel)
{
  if (returnModel)
  {
    # if full models are retained retain also the training-validation data splits for each model fit
    return(splitsAndFits)
  }
  # if only model statistics are kept then drop the data splits and bind together each fit's single 
  # row statistics tibble into one tibble with kr rows (k-folds x r-repetitions = kr fits)
  return(bind_rows(splitsAndFits$fit) %>% 
           mutate(repetition = as.numeric(str_replace(splitsAndFits$id, "Repeat", "")), 
                  fold = as.numeric(str_replace(splitsAndFits$id2, "Fold", ""))))
}

get_dbh_stats = function(name, model, trainingData, validationData, validationWeights = validationData$heightWeight, significant = TRUE, tDegreesOfFreedom = 8)
{
  dbhModelStats = create_model_stats(NULL, name = name, fittingMethod = class(model)[1]) %>% 
    mutate(name = name, 
           coefficients = get_model_coefficients(model), 
           isConverged = is_model_converged(model),
           adaptiveWeightFraction = get_adaptive_weighting(model))
  if (dbhModelStats$isConverged == FALSE)
  {
    warning(paste0(dbhModelStats$fitting, " ", formula(model)[2], " model using ", name, " is not converged."))
  }

  if (is(model, "gamm"))
  {
    predictedDbh = predict(model$gam, validationData) 
  } else if (is(model, "nlme"))
  {
    # see notes in get_height_stats()
    nlmeValidationData = left_join(bind_cols(validationData, bind_rows(model$coefficients$fixed)),
                                   bind_cols(StandID = as.numeric(rownames(model$coefficients$random$StandID), levels = levels(validationData$StandID)), model$coefficients$random$StandID),
                                   by = "StandID")
    randomEffectsColumns = setdiff(names(nlmeValidationData), names(validationData))
    nlmeValidationData %<>% replace_na(setNames(as.list(rep(0, length(randomEffectsColumns))), randomEffectsColumns))
    predictedDbh = eval(getCovariateFormula(model)[[2]], nlmeValidationData)
  } else {
    predictedDbh = predict(model, validationData)
  }
  validationResiduals = predictedDbh - validationData$DBH
  dbhByHeightClass = validationData %>% mutate(residuals = validationResiduals) %>%
    group_by(heightClass) %>%
    summarize(n = sum(TreeCount),
              nNaturalRegen = sum(TreeCount * (isPlantation == FALSE)),
              nPlantation = sum(TreeCount * isPlantation),
              meanBiasPerTree = sum(TreeCount * residuals) / n,
              meanBiasPerTreePct = 100 * sum(TreeCount * residuals / DBH) / n,
              meanDbh = sum(TreeCount * DBH) / n,
              meanNaturalRegenDbh = if_else(nNaturalRegen > 0, sum(TreeCount * DBH * (isPlantation == FALSE)) / nNaturalRegen, NA_real_),
              meanPlantationDbh = if_else(nPlantation > 0, sum(TreeCount * DBH * isPlantation) / nPlantation, NA_real_),
              minPlantationNaturalRegenN = min(nPlantation, nNaturalRegen),
              plantationEffect = meanPlantationDbh - meanNaturalRegenDbh,
              plantationEffectPct = 100 * plantationEffect / meanDbh,
              .groups = "drop") %>%
    filter(n > 0)

  nObservations = sum(trainingData$TreeCount)
  if (is(model, "gam"))
  {
    effectiveDegreesOfFreedom = sum(model$edf) + 1 # for GAMs use indicated effective degrees of freedom
  } else if(is(model, "gamm")) {
    effectiveDegreesOfFreedom = sum(model$gam$edf) + 1 # gamm effective degrees of freedom include random effect degrees of freedom from s(bs = "re) but not from random = list()
  } else {
    effectiveDegreesOfFreedom = length(coef(model)) + 1 # for linear and nonlinear regressions assume one degree of freedom per model parameter, coef(nlme()) is data.frame with length() returning the number of columns
  }
  residualDegreesOfFreedom = nObservations - effectiveDegreesOfFreedom # gamm() and nlme() don't implement df.residual(model), simplest just to calculate it here
  
  # logLik.lm() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/logLik.R)
  #  log likelihood = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w*res^2))))
  # logLik.nls() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/nls.R)
  #  log likelihood = -N/2 * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(regression$m$resid()^2)))
  #                 = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(regression$m$resid()^2)))) if no weights are zero
  standardDeviation = sqrt(1/residualDegreesOfFreedom * sum(validationWeights * validationResiduals^2)) / sqrt(validationWeights)
  logLikelihoodGaussian = sum(validationData$TreeCount * dnorm(validationResiduals, sd = standardDeviation, log = TRUE))
  logLikelihoodT = sum(validationData$TreeCount * (dt(validationResiduals / standardDeviation, df = tDegreesOfFreedom, log = TRUE) - log(standardDeviation)))

  heightDiameterRatio = validationData$TotalHt / (0.01 * predictedDbh)
  speciesLimits = get_species_limits(validationData)
  validationTreeCountTotal = sum(validationData$TreeCount)
  dbhModelStats$aic = -2*logLikelihoodGaussian + 2 * effectiveDegreesOfFreedom # calculate AIC and BIC manually because nlrob objects implement weighting differently from nls and gslnls
  dbhModelStats$aict = -2*logLikelihoodT + 2 * effectiveDegreesOfFreedom
  dbhModelStats$bias = sum(validationData$TreeCount * validationResiduals) / validationTreeCountTotal
  dbhModelStats$bic = -2*logLikelihoodGaussian + effectiveDegreesOfFreedom * log(nObservations)
  dbhModelStats$bict = -2*logLikelihoodT + effectiveDegreesOfFreedom * log(nObservations)
  dbhModelStats$mae = sum(validationData$TreeCount * abs(validationResiduals)) / validationTreeCountTotal
  dbhModelStats$mab = sum(dbhByHeightClass$n * abs(dbhByHeightClass$meanBiasPerTree)) / sum(dbhByHeightClass$n)
  dbhModelStats$mapb = sum(dbhByHeightClass$n * abs(dbhByHeightClass$meanBiasPerTreePct)) / sum(dbhByHeightClass$n)
  dbhModelStats$mape = 100 * sum(validationData$TreeCount * abs(validationResiduals / validationData$DBH)) / validationTreeCountTotal
  dbhModelStats$meanAbsolutePlantationEffect = sum(dbhByHeightClass$minPlantationNaturalRegenN * abs(dbhByHeightClass$plantationEffect), na.rm = TRUE) / sum(dbhByHeightClass$minPlantationNaturalRegenN * (is.na(dbhByHeightClass$plantationEffect) == FALSE), na.rm = TRUE)
  dbhModelStats$meanAbsolutePercentPlantationEffect = sum(dbhByHeightClass$minPlantationNaturalRegenN * abs(dbhByHeightClass$plantationEffectPct), na.rm = TRUE) / sum(dbhByHeightClass$minPlantationNaturalRegenN * (is.na(dbhByHeightClass$plantationEffect) == FALSE), na.rm = TRUE)
  dbhModelStats$nNonPhysical = sum(is.na(predictedDbh) | (predictedDbh < speciesLimits$dbhMin) | (predictedDbh > speciesLimits$dbhMax))
  dbhModelStats$nObservations = nObservations
  dbhModelStats$nTaperImplausible = sum((heightDiameterRatio < speciesLimits$heightDiameterRatioMin) | (heightDiameterRatio > speciesLimits$heightDiameterRatioMax), na.rm = TRUE)
  dbhModelStats$nValidation = validationTreeCountTotal
  dbhModelStats$nse = 1 - sum(validationData$TreeCount * validationResiduals^2) / sum(validationData$TreeCount * (validationData$DBH - sum(validationData$TreeCount * validationData$DBH) / validationTreeCountTotal)^2)
  dbhModelStats$rmse = sqrt(sum(validationData$TreeCount * validationResiduals^2) / validationTreeCountTotal)
  dbhModelStats$rmspe = 100 * sqrt(sum(validationData$TreeCount * (validationResiduals / validationData$DBH)^2) / validationTreeCountTotal)
  dbhModelStats$significant = significant
  
  naturalRegenIndices = which(validationData$isPlantation == FALSE)
  dbhNaturalRegen = validationData$DBH[naturalRegenIndices]
  predictedDbhNaturalRegen = predictedDbh[naturalRegenIndices]
  residualsNaturalRegen = predictedDbhNaturalRegen - dbhNaturalRegen
  naturalRegenTreeCount = validationData$TreeCount[naturalRegenIndices]
  naturalRegenTreeCountTotal = sum(naturalRegenTreeCount)
  dbhModelStats$biasNaturalRegen = sum(naturalRegenTreeCount * residualsNaturalRegen) / naturalRegenTreeCountTotal
  dbhModelStats$maeNaturalRegen = sum(naturalRegenTreeCount * abs(residualsNaturalRegen)) / naturalRegenTreeCountTotal
  dbhModelStats$mapeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / predictedDbhNaturalRegen)) / naturalRegenTreeCountTotal
  dbhModelStats$nseNaturalRegen = 1 - sum(naturalRegenTreeCount * residualsNaturalRegen^2) / sum(naturalRegenTreeCount * (dbhNaturalRegen - sum(naturalRegenTreeCount * dbhNaturalRegen) / naturalRegenTreeCountTotal)^2)
  dbhModelStats$paeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / dbhNaturalRegen)) / naturalRegenTreeCountTotal
  dbhModelStats$rmseNaturalRegen = sqrt(sum(naturalRegenTreeCount * residualsNaturalRegen^2) / naturalRegenTreeCountTotal)
  dbhModelStats$rmspeNaturalRegen = 100 * sqrt(sum(naturalRegenTreeCount * (residualsNaturalRegen / predictedDbhNaturalRegen)^2) / naturalRegenTreeCountTotal)
  
  plantationIndices = which(validationData$isPlantation)
  dbhPlantation = validationData$DBH[plantationIndices]
  predictedDbhPlantation = predictedDbh[plantationIndices]
  residualsPlantation = predictedDbhPlantation - dbhPlantation
  plantationTreeCount = validationData$TreeCount[plantationIndices]
  plantationTreeCountTotal = sum(plantationTreeCount)
  dbhModelStats$biasPlantation = sum(plantationTreeCount * residualsPlantation) / plantationTreeCountTotal
  dbhModelStats$maePlantation = sum(plantationTreeCount * abs(residualsPlantation)) / plantationTreeCountTotal
  dbhModelStats$mapePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / predictedDbhPlantation)) / plantationTreeCountTotal
  dbhModelStats$nsePlantation = 1 - sum(plantationTreeCount * residualsPlantation^2) / sum(plantationTreeCount * (dbhPlantation - sum(plantationTreeCount * dbhPlantation) / plantationTreeCountTotal)^2)
  dbhModelStats$paePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / dbhPlantation)) / plantationTreeCountTotal
  dbhModelStats$pearsonPlantation = cor(predictedDbhPlantation, dbhPlantation)
  dbhModelStats$rmsePlantation = sqrt(sum(plantationTreeCount * residualsPlantation^2) / plantationTreeCountTotal)
  dbhModelStats$rmspePlantation = 100 * sqrt(sum(plantationTreeCount * (residualsPlantation / predictedDbhPlantation)^2) / plantationTreeCountTotal)
  
  return(dbhModelStats)
}

get_elapsed_time = function(start)
{
  return(as.numeric(difftime(Sys.time(), start, units = "secs")))
}

get_fit_return_value = function(model, modelStats, returnModel)
{
  if (returnModel)
  {
    model$stats = modelStats
    return(model)
  } else {
    return(modelStats)
  }
}

get_height_stats = function(name, model, trainingData, validationData, validationWeights = validationData$dbhWeight, significant = TRUE, tDegreesOfFreedom = 8)
{
  heightModelStats = create_model_stats(name = name, fittingMethod = class(model)[1]) %>% 
    mutate(coefficients = get_model_coefficients(model), 
           isConverged = is_model_converged(model), 
           significant = significant, 
           adaptiveWeightFraction = get_adaptive_weighting(model))
  if (heightModelStats$isConverged == FALSE)
  {
    warning(paste0(heightModelStats$fitting, " ", formula(model)[2], " model using ", name, " is not converged."))
  }

  dbhClassSize = 10 # cm
  if (is(model, "gamm"))
  {
    predictedHeight = predict(model$gam, validationData) 
  } else if (is(model, "nlme"))
  {
    # predict.nlme(model, newdata) fails in nlme::asOneFormula() by default
    # predict.nlme(model, newdata) returns NAs for unknown random effects with do.call() workaround used by fit_nlme()
    # To avoid predicting NAs, this code path bypasses predict.nlme() and replaceds unknown random effects with zeros.
    # It relies on the do.call() workaround to be able to retrieve the model formula to evaluate.
    nlmeValidationData = left_join(bind_cols(validationData, bind_rows(model$coefficients$fixed)),
                                   bind_cols(StandID = as.numeric(rownames(model$coefficients$random$StandID), levels = levels(validationData$StandID)), model$coefficients$random$StandID),
                                   by = "StandID")
    randomEffectsColumns = setdiff(names(nlmeValidationData), names(validationData))
    nlmeValidationData %<>% replace_na(setNames(as.list(rep(0, length(randomEffectsColumns))), randomEffectsColumns))
    predictedHeight = eval(getCovariateFormula(model)[[2]], nlmeValidationData)
  } else {
    predictedHeight = predict(model, validationData)
  }
  validationResiduals = predictedHeight - validationData$TotalHt
  
  heightByDbhClass = validationData %>% mutate(residuals = validationResiduals) %>%
    group_by(dbhClass) %>%
    summarize(n = sum(TreeCount),
              nNaturalRegen = sum(TreeCount * (isPlantation == FALSE)),
              nPlantation = sum(TreeCount * isPlantation),
              meanBiasPerTree = sum(TreeCount * residuals) / n,
              meanBiasPerTreePct = 100 * sum(TreeCount * residuals / TotalHt) / n,
              meanHeight = sum(TreeCount * TotalHt) / n,
              meanNaturalRegenHeight = if_else(nNaturalRegen > 0, sum(TreeCount * TotalHt * (isPlantation == FALSE)) / nNaturalRegen, NA_real_),
              meanPlantationHeight = if_else(nPlantation > 0, sum(TreeCount * TotalHt * isPlantation) / nPlantation, NA_real_),
              minPlantationNaturalRegenN = min(nPlantation, nNaturalRegen),
              plantationEffect = meanPlantationHeight - meanNaturalRegenHeight,
              plantationEffectPct = 100 * plantationEffect / meanHeight,
              .groups = "drop") %>%
    filter(n > 0)

  nObservations = sum(trainingData$TreeCount)
  if (is(model, "gam"))
  {
    effectiveDegreesOfFreedom = sum(model$edf) + 1
  } else if (is(model, "gamm"))
  {
    effectiveDegreesOfFreedom = sum(model$gam$edf) + 1
  } else {
    effectiveDegreesOfFreedom = length(coef(model)) + 1
  }
  residualDegreesOfFreedom = nObservations - effectiveDegreesOfFreedom
  
  standardDeviation = sqrt(1/residualDegreesOfFreedom * sum(validationWeights * validationResiduals^2)) / sqrt(validationWeights)
  logLikelihoodGaussian = sum(validationData$TreeCount * dnorm(validationResiduals, sd = standardDeviation, log = TRUE))
  logLikelihoodT = sum(validationData$TreeCount * dt(validationResiduals / standardDeviation, df = tDegreesOfFreedom, log = TRUE) - log(standardDeviation))
  
  heightDiameterRatio = predictedHeight / (0.01 * validationData$DBH)
  speciesLimits = get_species_limits(validationData)
  validationTreeCountTotal = sum(validationData$TreeCount)
  heightModelStats$aic = -2*logLikelihoodGaussian + 2 * effectiveDegreesOfFreedom # see get_dbh_stats(): same nlrob issue
  heightModelStats$aict = -2*logLikelihoodT + 2 * effectiveDegreesOfFreedom
  heightModelStats$bias = sum(validationData$TreeCount * validationResiduals) / validationTreeCountTotal
  heightModelStats$bic = -2*logLikelihoodGaussian + effectiveDegreesOfFreedom * log(nObservations)
  heightModelStats$bict = -2*logLikelihoodT + effectiveDegreesOfFreedom * log(nObservations)
  heightModelStats$mab = sum(heightByDbhClass$n * abs(heightByDbhClass$meanBiasPerTree)) / sum(heightByDbhClass$n)
  heightModelStats$mapb = sum(heightByDbhClass$n * abs(heightByDbhClass$meanBiasPerTreePct)) / sum(heightByDbhClass$n)
  heightModelStats$mae = sum(validationData$TreeCount * abs(validationResiduals)) / validationTreeCountTotal
  heightModelStats$mape = 100 * sum(validationData$TreeCount * abs(validationResiduals / validationData$TotalHt)) / validationTreeCountTotal
  heightModelStats$meanAbsolutePlantationEffect = sum(heightByDbhClass$minPlantationNaturalRegenN * abs(heightByDbhClass$plantationEffect), na.rm = TRUE) / sum(heightByDbhClass$minPlantationNaturalRegenN * (is.na(heightByDbhClass$plantationEffect) == FALSE), na.rm = TRUE)
  heightModelStats$meanAbsolutePercentPlantationEffect = sum(heightByDbhClass$minPlantationNaturalRegenN * abs(heightByDbhClass$plantationEffectPct), na.rm = TRUE) / sum(heightByDbhClass$minPlantationNaturalRegenN * (is.na(heightByDbhClass$plantationEffect) == FALSE), na.rm = TRUE)
  heightModelStats$nNonPhysical = sum(is.na(predictedHeight) | (predictedHeight < speciesLimits$heightMin) | (predictedHeight > speciesLimits$heightMax))
  heightModelStats$nObservations = nObservations
  heightModelStats$nTaperImplausible = sum((heightDiameterRatio < speciesLimits$heightDiameterRatioMin) | (heightDiameterRatio > speciesLimits$heightDiameterRatioMax), na.rm = TRUE)
  heightModelStats$nValidation = validationTreeCountTotal
  heightModelStats$nse = 1 - sum(validationData$TreeCount * validationResiduals^2) / sum((validationData$TreeCount * validationData$TotalHt - sum(validationData$TreeCount * validationData$TotalHt) / validationTreeCountTotal)^2)
  heightModelStats$rmse = sqrt(sum(validationData$TreeCount * validationResiduals^2) / validationTreeCountTotal)
  heightModelStats$rmspe = 100 * sqrt(sum(validationData$TreeCount * (validationResiduals / validationData$TotalHt)^2) / validationTreeCountTotal)
  heightModelStats$significant = significant

  naturalRegenIndices = which(validationData$isPlantation == FALSE)
  heightNaturalRegen = validationData$TotalHt[naturalRegenIndices]
  predictedHeightNaturalRegen = predictedHeight[naturalRegenIndices]
  residualsNaturalRegen = predictedHeightNaturalRegen - heightNaturalRegen
  naturalRegenTreeCount = validationData$TreeCount[naturalRegenIndices]
  naturalRegenTreeCountTotal = sum(naturalRegenTreeCount)
  heightModelStats$biasNaturalRegen = sum(naturalRegenTreeCount * residualsNaturalRegen) / naturalRegenTreeCountTotal
  heightModelStats$maeNaturalRegen = sum(naturalRegenTreeCount * abs(residualsNaturalRegen)) /naturalRegenTreeCountTotal
  heightModelStats$mapeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / heightNaturalRegen)) / naturalRegenTreeCountTotal
  heightModelStats$nseNaturalRegen = 1 - sum(naturalRegenTreeCount * residualsNaturalRegen^2) / sum((heightNaturalRegen - sum(naturalRegenTreeCount * heightNaturalRegen) / naturalRegenTreeCountTotal)^2)
  heightModelStats$paeNaturalRegen = 100 * sum(naturalRegenTreeCount * abs(residualsNaturalRegen / heightNaturalRegen)) / naturalRegenTreeCountTotal
  heightModelStats$rmseNaturalRegen = sqrt(sum(naturalRegenTreeCount * residualsNaturalRegen^2) / naturalRegenTreeCountTotal)
  heightModelStats$rmspeNaturalRegen = 100 * sqrt(sum(naturalRegenTreeCount * (residualsNaturalRegen / heightNaturalRegen)^2) / naturalRegenTreeCountTotal)
  
  plantationIndices = which(validationData$isPlantation)
  heightPlantation = validationData$TotalHt[plantationIndices]
  predictedHeightPlantation = predictedHeight[plantationIndices]
  residualsPlantation = predictedHeightPlantation - heightPlantation
  plantationTreeCount = validationData$TreeCount[plantationIndices]
  plantationTreeCountTotal = sum(plantationTreeCount)
  heightModelStats$biasPlantation = sum(plantationTreeCount * residualsPlantation) / plantationTreeCountTotal
  heightModelStats$maePlantation = sum(plantationTreeCount * abs(residualsPlantation)) / plantationTreeCountTotal
  heightModelStats$mapePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / heightPlantation)) / plantationTreeCountTotal
  heightModelStats$nsePlantation = 1 - sum(plantationTreeCount * residualsPlantation^2) / sum((plantationTreeCount * heightPlantation - sum(plantationTreeCount * heightPlantation) / plantationTreeCountTotal)^2)
  heightModelStats$paePlantation = 100 * sum(plantationTreeCount * abs(residualsPlantation / heightPlantation) / plantationTreeCountTotal)
  heightModelStats$rmsePlantation = sqrt(sum(plantationTreeCount * residualsPlantation^2) / plantationTreeCountTotal)
  heightModelStats$rmspePlantation = 100 * sqrt(sum(plantationTreeCount * (residualsPlantation / heightPlantation)^2) / plantationTreeCountTotal)
  
  return(heightModelStats)
}

get_preferred_model_linetype_legend = function()
{
  return(ggplot() +
           geom_segment(aes(x = x, xend = xend, y = y, yend = yend, color = color, linetype = linetype), tibble(x = 0, xend = 1, y = 0, yend = 0, color = c(FALSE, TRUE, "reference curve"), linetype = c(FALSE, TRUE, "reference curve")), alpha = 0, show.legend = TRUE) +
           guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1, override.aes = list(alpha = 1))) +
           labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
           scale_color_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("grey25", "grey25", "grey70")) +
           theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.direction = "horizontal", legend.position = c(0.5, 0.5), panel.grid = element_blank()))
}

get_list_coefficients = function(modelCrossValidationListOrStatsTibble, fitSet = "primary", fixedWeight = NA_real_)
{
  # complete vfold_cv return tibble: vfold_cv, tbl_df, tbl, data.frame
  # vfold_cv() %>% select(-splits): tbl_df, tbl, data.frame
  if (is(modelCrossValidationListOrStatsTibble, "vfold_cv"))
  {
    # cross validation tibble
    coefficients = bind_rows(lapply(modelCrossValidationListOrStatsTibble$fit, function(fit) 
      { 
        return(fit$stats$coefficients %>% 
                 mutate(name = fit$stats$name,
                        fitting = fit$stats$fitting,
                        significant = fit$stats$significant))
      })) %>%
      mutate(repetition = as.numeric(str_replace(modelCrossValidationListOrStatsTibble$id, "Repeat", "")), 
             fold = as.numeric(str_replace(modelCrossValidationListOrStatsTibble$id2, "Fold", "")))
  } else if (is(modelCrossValidationListOrStatsTibble, "tbl_df"))
  {
    # model statistics tibble
    coefficients = modelCrossValidationListOrStatsTibble$coefficients %>% 
      mutate(name = modelCrossValidationListOrStatsTibble$name,
             fitting = modelCrossValidationListOrStatsTibble$fitting,
             significant = modelCrossValidationListOrStatsTibble$significant,
             repetition = modelCrossValidationListOrStatsTibble$repetition, 
             fold = modelCrossValidationListOrStatsTibble$fold)
  } else 
  {
    # single model
    coefficients = modelCrossValidationListOrStatsTibble$stats$coefficients %>% 
      mutate(name = modelCrossValidationListOrStatsTibble$stats$name, 
             fitting = modelCrossValidationListOrStatsTibble$stats$fitting, 
             significant = modelCrossValidationListOrStatsTibble$stats$significant,
             repetition = 1, fold = 1)
  }
  
  fitSetArgument = fitSet # change names for disambiguation (otherwise mutate() assigns the coefficient tibble's fitSet column to itself)
  fixedWeightArgument = fixedWeight
  coefficients %<>% mutate(fitSet = fitSetArgument,
                           fixedWeight = fixedWeightArgument) %>%
   relocate(fitSet, fixedWeight, name, repetition, fold)
  return(coefficients)
}

get_list_stats = function(modelCrossValidationListOrStatsTibble, fitSet = "primary", fixedWeight = NA_real_)
{
  if (is(modelCrossValidationListOrStatsTibble, "vfold_cv"))
  {
    # cross validation list
    stats = bind_rows(lapply(modelCrossValidationListOrStatsTibble$fit, get_model_stats)) %>% 
      mutate(repetition = as.numeric(str_replace(modelCrossValidationListOrStatsTibble$id, "Repeat", "")), 
             fold = as.numeric(str_replace(modelCrossValidationListOrStatsTibble$id2, "Fold", "")))
  } else if (is(modelCrossValidationListOrStatsTibble, "tbl_df"))
  {
    # model statistics tibble
    stats = modelCrossValidationListOrStatsTibble %>%
      mutate(fitSet = fitSet, fixedWeight = fixedWeight)
  } else {
    # single model
    stats = get_model_stats(modelCrossValidationListOrStatsTibble) %>% mutate(repetition = 1, fold = 1)
  }
  
  fitSetArgument = fitSet # change names for disambiguation (otherwise mutate() assigns the coefficient tibble's fitSet column to itself)
  fixedWeightArgument = fixedWeight
  stats %<>% mutate(fitSet = fitSetArgument,
                    fixedWeight = fixedWeightArgument) %>%
    select(-coefficients) %>%
    relocate(fitSet, fixedWeight, name, repetition, fold)
  return(stats)
}

get_model_coefficients = function(model)
{
  if (is.null(model$coefficients))
  {
    if (is(model, "gamm"))
    {
      # s(StandID, bs = "re")
      coefficients = bind_cols(!!!coefficients(model$gam),
                               !!!set_names(model$lme$coefficients$fixed, names(model$lme$coefficients$fixed)),
                               model$lme$coefficients$random$g,
                               model$lme$coefficients$random$g.0,
                               model$lme$coefficients$random$g.1)
      # random = list(StandID = ~1) is used
      #coefficients = bind_cols(!!!coefficients(model$gam),
      #                         !!!set_names(model$lme$coefficients$fixed, names(model$lme$coefficients$fixed)),
      #                         model$lme$coefficients$random$g,
      #                         model$lme$coefficients$random$g.0,
      #                         !!!set_names(model$lme$coefficients$random$StandID, str_replace(rownames(model$lme$coefficients$random$StandID), "1/1/", "XrStand")))
    } else if (is.null(model$m) == FALSE) {
      coefficients = tibble(!!!set_names(model$m$getPars(), names(model$m$getPars())))
    } else {
      stop(paste("Regression for", model$name, "lacks a both a coefficients property and an m property."))
    }
  }
  else
  {
    if (is(model, "nlme"))
    {
      coefficients = bind_cols(tibble(!!!set_names(model$coefficients$fixed, names(model$coefficients$fixed))),
                               tibble(!!!set_names(model$coefficients$random$StandID, str_c(colnames(model$coefficients$random$StandID), rownames(model$coefficients$random$StandID)))))
    } else {
      coefficients = tibble(!!!set_names(model$coefficients, names(model$coefficients)))
    }
  }

  if ("(Intercept)" %in% names(coefficients))
  {
    coefficients %<>% rename(a0 = `(Intercept)`)
  }
  
  if (is(model, "lm"))
  {
    coefficientNames = names(coefficients)
    # change linear regression nomenclature to match nonlinear
    if ("DBH" %in% coefficientNames) # regressions for height imputation
    {
      coefficients %<>% rename(a1 = DBH)
      if ("I(isPlantation * DBH)" %in% coefficientNames)
      {
        coefficients %<>% rename(a1p = `I(isPlantation * DBH)`)
      }
      if ("I(DBH^2)" %in% coefficientNames)
      {
        coefficients %<>% rename(a2 = `I(DBH^2)`)
      }
      if ("I(isPlantation * DBH^2)" %in% coefficientNames)
      {
        coefficients %<>% rename(a2p = `I(isPlantation * DBH^2)`)
      }
    }
    if ("I(TotalHt - 1.37)" %in% coefficientNames) # regressions for diameter imputation
    {
      coefficients %<>% rename(a1 = `I(TotalHt - 1.37)`)
      if ("I(isPlantation * (TotalHt - 1.37))" %in% coefficientNames)
      {
        coefficients %<>% rename(a1p = `I(isPlantation * (TotalHt - 1.37))`)
      }
      if ("I((TotalHt - 1.37)^2)" %in% coefficientNames)
      {
        coefficients %<>% rename(a2 = `I((TotalHt - 1.37)^2)`)
      }
      if ("I(isPlantation * (TotalHt - 1.37)^2)" %in% coefficientNames)
      {
        coefficients %<>% rename(a2p = `I(isPlantation * (TotalHt - 1.37)^2)`)
      }
    }
  }
  
  return(coefficients)
}

get_model_stats = function(modelOrStats)
{
  modelStats = modelOrStats
  if (is.null(modelOrStats$stats) == FALSE)
  {
    modelStats = modelOrStats$stats
  }
  return(modelStats)
}

get_species_limits = function(trees)
{
  return(tibble(speciesGroup = trees$speciesGroup, DBH = trees$DBH) %>% 
           mutate(heightMax = case_when(speciesGroup == "DF" ~ 85,
                                        speciesGroup == "RA" ~ 60,
                                        speciesGroup == "WH" ~ 80,
                                        speciesGroup == "BM" ~ 50,
                                        speciesGroup == "OM" ~ 45,
                                        speciesGroup == "RC" ~ 80,
                                        speciesGroup == "other" ~ 85), # m
                  heightMin = 1.37,
                  dbhMax = case_when(speciesGroup == "DF" ~ 250,
                                     speciesGroup == "RA" ~ 150,
                                     speciesGroup == "WH" ~ 200,
                                     speciesGroup == "BM" ~ 200,
                                     speciesGroup == "OM" ~ 250,
                                     speciesGroup == "RC" ~ 250,
                                     speciesGroup == "other" ~ 225), # cm
                  dbhMin = 0.3, # cm, for simplicity assume a minimum leader diameter of 3 mm and a maximum DBH of 3 m though, for some species, the minimum is larger and the maximum smaller (could also consider broken top increases in Sitka spruce, though those aren't predicted here
                  heightDiameterRatioMax = case_when(speciesGroup == "DF" ~ 1 + 1000 * DBH^-0.51,
                                                     speciesGroup == "RA" ~ 1 + 1400 * DBH^-0.65,
                                                     speciesGroup == "WH" ~ 5 + 1000 * DBH^-0.56,
                                                     speciesGroup == "BM" ~ 1 + 1000 * DBH^-0.63,
                                                     speciesGroup == "OM" ~ 1 + 800 * DBH^-0.60,
                                                     speciesGroup == "RC" ~ 5 + 1000 * DBH^-0.59,
                                                     speciesGroup == "other" ~ 10 + 1000 * DBH^-0.61),
                  heightDiameterRatioMin = 2))
}

impute_basal_area = function(Species, heightInM, isPlantation)
{
  # preferred fits from ends of PSME.R, ALRU2.R, TSHE.R, ...
  basalAreaInM2 = case_match(Species,
                             "DF" ~ 0.5138900 * (exp(0.0008023 * (heightInM - 1.37)^(1.7625191 - 0.0677262 * isPlantation)) - 1),
                             "RA" ~ 6.033354e+02 * (exp(8.625243e-07 * (heightInM - 1.37)^(1.865767e+00 - 1.471688e-01 * isPlantation)) - 1),
                             "WH" ~ 0.0002814 * (heightInM - 1.37)^1.9664137,
                             "BM" ~ 2.230479e+02 * (exp(1.331652e-06 * (heightInM - 1.37)^(2.212715e+00 - 9.147243e-01 * isPlantation)) - 1),
                             "OM" ~ (0.0006551 - 0.0003407 * isPlantation) * (heightInM - 1.37)^1.8674774,
                             "RC" ~ 0.0003335 * (heightInM - 1.37)^(2.1341090 - 0.0438769 * isPlantation),
                             .default = 6.313626e-04 * (heightInM - 1.37)^(-8.764352e-05 + 1.556025e-01 * isPlantation))
  basalAreaInM2 = if_else(basalAreaInM2 < 0.25 * pi * (0.01 * 2.54)^2, 0.25 * pi * (0.01 * 2.54)^2, basalAreaInM2) # clamp regressions to minimum cruised basal area: blocks implied negative DBHes
  return(basalAreaInM2)
}

impute_height = function(Species, DBH, isPlantation)
{
  # preferred base nonlinear regression form for height prediction by species
  return(case_match(Species,
                    "DF" ~ 1.37 + (0.001005 + 0.158058 * isPlantation)*DBH^((5.392558 - 3.115258  * isPlantation)*DBH^(-0.178754 + 0.045881 * isPlantation)), # Sibbesen
                    "RA" ~ 1.37 + (23.50353 + 4.61621 * isPlantation) * (1 - exp(-0.02751*DBH^(1.18019 - 0.06838 * isPlantation))), # Weibull
                    "WH" ~ 1.37 + (79.1071 - 31.3475 * isPlantation) / (1 + (191.3222 - 66.8086 * isPlantation) * DBH^(-1.2302 - 0.1349 * isPlantation)), # Hossfeld
                    "BM" ~ 1.37 + 24.732461 * (1 + ((1.37/24.732461)^(1 - 0.538452) - 1) * exp((-(0.021496 + 0.005265 * isPlantation) * DBH)/0.538452^(0.538452/(1 - 0.538452))))^(1/(1 - 0.538452)), # unified Richards
                    "OM" ~ 1.37 + (0.33786 + 0.26972 * isPlantation)*DBH^(1.78798*DBH^(-0.14864 - 0.05623 * isPlantation)), # Sibbesen
                    "RC" ~ 1.37 + (75.108 - 25.495 * isPlantation) * DBH^1.182 / (226.369 - 101.183 * isPlantation + DBH^1.182), # Michaelis-Menten
                    .default = (0.815747 + 0.250760 * isPlantation) * DBH + (-0.003644 - -0.005235 * isPlantation) * DBH^2)) # parabolic
}

is_model_converged = function(model)
{
  switch(class(model)[1], 
         "gam" = { isConverged = model$converged }, 
         "gamm" = { isConverged = TRUE }, # gamm() appears not to have a convergence flag
         "gnls" = { isConverged = (model$numIter < 500) }, # crude approximation since gnls() doesn't report convergence info
         "gsl_nls" = { isConverged = model$convInfo$isConv }, 
         "lm" = { isConverged = TRUE }, # linear models are deterministic, so no convergence to check
         "nlme" = { isConverged = TRUE }, # nlme() errors if not converged
         "nlrob" = { isConverged = model$status == "converged" },
         stop(paste0("Unhandled model fitting method ", class(model)[1], ".")))
  return(isConverged)
}

plot_auc_bank = function(aucs, fillLabel = "median\nAUC", omitMab = FALSE, xLimits = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"))
{
  if (omitMab)
  {
    aucBank = ggplot(aucs) +
        geom_raster(aes(x = species, y = name, fill = aucMae)) +
        labs(title = "a) MAE", x = NULL, y = NULL, fill = fillLabel) +
        scale_y_discrete(limits = rev)
    titles = c("b) RMSE", "c) AIC", "d) model efficiency")
  } else {
    aucBank = ggplot(aucs) +
        geom_raster(aes(x = species, y = name, fill = aucMab)) +
        labs(title = "a) MAB", x = NULL, y = NULL, fill = fillLabel) +
        scale_y_discrete(limits = rev) +
      ggplot(aucs) +
        geom_raster(aes(x = species, y = name, fill = aucMae)) +
        labs(title = "b) MAE", x = NULL, y = NULL, fill = fillLabel) +
        scale_y_discrete(labels = NULL, limits = rev)
    titles = c("c) RMSE", "d) AIC", "e) model efficiency")
  }

  aucBank = aucBank +
    ggplot(aucs) +
      geom_raster(aes(x = species, y = name, fill = aucRmse)) +
      labs(title = titles[1], x = NULL, y = NULL, fill = fillLabel) +
      scale_y_discrete(labels = NULL, limits = rev) +
    ggplot(aucs) +
      geom_raster(aes(x = species, y = name, fill = aucDeltaAicN)) +
      labs(title = titles[2], x = NULL, y = NULL, fill = fillLabel) +
      scale_y_discrete(labels = NULL, limits = rev) +
    ggplot(aucs) +
      geom_raster(aes(x = species, y = name, fill = aucNse)) +
      labs(title = titles[3], x = NULL, y = NULL, fill = fillLabel) +
      scale_y_discrete(labels = NULL, limits = rev) +
    plot_annotation(theme = theme(plot.margin =  margin())) +
    plot_layout(nrow = 1, guides = "collect") &
      scale_fill_scico(palette = "bam", limits = c(0, 1), na.value = rgb(0.9642, 0.9444, 0.9435)) &
      scale_x_discrete(limits = xLimits) &
      #scale_x_discrete(labels = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), limits = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")) &
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.spacing.y = unit(0.3, "line"), panel.grid = element_blank())
  
  return(aucBank)
}

plot_exploratory = function(liveUnbrokenTrees, plotLetters = c("a)", "b)", "c)"), speciesLabel = NULL, distributionLegendPositionY = 1, maxTreesMeasured = 400, omitLegends = FALSE, omitQuantiles = FALSE, omitXlabels = FALSE)
{
  dbhQuantiles = liveUnbrokenTrees %>% mutate(diameterClass = 2.5 * (ceiling(DBH / 2.5) - 0.5)) %>% group_by(diameterClass) %>%
    reframe(count = n(), quantiles = c("min", "q025", "q10", "q20", "q25", "q30", "q40", "median", "q60", "q70", "q75", "q80", "q90", "q975", "max"), height = quantile(TotalHt, probs = c(0, 0.025, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.975, 1), na.rm = TRUE), mean = mean(TotalHt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = quantiles, values_from = height)
  heightQuantiles = liveUnbrokenTrees %>% mutate(heightClass = 1 * (ceiling(TotalHt / 1) - 0.5)) %>% group_by(heightClass) %>%
    reframe(count = n(), quantiles = c("min", "q025", "q10", "q20", "q25", "q30", "q40", "median", "q60", "q70", "q75", "q80", "q90", "q975", "max"), dbh = quantile(DBH, probs = c(0, 0.025, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.975, 1), na.rm = TRUE), mean = mean(DBH, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = quantiles, values_from = dbh)

  distributionLegendPosition = c(1, distributionLegendPositionY)
  treeLegendPosition = c(1, 0.04)
  if (omitLegends)
  {
    distributionLegendPosition = "none"
    treeLegendPosition = "none"
  }
  
  dbhXlabel = "DBH, cm"
  heightXlabel = "height, m"
  if (omitXlabels)
  {
    dbhXlabel = NULL
    heightXlabel = NULL
  }
  
  heightPower = 1
  dbhPower = 1
  exploratoryPlots = ggplot() +
    geom_bin_2d(aes(x = DBH, y = TotalHt, fill = after_stat(count), weight = TreeCount), liveUnbrokenTrees %>% filter(is.na(TotalHt) == FALSE), binwidth = c(2.5, 1)) +
    geom_path(aes(x = diameterClass, y = mean, color = "mean height", linetype = "mean height"), dbhQuantiles %>% filter(count > 10), na.rm = TRUE) +
    #geom_path(aes(x = diameterClass, y = median, color = "median height", linetype = "median height"), dbhQuantiles %>% filter(count > 10), na.rm = TRUE) +
    geom_path(aes(x = mean, y = heightClass, color = "mean DBH", linetype = "mean DBH"), heightQuantiles %>% filter(count > 10), na.rm = TRUE) +
    #geom_path(aes(x = median, y = heightClass, color = "median DBH", linetype = "median DBH"), heightQuantiles %>% filter(count > 10), na.rm = TRUE) +
    annotate("text", x = 0, y = 81.5, label = paste(plotLetters[1], speciesLabel), hjust = 0, size = 3.5) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 80)) +
    labs(x = dbhXlabel, y = "height, m, of unbroken stem", color = NULL, fill = "trees\nmeasured", linetype = NULL) +
    guides(color = guide_legend(order = 1), fill = guide_colorbar(order = 2), linetype = guide_legend(order = 1)) +
    scale_color_manual(breaks = c("mean height", "median height", "mean DBH", "median DBH"), labels = c("mean\nheight", "median\nheight", "mean\nDBH", "median\nDBH"), values = c("green2", "green2", "burlywood2", "burlywood2")) +
    scale_fill_viridis_c(breaks = c(1, 3, 10, 33, 100, 330), limits = c(1, maxTreesMeasured), trans = "log10") +
    scale_linetype_manual(breaks = c("mean height", "median height", "mean DBH", "median DBH"), labels = c("mean\nheight", "median\nheight", "mean\nDBH", "median\nDBH"), values = c("solid", "longdash", "solid", "longdash")) +
    theme(legend.key.height = unit(0.95, "line"), legend.justification = c(1, 0), legend.position = treeLegendPosition, legend.title = element_text(size = 9.5), legend.spacing.y = unit(0.25, "line"))
  if (omitQuantiles)
  {
    return(exploratoryPlots)
  }
  
  exploratoryPlots = exploratoryPlots +
    ggplot(dbhQuantiles) +
      geom_ribbon(aes(x = diameterClass, ymin = 100 * (q025 - mean) / mean^heightPower, ymax = 100 * (q975 - mean) / mean^heightPower, alpha = "95% probability"), fill = "forestgreen") +
      geom_ribbon(aes(x = diameterClass, ymin = 100 * (q10 - mean) / mean^heightPower, ymax = 100 * (q90 - mean) / mean^heightPower, alpha = "80% probability"), fill = "forestgreen") +
      geom_ribbon(aes(x = diameterClass, ymin = 100 * (q25 - mean) / mean^heightPower, ymax = 100 * (q75 - mean) / mean^heightPower, alpha = "50% probability"), fill = "forestgreen") +
      geom_path(aes(x = diameterClass, y = 100 * (min - mean) / mean^heightPower, color = "max or min", linetype = "max or min"), na.rm = TRUE, linewidth = 0.3) +
      #geom_path(aes(x = diameterClass, y = 100 * (q10 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (q20 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (q30 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (q40 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_segment(x = 0, xend = 185, y = 0, yend = 0, color = "forestgreen", linewidth = 0.4) +
      geom_path(aes(x = diameterClass, y = 100 * (q60 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (q70 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (q80 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      #geom_path(aes(x = diameterClass, y = 100 * (q90 - mean) / mean^heightPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = diameterClass, y = 100 * (max - mean) / mean^heightPower, color = "max or min", linetype = "max or min"), na.rm = TRUE, linewidth = 0.3) +
      annotate("text", x = 0, y = 154, label = paste(plotLetters[2], speciesLabel), hjust = 0, size = 3.5) +
      coord_cartesian(xlim = c(0, 196), ylim = c(-50, 150)) +
      scale_alpha_manual(breaks = c("95% probability", "80% probability", "50% probability"), values = c(0.1, 0.2, 0.3)) +
      scale_color_manual(breaks = c("10% contour", "max or min"), values = c("grey50", "grey70")) +
      scale_linetype_manual(breaks = c("10% contour", "max or min"), values = c("dashed", "dotted")) +
      labs(x = dbhXlabel, y = "departure from mean height, %", alpha = NULL, color = NULL, linetype = NULL) +
      theme(legend.position = "none") +
    ggplot(heightQuantiles) +
      geom_ribbon(aes(x = heightClass, ymin = 100 * (q025 - mean) / mean^dbhPower, ymax = 100 * (q975 - mean) / mean^dbhPower, alpha = "95% probability"), fill = "burlywood4") +
      geom_ribbon(aes(x = heightClass, ymin = 100 * (q10 - mean) / mean^dbhPower, ymax = 100 * (q90 - mean) / mean^dbhPower, alpha = "80% probability"), fill = "burlywood4") +
      geom_ribbon(aes(x = heightClass, ymin = 100 * (q25 - mean) / mean^dbhPower, ymax = 100 * (q75 - mean) / mean^dbhPower, alpha = "50% probability"), fill = "burlywood4") +
      geom_path(aes(x = heightClass, y = 100 * (min - mean) / mean^dbhPower, color = "max or min", linetype = "max or min"), na.rm = TRUE, linewidth = 0.3) +
      #geom_path(aes(x = heightClass, y = 100 * (q10 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (q20 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (q30 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (q40 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_segment(x = 0, xend = 77.5, y = 0, yend = 0, color = "burlywood4", linewidth = 0.4) +
      geom_path(aes(x = heightClass, y = 100 * (q60 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (q70 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (q80 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      #geom_path(aes(x = heightClass, y = 100 * (q90 - mean) / mean^dbhPower, color = "10% contour", linetype = "10% contour"), na.rm = TRUE, linewidth = 0.3) +
      geom_path(aes(x = heightClass, y = 100 * (max - mean) / mean^dbhPower, color = "max or min", linetype = "max or min"), na.rm = TRUE, linewidth = 0.3) +
      annotate("text", x = 0, y = 154, label = paste(plotLetters[3], speciesLabel), hjust = 0, size = 3.5) +
      coord_cartesian(xlim = c(0, 80), ylim = c(-50, 150)) +
      guides(alpha = guide_legend(order = 1, override.aes = list(fill = "grey30")), color = guide_legend(order = 2), linetype = guide_legend(order = 2)) +
      scale_alpha_manual(breaks = c("95% probability", "80% probability", "50% probability"), values = c(0.1, 0.2, 0.3)) +
      scale_color_manual(breaks = c("10% contour", "max or min"), values = c("grey50", "grey70")) +
      scale_linetype_manual(breaks = c("10% contour", "max or min"), values = c("dashed", "dotted")) +
      labs(x = heightXlabel, y = "departure from mean DBH, %", alpha = NULL, color = NULL, linetype = NULL) +
      theme(legend.justification = c(1, 1), legend.position = distributionLegendPosition) +
    plot_layout(nrow = 1, ncol = 3, widths = c(260, 200, 200))
  return(exploratoryPlots)
}

plot_qq = function(diameterRegression1, diameterRegression2, diameterRegression3, diameterRegression4,
                   heightRegression1, heightRegression2, heightRegression3, heightRegression4,
                   speciesName, tDegreesOfFreedom = 7, tSkew = 2.25)
{
  #heightColors = c("#ff6961", "#ffb480", "#f8f38d", "#42d6a4")
  #dbhColors = c("#08cad1", "#59adf6", "#9d94ff", "#c780e8")
  heightColors = viridis::viridis_pal(option = "plasma", end = 0.9)(4)
  dbhColors = viridis::viridis_pal(end = 0.9)(4)
  qqPlot = ggplot() +
      geom_qq_line(aes(sample = -residuals(diameterRegression1), color = diameterRegression1$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(diameterRegression2), color = diameterRegression2$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(diameterRegression3), color = diameterRegression3$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(diameterRegression4), color = diameterRegression4$name), alpha = 0.4) +
      geom_qq(aes(sample = -residuals(diameterRegression1), color = diameterRegression1$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression2), color = diameterRegression2$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression3), color = diameterRegression3$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression4), color = diameterRegression4$name), alpha = 0.8, geom = "line") +
      annotate("text", x = -10.5, y = 160, label = paste0("'a) ", speciesName, " height, '*epsilon~'~'~'N(0, '*sigma*'²)'"), hjust = 0, parse = TRUE, size = 3.4) +
      coord_cartesian(xlim = c(-10, 13), ylim = c(-110, 160)) +
      labs(x = NULL, y = "sample quantile", color = NULL) +
      scale_color_manual(values = heightColors) +
      theme(legend.key.height = unit(0.8, "line"), legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
    ggplot() +
      geom_qq_line(aes(sample = -residuals(heightRegression1), color = heightRegression1$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(heightRegression2), color = heightRegression2$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(heightRegression3), color = heightRegression3$name), alpha = 0.4) +
      geom_qq_line(aes(sample = -residuals(heightRegression4), color = heightRegression4$name), alpha = 0.4) +
      geom_qq(aes(sample = -residuals(heightRegression1), color = heightRegression1$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression2), color = heightRegression2$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression3), color = heightRegression3$name), alpha = 0.8, geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression4), color = heightRegression4$name), alpha = 0.8, geom = "line") +
      annotate("text", x = -10.5, y = 160, label = paste0("'b) ", speciesName, " DBH, '*epsilon~'~'~'N(0, '*sigma*'²)'"), hjust = 0, parse = TRUE, size = 3.4) +
      coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
      labs(x = NULL, y = NULL, color = NULL) +
      scale_color_manual(values = dbhColors) +
      theme(legend.key.height = unit(0.8, "line"),legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
    ggplot() +
      geom_qq_line(aes(sample = -residuals(diameterRegression1), color = diameterRegression1$name), alpha = 0.4, distribution = qt, dparams = list(df = tDegreesOfFreedom)) +
      geom_qq_line(aes(sample = -residuals(diameterRegression2), color = diameterRegression2$name), alpha = 0.4, distribution = qt, dparams = list(df = tDegreesOfFreedom)) +
      geom_qq_line(aes(sample = -residuals(diameterRegression3), color = diameterRegression3$name), alpha = 0.4, distribution = qt, dparams = list(df = tDegreesOfFreedom)) +
      geom_qq_line(aes(sample = -residuals(diameterRegression4), color = diameterRegression4$name), alpha = 0.4, distribution = qt, dparams = list(df = tDegreesOfFreedom)) +
      geom_qq(aes(sample = -residuals(diameterRegression1), color = diameterRegression1$name), alpha = 0.8, distribution = qt, dparams = list(df = tDegreesOfFreedom), geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression2), color = diameterRegression2$name), alpha = 0.8, distribution = qt, dparams = list(df = tDegreesOfFreedom), geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression3), color = diameterRegression3$name), alpha = 0.8, distribution = qt, dparams = list(df = tDegreesOfFreedom), geom = "line") +
      geom_qq(aes(sample = -residuals(diameterRegression4), color = diameterRegression4$name), alpha = 0.8, distribution = qt, dparams = list(df = tDegreesOfFreedom), geom = "line") +
      annotate("text", x = -10.5, y = 160, label = paste0("'c) ", speciesName, " height, '*epsilon~'~'~'t(df = ", tDegreesOfFreedom, ")'"), hjust = 0, parse = TRUE, size = 3.4) +
      coord_cartesian(xlim = c(-10, 13), ylim = c(-110, 160)) +
      labs(x = "theoretical quantile", y = "sample quantile", color = NULL) +
      scale_color_manual(values = heightColors) +
      theme(legend.justification = c(1, 0), legend.position = "none") +
    ggplot() + # qst()'s omega (scale) parameter can be left as 1 as its only effect is rotation, xi (location) can be left as zero as its only effect is a translation in theoretical quantile
      geom_qq_line(aes(sample = -residuals(heightRegression1), color = heightRegression1$name), alpha = 0.4, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0)) +
      geom_qq_line(aes(sample = -residuals(heightRegression2), color = heightRegression2$name), alpha = 0.4, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0)) +
      geom_qq_line(aes(sample = -residuals(heightRegression3), color = heightRegression3$name), alpha = 0.4, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0)) +
      geom_qq_line(aes(sample = -residuals(heightRegression4), color = heightRegression4$name), alpha = 0.4, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0)) +
      geom_qq(aes(sample = -residuals(heightRegression1), color = heightRegression1$name), alpha = 0.8, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0), geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression2), color = heightRegression2$name), alpha = 0.8, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0), geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression3), color = heightRegression3$name), alpha = 0.8, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0), geom = "line") +
      geom_qq(aes(sample = -residuals(heightRegression4), color = heightRegression4$name), alpha = 0.8, distribution = sn::qst, dparams = list(nu = tDegreesOfFreedom, alpha = tSkew, omega = 1, xi = 0), geom = "line") +
      annotate("text", x = -10.5, y = 160, label = paste0("'d) ", speciesName, " DBH, '*epsilon~'~'~'t(df = ", tDegreesOfFreedom, ", '*alpha*' = ", tSkew, ")'"), hjust = 0, parse = TRUE, size = 3.4) +
      coord_cartesian(xlim = c(-10, 16.5), ylim = c(-110, 160)) +
      labs(x = "theoretical quantile", y = NULL, color = NULL) +
      scale_color_manual(values = dbhColors) +
      theme(legend.justification = c(1, 0), legend.position = "none") +
      plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
      plot_layout(nrow = 2, ncol = 2, widths = c(10 + 13, 10 + 16.5))
  return(qqPlot)
}


## load data
# Notable properties of cruise data loaded into trees2016
#  - Plots are either count plots, where trees aren't measured, or count plots, where all trees and snags are measured
#    for DBH and a subset measured for height. Thus, all stems with heights (TotalHt if unbroke, Ht2 if broken) are 
#    also have DBH measurements.
#  - Trees are on CO (count) and IP (measure) plots. Count plots are variable radius and count trees by species. 
#    Measure plots are nested variable radius and small tree fixed plots. Small tree plots count trees trees up to four
#    inches DBH (10 cm) by species with DBH to the nearest inch (excepting two records). CB and IB plots contain only 
#    records for other species with TreeCount = 0.
#  - Trees on variable radius measure plots generally have TreeCount = 1, as expected, but 462 records have TreeCount
#    = 2. Since it's very unlikely two trees are on the same plot with the same DBH and, often, the same height, these
#    records are assumed to be incorrect and the tree count is changed to one.
stands2022 = read_xlsx("GIS/Planning/Elliott Stand Data Feb2022.xlsx") %>% 
  mutate(Cruised_Si = na_if(Cruised_Si, 0),
         ODSL_Site_ = na_if(ODSL_Site_, 0),
         siteSpecies = if_else(startsWith(ODSL_VEG_L, "1W") | startsWith(ODSL_VEG_L, "WX"), "hemlock", 
                               if_else(startsWith(ODSL_VEG_L, "1H") | startsWith(ODSL_VEG_L, "HX"), "hardwood",
                                       if_else(startsWith(ODSL_VEG_L, "OT"), "other",
                                               "Douglas-fir"))),
         standAge2016 = pmax(if_else((Age_2020 - 4) > (Age_2015 + 1), Age_2015 + 1, Age_2020 - 4), 0),
         standArea = 0.404686 * GrossAc,  # ac to ha
         isPlantation = standAge2016 < 70)

plots2016 = read_xlsx("GIS/Trees/2015-16 cruise/CruisePlots_All_20151211.xlsx") # both 20151211 and 20160111 missing coordinates for 171 plots in stands 1661 and 2470

trees2016 = left_join(left_join(read_xlsx("trees/Elliott final cruise records 2015-16.xlsx", sheet = "CRUISERECS"),
                                stands2022 %>% select(StandID, standAge2016, standArea, isPlantation),
                                by = c("StandID")),
                      plots2016 %>% select(STAND, PltInteger, elevation, slope, aspect, topographicShelterIndex, x, y) %>% rename(PlotID = PltInteger),
                      by = c("PlotID")) %>%
  mutate(speciesGroup = factor(if_else(Species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), Species, "other"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other")),
         BHAge = na_if(BHAge, 0), # years
         DBH = na_if(2.54 * DBH, 0), # inches to cm
         Dia1 = na_if(2.54 * Dia1, 0),
         CrownRatio = na_if(CrownRatio, 0),
         Ht1 = na_if(0.3048 * Ht1, 0), # feet to m
         Ht2 = na_if(0.3048 * Ht2, 0),
         isConifer = Species %in% c("DF", "WH", "RC", "SS", "CX", "PC", "PY", "GF", "LP"),
         isLive = (CompCode %in% c("D.", "SN")) == FALSE,
         isLiveUnbroken = isLive & (CompCode != "BT"),
         plotRadius = if_else(SamplingMethod == "BAF", 100 / 2.54 * 0.3048 / (12 * sqrt(SampleFactor / 10890)), 0.3048 * sqrt(43560 / (pi * SampleFactor))), # m
         SampleFactor = 2.47105 * if_else(SamplingMethod == "BAF", 0.092903, 1) * SampleFactor, # convert BAF from ft²/ac to m²/ha and TPA to TPH, BAF conversion is BAF ft²/ac * 2.47105 ac/ha * 0.092903 m²/ft² = 0.229568 m²/ha / ft²/ac
         TotalHt = na_if(0.3048 * TotalHt, 0),
         TreeCount = if_else((PlotType == "IP") & (SamplingMethod == "BAF") & (TreeCount > 1), 1, TreeCount), # fix tree duplication per notes above
         basalArea = 0.25 * pi * (0.01*DBH)^2, # m² 
         breastHeight = 1.37, # m, used for offset in lm() height regressions
         heightDiameterRatio = TotalHt / (0.01 * DBH), # (DBH conversion from cm to m)
         imputedBasalArea = if_else(is.na(TotalHt) == FALSE, impute_basal_area(Species, TotalHt, isPlantation), basalArea), # m², imputed whenever height is available to impute
         imputedHeight = if_else(is.na(TotalHt) == FALSE, TotalHt, if_else(is.na(DBH) == FALSE, impute_height(Species, DBH, isPlantation), NA_real_)), # where possible, perform basic height imputation
         treeBasalAreaPerHectare = SampleFactor * TreeCount * if_else(SamplingMethod == "BAF", 1, basalArea), # m²/ha, measure plots have TreeCount = 1 for each tree, count plots have TreeCount = 0-41 depending on the number of trees present
         treeBasalAreaPerHectareApprox = SampleFactor * TreeCount * if_else(SamplingMethod == "BAF", if_else(is.na(basalArea) == FALSE, imputedBasalArea / basalArea, 1), imputedBasalArea)) %>% # m²/ha, stack basal area regression on height regression when possible; BAF has to be used with prism trees but imputed basal area is used where possible to introduce estimation error
  group_by(StandID) %>%
  arrange(desc(isLiveUnbroken), desc(DBH), .by_group = TRUE) %>% # put largest diameter live trees first in each stand for calculating BAL (numbers sort before NA)
  mutate(plotsInStand = length(unique(PlotID)), # nested fixed radius and BAF plots share same plot ID
         standBasalAreaPerHectare = sum(isLive * treeBasalAreaPerHectare) / plotsInStand, # m²/ha
         standBasalAreaApprox = sum(isLive * treeBasalAreaPerHectareApprox) / plotsInStand, # m²/ha
         basalAreaLarger = (cumsum(isLive * treeBasalAreaPerHectare) - treeBasalAreaPerHectare[1]) / plotsInStand, # m²/ha
         measurePlotsInStand = length(unique(PlotID * (PlotType == "IP"))) - any(PlotType %in% c("IB", "CO", "CB")), # if any IB, CO, or CB plots are present they'll introduce zero as a unique value which has to be subtracted from the plot count
         measureTreeTphContribution = if_else((TreeCount == 0) | is.na(DBH), NA_real_, SampleFactor * TreeCount * if_else(SamplingMethod == "BAF",  1 / basalArea, 1)), # trees per hectare, with trees not on measure plots NAed out
         meanTreesPerBafPlot = sum(TreeCount * (SamplingMethod == "BAF")) / plotsInStand,
         meanTreesPerBafMeasurePlot = sum(TreeCount * (SamplingMethod == "BAF") * (PlotType == "IP")) / measurePlotsInStand,
         tph = meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * sum(isLive * measureTreeTphContribution, na.rm = TRUE) / measurePlotsInStand, # stand's total trees per hectare
         qmd = sqrt(standBasalAreaPerHectare / (pi/4 * 0.01^2 * tph)), # quadratic mean diameter, cm
         relativeDiameter = DBH / qmd) %>%
  # top height by tallest trees in stand, regardless of plot
  #arrange(desc(isLiveUnbroken), desc(TotalHt), .by_group = TRUE) %>% # put tallest live trees without broken tops first in each stand
  #mutate(topHeightTph = pmin(cumsum(if_else(is.na(TotalHt), 0, measureTreeTphContribution)), 100), # TPH total towards the H100 definition of top height, trees not measured for TotalHt are skipped
  #       topHeightWeight = pmax((topHeightTph - lag(topHeightTph, default = 0)) / measureTreeTphContribution, 0), # clamp remaining fraction to [0, 1] to get individual trees' contributions to the top height average
  #       topHeight = sum(topHeightWeight * TotalHt, na.rm = TRUE) / sum(topHeightWeight, na.rm = TRUE), # m, tallest 100 trees per hectare
  #       relativeHeight = TotalHt / topHeight, # individual trees' heights as a fraction of top height, may be greater than 1, especially for retention trees (debatable if imputed heights should be included but, for now, trees not measured for height are left with NA relative height)
  #       tallerApproxBasalArea = (cumsum(isLive * treeBasalAreaPerHectareApprox) - treeBasalAreaPerHectareApprox[1]) / plotsInStand,
  #       tallerTph = cumsum(isLiveUnbroken * SampleFactor * TreeCount * if_else(SamplingMethod == "BAF",  1 / basalArea, 1)) / plotsInStand) %>% 
  # top height by estimating H100 on each plot and then averaging all plots (slower than pooling by stand)
  group_by(StandID, PlotID, isLiveUnbroken) %>%
  arrange(desc(TotalHt), .by_group = TRUE) %>% 
  mutate(topHeightTph = if_else(isLiveUnbroken, pmin(cumsum(if_else(is.na(TotalHt), 0, measureTreeTphContribution)), 100), NA_real_),
         topHeightWeight = if_else(isLiveUnbroken, pmax((topHeightTph - lag(topHeightTph, default = 0)) / measureTreeTphContribution, 0), NA_real_),
         topHeight = if_else(isLiveUnbroken, sum(topHeightWeight * TotalHt, na.rm = TRUE) / sum(topHeightWeight, na.rm = TRUE), NA_real_),
         topHeightMask = if_else(isLiveUnbroken & (is.na(topHeightWeight) == FALSE), row_number() == n(), NA_real_)) %>% # topHeightWeight is NA for count plots
  group_by(StandID) %>%
  arrange(desc(isLiveUnbroken), desc(TotalHt), .by_group = TRUE) %>%
  mutate(topHeight = sum(topHeightMask * topHeightTph * topHeight, na.rm = TRUE) / sum(topHeightMask * topHeightTph, na.rm = TRUE),
         relativeHeight = TotalHt / topHeight,
         tallerApproxBasalArea = (cumsum(isLive * treeBasalAreaPerHectareApprox) - treeBasalAreaPerHectareApprox[1]) / plotsInStand,
         tallerTph = cumsum(isLiveUnbroken * SampleFactor * TreeCount * if_else(SamplingMethod == "BAF",  1 / basalArea, 1)) / plotsInStand) %>%
  ungroup()

heightClassBreaks = trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE) %>%
  group_by(speciesGroup) %>%
  group_modify(~{
                  quantileBreaks = seq(0, 1, length.out = min(50, sum(.$TreeCount) / (5 * 10))) # constrain maximum number of classes based on data availability: setting the max to n / (meanClassN*k) classes averages meanClassN samples per class in validation folds => primarily affects low n species: Oregon myrtle, western redcedar, and other
                  return(tibble(heightBreaks = unique(ceiling(c(0, quantile(.$TotalHt, probs = quantileBreaks, na.rm = TRUE))))))
                }) %>%
  unstack(heightBreaks ~ speciesGroup) # list of height class breaks, named by species
dbhClassBreaks = trees2016 %>% filter(isLiveUnbroken, DBH > 2.54 * 3.5) %>%
  group_by(speciesGroup) %>%
  group_modify(~{
                  quantileBreaks = seq(0, 1, length.out = min(50, sum(.$TreeCount) / (5 * 10) - 3))
                  return(tibble(dbhBreaks = unique(c(2.5 * c(0, 1.5, 2.5, 3.5), 2.5 * ceiling(quantile(.$DBH, probs = quantileBreaks, na.rm = TRUE) / 2.5) + 0.5 * 2.5))))
                }) %>%
  unstack(dbhBreaks ~ speciesGroup) # list of DBH class breaks, named by species

trees2016 %<>% group_by(speciesGroup) %>%
  mutate(heightClass = cut(TotalHt, breaks = heightClassBreaks[[cur_group()$speciesGroup]], labels = 0.5 * (head(heightClassBreaks[[cur_group()$speciesGroup]], -1) + tail(heightClassBreaks[[cur_group()$speciesGroup]], -1))),
         dbhClass = cut(DBH, breaks = dbhClassBreaks[[cur_group()$speciesGroup]], labels = 0.5 * (head(dbhClassBreaks[[cur_group()$speciesGroup]], -1) + tail(dbhClassBreaks[[cur_group()$speciesGroup]], -1)))) %>%
  ungroup()

if (htDiaOptions$includeInvestigatory)
{
  # plots without spatial locations
  print(trees2016 %>% filter(is.na(elevation)) %>% group_by(PlotID) %>% summarize(trees = n(), .groups = "drop"), n = 51)
  # distribution of estimated stand basal areas
  ggplot(trees2016) + geom_histogram(aes(x = standBasalAreaApprox))
  # tree contribution to top height verification
  #print(tibble(treeTphContribution = rep(25.5, 30) / 5) %>% 
  #        mutate(remainingTph = 100 - cumsum(treeTphContribution),
  #               remainingFraction = (treeTphContribution + remainingTph) / treeTphContribution,
  #               weight = if_else(remainingFraction >= 1, 1, if_else(remainingFraction > 0, remainingFraction, 0))),
  #      n = 40)

  # check plots for calculated stand-level quantities: BA, TPH, QMD, H100
  standsFromTrees2016 = trees2016 %>% group_by(StandID) %>%
    summarize(plots = plotsInStand[1], measurePlots = measurePlotsInStand[1], meanTreesPerBafPlot = meanTreesPerBafPlot[1], meanTreesPerBafMeasurePlot = meanTreesPerBafMeasurePlot[1],
              tph = tph[1], topHeight = topHeight[1], standBasalAreaPerHectare = standBasalAreaPerHectare[1], standBasalAreaApprox = standBasalAreaApprox[1],
              qmd = sqrt(standBasalAreaPerHectare / (pi / (4 * 100^2) * tph)))
  standsFromTrees2016 %>% summarise(stands = n(), plotsOK = sum(measurePlots <= plots), trees = sum(is.na(meanTreesPerBafPlot) == FALSE), measureTrees = sum(is.na(meanTreesPerBafMeasurePlot) == FALSE),
                                    topHeight = sum(is.na(topHeight) == FALSE), basalArea = sum(is.na(standBasalAreaPerHectare) == FALSE), basalAreaApprox = sum(is.na(standBasalAreaApprox) == FALSE))
  reinekeSdi = crossing(sdi = c(100 * seq(1, 9), 1000 * seq(1, 10)),
                        tph = c(1, 10000)) %>%
    mutate(qmd = 25.4 * (sdi/tph)^(1/1.605)) # sdi = tph * (qmd/25.4)^1.605 => (sdi/tph)^(1/1.605) = (qmd/25.4)
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = plots, y = measurePlots), standsFromTrees2016, alpha = 0.2, color = "grey25", shape = 16) +
    coord_cartesian(xlim = c(0, 90)) +
    labs(x = "plots", y = "measure plots") +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 120, yend = 120), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = standBasalAreaPerHectare, y = standBasalAreaApprox), standsFromTrees2016, alpha = 0.2, color = "grey25", shape = 16) +
    coord_cartesian(xlim = c(0, 120)) +
    labs(x = bquote("basal area, m"^2*" ha"^-1), y = bquote("approximate basal area, m"^2*" ha"^-1)) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 10, yend = 10), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = meanTreesPerBafPlot, y = meanTreesPerBafMeasurePlot), standsFromTrees2016, alpha = 0.2, color = "grey25", shape = 16) +
    #coord_cartesian(xlim = c(0, 120)) +
    labs(x = "mean trees per BAF plot", y = "mean trees per BAF measure plot") +
  ggplot() +
    geom_line(aes(x = tph, y = qmd, group = sdi), reinekeSdi, color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = tph, y = qmd), standsFromTrees2016, alpha = 0.5, shape = 16) +
    coord_cartesian(xlim = c(50, 5000), ylim = c(5, 95)) +
    labs(x = "TPH", y = "QMD, cm") +
    scale_x_log10(breaks = c(50, 100, 200, 500, 1000, 5000), minor_breaks = c(60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 2000, 3000, 4000, 6000, 7000)) +
    scale_y_log10(breaks = c(5, 10, 20, 50, 100), minor_breaks = c(6, 7, 8, 9, 30, 40, 60, 70, 80, 90)) +
  ggplot() +
    geom_histogram(aes(y = topHeight), standsFromTrees2016, binwidth = 2) +
    labs(x = "stands", y = bquote("H"[100]*", m")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 3)

  # check plots for tree-level properties derived from stand-level properties
  ggplot() +
    geom_segment(aes(x = 1.5, y = 0, xend = 1.5, yend = 3000), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_histogram(aes(x = relativeHeight, fill = speciesGroup), trees2016, binwidth = 0.05, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, 2500)) +
    labs(x = "relative height", y = "trees measured", fill = NULL) +
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
  ggplot() +
    geom_segment(aes(x = 5.4, y = 0, xend = 5.4, yend = 3000), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_histogram(aes(x = relativeDiameter, fill = speciesGroup), trees2016, binwidth = 0.05, na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 7.5), ylim = c(0, 2500)) +
    labs(x = "relative diameter", y = "trees measured", fill = NULL) +
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    theme(legend.position = "none") +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 125, yend = 125), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = standBasalAreaPerHectare, y = basalAreaLarger, color = speciesGroup), trees2016, alpha = 0.2, shape = 16) +
    guides(color = "none") +
    labs(x = bquote("stand basal area, m"^2*" ha"^-1), y = bquote("basal area larger, m"^2*" ha"^-1), color = NULL) +
    scale_color_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 180, yend = 180), color = "grey80", linewidth = 0.3, linetype = "longdash") +
    geom_point(aes(x = standBasalAreaApprox, y = tallerApproxBasalArea, color = speciesGroup), trees2016, alpha = 0.2, shape = 16) +
    guides(color = "none") +
    labs(x = bquote("approximate stand basal area, m"^2*" ha"^-1), y = bquote("basal area taller, m"^2*" ha"^-1), color = NULL) +
    scale_color_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 2)
  
  # correlations among predictors
  predictorLabels = c("DBH", "height", "height:diameter", "stand age", "BA", "BAL", "ABA", "AAT", "elevation", "slope", "sin(aspect)", "cos(aspect)", "TSI", "H100", "RelHt", "QMD", "RelDbh")
  predictorLevels = c("DBH", "TotalHt", "heightDiameterRatio", "standAge", "standBasalAreaPerHectare", "basalAreaLarger", "standBasalAreaApprox", "tallerApproxBasalArea", "elevation", "slope", "sinAspect", "cosAspect", "topographicShelterIndex", "topHeight", "relativeHeight", "qmd", "relativeDiameter")
  predictorCorrelation = trees2016 %>% filter(isLiveUnbroken) %>% 
    mutate(heightDiameterRatio = TotalHt / (0.01 * DBH), sinAspect = sin(pi/180 * aspect), cosAspect = cos(pi/180 * aspect)) %>% 
    group_by(speciesGroup) %>%
    group_modify(~as_tibble(cor(as.matrix(.x %>% mutate(standAge = standAge2020 - 4) %>% select(DBH, TotalHt, heightDiameterRatio, standAge, standBasalAreaPerHectare, basalAreaLarger, standBasalAreaApprox, tallerApproxBasalArea, elevation, slope, sinAspect, cosAspect, topographicShelterIndex, topHeight, relativeHeight, qmd, relativeDiameter) %>% drop_na())), rownames = "predictor1")) %>%
    pivot_longer(-c(speciesGroup, predictor1), names_to = "predictor2", values_to = "correlation") %>%
    mutate(predictor1 = factor(predictor1, labels = predictorLabels, levels = predictorLevels),
           predictor2 = factor(predictor2, labels = predictorLabels, levels = predictorLevels))
  ggplot(predictorCorrelation %>% filter(predictor1 %in% c("DBH", "height", "height:diameter"))) +
    geom_raster(aes(x = predictor1, y = predictor2, fill = correlation)) +
    labs(x = NULL, y = NULL, fill = "correlation") +
    facet_wrap(vars(speciesGroup), nrow = 2, ncol = 4) +
    scale_fill_scico(palette = "vik", limits = c(-1, 1)) +
    scale_y_discrete(limits = rev) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.spacing.y = unit(0.5, "line"), strip.background = element_rect(fill = "grey95"))
  print(predictorCorrelation %>% filter(predictor1 %in% c("DBH", "height", "height:diameter"), predictor2 %in% c("DBH", "height", "height:diameter")) %>%
    pivot_wider(names_from = "predictor2", values_from = "correlation"), n = 21)
  
  # number of multitree height-diameter measurements by species
  trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE) %>% group_by(speciesGroup) %>%
    summarize(single = sum(TreeCount == 1), double = sum(TreeCount == 2), triple = sum(TreeCount == 3), multiple = sum(TreeCount > 3)) %>%
    mutate(pctSingle = 100 * single / (single + double + triple + multiple),
           pctDouble = 100 * double / (single + double + triple + multiple),
           pctTriple = 100 * triple / (single + double + triple + multiple),
           pctMulti = 100 * multiple / (single + double + triple + multiple))
  
  # export tree data joined with plots in R for joining in GIS
  plotTreeProperties = trees2016 %>% group_by(PlotID) %>%
    summarize(liveTrees = sum(isLive), 
              snags = sum(isLive == FALSE), 
              tph = sum(SampleFactor * isLive * TreeCount), 
              primarySpecies = unique(Species)[which.max(tabulate(match(Species, unique(Species))))], # mode of species
              stemsWithDbh = sum(is.na(DBH) == FALSE), 
              stemsWithHeight = sum((is.na(TotalHt) == FALSE) | (is.na(Ht2) == FALSE))) %>%
    mutate(stems = liveTrees + snags) %>%
    rename(PltInteger = PlotID) # for GIS joins
  plotTreeProperties %>% summarize(plots = n(), measure = sum(stemsWithDbh > 0), count = n() - measure)
  #writer::write_csv(plotTreeProperties, "GIS/Trees/2015-16 cruise/CruisePlots_All_20151211 treeProperties.csv")
  
  # export stand properties, including 2016 inventory where available
  # trees.R needs areas for all stands, so join cruised stands with all stands defined for 2016 inventory.
  write_xlsx(list(stands = left_join(stands2022 %>% mutate(StandID = as.integer(StandID)) %>% select(StandID, standAge2016, standArea, isPlantation),
                                     trees2016 %>% group_by(StandID) %>%
                                       summarize(standAge2016 = standAge2016[1], 
                                                 isPlantation = isPlantation[1], 
                                                 standArea = standArea[1], 
                                                 plotsInStand = plotsInStand[1], 
                                                 measurePlotsInStand = measurePlotsInStand[1],
                                                 tph = tph[1], 
                                                 topHeight = topHeight[1], 
                                                 qmd = qmd[1], 
                                                 standBasalAreaPerHectare = standBasalAreaPerHectare[1], 
                                                 standBasalAreaApprox = standBasalAreaApprox[1]) %>%
                                       mutate(StandID = as.integer(StandID)),
                                     by = "StandID") %>%
                            mutate(standAge2016 = if_else(is.na(standAge2016.x), standAge2016.y, standAge2016.x),
                                   standArea = if_else(is.na(standArea.x), standArea.y, standArea.x),
                                   isPlantation = if_else(is.na(isPlantation.x), isPlantation.y, isPlantation.x)) %>%
                            select(-standArea.x, -standArea.y, -standAge2016.x, -standAge2016.y, -isPlantation.x, -isPlantation.y) %>%
                            rename(standID2016 = StandID) %>%
                            relocate(standID2016, standAge2016, standArea,	isPlantation)),
             "GIS/Trees/2015-16 cruise.xlsx")
}


## data tabulation and basic plotting
if (htDiaOptions$includeInvestigatory)
{
  # Table 1
  trees2016summary = trees2016 %>%
    #group_by(Species) %>%
    group_by(speciesGroup) %>% 
    summarize(stands = n_distinct(StandID),
              stems = sum(TreeCount),
              live = sum(TreeCount * isLive),
              plantation = sum(TreeCount * isLive * isPlantation), 
              retention = sum(TreeCount * isLive * (CompCode == "RT")), 
              dbh = sum(TreeCount * isLive * (is.na(DBH) == FALSE), na.rm = TRUE), 
              height = sum(TreeCount * isLive * (is.na(TotalHt) == FALSE), na.rm = TRUE), 
              brokenTop = sum(TreeCount * isLive * (is.na(Ht2) == FALSE), na.rm = TRUE), 
              snag = sum(TreeCount * (isLive == FALSE)), 
              snagDbh = sum(TreeCount * (isLive == FALSE) * (is.na(DBH) == FALSE)), 
              age = sum(TreeCount * (is.na(BHAge) == FALSE), na.rm = TRUE), 
              crownRatio = sum(TreeCount * (is.na(CrownRatio) == FALSE), na.rm = TRUE), 
              taperDia = sum(TreeCount * (is.na(Dia1) == FALSE), na.rm = TRUE), 
              taperHt = sum(TreeCount * (is.na(Ht1) == FALSE), na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(pctStems = 100 * stems / sum(stems), 
           pctStands = 100 * stands / n_distinct(trees2016$StandID)) %>%
    relocate(speciesGroup, stands, pctStems) %>%
    #relocate(Species, stands, pctStems) %>%
    arrange(desc(stems)) %>%
    bind_rows(summarize(., across(where(is.numeric), ~if_else(is.integer(.x), max(.x), sum(.x))))) %>%
    mutate(speciesGroup = replace_na(as.character(speciesGroup), "total"))
    #mutate(Species = replace_na(as.character(Species), "total"))
  print(trees2016summary, n = 25)
  trees2016 %>% group_by(isConifer) %>% summarize(heightMeasureTrees = sum(TreeCount * isLive * (is.na(TotalHt) == FALSE)))

  # plot data summary
  trees2016 %>%
    group_by(PlotID) %>%
    mutate(isLive = CompCode %in% c("D.", "SN") == FALSE,
           isMeasurePlot = sum(is.na(DBH) == FALSE) > 0) %>% # exact check used is unimportant as all stems have DBH on measure plots
    ungroup() %>%
    summarize(standsSampled = n_distinct(StandID),
              plots = n_distinct(PlotID),
              measurePlots = n_distinct(isMeasurePlot * PlotID) - 1, # minus one since zero indicates count plot (plot IDs start with 1)
              countPlots = plots - measurePlots,
              liveTrees = sum(TreeCount * isLive), 
              measureTrees = sum(TreeCount * isLive * (DBH > 0), na.rm = TRUE), 
              countTrees = sum(TreeCount * isLive * is.na(DBH)),
              heightTrees = sum(TreeCount * (is.na(TotalHt) == FALSE) & is.na(Ht2), na.rm = TRUE),
              brokenTrees = sum(TreeCount * (Ht2 > 0), na.rm = TRUE),
              snags = sum(TreeCount * (isLive == FALSE)), 
              measureSnags = sum(TreeCount * (isLive == FALSE) * (DBH > 0), na.rm = TRUE), 
              countSnags = sum(TreeCount * (isLive == FALSE) * is.na(DBH)))
  # measured snags
  print(trees2016 %>% filter(CompCode == "RT", isPlantation == FALSE) %>% select(StandID, Species, DBH, standAge2020), n = 35)
  # height tree counts by species group  
  trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE) %>% 
    group_by(speciesGroup) %>% 
    summarize(trees = n(), conifers = sum(isConifer), broadleaves = sum(isConifer == FALSE), plantation = sum(isPlantation),
              .groups = "drop") %>%
    bind_rows(summarize(., across(where(is.factor), ~"total"), across(where(is.numeric), sum))) %>%
    mutate(pctPlantation = 100 * plantation/trees)
  # estimated total number of trees in inventoried stands
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
  
  # ranges of predictor variables
  print(liveUnbrokenTrees2016 %>% group_by(speciesGroup) %>% 
    summarize(quantile = c(0, 0.5, 1), 
              dbh = quantile(DBH, quantile, na.rm = TRUE),
              height = quantile(TotalHt, quantile, na.rm = TRUE),
              tph = quantile(tph, quantile, na.rm = TRUE),
              ba = quantile(standBasalAreaPerHectare, quantile, na.rm = TRUE),
              bal = quantile(basalAreaLarger, quantile, na.rm = TRUE),
              aa = quantile(standBasalAreaApprox, quantile, na.rm = TRUE),
              aat = quantile(tallerApproxBasalArea, quantile, na.rm = TRUE),
              elevation = quantile(elevation, quantile, na.rm = TRUE),
              slope = quantile(slope, quantile, na.rm = TRUE),
              aspect = quantile(aspect, quantile, na.rm = TRUE),
              tsi = quantile(topographicShelterIndex, quantile, na.rm = TRUE),
              topHt = quantile(topHeight, quantile, na.rm = TRUE),
              relHt = quantile(relativeHeight, quantile, na.rm = TRUE),
              .groups = "drop"),
    n = 25)
  
  ggplot(trees2016 %>% filter(isLiveUnbroken)) + # lower violin in pair is for plantations
    geom_violin(aes(x = relativeHeight, y = speciesGroup, color = speciesGroup), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 3)) +
    ggh4x::facet_nested(rows = vars(speciesGroup, factor(isPlantation, levels = c(FALSE, TRUE), labels = c("natural regen", "plantation"))), labeller = label_wrap_gen(width = 15), scales = "free_y", switch = "y") +
    guides(color = "none") +
    labs(x = "relative height", y = NULL, color = NULL) +
    scale_color_manual(breaks = levels(liveUnbrokenTrees2016$speciesGroup), limits = levels(liveUnbrokenTrees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_y_discrete(labels = NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
  
  # distribution of height measure trees by height and diameter class
  #trees2016 %>% filter(is.na(TotalHt) == FALSE, is.na(heightClass))
  #trees2016 %>% filter(is.na(DBH) == FALSE, is.na(dbhClass)) # one red alder snag > 116 cm
  trees2016classified = trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE) %>%
    group_by(speciesGroup, heightClass, dbhClass) %>%
    summarize(n = sum(TreeCount), .groups = "drop")
  ggplot(trees2016classified) +
    geom_bar(aes(x = as.numeric(levels(heightClass))[heightClass], y = 100 * after_stat(count / tapply(count, PANEL, sum)[PANEL]), weight = n)) + # https://stackoverflow.com/questions/68227541/ggplot-geom-bar-plot-percentages-by-group-and-facet-wrap
    facet_wrap(vars(speciesGroup)) +
    labs(x = "unbroken height, m", y = "fraction of stems, %", title = "a) height classes") +
  ggplot(trees2016classified) +
    geom_bar(aes(x = as.numeric(levels(dbhClass))[dbhClass], y = 100 * after_stat(count / tapply(count, PANEL, sum)[PANEL]), weight = n)) +
    facet_wrap(vars(speciesGroup)) +
    labs(x = "DBH, cm", y = "fraction of stems, %", title = "b) DBH classes") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() &
    coord_cartesian(ylim = c(0, 10))
  
  # sizes of trees and taper limits
  treeStemLimits = get_species_limits(trees2016)
  ggplot(trees2016) +
    geom_histogram(aes(x = TotalHt, y = 100 * after_stat(count / tapply(count, PANEL, sum)[PANEL]), fill = speciesGroup), binwidth = 2, na.rm = TRUE) +
    facet_wrap(vars(speciesGroup)) +
    labs(x = "height, m", y = "percentage of trees") +
  ggplot(trees2016) +
    geom_histogram(aes(x = DBH, y = 100 * after_stat(count / tapply(count, PANEL, sum)[PANEL]), fill = speciesGroup), binwidth = 2.5, na.rm = TRUE) +
    facet_wrap(vars(speciesGroup)) +
    labs(x = "DBH, cm", y = NULL) +
  ggplot() +
    geom_line(aes(x = DBH, y = heightDiameterRatioMin), treeStemLimits, color = "grey70", linetype = "longdash", na.rm = TRUE) +
    geom_point(aes(x = DBH, y = TotalHt / (0.01 * DBH), color = speciesGroup), trees2016, alpha = 0.1, na.rm = TRUE, shape = 16) +
    geom_line(aes(x = DBH, y = heightDiameterRatioMax), treeStemLimits, color = "grey70", linetype = "longdash", na.rm = TRUE) +
    facet_wrap(vars(speciesGroup)) +
    coord_cartesian(ylim = c(1, 400)) +
    labs(x = "DBH, cm", y = "height-diameter ratio") +
  plot_layout(guides = "collect") &
    guides(color = "none", fill = "none") &
    scale_color_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) &
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), limits = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65"))

  ## Figures A1-4: species level exploratory plots
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "DF"), speciesLabel = "Douglas-fir", maxTreesMeasured = 150, omitLegends = TRUE, omitXlabels = TRUE) /
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "RA"), speciesLabel = "red alder", maxTreesMeasured = 150, distributionLegendPositionY = 0.92, plotLetters = c("d)", "e)", "f)")) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
  #ggsave("trees/height-diameter/figures/Figure A1 PSME-ALRU2.png", height = 13, width = 20, units = "cm", dpi = 250)
  
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "WH"), speciesLabel = "western hemlock", maxTreesMeasured = 150, omitLegends = TRUE) /
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "BM"), speciesLabel = "bigleaf maple", maxTreesMeasured = 150, distributionLegendPositionY = 0.92, plotLetters = c("d)", "e)", "f)"), ) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
  #ggsave("trees/height-diameter/figures/Figure A2 TSHE-ACMA3.png", height = 13, width = 20, units = "cm", dpi = 250)
    
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "OM"), speciesLabel = "Oregon myrtle", maxTreesMeasured = 150, distributionLegendPositionY = 0.92, omitXlabels = TRUE) /
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "RC"), speciesLabel = "western redcedar", maxTreesMeasured = 150, plotLetters = c("d)", "e)", "f)"), omitLegends = TRUE) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
  #ggsave("trees/height-diameter/figures/Figure A3 UMCA-THPL.png", height = 13, width = 20, units = "cm", dpi = 250)
  
  plot_exploratory(trees2016 %>% filter(isLiveUnbroken, speciesGroup == "other"), speciesLabel = "other species ", distributionLegendPositionY = 0.92) +
  plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
  #ggsave("trees/height-diameter/figures/Figure A4 other species.png", height = 1/3*(18 - 1) + 1, width = 20, units = "cm", dpi = 250)
}


## stand-level summaries and clustering
if (htDiaOptions$includeInvestigatory)
{
  treesByStand2016 = trees2016 %>% 
    group_by(StandID) %>% 
    summarize(speciesGroup = names(sort(-table(speciesGroup)))[1], # one liner for mode of character vector (https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode/8189441#8189441)
              plots = n_distinct(PlotID),
              measurePlots = n_distinct(PlotID * (is.na(DBH) == FALSE)) - any(is.na(DBH)), # plot IDs start at 1 so multiplying by is.na(DBH) == FALSE introduces zero as a plot ID
              trees = sum(TreeCount),
              live = sum(isLive * TreeCount), plantation = sum(isPlantation), retention = sum(CompCode == "RT"), snag = sum(isLive == FALSE), 
              dbh = sum(isLive * TreeCount * (DBH > 0), na.rm = TRUE), 
              height = sum(isLive * TreeCount * (is.na(TotalHt) == FALSE), na.rm = TRUE), 
              age = sum(TreeCount * (BHAge > 0), na.rm = TRUE), 
              crownRatio = sum(TreeCount * (CrownRatio > 0), na.rm = TRUE), 
              dia1 = sum(TreeCount * (Dia1 > 0), na.rm = TRUE), 
              height1 = sum(TreeCount * (Ht1 > 0), na.rm = TRUE), 
              height2 = sum(TreeCount * (Ht2 > 0), na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(speciesGroup = factor(speciesGroup, levels = levels(trees2016$speciesGroup))) # restore factor levels lost in tabling
  
  ggplot(treesByStand2016) +
    geom_histogram(aes(x = measurePlots, fill = speciesGroup), binwidth = 1) + # 492 stands with 26 plots
    coord_cartesian(xlim = c(0, 45)) +
    labs(x = "measure plots", y = "number of stands", fill = "most\ncommon\nspecies") +
  ggplot(treesByStand2016) +
    geom_histogram(aes(x = trees, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "trees counted", y = NULL, fill = "most\ncommon\nspecies") +
  ggplot(treesByStand2016) +
    geom_histogram(aes(x = dbh, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "DBH measure trees", y = "number of stands", fill = "most\ncommon\nspecies") +
  ggplot(treesByStand2016) +
    geom_histogram(aes(x = height, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "height measure trees", y = NULL, fill = "most\ncommon\nspecies") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 2, guides = "collect") &
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) &
    theme(legend.spacing.y = unit(0.2, "line"))
  
  speciesCountByStand2016 = trees2016 %>% filter(isLiveUnbroken) %>% group_by(StandID) %>%
    summarize(psmeTph = sum((speciesGroup == "DF") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              alruTph = sum((speciesGroup == "RA") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              tsheTph = sum((speciesGroup == "WH") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              acmaTph = sum((speciesGroup == "BM") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              umcaTph = sum((speciesGroup == "OM") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              thplTph = sum((speciesGroup == "RC") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              otherTph = sum((speciesGroup == "other") * meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, na.rm = TRUE),
              totalTph = psmeTph + alruTph + tsheTph, acmaTph + umcaTph + thplTph + otherTph,
              psmeBA = sum((speciesGroup == "DF") * treeBasalAreaPerHectare) / plotsInStand[1], # m²/ha
              alruBA = sum((speciesGroup == "RA") * treeBasalAreaPerHectare) / plotsInStand[1],
              tsheBA = sum((speciesGroup == "WH") * treeBasalAreaPerHectare) / plotsInStand[1],
              acmaBA = sum((speciesGroup == "BM") * treeBasalAreaPerHectare) / plotsInStand[1],
              umcaBA = sum((speciesGroup == "OM") * treeBasalAreaPerHectare) / plotsInStand[1],
              thplBA = sum((speciesGroup == "RC") * treeBasalAreaPerHectare) / plotsInStand[1],
              otherBA = sum((speciesGroup == "other") * treeBasalAreaPerHectare) / plotsInStand[1],
              totalBA = psmeBA + alruBA + tsheBA + acmaBA + umcaBA + thplBA + otherBA,
              psmePct = 100 * psmeBA / totalBA,
              alruPct = 100 * alruBA / totalBA,
              tshePct = 100 * tsheBA / totalBA,
              acmaPct = 100 * acmaBA / totalBA,
              umcaPct = 100 * umcaBA / totalBA,
              thplPct = 100 * thplBA / totalBA,
              otherPct = 100 * otherBA / totalBA,
              primaryPct = pmax(psmePct, alruPct, tshePct, acmaPct, umcaPct, thplPct, otherPct),
              secondaryPct = rev(sort(c_across(c(psmePct, alruPct, tshePct, acmaPct, umcaPct, thplPct, otherPct))))[2],
              primarySpecies = if_else(primaryPct == psmePct, "PSME",
                                       if_else(primaryPct == alruPct, "ALRU",
                                               if_else(primaryPct == tshePct, "TSHE",
                                                       if_else(primaryPct == acmaPct, "ACMA",
                                                               if_else(primaryPct == umcaPct, "UMCA",
                                                                       if_else(primaryPct == thplPct, "THPL", "other")))))),
              secondarySpecies = if_else(secondaryPct == psmePct, "PSME",
                                         if_else(secondaryPct == alruPct, "ALRU",
                                                 if_else(secondaryPct == tshePct, "TSHE",
                                                         if_else(secondaryPct == acmaPct, "ACMA",
                                                                 if_else(secondaryPct == umcaPct, "UMCA",
                                                                         if_else(secondaryPct == thplPct, "THPL", "other")))))),
              isPlantation = isPlantation[1],
              standArea = standArea[1],
              standAge2016 = standAge2016[1],
              topHeight = topHeight[1],
              qmd = sqrt(totalBA / (pi / (4 * 100^2) * totalTph)))
  #speciesCountByStand2016 %>% group_by(primarySpecies, secondarySpecies) %>% summarize(n = n())
  
  reinekeSdi = crossing(sdi = c(100 * seq(1, 9), 1000 * seq(1, 9), 10000),
                        tph = c(1, 10000)) %>%
    mutate(qmd = 25.4 * (sdi/tph)^(1/1.605)) # sdi = tph * (qmd/25.4)^1.605 => (sdi/tph)^(1/1.605) = (qmd/25.4)
  
  # Reineke SDI
  ggplot() +
    geom_path(aes(x = tph, y = qmd, group = sdi), reinekeSdi, color = "grey80", linetype = "longdash", linewidth = 0.3) +
    geom_point(aes(x = totalTph, y = qmd, color = primarySpecies, shape = secondarySpecies), speciesCountByStand2016, alpha = 0.5) +
    coord_cartesian(xlim = c(30, 2500), ylim = c(5, 107)) +
    labs(x = "trees per hectare", y = "QMD, cm", color = "primary species", shape = "secondary species") +
    scale_color_manual(breaks = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_shape_manual(breaks = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), values = c(15, 16, 17, 22, 21, 24, 25)) +
    scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000), minor_breaks = c(30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 3000, 4000)) +
    scale_y_log10(breaks = c(seq(1, 9), 10 * seq(1, 10)), minor_breaks = c(15, 25, 35, 45, 55, 65, 75, 85, 95, 150)) +
    theme(legend.spacing.y = unit(0.3, "line"))
  # TPH distribution by species and stand
  ggplot(speciesCountByStand2016) + # could also pivot to longform and wrap_facet(speciesGroup)
    geom_histogram(aes(x = psmeTph), fill = "forestgreen", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "Douglas-fir TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = alruTph), fill = "red2", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "red alder TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = tsheTph), fill = "blue2", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "western hemlock TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = acmaTph), fill = "green3", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "bigleaf maple TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = umcaTph), fill = "mediumorchid", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "Oregon myrtle TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = thplTph), fill = "firebrick", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "western redcedar TPH", y = "stands") +
  ggplot(speciesCountByStand2016) +
    geom_histogram(aes(x = otherTph), fill = "grey65", binwidth = 50) +
    coord_cartesian(xlim = c(0, 2000), ylim = c(0, 700)) +
    labs(x = "other species TPH", y = "stands") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 4, guides = "collect")
  # alternate TPH breakdown
  ggplot(speciesCountByStand2016) +
    geom_point(aes(x = totalTph - psmeTph, y = psmeTph, color = primarySpecies, shape = secondarySpecies), alpha = 0.5) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    labs(x = bquote("total of all other species, trees ha"^-1), y = bquote("Douglas-fir, trees ha"^-1), color = "primary species", shape = "secondary species") +
    scale_color_manual(breaks = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_shape_manual(breaks = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), values = c(15, 16, 17, 22, 21, 24, 25)) +
    theme(legend.spacing.y = unit(0.3, "line"))

  # basic clustering of stand types
  # k-means and mean shift perform poorly here, presumably due to being asked to partition continuous data. Data visualization
  # here could use either scaled (normalized) or unscaled distances, the former emphasizing dissimilarity in species besides
  # Douglas-fir and the latter (presumably) being more directly representative of the distribution of trees in the ground.
  basalAreaDistances = dist(speciesCountByStand2016 %>% select(isPlantation, psmeBA, alruBA, tsheBA, acmaBA, umcaBA, thplBA, otherBA))
  standHierarchyBA = hclust(basalAreaDistances, method = "ward.D") # produces the most even area distribution among hclust()'s methods
  # ggdendro::ggdendrogram(standHierarchyBA)

  speciesBasalAreaByCluster = speciesCountByStand2016 %>% mutate(clusterID = cutree(standHierarchyBA, k = 20)) %>%
    group_by(clusterID) %>%
    summarize(PSME = sum(standArea * psmeBA) / sum(standArea), # could also join clustersBA$centers
              ALRU = sum(standArea * alruBA) / sum(standArea),
              TSHE = sum(standArea * tsheBA) / sum(standArea),
              ACMA = sum(standArea * acmaBA) / sum(standArea),
              UMCA = sum(standArea * umcaBA) / sum(standArea),
              THPL = sum(standArea * thplBA) / sum(standArea),
              other = sum(standArea * otherBA) / sum(standArea),
              meanAge2016 = sum(standArea * standAge2016) / sum(standArea),
              meanTopHeight = sum(standArea * topHeight) / sum(standArea),
              meanTotalBasalArea = sum(standArea * totalBA) / sum(standArea),
              plantationArea = sum(isPlantation * standArea),
              totalArea = sum(standArea)) %>%
    arrange(totalArea) %>%
    mutate(clusterID = factor(clusterID, levels = clusterID),
           ageLabel = if_else(totalArea == max(totalArea), sprintf('bar(age) == %.0f~"years"', meanAge2016), sprintf("%.0f", meanAge2016)),
           ageLabelX = if_else(totalArea > 270, 30, totalArea + 30),
           ageLabelColor = if_else(totalArea > 270, "white", "black"),
           topHeightLabel = if_else(totalArea == max(totalArea), sprintf('bar(H)[100] == "%.1f"~"m"', meanTopHeight), sprintf('"%.1f"', meanTopHeight)),
           topHeightLabelColor = if_else(meanTotalBasalArea > 10, "white", "black"),
           topHeightLabelX = if_else(meanTotalBasalArea > 10, 1, meanTotalBasalArea + 1))
  
  # clustering
  speciesBasalAreaOneRowPerClusterSlice = speciesBasalAreaByCluster %>% group_by(clusterID) %>% slice(1)

  ggplot() +
    geom_col(aes(x = area, y = clusterID, fill = isPlantation, group = fct_rev(isPlantation)), orientation = "y", speciesBasalAreaByCluster %>% mutate(naturalRegenArea = totalArea - plantationArea) %>% 
               select(clusterID, naturalRegenArea, plantationArea) %>%
               pivot_longer(cols = c("naturalRegenArea", "plantationArea"), names_to = "isPlantation", values_to = "area") %>%
               mutate(isPlantation = factor(isPlantation, labels = c("natural regeneration", "plantation"), levels = c("naturalRegenArea", "plantationArea")))) +
    geom_text(aes(x = ageLabelX, y = clusterID, label = ageLabel), speciesBasalAreaByCluster, color = speciesBasalAreaByCluster$ageLabelColor, hjust = 0, parse = TRUE, size = 3) +
    coord_cartesian(xlim = c(0, 2550)) +
    labs(x = "area inventoried, ha", y = "stand type classification (hierarchical cluster ID)", fill = NULL, title = "a) age and area") +
    scale_fill_manual(breaks = c("natural regeneration", "plantation"), values = c("grey10", "grey35")) +
    scale_x_continuous(expand = c(0.012, 0)) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.02)) +
  ggplot() +
    geom_bar(aes(y = clusterID, fill = fct_rev(speciesGroup), weight = basalArea), speciesBasalAreaByCluster %>% select(-meanAge2016, -meanTopHeight, -starts_with("age"), -starts_with("topHeight")) %>%
               pivot_longer(cols = -c("clusterID", "meanTotalBasalArea", "plantationArea", "totalArea"), names_to = "speciesGroup", values_to = "basalArea") %>%
               mutate(speciesGroup = factor(speciesGroup, levels = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other")))) +
    geom_text(aes(x = topHeightLabelX, y = clusterID, label = topHeightLabel), speciesBasalAreaOneRowPerClusterSlice, color = speciesBasalAreaOneRowPerClusterSlice$topHeightLabelColor, hjust = 0, parse = TRUE, size = 3) +
    coord_cartesian(xlim = c(0, 83)) +
    labs(x = bquote("mean basal area, m"^2*" ha"^-1), y = NULL, fill = NULL, title = "b) top height and density") +
    scale_fill_manual(breaks = c("PSME", "ALRU", "TSHE", "ACMA", "UMCA", "THPL", "other"), labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_x_continuous(expand = c(0.008, 0)) +
    scale_y_discrete(labels = NULL) +
    theme(legend.key.height = unit(1, "line"), legend.key.width = unit(1, "line")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 1, ncol = 2, widths = c(0.4, 0.6))
  #ggsave("trees/height-diameter/figures/Figure S90 Elliott stand clusters.png", height = 10.5, width = 22, units = "cm", dpi = 250)
  
  # tree counts
  trees2016 %>% group_by(StandID) %>% summarize(standArea = standArea[1], tph = tph[1]) %>%
    summarize(areaInventoried = sum(standArea), treeCountInventoried = sum(standArea * tph)) %>%
    mutate(totalArea = 33397.1, totalTrees = treeCountInventoried * totalArea / areaInventoried)
}


## Douglas-fir site index regression: not enough data for other species
if (htDiaOptions$includeInvestigatory)
{
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
  psmeSiteIndexModelNonlinear = gsl_nls(Cruised_Si ~ b0 + b1*Elev_Mean^b2 + b3*SlopeMeanPercent^b4 + b5*planted, psmeStands2022, start = list(b0 = 120, b1 = -1E-6, b2 = 2, b3 = -2E-3, b4 = 5, b5 = 20), control = gsl_nls_control(maxiter = 250))
  c(linear = AIC(psmeSiteIndexModelLinear), nonlinear = AIC(psmeSiteIndexModelNonlinear))
  
  ggplot(psmeStands2022) + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "longdash") + 
    geom_point(aes(x = Cruised_Si, y = predict(psmeSiteIndexModelLinear), color = planted), alpha = 0.3) +
    labs(x = "measured 50-year site index, feet", y = "linear model prediction, feet", color = NULL) +
    theme(legend.position = "none") +
  ggplot(psmeStands2022) + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "longdash") + 
    geom_point(aes(x = Cruised_Si, y = predict(psmeSiteIndexModelNonlinear, psmeStands2022), color = planted), alpha = 0.3) +
    labs(x = "measured 50-year site index, feet", y = "nonlinear model prediction, feet", color = "stand age") +
    scale_color_discrete(breaks = c(FALSE, TRUE), labels = c("≥100 years", "<100 years")) +
    theme(legend.justification = c(1, 0), legend.position = c(0.98, 0.02))
  
  ggplot(psmeStands2022) +
    geom_point(aes(x = Cruised_Si, y = -residuals(psmeSiteIndexModelLinear), color = planted), alpha = 0.3, shape = 16) +
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
  liveUnbrokenTrees2016 = trees2016 %>% filter(isLiveUnbroken)
  ggplot(liveUnbrokenTrees2016) +
    geom_histogram(aes(x = DBH, y = 100 * ..count../sum(..count..), fill = speciesGroup, alpha = isPlantation), binwidth = 2.5, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, 4.4)) +
    labs(x = "DBH, cm", y = "percentage of live, unbroken stems measured", alpha = NULL, fill = NULL) +
    scale_alpha_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c(1, 0.7)) +
    scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), values = c("green3", "red2", "blue2", "cyan2", "darkorchid3", "firebrick", "grey35")) +
    theme(legend.position = "none") +
  ggplot(liveUnbrokenTrees2016) +
    geom_histogram(aes(x = 100 * ..count../sum(..count..), y = TotalHt, fill = speciesGroup, alpha = isPlantation), binwidth = 1, na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 4.4)) +
    labs(x = "percentage of live stems measured", y = "height, m", alpha = NULL, fill = NULL) +
    scale_alpha_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c(1, 0.7)) +
    scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), labels = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other"), values = c("green3", "red2", "blue2", "cyan2", "darkorchid3", "firebrick", "grey35")) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.spacing.y = unit(0.3, "line"))
}

  
## site index plots
if (htDiaOptions$includeInvestigatory)
{
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
}


## stand tree pooling
if (htDiaOptions$includeInvestigatory)
{
  # stand merge by plantation age
  # for now, proxy HX as Pacific dogwood (Cornus nuttalli)
  treesOfAge = trees2016 %>% filter(standAge2016 == 30, isLive, is.na(DBH) == FALSE) %>% 
    mutate(year = 2016L, 
           expansionFactor = meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, 
           species = case_match(Species, "DF" ~ "PSME", "RA" ~ "ALRU", "WH" ~ "TSHE", "BM" ~ "ACMA", "OM" ~ "UMCA", "RC" ~ "THPL", "CH" ~ "CHCH", "HX" ~ "CONU", "TO" ~ "LIDE"),
           condition = as.integer(isLiveUnbroken == FALSE)) %>%
    rename(stand = StandID, plot = PlotID, tree = TreeID, age = standAge2016, dbh = DBH, height = imputedHeight) %>%
    select(stand, plot, tree, species, year, age, dbh, height, expansionFactor, condition)
  #treesOfAge %>% group_by(species) %>% summarize(n = n())
  #treesOfAge %>% group_by(stand) %>% summarize(uniqueTrees = n(), tph = sum(expansionFactor))
}
