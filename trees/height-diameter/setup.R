#.libPaths(.libPaths()[2])
#install.packages(c("cmocean", "dplyr", "ggplot2", "gslnls", "magrittr", "mgcv", "nlme", "nls.multstart", "nlstools", "patchwork", "readxl", "robustbase", "sn", "stringr", "tibble", "tidyr", "writexl"))
library(cmocean)
library(dplyr)
library(forcats)
library(ggplot2)
library(gslnls)
library(magrittr)
library(mgcv)
library(nlme)
library(nls.multstart)
library(nlstools)
library(patchwork)
library(readxl)
library(robustbase)
library(sn)
library(stringr)
library(tibble)
library(tidyr)
library(writexl)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             legend.key.height = unit(0.85, "line"),
                             legend.spacing.y = unit(0, "line"),
                             legend.title = element_text(size = 10),
                             panel.border = element_blank()))
includeInvestigatory = FALSE # default to excluding plotting and other add ons in species scripts

as_row = function(model = NULL, name = NULL, fitSet = "primary", fixedWeight = NA_real_)
{
  if (is.null(model))
  {
    if (is.null(name))
    {
      stop("Name must be specified if regression is null.")
    }
    return(tibble(fitSet = fitSet, name = name, 
                  fitting = NA_character_, n = NA_real_, significant = NA_real_, fixedWeight = fixedWeight,
                  bias = NA_real_, biasNaturalRegen = NA_real_, biasPlantation = NA_real_,
                  aic = NA_real_, aicN = NA_real_, aict = NA_real_,
                  bic = NA_real_, bict = NA_real_,
                  mab = NA_real_, mapb = NA_real_,
                  mape = NA_real_, mapeNaturalRegen = NA_real_, mapePlantation = NA_real_,
                  nse = NA_real_, nseNaturalRegen = NA_real_, nsePlantation = NA_real_,
                  pae = NA_real_, paeNaturalRegen = NA_real_, paePlantation = NA_real_,
                  pearson = NA_real_, pearsonNaturalRegen = NA_real_, pearsonPlantation = NA_real_,
                  rmse = NA_real_, rmseNaturalRegen = NA_real_, rmsePlantation = NA_real_,
                  rmspe = NA_real_,
                  power = NA_real_, powerPlantation = NA_real_,
                  adaptiveWeightFraction = NA_real_, maxResidual = NA_real_))
  }
  if (("mape" %in% names(model)) == FALSE)
  {
    stop(paste("Regression for", model$name, " is missing summary statistics."))
  }
  if ((is.null(model$convInfo) == FALSE) && (model$convInfo$isConv == FALSE))
  {
    stop(paste("Regression for", model$name, " is not converged."))
  }
  
  power = NA_real_
  powerPlantation = NA_real_
  if (is.null(model$modelStruct$varStruct) == FALSE)
  {
    power = model$modelStruct$varStruct[1]
    if (length(model$modelStruct$varStruct) > 1)
    {
      powerPlantation = model$modelStruct$varStruct[2]
    }
  }
  
  residuals = residuals(model)
  return(tibble(fitSet = fitSet, name = model$name, 
                fitting = class(model)[1], n = sum(is.na(residuals) == FALSE), significant = model$significant, fixedWeight = fixedWeight,
                aic = model$aic, aicN = model$aic / nobs(model), aict = model$aict,
                bic = model$bic, bict = model$bict,
                bias = model$bias, biasNaturalRegen = model$biasNaturalRegen, biasPlantation = model$biasPlantation,
                mab = model$mab, mapb = model$mapb,
                mae = model$mae, maeNaturalRegen = model$mapeNaturalRegen, maePlantation = model$mapePlantation,
                mape = model$mape, mapeNaturalRegen = model$mapeNaturalRegen, mapePlantation = model$mapePlantation,
                nse = model$nse, nseNaturalRegen = model$nseNaturalRegen, nsePlantation = model$nsePlantation,
                pearson = model$pearson, pearsonNaturalRegen = model$pearsonNaturalRegen, pearsonPlantation = model$pearsonPlantation,
                rmse = model$rmse, rmseNaturalRegen = model$rmseNaturalRegen, rmsePlantation = model$rmsePlantation,
                rmspe = model$rmspe,
                power = power, powerPlantation = powerPlantation,
                adaptiveWeightFraction = model$adaptiveWeightFraction, maxResidual = residuals[which.max(abs(residuals))]))
}

confint_nlrob = function(regression, level = 0.99, df = df.residual(regression), weights = regression$weights)
{
  if (is.null(regression$rweights))
  {
    stop("regression$rweights is not set. Is the regression of type nlrob?")
  }
  if (is.null(weights))
  {
    stop("Either regression$weights is not set or weights were not specified.")
  }
  squaredDeviation = sum(weights * regression$rweights * residuals(regression)^2) / df.residual(regression)
  gradient = regression$m$gradient()
  levels = c((1 - level)/2, 1 - (1 - level)/2)
  parameterValues = regression$m$getPars()
  confidenceInterval = parameterValues + sqrt(diag(squaredDeviation * solve(t(gradient) %*% gradient))) %o% qt(p = levels, df = df)
  colnames(confidenceInterval) = sprintf("%g%%", 100*levels)
  rownames(confidenceInterval) = names(parameterValues)
  return(confidenceInterval)
}

# wrap calls to fitting functions for consistency of arguments and use of get_*_error()
# This is a little fragile from a code maintenance perspective as the weights need to be column in data for R to flow them
# correctly but the risk appears low and worth the simplification elsewhere.
fit_gam = function(name, formula, data, start, tDegreesOfFreedom = 8, control = gsl_nls_control(), bam = FALSE, significant = TRUE, nthreads = 1)
{
  if (formula[2] == "TotalHt()")
  {
    if (bam)
    {
      model = bam(formula = formula, data = data, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
    }
    else
    {
      model = gam(formula = formula, data = data, method = "REML", select = TRUE, weights = dbhWeight, nthreads = nthreads)
    }
    model = get_height_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  else
  {
    stopifnot("Expected response variable to be DBH." = formula[2] == "DBH()")
    if (bam)
    {
      model = bam(formula = formula, data = data, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
    }
    else
    {
      model = gam(formula = formula, data = data, method = "REML", select = TRUE, weights = heightWeight, nthreads = nthreads)
    }
    model = get_dbh_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  
  return(model)
}

fit_gsl_nls = function(name, formula, data, start, tDegreesOfFreedom = 8, control = gsl_nls_control(), significant = TRUE)
{
  if (formula[2] == "TotalHt()")
  {
    model = gsl_nls(fn = formula, data = data, start = start, weights = dbhWeight, control = control)
    model = get_height_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  else
  {
    stopifnot("Expected response variable to be DBH." = formula[2] == "DBH()")
    model = gsl_nls(fn = formula, data = data, start = start, weights = heightWeight, control = control)
    model = get_dbh_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  
  return(model)
}

fit_lm = function(name, formula, data, tDegreesOfFreedom = 8, significant = TRUE)
{
  if (formula[2] == "TotalHt()")
  {
    model = lm(formula = formula, data = data, offset = breastHeight, weights = dbhWeight)
    model = get_height_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  else
  {
    stopifnot("Expected response variable to be DBH." = formula[2] == "DBH()")
    model = lm(formula = formula, data = data, weights = heightWeight)
    model = get_dbh_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  
  return(model)
}

fit_nlrob = function(name, formula, data, start, tDegreesOfFreedom = 8, control = nls.control(), maxit = 20, significant = TRUE)
{
  if (formula[2] == "TotalHt()")
  {
    model = nlrob(formula = formula, data = data, maxit = maxit, start = start, weights = dbhWeight, control = control)
    model = get_height_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  else
  {
    stopifnot("Expected response variable to be DBH." = formula[2] == "DBH()")
    model = nlrob(formula = formula, data = data, maxit = maxit, start = start, weights = heightWeight, control = control)
    model = get_dbh_error(name = name, model = model, data = data, significant = significant, tDegreesOfFreedom = tDegreesOfFreedom)
  }
  
  return(model)
}

get_coefficients = function(regression)
{
  if (is.null(regression$coefficients))
  {
    if (is.null(regression$m))
    {
      stop(paste("Regression for", regression$name, "lacks a both a coefficients and an m property."))
    }
    coefficients = tibble(!!!set_names(regression$m$getPars(), names(regression$m$getPars())))
  }
  else
  {
    coefficients = tibble(!!!set_names(regression$coefficients, names(regression$coefficients)))
  }
  coefficients %<>% mutate(name = regression$name, fitting = class(regression)[1]) %>%
    relocate(name, fitting)
  
  if (class(regression)[1] == "lm")
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

get_dbh_error = function(name, model, data, significant = TRUE, tDegreesOfFreedom = 8, weights = data$heightWeight)
{
  model$name = name
  model$adaptiveWeightFraction = 0
  if (class(model)[1] == "nlrob")
  {
    model$adaptiveWeightFraction = sum(model$rweights != 1) / nobs(model)
  }
  if (is.null(model$weights))
  {
    # if weights aren't set, capture them (needed for nlrob(), not needed by lm(), gam(), gsl_nls()), nls())
    model$weights = weights
  }

  fittedValues = predict(model, data)
  residuals = fittedValues - data$DBH
  errorByHeightClass = tibble(heightClass = 1*floor(data$TotalHt/1) + 0.5*1, dbh = data$DBH, residual = residuals) %>%
    group_by(heightClass) %>%
    summarize(n = n(),
              meanBiasPerTree = sum(residual) / n,
              meanBiasPerTreePct = 100 * sum(residual / dbh) / n,
              .groups = "drop") %>%
    filter(n >= 10)
  
  effectiveDegreesOfFreedom = length(coef(model)) + 1 # for linear and nonlinear regressions assume one degree of freedom per model parameter
  if (is.null(model$edf) == FALSE)
  {
    effectiveDegreesOfFreedom = sum(model$edf) + 1 # for GAMs use indicated effective degrees of freedom
  }

  standardDeviation = sqrt(1/df.residual(model) * sum(weights * residuals^2)) / sqrt(weights)
  logLikelihoodGaussian = sum(dnorm(residuals, sd = standardDeviation, log = TRUE))
  logLikelihoodT = sum(dt(residuals / standardDeviation, df = tDegreesOfFreedom, log = TRUE) - log(standardDeviation))
  
  # logLik.lm() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/logLik.R)
  #  log likelihood = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w*res^2))))
  # logLik.nls() (https://github.com/wch/r-source/blob/trunk/src/library/stats/R/nls.R)
  #  log likelihood = -N/2 * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(regression$m$resid()^2)))
  #                 = 1/2 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(regression$m$resid()^2)))) if no weights are zero
  model$aic = -2*logLikelihoodGaussian + 2 * effectiveDegreesOfFreedom # calculate AIC and BIC manually because nlrob objects implement weighting differently from nls and gslnls
  model$aict = -2*logLikelihoodT + 2 * effectiveDegreesOfFreedom
  model$bias = mean(residuals)
  model$bic = -2*logLikelihoodGaussian + effectiveDegreesOfFreedom * log(nobs(model))
  model$bict = -2*logLikelihoodT + effectiveDegreesOfFreedom * log(nobs(model))
  model$mae = mean(abs(residuals))
  model$mab = sum(errorByHeightClass$n * abs(errorByHeightClass$meanBiasPerTree)) / sum(errorByHeightClass$n)
  model$mapb = sum(errorByHeightClass$n * abs(errorByHeightClass$meanBiasPerTreePct)) / sum(errorByHeightClass$n)
  model$mape = 100 * mean(abs(residuals / data$DBH))
  model$nse = 1 - sum(residuals^2) / sum((data$DBH - mean(data$DBH))^2)
  model$pearson = cor(fittedValues, data$DBH)
  model$rmse = sqrt(mean(residuals^2))
  model$rmspe = 100 * sqrt(mean((residuals / data$DBH)^2))
  model$significant = significant
  
  dataNaturalRegen = data %>% filter(isPlantation == FALSE)
  dbhNaturalRegen = predict(model, dataNaturalRegen)
  residualsNaturalRegen = dbhNaturalRegen - dataNaturalRegen$DBH
  model$biasNaturalRegen = mean(residualsNaturalRegen)
  model$maeNaturalRegen = mean(abs(residualsNaturalRegen))
  model$nseNaturalRegen = 1 - sum(residualsNaturalRegen^2) / sum((dataNaturalRegen$DBH - mean(dataNaturalRegen$DBH))^2)
  model$paeNaturalRegen = 100 * mean(abs(residualsNaturalRegen / dataNaturalRegen$DBH))
  model$pearsonNaturalRegen = cor(residualsNaturalRegen, dataNaturalRegen$DBH)
  model$rmseNaturalRegen = sqrt(mean(residualsNaturalRegen^2))
  
  dataPlantation = data %>% filter(isPlantation)
  dbhPlantation = predict(model, dataPlantation)
  residualsPlantation = dbhPlantation - dataPlantation$DBH
  model$biasPlantation = mean(residualsPlantation)
  model$maePlantation = mean(abs(residualsPlantation))
  model$nsePlantation = 1 - sum(residualsPlantation^2) / sum((dataPlantation$DBH - mean(dataPlantation$DBH))^2)
  model$paePlantation = 100 * mean(abs(residualsPlantation / dataPlantation$DBH))
  model$pearsonPlantation = cor(residualsPlantation, dataPlantation$DBH)
  model$rmsePlantation = sqrt(mean(residualsPlantation^2))
  
  return(model)
}

get_height_error = function(name, model, data, significant = TRUE, tDegreesOfFreedom = 8, weights = data$dbhWeight)
{
  model$name = name
  model$adaptiveWeightFraction = 0
  if (class(model)[1] == "nlrob")
  {
    model$adaptiveWeightFraction = sum(model$rweights != 1) / nobs(model)
  }
  if (is.null(model$weights))
  {
    model$weights = weights
  }

  dbhClassSize = 10 # cm
  fittedValues = predict(model, data)
  residuals = fittedValues - data$TotalHt
  errorByDbhClass = tibble(dbhClass = dbhClassSize*floor(data$DBH/dbhClassSize) + 0.5*dbhClassSize, height = data$TotalHt, residual = residuals) %>%
    group_by(dbhClass) %>%
    summarize(n = n(),
              meanBiasPerTree = sum(residual) / n,
              meanBiasPerTreePct = 100 * sum(residual / height) / n,
              .groups = "drop") %>%
    filter(n >= 10)

  effectiveDegreesOfFreedom = length(coef(model)) + 1
  if (is.null(model$edf) == FALSE)
  {
    effectiveDegreesOfFreedom = sum(model$edf) + 1
  }
  standardDeviation = sqrt(1/df.residual(model) * sum(weights * residuals^2)) / sqrt(weights)
  logLikelihoodGaussian = sum(dnorm(residuals, sd = standardDeviation, log = TRUE))
  logLikelihoodT = sum(dt(residuals / standardDeviation, df = tDegreesOfFreedom, log = TRUE) - log(standardDeviation))
  
  model$aic = -2*logLikelihoodGaussian + 2 * effectiveDegreesOfFreedom # see get_dbh_error(): same nlrob issue
  model$aict = -2*logLikelihoodT + 2 * effectiveDegreesOfFreedom
  model$bias = mean(residuals)
  model$bic = -2*logLikelihoodGaussian + effectiveDegreesOfFreedom * log(nobs(model))
  model$bict = -2*logLikelihoodT + effectiveDegreesOfFreedom * log(nobs(model))
  model$mab = sum(errorByDbhClass$n * abs(errorByDbhClass$meanBiasPerTree)) / sum(errorByDbhClass$n)
  model$mapb = sum(errorByDbhClass$n * abs(errorByDbhClass$meanBiasPerTreePct)) / sum(errorByDbhClass$n)
  model$mae = mean(abs(residuals))
  model$mape = 100 * mean(abs(residuals / data$TotalHt))
  model$nse = 1 - sum(residuals^2) / sum((data$TotalHt - mean(data$TotalHt))^2)
  model$pearson = cor(fittedValues, data$TotalHt)
  model$rmse = sqrt(mean(residuals^2))
  model$rmspe = 100 * sqrt(mean((residuals / data$TotalHt)^2))
  model$significant = significant
                          
  dataNaturalRegen = data %>% filter(isPlantation == FALSE)
  heightNaturalRegen = predict(model, dataNaturalRegen)
  residualsNaturalRegen = heightNaturalRegen - dataNaturalRegen$TotalHt
  model$biasNaturalRegen = mean(residualsNaturalRegen)
  model$maeNaturalRegen = mean(abs(residualsNaturalRegen))
  model$nseNaturalRegen = 1 - sum(residualsNaturalRegen^2) / sum((dataNaturalRegen$TotalHt - mean(dataNaturalRegen$TotalHt))^2)
  model$paeNaturalRegen = 100 * mean(abs(residualsNaturalRegen / dataNaturalRegen$TotalHt))
  model$pearsonNaturalRegen = cor(residualsNaturalRegen, dataNaturalRegen$TotalHt)
  model$rmseNaturalRegen = sqrt(mean(residualsNaturalRegen^2))
  
  dataPlantation = data %>% filter(isPlantation)
  heightPlantation = predict(model, dataPlantation)
  residualsPlantation = heightPlantation - dataPlantation$TotalHt
  model$biasPlantation = mean(residualsPlantation)
  model$maePlantation = mean(abs(residualsPlantation))
  model$nsePlantation = 1 - sum(residualsPlantation^2) / sum((dataPlantation$TotalHt - mean(dataPlantation$TotalHt))^2)
  model$paePlantation = 100 * mean(abs(residualsPlantation / dataPlantation$TotalHt))
  model$pearsonPlantation = cor(residualsPlantation, dataPlantation$TotalHt)
  model$rmsePlantation = sqrt(mean(residualsPlantation^2))
  
  return(model)
}

impute_basal_area = function(Species, heightInM, isPlantation)
{
  # preferred fits from ends of HtDia PSME.R, ALRU2.R, THSE.R, ...
  basalAreaInM2 = recode(Species,
                         DF = (0.6895239663 - 0.4132029239 * isPlantation) * (exp((0.0003122770 - 0.0005283067 * isPlantation) * (heightInM - 1.37)^(1.9075142639 - 0.1005401747 * isPlantation)) - 1),
                         RA = 7.159320e+02 * (exp(4.270042e-07 * (heightInM - 1.37)^(1.904561e+00 - 1.121828e-01 * isPlantation)) - 1),
                         WH = (1.068577e-04 + 5.234021e-05 * isPlantation) * (heightInM - 1.37)^(2.178129e+00 - 1.778726e-01 * isPlantation),
                         BM = 5.348649e+02 * (exp(3.583647e-07 * (heightInM - 1.37)^(2.216220e+00 - 1.465682e-01 * isPlantation)) - 1),
                         OM = 1.3558236 * (exp(0.0001969 * (heightInM - 1.37)^(2.0561921 - 0.2726091 * isPlantation)) - 1),
                         RC = 3.445860e+02 * (exp(0.0001969 * (heightInM - 1.37)^(6.923252e-07 - 2.181920e+00 * isPlantation)) - 1),
                         .default = 1.828744e+02 * (exp(5.045842e-07 * (heightInM - 1.37)^(2.306458e+00 - 1.884820e-01 * isPlantation)) - 1))
  basalAreaInM2 = if_else(basalAreaInM2 < 0.25 * pi * (0.01 * 2.54)^2, 0.25 * pi * (0.01 * 2.54)^2, basalAreaInM2) # clamp regressions to minimum cruised basal area: blocks implied negative DBHes
  return(replace_na(basalAreaInM2, 0.25 * pi * (0.01 * 10)^2)) # assume any tree without a height is 10 cm DBH
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

plot_exploratory = function(liveUnbrokenTrees, plotLetters = c("a)", "b)", "c)"), speciesLabel = NULL, distributionLegendPositionY = 1, maxTreesMeasured = 400, omitLegends = FALSE, omitXlabels = FALSE)
{
  dbhQuantiles = liveUnbrokenTrees %>% mutate(diameterClass = 2.5 * (ceiling(DBH / 2.5) - 0.5)) %>% group_by(diameterClass) %>%
    summarize(count = n(), quantiles = c("min", "q025", "q10", "q20", "q25", "q30", "q40", "median", "q60", "q70", "q75", "q80", "q90", "q975", "max"), height = quantile(TotalHt, probs = c(0, 0.025, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.975, 1), na.rm = TRUE), mean = mean(TotalHt, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = quantiles, values_from = height)
  heightQuantiles = liveUnbrokenTrees %>% mutate(heightClass = 1 * (ceiling(TotalHt / 1) - 0.5)) %>% group_by(heightClass) %>%
    summarize(count = n(), quantiles = c("min", "q025", "q10", "q20", "q25", "q30", "q40", "median", "q60", "q70", "q75", "q80", "q90", "q975", "max"), dbh = quantile(DBH, probs = c(0, 0.025, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.975, 1), na.rm = TRUE), mean = mean(DBH, na.rm = TRUE), .groups = "drop") %>%
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
    geom_bin_2d(aes(x = DBH, y = TotalHt, fill = after_stat(count)), liveUnbrokenTrees %>% filter(is.na(TotalHt) == FALSE), binwidth = c(2.5, 1)) +
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
    theme(legend.key.height = unit(0.95, "line"), legend.justification = c(1, 0), legend.position = treeLegendPosition, legend.title = element_text(size = 9.5), legend.spacing.y = unit(0.25, "line")) +
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
                                plots2016 %>% select(STAND, PltInteger, elevation, slope, aspect, topographicShelterIndex, x, y) %>% rename(PlotID = PltInteger),
                                by = c("PlotID")) %>%
  rename(standAge2020 = Age_2020) %>%
  mutate(speciesGroup = factor(if_else(Species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), Species, "other"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other")),
         isConifer = Species %in% c("DF", "WH", "RC", "SS", "CX", "PC", "PY", "GF", "LP"),
         isLiveUnbroken = (CompCode %in% c("BT", "D.", "SN")) == FALSE,
         BHAge = na_if(BHAge, 0), # years
         DBH = na_if(2.54 * DBH, 0), # inches to cm
         Dia1 = na_if(2.54 * Dia1, 0),
         Ht1 = na_if(0.3048 * Ht1, 0), # feet to m
         Ht2 = na_if(0.3048 * Ht2, 0),
         TotalHt = na_if(0.3048 * TotalHt, 0),
         basalArea = 0.25 * pi * (0.01*DBH)^2, # m² 
         breastHeight = 1.37, # m, used for offset in lm() height regressions
         dbhWeightDefault = DBH^-1,
         heightWeightDefault = TotalHt^-2,
         isPlantation = standAge2020 < 75,
         SampleFactor = 2.47105 * SampleFactor, # trees per acre to trees per hectare
         standArea = 0.404686 * GrossAc,  # ac to ha
         standSampleFactor = mean(SampleFactor),
         treeBasalAreaPerHectare = SampleFactor * if_else(SamplingMethod == "BAF", 0.092903, basalArea), # m²/ha, conversion factor is either 2.47105 * 0.092903 = 0.229568 m²/ha / ft²/ac or 2.47105 ac/ha
         heightDiameterRatio = TotalHt / (0.01 * DBH), # (DBH conversion from cm to m)
         imputedHeight = if_else(is.na(TotalHt) == FALSE, TotalHt, if_else(is.na(DBH) == FALSE, impute_height(Species, DBH, isPlantation), NA_real_)), # where possible, perform basic height imputation
         treeBasalAreaPerHectareApprox = SampleFactor * if_else(SamplingMethod == "BAF", 0.092903, impute_basal_area(Species, TotalHt, isPlantation))) %>% # m²/ha, stack basal area regression on height regression when possible; BAF has to be used with prism trees
  select(-GrossAc) %>%
  group_by(PlotID) %>%
  mutate(plotTrees = sum((CompCode %in% c("D.", "SN") == FALSE) * SampleFactor),
         plotTreesWithDbh = sum(if_else(is.na(DBH), 0, SampleFactor)),
         plotContributionToStandBasalArea = 0,
         plotContributionToStandBasalArea = replace(plotContributionToStandBasalArea, 1, sum(treeBasalAreaPerHectare)),
         plotContributionToStandApproxBasalArea = 0,
         plotContributionToStandApproxBasalArea = replace(plotContributionToStandApproxBasalArea, 1, sum(treeBasalAreaPerHectareApprox))) %>%
  group_by(StandID) %>%
  arrange(desc(isLiveUnbroken), desc(DBH), .by_group = TRUE) %>% # put largest diameter live trees first in each stand for calculating BAL (numbers sort before NA)
  mutate(plotsInStand = length(unique(PlotID)),
         standBasalAreaPerHectare = sum(plotContributionToStandBasalArea) / plotsInStand, # m²/ha
         standBasalAreaApprox = sum(plotContributionToStandApproxBasalArea) / plotsInStand, # m²/ha
         basalAreaLarger = (cumsum(isLiveUnbroken * treeBasalAreaPerHectare) - treeBasalAreaPerHectare[1]) / plotsInStand, # m²/ha
         treeTphContribution = SampleFactor / plotsInStand, # trees per hectare
         tph = sum(treeTphContribution)) %>% # stand trees per hectare
  arrange(desc(isLiveUnbroken), desc(imputedHeight), .by_group = TRUE) %>% # put tallest live trees without broken tops first in each stand
  mutate(remainingTopHeightTph = pmax(100 - cumsum(if_else(is.na(TotalHt), 0, treeTphContribution)), 0), # remaining TPH contribution to H100 definition of top height, trees not measured for TotalHt are skipped
         remainingTopHeightFraction = remainingTopHeightTph / treeTphContribution,
         topHeightWeight = if_else(remainingTopHeightFraction >= 1, 1, if_else(remainingTopHeightFraction > 0, remainingTopHeightFraction, 0)), # clamp remaining fraction to [0, 1] to get individual trees' contributions to the top height average
         topHeight = sum(topHeightWeight * TotalHt, na.rm = TRUE) / sum((is.na(TotalHt) == FALSE) * topHeightWeight), # m, tallest 100 trees per hectare
         relativeHeight = TotalHt / topHeight, # individual trees' heights as a fraction of top height, may be greater than 1, especially for retention trees (debatable if imputed heights should be included but, for now, trees not measured for height are left with NA relative height)
         tallerApproxBasalArea = (cumsum(isLiveUnbroken * treeBasalAreaPerHectareApprox) - treeBasalAreaPerHectareApprox[1]) / plotsInStand,
         tallerTph = cumsum(isLiveUnbroken * SampleFactor) / plotsInStand) %>% 
  ungroup()

if (includeInvestigatory)
{
  print(trees2016 %>% filter(is.na(elevation)) %>% group_by(PlotID) %>% summarize(trees = n(), .groups = "drop"), n = 51)
  ggplot(trees2016) + geom_histogram(aes(x = standBasalAreaApprox))
  
  print(tibble(treeTphContribution = rep(25.5, 30) / 5) %>% 
          mutate(remainingTph = 100 - cumsum(treeTphContribution),
                 remainingFraction = (treeTphContribution + remainingTph) / treeTphContribution,
                 weight = if_else(remainingFraction >= 1, 1, if_else(remainingFraction > 0, remainingFraction, 0))),
        n = 40)
  
  plotTreeProperties = trees2016 %>% group_by(PlotID) %>%
    summarize(liveTrees = sum(CompCode %in% c("D.", "SN") == FALSE), snags = sum(CompCode %in% c("D.", "SN")), tph = plotTrees[1], primarySpecies = unique(Species)[which.max(tabulate(match(Species, unique(Species))))], stemsWithDbh = sum(is.na(DBH) == FALSE), stemsWithHeight = sum((is.na(TotalHt) == FALSE) | (is.na(Ht2) == FALSE))) %>% # mode of species
    mutate(stems = liveTrees + snags) %>%
    rename(PltInteger = PlotID) # for GIS joins
  plotTreeProperties %>% summarize(plots = n(), measure = sum(stemsWithDbh > 0), count = n() - measure)
  #write_xlsx(plotTreeProperties, "GIS/Trees/2015-16 cruise/CruisePlots_All_treeProperties.xlsx")
}


## data tabulation and basic plotting
if (includeInvestigatory)
{
  trees2016summary = trees2016 %>% 
    mutate(isLive = CompCode %in% c("D.", "SN") == FALSE) %>%
    #mutate(speciesClassification = if_else(Species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), Species, "other")) %>%
    #group_by(speciesClassification) %>%
    group_by(Species) %>% 
    summarize(stands = length(unique(StandID)),
              pctStems = 100 * n() / nrow(trees2016), 
              trees = n(),
              live = sum(isLive), plantation = sum(isPlantation), retention = sum(CompCode == "RT"), snag = sum(isLive == FALSE), 
              dbh = sum(isLive & (DBH > 0), na.rm = TRUE), 
              height = sum(isLive & (TotalHt > 0), na.rm = TRUE), 
              age = sum(BHAge > 0, na.rm = TRUE), 
              crownRatio = sum(CrownRatio > 0, na.rm = TRUE), 
              dia1 = sum(Dia1 > 0, na.rm = TRUE), height1 = sum(Ht1 > 0, na.rm = TRUE), 
              height2 = sum(Ht2 > 0, na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(pctStands = 100 * stands / length(unique(trees2016$StandID))) %>%
    arrange(desc(trees))
  print(trees2016summary, n = 25)
  
  trees2016 %>% summarize(live = sum(CompCode %in% c("D.", "SN") == FALSE), measureTrees = sum((CompCode %in% c("D.", "SN") == FALSE) * (DBH > 0), na.rm = TRUE), countTrees = sum((CompCode %in% c("D.", "SN") == FALSE) * is.na(DBH)),
                          snag = sum(CompCode %in% c("D.", "SN")), measureSnag = sum(CompCode %in% c("D.", "SN") * (DBH > 0), na.rm = TRUE), countSnag = sum(CompCode %in% c("D.", "SN") * is.na(DBH)))
  
  trees2016 %>% filter((CompCode %in% c("D.", "SN")) == FALSE, TotalHt > 0 | Ht2 > 0) %>% summarize(n = n()) # measured snags
  print(trees2016 %>% filter(CompCode == "RT", isPlantation == FALSE) %>% select(StandID, Species, DBH, standAge2020), n = 35)
  
  trees2016 %>% group_by(StandID) %>% summarize(standArea = standArea[1], tph = tph[1]) %>%
    summarize(trees = sum(standArea * tph))
  
  liveUnbrokenTrees2016 %>% filter(is.na(TotalHt) == FALSE) %>% 
    summarize(conifers = sum(isConifer), broadleaves = sum(isConifer == FALSE))
  
  liveUnbrokenTrees2016 %>% filter(is.na(TotalHt) == FALSE) %>% 
    group_by(speciesGroup) %>% 
    summarize(trees = n(), plantation = sum(isPlantation), pctPlantation = 100 * plantation/trees)
  
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
  
  ggplot(liveUnbrokenTrees2016) + # lower violin in pair is for plantations
    geom_violin(aes(x = relativeHeight, y = speciesGroup, color = speciesGroup), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 3)) +
    ggh4x::facet_nested(rows = vars(speciesGroup, factor(isPlantation, levels = c(FALSE, TRUE), labels = c("natural regen", "plantation"))), labeller = label_wrap_gen(width = 15), scales = "free_y", switch = "y") +
    guides(color = "none") +
    labs(x = "relative height", y = NULL, color = NULL) +
    scale_color_manual(breaks = levels(liveUnbrokenTrees2016$speciesGroup), limits = levels(liveUnbrokenTrees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
    scale_y_discrete(labels = NULL) +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
}


## stand-level summaries
if (includeInvestigatory)
{
  standSummary2016 = trees2016 %>% 
    mutate(isLive = CompCode %in% c("D.", "SN") == FALSE) %>%
    group_by(StandID) %>% 
    summarize(speciesGroup = names(sort(-table(speciesGroup)))[1], # one liner for mode of character vector (https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode/8189441#8189441)
              plots = length(unique(PlotID)),
              measurePlots = length(unique(PlotID * (is.na(DBH) == FALSE))) - any(is.na(DBH)), # plot IDs start at 1 so multiplying by is.na(DBH) == FALSE introduces zero as a plot ID
              trees = n(),
              live = sum(isLive), plantation = sum(isPlantation), retention = sum(CompCode == "RT"), snag = sum(isLive == FALSE), 
              dbh = sum(isLive & (DBH > 0), na.rm = TRUE), 
              height = sum(isLive & (TotalHt > 0), na.rm = TRUE), 
              age = sum(BHAge > 0, na.rm = TRUE), 
              crownRatio = sum(CrownRatio > 0, na.rm = TRUE), 
              dia1 = sum(Dia1 > 0, na.rm = TRUE), height1 = sum(Ht1 > 0, na.rm = TRUE), 
              height2 = sum(Ht2 > 0, na.rm = TRUE), 
              .groups = "drop") %>%
    mutate(speciesGroup = factor(speciesGroup, levels = levels(trees2016$speciesGroup))) # restore factor levels lost in tabling
  
  ggplot(standSummary2016) +
    geom_histogram(aes(x = measurePlots, fill = speciesGroup), binwidth = 1) + # 492 stands with 26 plots
    coord_cartesian(xlim = c(0, 45)) +
    labs(x = "measure plots", y = "number of stands", fill = "most\ncommon\nspecies") +
  ggplot(standSummary2016) +
    geom_histogram(aes(x = trees, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "trees counted", y = NULL, fill = "most\ncommon\nspecies") +
  ggplot(standSummary2016) +
    geom_histogram(aes(x = dbh, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "DBH measure trees", y = "number of stands", fill = "most\ncommon\nspecies") +
  ggplot(standSummary2016) +
    geom_histogram(aes(x = height, fill = speciesGroup), binwidth = 5) +
    coord_cartesian(xlim = c(0, 260), ylim = c(0, 170)) +
    labs(x = "height measure trees", y = NULL, fill = "most\ncommon\nspecies") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(nrow = 2, ncol = 2, guides = "collect") &
    scale_fill_manual(breaks = levels(trees2016$speciesGroup), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) &
    theme(legend.spacing.y = unit(0.2, "line"))
}


## Douglas-fir site index regression: not enough data for other species
if (includeInvestigatory)
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
  psmeSiteIndexModelNonlinear = nls(Cruised_Si ~ b0 + b1*Elev_Mean^b2 + b3*SlopeMeanPercent^b4 + b5*planted, psmeStands2022, start = list(b0 = 120, b1 = -1E-6, b2 = 2, b3 = -2E-3, b4 = 5, b5 = 20), control = list(maxiter = 250))
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
if (includeInvestigatory)
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


## persist workspace to disk
#save.image("trees/height-diameter/height-diameter.Rdata")
#load("trees/height-diameter/height-diameter.Rdata")