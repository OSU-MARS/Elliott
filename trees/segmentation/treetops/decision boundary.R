# assumes library()s and treetopDataDsm, treetopDataChm, and treetopDataCmm from treetops.R setup
library(rsample) # blocking by merge cluster not needed for radius accuracy assessment so vfold_cv() is ok

handlers(global = TRUE)
handlers("cli")
plan(multisession, workers = 0.5 * future::availableCores()) # no gain for form selection with vfold_cv() but effective for best fit searches and classifying treetops in tiles

fit_radius_power = function(maximaData, startingParameters, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
{
  progressBar = progressor(steps = folds * repetitions)
  
  fitFunction = function(dataFold)
  {
    #dataFold = splitsAndFits$splits[[1]]
    foldTrainingData = analysis(dataFold)
    get_radius_accuracy_power_cv = function(coefficients)
    {
      return(get_radius_accuracy_power(foldTrainingData, coefficients))
    }
    
    fit = optim(par = startingParameters, get_radius_accuracy_power_cv, method = "BFGS", control = list(fnscale = -1))
    names(fit$par) = c("a0", "a1", "b1")
    a0 = fit$par[1]
    a1 = fit$par[2]
    b1 = fit$par[3]

    validationData = assessment(dataFold)
    treetopPrediction = factor(validationData$radius > (a0 + a1 * validationData$height^b1), levels = c(FALSE, TRUE))
    
    confusion = caret::confusionMatrix(validationData$isTreetopRadius, treetopPrediction)
    overallAccuracyByHeightClass = validationData %>% select(isTreetopRadius, height) %>% 
      mutate(heightClass = round(height),
             treetopPrediction = treetopPrediction) %>%
      group_by(heightClass) %>%
      summarize(n = n(), overallAccuracy = sum(isTreetopRadius == treetopPrediction) / n(), .groups = "drop")
    
    progressBar()
    return(tibble(fit = list(fit),
                  overallAccuracy = confusion$overall["Accuracy"],
                  confusionMatrix = list(confusion),
                  overallAccuracyByHeight = list(overallAccuracyByHeightClass)))
  }
  
  splitsAndFits = vfold_cv(maximaData, v = folds, repeats = repetitions) %>%
    mutate(fit = purrr::map(splits, fitFunction)) %>% 
    select(-splits) %>% 
    unnest(fit)
  return(rename_vfold_cv_ids(splitsAndFits))
}

fit_radius_quadratic = function(maximaData, startingParameters, folds = treetopOptions$folds, repetitions = treetopOptions$repetitions)
{
  progressBar = progressor(steps = folds * repetitions)
  
  fitFunction = function(dataFold)
  {
    #dataFold = splitsAndFits$splits[[1]]
    foldTrainingData = analysis(dataFold)
    get_radius_accuracy_dsm_quadratic_cv = function(coefficients)
    {
      return(get_radius_accuracy_quadratic(foldTrainingData, coefficients))
    }
    
    fit = optim(par = startingParameters, get_radius_accuracy_dsm_quadratic_cv, method = "BFGS", control = list(fnscale = -1))
    names(fit$par) = c("a0", "a1", "a2")
    a0 = fit$par[1]
    a1 = fit$par[2]
    a2 = fit$par[3]
    
    validationData = assessment(dataFold)
    treetopPrediction = factor(validationData$radius > (a0 + a1 * validationData$height + a2 * validationData$height^2), levels = c(FALSE, TRUE))
    
    confusion = caret::confusionMatrix(validationData$isTreetopRadius, treetopPrediction)
    overallAccuracyByHeightClass = validationData %>% select(isTreetopRadius, height) %>% 
      mutate(heightClass = round(height),
             treetopPrediction = treetopPrediction) %>%
      group_by(heightClass) %>%
      summarize(n = n(), overallAccuracy = sum(isTreetopRadius == treetopPrediction) / n(), .groups = "drop")
    
    progressBar()
    return(tibble(fit = list(fit),
                  overallAccuracy = confusion$overall["Accuracy"],
                  confusionMatrix = list(confusion),
                  overallAccuracyByHeight = list(overallAccuracyByHeightClass)))
  }
  
  splitsAndFits = vfold_cv(maximaData, v = folds, repeats = repetitions) %>%
    mutate(fit = purrr::map(splits, fitFunction)) %>% 
    select(-splits) %>% 
    unnest(fit)
  return(rename_vfold_cv_ids(splitsAndFits))
}

get_radius_accuracy_quadratic = function(data, coefficients)
{
  a0 = coefficients[1]
  a1 = coefficients[2]
  a2 = coefficients[3]
  isTreetopRadiusParameterSearch = factor(data$radius > (a0 + a1 * data$height + a2 * data$height^2), levels = c(FALSE, TRUE)) 
  return(caret::confusionMatrix(isTreetopRadiusParameterSearch, data$isTreetopRadius)$overall["Accuracy"])
}

get_radius_accuracy_power = function(data, coefficients)
{
  a0 = coefficients[1]
  a1 = coefficients[2]
  b1 = coefficients[3]
  isTreetopRadiusParameterSearch = factor(data$radius > (a0 + a1 * data$height^b1), levels = c(FALSE, TRUE))
  return(caret::confusionMatrix(isTreetopRadiusParameterSearch, data$isTreetopRadius)$overall["Accuracy"])
}

print_fits_power = function(fits)
{
  tibble(fits = fits) %>% unnest_wider(col = "fits") %>% unnest_wider(col = "par", names_sep = "") %>%
    rename(a0 = par1, a1 = par2, b1 = par3, overallAccuracy = value) %>%
    mutate(a0 = tibble::num(a0, digits = 8), a1 = tibble::num(a1, digits = 8), b1 = tibble::num(b1, digits = 8), overallAccuracy = tibble::num(overallAccuracy, digits = 6)) %>%
    arrange(desc(overallAccuracy))
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
# power       0.42644864 + 0.03226514 h^1.02543385                        0.929408 BFGS -> 2x25 cross validation + Kolmogorov-Smirnov selected
#             0.42497193 + 0.03243895 h^1.02453600                        0.929406 
#             0.42495405 + 0.03256494 h^1.02310747                        0.929402
#             0.42666305 + 0.03225818 h^1.02559025                        0.929400 
#             0.42523727 + 0.03232708 h^1.02523438                        0.929397
# cubic       0.437904 + 0.0468697 h - 5.203931e-04 h² + 4.891224e-06 h³  0.928593 Nelder-Mead but of questionable utility
#             CHM
# power       0.51709010 + 0.11263928 h^0.67990521                        0.904037 BFGS -> 2x25 cross validation + Kolmogorov-Smirnov selected
#             0.49809213 + 0.12173326 h^0.66239850                        0.904026
#             0.49794017 + 0.12201314 h^0.66189253                        0.904017
#             0.48241636 + 0.11614377 h^0.68408158                        0.904012
#             0.53533571 + 0.11274569 h^0.67705943                        0.904010
#             CMM
# power      -0.22250373 + 0.05153255 h^0.93888880                        0.866874 BFGS -> 2x25 cross validation + Kolmogorov-Smirnov selected
#            -0.17995372 + 0.04011842 h^0.99411243                        0.866528
#            -0.15315471 + 0.04290902 h^0.97412857                        0.866461
#            -0.14833223 + 0.04001325 h^0.98906210                        0.866461
#            -0.14567118 + 0.03973732 h^0.99257187                        0.866372
#
#             0.11816692 + 0.03506441 h^0.97880166                        0.861297
#             0.11812274 + 0.04008985 h^0.93729610                        0.861040
#             0.15394049 + 0.03312562 h^0.98698206                        0.860773
#             0.15713721 + 0.03580347 h^0.95656241                        0.860717
#             0.17222569 + 0.03045372 h^1.00218924                        0.860650
# PRELIMINARY: change from treetopDataChm and treetopDataCmm to acceptedTreetops46chm and acceptedTreetops46cmm if/when those layers are manually reviewed and editied.
with_progress({
  progressBar = progressor(steps = 100) # 6m for 100 iterations @ 8 workers
  dsmFits = future_map(1:100, function(iteration) # future_map() requires an argument, even if unused
    {
      dsmFitPower = optim(par = c(0.420, 0.033, 1.02) + c(0.01, 0.001, 0.01) * runif(3), method = "BFGS", control = list(fnscale = -1, trace = 0), fn = function(coefficients)
      {
        return(get_radius_accuracy_power(treetopDataDsm, coefficients))
      })
      progressBar()
      return(dsmFitPower)
    }, .options = furrr_options(seed = TRUE))
})
  
print_fits_power(dsmFits)

with_progress({ # ~1 m
  progressBar = progressor(steps = 100)
  chmFits = future_map(1:100, function(iteration)
  {
    chmFitPower = optim(par = c(0.48, 0.11, 0.65) + c(0.01, 0.001, 0.01) * runif(3), method = "BFGS", control = list(fnscale = -1, trace = 0), fn = function(coefficients)
    {
      return(get_radius_accuracy_power(treetopDataChm, coefficients))
    })
    progressBar()
    return(chmFitPower)
  }, .options = furrr_options(seed = TRUE))
})

print_fits_power(chmFits)

with_progress({
  progressBar = progressor(steps = 100)
  cmmFits = future_map(1:100, function(iteration)
  {
    #cmmFitPower = optim(par = c(-0.15, 0.04, 0.98) + c(0.01, 0.001, 0.01) * runif(3), method = "BFGS", control = list(fnscale = -1, trace = 0), fn = function(coefficients)
    cmmFitPower = optim(par = c(0.2, 0.03, 0.98) + c(0.01, 0.001, 0.01) * runif(3), method = "BFGS", control = list(fnscale = -1, trace = 0), fn = function(coefficients)
    {
      return(get_radius_accuracy_power(treetopDataCmm, coefficients))
    })
    progressBar()
    return(cmmFitPower)
  }, .options = furrr_options(seed = TRUE))
})

print_fits_power(cmmFits)

## DSM
if (treetopOptions$includeSetup)
{
  powerStart = Sys.time() # ~63s, 9900X
  radiusDsmAccuracyPower = fit_radius_power(treetopDataDsm, c(0.420, 0.033, 1.02))
  Sys.time() - powerStart
  saveRDS(radiusDsmAccuracyPower, "trees/segmentation/treetops/radius DSM power s4268 458k 2x25.Rds")
  
  quadraticStart = Sys.time()
  radiusDsmAccuracyQuadratic = fit_radius_quadratic(treetopDataDsm, c(0.420, 0.0374, -0.00005))
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
    treetopDataDsm %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height), levels = c(FALSE, TRUE))) 
    confusion = caret::confusionMatrix(treetopDataDsm$isTreetopRadiusParameterSearch, treetopDataDsm$isTreetopRadius)
    return(confusion$overall["Accuracy"])
  }
  
  get_radius_accuracy_dsm_quadratic = function(coefficients)
  {
    return(get_radius_accuracy_quadratic(treetopDataDsm, coefficients))
  }
  
  get_radius_accuracy_dsm_cubic = function(coefficients)
  {
    a0 = coefficients[1]
    a1 = coefficients[2]
    a2 = coefficients[3]
    a3 = coefficients[4]
    treetopDataDsm %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height + a2 * height^2 + a3 * height^3), levels = c(FALSE, TRUE))) 
    confusion = caret::confusionMatrix(treetopDataDsm$isTreetopRadiusParameterSearch, treetopDataDsm$isTreetopRadius)
    return(confusion$overall["Accuracy"])
  }
  
  dsmFitQuadratic = optim(par = c(0.420, 0.0374, -0.00005), get_radius_accuracy_dsm_quadratic, method = "CG", control = list(fnscale = -1, trace = 1)) # negative fnscale maximizes
  dsmFitQuadratic
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
  #   treetopDataDsm %<>% mutate(isTreetopRadiusParameterSearch = factor(radius > (a0 + a1 * height + a2 * height^2), levels = c(FALSE, TRUE))) 
  #   #confusion = caret::confusionMatrix(treetopDataDsm$isTreetopRadiusParameterSearch, treetopDataDsm$isTreetopRadius)
  #   confusion = caret::confusionMatrix(treetopDataDsm$isTreetopRadiusParameterSearch, treetopDataDsm$isTreetopRadius)
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


## CHM
if (treetopOptions$includeSetup)
{
  powerStart = Sys.time() # ~6.5m, 9900X
  radiusChmAccuracyPower = fit_radius_power(treetopDataChm, c(0.48, 0.11, 0.65))
  Sys.time() - powerStart
  #saveRDS(radiusChmAccuracyPower, "trees/segmentation/treetops/radius CHM power s4268 458k 2x25.Rds")
  
  quadraticStart = Sys.time()
  radiusChmAccuracyQuadratic = fit_radius_quadratic(treetopDataChm, c(0.420, 0.039, -0.0001))
  Sys.time() - quadraticStart
  
  # power preferred @ p = 0.0008
  tibble(power = tibble::num(mean(radiusChmAccuracyPower$overallAccuracy), digits = 5), 
         quad = tibble::num(mean(radiusChmAccuracyQuadratic$overallAccuracy), digits = 5), 
         kolmoP = ks.test(radiusChmAccuracyPower$overallAccuracy, radiusChmAccuracyQuadratic$overallAccuracy)$p.value)

  radiusChmAccuracyPower %>% unnest_wider(col = "fit") %>% unnest_wider(col = "par", names_sep = "") %>%
    reframe(a0 = range(para0), a1 = range(para1), b1 = range(parb1))
  radiusChmAccuracyQuadratic %>% unnest_wider(col = "fit") %>% unnest_wider(col = "par", names_sep = "") %>%
    reframe(a0 = range(para0), a1 = range(para1), a2 = range(para2))
}


## CMM
if (treetopOptions$includeSetup)
{
  powerStart = Sys.time() # ~58s, 9900X
  radiusCmmAccuracyPower = fit_radius_power(treetopDataCmm, c(-0.15, 0.04, 0.98))
  Sys.time() - powerStart
  saveRDS(radiusCmmAccuracyPower, "trees/segmentation/treetops/radius CMM power s4268 458k 2x25.Rds")
  
  quadraticStart = Sys.time() # ~1.9m, 9900X
  radiusCmmAccuracyQuadratic = fit_radius_quadratic(treetopDataCmm, c(0.420, 0.0374, -0.00005))
  Sys.time() - quadraticStart
  
  tibble(power = tibble::num(mean(radiusCmmAccuracyPower$overallAccuracy), digits = 5), 
         quad = tibble::num(mean(radiusCmmAccuracyQuadratic$overallAccuracy), digits = 5), 
         kolmoP = ks.test(radiusCmmAccuracyPower$overallAccuracy, radiusCmmAccuracyQuadratic$overallAccuracy)$p.value)
  
  radiusCmmAccuracyPower %>% unnest_wider(col = "fit") %>% unnest_wider(col = "par", names_sep = "") %>%
    reframe(a0 = range(para0), a1 = range(para1), b1 = range(parb1))
  radiusCmmAccuracyQuadratic %>% unnest_wider(col = "fit") %>% unnest_wider(col = "par", names_sep = "") %>%
    reframe(a0 = range(para0), a1 = range(para1), a2 = range(para2))
}


## classify local maxima
localMaximaTilePaths = list.files(localMaximaChmCmmPath, "\\.gpkg$", full.names = TRUE)
radiusTreetopsPath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/treetops/radius"

minimumHeightInM = 1

stands2016 = st_transform(st_read("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg", quiet = TRUE, layer = "unified stands 2016"),
                          make_compound_crs(6557, 8228)) %>% # keep in sync with same code in treetopJob.R
  select(standID2016)

treetopStartTime = Sys.time() # 7m21s @ 8 workers, 5m56s @ 12 workers, 9900X
with_progress({
  progressBar = progressor(steps = length(localMaximaTilePaths))
  radiusTreetops = bind_rows(future_map(localMaximaTilePaths, function(localMaximaTilePath)
  {
    treetopsDsm = st_join(st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaDsm") %>% select(tile, id, height, radius) %>%
                            filter(is.na(height) == FALSE, # exclude local maxima along north and west edges of DTM
                                   0.3048 * height >= minimumHeightInM, 
                                   0.3048 * radius > 0.42644864 + 0.03226514 * (0.3048 * height)^1.02543385),
                          stands2016, left = TRUE) %>%
      rename(treeID = id)
    
    treetopsChm = st_join(st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaChm") %>% select(tile, id, height, radius) %>%
                            filter(is.na(height) == FALSE,
                                   0.3048 * height >= minimumHeightInM, 
                                   0.3048 * radius > 0.51709010 + 0.11263928 * (0.3048 * height)^0.67990521),
                          stands2016, left = TRUE) %>%
      rename(treeID = id)
    
    treetopsCmm = st_join(st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaCmm") %>% 
                            mutate(is.na(height) == FALSE, # exclude local maxima along north and west edges of DTM
                                   dsmHeight = height,
                                   height = dsmHeight + cmmZ - dsmZ) %>%
                            select(tile, id, height, radius, dsmHeight) %>%
                            filter(is.na(height) == FALSE,
                                   0.3048 * height >= minimumHeightInM, 
                                   0.3048 * radius > -0.22250373 + 0.05153255 * (0.3048 * height)^0.93888880),
                                   #0.3048 * radius > 0.15370334 + 0.03378425 * (0.3048 * height)^0.97994540),
                          stands2016, left = TRUE) %>%
      rename(treeID = id)
    #tibble(dsm = nrow(treetopsDsm), chm = nrow(treetopsChm), cmm = nrow(treetopsCmm))
    
    tileFile = basename(localMaximaTilePath)
    treetopTilePath = file.path(radiusTreetopsPath, tileFile)
    st_write(treetopsDsm, treetopTilePath, layer = "treetopsDsm", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(treetopsChm, treetopTilePath, layer = "treetopsChm", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    st_write(treetopsCmm, treetopTilePath, layer = "treetopsCmm", delete_dsn = FALSE, delete_layer = TRUE, quiet = TRUE)
    
    progressBar()
    
    tileName = tools::file_path_sans_ext(tileFile)
    return(bind_rows(st_drop_geometry(treetopsDsm) %>% mutate(heightClassInM = round(0.3048 * height)) %>% group_by(standID2016, heightClassInM) %>%
                       summarize(method = "DSM radius", tile = tile[1], treetops = n(), .groups = "drop"),
                     st_drop_geometry(treetopsChm) %>% mutate(heightClassInM = round(0.3048 * height)) %>% group_by(standID2016, heightClassInM) %>%
                       summarize(method = "CHM radius", tile = tile[1], treetops = n(), .groups = "drop"),
                     st_drop_geometry(treetopsCmm) %>% mutate(heightClassInM = round(0.3048 * height)) %>% group_by(standID2016, heightClassInM) %>%
                       summarize(method = "CMM radius", tile = tile[1], treetops = n(), .groups = "drop")) %>%
            relocate(method, tile, standID2016, heightClassInM))
  }, .options = furrr_options(seed = TRUE))) # unclear why seed warnings occur as nothing here's randomized but it's handled just in case
})
Sys.time() - treetopStartTime
radiusTreetops %<>% mutate(method = factor(method, levels = c("DSM", "CHM", "CMM")))

#writexl::write_xlsx(radiusTreetops, file.path(radiusTreetopsPath, "standsByHeightClass.xlsx")) # 4.8 MB DSM only, 17.7 MB all surfaces

radiusTreetops %>% group_by(method) %>%
  summarize(tiles = length(unique(tile)), stands = length(unique(standID2016)), maxHeightInM = max(heightClassInM), 
            elliottTreetops1m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000), 0, treetops)),
            elliottTreetops5m = sum(if_else(is.na(standID2016) | (standID2016 >= 4000) | (heightClassInM < 5), 0, treetops)), 
            totalTreetops1m = sum(treetops), totalTreetops5m = sum(if_else(heightClassInM >= 5, treetops, 0)))

if (treetopOptions$includeInvestigatory)
{
  # full area distribution: all 561 tiles
  ggplot() +
    geom_histogram(aes(y = heightClassInM, weight = treetops, alpha = heightClassInM > 5), radiusTreetops %>% filter(surface == "DSM"), binwidth = 1, width = 1) +
    labs(x = "treetops", y = "height above ground, m", alpha = NULL, title = paste(plotLetters[1], "DSM radius")) +
  ggplot() +
    geom_histogram(aes(y = heightClassInM, weight = treetops, alpha = heightClassInM > 5), radiusTreetops %>% filter(surface == "CHM"), binwidth = 1, width = 1) +
    labs(x = "treetops", y = NULL, alpha = NULL, title = paste(plotLetters[2], "CHM radius")) +
  ggplot() +
    geom_histogram(aes(y = heightClassInM, weight = treetops, alpha = heightClassInM > 5), radiusTreetops %>% filter(surface == "CMM"), binwidth = 1, width = 1) +
    labs(x = "treetops", y = NULL, alpha = NULL, title = paste(plotLetters[3], "CMM radius")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() &
    guides(alpha = "none") &
    scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.5)) &
    coord_trans(x = scales::pseudo_log_trans(sigma = 10), xlim = c(0, 900000), ylim = c(0, 110)) &
    scale_x_continuous(breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), labels = scales::comma, minor_breaks = c(10 * 2:9, 100 * 2:9, 1000 * 2:9, 10000 * 2:9, 100000 * 2:9)) &
    scale_y_continuous(breaks = seq(0, 110, by = 10), expand = c(0, 1))
  
  # local maxima distribution on current tile
  dsm = st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaDsm")
  chm = st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaChm")
  cmm = st_read(localMaximaTilePath, quiet = TRUE, layer = "localMaximaCmm") %>% mutate(dsmHeight = height, height = dsmHeight + cmmZ - dsmZ)
  tibble(dsm = nrow(dsm), chm = nrow(chm), cmm = nrow(cmm))
  
  ggplot() +
    geom_histogram(aes(y = heightClassInM, alpha = heightClassInM > 5), dsm %>% mutate(heightClassInM = round(0.3048 * height)), binwidth = 1) +
    labs(x = "local maxima", y = "height above ground, m", alpha = NULL, title = paste(plotLetters[1], "DSM")) +
  ggplot() +
    geom_histogram(aes(y = heightClassInM, alpha = heightClassInM > 5), chm %>% mutate(heightClassInM = round(0.3048 * height)), binwidth = 1) +
    labs(x = "local maxima", y = NULL, alpha = NULL, title = paste(plotLetters[2], "CHM")) +
  ggplot() +
    geom_histogram(aes(y = heightClassInM, alpha = heightClassInM > 5), cmm %>% mutate(heightClassInM = round(0.3048 * height)), binwidth = 1) +
    labs(x = "local maxima", y = NULL, alpha = NULL, title = paste(plotLetters[3], "CMM")) +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() &
    guides(alpha = "none") &
    scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.5)) &
    scale_x_continuous(labels = scales::comma) &
    scale_y_continuous(breaks = seq(0, 100, by = 10))
}