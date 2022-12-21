library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(rsample)
library(stringr)

fitFunction = function(data)
{
  return(gsl_nls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, 
                 analysis(data), 
                 start = thplHeightFromDiameter$sharmaParton$m$getPars(), 
                 weights = dbhWeight))
}

fits = vfold_cv(thpl2016, v = 10, repeats = 10) %>% mutate(fit = map(splits, fitFunction))

fitParameters = bind_rows(pmap(fits, function(splits, id, id2, fit)
                          {
                            return(c(repetition = as.numeric(str_replace(id, "Repeat", "")),
                                     fold = as.numeric(str_replace(id2, "Fold", "")),
                                     fit$m$getPars()))
                          }))
fitStatistics = bind_rows(pmap(fits, function(splits, id, id2, fit)
                          {
                            validate = assessment(splits)
                            predictedHeight = predict(fit, validate)
                            residuals = validate$TotalHt - predictedHeight
                            return(tibble(repetition = as.numeric(str_replace(id, "Repeat", "")),
                                          fold = as.numeric(str_replace(id2, "Fold", "")),
                                          mae = mean(abs(residuals)), 
                                          rmse = sqrt(mean(residuals^2))))
                          }))
predictions = bind_rows(pmap(fits, function(splits, id, id2, fit)
                        {
                          return(tibble(repetition = as.numeric(str_replace(id, "Repeat", "")),
                                        fold = as.numeric(str_replace(id2, "Fold", "")),
                                        DBH = thpl2016$DBH, 
                                        predictedHeight = predict(fit, thpl2016)))
                        })) %>%
  mutate(groupID = str_c("repetition ", repetition, ", fold ", fold))
#validations = bind_rows(pmap(fits, function(splits, id, id2, fit) 
#                             {
#                               validate = assessment(splits)
#                               return(tibble(repetition = as.numeric(str_replace(id, "Repeat", "")),
#                                             fold = as.numeric(str_replace(id2, "Fold", "")),
#                                             DBH = validate$DBH, TotalHt = validate$TotalHt, 
#                                             predictedHeight = predict(fit, validate)))
#                             })) %>%
#  mutate(residual = TotalHt - predictedHeight)

ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), thpl2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = DBH, y = predictedHeight, group = groupID), predictions, alpha = 0.01, color = "firebrick") +
  guides(color = "none") +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))

ggplot() +
  geom_violin(aes(x = statistic, y = value, color = statistic), fitStatistics %>% pivot_longer(cols = c("mae", "rmse"), names_to = "statistic", values_to = "value"), draw_quantiles = c(0.25, 0.5, 0.75)) +
  guides(color = "none") +
  labs(x = NULL, y = "height error statistic value") +
  scale_x_discrete(breaks = c("mae", "rmse"), labels = c("MAE, m", "RMSE, m"))

ggplot() +
  geom_violin(aes(x = value, y = parameter, color = parameter), fitParameters %>% pivot_longer(cols = c("a1", "b1", "b1p", "b2", "b2p", "b3", "b4"), names_to = "parameter", values_to = "value"), draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
  guides(color = "none") +
  labs(x = "Sharma-Parton height coefficient value", y = NULL) +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50), minor_breaks = c(0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9, 30, 40), trans = scales::pseudo_log_trans(sigma = 0.03)) +
  scale_y_discrete(limits = rev)

# histogram version of coefficient violins: not as effective a presentation
#ggplot() +
#  geom_histogram(aes(x = value, y = after_stat(density), fill = parameter), fitParameters %>% pivot_longer(cols = c("a1", "b1", "b1p", "b2", "b2p", "b3", "b4"), names_to = "parameter", values_to = "value"), binwidth = 0.1) +
#  facet_wrap(vars(parameter)) +
#  guides(fill = "none") +
#  labs(x = "Sharma-Parton height coefficients", y = "probability density") +
#  scale_x_continuous(breaks = c(0, 1, 5, 10, 50), minor_breaks = c(2, 3, 4, 6, 7, 8, 9, 20, 30, 40), trans = scales::pseudo_log_trans(sigma = 0.03)) +
#  theme(strip.background = element_rect(fill = "grey95"))
  