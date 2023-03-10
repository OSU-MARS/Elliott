library(dplyr)
library(ggplot2)
library(furrr)
library(progressr)
library(rsample)
library(tidyr)
library(WeightedROC)

handlers(global = TRUE)
handlers("progress")
plan(multisession, workers = 1) # trivial unpack and return is fastest single threaded

iterations = 100
folds = 10
repetitions = 10

auc = rep(NA_real_, iterations)
with_progress({
  progressBar = progressor(steps = iterations)
  for (iteration in 1:iterations)
  {
    data = tibble(distribution1 = rt(folds^2 * repetitions, df = 10),
                  distribution2 = rt(folds^2 * repetitions, df = 10) - 0.2)
  
    samples = vfold_cv(data, v = folds, repeats = repetitions) %>% mutate(fit = future_map(splits, function(fold) 
      {
        return(mutate(assessment(fold), distribution1 = mean(distribution1), distribution2 = mean(distribution2))) 
      }))
    observations = bind_rows(samples$fit) %>% 
      pivot_longer(cols = c("distribution1", "distribution2"), names_to = "distribution", values_to = "sample") %>%
      mutate(distribution = as.factor(if_else(distribution == "distribution1", 0, 1)))
    # AUC estimates the probability a sample from distribution 2 (distribution = 1) is greater than one from 
    # distribution 1 (distribution = 0)
    # The error rates of binary classification---deciding which distribution a sample is from based on its value---
    # are empirical probability measures. See, for example, https://mlu-explain.github.io/roc-auc/.
    roc = WeightedROC(guess = observations$sample, label = observations$distribution)
    auc[iteration] = WeightedAUC(roc)
    
    progressBar(message = sprintf("AUC = %.3f", auc[iteration]))
  }
})

ggplot() +
  geom_histogram(aes(x = auc), binwidth = 0.01) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "AUC", y = paste0(folds, "x", repetitions))

quantile(auc, probs = c(0.015, 0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975, 0.995))

# generalization to more complex tibbles
#foo %>% group_by(key1, key2, key3) %>%
#  group_modify(~{
#    fooSubset = foo %>% filter(key1 == .y$key1, key2 == .y$key2, key3 == "special important value of key3")
#    return(tibble(auc = WeightedAUC(WeightedROC(guess = c(.x$bar, fooSubset$bar),
#                                                label = factor(c(rep(0, nrow(.x)), rep(1, nrow(fooSubset))), levels = c(0, 1))))))
#  })

