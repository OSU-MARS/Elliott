library(rsample)
library(tidyverse)

foo = vfold_cv(thpl2016) %>%
  mutate(fit = map(splits, function(split, ...) { return(gsl_nls(TotalHt ~ 1.37 + a1*topHeight^(b1 + b1p * isPlantation)*(1 - exp((b2 + b2p * isPlantation)*(tph/standBasalAreaPerHectare)^b3*DBH))^b4, assessment(split), start = thplHeightFromDiameterSharmaParton$m$getPars(), weights = pmin(DBH^-1.2, 1))) }))

foo %>% pmap_dfr(~tibble(fold = ..2, TotalHt = predict(..3, thpl2016)))
bar = map2_dfr(foo$id, foo$fit, function(id, fit) { return(tibble(fold = id, DBH = thpl2016$DBH, TotalHt = predict(fit, thpl2016))) })

ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), thpl2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = DBH, y = TotalHt, color = fold, group = fold), bar, alpha = 0.5) +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))
