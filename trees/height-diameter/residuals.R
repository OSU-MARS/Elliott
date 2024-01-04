## GAM-based variance estimation
# Using by = interaction(speciesGroup, isPlantation) takes ~61 seconds for height and 160 for DBH (Zen 3, 4.8 GHz).
# furring by species group takes ~1 s but means seven models to analyze.
plan(multisession, workers = 7)
liveUnbrokenHeightTrees2016 = trees2016 %>% filter(isLiveUnbroken, is.na(TotalHt) == FALSE)

#baseFormsForHeight = liveUnbrokenHeightTrees2016 %>% split(liveUnbrokenHeightTrees2016$speciesGroup) %>%
#  future_map(~gam(TotalHt ~ s(DBH, bs = "ts", by = as.factor(isPlantation), k = 15, pc = c(DBH = -2)), data = ., method = "REML", select = TRUE, nthreads = 4))
#baseFormsForDbh = liveUnbrokenHeightTrees2016 %>% split(liveUnbrokenHeightTrees2016$speciesGroup) %>%
#  future_map(~gam(DBH ~ s(TotalHt, bs = "ts", by = as.factor(isPlantation), k = 15, pc = c(TotalHt = 1.37)), data = liveUnbrokenHeightTrees2016, method = "REML", select = TRUE, nthreads = 8))
baseFormForHeight = gam(TotalHt ~ s(DBH, bs = "ts", by = interaction(speciesGroup, as.factor(isPlantation)), k = 16, pc = c(DBH = -2)), data = liveUnbrokenHeightTrees2016, method = "REML", select = TRUE, weights = pmin(TreeCount/(1.2*DBH^0.8), 5*TreeCount), nthreads = 8)
baseFormForDbh = gam(DBH ~ s(TotalHt, bs = "ts", by = interaction(speciesGroup, as.factor(isPlantation)), k = 12, pc = c(TotalHt = 1.37)), data = liveUnbrokenHeightTrees2016, method = "REML", select = TRUE, weights = pmin(TreeCount/(7.2*(TotalHt - 1.3)^1.3), 5*TreeCount), nthreads = 8)
#summary(baseFormForHeight)
#summary(baseFormForDbh)

predictionError = tibble(speciesGroup = liveUnbrokenHeightTrees2016$speciesGroup,
                         TreeCount = liveUnbrokenHeightTrees2016$TreeCount,
                         TotalHt = liveUnbrokenHeightTrees2016$TotalHt, 
                         DBH = liveUnbrokenHeightTrees2016$DBH, 
                         isPlantation = liveUnbrokenHeightTrees2016$isPlantation,
                         predictedHeight = predict(baseFormForHeight, liveUnbrokenHeightTrees2016), # a couple seconds
                         predictedDbh = predict(baseFormForDbh, liveUnbrokenHeightTrees2016)) %>%
  mutate(speciesGroup = recode_factor(speciesGroup, "DF" = "Psme", "RA" = "Alru", "WH" = "Tshe", "BM" = "Acma", "OM" = "Umca", "RC" = "Thpl"),
         heightResidual = predictedHeight - TotalHt, dbhResidual = predictedDbh - DBH,
         heightResidualSquared = heightResidual^2, dbhResidualSquared = dbhResidual^2)

empiricalHeightVariance = predictionError %>% 
  mutate(dbhClass = if_else(DBH < 100, 5 * floor(DBH / 5) + 0.5 * 5, 25 * floor(DBH / 25) + 0.5 * 25)) %>%
  group_by(speciesGroup, isPlantation, dbhClass) %>%
  summarize(n = n(),
            variance = if_else(n() > 1, sum(heightResidualSquared) / (n() - 1), NA_real_), .groups = "drop") %>%
  mutate(varianceQ025 = (n - 1) * variance / qchisq(0.975, n - 1),
         varianceQ975 = (n - 1) * variance / qchisq(0.025, n - 1)) %>%
  relocate(speciesGroup, isPlantation, dbhClass, n, varianceQ025, variance, varianceQ975)
empiricalDbhVariance = predictionError %>% 
  mutate(heightClass = if_else(TotalHt < 70, floor(TotalHt) + 0.5 * 1, 2.5 * floor(TotalHt / 2.5) + 0.5 * 2.5)) %>%
  group_by(speciesGroup, isPlantation, heightClass) %>%
  summarize(n = n(),
            variance = if_else(n() > 1, sum(dbhResidualSquared) / (n() - 1), NA_real_), .groups = "drop") %>%
  mutate(varianceQ025 = (n - 1) * variance / qchisq(0.975, n - 1),
         varianceQ975 = (n - 1) * variance / qchisq(0.025, n - 1)) %>%
  relocate(speciesGroup, isPlantation, heightClass, n, varianceQ025, variance, varianceQ975)

varianceForHeightWeightMultiplier = list(psme = 1.8, alru = 1.6, tshe = 1.3, acma = 1.3, umca = 1.3, thpl = 0.2, other = 1.1)
varianceForHeightWeightPower = list(psme = 0.7, alru = 0.8, tshe = 0.8, acma = 0.8, umca = 0.7, thpl = 1.2, other = 0.9)
varianceForHeight = list(psme = gsl_nls(heightResidualSquared ~ a1*DBH^(b1 + b1p*isPlantation), predictionError %>% filter(speciesGroup == "Psme"), start = list(a1 = 1, b1 = 1, b1p = 0), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$psme*DBH^varianceForHeightWeightPower$psme), 5*TreeCount))) # a1p statistically significant but fit is implausible on a1-b1 divergence and cancellation
varianceForHeight$alru = gsl_nls(heightResidualSquared ~ a1*DBH^(b1 + b1p*isPlantation), predictionError %>% filter(speciesGroup == "Alru"), start = list(a1 = 1, b1 = 1, b1p = 0), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$alru*DBH^varianceForHeightWeightPower$alru), 5*TreeCount))
varianceForHeight$tshe = gsl_nls(heightResidualSquared ~ a1*DBH^b1, predictionError %>% filter(speciesGroup == "Tshe"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$tshe*DBH^varianceForHeightWeightPower$tshe), 5*TreeCount))
varianceForHeight$acma = gsl_nls(heightResidualSquared ~ a1*DBH^b1, predictionError %>% filter(speciesGroup == "Acma"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$acma*DBH^varianceForHeightWeightPower$acma), 5*TreeCount))
varianceForHeight$umca = gsl_nls(heightResidualSquared ~ a1*DBH^b1, predictionError %>% filter(speciesGroup == "Umca"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$umca*DBH^varianceForHeightWeightPower$umca), 5*TreeCount))
varianceForHeight$thpl = gsl_nls(heightResidualSquared ~ a1*DBH^b1, predictionError %>% filter(speciesGroup == "Thpl"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$thpl*DBH^varianceForHeightWeightPower$thpl), 5*TreeCount))
varianceForHeight$other = gsl_nls(heightResidualSquared ~ a1*DBH^b1, predictionError %>% filter(speciesGroup == "other"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForHeightWeightMultiplier$other*DBH^varianceForHeightWeightPower$other), 5*TreeCount))
bind_rows(imap(varianceForHeight, function(model, name) { return(c(responseVariable = "height", name = name, converged = model$convInfo$isConv, round(model$m$getPars(), 3))) })) %>% mutate(a1 = round(as.numeric(a1), 2), b1 = round(as.numeric(b1), 2), b1p = round(as.numeric(b1p), 3))
#bind_rows(imap(varianceForHeight, function(model, name) { return(c(responseVariable = "height", name = name, converged = model$convInfo$isConv, round(model$m$getPars(), 3))) })) %>% summarize(a1 = mean(as.numeric(a1)), b1 = mean(as.numeric(b1)))
#lapply(varianceForHeight, confint2, level = 0.99)

varianceForDbhWeightMultiplier = list(psme = 5.7, alru = 27, tshe = 1.1, acma = 9.7, umca = 3.5, thpl = 2.4, other = 1.3)
varianceForDbhWeightPower = list(psme = 1.1, alru = 0.6, tshe = 1.5, acma = 1.2, umca = 1.6, thpl = 1.4, other = 1.8)
varianceForDbh = list(psme = gsl_nls(dbhResidualSquared ~ (a1 + a1p*isPlantation)*(TotalHt - 1.37)^(b1 + b1p*isPlantation), predictionError %>% filter(speciesGroup == "Psme"), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$psme*(TotalHt - 1.37)^varianceForDbhWeightPower$psme), 5*TreeCount)))
varianceForDbh$alru = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^(b1 + b1p*isPlantation), predictionError %>% filter(speciesGroup == "Alru"), start = list(a1 = 1, b1 = 1, b1p = 0), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$alru*(TotalHt - 1.37)^varianceForDbhWeightPower$alru), 5*TreeCount))
varianceForDbh$tshe = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^b1, predictionError %>% filter(speciesGroup == "Tshe"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$tshe*(TotalHt - 1.37)^varianceForDbhWeightPower$tshe), 5*TreeCount)) # b1p makes a1 not significant
varianceForDbh$acma = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^b1, predictionError %>% filter(speciesGroup == "Acma"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$acma*(TotalHt - 1.37)^varianceForDbhWeightPower$acma), 5*TreeCount))
varianceForDbh$umca = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^b1, predictionError %>% filter(speciesGroup == "Umca"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$umca*(TotalHt - 1.37)^varianceForDbhWeightPower$umca), 5*TreeCount))
varianceForDbh$thpl = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^b1, predictionError %>% filter(speciesGroup == "Thpl"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$thpl*(TotalHt - 1.37)^varianceForDbhWeightPower$thpl), 5*TreeCount))
varianceForDbh$other = gsl_nls(dbhResidualSquared ~ a1*(TotalHt - 1.37)^b1, predictionError %>% filter(speciesGroup == "other"), start = list(a1 = 1, b1 = 1), weight = pmin(TreeCount/(varianceForDbhWeightMultiplier$other*(TotalHt - 1.37)^varianceForDbhWeightPower$other), 5*TreeCount))
bind_rows(imap(varianceForDbh, function(model, name) { return(c(responseVariable = "DBH", name = name, converged = model$convInfo$isConv, round(model$m$getPars(), 3))) })) %>% mutate(a1 = round(as.numeric(a1), 2), b1 = round(as.numeric(b1), 2), b1p = round(as.numeric(b1p), 2))
#bind_rows(imap(varianceForDbh, function(model, name) { return(c(responseVariable = "DBH", name = name, converged = model$convInfo$isConv, round(model$m$getPars(), 3))) })) %>% summarize(a1 = mean(as.numeric(a1)), b1 = mean(as.numeric(b1)))
#lapply(varianceForDbh, confint2, level = 0.99)

predict_gsl_nls = function(model, newData, weights, level = 0.95)
{
  # modification of gslnls::predict.gsl_nls() to include weights (https://github.com/JorisChau/gslnls/blob/master/R/nls_methods.R)
  fit = model$m$predict(newData)
  Fdot = diag(pmin(weights, 1)) %*% model$m$gradient1(newData)
  Rmat = qr.R(qr(Fdot))
  scale = sigma(model)
  a = c((1 - level) / 2, (1 + level) / 2)
  ses = scale * sqrt(1 + rowSums(Fdot %*% chol2inv(Rmat) * Fdot))
  ci = fit + ses %o% qt(a, df.residual(model))
  return(cbind(fit = fit, lwr = ci[, 1], upr = ci[, 2]))
}

modeledVariance = crossing(tibble(DBH = seq(0, 250), TotalHt = seq(1.37, 85, length.out = 251)),
                           tibble(isPlantation = c(TRUE, FALSE))) %>%
  mutate(estHeightVarPsme = predict_gsl_nls(varianceForHeight$psme, ., pmin(1/(varianceForHeightWeightMultiplier$psme*.$DBH^varianceForHeightWeightPower$psme), 5)), # column order is fit, lower, upper
         estHeightVarAlru = predict_gsl_nls(varianceForHeight$alru, ., pmin(1/(varianceForHeightWeightMultiplier$alru*.$DBH^varianceForHeightWeightPower$alru), 5)),
         estHeightVarTshe = predict_gsl_nls(varianceForHeight$tshe, ., pmin(1/(varianceForHeightWeightMultiplier$tshe*.$DBH^varianceForHeightWeightPower$tshe), 5)),
         estHeightVarAcma = predict_gsl_nls(varianceForHeight$acma, ., pmin(1/(varianceForHeightWeightMultiplier$acma*.$DBH^varianceForHeightWeightPower$acma), 5)),
         estHeightVarUmca = predict_gsl_nls(varianceForHeight$umca, ., pmin(1/(varianceForHeightWeightMultiplier$umca*.$DBH^varianceForHeightWeightPower$umca), 5)),
         estHeightVarThpl = predict_gsl_nls(varianceForHeight$thpl, ., pmin(1/(varianceForHeightWeightMultiplier$thpl*.$DBH^varianceForHeightWeightPower$thpl), 5)),
         estHeightVarOther = predict_gsl_nls(varianceForHeight$other, ., pmin(1/(varianceForHeightWeightMultiplier$other*.$DBH^varianceForHeightWeightPower$other), 5)),
         estDbhVarPsme = predict_gsl_nls(varianceForDbh$psme, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$psme), 5)),
         estDbhVarAlru = predict_gsl_nls(varianceForDbh$alru, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$alru), 5)),
         estDbhVarTshe = predict_gsl_nls(varianceForDbh$tshe, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$tshe), 5)),
         estDbhVarAcma = predict_gsl_nls(varianceForDbh$acma, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$acma), 5)),
         estDbhVarUmca = predict_gsl_nls(varianceForDbh$umca, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$umca), 5)),
         estDbhVarThpl = predict_gsl_nls(varianceForDbh$thpl, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$thpl), 5)),
         estDbhVarOther = predict_gsl_nls(varianceForDbh$other, ., pmin(1/(varianceForDbhWeightMultiplier$psme*(.$TotalHt - 1.37)^varianceForDbhWeightPower$other), 5))) %>%
  pivot_longer(cols = c(starts_with("estHeightVar"), starts_with("estDbhVar")), names_to = c("estimate", "speciesGroup"), names_pattern = "(estHeightVar|estDbhVar)(\\w{4,5})") %>%
  filter((speciesGroup == "Psme") |
         ((speciesGroup == "Alru") & (((estimate == "estHeightVar") & (DBH < 125)) | ((estimate == "estDbhVar") & (TotalHt < 50)))) |
         ((speciesGroup == "Tshe") & (((estimate == "estHeightVar") & (DBH < 175)) | ((estimate == "estDbhVar") & (TotalHt < 75)))) |
         ((speciesGroup == "Acma") & (((estimate == "estHeightVar") & (DBH < 175)) | ((estimate == "estDbhVar") & (TotalHt < 50)))) |
         ((speciesGroup == "Umca") & (((estimate == "estHeightVar") & (DBH < 125)) | ((estimate == "estDbhVar") & (TotalHt < 40)))) |
         ((speciesGroup == "Thpl") & (((estimate == "estHeightVar") & (DBH < 200)) | ((estimate == "estDbhVar") & (TotalHt < 70)))) |
         ((speciesGroup == "Other") & (((estimate == "estHeightVar") & (DBH < 175)) | ((estimate == "estDbhVar") & (TotalHt < 60))))) %>%
  pivot_wider(id_cols = c("DBH", "TotalHt", "speciesGroup", "isPlantation"), names_from = "estimate", values_from = "value") %>%
  mutate(speciesGroup = factor(if_else(speciesGroup != "Other", speciesGroup, "other"), levels = levels(predictionError$speciesGroup)))
#print(modeledVariance, n = 100)

## Figure S4: residual variance model
speciesGroupColors = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")
speciesGroupNames = c("Douglas-fir", "red alder", "western hemlock", "bigleaf maple", "Oregon myrtle", "western redcedar", "other species")
ggplot() +
  geom_errorbar(aes(x = dbhClass, ymin = varianceQ025, ymax = varianceQ975, alpha = n), empiricalHeightVariance, color = "grey50", linewidth = 0.3, na.rm = TRUE) +
  geom_point(aes(x = dbhClass, y = variance, alpha = n, shape = isPlantation, size = isPlantation), empiricalHeightVariance, color = "grey40", na.rm = TRUE) +
  geom_ribbon(aes(x = DBH, ymin = pmax(estHeightVar[,2], 0), ymax = estHeightVar[,3], fill = speciesGroup, group = paste(speciesGroup, isPlantation), linetype = isPlantation), modeledVariance, alpha = 0.1) +
  geom_path(aes(x = DBH, y = estHeightVar[,1], color = speciesGroup, group = paste(speciesGroup, isPlantation), linetype = isPlantation), modeledVariance, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 220), ylim = c(0, 125)) +
  facet_wrap(vars(speciesGroup), labeller = labeller(speciesGroup = c("Psme" = "Douglas-fir", "Alru" = "red alder", "Tshe" = "western hemlock", "Acma" = "bigleaf maple", "Umca" = "Oregon myrtle", "Thpl" = "western redcedar", "other" = "other species"))) +
  labs(x = "DBH, cm", y = bquote(widehat(var)["height"]*", m²"), alpha = "variance\nestimate", color = "variance\nmodel", fill = "variance\nmodel", linetype = NULL, shape = NULL, size = NULL, title = "a) height prediction") +
  scale_x_continuous(breaks = seq(0, 300, by = 100)) +
ggplot() +
  geom_errorbar(aes(xmin = varianceQ025, xmax = varianceQ975, y = heightClass, alpha = n), empiricalDbhVariance, color = "grey50", linewidth = 0.3, na.rm = TRUE, orientation = "y") +
  geom_point(aes(x = variance, y = heightClass, alpha = n, shape = isPlantation, size = isPlantation), empiricalDbhVariance, color = "grey40", na.rm = TRUE) +
  geom_ribbon(aes(xmin = pmax(estDbhVar[,2], 0), xmax = estDbhVar[,3], y = TotalHt, fill = speciesGroup, group = paste(speciesGroup, isPlantation), linetype = isPlantation), modeledVariance, alpha = 0.1) +
  geom_path(aes(x = estDbhVar[,1], y = TotalHt, color = speciesGroup, group = paste(speciesGroup, isPlantation), linetype = isPlantation), modeledVariance, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1550), ylim = c(0, 83)) +
  facet_wrap(vars(speciesGroup), labeller = labeller(speciesGroup = c("Psme" = "Douglas-fir", "Alru" = "red alder", "Tshe" = "western hemlock", "Acma" = "bigleaf maple", "Umca" = "Oregon myrtle", "Thpl" = "western redcedar", "other" = "other species"))) +
  labs(x = bquote(widehat(var)["DBH"]*", cm²"), y = "height, m", alpha = "variance\nestimate", color = "variance\nmodel", fill = "variance\nmodel", linetype = NULL, shape = NULL, size = NULL, title = "b) DBH prediction") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  guides(alpha = guide_legend(order = 3, reverse = TRUE), color = guide_legend(order = 1), fill = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 2), size = guide_legend(order = 2)) &
  scale_alpha_continuous(breaks = c(2, 10, 100, 1000), labels = c("nᵢ = 2", "nᵢ = 10", "nᵢ = 100", "nᵢ = 1000"), limits = c(2, 1000), range = c(0.1, 1.0), trans = "log10") &
  scale_color_manual(breaks = levels(modeledVariance$speciesGroup), labels = speciesGroupNames, limits = levels(modeledVariance$speciesGroup), values = speciesGroupColors) &
  scale_fill_manual(breaks = levels(modeledVariance$speciesGroup), labels = speciesGroupNames, limits = levels(modeledVariance$speciesGroup), values = speciesGroupColors) &
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural\nregeneration", "plantation\n(if significant)"), values = c("solid", "longdash")) &
  scale_shape_manual(breaks = c(FALSE, TRUE), labels = c("natural\nregeneration", "plantation\n(if significant)"), values = c(16, 18)) &
  scale_size_manual(breaks = c(FALSE, TRUE), labels = c("natural\nregeneration", "plantation\n(if significant)"), values = c(1.5, 1.9)) &
  theme(legend.spacing.y = unit(0.5, "line"), strip.background = element_rect(fill = "grey95"), strip.text = element_text(size = 8))
#ggsave("trees/height-diameter/figures/Figure S04 error model gsl_nls.png", height = 15, width = 22, units = "cm", dpi = 200)


## height interquartile ranges
library(quantreg)
psmeHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 70.4, a1p = -25.3, a2 = 660, a2p = -348, b1 = 1.61), tau = 0.25)
psmeHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), psme2016, start = list(a1 = 93.9, a1p = -10.3, a2 = 108, a2p = 4.8, b1 = 1.13), tau = 0.75)
#psmeHeightFromDiameterSharmaPartonQ75 = nlrq(TotalHt ~ 1.37 + a1*topHeight^a2*(1 - exp(b1*(tph/standBasalAreaPerHectare)^b2*DBH))^b3, psme2016, start = list(a1 = 52.41, a2 = 0.135, b1 = -0.017, b2 = -0.071, b3 = 1.02), tau = 0.75)
#psmeHeightFromDiameterSharmaPartonQ75 = nlrq(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 52.41, a2 = 0.135, a2p = -0, b1 = -0.017, b1p = -0, b2 = -0.071, b2p = -0, b3 = 1.02, b3p = -0), tau = 0.75)
psmeHeightFromDiameterSharmaPartonQ25 = nlrq(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 42.8, a2 = 0.114, a2p = -0.045, b1 = -0.023, b1p = -0.026, b2 = 0.039, b2p = -0.258, b3 = 1.94, b3p = -0.94), tau = 0.25)
psmeHeightFromDiameterSharmaPartonQ75 = nlrq(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp((b1 + b1p * isPlantation)*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), psme2016, start = list(a1 = 57.5, a2 = 0.065, a2p = 0.071, b1 = -0.021, b1p = 0.010, b2 = 0.028, b2p = -0.115, b3 = 1.39, b3p = -0.573), tau = 0.75)

alruHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 22.7, a1p = 5.55, a2 = 101, b1 = 1.54, b1p = -0.067), tau = 0.25)
alruHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^(b1 + b1p * isPlantation) / (a2 + + DBH^(b1 + b1p * isPlantation)), alru2016, start = list(a1 = 31.6, a1p = 5.53, a2 = 49.6, b1 = 1.37, b1p = -0.099), tau = 0.75)
alruHeightFromDiameterSharmaPartonQ25 = nlrq(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 21.3, a2 = 0.038, a2p = 0.068, b1 = -0.039, b2 = 0.084, b2p = -0.127, b3 = 1.19, b3p = -0.14), tau = 0.25)
alruHeightFromDiameterSharmaPartonQ75 = nlrq(TotalHt ~ 1.37 + a1*topHeight^(a2 + a2p * isPlantation)*(1 - exp(b1*(tph/standBasalAreaPerHectare)^(b2 + b2p * isPlantation)*DBH))^(b3 + b3p * isPlantation), alru2016, start = list(a1 = 21.3, a2 = 0.038, a2p = 0.068, b1 = -0.039, b2 = 0.084, b2p = -0.127, b3 = 1.19, b3p = -0.14), tau = 0.75)

tsheHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264), tau = 0.25)
tsheHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), tshe2016, start = list(a1 = 74.8, a1p = -19.0, a2 = 200, a2p = -77.4, b1 = 1.264), tau = 0.75)

acmaHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), tau = 0.25)
acmaHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), acma2016, start = list(a1 = 41.2, a2 = 49.0, a2p = -9.29, b1 = 0.986), tau = 0.75)

umcaHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27), tau = 0.25)
umcaHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + (a1 + a1p*isPlantation)*DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), umca2016, start = list(a1 = 21.4, a1p = -4.37, a2 = 43.8, a2p = -20.0, b1 = 1.27), tau = 0.75)

thplHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176), tau = 0.25)
thplHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + (a1 + a1p * isPlantation) * DBH^b1 / (a2 + a2p * isPlantation + DBH^b1), thpl2016, start = list(a1 = 70.3, a1p = -18.7, a2 = 200, a2p = -68.2, b1 = 1.176), tau = 0.75)

otherHeightFromDiameterMichaelisMentenQ25 = nlrq(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 34.6, a2 = 40.2, b1 = 1.03), tau = 0.25)
otherHeightFromDiameterMichaelisMentenQ75 = nlrq(TotalHt ~ 1.37 + a1*DBH^b1 / (a2 + DBH^b1), other2016, start = list(a1 = 34.6, a2 = 40.2, b1 = 1.03), tau = 0.75)

heightIqr = bind_rows(tibble(species = "PSME", name = "Michaelis-Menten", fitting = "nlrq", DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, q25 = predict(psmeHeightFromDiameterMichaelisMentenQ25), q75 = predict(psmeHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "PSME", name = "Sharma-Parton", fitting = "nlrq", DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, q25 = predict(psmeHeightFromDiameterSharmaPartonQ25), q75 = predict(psmeHeightFromDiameterSharmaPartonQ75)),
                      tibble(species = "ALRU2", name = "Michaelis-Menten", fitting = "nlrq", DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, q25 = predict(alruHeightFromDiameterMichaelisMentenQ25), q75 = predict(alruHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "ALRU2", name = "Sharma-Parton", fitting = "nlrq", DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, q25 = predict(alruHeightFromDiameterSharmaPartonQ25), q75 = predict(alruHeightFromDiameterSharmaPartonQ75)),
                      tibble(species = "TSHE", name = "Michaelis-Menten", fitting = "nlrq", DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, q25 = predict(tsheHeightFromDiameterMichaelisMentenQ25), q75 = predict(tsheHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "ACMA3", name = "Michaelis-Menten", fitting = "nlrq", DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, q25 = predict(acmaHeightFromDiameterMichaelisMentenQ25), q75 = predict(acmaHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "UMCA", name = "Michaelis-Menten", fitting = "nlrq", DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, q25 = predict(umcaHeightFromDiameterMichaelisMentenQ25), q75 = predict(umcaHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "THPL", name = "Michaelis-Menten", fitting = "nlrq", DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, q25 = predict(thplHeightFromDiameterMichaelisMentenQ25), q75 = predict(thplHeightFromDiameterMichaelisMentenQ75)),
                      tibble(species = "other", name = "Michaelis-Menten", fitting = "nlrq", DBH = other2016$DBH, isPlantation = other2016$isPlantation, q25 = predict(otherHeightFromDiameterMichaelisMentenQ25), q75 = predict(otherHeightFromDiameterMichaelisMentenQ75))) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"),  levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         iqr = q75 - q25)

ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), psme2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterSharmaPartonQ25), color = "q25", group = psme2016$isPlantation, linetype = psme2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterSharmaPartonQ75), color = "q75", group = psme2016$isPlantation, linetype = psme2016$isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 83, label = "a) Douglas-fir", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250)) +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL, linetype = NULL) +
  scale_color_discrete(breaks = c("q25", "q75"), labels = c("lower quartile", "upper quartile")) +
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c("solid", "dashed")) +
  theme(legend.justification = c(0, 1), legend.position = "none") +
ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), alru2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterMichaelisMentenQ25), color = "q25", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterMichaelisMentenQ75), color = "q75", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  #geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterSharmaPartonQ25), color = "q25", group = alru2016$isPlantation, linetype = alru2016$isPlantation), alpha = 0.5) +
  #geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterSharmaPartonQ75), color = "q75", group = alru2016$isPlantation, linetype = alru2016$isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 83, label = "b) red alder", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 120)) +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL, linetype = NULL) +
  scale_color_discrete(breaks = c("q25", "q75"), labels = c("lower quartile", "upper quartile")) +
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c("solid", "dashed")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
plot_layout(nrow = 2, widths = c(250, 150))


## diameter interquartile range
psmeDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), tau = 0.25)
psmeDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), psme2016, start = list(a1 = 2.67, b1 = 0.813, b1p = -0.126, b2 = 0.0067, b2p = 0.0096), tau = 0.75)
alruDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.35, b1 = 1.37, b1p = -0.24, b2 = -0.033, b2p = 0.022), tau = 0.25)
alruDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), alru2016, start = list(a1 = 1.35, b1 = 1.37, b1p = -0.24, b2 = -0.033, b2p = 0.022), tau = 0.75)
tsheDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.31, b1 = 0.818, b2 = 0.010), tau = 0.25)
tsheDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^b1 * exp(b2 * (TotalHt - 1.37)), tshe2016, start = list(a1 = 2.31, b1 = 0.818, b2 = 0.010), tau = 0.75)
acmaDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), tau = 0.25)
acmaDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), acma2016, start = list(a1 = 1.20, b1 = 1.52, b1p = -0.32, b2 = -0.038, b2p = 0.037), tau = 0.75)
umcaDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 2.48, b1 = 1.14, b1p = -0.57, b2 = -0.019, b2p = 0.092), tau = 0.25)
umcaDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), umca2016, start = list(a1 = 2.48, b1 = 1.14, b1p = -0.57, b2 = -0.019, b2p = 0.092), tau = 0.75)
thplDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^b1 * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.12, b1 = 1.01, b2 = 0.0038, b2p = 0.0025), tau = 0.25)
thplDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^b1 * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), thpl2016, start = list(a1 = 2.12, b1 = 1.01, b2 = 0.0038, b2p = 0.0025), tau = 0.75)
otherDiameterFromHeightRuarkQ25 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 2.54, b1 = 0.56, b1p = -0.28, b2 = 0.040, b2p = 0.043), tau = 0.25)
otherDiameterFromHeightRuarkQ75 = nlrq(DBH ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation) * exp((b2 + b2p * isPlantation) * (TotalHt - 1.37)), other2016, start = list(a1 = 2.54, b1 = 0.56, b1p = -0.28, b2 = 0.040, b2p = 0.043), tau = 0.75)


diameterIqr = bind_rows(tibble(species = "PSME", name = "Ruark", fitting = "nlrq", TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, q25 = predict(psmeDiameterFromHeightRuarkQ25), q75 = predict(psmeDiameterFromHeightRuarkQ75)),
                        tibble(species = "ALRU2", name = "Ruark", fitting = "nlrq", TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, q25 = predict(alruDiameterFromHeightRuarkQ25), q75 = predict(alruDiameterFromHeightRuarkQ75)),
                        tibble(species = "TSHE", name = "Ruark", fitting = "nlrq", TotalHt = tshe2016$TotalHt, isPlantation = tshe2016$isPlantation, q25 = predict(tsheDiameterFromHeightRuarkQ25), q75 = predict(tsheDiameterFromHeightRuarkQ75)),
                        tibble(species = "ACMA3", name = "Ruark", fitting = "nlrq", TotalHt = acma2016$TotalHt, isPlantation = acma2016$isPlantation, q25 = predict(acmaDiameterFromHeightRuarkQ25), q75 = predict(acmaDiameterFromHeightRuarkQ75)),
                        tibble(species = "UMCA", name = "Ruark", fitting = "nlrq", TotalHt = umca2016$TotalHt, isPlantation = umca2016$isPlantation, q25 = predict(umcaDiameterFromHeightRuarkQ25), q75 = predict(umcaDiameterFromHeightRuarkQ75)),
                        tibble(species = "THPL", name = "Ruark", fitting = "nlrq", TotalHt = thpl2016$TotalHt, isPlantation = thpl2016$isPlantation, q25 = predict(thplDiameterFromHeightRuarkQ25), q75 = predict(thplDiameterFromHeightRuarkQ75)),
                        tibble(species = "other", name = "Ruark", fitting = "nlrq", TotalHt = other2016$TotalHt, isPlantation = other2016$isPlantation, q25 = predict(otherDiameterFromHeightRuarkQ25), q75 = predict(otherDiameterFromHeightRuarkQ75))) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"),  levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         iqr = q75 - q25)

ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), psme2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = predict(psmeDiameterFromHeightRuarkQ25), y = psme2016$TotalHt, color = "q25", group = psme2016$isPlantation, linetype = psme2016$isPlantation), alpha = 0.5) +
  geom_line(aes(x = predict(psmeDiameterFromHeightRuarkQ75), y = psme2016$TotalHt, color = "q75", group = psme2016$isPlantation, linetype = psme2016$isPlantation), alpha = 0.5) +
  annotate("text", x = 0, y = 83, label = "a) Douglas-fir", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250)) +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL, linetype = NULL) +
  scale_color_discrete(breaks = c("q25", "q75"), labels = c("lower quartile", "upper quartile")) +
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c("solid", "dashed")) +
  theme(legend.justification = c(0, 1), legend.position = "none") +
ggplot() +
  geom_point(aes(x = DBH, y = TotalHt), alru2016, alpha = 0.1, color = "grey25", shape = 16) +
  geom_line(aes(x = predict(alruDiameterFromHeightRuarkQ25), y = alru2016$TotalHt, color = "q25", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = predict(alruDiameterFromHeightRuarkQ75), y = alru2016$TotalHt, color = "q75", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  annotate("text", x = 0, y = 83, label = "b) red alder", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 120)) +
  labs(x = "DBH, cm", y = "height of unbroken stem, m", color = NULL, linetype = NULL) +
  scale_color_discrete(breaks = c("q25", "q75"), labels = c("lower quartile", "upper quartile")) +
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c("solid", "dashed")) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
plot_layout(nrow = 2, widths = c(250, 150))


## direct fitting of height residuals
#acmaHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#acmaHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#tsheHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#tsheHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#umcaHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#umcaHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#thplHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#thplHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p*isPlantation)*DBH^(b1 + b1p*isPlantation), tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
#confint_nlrob(thplHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, nrow(thpl2016)))
psmeHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
alruHeightChapmanRichardsResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightFromDiameterChapmanRichards))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
tsheHeightChapmanRichardsResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterChapmanRichards))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
acmaHeightChapmanRichardsResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterChapmanRichards))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
umcaHeightChapmanRichardsResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterChapmanRichards))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
thplHeightChapmanRichardsResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterChapmanRichards))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
otherHeightChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = other2016$DBH, isPlantation = other2016$isPlantation, residual = abs(residuals(otherHeightFromDiameterChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))

psmeHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
alruHeightMichaelisMentenResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
tsheHeightMichaelisMentenResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
acmaHeightMichaelisMentenResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
umcaHeightMichaelisMentenResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
thplHeightMichaelisMentenResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
otherHeightMichaelisMentenResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = other2016$DBH, isPlantation = other2016$isPlantation, residual = abs(residuals(otherHeightFromDiameterMichaelisMenten))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))

psmeHeightSharmaPartonResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeHeightFromDiameterSharmaParton))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
alruHeightSharmaPartonResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightFromDiameterSharmaParton))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
tsheHeightSharmaPartonResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterSharmaParton))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
acmaHeightSharmaPartonResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterSharmaParton))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
umcaHeightSharmaPartonResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterSharmaParton))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
thplHeightSharmaPartonResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterSharmaParton))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
otherHeightSharmaPartonResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = other2016$DBH, isPlantation = other2016$isPlantation, residual = abs(residuals(otherHeightFromDiameterSharmaParton))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))

psmeHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1, b1p = 0)) # singular gradient with a1p
alruHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
tsheHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
acmaHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
umcaHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
thplHeightSharmaZhangResidualPower = nlrob(residual ~ a1*DBH^b1, tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterSharmaZhang))), start = list(a1 = 1, b1 = 1)) # a1p, b1p not significant
otherHeightSharmaZhangResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = other2016$DBH, isPlantation = other2016$isPlantation, residual = abs(residuals(otherHeightFromDiameterSharmaZhang))), maxit = 50, start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))

psmeHeightChapmanRichardsResidualPower$confint = confint_nlrob(psmeHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, psmeHeightChapmanRichardsResidualPower$nobs))
alruHeightChapmanRichardsResidualPower$confint = confint_nlrob(alruHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, alruHeightChapmanRichardsResidualPower$nobs))
tsheHeightChapmanRichardsResidualPower$confint = confint_nlrob(tsheHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, tsheHeightChapmanRichardsResidualPower$nobs))
acmaHeightChapmanRichardsResidualPower$confint = confint_nlrob(acmaHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, acmaHeightChapmanRichardsResidualPower$nobs))
umcaHeightChapmanRichardsResidualPower$confint = confint_nlrob(umcaHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, umcaHeightChapmanRichardsResidualPower$nobs))
thplHeightChapmanRichardsResidualPower$confint = confint_nlrob(thplHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, thplHeightChapmanRichardsResidualPower$nobs))
otherHeightChapmanRichardsResidualPower$confint = confint_nlrob(otherHeightChapmanRichardsResidualPower, level = 0.99, weights = rep(1, otherHeightChapmanRichardsResidualPower$nobs))

psmeHeightMichaelisMentenResidualPower$confint = confint_nlrob(psmeHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, psmeHeightMichaelisMentenResidualPower$nobs))
alruHeightMichaelisMentenResidualPower$confint = confint_nlrob(alruHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, alruHeightMichaelisMentenResidualPower$nobs))
tsheHeightMichaelisMentenResidualPower$confint = confint_nlrob(tsheHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, tsheHeightMichaelisMentenResidualPower$nobs))
acmaHeightMichaelisMentenResidualPower$confint = confint_nlrob(acmaHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, acmaHeightMichaelisMentenResidualPower$nobs))
umcaHeightMichaelisMentenResidualPower$confint = confint_nlrob(umcaHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, umcaHeightMichaelisMentenResidualPower$nobs))
thplHeightMichaelisMentenResidualPower$confint = confint_nlrob(thplHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, thplHeightMichaelisMentenResidualPower$nobs))
otherHeightMichaelisMentenResidualPower$confint = confint_nlrob(otherHeightMichaelisMentenResidualPower, level = 0.99, weights = rep(1, otherHeightMichaelisMentenResidualPower$nobs))

psmeHeightSharmaPartonResidualPower$confint = confint_nlrob(psmeHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, psmeHeightSharmaPartonResidualPower$nobs))
alruHeightSharmaPartonResidualPower$confint = confint_nlrob(alruHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, alruHeightSharmaPartonResidualPower$nobs))
tsheHeightSharmaPartonResidualPower$confint = confint_nlrob(tsheHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, tsheHeightSharmaPartonResidualPower$nobs))
acmaHeightSharmaPartonResidualPower$confint = confint_nlrob(acmaHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, acmaHeightSharmaPartonResidualPower$nobs))
umcaHeightSharmaPartonResidualPower$confint = confint_nlrob(umcaHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, umcaHeightSharmaPartonResidualPower$nobs))
thplHeightSharmaPartonResidualPower$confint = confint_nlrob(thplHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, thplHeightSharmaPartonResidualPower$nobs))
otherHeightSharmaPartonResidualPower$confint = confint_nlrob(otherHeightSharmaPartonResidualPower, level = 0.99, weights = rep(1, otherHeightSharmaPartonResidualPower$nobs))

psmeHeightSharmaZhangResidualPower$confint = confint_nlrob(psmeHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, psmeHeightSharmaZhangResidualPower$nobs))
alruHeightSharmaZhangResidualPower$confint = confint_nlrob(alruHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, alruHeightSharmaZhangResidualPower$nobs))
tsheHeightSharmaZhangResidualPower$confint = confint_nlrob(tsheHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, tsheHeightSharmaZhangResidualPower$nobs))
acmaHeightSharmaZhangResidualPower$confint = confint_nlrob(acmaHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, acmaHeightSharmaZhangResidualPower$nobs))
umcaHeightSharmaZhangResidualPower$confint = confint_nlrob(umcaHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, umcaHeightSharmaZhangResidualPower$nobs))
thplHeightSharmaZhangResidualPower$confint = confint_nlrob(thplHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, thplHeightSharmaZhangResidualPower$nobs))
otherHeightSharmaZhangResidualPower$confint = confint_nlrob(otherHeightSharmaZhangResidualPower, level = 0.99, weights = rep(1, otherHeightSharmaZhangResidualPower$nobs))

psmeHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = psme2016$DBH, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 2.1, a1p = -1.6, b1 = 0.17, b1p = 0.34)) # singular gradient
#alruHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 0.94, b1 = 0.38, b1p = -0.04), control = nls.control(maxiter = 500)) # step factor
tsheHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + a1*DBH^b1, tibble(DBH = tshe2016$DBH, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 1, b1 = 1))
acmaHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + a1*DBH^b1, tibble(DBH = acma2016$DBH, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 1, b1 = 1))
umcaHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + a1*DBH^b1, tibble(DBH = umca2016$DBH, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 1, b1 = 1))
thplHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + a1*DBH^b1, tibble(DBH = thpl2016$DBH, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 1, b1 = 1))
otherHeightMichaelisMentenResidualConstPower = nlrob(residual ~ a0 + (a1 + a1p * isPlantation) * DBH^(b1 + b1p * isPlantation), tibble(DBH = other2016$DBH, isPlantation = other2016$isPlantation, residual = abs(residuals(otherHeightFromDiameterMichaelisMenten))), start = list(a0 = 0, a1 = 1, a1p = 0, b1 = 1, b1p = 0))

tribble(~species, ~name, ~mae, ~aic,
        "PSME", "power", mean(abs(residuals(psmeHeightMichaelisMentenResidualPower))), AIC(psmeHeightMichaelisMentenResidualPower),
        "PSME", "const+power", mean(abs(residuals(psmeHeightMichaelisMentenResidualConstPower))), AIC(psmeHeightMichaelisMentenResidualConstPower),
        "TSHE", "power", mean(abs(residuals(tsheHeightMichaelisMentenResidualPower))), AIC(tsheHeightMichaelisMentenResidualPower),
        "TSHE", "const+power", mean(abs(residuals(tsheHeightMichaelisMentenResidualConstPower))), AIC(tsheHeightMichaelisMentenResidualConstPower),
        "ACMA", "power", mean(abs(residuals(acmaHeightMichaelisMentenResidualPower))), AIC(acmaHeightMichaelisMentenResidualPower),
        "ACMA", "const+power", mean(abs(residuals(acmaHeightMichaelisMentenResidualConstPower))), AIC(acmaHeightMichaelisMentenResidualConstPower),
        "UMCA", "power", mean(abs(residuals(umcaHeightMichaelisMentenResidualPower))), AIC(umcaHeightMichaelisMentenResidualPower),
        "UMCA", "const+power", mean(abs(residuals(umcaHeightMichaelisMentenResidualConstPower))), AIC(umcaHeightMichaelisMentenResidualConstPower),
        "THPL", "power", mean(abs(residuals(thplHeightMichaelisMentenResidualPower))), AIC(thplHeightMichaelisMentenResidualPower),
        "THPL", "const+power", mean(abs(residuals(thplHeightMichaelisMentenResidualConstPower))), AIC(thplHeightMichaelisMentenResidualConstPower),
        "other", "power", mean(abs(residuals(otherHeightMichaelisMentenResidualPower))), AIC(otherHeightMichaelisMentenResidualPower),
        "other", "const+power", mean(abs(residuals(otherHeightMichaelisMentenResidualConstPower))), AIC(otherHeightMichaelisMentenResidualConstPower))


## direct fitting of diameter residuals
psmeDiameterGslNlsResidualConstPower = gsl_nls(residual ~ a0 + a1*(TotalHt - 1.37)^b1, tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))), start = list(a0 = 0, a1 = 1, b1 = 1), control = gsl_nls_control(maxiter = 250))
psmeDiameterGslNlsResidualPower = gsl_nls(residual ~ a1*(TotalHt - 1.37)^b1, tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1))
psmeDiameterResidualLm01 = lm(residual ~ 0 + I(TotalHt - 1.37), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))))
psmeDiameterResidualLm0.5 = lm(residual ~ I(sqrt(TotalHt - 1.37)), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))))
psmeDiameterResidualLm1 = lm(residual ~ I(TotalHt - 1.37), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))))
psmeDiameterResidualLm2 = lm(residual ~ I(TotalHt - 1.37) + I((TotalHt - 1.37)^2), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))))
psmeDiameterResidualLm3 = lm(residual ~ I(TotalHt - 1.37) + I((TotalHt - 1.37)^2) + I((TotalHt - 1.37)^3), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))))
psmeDiameterRuarkResidualConstPower = gsl_nls(residual ~ a0 + a1*(TotalHt - 1.37)^b1, tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))), start = list(a0 = 0, a1 = 1, b1 = 1), control = gsl_nls_control(maxiter = 250)) # parameter evaporation
psmeDiameterRuarkResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1, b1p = 0)) # NaN-inf with a1p
psmeDiameterRuarkResidualPower$confint = confint_nlrob(psmeDiameterRuarkResidualPower, level = 0.99, weights = rep(1, psmeDiameterRuarkResidualPower$nobs))

tribble(~name, ~mae, ~aic,
        "gsl_nls const", mean(abs(residuals(psmeDiameterGslNlsResidualConstPower))), AIC(psmeDiameterGslNlsResidualConstPower),
        "gsl_nls", mean(abs(residuals(psmeDiameterGslNlsResidualPower))), AIC(psmeDiameterGslNlsResidualPower),
        "nlrob const", mean(abs(residuals(psmeDiameterRuarkResidualConstPower))), AIC(psmeDiameterRuarkResidualConstPower),
        "nlrob", mean(abs(residuals(psmeDiameterRuarkResidualPower))), AIC(psmeDiameterRuarkResidualPower),
        "lm01", mean(abs(residuals(psmeDiameterResidualLm01))), AIC(psmeDiameterResidualLm01),
        "lm0.5", mean(abs(residuals(psmeDiameterResidualLm0.5))), AIC(psmeDiameterResidualLm0.5),
        "lm1", mean(abs(residuals(psmeDiameterResidualLm1))), AIC(psmeDiameterResidualLm1),
        "lm2", mean(abs(residuals(psmeDiameterResidualLm2))), AIC(psmeDiameterResidualLm2),
        "lm3", mean(abs(residuals(psmeDiameterResidualLm3))), AIC(psmeDiameterResidualLm3)) %>%
  mutate(deltaAic = aic - min(aic)) %>%
  arrange(desc(deltaAic))

tibble(gsl_nls = mean(abs(residuals(psmeDiameterGslNlsResidualPower))), nlrob = mean(abs(residuals(psmeDiameterRuarkResidualPower))), lm0.5 = mean(abs(psmeDiameterResidualLm0.5$residuals)), lm01 = mean(abs(psmeDiameterResidualLm01$residuals)), lm1 = mean(abs(psmeDiameterResidualLm1$residuals)), lm2 = mean(abs(psmeDiameterResidualLm2$residuals)))

ggplot() +
  geom_segment(aes(x = 1.37, xend = 1.37, y = 0, yend = 130), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = psme2016$TotalHt, y = abs(residuals(psmeDiameterFromHeightRuark))), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = psme2016$TotalHt, y = abs(residuals(psmeDiameterFromHeightRuark)), color = "Ruark GAM", fill = "Ruark GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  #geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterGslNlsResidualConstPower), color = "gsl_nls() const", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterGslNlsResidualPower), color = "gsl_nls()", group = psme2016$isPlantation)) +
  #geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterResidualLm01), color = "lm01", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterResidualLm3), color = "lm3", group = psme2016$isPlantation)) +
  #geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterRuarkResidualConstPower), color = "nlrob() const", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$TotalHt, y = predict(psmeDiameterRuarkResidualPower), color = "nlrob()", group = psme2016$isPlantation)) +
  labs(x = "height, m", y = "diameter residual, cm", color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))

alruDiameterGslNlsResidualConstPower = gsl_nls(residual ~ a0 + a1*(TotalHt - 1.37)^b1, tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))), start = list(a0 = 0, a1 = 1, b1 = 1), control = gsl_nls_control(maxiter = 250))
alruDiameterGslNlsResidualPower = gsl_nls(residual ~ a1*(TotalHt - 1.37)^b1, tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1))
alruDiameterResidualLm01 = lm(residual ~ 0 + I(TotalHt - 1.37), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))))
alruDiameterResidualLm0.5 = lm(residual ~ I(sqrt(TotalHt - 1.37)), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))))
alruDiameterResidualLm1 = lm(residual ~ I(TotalHt - 1.37), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))))
alruDiameterResidualLm2 = lm(residual ~ I(TotalHt - 1.37) + I((TotalHt - 1.37)^2), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))))
alruDiameterResidualLm3 = lm(residual ~ I(TotalHt - 1.37) + I((TotalHt - 1.37)^2) + I((TotalHt - 1.37)^3), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))))
alruDiameterRuarkResidualConstPower = gsl_nls(residual ~ a0 + a1*(TotalHt - 1.37)^b1, tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))), start = list(a0 = 0, a1 = 1, b1 = 1), control = gsl_nls_control(maxiter = 250)) # parameter evaporation
alruDiameterRuarkResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightRuark))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
alruDiameterRuarkResidualPower$confint = confint_nlrob(alruDiameterRuarkResidualPower, level = 0.99, weights = rep(1, alruDiameterRuarkResidualPower$nobs))

tribble(~name, ~mae, ~aic,
        "gsl_nls const", mean(abs(residuals(alruDiameterGslNlsResidualConstPower))), AIC(alruDiameterGslNlsResidualConstPower),
        "gsl_nls", mean(abs(residuals(alruDiameterGslNlsResidualPower))), AIC(alruDiameterGslNlsResidualPower),
        "nlrob const", mean(abs(residuals(alruDiameterRuarkResidualConstPower))), AIC(alruDiameterRuarkResidualConstPower),
        "nlrob", mean(abs(residuals(alruDiameterRuarkResidualPower))), AIC(alruDiameterRuarkResidualPower),
        "lm01", mean(abs(residuals(alruDiameterResidualLm01))), AIC(alruDiameterResidualLm01),
        "lm0.5", mean(abs(residuals(alruDiameterResidualLm0.5))), AIC(alruDiameterResidualLm0.5),
        "lm1", mean(abs(residuals(alruDiameterResidualLm1))), AIC(alruDiameterResidualLm1),
        "lm2", mean(abs(residuals(alruDiameterResidualLm2))), AIC(alruDiameterResidualLm2),
        "lm3", mean(abs(residuals(alruDiameterResidualLm3))), AIC(alruDiameterResidualLm3)) %>%
  mutate(deltaAic = aic - min(aic)) %>%
  arrange(desc(deltaAic))

tibble(gsl_nls = mean(abs(residuals(alruDiameterGslNlsResidualPower))), nlrob = mean(abs(residuals(alruDiameterRuarkResidualPower))), lm0.5 = mean(abs(alruDiameterResidualLm0.5$residuals)), lm01 = mean(abs(alruDiameterResidualLm01$residuals)), lm1 = mean(abs(alruDiameterResidualLm1$residuals)), lm2 = mean(abs(alruDiameterResidualLm2$residuals)))

ggplot() +
  geom_segment(aes(x = 1.37, xend = 1.37, y = 0, yend = 70), color = "grey70", linetype = "longdash") +
  geom_point(aes(x = alru2016$TotalHt, y = abs(residuals(alruDiameterFromHeightRuark)), color = isPlantation), alpha = 0.1, color = "grey25", shape = 16) +
  geom_smooth(aes(x = alru2016$TotalHt, y = abs(residuals(alruDiameterFromHeightRuark)), color = "Ruark GAM", fill = "Ruark GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  #geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterGslNlsResidualConstPower), color = "gsl_nls() const", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterGslNlsResidualPower), color = "gsl_nls()", group = alru2016$isPlantation)) +
  #geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterResidualLm01), color = "lm01", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterResidualLm3), color = "lm3", group = alru2016$isPlantation)) +
  #geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterRuarkResidualConstPower), color = "nlrob() const", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$TotalHt, y = predict(alruDiameterRuarkResidualPower), color = "nlrob()", group = alru2016$isPlantation)) +
  labs(x = "height, m", y = "diameter residual, cm", color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))

psmeDiameterChapmanFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightChapmanForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
alruDiameterChapmanFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightChapmanForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
tsheDiameterChapmanFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = tshe2016$TotalHt, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheDiameterFromHeightChapmanForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
acmaDiameterChapmanFormResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = acma2016$TotalHt, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaDiameterFromHeightChapmanForm))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
umcaDiameterChapmanFormResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = umca2016$TotalHt, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaDiameterFromHeightChapmanForm))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
thplDiameterChapmanFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = thpl2016$TotalHt, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplDiameterFromHeightChapmanForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
otherDiameterChapmanFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = other2016$TotalHt, isPlantation = other2016$isPlantation, residual = abs(pmin(-residuals(otherDiameterFromHeightChapmanForm), 125))), start = list(a1 = 1, b1 = 1, b1p = 0), maxit = 150, control = nls.control(maxiter = 500)) # a1p not significant, constrain maximum error to avoid polynomial runaway

psmeDiameterChapmanRichardsResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightChapmanRichards))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
alruDiameterChapmanRichardsResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightChapmanRichards))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
tsheDiameterChapmanRichardsResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = tshe2016$TotalHt, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheDiameterFromHeightChapmanRichards))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
acmaDiameterChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = acma2016$TotalHt, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaDiameterFromHeightChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
umcaDiameterChapmanRichardsResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = umca2016$TotalHt, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaDiameterFromHeightChapmanRichards))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
thplDiameterChapmanRichardsResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = thpl2016$TotalHt, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplDiameterFromHeightChapmanRichards))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
otherDiameterChapmanRichardsResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = other2016$TotalHt, isPlantation = other2016$isPlantation, residual = abs(residuals(otherDiameterFromHeightChapmanRichards))), maxit = 100, start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant

tsheDiameterRuarkResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = tshe2016$TotalHt, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
acmaDiameterRuarkResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = acma2016$TotalHt, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaDiameterFromHeightRuark))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
umcaDiameterRuarkResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = umca2016$TotalHt, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaDiameterFromHeightRuark))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
thplDiameterRuarkResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = thpl2016$TotalHt, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
otherDiameterRuarkResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = other2016$TotalHt, isPlantation = other2016$isPlantation, residual = abs(residuals(otherDiameterFromHeightRuark))), start = list(a1 = 1, b1 = 1, b1p = 0), maxit = 100, control = nls.control(maxiter = 100)) # a1p not significant

psmeDiameterSibbesenFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = psme2016$TotalHt, isPlantation = psme2016$isPlantation, residual = abs(residuals(psmeDiameterFromHeightSibbesenForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
alruDiameterSibbesenFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = alru2016$TotalHt, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruDiameterFromHeightSibbesenForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
tsheDiameterSibbesenFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = tshe2016$TotalHt, isPlantation = tshe2016$isPlantation, residual = abs(residuals(tsheDiameterFromHeightSibbesenForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant
acmaDiameterSibbesenFormResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = acma2016$TotalHt, isPlantation = acma2016$isPlantation, residual = abs(residuals(acmaDiameterFromHeightSibbesenForm))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
umcaDiameterSibbesenFormResidualPower = nlrob(residual ~ (a1 + a1p * isPlantation) * (TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = umca2016$TotalHt, isPlantation = umca2016$isPlantation, residual = abs(residuals(umcaDiameterFromHeightSibbesenForm))), start = list(a1 = 1, a1p = 0, b1 = 1, b1p = 0))
thplDiameterSibbesenFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = thpl2016$TotalHt, isPlantation = thpl2016$isPlantation, residual = abs(residuals(thplDiameterFromHeightSibbesenForm))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
otherDiameterSibbesenFormResidualPower = nlrob(residual ~ a1*(TotalHt - 1.37)^(b1 + b1p * isPlantation), tibble(TotalHt = other2016$TotalHt, isPlantation = other2016$isPlantation, residual = abs(residuals(otherDiameterFromHeightSibbesenForm))), maxit = 100, start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p not significant

psmeDiameterChapmanFormResidualPower$confint = confint_nlrob(psmeDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, psmeDiameterChapmanFormResidualPower$nobs))
alruDiameterChapmanFormResidualPower$confint = confint_nlrob(alruDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, alruDiameterChapmanFormResidualPower$nobs))
tsheDiameterChapmanFormResidualPower$confint = confint_nlrob(tsheDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, tsheDiameterChapmanFormResidualPower$nobs))
acmaDiameterChapmanFormResidualPower$confint = confint_nlrob(acmaDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, acmaDiameterChapmanFormResidualPower$nobs))
umcaDiameterChapmanFormResidualPower$confint = confint_nlrob(umcaDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, umcaDiameterChapmanFormResidualPower$nobs))
thplDiameterChapmanFormResidualPower$confint = confint_nlrob(thplDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, thplDiameterChapmanFormResidualPower$nobs))
otherDiameterChapmanFormResidualPower$confint = confint_nlrob(otherDiameterChapmanFormResidualPower, level = 0.99, weights = rep(1, otherDiameterChapmanFormResidualPower$nobs))

psmeDiameterChapmanRichardsResidualPower$confint = confint_nlrob(psmeDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, psmeDiameterChapmanRichardsResidualPower$nobs))
alruDiameterChapmanRichardsResidualPower$confint = confint_nlrob(alruDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, alruDiameterChapmanRichardsResidualPower$nobs))
tsheDiameterChapmanRichardsResidualPower$confint = confint_nlrob(tsheDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, tsheDiameterChapmanRichardsResidualPower$nobs))
acmaDiameterChapmanRichardsResidualPower$confint = confint_nlrob(acmaDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, acmaDiameterChapmanRichardsResidualPower$nobs))
umcaDiameterChapmanRichardsResidualPower$confint = confint_nlrob(umcaDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, umcaDiameterChapmanRichardsResidualPower$nobs))
thplDiameterChapmanRichardsResidualPower$confint = confint_nlrob(thplDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, thplDiameterChapmanRichardsResidualPower$nobs))
otherDiameterChapmanRichardsResidualPower$confint = confint_nlrob(otherDiameterChapmanRichardsResidualPower, level = 0.99, weights = rep(1, otherDiameterChapmanRichardsResidualPower$nobs))

tsheDiameterRuarkResidualPower$confint = confint_nlrob(tsheDiameterRuarkResidualPower, level = 0.99, weights = rep(1, tsheDiameterRuarkResidualPower$nobs))
acmaDiameterRuarkResidualPower$confint = confint_nlrob(acmaDiameterRuarkResidualPower, level = 0.99, weights = rep(1, acmaDiameterRuarkResidualPower$nobs))
umcaDiameterRuarkResidualPower$confint = confint_nlrob(umcaDiameterRuarkResidualPower, level = 0.99, weights = rep(1, umcaDiameterRuarkResidualPower$nobs))
thplDiameterRuarkResidualPower$confint = confint_nlrob(thplDiameterRuarkResidualPower, level = 0.99, weights = rep(1, thplDiameterRuarkResidualPower$nobs))
#otherDiameterRuarkResidualPower$confint = confint_nlrob(otherDiameterRuarkResidualPower, level = 0.99, weights = rep(1, otherDiameterRuarkResidualPower$nobs))
otherDiameterRuarkResidualPower$confint = matrix(NA_real_, nrow = 3, ncol = 2)
colnames(otherDiameterRuarkResidualPower$confint) = c("0.5%", "99.5%")
rownames(otherDiameterRuarkResidualPower$confint) = c("a1", "b1", "b1p")

psmeDiameterSibbesenFormResidualPower$confint = confint_nlrob(psmeDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, psmeDiameterSibbesenFormResidualPower$nobs))
alruDiameterSibbesenFormResidualPower$confint = confint_nlrob(alruDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, alruDiameterSibbesenFormResidualPower$nobs))
tsheDiameterSibbesenFormResidualPower$confint = confint_nlrob(tsheDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, tsheDiameterSibbesenFormResidualPower$nobs))
acmaDiameterSibbesenFormResidualPower$confint = confint_nlrob(acmaDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, acmaDiameterSibbesenFormResidualPower$nobs))
umcaDiameterSibbesenFormResidualPower$confint = confint_nlrob(umcaDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, umcaDiameterSibbesenFormResidualPower$nobs))
thplDiameterSibbesenFormResidualPower$confint = confint_nlrob(thplDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, thplDiameterSibbesenFormResidualPower$nobs))
otherDiameterSibbesenFormResidualPower$confint = confint_nlrob(otherDiameterSibbesenFormResidualPower, level = 0.99, weights = rep(1, otherDiameterSibbesenFormResidualPower$nobs))

residualPower = tribble(~species, ~htName, ~ht_b1, ~ht_b1p, ~diaName, ~dia_b1, ~dia_b1p, ~ht_b1_min, ~ht_b1_max, ~ht_b1p_min, ~ht_b1p_max, ~dia_b1_min, ~dia_b1_max, ~dia_b1p_min, ~dia_b1p_max,
                        "PSME", "Chapman-Richards", psmeHeightChapmanRichardsResidualPower$m$getPars()["b1"], psmeHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", psmeDiameterChapmanFormResidualPower$m$getPars()["b1"], psmeDiameterChapmanFormResidualPower$m$getPars()["b1p"], psmeHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], psmeHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], psmeHeightChapmanRichardsResidualPower$confint["b1p", "0.5%"], psmeHeightChapmanRichardsResidualPower$confint["b1p", "99.5%"], psmeDiameterChapmanFormResidualPower$confint["b1", "0.5%"], psmeDiameterChapmanFormResidualPower$confint["b1", "99.5%"], psmeDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], psmeDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "PSME", "Michaelis-Menten", psmeHeightMichaelisMentenResidualPower$m$getPars()["b1"], psmeHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", psmeDiameterChapmanRichardsResidualPower$m$getPars()["b1"], psmeDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], psmeHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], psmeHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], psmeHeightMichaelisMentenResidualPower$confint["b1p", "0.5%"], psmeHeightMichaelisMentenResidualPower$confint["b1p", "99.5%"], psmeDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], psmeDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], psmeDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], psmeDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "PSME", "Sharma-Parton", psmeHeightSharmaPartonResidualPower$m$getPars()["b1"], psmeHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", psmeDiameterRuarkResidualPower$m$getPars()["b1"], psmeDiameterRuarkResidualPower$m$getPars()["b1p"], psmeHeightSharmaPartonResidualPower$confint["b1", "0.5%"], psmeHeightSharmaPartonResidualPower$confint["b1", "99.5%"], psmeHeightSharmaPartonResidualPower$confint["b1p", "0.5%"], psmeHeightSharmaPartonResidualPower$confint["b1p", "99.5%"], psmeDiameterRuarkResidualPower$confint["b1", "0.5%"], psmeDiameterRuarkResidualPower$confint["b1", "99.5%"], psmeDiameterRuarkResidualPower$confint["b1p", "0.5%"], psmeDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "PSME", "Sharma-Zhang", psmeHeightSharmaZhangResidualPower$m$getPars()["b1"], psmeHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", psmeDiameterSibbesenFormResidualPower$m$getPars()["b1"], psmeDiameterSibbesenFormResidualPower$m$getPars()["b1p"], psmeHeightSharmaZhangResidualPower$confint["b1", "0.5%"], psmeHeightSharmaZhangResidualPower$confint["b1", "99.5%"], psmeHeightSharmaZhangResidualPower$confint["b1p", "0.5%"], psmeHeightSharmaZhangResidualPower$confint["b1p", "99.5%"], psmeDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], psmeDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], psmeDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], psmeDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "ALRU2", "Chapman-Richards", alruHeightChapmanRichardsResidualPower$m$getPars()["b1"], alruHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", alruDiameterChapmanFormResidualPower$m$getPars()["b1"], alruDiameterChapmanFormResidualPower$m$getPars()["b1p"], alruHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], alruHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], alruHeightChapmanRichardsResidualPower$confint["b1p", "0.5%"], alruHeightChapmanRichardsResidualPower$confint["b1p", "99.5%"], alruDiameterChapmanFormResidualPower$confint["b1", "0.5%"], alruDiameterChapmanFormResidualPower$confint["b1", "99.5%"], alruDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], alruDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "ALRU2", "Michaelis-Menten", alruHeightMichaelisMentenResidualPower$m$getPars()["b1"], alruHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", alruDiameterChapmanRichardsResidualPower$m$getPars()["b1"], alruDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], alruHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], alruHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], alruHeightMichaelisMentenResidualPower$confint["b1p", "0.5%"], alruHeightMichaelisMentenResidualPower$confint["b1p", "99.5%"], alruDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], alruDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], alruDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], alruDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "ALRU2", "Sharma-Parton", alruHeightSharmaPartonResidualPower$m$getPars()["b1"], alruHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", alruDiameterRuarkResidualPower$m$getPars()["b1"], alruDiameterRuarkResidualPower$m$getPars()["b1p"], alruHeightSharmaPartonResidualPower$confint["b1", "0.5%"], alruHeightSharmaPartonResidualPower$confint["b1", "99.5%"], alruHeightSharmaPartonResidualPower$confint["b1p", "0.5%"], alruHeightSharmaPartonResidualPower$confint["b1p", "99.5%"], alruDiameterRuarkResidualPower$confint["b1", "0.5%"], alruDiameterRuarkResidualPower$confint["b1", "99.5%"], alruDiameterRuarkResidualPower$confint["b1p", "0.5%"], alruDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "ALRU2", "Sharma-Zhang", alruHeightSharmaZhangResidualPower$m$getPars()["b1"], alruHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", alruDiameterSibbesenFormResidualPower$m$getPars()["b1"], alruDiameterSibbesenFormResidualPower$m$getPars()["b1p"], alruHeightSharmaZhangResidualPower$confint["b1", "0.5%"], alruHeightSharmaZhangResidualPower$confint["b1", "99.5%"], alruHeightSharmaZhangResidualPower$confint["b1p", "0.5%"], alruHeightSharmaZhangResidualPower$confint["b1p", "99.5%"], alruDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], alruDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], alruDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], alruDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "TSHE", "Chapman-Richards", tsheHeightChapmanRichardsResidualPower$m$getPars()["b1"], tsheHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", tsheDiameterChapmanFormResidualPower$m$getPars()["b1"], tsheDiameterChapmanFormResidualPower$m$getPars()["b1p"], tsheHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], tsheHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, tsheDiameterChapmanFormResidualPower$confint["b1", "0.5%"], tsheDiameterChapmanFormResidualPower$confint["b1", "99.5%"], tsheDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], tsheDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "TSHE", "Michaelis-Menten", tsheHeightMichaelisMentenResidualPower$m$getPars()["b1"], tsheHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", tsheDiameterChapmanRichardsResidualPower$m$getPars()["b1"], tsheDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], tsheHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], tsheHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, tsheDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], tsheDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], tsheDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], tsheDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "TSHE", "Sharma-Parton", tsheHeightSharmaPartonResidualPower$m$getPars()["b1"], tsheHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", tsheDiameterRuarkResidualPower$m$getPars()["b1"], tsheDiameterRuarkResidualPower$m$getPars()["b1p"], tsheHeightSharmaPartonResidualPower$confint["b1", "0.5%"], tsheHeightSharmaPartonResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, tsheDiameterRuarkResidualPower$confint["b1", "0.5%"], tsheDiameterRuarkResidualPower$confint["b1", "99.5%"], tsheDiameterRuarkResidualPower$confint["b1p", "0.5%"], tsheDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "TSHE", "Sharma-Zhang", tsheHeightSharmaZhangResidualPower$m$getPars()["b1"], tsheHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", tsheDiameterSibbesenFormResidualPower$m$getPars()["b1"], tsheDiameterSibbesenFormResidualPower$m$getPars()["b1p"], tsheHeightSharmaZhangResidualPower$confint["b1", "0.5%"], tsheHeightSharmaZhangResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, tsheDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], tsheDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], tsheDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], tsheDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "ACMA3", "Chapman-Richards", acmaHeightChapmanRichardsResidualPower$m$getPars()["b1"], acmaHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", acmaDiameterChapmanFormResidualPower$m$getPars()["b1"], acmaDiameterChapmanFormResidualPower$m$getPars()["b1p"], acmaHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], acmaHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, acmaDiameterChapmanFormResidualPower$confint["b1", "0.5%"], acmaDiameterChapmanFormResidualPower$confint["b1", "99.5%"], acmaDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], acmaDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "ACMA3", "Michaelis-Menten", acmaHeightMichaelisMentenResidualPower$m$getPars()["b1"], acmaHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", acmaDiameterChapmanRichardsResidualPower$m$getPars()["b1"], acmaDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], acmaHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], acmaHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, acmaDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], acmaDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], acmaDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], acmaDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "ACMA3", "Sharma-Parton", acmaHeightSharmaPartonResidualPower$m$getPars()["b1"], acmaHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", acmaDiameterRuarkResidualPower$m$getPars()["b1"], acmaDiameterRuarkResidualPower$m$getPars()["b1p"], acmaHeightSharmaPartonResidualPower$confint["b1", "0.5%"], acmaHeightSharmaPartonResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, acmaDiameterRuarkResidualPower$confint["b1", "0.5%"], acmaDiameterRuarkResidualPower$confint["b1", "99.5%"], acmaDiameterRuarkResidualPower$confint["b1p", "0.5%"], acmaDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "ACMA3", "Sharma-Zhang", acmaHeightSharmaZhangResidualPower$m$getPars()["b1"], acmaHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", acmaDiameterSibbesenFormResidualPower$m$getPars()["b1"], acmaDiameterSibbesenFormResidualPower$m$getPars()["b1p"], acmaHeightSharmaZhangResidualPower$confint["b1", "0.5%"], acmaHeightSharmaZhangResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, acmaDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], acmaDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], acmaDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], acmaDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "UMCA", "Chapman-Richards", umcaHeightChapmanRichardsResidualPower$m$getPars()["b1"], umcaHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", umcaDiameterChapmanFormResidualPower$m$getPars()["b1"], umcaDiameterChapmanFormResidualPower$m$getPars()["b1p"], umcaHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], umcaHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, umcaDiameterChapmanFormResidualPower$confint["b1", "0.5%"], umcaDiameterChapmanFormResidualPower$confint["b1", "99.5%"], umcaDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], umcaDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "UMCA", "Michaelis-Menten", umcaHeightMichaelisMentenResidualPower$m$getPars()["b1"], umcaHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", umcaDiameterChapmanRichardsResidualPower$m$getPars()["b1"], umcaDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], umcaHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], umcaHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, umcaDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], umcaDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], umcaDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], umcaDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "UMCA", "Sharma-Parton", umcaHeightSharmaPartonResidualPower$m$getPars()["b1"], umcaHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", umcaDiameterRuarkResidualPower$m$getPars()["b1"], umcaDiameterRuarkResidualPower$m$getPars()["b1p"], umcaHeightSharmaPartonResidualPower$confint["b1", "0.5%"], umcaHeightSharmaPartonResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, umcaDiameterRuarkResidualPower$confint["b1", "0.5%"], umcaDiameterRuarkResidualPower$confint["b1", "99.5%"], umcaDiameterRuarkResidualPower$confint["b1p", "0.5%"], umcaDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "UMCA", "Sharma-Zhang", umcaHeightSharmaZhangResidualPower$m$getPars()["b1"], umcaHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", umcaDiameterSibbesenFormResidualPower$m$getPars()["b1"], umcaDiameterSibbesenFormResidualPower$m$getPars()["b1p"], umcaHeightSharmaZhangResidualPower$confint["b1", "0.5%"], umcaHeightSharmaZhangResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, umcaDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], umcaDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], umcaDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], umcaDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "THPL", "Chapman-Richards", thplHeightChapmanRichardsResidualPower$m$getPars()["b1"], thplHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", thplDiameterChapmanFormResidualPower$m$getPars()["b1"], thplDiameterChapmanFormResidualPower$m$getPars()["b1p"], thplHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], thplHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, thplDiameterChapmanFormResidualPower$confint["b1", "0.5%"], thplDiameterChapmanFormResidualPower$confint["b1", "99.5%"], thplDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], thplDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "THPL", "Michaelis-Menten", thplHeightMichaelisMentenResidualPower$m$getPars()["b1"], thplHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", thplDiameterChapmanRichardsResidualPower$m$getPars()["b1"], thplDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], thplHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], thplHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, thplDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], thplDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], thplDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], thplDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "THPL", "Sharma-Parton", thplHeightSharmaPartonResidualPower$m$getPars()["b1"], thplHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", thplDiameterRuarkResidualPower$m$getPars()["b1"], thplDiameterRuarkResidualPower$m$getPars()["b1p"], thplHeightSharmaPartonResidualPower$confint["b1", "0.5%"], thplHeightSharmaPartonResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, thplDiameterRuarkResidualPower$confint["b1", "0.5%"], thplDiameterRuarkResidualPower$confint["b1", "99.5%"], thplDiameterRuarkResidualPower$confint["b1p", "0.5%"], thplDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "THPL", "Sharma-Zhang", thplHeightSharmaZhangResidualPower$m$getPars()["b1"], thplHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", thplDiameterSibbesenFormResidualPower$m$getPars()["b1"], thplDiameterSibbesenFormResidualPower$m$getPars()["b1p"], thplHeightSharmaZhangResidualPower$confint["b1", "0.5%"], thplHeightSharmaZhangResidualPower$confint["b1", "99.5%"], NA_real_, NA_real_, thplDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], thplDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], thplDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], thplDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"],
                        "other", "Chapman-Richards", otherHeightChapmanRichardsResidualPower$m$getPars()["b1"], otherHeightChapmanRichardsResidualPower$m$getPars()["b1p"], "Chapman-Richards form", otherDiameterChapmanFormResidualPower$m$getPars()["b1"], otherDiameterChapmanFormResidualPower$m$getPars()["b1p"], otherHeightChapmanRichardsResidualPower$confint["b1", "0.5%"], otherHeightChapmanRichardsResidualPower$confint["b1", "99.5%"], otherHeightChapmanRichardsResidualPower$confint["b1p", "0.5%"], otherHeightChapmanRichardsResidualPower$confint["b1p", "99.5%"], otherDiameterChapmanFormResidualPower$confint["b1", "0.5%"], otherDiameterChapmanFormResidualPower$confint["b1", "99.5%"], otherDiameterChapmanFormResidualPower$confint["b1p", "0.5%"], otherDiameterChapmanFormResidualPower$confint["b1p", "99.5%"],
                        "other", "Michaelis-Menten", otherHeightMichaelisMentenResidualPower$m$getPars()["b1"], otherHeightMichaelisMentenResidualPower$m$getPars()["b1p"], "Chapman-Richards", otherDiameterChapmanRichardsResidualPower$m$getPars()["b1"], otherDiameterChapmanRichardsResidualPower$m$getPars()["b1p"], otherHeightMichaelisMentenResidualPower$confint["b1", "0.5%"], otherHeightMichaelisMentenResidualPower$confint["b1", "99.5%"], otherHeightMichaelisMentenResidualPower$confint["b1p", "0.5%"], otherHeightMichaelisMentenResidualPower$confint["b1p", "99.5%"], otherDiameterChapmanRichardsResidualPower$confint["b1", "0.5%"], otherDiameterChapmanRichardsResidualPower$confint["b1", "99.5%"], otherDiameterChapmanRichardsResidualPower$confint["b1p", "0.5%"], otherDiameterChapmanRichardsResidualPower$confint["b1p", "99.5%"],
                        "other", "Sharma-Parton", otherHeightSharmaPartonResidualPower$m$getPars()["b1"], otherHeightSharmaPartonResidualPower$m$getPars()["b1p"], "Ruark", otherDiameterRuarkResidualPower$m$getPars()["b1"], otherDiameterRuarkResidualPower$m$getPars()["b1p"], otherHeightSharmaPartonResidualPower$confint["b1", "0.5%"], otherHeightSharmaPartonResidualPower$confint["b1", "99.5%"], otherHeightSharmaPartonResidualPower$confint["b1p", "0.5%"], otherHeightSharmaPartonResidualPower$confint["b1p", "99.5%"], otherDiameterRuarkResidualPower$confint["b1", "0.5%"], otherDiameterRuarkResidualPower$confint["b1", "99.5%"], otherDiameterRuarkResidualPower$confint["b1p", "0.5%"], otherDiameterRuarkResidualPower$confint["b1p", "99.5%"],
                        "other", "Sharma-Zhang", otherHeightSharmaZhangResidualPower$m$getPars()["b1"], otherHeightSharmaZhangResidualPower$m$getPars()["b1p"], "Sibbesen form", otherDiameterSibbesenFormResidualPower$m$getPars()["b1"], otherDiameterSibbesenFormResidualPower$m$getPars()["b1p"], otherHeightSharmaZhangResidualPower$confint["b1", "0.5%"], otherHeightSharmaZhangResidualPower$confint["b1", "99.5%"], otherHeightSharmaZhangResidualPower$confint["b1p", "0.5%"], otherHeightSharmaZhangResidualPower$confint["b1p", "99.5%"], otherDiameterSibbesenFormResidualPower$confint["b1", "0.5%"], otherDiameterSibbesenFormResidualPower$confint["b1", "99.5%"], otherDiameterSibbesenFormResidualPower$confint["b1p", "0.5%"], otherDiameterSibbesenFormResidualPower$confint["b1p", "99.5%"]) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"),  levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         ht_b1p_min = ht_b1_min + ht_b1p_min + sqrt(2 * (ht_b1 - ht_b1_min) * (ht_b1p - ht_b1p_min)), # shift confidence intervals from plantation effects to combined effects: the covariance is equal and opposite the main effect variance
         ht_b1p_max = ht_b1_max + ht_b1p_max - sqrt(2 * (ht_b1 - ht_b1_max) * (ht_b1p - ht_b1p_max)),
         dia_b1p_min = dia_b1_min + dia_b1p_min + sqrt(2 * (dia_b1 - dia_b1_min) * (dia_b1p - dia_b1p_min)), 
         dia_b1p_max = dia_b1_max + dia_b1p_max - sqrt(2 * (dia_b1 - dia_b1_max) * (dia_b1p - dia_b1p_max)))
#residualPower %>% mutate(sum_center = ht_b1 + ht_b1p, sum_min = ht_b1_min + ht_b1p_min, ht_b1p_min2 = ht_b1_min + ht_b1p_min + sqrt(2*(ht_b1 - ht_b1_min) * (ht_b1p - ht_b1p_min))) %>% select(species, ht_b1_min, ht_b1, ht_b1p_min, ht_b1p, sum_min, ht_b1p_min2, sum_center)
#residualPower %>% filter(species == "other species") %>% select(-starts_with("ht_"))
residualPower %>% select(species, ht_b1, ht_b1p, dia_b1, dia_b1p) %>%
  mutate(ht_b1p = ht_b1 + ht_b1p, dia_b1p = dia_b1 + dia_b1p) %>% 
  group_by(species) %>%
  summarize(ht_b1 = round(2 * mean(ht_b1), 1), ht_b1p = round(2 * mean(ht_b1p), 1), # Table S3
            dia_b1 = round(2 * mean(dia_b1), 1), dia_b1p = round(2 * mean(dia_b1p), 1))
residualPower %>% summarize(ht_b1 = mean(ht_b1), ht_b1p = mean(ht_b1 + if_else(is.na(ht_b1p), 0, ht_b1p)), dia_b1 = mean(dia_b1), dia_b1p = mean(dia_b1 + dia_b1p))


## Figures S11-13: estimates of residual variance
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
ggsave("trees/height-diameter/figures/Figure S11 height interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

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
ggsave("trees/height-diameter/figures/Figure S12 diameter interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

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
ggsave("trees/height-diameter/figures/Figure S13 residual power estimates.png", height = 14, width = 20, units = "cm", dpi = 150)


## Figures S14-20: Q-Q plots
plot_qq(psmeHeightFromDiameterPreferred$chapmanRichards, psmeHeightFromDiameterPreferred$michaelisMenten, psmeHeightFromDiameterPreferred$sharmaParton, psmeHeightFromDiameterPreferred$sharmaZhang,
        psmeDiameterFromHeightPreferred$ChapmanRichards, psmeDiameterFromHeightPreferred$chapmanForm, psmeDiameterFromHeightPreferred$ruark, psmeDiameterFromHeightPreferred$sibbesenForm,
        "Douglas-fir")
ggsave("trees/height-diameter/figures/Figure S14 PSME Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(alruHeightFromDiameterPreferred$chapmanRichards, alruHeightFromDiameterPreferred$michaelisMenten, alruHeightFromDiameterPreferred$sharmaParton, alruHeightFromDiameterPreferred$sharmaZhang,
        alruDiameterFromHeightPreferred$chapmanRichards, alruDiameterFromHeightPreferred$chapmanForm, alruDiameterFromHeightPreferred$ruark, alruDiameterFromHeightPreferred$sibbesenForm,
        "red alder", tDegreesOfFreedom = 8, tSkew = 2.1)
ggsave("trees/height-diameter/figures/Figure S15 ALRU2 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(tsheHeightFromDiameterPreferred$chapmanRichards, tsheHeightFromDiameterPreferred$michaelisMenten, tsheHeightFromDiameterPreferred$sharmaParton, tsheHeightFromDiameterPreferred$sharmaZhang,
        tsheDiameterFromHeightPreferred$chapmanRichards, tsheDiameterFromHeightPreferred$chapmanForm, tsheDiameterFromHeightPreferred$ruark, tsheDiameterFromHeightPreferred$sibbesenForm,
        "western hemlock", tDegreesOfFreedom = 7, tSkew = 2)
ggsave("trees/height-diameter/figures/Figure S16 TSHE Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(acmaHeightFromDiameterPreferred$chapmanRichards, acmaHeightFromDiameterPreferred$michaelisMenten, acmaHeightFromDiameterPreferred$sharmaParton, acmaHeightFromDiameterPreferred$sharmaZhang,
        acmaDiameterFromHeightPreferred$chapmanRichards, acmaDiameterFromHeightPreferred$chapmanForm, acmaDiameterFromHeightPreferred$ruark, acmaDiameterFromHeightPreferred$sibbesenForm,
        "bigleaf maple", tDegreesOfFreedom = 10, tSkew = 4)
ggsave("trees/height-diameter/figures/Figure S17 ACMA3 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(umcaHeightFromDiameterPreferred$chapmanRichards, umcaHeightFromDiameterPreferred$michaelisMenten, umcaHeightFromDiameterPreferred$sharmaParton, umcaHeightFromDiameterPreferred$sharmaZhang,
        umcaDiameterFromHeightPreferred$chapmanRichards, umcaDiameterFromHeightPreferred$chapmanForm, umcaDiameterFromHeightPreferred$ruark, umcaDiameterFromHeightPreferred$sibbesenForm,
        "Oregon myrtle", tDegreesOfFreedom = 8, tSkew = 4)
ggsave("trees/height-diameter/figures/Figure S18 UMCA Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(thplHeightFromDiameterPreferred$chapmanRichards, thplHeightFromDiameterPreferred$michaelisMenten, thplHeightFromDiameterPreferred$sharmaParton, thplHeightFromDiameterPreferred$sharmaZhang,
        thplDiameterFromHeightPreferred$chapmanRichards, thplDiameterFromHeightPreferred$chapmanForm, thplDiameterFromHeightPreferred$ruark, thplDiameterFromHeightPreferred$sibbesenForm,
        "western redcedar", tDegreesOfFreedom = 8, tSkew = 3)
ggsave("trees/height-diameter/figures/Figure S19 THPL Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(otherHeightFromDiameterPreferred$chapmanRichards, otherHeightFromDiameterPreferred$michaelisMenten, otherHeightFromDiameterPreferred$sharmaParton, otherHeightFromDiameterPreferred$sharmaZhang,
        otherDiameterFromHeightPreferred$chapmanRichards, otherDiameterFromHeightPreferred$chapmanForm, otherDiameterFromHeightPreferred$ruark, otherDiameterFromHeightPreferred$sibbesenForm,
        "other species", tDegreesOfFreedom = 3, tSkew = 1)
ggsave("trees/height-diameter/figures/Figure S20 other Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)


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