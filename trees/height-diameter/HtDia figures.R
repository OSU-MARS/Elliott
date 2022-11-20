# load data from Elliott Stand Data Feb2022.R, get regressions and summaries from HtDia *.R
# form regression summary
heightDiameterParameters = bind_rows(psmeParameters, alruParameters, tsheParameters, acmaParameters,
                                     umcaParameters, thplParameters, otherParameters)
#write_xlsx(heightDiameterParameters, "trees/height-diameter/HtDia parameters.xlsx")

heightDiameterResults = bind_rows(psmeHeightFromDiameterResults,
                                  psmeHeightFromDiameterResultsGnls,
                                  alruHeightFromDiameterResults,
                                  alruHeightFromDiameterResultsGnls,
                                  tsheHeightFromDiameterResults,
                                  tsheHeightFromDiameterResultsGnls,
                                  acmaHeightFromDiameterResults,
                                  acmaHeightFromDiameterResultsGnls,
                                  umcaHeightFromDiameterResults,
                                  umcaHeightFromDiameterResultsGnls,
                                  thplHeightFromDiameterResults,
                                  thplHeightFromDiameterResultsGnls,
                                  otherHeightFromDiameterResults,
                                  otherHeightFromDiameterResultsGnls,
                                  psmeDiameterFromHeightResults,
                                  alruDiameterFromHeightResults,
                                  tsheDiameterFromHeightResults,
                                  acmaDiameterFromHeightResults,
                                  umcaDiameterFromHeightResults,
                                  thplDiameterFromHeightResults,
                                  otherDiameterFromHeightResults) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"), levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         speciesFraction = recode(species, "Douglas-fir" = 0.750, "red alder" = 0.101, "western hemlock" = 0.056, "bigleaf maple" = 0.029, "Oregon myrtle" = 0.025, "western redcedar" = 0.013, "other species" = 0.017),
         weighting = if_else(fitting %in% c("gnls", "nlrob"), "reweighted", "fixed weights"),
         sizeShapeAlpha = as.factor(if_else(significant == TRUE, weighting, "not significant"))) %>%
  group_by(responseVariable, species) %>%
  mutate(deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) %>%
  group_by(responseVariable, name) %>%
  mutate(weightedMae = sum(speciesFraction * if_else(is.na(mae), 100, mae)))
heightDiameterAccuracyLevels = heightDiameterResults %>% summarize(weightedMae = weightedMae[1], .groups = "keep") %>%
  arrange(weightedMae)

#heightDiameterResults %>% filter(name == "linear") %>% select(responseVariable, species, nse, pearson)

print(heightDiameterResults %>% filter(str_detect(name, "GNLS") == FALSE, str_detect(name, "RelHt") == FALSE) %>% 
        group_by(responseVariable, species) %>% 
        slice_min(pae, n = 4) %>% 
        select(-responseVariable, -biasNaturalRegen, -biasPlantation, -rmseNaturalRegen, -rmsePlantation, -nseNaturalRegen, -nsePlantation, -pearson, -pearsonNaturalRegen, -pearsonPlantation, -aic, -bic), n = 60)
#write_xlsx(heightDiameterResults, "trees/height-diameter/HtDia results.xlsx")

# summarize fitting success for Section 2.4
heightDiameterResults %>% ungroup() %>%
  summarize(nlrob = sum(fitting == "nlrob", na.rm = TRUE), gslNls = sum(fitting == "gsl_nls", na.rm = TRUE), gnls = sum((fitting == "gnls") | (is.na(fitting) & str_detect(name, "GNLS")), na.rm = TRUE), lm = sum(fitting == "lm", na.rm = TRUE), fail = sum(is.na(fitting)), uniqueNonlinear = n() - gnls - lm) %>%
  mutate(nlrobPct = 100 * nlrob / uniqueNonlinear, gslNlsPct = nlrobPct + 100 * gslNls / uniqueNonlinear, totalFailPct = 100 * fail / uniqueNonlinear)
heightDiameterResults %>% filter(is.na(fitting)) %>% select(responseVariable, species, name)

## Figure 1: overall dataset summary
plot_exploratory(liveUnbrokenTrees2016 %>% filter(isConifer), speciesLabel = "conifer", maxTreesMeasured = 170, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(isConifer == FALSE), speciesLabel = "broadleaf", maxTreesMeasured = 170, plotLetters = c("d)", "e)", "f)")) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/Figure 1 height-diameter distribution.png", height = 13, width = 22, units = "cm", dpi = 150)
#plot_exploratory(liveUnbrokenTrees2016) + plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/Figure 1 height-diameter distribution.png", height = 10, width = 22, units = "cm", dpi = 150)


## Figure 2: height-diameter error summary
heightFromDiameterResults = heightDiameterResults %>% filter(responseVariable == "height", str_detect(name, "GNLS") == FALSE, str_detect(name, "RelHt") == FALSE) %>%
  group_by(responseVariable, species) %>%
  mutate(name = factor(name, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "height"))$name),
         deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) # recalculate ΔAICn with relative height forms excluded, requires regrouping by species

ggplot(heightFromDiameterResults) +
  geom_point(aes(x = pae, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(heightFromDiameterResults) +
  geom_point(aes(x = rmse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 11)) +
  labs(x = "RMSE, m", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_segment(x = 0.24, xend = 0.26, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.24, xend = 0.26, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_point(aes(x = deltaAicN, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  labs(x = bquote("normalized "*Delta*"AIC"), y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  coord_cartesian(xlim = c(0, 0.25)) + # exclude high AIC of linear, parabolic, and Douglas-fir power+Curtis fits to avoid squashing of nonlinear differences
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_point(aes(x = nse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_point(aes(x = nse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Pearson's R", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 2), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/Figure 2 height accuracy.png", height = 11, width = 22, units = "cm", dpi = 150)


## Figure 3: diameter-height error summary
diameterFromHeightResults = heightDiameterResults %>% filter(responseVariable == "DBH", str_detect(name, "BAL") == FALSE) %>%
  group_by(responseVariable, species) %>%
  mutate(name = factor(name, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "DBH"))$name),
         deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) # recalculate ΔAICn with BAL names excluded

ggplot(diameterFromHeightResults) +
  geom_point(aes(x = pae, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(diameterFromHeightResults) +
  geom_point(aes(x = rmse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "RMSE, cm", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_segment(x = 0.245, xend = 0.265, y = "Schnute", yend = "Schnute", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", size = 0.4) +
  geom_segment(x = 0.245, xend = 0.265, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", size = 0.4) +
  geom_segment(x = 0.245, xend = 0.265, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", size = 0.4) +
  geom_point(aes(x = deltaAicN, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 0.255)) +
  labs(x = bquote("normalized "*Delta*"AIC"), y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  #scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), trans = scales::pseudo_log_trans()) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_point(aes(x = nse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "model efficiency", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_point(aes(x = pearson, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Pearson's R", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 5, guides = "collect") &
  guides(color = guide_legend(byrow = TRUE, order = 1, ncol = 4), alpha = guide_legend(byrow = TRUE, order = 2, ncol = 2), shape = guide_legend(byrow = TRUE, order = 2, ncol = 2), size = guide_legend(byrow = TRUE, order = 2, ncol = 2)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/Figure 3 diameter accuracy.png", height = 11, width = 22, units = "cm", dpi = 150)

## Figure 4: Douglas-fir, red alder, and western hemlock regression comparision
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = psme2016physio$DBH, y = psmeHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007")) + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) Douglas-fir height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  #scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = psmeDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = psme2016physio$TotalHt, color = "Chapman-Richards physiological", group = psme2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = psmeDiameterFromHeightChapmanRichards$fitted.values, y = psme2016$TotalHt, color = "Chapman-Richards", group = psme2016$isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightSibbesenForm$fitted.values, y = psme2016$TotalHt, color = "Sibbesen", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007")) + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) Douglas-fir diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alru2016physio$DBH, y = alruHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) red alder height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  #scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alruDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = alru2016physio$TotalHt, color = "Chapman-Richards physiological", group = alru2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = alruDiameterFromHeightChapmanRichards$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation)) +
  geom_line(aes(x = alruDiameterFromHeightSibbesenForm$fitted.values, y = alru2016$TotalHt, color = "Sibbesen", group = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) red alder diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = tshe2016physio$DBH, y = tsheHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = tshe2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "e) western hemlock height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = tsheDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = tshe2016physio$TotalHt, color = "Chapman-Richards physiological", group = tshe2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = tsheDiameterFromHeightChapmanRichards$fitted.values, y = tshe2016$TotalHt, color = "Chapman-Richards", group = tshe2016$isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightSibbesenForm$fitted.values, y = tshe2016$TotalHt, color = "Sibbesen", group = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "f) western hemlock diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
plot_layout(nrow = 3, ncol = 2)
#ggsave("trees/height-diameter/Figure 4 PSME-ALRU2-TSHE curves.png", height = 20, width = 17, units = "cm")

## Figure 5: bigleaf maple, Oregon myrtle, and western redcedar regression comparison
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = acma2016physio$DBH, y = acmaHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = acma2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = acma2016$isPlantation)) +
  geom_line(aes(x = acma2016$DBH, y = acmaHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = acma2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) bigleaf maple height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = acma2016$DBH, y = acma2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = acmaDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = acma2016physio$TotalHt, color = "Chapman-Richards physiological", group = acma2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = acmaDiameterFromHeightChapmanRichards$fitted.values, y = acma2016$TotalHt, color = "Chapman-Richards", group = acma2016$isPlantation)) +
  geom_line(aes(x = acmaDiameterFromHeightSibbesenForm$fitted.values, y = acma2016$TotalHt, color = "Sibbesen", group = acma2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) bigleaf maple diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = umca2016physio$DBH, y = umcaHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = umca2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = umca2016$isPlantation)) +
  geom_line(aes(x = umca2016$DBH, y = umcaHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = umca2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) Oregon myrtle height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  #scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
ggplot() +
  geom_point(aes(x = umca2016$DBH, y = umca2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = umcaDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = umca2016physio$TotalHt, color = "Chapman-Richards physiological", group = umca2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = umcaDiameterFromHeightChapmanRichards$fitted.values, y = umca2016$TotalHt, color = "Chapman-Richards", group = umca2016$isPlantation)) +
  geom_line(aes(x = umcaDiameterFromHeightSibbesenForm$fitted.values, y = umca2016$TotalHt, color = "Sibbesen", group = umca2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) Oregon myrtle diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = thpl2016physio$DBH, y = thplHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = thpl2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thpl2016$DBH, y = thplHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = thpl2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "e) western redcedar height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = thpl2016$DBH, y = thpl2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = thplDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = thpl2016physio$TotalHt, color = "Chapman-Richards physiological", group = thpl2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = thplDiameterFromHeightChapmanRichards$fitted.values, y = thpl2016$TotalHt, color = "Chapman-Richards", group = thpl2016$isPlantation)) +
  geom_line(aes(x = thplDiameterFromHeightSibbesenForm$fitted.values, y = thpl2016$TotalHt, color = "Sibbesen", group = thpl2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "f) western redcedar diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "none") +
plot_layout(nrow = 3, ncol = 2)
#ggsave("trees/height-diameter/Figure 5 ACMA3-UMCA-THPL curves.png", height = 20, width = 17, units = "cm")

## Figure 6: minority species merged regression comparison
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = other2016physio$DBH, y = otherHeightFromDiameterSharmaPartonBalPhysio$fitted.values, color = "Sharma-Parton BAL physiological", group = other2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = other2016$isPlantation)) +
  geom_line(aes(x = other2016$DBH, y = otherHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = other2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) other species height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  #scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) +
ggplot() +
  geom_point(aes(x = other2016$DBH, y = other2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = otherDiameterFromHeightChapmanRichardsPhysio$fitted.values, y = other2016physio$TotalHt, color = "Chapman-Richards physiological", group = other2016physio$isPlantation), alpha = 0.5) +
  geom_line(aes(x = otherDiameterFromHeightChapmanRichards$fitted.values, y = other2016$TotalHt, color = "Chapman-Richards", group = other2016$isPlantation)) +
  geom_line(aes(x = otherDiameterFromHeightSibbesenForm$fitted.values, y = other2016$TotalHt, color = "Sibbesen", group = other2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) other species diameter", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Chapman-Richards physiological", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 0.9))
#ggsave("trees/height-diameter/Figure 6 other species curves.png", height = 20/3 + 0.5, width = 17, units = "cm")


## Figures S1-3: species level exploratory plots
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "DF"), speciesLabel = "Douglas-fir", maxTreesMeasured = 150, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "RA"), speciesLabel = "red alder", maxTreesMeasured = 150, plotLetters = c("d)", "e)", "f)"), omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "WH"), speciesLabel = "western hemlock", maxTreesMeasured = 150, plotLetters = c("g)", "h)", "i)"), omitLegends = TRUE) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/Figure S01 PSME-ALRU2-TSHE.png", height = 18, width = 22, units = "cm", dpi = 150)

plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "BM"), speciesLabel = "bigleaf maple", maxTreesMeasured = 150, omitLegends = TRUE, omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "OM"), speciesLabel = "Oregon myrtle", maxTreesMeasured = 150, distributionLegendPositionY = 0.92, plotLetters = c("d)", "e)", "f)"), omitXlabels = TRUE) /
plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "RC"), speciesLabel = "western redcedar", maxTreesMeasured = 150, plotLetters = c("g)", "h)", "i)"), omitLegends = TRUE) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
#ggsave("trees/height-diameter/Figure S02 ACMA3-UMCA-THPL.png", height = 18, width = 22, units = "cm", dpi = 150)

plot_exploratory(liveUnbrokenTrees2016 %>% filter(speciesGroup == "other"), speciesLabel = "other species ", distributionLegendPositionY = 0.92) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
ggsave("trees/height-diameter/Figure S03 other species.png", height = 1/3*(18 - 1) + 1, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 1 Douglas-fir.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 2 red alder.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 3 western hemlock.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 4 bigleaf maple.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 5 Oregon myrtle.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 6 western redcedar.png", height = 10, width = 22, units = "cm", dpi = 150)
#ggsave("trees/height-diameter/exploratory 7 minority species.png", height = 10, width = 22, units = "cm", dpi = 150)


## Figures S4-6: estimates of residual variance
# heightIqr, diameterIqr, and residualPower are calculated in residuals.R
ggplot() +
  geom_line(aes(x = seq(0, 250), y = seq(0, 250)^0.5), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = DBH, y = iqr, color = species, group = paste(species, name, isPlantation)), heightIqr %>% filter(name == "Michaelis-Menten", isPlantation == FALSE)) +
  annotate("text", x = 0, y = 24, label = "a) natural regeneration", hjust = 0, size = 3.5) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "DBH, cm", y = "Michaelis-Menten interquartile height range, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_line(aes(x = seq(0, 250), y = seq(0, 250)^0.5), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = DBH, y = iqr, color = species, group = paste(species, name, isPlantation)), heightIqr %>% filter(name == "Michaelis-Menten", isPlantation == TRUE)) +
  annotate("text", x = 0, y = 24, label = "b) plantations", hjust = 0, size = 3.5) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.justification = c(1, 0), legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
ggsave("trees/height-diameter/Figure S04 height interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

ggplot() +
  geom_line(aes(x = seq(0, 100), y = 80/80*seq(0, 100)^0.9), color = "grey70", linetype = "longdash") +
  geom_line(aes(x = TotalHt, y = iqr, color = species, group = paste(species, name, isPlantation)), diameterIqr %>% filter(isPlantation == FALSE)) +
  annotate("text", x = 0, y = 62, label = "a) natural regeneration", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 82), ylim = c(0, 61)) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "height, m", y = "Ruark interquartile DBH range, cm", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_line(aes(x = TotalHt, y = iqr, color = species, group = paste(species, name, isPlantation)), diameterIqr %>% filter(isPlantation == TRUE)) +
  geom_line(aes(x = seq(0, 100), y = 30/80*seq(0, 100)^1.2), color = "grey70", linetype = "longdash") +
  annotate("text", x = 0, y = 62, label = "b) plantations", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 82), ylim = c(0, 61)) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "height, m", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.justification = c(1, 0), legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt")))
ggsave("trees/height-diameter/Figure S05 diameter interquartile range.png", height = 10, width = 18, units = "cm", dpi = 150)

ggplot(residualPower) +
  geom_errorbar(aes(xmin = ht_b1_min, xmax = ht_b1_max, y = htName, color = species), width = 0.25) +
  geom_errorbar(aes(xmin = ht_b1p_min, xmax = ht_b1p_max, y = htName, color = species), width = 0.25, na.rm = TRUE) +
  geom_point(aes(x = ht_b1, y = htName, color = species, shape = "b1"), size = 1.7) +
  geom_point(aes(x = ht_b1 + ht_b1p, y = htName, color = species, shape = "b1p"), size = 1.7, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 2)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  labs(x = "power fitted to height residuals", y = NULL, color = NULL, shape = NULL) +
  scale_y_discrete(limits = rev) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
ggplot(residualPower) +
  geom_errorbar(aes(xmin = dia_b1_min, xmax = dia_b1_max, y = diaName, color = species), width = 0.25) +
  geom_errorbar(aes(xmin = dia_b1p_min, xmax = dia_b1p_max, y = diaName, color = species), width = 0.25) +
  geom_point(aes(x = dia_b1, y = diaName, color = species, shape = "b1"), size = 1.7) +
  geom_point(aes(x = dia_b1 + dia_b1p, y = diaName, color = species, shape = "b1p"), size = 1.7) +
  coord_cartesian(xlim = c(0, 2)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  labs(x = "power fitted to DBH residuals", y = NULL, color = NULL, shape = NULL) +
  scale_y_discrete(limits = rev) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_blank()) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, guides = "collect") &
  guides(color = "none") &
  scale_color_manual(breaks = levels(heightIqr$species), limits = levels(heightIqr$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  scale_shape_manual(breaks = c("b1", "b1p"), labels = c("base power", "plantations"), values = c(15, 17)) &
  theme(legend.position = "bottom")
ggsave("trees/height-diameter/Figure S06 residual power estimates.png", height = 14, width = 20, units = "cm", dpi = 150)


## Figures S7-: Q-Q plots
plot_qq(psmeHeightFromDiameterChapmanRichards, psmeHeightFromDiameterMichaelisMenten, psmeHeightFromDiameterSharmaParton, psmeHeightFromDiameterSharmaZhang,
        psmeDiameterFromHeightChapmanRichards, psmeDiameterFromHeightChapmanForm, psmeDiameterFromHeightRuark, psmeDiameterFromHeightSibbesenForm,
        "Douglas-fir")
ggsave("trees/height-diameter/Figure S07 PSME Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(alruHeightFromDiameterChapmanRichards, alruHeightFromDiameterMichaelisMenten, alruHeightFromDiameterSharmaParton, alruHeightFromDiameterSharmaZhang,
        alruDiameterFromHeightChapmanRichards, alruDiameterFromHeightChapmanForm, alruDiameterFromHeightRuark, alruDiameterFromHeightSibbesenForm,
        "red alder", tDegreesOfFreedom = 8, tSkew = 2.1)
ggsave("trees/height-diameter/Figure S08 ALRU2 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(tsheHeightFromDiameterChapmanRichards, tsheHeightFromDiameterMichaelisMenten, tsheHeightFromDiameterSharmaParton, tsheHeightFromDiameterSharmaZhang,
        tsheDiameterFromHeightChapmanRichards, tsheDiameterFromHeightChapmanForm, tsheDiameterFromHeightRuark, tsheDiameterFromHeightSibbesenForm,
        "western hemlock", tDegreesOfFreedom = 7, tSkew = 2)
ggsave("trees/height-diameter/Figure S08 TSHE Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(acmaHeightFromDiameterChapmanRichards, acmaHeightFromDiameterMichaelisMenten, acmaHeightFromDiameterSharmaParton, acmaHeightFromDiameterSharmaZhang,
        acmaDiameterFromHeightChapmanRichards, acmaDiameterFromHeightChapmanForm, acmaDiameterFromHeightRuark, acmaDiameterFromHeightSibbesenForm,
        "bigleaf maple", tDegreesOfFreedom = 10, tSkew = 4)
ggsave("trees/height-diameter/Figure S09 ACMA3 Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(umcaHeightFromDiameterChapmanRichards, umcaHeightFromDiameterMichaelisMenten, umcaHeightFromDiameterSharmaParton, umcaHeightFromDiameterSharmaZhang,
        umcaDiameterFromHeightChapmanRichards, umcaDiameterFromHeightChapmanForm, umcaDiameterFromHeightRuark, umcaDiameterFromHeightSibbesenForm,
        "Oregon myrtle", tDegreesOfFreedom = 8, tSkew = 4)
ggsave("trees/height-diameter/Figure S10 UMCA Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(thplHeightFromDiameterChapmanRichards, thplHeightFromDiameterMichaelisMenten, thplHeightFromDiameterSharmaParton, thplHeightFromDiameterSharmaZhang,
        thplDiameterFromHeightChapmanRichards, thplDiameterFromHeightChapmanForm, thplDiameterFromHeightRuark, thplDiameterFromHeightSibbesenForm,
        "western redcedar", tDegreesOfFreedom = 8, tSkew = 3)
ggsave("trees/height-diameter/Figure S11 THPL Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)

plot_qq(otherHeightFromDiameterChapmanRichards, otherHeightFromDiameterMichaelisMenten, otherHeightFromDiameterSharmaParton, otherHeightFromDiameterSharmaZhang,
        otherDiameterFromHeightChapmanRichards, otherDiameterFromHeightChapmanForm, otherDiameterFromHeightRuark, otherDiameterFromHeightSibbesenForm,
        "other species", tDegreesOfFreedom = 3, tSkew = 1)
ggsave("trees/height-diameter/Figure S12 other Q-Q.png", height = 11, width = 16, units = "cm", dpi = 150)


## companion variance structure summary from GNLS
ggplot(heightDiameterResults %>% filter(is.na(power) == FALSE) %>% mutate(name = str_replace(name, " GNLS", ""))) +
  geom_point(aes(x = power, y = name, color = species)) +
  coord_cartesian(xlim = c(0, 2)) +
  facet_grid(rows = vars(species), labeller = label_wrap_gen(width = 10), switch = "y") +
  guides(color = guide_legend(ncol = 7)) +
  labs(x = "power fitted to height residuals", y = NULL, color = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  scale_y_discrete(limits = rev) +
  theme(legend.position = "none", strip.background = element_blank(), strip.placement = "outside")

ggplot(heightDiameterResults %>% filter(responseVariable == "height", fitting != "gnls")) +
  geom_point(aes(x = bias, y = species, color = species)) +
  coord_cartesian(xlim = c(-3, 2.5)) +
  guides(color = "none") +
  labs(x = "nlrob() or gsl_nls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.position = "none") +
ggplot(heightDiameterResults %>% filter(responseVariable == "height", fitting == "gnls")) +
  geom_point(aes(x = bias, y = species, color = species)) +
  coord_cartesian(xlim = c(-3, 2.5)) +
  guides(color = "none") +
  labs(x = "gnls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.position = "none")
  