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
                                  otherDiameterFromHeightResults,
                                  otherDiameterFromHeightResultsGnls) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other species"), levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         speciesFraction = recode(species, "Douglas-fir" = 0.750, "red alder" = 0.101, "western hemlock" = 0.056, "bigleaf maple" = 0.029, "Oregon myrtle" = 0.025, "western redcedar" = 0.013, "other species" = 0.017),
         isBaseForm = (str_detect(name, "AA\\+T") == FALSE) & (str_detect(name, "BA\\+L") == FALSE) & (str_detect(name, "physio") == FALSE) & (str_detect(name, "RelHt") == FALSE) & (fitting != "gnls"),
         weighting = if_else(fitting %in% c("gnls", "nlrob"), "reweighted", "fixed weights"),
         sizeShapeAlpha = as.factor(if_else(significant == TRUE, weighting, "not significant"))) %>%
  group_by(responseVariable, species) %>%
  mutate(deltaAicN = aic/n - min(aic / n, na.rm = TRUE),
         deltaAicNrank = dense_rank(deltaAicN),
         maeRank = dense_rank(mae),
         nseRank = dense_rank(desc(nse)),
         paeRank = dense_rank(pae),
         pearsonRank = dense_rank(desc(pearson)),
         rmseRank = dense_rank(rmse)) %>%
  group_by(responseVariable, name) %>%
  mutate(weightedMae = sum(speciesFraction * if_else(is.na(mae), 100, mae))) %>%
  ungroup()
heightDiameterAccuracyLevels = heightDiameterResults %>% group_by(responseVariable, name) %>%
  summarize(weightedMae = weightedMae[1], .groups = "keep") %>%
  arrange(weightedMae)
#write_xlsx(heightDiameterResults, "trees/height-diameter/HtDia results.xlsx")


## summaries for Results
#unique((heightDiameterResults %>% filter(isBaseForm))$name)
heightDiameterResults %>% filter(isBaseForm) %>% group_by(responseVariable, species) %>% 
  mutate(deltaAicNgamRank = deltaAicNrank[which(name == "GCV GAM")],
         maeGamRank = maeRank[which(name == "GCV GAM")],
         nseGamRank = nseRank[which(name == "GCV GAM")],
         pearsonGamRank = pearsonRank[which(name == "GCV GAM")],
         rmseGamRank = rmseRank[which(name == "GCV GAM")]) %>%
  summarize(n = n(),
            deltaAicN = sum(deltaAicNrank < deltaAicNgamRank, na.rm = TRUE),
            mae = sum(maeRank < maeGamRank, na.rm = TRUE),
            nse = sum(nseRank < nseGamRank, na.rm = TRUE),
            pearson = sum(pearsonRank < pearsonGamRank, na.rm = TRUE),
            rmse = sum(rmseRank < rmseGamRank, na.rm = TRUE),
            .groups = "drop_last") %>%
  summarize(deltaAicN = 100 * sum(deltaAicN) / sum(n),
            mae = 100 * sum(mae) / sum(n),
            nse = 100 * sum(nse) / sum(n),
            pearson = 100 * sum(pearson) / sum(n),
            rmse = 100 * sum(rmse) / sum(n))
#  summarize(deltaAicN = 100 * sum(deltaAicNrank < deltaAicNgamRank, na.rm = TRUE) / n(),
#            mae = 100 * sum(maeRank < maeGamRank, na.rm = TRUE) / n(),
#            nse = 100 * sum(nseRank < nseGamRank, na.rm = TRUE) / n(),
#            pearson = 100 * sum(pearsonRank < pearsonGamRank, na.rm = TRUE) / n(),
#            rmse = 100 * sum(rmseRank < rmseGamRank, na.rm = TRUE) / n(),
#            .groups = "drop")

heightDiameterParameters %>% filter(fitting %in% c("lm", "nlrob", "gsl_nls")) %>%
  #group_by(responseVariable) %>%
  summarise(anyP = sum((is.na(a1p) == FALSE) | (is.na(a2p) == FALSE) | (is.na(a3p) == FALSE) | (is.na(a4p) == FALSE) | (is.na(b1p) == FALSE) | (is.na(b2p) == FALSE) | (is.na(b3p) == FALSE)),
            a1 = sum(is.na(a1) == FALSE), a1p = sum(is.na(a1p) == FALSE),
            a2 = sum(is.na(a2) == FALSE), a2p = sum(is.na(a2p) == FALSE),
            a3 = sum(is.na(a3) == FALSE), a3p = sum(is.na(a3p) == FALSE),
            a4 = sum(is.na(a4) == FALSE), a4p = sum(is.na(a4p) == FALSE), # plantation effects not tested for a5 to a7
            b1 = sum(is.na(b1) == FALSE), b1p = sum(is.na(b1p) == FALSE),
            b2 = sum(is.na(b2) == FALSE), b2p = sum(is.na(b2p) == FALSE),
            b3 = sum(is.na(b3) == FALSE), b3p = sum(is.na(a3p) == FALSE)) %>%
  mutate(total = a1 + a2 + a3 + a4 + b1 + b2 + b3, 
         totalP = a1p + a2p + a3p + a4p + b1p + b2p + b3p,
         hasPpct = 100 * anyP/a1)

# AIC selection of generalizing predictors
heightDiameterResults %>% filter(fitting %in% c("nlrob", "gsl_nls")) %>% 
  mutate(baseName = str_remove(name, " AA\\+T| BA\\+L| physio| RelHt")) %>%
  group_by(responseVariable, species, baseName) %>%
  mutate(n = n(), 
         deltaAicNgeneralized = deltaAicN - max(isBaseForm * deltaAicN),
         isPhysio = str_detect(name, "physio")) %>%
  filter(n > 1) %>%
  #select(responseVariable, species, isBaseForm, baseName, name, deltaAicNgeneralized) %>%
  #arrange(desc(responseVariable), species, baseName, name) %>%
  group_by(responseVariable) %>%
  filter(isBaseForm == FALSE) %>%
  summarize(n = n(), 
            base = sum(deltaAicNgeneralized > 0), 
            generalized = sum(deltaAicNgeneralized < 0), 
            physioDeltaAicNmed = median(if_else(isPhysio, deltaAicNgeneralized, NA_real_), na.rm = TRUE), 
            standDeltaAicNmed = median(if_else(isPhysio, NA_real_, deltaAicNgeneralized), na.rm = TRUE)) %>%
  mutate(basePct = 100 * base / n, generalizedPct = 100 * generalized / n) %>%
  arrange(desc(responseVariable))

# parameter counts
#heightDiameterParameterCounts = heightDiameterParameters %>% rowwise() %>% mutate(parameters = ncol(heightDiameterParameters) - sum(is.na(c_across(is.numeric)))) %>%
#  select(responseVariable, species, name, fitting, parameters) # slow
#print(heightDiameterParameterCounts %>% filter(fitting == "gam"), n = 50)


# summarize fitting success for Section S1
heightDiameterResults %>% 
  summarize(nlrob = sum(fitting == "nlrob", na.rm = TRUE), 
            gslNls = sum(fitting == "gsl_nls", na.rm = TRUE), 
            gnls = sum((fitting == "gnls") | (is.na(fitting) & str_detect(name, "GNLS")), na.rm = TRUE), 
            lm = sum(fitting == "lm", na.rm = TRUE), 
            gam = sum(fitting == "gam", na.rm = TRUE),
            fail = sum(is.na(fitting)), 
            uniqueNonlinear = n() - gnls - lm - gam) %>%
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
heightFromDiameterResults = heightDiameterResults %>% filter(responseVariable == "height", str_detect(name, "GNLS") == FALSE, str_detect(name, "RelHt") == FALSE, name != "GCV GAM BA+L physio") %>%
  group_by(responseVariable, species) %>%
  mutate(name = factor(name, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "height"))$name),
         deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) # recalculate ΔAICn with relative height forms excluded, requires regrouping by species

ggplot(heightFromDiameterResults) +
  geom_point(aes(x = pae, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(heightFromDiameterResults) +
  geom_point(aes(x = rmse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 11)) +
  labs(x = "RMSE, m", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_discrete(labels = NULL) +
ggplot(heightFromDiameterResults) +
  geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM", yend = "GCV GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM BA+L", yend = "GCV GAM BA+L", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM physio", yend = "GCV GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_point(aes(x = deltaAicN, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  labs(x = bquote("normalized "*Delta*"AIC"), y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  coord_cartesian(xlim = c(0, 0.265)) + # exclude high AIC of linear, parabolic, and Douglas-fir power+Curtis fits to avoid squashing of nonlinear differences
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
#ggsave("trees/height-diameter/Figure 2 height accuracy.png", height = 12, width = 20, units = "cm", dpi = 150)


## Figure 3: diameter-height error summary
diameterFromHeightResults = heightDiameterResults %>% filter(responseVariable == "DBH", str_detect(name, "BA\\+L") == FALSE, str_detect(name, "GNLS") == FALSE, name != "GCV GAM AA+T physio") %>%
  group_by(responseVariable, species) %>%
  mutate(name = factor(name, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "DBH"))$name),
         deltaAicN = aic/n - min(aic / n, na.rm = TRUE)) # recalculate ΔAICn with BAL names excluded

ggplot(diameterFromHeightResults) +
  geom_point(aes(x = pae, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "MAE, %", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(diameterFromHeightResults) +
  geom_point(aes(x = rmse, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "RMSE, cm", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
  scale_y_discrete(labels = NULL) +
ggplot(diameterFromHeightResults) +
  geom_segment(x = 0.255, xend = 0.275, y = "Schnute", yend = "Schnute", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "linear", yend = "linear", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "parabolic", yend = "parabolic", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM", yend = "GCV GAM", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM AA+T", yend = "GCV GAM AA+T", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_segment(x = 0.255, xend = 0.275, y = "GCV GAM physio", yend = "GCV GAM physio", arrow = arrow(length = unit(0.2, "line"), type = "closed"), color = "grey70", linewidth = 0.4) +
  geom_point(aes(x = deltaAicN, y = name, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 0.265)) +
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
#ggsave("trees/height-diameter/Figure 3 diameter accuracy.png", height = 12, width = 20, units = "cm", dpi = 150)


## Figure 4: comparison of accuracy metrics
accuracyCorrelation = bind_rows(as.data.frame(cor(heightDiameterResults %>% filter(responseVariable == "height") %>% select(maeRank, rmseRank, deltaAicNrank, nseRank, pearsonRank), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "height"),
                                as.data.frame(cor(heightDiameterResults %>% filter(responseVariable == "DBH") %>% select(maeRank, rmseRank, deltaAicNrank, nseRank, pearsonRank), use = "pairwise.complete.obs")) %>%
                                  rownames_to_column("metricX") %>% gather("metricY", "correlation", -metricX) %>%
                                  mutate(responseVariable = "DBH")) %>%
  mutate(metricX = factor(metricX, levels = c("maeRank", "rmseRank", "deltaAicNrank", "nseRank", "pearsonRank"), labels = c("MAE", "RMSE", "normalized\nΔAIC", "model\nefficiency", "Pearson's\nR")),
         metricY = factor(metricY, levels = c("maeRank", "rmseRank", "deltaAicNrank", "nseRank", "pearsonRank"), labels = c("MAE", "RMSE", "normalized\nΔAIC", "model\nefficiency", "Pearson's\nR")))

ggplot(accuracyCorrelation %>% filter(responseVariable == "height")) + 
  coord_equal() +
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  labs(x = NULL, y = NULL, fill = "correlation", title = "a) height imputation") +
ggplot(accuracyCorrelation %>% filter(responseVariable == "DBH")) + 
  geom_raster(aes(x = metricX, y = metricY, fill = correlation)) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "correlation", title = "b) DBH imputation") +
  scale_y_discrete(labels = NULL) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  scale_fill_viridis_c(limits = c(-0.1, 1)) &
  theme(legend.spacing.y = unit(0.5, "line"), title = element_text(size = 8))
ggsave("trees/height-diameter/Figure 4 accuracy metric correlation.png", height = 9, width = 20, units = "cm")


## Figure 5: preferred model forms
nonlinearForms = heightDiameterResults %>% filter(str_detect(name, "GNLS") == FALSE, (responseVariable != "DBH") | (str_detect(name, "BA\\+L") == FALSE), (responseVariable == "DBH") | (str_detect(name, "RelHt") == FALSE)) %>%
  group_by(responseVariable, species)

preferredForms = bind_rows(nonlinearForms %>% filter(responseVariable == "height") %>% slice_min(maeRank, n = 4) %>% arrange(desc(maeRank)) %>% mutate(criteria = "mae", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name, speciesFraction),
                           nonlinearForms %>% filter(responseVariable == "DBH") %>% slice_min(rmseRank, n = 4) %>% arrange(desc(rmseRank)) %>% mutate(criteria = "rmse", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name, speciesFraction),
                           nonlinearForms %>% slice_min(deltaAicNrank, n = 4) %>% arrange(desc(deltaAicNrank)) %>% mutate(criteria = "aicN", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name, speciesFraction),
                           nonlinearForms %>% slice_min(nseRank, n = 4) %>% arrange(desc(nseRank)) %>% mutate(criteria = "nse", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name, speciesFraction),
                           nonlinearForms %>% slice_min(pearsonRank, n = 4) %>% arrange(desc(pearsonRank)) %>% mutate(criteria = "pearson", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name, speciesFraction))

ggplot(preferredForms %>% filter(responseVariable == "height", criteria == "mae")) +
  geom_bar(aes(x = 100 * after_stat(count / sum(count)), y = fct_reorder(name, speciesFraction, sum), fill = species, group = species, weight = speciesFraction)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "a) height models by MAE") +
ggplot(preferredForms %>% filter(responseVariable == "height", criteria == "aicN")) +
  geom_bar(aes(x = 100 * after_stat(count / sum(count)), y = fct_reorder(name, speciesFraction, sum), fill = species, group = species, weight = speciesFraction)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "b) height models by ΔAIC") +
ggplot(preferredForms %>% filter(responseVariable == "DBH", criteria == "rmse")) +
  geom_bar(aes(x = 100 * after_stat(count / sum(count)), y = fct_reorder(name, speciesFraction, sum), fill = species, group = species, weight = speciesFraction)) +
  labs(x = "frequency of selection, %", y = NULL, fill = NULL, title = "c) DBH models by RMSE") +
ggplot(preferredForms %>% filter(responseVariable == "DBH", criteria == "aicN")) +
  geom_bar(aes(x = 100 * after_stat(count / sum(count)), y = fct_reorder(name, speciesFraction, sum), fill = species, group = species, weight = speciesFraction)) +
  labs(x = "frequency of selection, %", y = NULL, fill = NULL, title = "d) DBH models by ΔAIC") +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 2, ncol = 2, guides = "collect") &
  coord_cartesian(xlim = c(0, 25)) &
  guides(fill = guide_legend(byrow = TRUE)) &
  scale_fill_manual(breaks = levels(preferredForms$species), limits = levels(preferredForms$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  theme(legend.position = "bottom", title = element_text(size = 8))
ggsave("trees/height-diameter/Figure 5 model selection.png", height = 12, width = 20, units = "cm")


## Figure 6: Douglas-fir, red alder, and western hemlock regression comparison
# Generalized fits shown in Figures 6, 7, and 8 follow from selections in Figure 5. Overlaid base forms are
# from the below query.
preferredBaseForms = bind_rows(nonlinearForms %>% filter(isBaseForm, responseVariable == "height") %>% slice_min(maeRank, n = 4) %>% arrange(desc(maeRank)) %>% mutate(criteria = "mae", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name),
                               nonlinearForms %>% filter(isBaseForm, responseVariable == "DBH") %>% slice_min(rmseRank, n = 4, na.rm = TRUE) %>% arrange(desc(rmseRank)) %>% mutate(criteria = "rmse", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name),
                               nonlinearForms %>% filter(isBaseForm) %>% slice_min(deltaAicNrank, n = 4) %>% arrange(desc(deltaAicNrank)) %>% mutate(criteria = "aicN", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name),
                               nonlinearForms %>% filter(isBaseForm) %>% slice_min(nseRank, n = 4) %>% arrange(desc(nseRank)) %>% mutate(criteria = "nse", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name),
                               nonlinearForms %>% filter(isBaseForm) %>% slice_min(pearsonRank, n = 4) %>% arrange(desc(pearsonRank)) %>% mutate(criteria = "pearson", rank = row_number()) %>% select(responseVariable, species, criteria, rank, name))
preferredBaseForms %>% filter(species == "western hemlock") %>% pivot_wider(names_from = criteria, values_from = name) %>%
  arrange(desc(responseVariable), species) 

ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = psme2016physio$DBH, y = psmeHeightFromDiameterChapmanRichardsBalPhysio$fitted.values, color = "Chapman-Richards BAL physio", group = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterGam$fitted.values, color = "GCV GAM", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterProdan$fitted.values, color = "Prodan", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  #geom_line(aes(x = psme2016$DBH, y = psmeHeightFromDiameterSibbesen$fitted.values, color = "Sibbesen", group = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007")) + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) Douglas-fir height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("GCV GAM", "Prodan", "Chapman-Richards BAL physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  theme(legend.key = element_rect(fill = alpha("white", 0.5)), legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = psmeDiameterFromHeightSharmaParton$fitted.values, y = psme2016$TotalHt, color = "modified Sharma-Parton", group = psme2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psmeDiameterFromHeightGam$fitted.values, y = psme2016$TotalHt, color = "GCV GAM", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psmeDiameterFromHeightRuark$fitted.values, y = psme2016$TotalHt, color = "Ruark", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007")) + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) Douglas-fir DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("GCV GAM", "Ruark", "modified Sharma-Parton", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.10, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alru2016physio$DBH, y = alruHeightFromDiameterChapmanRichardsBalPhysio$fitted.values, color = "Chapman-Richards BAL physio", group = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterChapmanRichards$fitted.values, color = "Chapman-Richards", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = alruHeightFromDiameterWeibull$fitted.values, color = "Weibull", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) red alder height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Weibull", "Chapman-Richards BAL physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  #scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen", "Sharma-Parton BAL physiological", "Temesgen et al. 2007"), values = c("#ac92eb", "#4dc1e8", "#a0d568", "#ffce54", "#ed5564", "grey65")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = alruDiameterFromHeightSibbesenFormPhysio$fitted.values, y = alru2016physio$TotalHt, color = "Chapman-Richards physio", group = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alruDiameterFromHeightChapmanRichards$fitted.values, y = alru2016$TotalHt, color = "Chapman-Richards", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alruDiameterFromHeightSibbesenForm$fitted.values, y = alru2016$TotalHt, color = "Sibbesen form", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) red alder DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen form", "Chapman-Richards physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = tshe2016physio$DBH, y = tsheHeightFromDiameterChapmanRichardsBalPhysio$fitted.values, color = "Chapman-Richards BAL physio", group = tshe2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterGam$fitted.values, color = "GCV GAM", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = tshe2016$DBH, y = tsheHeightFromDiameterProdan$fitted.values, color = "Prodan", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "e) western hemlock height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("GCV GAM", "Prodan", "Chapman-Richards BAL physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey50")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = tshe2016$DBH, y = tshe2016$TotalHt), alpha = 0.15, color = "grey25", na.rm = TRUE, shape = 16) +
  geom_line(aes(x = tsheDiameterFromHeightSibbesenFormAat$fitted.values, y = tshe2016$TotalHt, color = "Sibbesen form AA+T", group = tshe2016$isPlantation), alpha = 0.4) +
  geom_line(aes(x = tsheDiameterFromHeightChapmanRichards$fitted.values, y = tshe2016$TotalHt, color = "Chapman-Richards", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = tsheDiameterFromHeightSibbesenForm$fitted.values, y = tshe2016$TotalHt, color = "Sibbesen form", group = tshe2016$isPlantation, linetype = tshe2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), color = "Temesgen et al. 2007"), linetype = "longdash") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "f) western hemlock DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("Chapman-Richards", "Sibbesen form", "Sibbesen form AA+T", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey50")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(nrow = 3, ncol = 2, ) &
  guides(color = guide_legend(order = 1), linetype = "none") &
  scale_linetype_manual(breaks = c(FALSE, TRUE), labels = c("natural regeneration", "plantation"), values = c("solid", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/Figure 6 PSME-ALRU2-TSHE curves.png", height = 17, width = 20, units = "cm")

## Figure 7: bigleaf maple, Oregon myrtle, and western redcedar regression comparison
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
#ggsave("trees/height-diameter/Figure 7 ACMA3-UMCA-THPL curves.png", height = 20, width = 17, units = "cm")

## Figure 8: minority species merged regression comparison
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
#ggsave("trees/height-diameter/Figure 8 other species curves.png", height = 20/3 + 0.5, width = 17, units = "cm")


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
  scale_color_manual(breaks = levels(residualPower$species), limits = levels(residualPower$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  scale_shape_manual(breaks = c("b1", "b1p"), labels = c("base power", "plantations"), values = c(15, 17)) &
  theme(legend.position = "bottom")
ggsave("trees/height-diameter/Figure S06 residual power estimates.png", height = 14, width = 20, units = "cm", dpi = 150)


## Figures S7-12: Q-Q plots
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
  coord_cartesian(xlim = c(-2, 2)) +
  guides(color = "none") +
  labs(x = "nlrob() or gsl_nls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.position = "none") +
ggplot(heightDiameterResults %>% filter(responseVariable == "height", fitting == "gnls")) +
  geom_point(aes(x = bias, y = species, color = species)) +
  coord_cartesian(xlim = c(-2, 2)) +
  guides(color = "none") +
  labs(x = "gnls() bias, m height", y = NULL, color = NULL, shape = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  scale_y_discrete(labels = NULL) +
  theme(legend.position = "none")
  