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
  group_by(responseVariable, form) %>%
  mutate(species = factor(species, labels = c("Douglas-fir", "western redcedar", "western hemlock", "red alder", "bigleaf maple", "Oregon myrtle", "other"),  levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         speciesFraction = recode(species, "Douglas-fir" = 0.750, "red alder" = 0.101, "western hemlock" = 0.056, "bigleaf maple" = 0.029, "Oregon myrtle" = 0.025, "western redcedar" = 0.013, "other" = 0.017),
         weighting = if_else(fitting == "nlrob", "reweighted", "fixed weights"),
         weightedMae = sum(speciesFraction * if_else(is.na(mae), 100, mae)),
         sizeShapeAlpha = as.factor(if_else(significant == TRUE, weighting, "not significant")))
heightDiameterAccuracyLevels = heightDiameterResults %>% summarize(weightedMae = weightedMae[1], .groups = "keep") %>%
  arrange(weightedMae)

print(heightDiameterResults %>% filter(str_detect(form, "Gnls") == FALSE, str_detect(form, "RelHt") == FALSE) %>% 
        group_by(responseVariable, species) %>% 
        slice_min(pae, n = 4) %>% 
        select(-responseVariable, -biasNR, -biasPl, -rmseNR, -rmsePl, -nseNR, -nsePl, -pearson, -pearsonNR, -pearsonPl, -aic, -bic), n = 60)
#write_xlsx(heightDiameterResults, "trees/height-diameter/HtDia results.xlsx")

## Figure 1: dataset summary
dbhQuantiles = liveTrees2016 %>% mutate(diameterClass = 2.5 * (ceiling(DBH / 2.5) - 0.5)) %>% group_by(diameterClass) %>%
  summarize(count = n(), quantiles = c("min", "q025", "q10", "q25", "median", "q75", "q90", "q975", "max"), height = quantile(TotalHt, probs = c(0, 0.025, 0.10, 0.25, 0.5, 0.75, 0.90, 0.975, 1), na.rm = TRUE), mean = mean(TotalHt, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = quantiles, values_from = height)
heightQuantiles = liveTrees2016 %>% mutate(heightClass = 1 * (ceiling(TotalHt / 1) - 0.5)) %>% group_by(heightClass) %>%
  summarize(count = n(), quantiles = c("min", "q025", "q10", "q25", "median", "q75", "q90", "q975", "max"), dbh = quantile(DBH, probs = c(0, 0.025, 0.10, 0.25, 0.5, 0.75, 0.90, 0.975, 1), na.rm = TRUE), mean = mean(DBH, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = quantiles, values_from = dbh)

ggplot() +
  geom_bin_2d(aes(x = DBH, y = TotalHt, fill = ..count..), liveTrees2016 %>% filter(is.na(TotalHt) == FALSE), binwidth = c(2.5, 1)) +
  geom_path(aes(x = diameterClass, y = mean, color = "mean height", linetype = "mean height"), dbhQuantiles %>% filter(count > 10), na.rm = TRUE) +
  geom_path(aes(x = diameterClass, y = median, color = "median height", linetype = "median height"), dbhQuantiles %>% filter(count > 10), na.rm = TRUE) +
  geom_path(aes(x = mean, y = heightClass, color = "mean DBH", linetype = "mean DBH"), heightQuantiles %>% filter(count > 10), na.rm = TRUE) +
  geom_path(aes(x = median, y = heightClass, color = "median DBH", linetype = "median DBH"), heightQuantiles %>% filter(count > 10), na.rm = TRUE) +
  annotate("text", x = 0, y = 82.13, label = "a)", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 80)) +
  labs(x = "DBH, cm", y = "height, m, of unbroken stem", color = NULL, fill = "trees\nmeasured", linetype = NULL) +
  guides(color = guide_legend(order = 1), fill = guide_colorbar(order = 2), linetype = guide_legend(order = 1)) +
  scale_color_manual(breaks = c("mean height", "median height", "mean DBH", "median DBH"), labels = c("mean\nheight", "median\nheight", "mean\nDBH", "median\nDBH"), values = c("green2", "green2", "burlywood2", "burlywood2")) +
  scale_fill_viridis_c(breaks = c(1, 10, 100, 330), trans = "log10") +
  scale_linetype_manual(breaks = c("mean height", "median height", "mean DBH", "median DBH"), labels = c("mean\nheight", "median\nheight", "mean\nDBH", "median\nDBH"), values = c("solid", "longdash", "solid", "longdash")) +
  theme(legend.key.height = unit(0.95, "line"), legend.justification = c(1, 0), legend.position = c(1, 0.02), legend.spacing.y = unit(0.25, "line")) +
ggplot(dbhQuantiles) +
  geom_ribbon(aes(x = diameterClass, ymin = 100 * (q025 - mean) / mean, ymax = 100 * (q975 - mean) / mean, alpha = "95% probability"), fill = "forestgreen") +
  geom_ribbon(aes(x = diameterClass, ymin = 100 * (q10 - mean) / mean, ymax = 100 * (q90 - mean) / mean, alpha = "80% probability"), fill = "forestgreen") +
  geom_ribbon(aes(x = diameterClass, ymin = 100 * (q25 - mean) / mean, ymax = 100 * (q75 - mean) / mean, alpha = "50% probability"), fill = "forestgreen") +
  geom_segment(x = 0, xend = 185, y = 0, yend = 0, color = "grey35") +
  annotate("text", x = 0, y = 154, label = "b)", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 196), ylim = c(-50, 150)) +
  scale_alpha_manual(breaks = c("95% probability", "80% probability", "50% probability"), values = c(0.1, 0.2, 0.5)) +
  labs(x = "DBH, cm", y = "departure from mean height, %", alpha = NULL) +
  theme(legend.position = "none") +
ggplot(heightQuantiles) +
  geom_ribbon(aes(x = heightClass, ymin = 100 * (q025 - mean) / mean, ymax = 100 * (q975 - mean) / mean, alpha = "95% probability"), fill = "burlywood4") +
  geom_ribbon(aes(x = heightClass, ymin = 100 * (q10 - mean) / mean, ymax = 100 * (q90 - mean) / mean, alpha = "80% probability"), fill = "burlywood4") +
  geom_ribbon(aes(x = heightClass, ymin = 100 * (q25 - mean) / mean, ymax = 100 * (q75 - mean) / mean, alpha = "50% probability"), fill = "burlywood4") +
  geom_segment(x = 0, xend = 77.5, y = 0, yend = 0, color = "grey35") +
  annotate("text", x = 0, y = 154, label = "c)", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 80), ylim = c(-50, 150)) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  scale_alpha_manual(breaks = c("95% probability", "80% probability", "50% probability"), values = c(0.1, 0.2, 0.5)) +
  labs(x = "height, m", y = "departure from mean DBH, %", alpha = NULL) +
  theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.98)) +
plot_layout(nrow = 1, ncol = 3, widths = c(260, 200, 200))
#ggsave("trees/height-diameter/Figure 1 height-diameter distribution.png", height = 10, width = 22, units = "cm", dpi = 150)


## Figure 2: height-diameter error summary
heightFromDiameterResults = heightDiameterResults %>% filter(responseVariable == "DBH") %>%
  mutate(form = factor(form, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "DBH"))$form))
diameterFromHeightResults = heightDiameterResults %>% filter(responseVariable == "height") %>%
  mutate(form = factor(form, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "height"))$form))

ggplot(heightFromDiameterResults %>% filter(str_detect(form, "RelHt") == FALSE)) +
  geom_point(aes(x = pae, y = form, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "mean absolute error, % DBH", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
ggplot(diameterFromHeightResults %>% filter(str_detect(form, "BAL") == FALSE)) +
  geom_point(aes(x = pae, y = form, color = species, alpha = sizeShapeAlpha, shape = sizeShapeAlpha, size = sizeShapeAlpha), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "mean absolute error, % tree height", y = NULL, color = NULL, alpha = NULL, shape = NULL, size = NULL) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  guides(color = guide_legend(order = 1, ncol = 4), alpha = guide_legend(order = 2, ncol = 2), shape = guide_legend(order = 2, ncol = 2), size = guide_legend(order = 2, ncol = 2)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  scale_alpha_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(0.75, 0.75, 0.6), drop = FALSE) &
  scale_shape_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(16, 18, 3), drop = FALSE) &
  scale_size_manual(breaks = c("reweighted", "fixed weights", "not significant"), values = c(1.5, 1.9, 1.4), drop = FALSE) &
  theme(legend.key.size = unit(0.2, "line"), legend.justification = "left", legend.position = "bottom")
#ggsave("trees/height-diameter/Figure 2 MAE.png", height = 11, width = 22, units = "cm", dpi = 150)

## Figure 3: Douglas-fir, red alder, and western hemlock regression comparision
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
#ggsave("trees/height-diameter/Figure 3 PSME-ALRU2-TSHE curves.png", height = 20, width = 17, units = "cm")

## Figure 4: bigleaf maple, Oregon myrtle, and western redcedar regression comparison
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
#ggsave("trees/height-diameter/Figure 4 ACMA3-UMCA-THPL curves.png", height = 20, width = 17, units = "cm")

## Figure 5: minority species merged regression comparison
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


## variance structure summary from GNLS
ggplot(heightDiameterResults %>% filter(is.na(power) == FALSE)) +
  geom_point(aes(x = 2*power, y = form, color = species)) +
  guides(color = guide_legend(ncol = 7)) +
  labs(x = "estimated variance power", y = NULL, color = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.position = "bottom")
  
heightDiameterResults %>% group_by(species) %>% 
  summarize(minPower = 2 * min(power, na.rm = TRUE), meanPower = 2 * mean(power, na.rm = TRUE), maxPower = 2 * max(power, na.rm = TRUE))
