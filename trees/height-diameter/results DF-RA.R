# extension of results.R with figures for just Douglas-fir and red alder

# Douglas-fir and red alder height-diameter distribution
plot_exploratory(trees2016 %>% filter(speciesGroup == "DF", isLiveUnbroken), speciesLabel = "Douglas-fir", omitQuantiles = TRUE) +
plot_exploratory(trees2016 %>% filter(speciesGroup == "RA", isLiveUnbroken), plotLetters = c("b)"), speciesLabel = "red alder", omitQuantiles = TRUE) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, ncol = 2, guides = "collect")
#ggsave("trees/height-diameter/figures/DF-RA height-DBH distribution.png", height = 6, width = 16, units = "cm", dpi = 275)

# DBH AUCs
diameterFromHeightModelComparison = heightDiameterModelRanking %>% 
  filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "DBH"))$name)))
plot_auc_bank(diameterFromHeightModelComparison, omitMab = TRUE, xLimits = c("Douglas-fir", "red alder"))
#ggsave("trees/height-diameter/figures/DF-RA DBH AUCs.png", height = 14, width = 28, units = "cm", dpi = 175)

# height AUCs
heightFromDiameterModelComparison = heightDiameterModelRanking %>% 
  filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height") %>%
  mutate(name = factor(name, levels = rev((heightDiameterModelDisplaySort %>% filter(responseVariable == "height"))$name)))
plot_auc_bank(heightFromDiameterModelComparison, omitMab = TRUE, xLimits = c("Douglas-fir", "red alder"))
#ggsave("trees/height-diameter/figures/DF-RA height AUCs.png", height = 14, width = 28, units = "cm", dpi = 175)

# model fits
modelFitStats = primaryResults %>% filter(species %in% c("Douglas-fir", "red alder")) %>%
  group_by(responseVariable, species, name) %>%
  summarize(modelFit = is.na(nse[1]) == FALSE, .groups = "drop_last") %>% 
  summarize(fitPct = 100 * mean(modelFit), .groups = "drop") %>%
  mutate(species = factor(species, levels = levels(predictorVariableResults$species)))

ggplot(modelFitStats %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.4) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "a) DBH fits") +
  scale_y_discrete(limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableStats %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH")) +
  geom_col(aes(x = significantPct, y = species, fill = species), width = 0.4) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "b) DBH forms") +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", significant)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = NULL, y = NULL, color = NULL, title = "c) DBH") +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "d) ABA+T") +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", hasPhysio)) +
  geom_count(aes(x = nPhysioSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "e) DBH physiographic") +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", hasRelHt)) +
  geom_count(aes(x = significantRelHt, y = species, color = species)) + 
  labs(x = NULL, y = NULL, fill = NULL, title = "f) RelHt") +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = c(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(modelFitStats %>% filter(responseVariable == "height")) +
  geom_col(aes(x = fitPct, y = species, fill = species), width = 0.4) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "model forms\nfit, %", y = NULL, fill = NULL, title = "g) height fits") +
  scale_y_discrete(limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableStats %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height")) +
  geom_col(aes(x = significantPct, y = species, fill = species), width = 0.4) + 
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "significant\nfits, %", y = NULL, fill = NULL, title = "h) height forms") +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", significant)) +
  geom_violin(aes(x = meanAbsolutePercentPlantationEffect, y = species, color = species, group = species), draw_quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE, width = 1.2) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "mean absolute\nplantation effect, %", y = NULL, color = NULL, title = "i) height") +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", hasBasalArea)) +
  geom_count(aes(x = nBasalAreaSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.2, 2.2)) +
  labs(x = "predictors\nretained", y = NULL, fill = NULL, title = "j) BA+L") +
  scale_x_continuous(breaks = seq(0, 2), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", hasPhysio)) +
  geom_count(aes(x = nPhysioSignificant, y = species, color = species)) + 
  coord_cartesian(xlim = c(-0.3, 5.3)) +
  labs(x = "predictors\nretained", y = NULL, fill = NULL, title = "k) height physiographic") +
  scale_x_continuous(breaks = seq(0, 5), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
ggplot(predictorVariableResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", hasRelDbh)) +
  geom_count(aes(x = significantRelDbh, y = species, color = species)) + 
  labs(x = "predictors\nretained", y = NULL, fill = NULL, title = "l) RelDbh") +
  coord_cartesian(xlim = c(-0.3, 1.3)) +
  scale_x_continuous(breaks = seq(0, 1), minor_breaks = NULL) +
  scale_y_discrete(labels = NULL, limits = rev(c("Douglas-fir", "red alder"))) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 6, widths = c(1.7, 1.7, 2.2, 1.7, 2.7, 1.2), guides = "collect") &
  scale_color_manual(breaks = c("Douglas-fir", "red alder"), limits = c("Douglas-fir", "red alder"), values = c("forestgreen", "red2")) &
  scale_fill_manual(breaks = c("Douglas-fir", "red alder"), limits = c("Douglas-fir", "red alder"), values = c("forestgreen", "red2")) &
  scale_size_area(max_size = 5.0) &
  theme(legend.position = "none")
#ggsave("trees/height-diameter/figures/DF-RA predictor variables.png", height = 8, width = 19, units = "cm", dpi = 200)

# model efficiency
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 40)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "a) base DBH prediction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", isBaseForm)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 40)) +
  labs(x = NULL, y = "fraction of fits, %", fill = NULL, title = "b) base height prediction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "DBH", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 40)) +
  labs(x = "model efficiency", y = NULL, fill = NULL, title = "c) generalized DBH prediction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
ggplot(primaryResults %>% filter(species %in% c("Douglas-fir", "red alder"), responseVariable == "height", isBaseForm == FALSE)) +
  geom_histogram(aes(x = if_else(nse > -0.05, nse, NA_real_), y = 100 * after_stat(count / sum(count)), fill = species), binwidth = 0.025, na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 40)) +
  labs(x = "model efficiency", y = "fraction of fits, %", fill = NULL, title = "d) generalized height prediction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 2, guides = "collect") &
  guides(color = "none", fill = guide_legend(ncol = 2)) &
  scale_fill_manual(breaks = c("Douglas-fir", "red alder"), limits = c("Douglas-fir", "red alder"), values = c("forestgreen", "red2")) &
  theme(legend.key.size = unit(0.6, "line"), legend.position = "bottom", legend.text = element_text(margin = margin(l = -2.5, r = 6.5)))
#ggsave("trees/height-diameter/figures/DF-RA model efficiency.png", height = 9, width = 19, units = "cm", dpi = 200)

# model fits
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gamAbatPhysio), y = psme2016physio$TotalHt, color = "REML GAM ABA+T physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$gam), y = psme2016$TotalHt, color = "REML GAM", group = psme2016$isPlantation, linetype = psme2016$isPlantation), orientation = "y") +
  geom_line(aes(x = predict(psmeDiameterFromHeightPreferred$ruark), y = psme2016$TotalHt, color = "Ruark", group = psme2016$isPlantation, linetype = psme2016$isPlantation), orientation = "y") +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "a) Douglas-fir DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Ruark", "REML GAM ABA+T physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = psme2016$DBH, y = psme2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = psme2016physio$DBH, y = predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh), color = "Sharma-Parton BA+L RelDbh physio", group = psme2016physio$isPlantation, linetype = psme2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$gam), color = "REML GAM", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = predict(psmeHeightFromDiameterPreferred$sibbesen), color = "Sibbesen", group = psme2016$isPlantation, linetype = psme2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "b) Douglas-fir height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = NULL, y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Sibbesen", "Sharma-Parton BA+L RelDbh physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey60")) +
  theme(legend.key = element_rect(fill = alpha("white", 0.5)), legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$gamAbatPhysioRelHt), y = alru2016physio$TotalHt, color = "REML GAM ABA+T RelHt physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4, orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$ruark), y = alru2016$TotalHt, color = "Ruark", group = alru2016$isPlantation, linetype = alru2016$isPlantation), orientation = "y") +
  geom_line(aes(x = predict(alruDiameterFromHeightPreferred$gam), y = alru2016$TotalHt, color = "REML GAM", group = alru2016$isPlantation, linetype = alru2016$isPlantation), orientation = "y") +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "c) red alder DBH", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = "height, m", color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Ruark", "REML GAM ABA+T RelHt physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "cyan", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
ggplot() +
  geom_point(aes(x = alru2016$DBH, y = alru2016$TotalHt), alpha = 0.08, color = "grey25", na.rm = TRUE, shape = 16, size = 1.2) +
  geom_line(aes(x = alru2016physio$DBH, y = predict(alruHeightFromDiameterPreferred$sharmaPartonBalPhysio), color = "Sharma-Parton BA+L physio", group = alru2016physio$isPlantation, linetype = alru2016physio$isPlantation), alpha = 0.4) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$gam), color = "REML GAM", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightFromDiameterPreferred$ratkowsky), color = "Ratkowsky", group = alru2016$isPlantation, linetype = alru2016$isPlantation)) +
  geom_line(aes(x = psme2016$DBH, y = 1.3 + exp(5.7567 - 6.7792*psme2016$DBH^-0.2795), linetype = "reference curve"), color = "grey70") + # Temesgen et al. 2007, Eq. 4
  annotate("text", x = 0, y = 85, label = "d) red alder height", hjust = 0, size = 3.5) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 85)) +
  guides(linetype = "none") +
  labs(x = "DBH, cm", y = NULL, color = NULL, linetype = NULL) +
  scale_color_manual(breaks = c("REML GAM", "Ratkowsky", "Sharma-Parton BA+L physio", "Temesgen et al. 2007"), values = c("dodgerblue2", "red2", "green2", "grey70")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03)) +
get_preferred_model_linetype_legend() +
plot_annotation(theme = theme(plot.margin = margin(1, 1, 1, 1, "pt"))) +
plot_layout(design = "12\n34\n55", heights = c(1, 1, 0)) &
  scale_linetype_manual(breaks = c(FALSE, TRUE, "reference curve"), labels = c("natural regeneration", "plantation", "reference curve"), values = c("solid", "longdash", "dashed")) &
  scale_y_continuous(breaks = seq(0, 100, by = 20))
#ggsave("trees/height-diameter/figures/DF-RA model fits.png", height = 9/16 * 22, width = 22, units = "cm", dpi = 200)