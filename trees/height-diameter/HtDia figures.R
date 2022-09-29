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
  group_by(responseVariable, method) %>%
  mutate(species = factor(species, levels = c("PSME", "THPL", "TSHE", "ALRU2", "ACMA3", "UMCA", "other")),
         weightedMae = sum(recode(species, PSME = 0.750, ALRU2 = 0.101, TSHE = 0.056, ACMA3 = 0.029, UMCA = 0.025, THPL = 0.013, other = 0.017) * if_else(is.na(mae), 100, mae)))
heightDiameterAccuracyLevels = heightDiameterResults %>% summarize(weightedMae = weightedMae[1], .groups = "keep") %>%
  arrange(weightedMae)

print(heightDiameterResults %>% filter(str_detect(method, "Gnls") == FALSE, str_detect(method, "RelHt") == FALSE) %>% 
        group_by(responseVariable, species) %>% 
        slice_min(pae, n = 4) %>% 
        select(-responseVariable, -biasNR, -biasPl, -rmseNR, -rmsePl, -nseNR, -nsePl, -pearson, -pearsonNR, -pearsonPl, -aic, -bic), n = 60)
#write_xlsx(heightDiameterResults, "trees/height-diameter/HtDia results.xlsx")

## Figure 1: height-diameter error summary
heightFromDiameterResults = heightDiameterResults %>% filter(responseVariable == "DBH") %>%
  mutate(method = factor(method, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "DBH"))$method))
diameterFromHeightResults = heightDiameterResults %>% filter(responseVariable == "height") %>%
  mutate(method = factor(method, levels = (heightDiameterAccuracyLevels %>% filter(responseVariable == "height"))$method))

ggplot(heightFromDiameterResults %>% filter(str_detect(method, "RelHt") == FALSE)) +
  geom_point(aes(x = pae, y = method, color = species), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "mean absolute error, % DBH", y = NULL, color = NULL, shape = NULL) +
ggplot(diameterFromHeightResults %>% filter(str_detect(method, "BAL") == FALSE)) +
  geom_point(aes(x = pae, y = method, color = species), na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "mean absolute error, % tree height", y = NULL, color = NULL, shape = NULL) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  guides(color = guide_legend(ncol = 7)) &
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) &
  theme(legend.justification = "center", legend.position = "bottom")
#ggsave("trees/height-diameter/Figure 1 MAE.png", height = 10, width = 20, units = "cm", dpi = 200)

## Figure 2: Douglas-fir, red alder, and western hemlock regression comparision
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
#ggsave("trees/height-diameter/Figure 2 PSME-ALRU2-TSHE curves.png", height = 20, width = 17, units = "cm")

## Figure 3: bigleaf maple, Oregon myrtle, and western redcedar regression comparison
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
#ggsave("trees/height-diameter/Figure 3 ACMA3-UMCA-THPL curves.png", height = 20, width = 17, units = "cm")

## Figure 3: minority species merged regression comparison
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
  geom_point(aes(x = 2*power, y = method, color = species)) +
  guides(color = guide_legend(ncol = 7)) +
  labs(x = "estimated variance power", y = NULL, color = NULL) +
  scale_color_manual(breaks = levels(heightFromDiameterResults$species), limits = levels(heightFromDiameterResults$species), values = c("forestgreen", "firebrick", "blue2", "red2", "green3", "mediumorchid1", "grey65")) +
  theme(legend.position = "bottom")
  
heightDiameterResults %>% group_by(species) %>% 
  summarize(minPower = 2 * min(power, na.rm = TRUE), meanPower = 2 * mean(power, na.rm = TRUE), maxPower = 2 * max(power, na.rm = TRUE))
