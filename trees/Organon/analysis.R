library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(readxl)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             legend.key.height = unit(0.85, "line"),
                             legend.spacing.y = unit(0, "line"),
                             legend.title = element_text(size = 10),
                             panel.border = element_blank(),
                             plot.title = element_text(size = 10)))

standTrajectoriesFile ="trees/Organon/Elliott stand trajectories 2016-2116.csv"
#standTrajectoriesFile ="trees/Organon/Elliott stand trajectories 80 percent site index.csv"
standTrajectories = left_join(read_csv(standTrajectoriesFile, col_types = cols(.default = "d", financialScenario = "c", regenMinCostSystem = "c", thinChainsawCrewWithWheeledHarvester = "c", thinForwardingMethod = "c", regenChainsawCrewWithFellerBuncherAndGrappleSwingYarder = "c", regenChainsawCrewWithFellerBuncherAndGrappleYoader = "c", regenChainsawCrewWithTrackedHarvester = "c", regenChainsawCrewWithWheeledHarvester = "c")) %>%
                                mutate(mai = standingMbfh / standAge), # MBF/ha
                              read_xlsx("trees/Organon/Elliott Organon cruise records 2015-16.xlsx", sheet = "stands") %>% rename(stand = id) %>% select(-age),
                              by = c("stand"))

stands2022 %>% group_by(is.na(plotsInStand)) %>% summarize(area = sum(standArea))

# MAI
ggplot() +
  geom_segment(aes(x = c(40, 65), y = 0, xend = c(40, 65), yend = 20000), color = "grey50", linetype = "longdash", linewidth = 0.4) +
  geom_bin_2d(aes(x = standAge, y = 1000 * mai, weight = area), standTrajectories %>% group_by(stand) %>% slice_min(standAge, n = 1), binwidth = c(2, 50)) +
  labs(x = "stand age, years", y = bquote("cruised MAI in 2016, BF ha"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = 1000 * mai, weight = area), standTrajectories %>% group_by(stand) %>% slice_max(standAge, n = -1), binwidth = c(2, 50)) +
  geom_segment(aes(x = c(40, 65), y = 0, xend = c(40, 65), yend = 20000), color = "grey50", linetype = "longdash", linewidth = 0.4) +
  labs(x = "stand age, years", y = bquote("Organon SWO predicted MAI, BF ha"^-1), fill = "area, ha") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(widths = c(200, 300), guides = "collect") &
  coord_cartesian(ylim = c(0, 5000)) &
  scale_fill_viridis_c(breaks = c(1, 10, 100, 500), limits = c(1, 500), trans = "log10") &
  scale_y_continuous(labels = scales::label_comma()) &
  theme(legend.spacing.y = unit(0.4, "line"))
#ggsave("trees/Organon/MAI measured+Organon SWO predicted - preliminary.png", height = 13, width = 20, units = "cm", dpi = 150)


# stand trajectories
logBreaks = c(1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000, 4000)
logMinorBreaks = c(4, 6, 7, 8, 9, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 2000, 3000)
sdi = crossing(tph = c(1, 1000, 4000), sdi = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000)) %>% mutate(qmd = 25.4 * (sdi / tph)^(1/1.605))

ggplot() +
  geom_path(aes(x = tph, y = qmd, group = sdi), sdi, color = "grey70", linetype = "longdash") +
  geom_path(aes(x = TPH, y = QMD, group = stand, linewidth = 2.47 * area), standTrajectories, arrow = arrow(length = unit(0.25, "line"), type = "closed")) +
  geom_path(aes(x = TPH, y = QMD, color = standingCmh, group = stand, linewidth = 2.47 * area), standTrajectories %>% group_by(stand) %>% slice_max(standAge, n = -1) %>% slice_min(standAge, n = 16)) +
  geom_label(aes(x = 4, y = 200, label = "missing stands"), color = "red", fill = alpha("white", 0.7), label.padding = unit(0.15, "line"), label.size = NA, hjust = 0, size = 3.0, vjust = -0.2) +
  geom_label(aes(x = tph, y = qmd, label = sdi), sdi %>% filter(tph == 4000, sdi %in% c(1000, 2000)), color = "grey70", fill = alpha("white", 0.7), label.padding = unit(0.15, "line"), label.size = NA, size = 3.0) +
  labs(x = "stand density, mean trees per hectare", y = "stand QMD, cm", color = "live stem\nvolume,\nm³ ha ¹", linewidth = "stand\narea, ac", title = "a) Organon SWO: 738 stands, 2016 ground data") +
  coord_cartesian(xlim = c(4, 4000), ylim = c(10, 200)) +
  scale_color_viridis_c(breaks = seq(0, 4500, by = 1500), limits = c(0, 4500)) +
  scale_linewidth_continuous(range = c(0.1, 3.0)) +
  scale_x_log10(breaks = logBreaks, minor_breaks = logMinorBreaks) +
  scale_y_log10(breaks = logBreaks, minor_breaks = logMinorBreaks) +
  theme(legend.position = "none", legend.spacing.y = unit(0.4, "line"))
#ggsave("Presentation/stand trajectories 2023-08-07 Organon SWO.png", units = "cm", width = 9.1, height = 9, dpi = 300)