library(arrow)
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

as_harvest_system = function(harvestSystem)
{
  return(factor(harvestSystem, levels = seq(0, 10), labels = c("None", "FallersGrappleSwingYarderProcessorLoader", "FallersGrappleYoaderProcessorLoader", "FellerBuncherGrappleSwingYarderProcessorLoader", "FellerBuncherGrappleYoaderProcessorLoader", "TrackedHarvesterForwarder", "TrackedHarvesterGrappleSwingYarderLoader", "TrackedHarvesterGrappleYoaderLoader", "WheeledHarvesterForwarder", "WheeledHarvesterGrappleSwingYarderLoader", "WheeledHarvesterGrappleYoaderLoader")))
}


## no management option
standTrajectoriesFile ="trees/Organon/Elliott stand trajectories 2016-2116.feather"
standTrajectories = left_join(read_feather(standTrajectoriesFile, mmap = FALSE) %>%
                                mutate(mai = standingMbfh / standAge), # MBF/ha
                              read_xlsx("trees/Organon/Elliott Organon cruise records 2015-16.xlsx", sheet = "stands") %>% rename(stand = id, cruiseAge2016 = age),
                              by = c("stand"))

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


# net annualized revenue and components, undiscounted
ggplot() +
  geom_bin_2d(aes(x = standAge, y = (regenPond2S + regenPond3S + regenPond4S - regenFallerGrappleSwingYarderCost) / standAge, weight = area), standTrajectories %>% filter(cruiseAge2016 < (2016 - 1950)), binwidth = c(2, 20)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 2000)) +
  labs(x = "stand age, years", y = bquote("net annualized revenue, US$ ha"^-1~"year"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = regenPond2S + regenPond3S + regenPond4S, weight = area), standTrajectories %>% filter(cruiseAge2016 < (2016 - 1950)), binwidth = c(2, 1000)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 325000)) +
  labs(x = "stand age, years", y = bquote("pond value, US$ ha"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = regenFallerGrappleSwingYarderCost, weight = area), standTrajectories %>% filter(cruiseAge2016 < (2016 - 1950)), binwidth = c(2, 1000)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 325000)) +
  labs(x = "stand age, years", y = bquote("harvest cost, US$ ha"^-1), fill = "area, ha") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_fill_viridis_c(breaks = c(1, 10, 100, 500), limits = c(1, 500), trans = "log10") &
  scale_y_continuous(labels = scales::label_comma()) &
  theme(legend.spacing.y = unit(0.4, "line"))
              
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



## intensive stands
standTrajectoriesFile ="trees/Organon/Elliott intensive prescriptions max LEV.feather"
#standTrajectoriesFile ="trees/Organon/Elliott intensive prescriptions max NPV.feather"
intensiveTrajectories = left_join(read_feather(standTrajectoriesFile, mmap = FALSE) %>%
                                    mutate(mai = standingMbfh / standAge), # MBF/ha
                                  read_xlsx("trees/Organon/Elliott Organon cruise records 2015-16.xlsx", sheet = "stands") %>% rename(stand = id, cruiseAge2016 = age),
                                  by = c("stand")) %>%
  mutate(thinMinCostSystem = as_harvest_system(thinMinCostSystem),
         regenMinCostSystem = as_harvest_system(regenMinCostSystem),
         thinCost = case_match(thinMinCostSystem, "None" ~ 0,
                                                  "FallersGrappleSwingYarderProcessorLoader" ~ thinFallerGrappleSwingYarderCost,
                                                  "FallersGrappleYoaderProcessorLoader" ~ thinFallerGrappleYoaderCost,
                                                  "FellerBuncherGrappleSwingYarderProcessorLoader" ~ thinFellerBuncherGrappleSwingYarderCost,
                                                  "FellerBuncherGrappleYoaderProcessorLoader" ~ thinChainsawCmhWithFellerBuncherAndGrappleYoader,
                                                  "TrackedHarvesterForwarder" ~ thinTrackedHarvesterForwarderCost,
                                                  "TrackedHarvesterGrappleSwingYarderLoader" ~ thinTrackedHarvesterGrappleSwingYarderCost,
                                                  "TrackedHarvesterGrappleYoaderLoader" ~ thinTrackedHarvesterGrappleYoaderCost,
                                                  "WheeledHarvesterForwarder" ~ thinWheeledHarvesterForwarderCost, 
                                                  "WheeledHarvesterGrappleSwingYarderLoader" ~ thinWheeledHarvesterGrappleSwingYarderCost,
                                                  "WheeledHarvesterGrappleYoaderLoader" ~ thinWheeledHarvesterGrappleYoaderCost),
         regenHarvestCost = case_match(regenMinCostSystem, "None" ~ 0,
                                       "FallersGrappleSwingYarderProcessorLoader" ~ regenFallerGrappleSwingYarderCost,
                                       "FallersGrappleYoaderProcessorLoader" ~ regenFallerGrappleYoaderCost,
                                       "FellerBuncherGrappleSwingYarderProcessorLoader" ~ regenFellerBuncherGrappleSwingYarderCost,
                                       "FellerBuncherGrappleYoaderProcessorLoader" ~ regenChainsawCmhWithFellerBuncherAndGrappleYoader,
                                       "TrackedHarvesterForwarder" ~ NA,
                                       "TrackedHarvesterGrappleSwingYarderLoader" ~ regenTrackedHarvesterGrappleSwingYarderCost,
                                       "TrackedHarvesterGrappleYoaderLoader" ~ regenTrackedHarvesterGrappleYoaderCost,
                                       "WheeledHarvesterForwarder" ~ NA, 
                                       "WheeledHarvesterGrappleSwingYarderLoader" ~ regenWheeledHarvesterGrappleSwingYarderCost,
                                       "WheeledHarvesterGrappleYoaderLoader" ~ regenWheeledHarvesterGrappleYoaderCost),
         thinNetRevenue = if_else(standAge == thin1, replace_na(thinPond2S, 0) + replace_na(thinPond3S, 0) + replace_na(thinPond4S, 0) - thinCost, 0),
         regenHarvestNetRevenue = if_else(standAge == rotation, replace_na(regenPond2S, 0) + replace_na(regenPond3S, 0) + replace_na(regenPond4S, 0) - regenHarvestCost, 0)) %>%
  group_by(stand, thin1, thin2, thin3, rotation) %>%
  mutate(netRevenue = sum(thinNetRevenue + regenHarvestNetRevenue),
         npvRotation = sum(if_else(standAge == rotation, NPV, 0)),
         levRotation = sum(if_else(standAge == rotation, LEV, 0))) %>%
  group_by(stand) %>%
  mutate(maxNetRevenue = max(netRevenue),
         maxNpv = max(npvRotation),
         maxLev = max(levRotation)) %>%
  ungroup()


# net harvest revenues by stand age
# TODO: filter out pre-2021 harvest
ggplot() +
  geom_point(aes(x = standAge, y = thinNetRevenue, color = "thin", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(netRevenue == maxNetRevenue, thinNetRevenue > 0), alpha = 0.3, shape = 16) +
  geom_point(aes(x = standAge, y = regenHarvestNetRevenue, color = "regeneration harvest", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(netRevenue == maxNetRevenue, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 80), ylim = c(-750 * 150/15, 150000)) +
  labs(x = "stand age, years", y = bquote("net harvest revenue, US$ ha"^-1), color = NULL, size = "stand area, ha") +
ggplot() +
  geom_point(aes(x = standAge, y = NPV, color = "thin", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(npvRotation == maxNpv, thinNetRevenue > 0), alpha = 0.3, shape = 16) +
  geom_point(aes(x = standAge, y = NPV, color = "regeneration harvest", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(npvRotation == maxNpv, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 80), ylim = c(-750 * 40/15, 40000)) +
  labs(x = "stand age, years", y = bquote("net present value, US$ ha"^-1), color = NULL, size = "stand area, ha") +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = LEV, size = area, color = "infinite rotations"), intensiveTrajectories %>% group_by(stand) %>% filter(levRotation == maxLev, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 80), ylim = c(-750, 15000)) +
  labs(x = "rotation age, years", y = bquote("land expectation value, US$ ha"^-1), color = NULL, size = "stand area, ha") +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(netRevenue == maxNetRevenue, thinNetRevenue > 0), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = "thinning\nintensity, %BA", color = NULL, size = "stand area, ha") +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(npvRotation == maxNpv, thinNetRevenue > 0), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = NULL, color = NULL, size = "stand area, ha") +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), intensiveTrajectories %>% group_by(stand) %>% filter(levRotation == maxLev, thinNetRevenue > 0), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = NULL, color = NULL, size = "stand area, ha") +
  theme(legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 3, guides = "collect", heights = c(1, 0.15)) &
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2, override.aes = list(alpha = 0.5))) &
  scale_color_manual(breaks = c("thin", "regeneration harvest", "infinite rotations"), limits = c("thin", "regeneration harvest", "infinite rotations"), values = c("darkorchid", "firebrick", "grey10")) &
  scale_size_area(limits = c(0, 50), max_size = 4) &
  scale_y_continuous(labels = scales::label_comma()) &
  theme(legend.spacing.y = unit(0.4, "line"))
#ggsave("Presentation/intensive management 2023-08-30 Organon SWO max LEV.png", units = "cm", width = 22, height = 12, dpi = 150)
#ggsave("Presentation/intensive management 2023-08-30 Organon SWO max NPV.png", units = "cm", width = 22, height = 12, dpi = 150)

# point checks and exploration
standSummary = intensiveTrajectories %>% filter(stand %in% c(552, 2359, 2414, 2445, 2449, 2461, 2509)) %>% group_by(stand, thin1, thin2, thin3, rotation) %>%
  summarize(netRevenue = netRevenue[1], npvRotation = npvRotation[1], levRotation = levRotation[1], .groups = "drop")
ggplot(standSummary) +
  geom_line(aes(x = rotation, y = netRevenue / rotation, color = thin1, group = paste(stand, thin1))) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "rotation length, years", y = bquote("net revenue, US$ ha"^-1~"year"^-1), color = "thin age,\nyears") +
ggplot(standSummary) +
  geom_line(aes(x = rotation, y = npvRotation, color = thin1, group = paste(stand, thin1))) +
  coord_cartesian(ylim = c(0, 6500)) +
  labs(x = "rotation length, years", y = bquote("net present value, US$ ha"^-1), color = "thin age,\nyears") +
ggplot(standSummary) +
  geom_line(aes(x = rotation, y = levRotation, color = thin1, group = paste(stand, thin1))) +
  coord_cartesian(ylim = c(0, 6500)) +
  labs(x = "rotation length, years", y = bquote("land expectation value, US$ ha"^-1), color = "thin age,\nyears") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_color_viridis_c(limits = c(34, 56)) &
  scale_y_continuous(labels = scales::label_comma())

ggplot(standSummary) +
  geom_line(aes(x = rotation, y = levRotation, color = thin1, group = paste(stand, thin1))) +
  coord_cartesian(ylim = c(0, 6500)) +
  facet_wrap(vars(stand)) +
  labs(x = "rotation length, years", y = bquote("net present value, US$ ha"^-1), color = "thin age,\nyears") +
  scale_color_viridis_c(limits = c(34, 56)) +
  scale_y_continuous(labels = scales::label_comma())
