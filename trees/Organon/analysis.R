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

read_stand_trajectories = function(trajectoryFilePath)
{
  return(read_feather(trajectoryFilePath, mmap = FALSE) %>%
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
    group_by(stand, thin1, thin2, thin3, rotation, financialScenario) %>%
    mutate(netRevenue = sum(thinNetRevenue + regenHarvestNetRevenue),
           npvRotation = sum(if_else(standAge == rotation, NPV, 0)),
           levRotation = sum(if_else(standAge == rotation, LEV, 0))) %>%
    group_by(stand, financialScenario) %>%
    mutate(maxNetRevenue = max(netRevenue),
           maxNpv = max(npvRotation),
           maxLev = max(levRotation)) %>%
    ungroup())
}


## unthinned trajectories (no management or clearcut only)
standTrajectoriesFile ="trees/Organon/Elliott stand trajectories 2016-2116.feather"
standTrajectories = left_join(read_feather(standTrajectoriesFile, mmap = FALSE) %>%
                                mutate(mai = standingMbfh / standAge), # MBF/ha
                              read_xlsx("trees/Organon/Elliott Organon cruise records 2015-16.xlsx", sheet = "stands") %>% rename(stand = id, standAge2016 = age),
                              by = c("stand")) %>%
  mutate(discountRatePct = if_else(financialScenario < 2, 0.5 * financialScenario, financialScenario - 1),
         isPlantation = standAge2016 < (2016 - 1950), 
         # currently considers only hand falling, ignoring steep slope harvest systems with tethered mechanized falling (feller buncher, harvester)
         regenNetRevenue = regenPond2S + regenPond3S + regenPond4S - pmin(regenFallerGrappleSwingYarderCost, regenFallerGrappleYoaderCost) - regenTaskCost)

optimumRotations = standTrajectories %>% filter(isPlantation, discountRatePct > 0) %>% group_by(stand, discountRatePct) %>%
  slice_max(LEV) %>%
  summarize(standAge2016 = standAge2016, rotationLength = standAge, .groups = "drop")
print(optimumRotations, n = 100)

print(standTrajectories %>% filter(stand == 16, discountRatePct > 0) %>% select(stand, discountRatePct, year, standAge, LEV), n = 140)

# MAI
ggplot() +
  geom_segment(aes(x = c(40, 65), y = 0, xend = c(40, 65), yend = 20000), color = "grey50", linetype = "longdash", linewidth = 0.4) +
  geom_bin_2d(aes(x = standAge, y = 1000 * mai, weight = area), standTrajectories %>% filter(discountRatePct == 0) %>% group_by(stand) %>% slice_min(standAge, n = 1), binwidth = c(2, 50)) +
  labs(x = "stand age, years", y = bquote("cruised MAI in 2016, BF ha"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = 1000 * mai, weight = area), standTrajectories %>% filter(discountRatePct == 0) %>% group_by(stand) %>% slice_max(standAge, n = -1), binwidth = c(2, 50)) +
  geom_segment(aes(x = c(40, 65), y = 0, xend = c(40, 65), yend = 20000), color = "grey50", linetype = "longdash", linewidth = 0.4) +
  labs(x = "stand age, years", y = bquote("Organon SWO predicted MAI, BF ha"^-1), fill = "area, ha") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(widths = c(200, 300), guides = "collect") &
  coord_cartesian(ylim = c(0, 5000)) &
  scale_fill_viridis_c(breaks = c(1, 10, 100, 500), limits = c(1, 500), trans = "log10") &
  scale_y_continuous(labels = scales::label_comma()) &
  theme(legend.spacing.y = unit(0.4, "line"))
#ggsave("trees/Organon/MAI measured+Organon SWO predicted - preliminary.png", height = 13, width = 20, units = "cm", dpi = 150)


# net annualized harvest revenue and components from intensive silviculture, undiscounted
ggplot() +
  geom_bin_2d(aes(x = standAge, y = regenNetRevenue / standAge, weight = area), standTrajectories %>% filter(isPlantation, discountRatePct == 0), binwidth = c(2, 20)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 2000)) +
  labs(x = "stand age, years", y = bquote("net annualized revenue, US$ ha"^-1~"year"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = regenPond2S + regenPond3S + regenPond4S, weight = area), standTrajectories %>% filter(isPlantation, discountRatePct == 0), binwidth = c(2, 1000)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 325000)) +
  labs(x = "stand age, years", y = bquote("pond value, US$ ha"^-1), fill = "area, ha") +
ggplot() +
  geom_bin_2d(aes(x = standAge, y = pmin(regenFallerGrappleSwingYarderCost, regenFallerGrappleYoaderCost) + regenTaskCost, weight = area), standTrajectories %>% filter(isPlantation, discountRatePct == 0), binwidth = c(2, 1000)) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 325000)) +
  labs(x = "stand age, years", y = bquote("harvest and reforestation cost, US$ ha"^-1), fill = "area, ha") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  scale_fill_viridis_c(breaks = c(1, 10, 100, 500), limits = c(1, 500), trans = "log10") &
  scale_y_continuous(labels = scales::label_comma()) &
  theme(legend.spacing.y = unit(0.4, "line"))

ggplot() +
  geom_line(aes(x = standAge, y = regenNetRevenue / standAge, group = stand), standTrajectories %>% filter(discountRatePct == 0, isPlantation), alpha = 0.1) +
  coord_cartesian(xlim = c(0, 175), ylim = c(0, 2000)) +
  labs(x = "stand age, years", y = bquote("net annualized revenue, US$ ha"^-1~"year"^-1))

ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 0, isPlantation), alpha = 0.1) +
  labs(x = NULL, y = bquote("LEV, US$ ha"^-1~"year"^-1), title = "(a) 0.0% discount rate") +
ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 0.5, isPlantation), alpha = 0.1) +
  labs(x = NULL, y = NULL, title = "(b) 0.5% discount rate") +
ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 1.0, isPlantation), alpha = 0.1) +
  labs(x = NULL, y = NULL, title = "(c) 1.0% discount rate") +
ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 2.0, isPlantation), alpha = 0.1) +
  labs(x = "stand age, years", y = bquote("LEV, US$ ha"^-1~"year"^-1), title = "(d) 2.0% discount rate") +
ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 3.0, isPlantation), alpha = 0.1) +
  labs(x = "stand age, years", y = NULL, title = "(e) 3.0% discount rate") +
ggplot() +
  geom_line(aes(x = standAge, y = LEV, group = stand), standTrajectories %>% filter(discountRatePct == 4.0, isPlantation), alpha = 0.1) +
  labs(x = "stand age, years", y = NULL, title = "(f) 4.0% discount rate") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout() &
  coord_cartesian(xlim = c(0, 175), ylim = c(0, NA))

# stand trajectories
logBreaks = c(1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000, 4000)
logMinorBreaks = c(4, 6, 7, 8, 9, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900, 2000, 3000)
sdi = crossing(tph = c(1, 1000, 4000), sdi = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000)) %>% mutate(qmd = 25.4 * (sdi / tph)^(1/1.605))

ggplot() +
  geom_path(aes(x = tph, y = qmd, group = sdi), sdi, color = "grey70", linetype = "longdash") +
  geom_path(aes(x = TPH, y = QMD, group = stand, linewidth = 2.47 * area), standTrajectories %>% filter(discountRatePct == 0), arrow = arrow(length = unit(0.25, "line"), type = "closed")) +
  geom_path(aes(x = TPH, y = QMD, color = standingCmh, group = stand, linewidth = 2.47 * area), standTrajectories %>% filter(discountRatePct == 0) %>% group_by(stand) %>% slice_max(standAge, n = -1) %>% slice_min(standAge, n = 16)) +
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



## plantation stands
standTrajectoriesFile ="trees/Organon/Elliott plantation prescriptions max LEV.feather"
#standTrajectoriesFile ="trees/Organon/Elliott plantation prescriptions max NPV.feather"
plantationTrajectoryReadStart = Sys.time() # ~22 s @ 4GB, 9900X
plantationTrajectories = left_join(read_stand_trajectories(standTrajectoriesFile) %>%
                                    mutate(mai = standingMbfh / standAge), # MBF/ha
                                  read_xlsx("trees/Organon/Elliott Organon cruise records 2015-16.xlsx", sheet = "stands") %>% rename(stand = id, standAge2016 = age),
                                  by = c("stand")) %>%
  mutate(discountRatePct = if_else(financialScenario < 2, 0.5 * financialScenario, financialScenario - 1),
         isPlantation = standAge2016 < (2016 - 1950))
(Sys.time() - plantationTrajectoryReadStart)

plantationTrajectories %>% filter(levRotation == maxLev, (standAge == thin1) | (standAge == rotation), discountRatePct> 0) %>% 
  group_by(discountRatePct, stand) %>% 
  summarize(area = area[1], isThinned = standAge[1] == thin1[1], .groups = "drop_last") %>%
  summarize(stands = n(), areaHa = sum(area), thinned = sum(isThinned), thinnedAreaHa = sum(isThinned * area)) %>%
  mutate(thinAreaPct = 100 * thinnedAreaHa / areaHa)

# LEV by discount rate
ggplot() +
  geom_point(aes(x = standAge, y = LEV, size = area, color = if_else(thin1 != -1, "with thinning", "clearcut only")), plantationTrajectories %>% filter(discountRatePct == 1.0, levRotation == maxLev, standAge == rotation), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, NA)) +
  labs(x = "rotation age, years", y = bquote("land expectation value, US$ ha"^-1), color = "preferred rotation", size = "stand area, ha", title = "(a) 1% discount rate") +
  scale_y_continuous(labels = scales::label_comma()) +
ggplot() +
  geom_point(aes(x = standAge, y = LEV, size = area, color = if_else(thin1 != -1, "with thinning", "clearcut only")), plantationTrajectories %>% filter(discountRatePct == 3.0, levRotation == maxLev, standAge == rotation), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, NA)) +
  labs(x = "rotation age, years", y = NULL, color = "preferred rotation", size = "stand area, ha", title = "(b) 3% discount rate") +
  scale_y_continuous(labels = scales::label_comma()) +
ggplot() +
  geom_point(aes(x = standAge, y = LEV, size = area, color = if_else(thin1 != -1, "with thinning", "clearcut only")), plantationTrajectories %>% filter(discountRatePct == 5.0, levRotation == maxLev, standAge == rotation), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, NA)) +
  labs(x = "rotation age, years", y = NULL, color = "preferred rotation", size = "stand area, ha", title = "(c) 5% discount rate") +
  scale_y_continuous(labels = scales::label_comma()) +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "with thinning", size = area), plantationTrajectories %>% filter(discountRatePct == 1.0, netRevenue == maxNetRevenue, standAge == thin1), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 40)) +
  labs(x = "thin age, years", y = "thinning\nintensity, % BA", color = "preferred rotation", size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "with thinning", size = area), plantationTrajectories %>% filter(discountRatePct == 3.0, netRevenue == maxNetRevenue, standAge == thin1), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 40)) +
  labs(x = "thin age, years", y = NULL, color = "preferred rotation", size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "with thinning", size = area), plantationTrajectories %>% filter(discountRatePct == 5.0, netRevenue == maxNetRevenue, standAge == thin1), alpha = 0.2, shape = 16) +
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 40)) +
  labs(x = "thin age, years", y = NULL, color = "preferred rotation", size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 3, heights = c(1, 0.25), guides = "collect") &
  guides(color = guide_legend(order = 1, override.aes = list(alpha = 0.7)), size = guide_legend(order = 2, override.aes = list(alpha = 0.4))) &
  scale_color_manual(breaks = c("with thinning", "clearcut only"), values = c("darkorchid", "firebrick")) &
  scale_size_area(limits = c(0, 50), max_size = 4) &
  theme(legend.spacing.y = unit(0.4, "line"))
#ggsave("trees/Organon/figures/plantation management Organon SWO discount rates.png", units = "cm", width = 22, height = 12, dpi = 200)

# net harvest revenues by stand age
# TODO: filter out pre-2021 harvest
discountRate = 1.0 # %
ggplot() +
  geom_point(aes(x = standAge, y = thinNetRevenue, color = "thin", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, netRevenue == maxNetRevenue, standAge == thin1), alpha = 0.3, shape = 16) +
  geom_point(aes(x = standAge, y = regenHarvestNetRevenue, color = "clearcut", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, netRevenue == maxNetRevenue, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, NA)) +
  labs(x = "stand age, years", y = bquote("net harvest revenue, US$ ha"^-1), color = NULL, size = "stand area, ha", title = "(a) cashflow") +
  scale_y_continuous(labels = scales::label_comma()) +
ggplot() +
  geom_point(aes(x = standAge, y = NPV, color = "thin", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, npvRotation == maxNpv, standAge == thin1), alpha = 0.3, shape = 16) +
  geom_point(aes(x = standAge, y = NPV, color = "clearcut", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, npvRotation == maxNpv, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, NA)) +
  labs(x = "stand age, years", y = bquote("net present value, US$ ha"^-1), color = NULL, size = "stand area, ha", title = "(b) single rotation") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = LEV, size = area, color = if_else(thin1 != -1, "thin", "clearcut")), plantationTrajectories %>% filter(discountRatePct == discountRate, levRotation == maxLev, standAge == rotation), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, NA)) +
  labs(x = "rotation age, years", y = bquote("land expectation value, US$ ha"^-1), color = NULL, size = "stand area, ha", title = "(c) multi-rotation") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, netRevenue == maxNetRevenue, standAge == thin1), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = "thinning\nintensity, %BA", color = NULL, size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, npvRotation == maxNpv, standAge == thin1), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = NULL, color = NULL, size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
ggplot() +
  geom_point(aes(x = standAge, y = 100 * BAintensity, color = "thin", size = area), plantationTrajectories %>% filter(discountRatePct == discountRate, levRotation == maxLev, standAge == thin1), alpha = 0.3, shape = 16) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 50)) +
  labs(x = "thin age, years", y = NULL, color = NULL, size = "stand area, ha") +
  scale_y_continuous(breaks = seq(0, 50, by = 10)) +
  theme(legend.position = "none") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 2, ncol = 3, guides = "collect", heights = c(1, 0.15)) &
  guides(color = guide_legend(order = 1, override.aes = list(alpha = 0.8)), size = guide_legend(order = 2, override.aes = list(alpha = 0.3))) &
  scale_color_manual(breaks = c("thin", "clearcut"), limits = c("thin", "clearcut"), values = c("darkorchid", "firebrick")) &
  scale_size_area(limits = c(0, 50), max_size = 4) &
  theme(legend.spacing.y = unit(0.4, "line"))
#ggsave("trees/Organon/figures/plantation management Organon SWO max LEV 6.png", units = "cm", width = 22, height = 12, dpi = 150)

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

# mortality
trees2021organon = read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE)
                             
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1000, yend = 1000), color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_bin_2d(aes(x = liveExpansionFactor2016, y = liveExpansionFactor2021), trees2021organon %>% filter(year <= 2021) %>% select(stand, plot, tag, year, liveExpansionFactor) %>% pivot_wider(id_cols = c("stand", "plot", "tag"), names_from = year, names_prefix = "liveExpansionFactor", values_from = liveExpansionFactor), binwidth = c(0.1, 0.1)) +
  labs(x = "trees per hectare, 2016", y = "trees per hectare, 2021", fill = "tree records") +
  scale_fill_viridis_c(labels = scales::label_comma(), trans = "log10") +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10000), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900), trans = scales::transform_pseudo_log()) +
  scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10000), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900), trans = scales::transform_pseudo_log())
