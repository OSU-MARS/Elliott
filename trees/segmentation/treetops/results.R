# assumes library()s, functions, and dataset from treetops.R setup
pointsOfInterest46 = bind_rows(st_read(file.path(acceptedTreetopsPath, "s04200w06840.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(acceptedTreetopsPath, "s04200w06810.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(acceptedTreetopsPath, "s04230w06810.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04230w06810"))
#st_drop_geometry(pointsOfInterest46) %>% filter(str_starts(notes, "broken top") | (notes %in% c("broken top", "point cloud ambiguous", "reiterated leader", "reiterated leader obscured by branch", "snag lacking observable top", "top obscured by branch", "top obscured by noise", "tree lacking observable top")))
#print(st_drop_geometry(pointsOfInterest46) %>% group_by(notes) %>% summarize(n = n()), n = 100)
#pointsOfInterest46 %>% filter(notes == "leanging hardwood snag")

chm46 = rast(file.path(dsmPath, "chm.vrt"))
missingOrAmbiguous46 = pointsOfInterest46 %>% filter(notes %in% c("point cloud ambiguous", "broken top obscured by branch", "reiterated leader obscured by branch", "snag lacking observable top", "top obscured by branch", "top obscured by noise", "top obscured by snag", "tree lacking observable top"))
missingOrAmbiguous46$height = 0.3048 * terra::extract(chm46, st_coordinates(missingOrAmbiguous46))[, 1] # convert to metric
missingOrAmbiguous46 = st_drop_geometry(missingOrAmbiguous46) %>% mutate(notes = factor(forcats::fct_collapse(factor(notes), `obscured by branch` = c("broken top obscured by branch", "reiterated leader obscured by branch", "top obscured by branch"), `obscured by noise` = c("top obscured by noise"), `local maxima absent` = c("snag lacking observable top", "top obscured by snag", "tree lacking observable top")), levels = c("obscured by branch", "obscured by noise", "local maxima absent", "point cloud ambiguous")))
missingOrAmbiguous46 %>% group_by(notes) %>% summarize(n = n())

# Figure TBD: dataset distribution by height
totalTileAreaHa = length(unique(treetopData$tile)) * (0.3048 * treetopOptions$tileSize)^2 / 10000

localMaximaHistogram = treetopData %>% mutate(heightClass = 0.5 * round(height / 0.5)) %>% group_by(heightClass, treetop) %>%
  summarize(maxima = n(), .groups = "drop_last") %>%
  mutate(maximaInHeightClass = sum(maxima))
missingOrAmbiguousHistogram = missingOrAmbiguous46 %>% mutate(heightClass = 0.5 * round(height / 0.5)) %>% group_by(heightClass, notes) %>%
  summarize(trees = n(), .groups = "drop_last")
treetopProbability = treetopData %>% mutate(heightClass = 0.5 * round(height / 0.5)) %>%
  group_by(heightClass, radius) %>%
  summarize(probability = sum(((treetop == "yes") | (treetop == "merge")) & as.logical(isTreetopRadiusDsm)) / n(), .groups = "drop") 

ggplot() +
  geom_col(aes(x = maxima / totalTileAreaHa, y = heightClass, alpha = heightClass >= 5, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 0.5) + # geom_col() does not stack reliably if width exceeds height class size
  geom_segment(aes(x = 0, y = 4.5, xend = 20, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  labs(x = bquote("local maxima ha"^-1), y = "height above ground, m",alpha = NULL, fill = "local maxima\nclassification", title = paste(plotLetters[1], "local maxima density")) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 90)) +
  scale_fill_manual(breaks = c("yes", "merge", "noise", "maybe noise", "no"), labels = c("single point treetop", "merge point", "residual noise", "processing artifact", "other"), values = c("purple", "green3", "red", "dodgerblue3", "grey90")) +
  scale_x_continuous(labels = scales::comma) +
ggplot() +
  geom_col(aes(x = maxima / maximaInHeightClass, y = heightClass, alpha = heightClass >= 5, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 0.5) +
  geom_segment(aes(x = 0, y = 4.5, xend = 1, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 90)) +
  labs(x = "probability", y = NULL, alpha = NULL, fill = "local maxima\nclassification", title = paste(plotLetters[2], "class probabilities")) +
  scale_fill_manual(breaks = c("yes", "merge", "noise", "maybe noise", "no"), labels = c("single point treetop", "merge point", "residual noise", "processing artifact", "other"), values = c("purple", "green3", "red", "dodgerblue3", "grey90")) +
  scale_x_continuous(labels = scales::percent) +
  theme(axis.title.x = element_text(vjust = 0.3)) +
ggplot() +
  geom_col(aes(x = trees / totalTileAreaHa, y = heightClass, alpha = heightClass >= 5, fill = notes, group = heightClass), missingOrAmbiguousHistogram, orientation = "y", width = 0.5) +
  geom_segment(aes(x = 0, y = 4.5, xend = 1.5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 90)) +
  labs(x = bquote("trees ha"^-1), y = NULL, alpha = NULL, fill = "missing and uncertain\ntrees", title = paste(plotLetters[3], "omissions")) +
  scale_fill_manual(breaks = c("obscured by branch", "obscured by noise", "local maxima absent", "point cloud ambiguous"), values = c("forestgreen", "red", "cyan", "grey80")) +
ggplot() +
  geom_tile(aes(x = radius, y = heightClass - 0.25, alpha = heightClass >= 5, fill = probability), treetopProbability) + 
  geom_line(aes(x = radius, y = height, color = "DSM"), tibble(height = seq(0, 90), radius = 0.42184646 + 0.03624721 * height^0.99839497)) +
  geom_segment(aes(x = 0, y = 4.5, xend = 5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_fixed(ratio = 0.5, xlim = c(0, NA), ylim = c(0, 90)) +
  labs(x = "dominance radius, m", y = NULL, alpha = NULL, color = "decision boundary", fill = "treetop probability", title = paste("   ", plotLetters[4], "treetop probabilities")) + # alignment hack for more consistent tile positioning
  scale_color_manual(breaks = c("DSM"), labels = c(bquote(f(h) == a[0] + a[1]*h^b[1])), values = c("cyan")) +
  scale_fill_viridis_c() + 
  theme(axis.title.x = element_text(vjust = -0.3)) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, widths = c(1, 0.7, 0.5, 0.6), guides = "collect") &
  guides(alpha = "none") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1))

# Figure TBD: accuracy distribution by height
radiusDsmAccuracyPower = readRDS("trees/segmentation/treetops/radius DSM power s4268 458k 2x25.Rds")
radiusDsmAccuracyPowerByHeight = unnest_cv_accuracy_by_height(radiusDsmAccuracyPower)

randomForestAccuracy = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde 2x25 m9n3.Rds")
randomForestAccuracyByHeight = unnest_cv_accuracy_by_height(randomForestAccuracy)

aucByHeight = left_join(radiusDsmAccuracyPowerByHeight %>% rename(nRadius = n, radiusAccuracy = overallAccuracy) %>% select(-meanN),
                        randomForestAccuracyByHeight %>% rename(nRandomForest = n, randomForestAccuracy = treetopAccuracy) %>% select(-meanN),
                        by = join_by(repetition, fold, heightClass)) %>%
  filter(is.na(nRadius) == FALSE, is.na(nRandomForest) == FALSE) %>% # WeightedROC() does not support incomplete cases
  group_by(heightClass) %>%
  summarize(auc = WeightedAUC(WeightedROC(guess = c(radiusAccuracy, randomForestAccuracy), 
                                          label = factor(c(rep(0, n()), rep(1, n())), levels = c(0, 1)), # lowerIsBetter = factor(c(rep(0, n()), rep(1, n())), levels = c(0, 1)) -> probability random forest is more accurate
                                          weight = c(nRadius, nRandomForest))),
            n = 0.5 * (mean(nRadius) + mean(nRandomForest)))

accuracyDeltaByHeight = left_join(radiusDsmAccuracyPowerByHeight %>% rename(nRadius = n, radiusAccuracy = overallAccuracy) %>% select(-meanN),
                                  randomForestAccuracyByHeight %>% rename(nRandomForest = n, randomForestAccuracy = treetopAccuracy) %>% select(-meanN),
                                  by = join_by(repetition, fold, heightClass)) %>%
  group_by(heightClass) %>%
  summarize(deltaAccuracyMedian = median(randomForestAccuracy, na.rm = TRUE) - median(radiusAccuracy, na.rm = TRUE))

ggplot() +
  geom_violin(aes(x = overallAccuracy, y = heightClass, color = after_stat(count), group = heightClass, weight = meanN), radiusDsmAccuracyPowerByHeight, draw_quantiles = c(0.5), linewidth = 0.3, width = 1.75) +
  geom_segment(aes(x = 0, y = 4.5, xend = 5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "treetop detection\naccuracy", y = "height above ground, m", color = "local\nmaxima", title = bquote(.(plotLetters[1])~"radius, "*f(h) == a[0] + a[1]*h^b[1])) + # paste(plotLetters[1], "f(h) = power function")
  scale_x_continuous(labels = scales::percent) +
  theme(plot.title = element_text(vjust = 0.66)) +
ggplot() +
  geom_violin(aes(x = treetopAccuracy, y = heightClass, color = after_stat(count), group = heightClass, weight = meanN), randomForestAccuracyByHeight, draw_quantiles = c(0.5), linewidth = 0.3, width = 1.75) +
  geom_segment(aes(x = 0, y = 4.5, xend = 5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  guides(color = "none") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "treetop detection\naccuracy", y = NULL, color = "local\nmaxima", title = paste(plotLetters[2], "clustered random forest")) +
  scale_x_continuous(labels = scales::percent) +
ggplot() +
  geom_col(aes(x = auc, y = heightClass, alpha = heightClass >= 5, fill = n), aucByHeight, orientation = "y") +
  geom_segment(aes(x = 0, y = 4.5, xend = 1, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  guides(alpha = "none", fill = "none") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "improvement\nprobability", y = NULL, alpha = NULL, fill = "local\nmaxima", title = paste(plotLetters[3], "learning success")) +
  scale_fill_gradient(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c(1, 10, 100, 1000, 10000, 100000), limits = c(1, NA), high = "#132B43", low = "#96F1FF", transform = "log10") +
  scale_x_continuous(labels = scales::percent) +
ggplot() +
  geom_col(aes(x = deltaAccuracyMedian, y = heightClass, alpha = heightClass >= 5, fill = 100 * deltaAccuracyMedian), accuracyDeltaByHeight, orientation = "y") +
  geom_segment(aes(x = -1, y = 4.5, xend = 1, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  labs(x = "median\naccuracy increase", y = NULL, alpha = NULL, fill = "learning\ngain", title = paste(plotLetters[4], "learning gain")) +
  coord_cartesian(xlim = c(-1, 1)) +
  guides(alpha = "none", fill = guide_colorbar()) +
  scale_fill_gradientn(colors = c("#FF0000", "#FF0000", "#FF0000", "grey80", "#009000", "#009000", "#009000"), values = c(-100, -10, -2, 0, 2, 10, 100), breaks = c(-100, -10, -1, 0, 1, 10, 100), labels = scales::percent(c(-100, -10, -1, 0, 1, 10, 100), scale = 1), limits = c(-100, 100), rescaler = scales::rescale_none, transform = scales::pseudo_log_trans(sigma = 0.25)) +
  scale_x_continuous(breaks = c(-1, -0.1, 0, 0.1, 1), labels = scales::percent, minor_breaks = c(-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), transform = scales::pseudo_log_trans(sigma = 0.02, base = 10)) +
plot_annotation(theme = theme()) +
plot_layout(nrow = 1, widths = c(1, 1, 0.7, 0.85), guides = "collect") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_color_gradient(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c(1, 10, 100, 1000, 10000, 100000), limits = c(1, NA), high = "#132B43", low = "#96F1FF", transform = "log10") &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1))

# Figure TBD: random forest variable selection and importance

# Figure TBD: forest-wide classifications by tile

# Figure TBD: density and top height in 2015-16 cruised stands by type


# Figure TBD: runtimes
# Get-Dsm: 2,244,000,000 DSM cells from 561 tiles (49364.3 Mpoints) in 06:01: 1739.64 GB at 1.55 tiles/s (88.0 Mpoints/tile, 4.8 GB/s).
# Get-LocalMaxima: Found 99,598,584 DSM, 20,299,060 CMM, and 96,860,358 CHM maxima within 561 tiles in 19:26 (177,537 DSM maxima/tile).
# DSM optim 14.3329 secs + 1.435918 m 2x25 cross validation 
# clustered random forest: Boruta 36.22h + VSURF + 39.21h + tune 8.823h + cross validation 1.528h + 4.409m fit

# Table TBD: dataset content
local46maxima = bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810"),
                          get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840"),
                          get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810"))
#local46maxima %>% group_by(tile) %>% summarize(`single top` = sum(treetop == "yes"), `merge point` = sum(treetop == "merge"), `residual noise` = sum(treetop == "noise"), other = sum(treetop == "no"))
acceptedTreetops46 = bind_rows(st_drop_geometry(st_read(file.path(acceptedTreetopsPath, "s04200w06840.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06840"),
                               st_drop_geometry(st_read(file.path(acceptedTreetopsPath, "s04200w06810.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06810"), # has vertical CRS which causes bind_rows() to fail on CRS mismatch
                               st_drop_geometry(st_read(file.path(acceptedTreetopsPath, "s04230w06810.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04230w06810"))

left_join(left_join(local46maxima %>% group_by(tile, treetop) %>% summarize(maxima = n(), .groups = "drop") %>% 
                    pivot_wider(id_cols = "tile", names_from = "treetop", values_from = "maxima") %>%
                    mutate(`total maxima` = no + yes + merge + noise + `maybe noise`,
                           `processing artifact` = c(0, 2, 4), # s04200w06810, s04230w06810, s04200w06840 
                           no = no + `maybe noise` - `processing artifact`),
                  acceptedTreetops46 %>% group_by(tile) %>% summarize(treetops = n()),
                  by = join_by(tile)),
          missingOrAmbiguous46 %>% group_by(tile, notes) %>% summarize(trees = n(), .groups = "drop") %>%
            pivot_wider(id_cols = "tile", names_from = "notes", values_from = "trees") %>%
            mutate(`obscured by noise` = replace_na(`obscured by noise`, 0)),
          by = join_by(tile)) %>%
  bind_rows(summarize(., across(where(is.numeric), sum))) %>%
  relocate(tile, yes, merge, noise, `processing artifact`, no, `total maxima`, treetops, `maybe noise`) %>%
  rename(`single top` = yes, `merge point` = merge, `residual noise` = noise, other = no) %>%
  mutate(tile = replace_na(tile, "total"))
