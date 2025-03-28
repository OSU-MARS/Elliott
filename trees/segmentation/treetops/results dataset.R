# assumes library()s, functions, and treetopDataDsm from treetops.R setup
library(writexl)

get_local_maxima_histogram = function(treetopData)
{
  return(treetopData %>% mutate(heightClass = round(height)) %>% group_by(heightClass, treetop) %>%
           summarize(maxima = n(), .groups = "drop_last") %>%
           mutate(maximaInHeightClass = sum(maxima)))
}

get_treetop_probability = function(treetopData)
{
  return(treetopData %>% mutate(heightClass = round(height)) %>%
           group_by(heightClass, radius) %>%
           summarize(probability = sum(((treetop == "yes") | (treetop == "merge")) & as.logical(isTreetopRadius)) / n(), .groups = "drop"))
}

plot_local_maxima_distribution = function(localMaximaHistogram, plotLetter = plotLetters[2], plotTitle = "class distribution")
{
  return(ggplot() +
           geom_col(aes(x = maxima / maximaInHeightClass, y = heightClass, alpha = heightClass >= 5, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 1) +
           geom_segment(aes(x = 0, y = 4.5, xend = 1, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
           coord_cartesian(xlim = c(0, NA), ylim = c(0, 90)) +
           labs(x = "probability", y = NULL, alpha = NULL, fill = "local maxima\nclassification", title = paste(plotLetter, plotTitle)) +
           scale_fill_manual(breaks = c("yes", "merge", "noise", "maybe noise", "no"), labels = c("single point treetop", "merge point", "residual noise", "processing artifact", "other"), values = c("purple", "green3", "red", "dodgerblue3", "grey90")) +
           scale_x_continuous(labels = scales::percent) +
           theme(axis.title.x = element_text(vjust = 0.3)))
}

plot_local_maxima_density = function(localMaximaHistogram, plotLetter = paste("   ", plotLetters[1]), plotTitle = "DSM local maxima density", yLabel = "height above ground, m")
{
  return(ggplot() +
           geom_col(aes(x = maxima / totalTileAreaHa, y = heightClass, alpha = heightClass >= 5, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 1) + # geom_col() does not stack reliably if width exceeds height class size
           geom_segment(aes(x = 0, y = 4.5, xend = 40, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
           labs(x = bquote("local maxima ha"^-1), y = yLabel, alpha = NULL, fill = "local maxima\nclassification", title = paste(plotLetter, plotTitle)) +
           coord_cartesian(xlim = c(0, 40), ylim = c(0, 90)) +
           scale_fill_manual(breaks = c("yes", "merge", "noise", "maybe noise", "no"), labels = c("single point treetop", "merge point", "residual noise", "processing artifact", "other"), values = c("purple", "green3", "red", "dodgerblue3", "grey90")) +
           scale_x_continuous(labels = scales::comma))
}

plot_treetop_probability = function(treetopProbability, fAtH = function(h) { return(0.42644864 + 0.03226514 * h^1.02543385) }, plotLetter = paste("  ", plotLetters[4]), plotTitle = "DSM treetop distribution")
{
  return(ggplot() +
           geom_tile(aes(x = radius, y = heightClass, alpha = heightClass >= 5, fill = probability), treetopProbabilityDsm %>% filter(heightClass >= minimumHeightClass)) + # near zero probabilities @ 1 m look like no data grey at 0.5 alpha
           geom_line(aes(x = radius, y = height, color = "DSM"), tibble(height = seq(0, 90), radius = fAtH(height))) +
           geom_segment(aes(x = 0, y = 4.5, xend = 5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
           coord_fixed(ratio = 0.5, xlim = c(0, NA), ylim = c(0, 90)) +
           labs(x = "dominance radius, m", y = NULL, alpha = NULL, color = "decision boundary", fill = "treetop probability", title = paste(plotLetter, plotTitle)) + # alignment hack for more consistent tile positioning
           scale_color_manual(breaks = c("DSM"), labels = c(bquote(f(h) == a[0] + a[1]*h^b[1])), values = c("cyan")) +
           scale_fill_viridis_c() + 
           theme(axis.title.x = element_text(vjust = -0.3)))
}

unnest_binary_confusion_median = function(crossValidatedAccuracy)
{
  confusion = crossValidatedAccuracy %>% select(repetition, fold) %>% mutate(confusion = vector(mode = "list", length = n()))
  for (row in 1:nrow(crossValidatedAccuracy))
  {
    confusion$confusion[[row]] = crossValidatedAccuracy$confusionMatrix[[row]]$table
  }
  confusion %<>% unnest_wider(col = confusion)
  if (("no" %in% names(confusion)) == FALSE)
  {
    # binary matrix from radius classification
    confusion %<>% rename(no = `FALSE`, yes = `TRUE`)
  }
  confusion %<>% pivot_longer(cols = c("no", "yes"), names_to = "prediction", values_to = "reference") %>%
    mutate(reference_no = reference[, 1], reference_yes = reference[, 2]) %>%
    select(-reference) %>%
    pivot_longer(cols = c("reference_no", "reference_yes"), names_prefix = "reference_", names_to = "reference", values_to = "n") %>%
    group_by(repetition, fold) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup() %>%
    mutate(prediction = forcats::fct_recode(factor(prediction, levels = c("yes", "no")), treetop = "yes", other = "no"),
           reference = forcats::fct_recode(factor(reference, levels = c("yes", "no")), treetop = "yes", other = "no"))
  
  confusionMedian = confusion %>% group_by(prediction, reference) %>%
    summarize(fraction = median(fraction), .groups = "drop")
  return(confusionMedian)
}

unnest_cv_accuracy_by_height = function(crossValidatedAccuracy)
{
  return(unnest(crossValidatedAccuracy %>% select(repetition, fold, overallAccuracyByHeight), cols = overallAccuracyByHeight) %>% 
           mutate(meanN = n / max(repetition)))
}

unnest_quinary_confusion_median = function(crossValidatedAccuracy)
{
  confusion = crossValidatedAccuracy %>% select(repetition, fold) %>% mutate(confusion = vector(mode = "list", length = n()))
  for (row in 1:nrow(crossValidatedAccuracy))
  {
    confusion$confusion[[row]] = crossValidatedAccuracy$confusionSubmatrix[[row]]$table
  }
  confusion %<>% unnest_wider(col = confusion) %>% 
    pivot_longer(cols = c("no", "yes", "merge", "noise", "maybe noise"), names_to = "prediction", values_to = "reference") %>% 
    mutate(reference_no = reference[, 1], reference_yes = reference[, 2], reference_merge = reference[, 3], reference_noise = reference[, 4], `reference_maybe noise` = reference[, 5]) %>%
    select(-reference) %>%
    pivot_longer(cols = c("reference_no", "reference_yes", "reference_merge", "reference_noise", "reference_maybe noise"), names_prefix = "reference_", names_to = "reference", values_to = "n") %>%
    group_by(repetition, fold) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup() %>%
    mutate(prediction = forcats::fct_recode(factor(prediction, levels = c("yes", "merge", "no", "noise", "maybe noise")), treetop = "yes", `merge point` = "merge", `residual noise` = "noise", other = "no", `processing artifact` = "maybe noise"),
           reference = forcats::fct_recode(factor(reference, levels = c("yes", "merge", "no", "noise", "maybe noise")), treetop = "yes", `merge point` = "merge", `residual noise` = "noise", other = "no", `processing artifact` = "maybe noise"))
  
  confusionMedian = confusion %>% group_by(prediction, reference) %>%
    summarize(fraction = median(fraction), .groups = "drop")
  return(confusionMedian)
}

acceptedTreetops46dsm = bind_rows(st_drop_geometry(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06840"),
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06810"), # has vertical CRS which causes bind_rows() to fail on CRS mismatch
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04230w06810"))
acceptedTreetops46chm = bind_rows(st_drop_geometry(st_read(file.path(acceptedTreetopsChmPath, "s04200w06840 chm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06840"),
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsChmPath, "s04200w06810 chm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06810"), # has vertical CRS which causes bind_rows() to fail on CRS mismatch
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsChmPath, "s04230w06810 chm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04230w06810"))
acceptedTreetops46cmm = bind_rows(st_drop_geometry(st_read(file.path(acceptedTreetopsCmmPath, "s04200w06840 cmm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06840"),
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsCmmPath, "s04200w06810 cmm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04200w06810"), # has vertical CRS which causes bind_rows() to fail on CRS mismatch
                                  st_drop_geometry(st_read(file.path(acceptedTreetopsCmmPath, "s04230w06810 cmm.gpkg"), layer = "treetops", quiet = TRUE)) %>% mutate(tile = "s04230w06810"))

chm46 = rast(file.path(dsmPath, "chm.vrt"))
pointsOfInterest46 = bind_rows(st_read(file.path(acceptedTreetopsDsmPath, "s04200w06840.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04200w06840"),
                               st_read(file.path(acceptedTreetopsDsmPath, "s04200w06810.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04200w06810"),
                               st_read(file.path(acceptedTreetopsDsmPath, "s04230w06810.gpkg"), layer = "points of interest", quiet = TRUE) %>% mutate(tile = "s04230w06810"))

missingOrAmbiguous46dsm = pointsOfInterest46 %>% filter(notes %in% c("point cloud ambiguous", "broken top obscured by branch", "reiterated leader obscured by branch", "snag lacking observable top", "top obscured by branch", "top obscured by noise", "top obscured by snag", "tree lacking observable top"))
missingOrAmbiguous46dsm$height = 0.3048 * terra::extract(chm46, st_coordinates(missingOrAmbiguous46dsm))[, 1] # convert to metric
missingOrAmbiguous46dsm = st_drop_geometry(missingOrAmbiguous46dsm) %>% mutate(notes = factor(forcats::fct_collapse(factor(notes), `obscured by branch` = c("broken top obscured by branch", "reiterated leader obscured by branch", "top obscured by branch"), `obscured by noise` = c("top obscured by noise"), `local maxima absent` = c("snag lacking observable top", "top obscured by snag", "tree lacking observable top")), levels = c("obscured by branch", "obscured by noise", "local maxima absent", "point cloud ambiguous")))

radiusDsmAccuracyPower = readRDS("trees/segmentation/treetops/radius DSM power s4268 458k 2x25.Rds")
radiusChmAccuracyPower = readRDS("trees/segmentation/treetops/radius CHM power s4268 458k 2x25.Rds")
radiusCmmAccuracyPower = readRDS("trees/segmentation/treetops/radius CMM power s4268 458k 2x25.Rds")
randomForestAccuracy = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde 2x25 m9n3.Rds")

#missingOrAmbiguous46dsm %>% group_by(notes) %>% summarize(n = n())
#st_drop_geometry(pointsOfInterest46) %>% filter(str_starts(notes, "broken top") | (notes %in% c("broken top", "point cloud ambiguous", "reiterated leader", "reiterated leader obscured by branch", "snag lacking observable top", "top obscured by branch", "top obscured by noise", "tree lacking observable top")))
#print(st_drop_geometry(pointsOfInterest46) %>% group_by(notes) %>% summarize(n = n()), n = 100)
#pointsOfInterest46 %>% filter(notes == "leanging hardwood snag")

figureDpi = 300
minimumHeightClass = 2 # m
totalTileAreaHa = length(unique(treetopDataDsm$tile)) * (0.3048 * treetopOptions$tileSize)^2 / 10000

## Figure 01: distribution of naturally regenerated and plantation stands + treetop dataset tiles
# .jpegs sourced from Elliott.qgz layouts
standStratification2022 = jpeg::readJPEG("GIS/Trees/2015-16 cruise/natural regen and plantation stands.jpeg")
cruiseStands2016 = jpeg::readJPEG("GIS/Trees/2015-16 cruise/2015-16 cruise stands.jpeg")

ggplot() +
  ggpubr::background_image(standStratification2022) +
  coord_fixed(ratio = dim(standStratification2022)[1] / dim(standStratification2022)[2]) +
  labs(title = paste(plotLetters[1], "primary distribution of stand structure")) +
ggplot() +
  ggpubr::background_image(cruiseStands2016) +
  coord_fixed(ratio = dim(cruiseStands2016)[1] / dim(cruiseStands2016)[2]) +
  labs(title = paste(plotLetters[2], "winter 2015–16 ground inventory and treetop dataset")) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1)
#ggsave("trees/segmentation/treetops/figures/Figure 01 stands.jpg", quality = 90, height = 11, width = 20, units = "cm", dpi = figureDpi)


## Figure 02: dataset tiles from QT Modeler
s04200w06810 = jpeg::readJPEG("GIS/DOGAMI/2021 OLC Coos County/images/s04200w06810 with noise grey 4.3.jpg")
s04200w06840 = jpeg::readJPEG("GIS/DOGAMI/2021 OLC Coos County/images/s04200w06840 with noise grey 4.3.jpg")
s04230w06810 = jpeg::readJPEG("GIS/DOGAMI/2021 OLC Coos County/images/s04230w06810 with noise grey 4.3.jpg")

ggplot() +
  ggpubr::background_image(s04200w06810) +
  coord_fixed(ratio = dim(s04200w06810)[1] / dim(s04200w06810)[2]) +
  labs(title = paste(plotLetters[1], "s04200w06810")) +
ggplot() +
  ggpubr::background_image(s04200w06840) +
  coord_fixed(ratio = dim(s04230w06810)[1] / dim(s04230w06810)[2]) +
  labs(title = paste(plotLetters[2], "s04200w06840")) +
ggplot() +
  ggpubr::background_image(s04230w06810) +
  coord_fixed(ratio = dim(s04200w06840)[1] / dim(s04200w06840)[2]) +
  labs(title = paste(plotLetters[3], "s04230w06810")) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 3)
#ggsave("trees/segmentation/treetops/figures/Figure 02 dataset tiles with noise.jpg", quality = 90, height = 23, width = 9.5, units = "cm", dpi = figureDpi)


## Figure 03: rings for DSM forest
distance = crossing(x = seq(-10, 10), y = seq(-10, 10)) %>% 
  mutate(distance = sqrt(x^2 + y^2),
         ring = round(distance, 0),
         angle = if_else(y >= 0, 0, 360) + 180 / pi * atan2(y, x),
         octant = 45 * if_else(angle <= 360 - 0.5 * 45, round(angle / 45), 0)) %>%
  filter(distance < 10.5) %>%
  relocate(ring) %>%
  arrange(ring, y, x)

ggplot() +
  geom_raster(aes(x = x, y = y, fill = as.factor(ring)), distance) + # geom_contour() not helpful here
  #geom_text(aes(x = x, y = y, label = ring, color = ring < 5), distance %>% filter(x >= 0, y == 0), size = 2.9) +
  coord_equal() +
  guides(color = "none") +
  labs(x = "x, surface model cells", y = "y, surface model cells", fill = "ring") +
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("black", "white")) +
  scale_fill_viridis_d(breaks = seq(0, 10), limits = rev, begin = 0.1)
#ggsave("trees/segmentation/treetops/figures/Figure 03 rings.png", height = 9, width = 11, units = "cm", dpi = figureDpi)
#write_xlsx(distance %>% filter(ring > 0), "trees/segmentation/rings.xlsx")


## Figure 04: DSM, CHM, and CMM surface comparison
# .jpegs sourced from Elliott ABA.qgz layouts
s04200w06810dsm = jpeg::readJPEG("trees/segmentation/treetops/figures/s04200w06810 DSM.jpeg")
s04200w06810chm = jpeg::readJPEG("trees/segmentation/treetops/figures/s04200w06810 CHM.jpeg")
s04200w06810cmm = jpeg::readJPEG("trees/segmentation/treetops/figures/s04200w06810 CMM.jpeg")

ggplot() +
  ggpubr::background_image(s04200w06810dsm) +
  coord_fixed(ratio = dim(s04200w06810dsm)[1] / dim(s04200w06810dsm)[2]) +
  labs(title = paste(plotLetters[1], "s04200w06810, 62 × 60 m DSM patch")) +
ggplot() +
  ggpubr::background_image(s04200w06810chm) +
  coord_fixed(ratio = dim(s04200w06810chm)[1] / dim(s04200w06810chm)[2]) +
  labs(title = paste(plotLetters[2], "s04200w06810, matching CHM patch")) +
ggplot() +
  ggpubr::background_image(s04200w06810cmm) +
  coord_fixed(ratio = dim(s04200w06810cmm)[1] / dim(s04200w06810cmm)[2]) +
  labs(title = paste(plotLetters[3], "s04230w06810, matching CMM patch")) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 3)
#ggsave("trees/segmentation/treetops/figures/Figure 04 s04200w06810 patch comparison.jpg", quality = 90, height = 23, width = 9.5, units = "cm", dpi = figureDpi)


## Figure 05: DSM distribution by height
localMaximaHistogramDsm = get_local_maxima_histogram(treetopDataDsm)
missingOrAmbiguousHistogramDsm = missingOrAmbiguous46dsm %>% mutate(heightClass = round(height)) %>% group_by(heightClass, notes) %>%
  summarize(trees = n(), .groups = "drop_last")
treetopProbabilityDsm = get_treetop_probability(treetopDataDsm)

overallOmmissionRate = missingOrAmbiguous46dsm %>% group_by(notes) %>% summarize(n = n()) %>% summarize(nOmitted = sum(if_else(notes != "point cloud ambiguous", n, 0.5 * n))) %>% # absent any more accurate indication, assume trees in half of ambiguous locations
  mutate(nTreetops = nrow(acceptedTreetops46dsm), omissionRatePct = 100 * nOmitted / (nTreetops + nOmitted))

plot_local_maxima_density(localMaximaHistogramDsm) +
plot_local_maxima_distribution(localMaximaHistogramDsm) +
ggplot() +
  geom_col(aes(x = trees / totalTileAreaHa, y = heightClass, alpha = heightClass >= 5, fill = notes, group = heightClass), missingOrAmbiguousHistogramDsm, orientation = "y", width = 1) +
  geom_segment(aes(x = 0, y = 4.5, xend = 1.5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, 90)) +
  labs(x = bquote("trees ha"^-1), y = NULL, alpha = NULL, fill = "missing and uncertain\ntrees", title = paste(plotLetters[3], "omissions")) +
  scale_fill_manual(breaks = c("obscured by branch", "obscured by noise", "local maxima absent", "point cloud ambiguous"), values = c("forestgreen", "red", "cyan", "grey80")) +
plot_treetop_probability(treetopProbabilityDsm) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, widths = c(1, 0.7, 0.5, 0.6), guides = "collect") &
  guides(alpha = "none") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1))
#ggsave("trees/segmentation/treetops/figures/Figure 05 DSM distribution by height class.png", height = 17, width = 20, units = "cm", dpi = figureDpi)


## Figure 06: CHM and CMM local maxima distribution and treetop decision boundaries
localMaximaHistogramChm = get_local_maxima_histogram(treetopDataChm)
localMaximaHistogramCmm = get_local_maxima_histogram(treetopDataCmm)

treetopProbabilityChm = get_treetop_probability(treetopDataChm)
treetopProbabilityCmm = get_treetop_probability(treetopDataCmm)

plot_local_maxima_density(localMaximaHistogramChm, plotLetter = paste("   ", plotLetters[1]), plotTitle = "CHM distribution") +
plot_treetop_probability(treetopProbabilityChm, fAtH = function(h) { return(0.51709010 + 0.11263928 * h^0.67990521) }, plotLetter = paste("  ", plotLetters[2]), plotTitle = "CHM treetops") +
plot_local_maxima_density(localMaximaHistogramCmm, plotLetter = plotLetters[3], plotTitle = "CMM distribution", yLabel = NULL) +
plot_treetop_probability(treetopProbabilityCmm, fAtH = function(h) { return(-0.22250373 + 0.05153255 * h^0.93888880) }, plotLetter = paste("  ", plotLetters[4]), plotTitle = "CMM treetops") +  
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, widths = c(1, 0.7, 1, 0.7), guides = "collect") &
  guides(alpha = "none") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1))
#ggsave("trees/segmentation/treetops/figures/Figure 06 CHM and CMM distribution by height class.png", height = 17, width = 20, units = "cm", dpi = figureDpi)


## Figure 07: treetop dataset density by height class
# PRELIMINARY: update to acceptedTreetops46chm and acceptedTreetops46cmm if/when those layers are manually reviewed and edited?
treetopsByHeight = left_join(left_join(acceptedTreetops46dsm %>% mutate(heightClass = round(0.3048 * height)) %>% group_by(heightClass) %>%
                                         summarize(treetopsPerHectareDsm = n() / totalTileAreaHa),
                                       treetopDataChm %>% filter(isTreetopRadius == TRUE) %>% mutate(heightClass = round(height)) %>% group_by(heightClass) %>%
                                         summarize(treetopsPerHectareChm = n() / totalTileAreaHa),
                                       by = join_by(heightClass)),
                             treetopDataCmm %>% filter(isTreetopRadius == TRUE) %>% mutate(heightClass = round(height)) %>% group_by(heightClass) %>%
                               summarize(treetopsPerHectareCmm = n() / totalTileAreaHa),
                             by = join_by(heightClass))

ggplot() +
  geom_col(aes(x = treetopsPerHectareDsm, y = heightClass, alpha = heightClass >= 5, fill = "DSM"), treetopsByHeight, orientation = "y", width = 1) +
  geom_segment(aes(x = 0, y = 4.5, xend = 25, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_trans(x = scales::pseudo_log_trans(), xlim = c(0, 25), ylim = c(0, 90)) +
  labs(x = bquote("treetops ha"^-1), y = "height above ground, m", alpha = NULL, color = NULL, fill = NULL, title = paste("    ", plotLetters[1], "DSM treetops")) +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20), minor_breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9, 30)) +
ggplot() +
  geom_col(aes(x = treetopsPerHectareChm, y = heightClass, alpha = heightClass >= 5, fill = "CHM"), treetopsByHeight, orientation = "y", width = 1) +
  geom_segment(aes(x = 0, y = 4.5, xend = 25, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_trans(x = scales::pseudo_log_trans(), xlim = c(0, 25), ylim = c(0, 90)) +
  labs(x = bquote("treetops ha"^-1), y = NULL, alpha = NULL, color = NULL, fill = NULL, title = paste(plotLetters[2], "CHM treetops")) +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20), minor_breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9, 30)) +
ggplot() +
  geom_col(aes(x = treetopsPerHectareCmm, y = heightClass, alpha = heightClass >= 5, fill = "CMM"), treetopsByHeight, orientation = "y", width = 1) +
  geom_segment(aes(x = 0, y = 4.5, xend = 25, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  coord_trans(x = scales::pseudo_log_trans(), xlim = c(0, 25), ylim = c(0, 90)) +
  labs(x = bquote("treetops ha"^-1), y = NULL, alpha = NULL, color = NULL, fill = NULL, title = paste(plotLetters[3], "CMM treetops")) +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20), minor_breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9, 30)) +
ggplot() +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.3) +
  geom_segment(aes(x = 0, y = 4.5, xend = 2, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  geom_line(aes(x = treetopsPerHectareCmm / treetopsPerHectareDsm, y = heightClass, alpha = heightClass >= 5, color = "CMM"), treetopsByHeight %>% filter(heightClass >= minimumHeightClass), orientation = "y") +
  geom_line(aes(x = treetopsPerHectareChm / treetopsPerHectareDsm, y = heightClass, alpha = heightClass >= 5, color = "CHM"), treetopsByHeight %>% filter(heightClass >= minimumHeightClass), orientation = "y") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 90)) +
  labs(x = bquote("treetops ha"^-1), y = NULL, alpha = NULL, color = NULL, fill = NULL, title = paste(plotLetters[4], "fraction of DSM tops")) +
  scale_x_continuous(labels = scales::percent) +
ggplot() + 
  geom_vline(xintercept = 1, color = "black", linewidth = 0.3) +
  geom_segment(aes(x = 0.8, y = 4.5, xend = 1.2, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  geom_line(aes(x = treetopsPerHectareCmm / treetopsPerHectareDsm, y = heightClass, alpha = heightClass >= 5, color = "CMM"), treetopsByHeight %>% filter(heightClass >= minimumHeightClass), orientation = "y") +
  geom_line(aes(x = treetopsPerHectareChm / treetopsPerHectareDsm, y = heightClass, alpha = heightClass >= 5, color = "CHM"), treetopsByHeight %>% filter(heightClass >= minimumHeightClass), orientation = "y") +
  coord_cartesian(xlim = c(0.85, 1.15), ylim = c(0, 90)) +
  labs(x = bquote("treetops ha"^-1), y = NULL, alpha = NULL, color = NULL, fill = NULL, title = paste(plotLetters[5], "expanded view")) +
  scale_x_continuous(breaks = seq(0.8, 1.2, by = 0.1), labels = scales::percent) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, widths = c(0.7, 0.7, 0.7, 1, 0.65), guides = "collect") &
  guides(alpha = "none", fill = "none") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_color_manual(breaks = c("DSM", "CHM", "CMM"), values = c("#609048FF", "#90A860FF", "#486030FF")) & # calecopal::redwood1 and 2, https://pmassicotte.github.io/paletteer_gallery/#qualitative
  scale_fill_manual(breaks = c("DSM", "CHM", "CMM"), values = c("#609048FF", "#90A860FF", "#486030FF")) &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1)) &
  theme(legend.margin = margin())
#ggsave("trees/segmentation/treetops/figures/Figure 07 dataset density by height class.png", height = 17, width = 20, units = "cm", dpi = figureDpi)


## Figure 08: median confusion matrices
# method       net tree count error, median %  treetops, M
# DSM forest  -1.63                             9.70 -> 9.86
# DSM radius   2.53                            10.19 -> 9.93
# CHM radius   2.42                             9.56 -> 9.33
# CMM radius  -4.32                             9.30 -> 9.70
rfDsmConfusionBinaryMedian = unnest_binary_confusion_median(randomForestAccuracy)
rfDsmConfusionQuinaryMedian = unnest_quinary_confusion_median(randomForestAccuracy)
radiusDsmConfusionMedian = unnest_binary_confusion_median(radiusDsmAccuracyPower)
radiusChmConfusionMedian = unnest_binary_confusion_median(radiusChmAccuracyPower)
radiusCmmConfusionMedian = unnest_binary_confusion_median(radiusCmmAccuracyPower)

medianErrorBySurface = bind_rows(rfDsmConfusionBinaryMedian %>% mutate(method = "DSM forest"),
                                 radiusDsmConfusionMedian %>% mutate(method = "DSM radius"),
                                 radiusChmConfusionMedian %>% mutate(method = "CHM radius"),
                                 radiusCmmConfusionMedian %>% mutate(method = "CMM radius")) %>% 
  mutate(surface = factor(method, levels = c("DSM forest", "DSM radius", "CHM radius", "CMM radius"))) %>%
  filter(prediction != reference) %>% 
  group_by(surface) %>%
  summarize(overallErrorPct = 100 * sum(fraction), netTreeCountErrorPct = -100 * diff(fraction))

ggplot() +
  geom_tile(aes(x = reference, y = prediction, fill = fraction), radiusDsmConfusionMedian) +
  geom_text(aes(x = reference, y = prediction, label = sprintf("%.1f%%", 100 * fraction), color = fraction > 0.60), rfDsmConfusionBinaryMedian, size = 2.7) +
  guides(fill = "none") +
  labs(x = "actual class", y = "predicted class", color = NULL, fill = "fraction of\nlocal maxima", title = paste("                    ", plotLetters[1], "DSM forest, merged"), subtitle = sprintf("                              %0.1f%% overall accuracy", 100 - medianErrorBySurface$overallErrorPct[which(medianErrorBySurface$surface == "DSM forest")])) +
ggplot() +
  geom_tile(aes(x = reference, y = prediction, fill = fraction), radiusDsmConfusionMedian) +
  geom_text(aes(x = reference, y = prediction, label = sprintf("%.1f%%", 100 * fraction), color = fraction > 0.60), radiusDsmConfusionMedian, size = 2.7) +
  guides(fill = "none") +
  labs(x = "actual class", y = NULL, color = NULL, fill = "fraction of\nlocal maxima", title = paste(plotLetters[2], "DSM radius"), subtitle = sprintf("      %0.1f%% overall accuracy", 100 - medianErrorBySurface$overallErrorPct[which(medianErrorBySurface$surface == "DSM radius")])) +
ggplot() +
  geom_tile(aes(x = reference, y = prediction, fill = fraction), radiusChmConfusionMedian) +
  geom_text(aes(x = reference, y = prediction, label = sprintf("%.1f%%", 100 * fraction), color = fraction > 0.60), radiusChmConfusionMedian, size = 2.7) +
  guides(fill = "none") +
  labs(x = "actual class", y = NULL, color = NULL, fill = "fraction of\nlocal maxima", title = paste(plotLetters[3], "CHM radius"), subtitle = sprintf("      %0.1f%% overall accuracy", 100 - medianErrorBySurface$overallErrorPct[which(medianErrorBySurface$surface == "CHM radius")])) +
ggplot() +
  geom_tile(aes(x = reference, y = prediction, fill = fraction), radiusCmmConfusionMedian) +
  geom_text(aes(x = reference, y = prediction, label = sprintf("%.1f%%", 100 * fraction), color = fraction > 0.60), radiusCmmConfusionMedian, size = 2.7) +
  guides(fill = "none") +
  labs(x = "actual class", y = NULL, color = NULL, fill = "fraction of\nlocal maxima", title = paste(plotLetters[4], "CMM radius"), subtitle = sprintf("      %0.1f%% overall accuracy", 100 - medianErrorBySurface$overallErrorPct[which(medianErrorBySurface$surface == "CMM radius")])) +
ggplot() +
  geom_tile(aes(x = reference, y = prediction, fill = if_else(fraction > 0, fraction, NA_real_)), rfDsmConfusionQuinaryMedian) +
  geom_text(aes(x = reference, y = prediction, label = if_else(fraction > 0, sprintf(if_else(fraction > 0.005, "%.1f%%", "%.1g%%"), 100 * fraction), "0%"), color = fraction > 0.60), rfDsmConfusionQuinaryMedian, size = 2.5) +
  labs(x = "actual class", y = "predicted class", color = NULL, fill = "fraction of\nlocal maxima", title = paste(plotLetters[5], "DSM forest, unmerged"), subtitle = sprintf("     %.1f%% overall accuracy", 100 * sum((rfDsmConfusionQuinaryMedian %>% filter(prediction == reference))$fraction))) +
  scale_x_discrete(breaks = c("treetop", "merge point", "other", "residual noise", "processing artifact"), labels = c("treetop", "merge\npoint", "other", "residual\nnoise", "proc.\nartifact")) +
  theme(legend.margin = margin(l = -105), plot.subtitle = element_text(hjust = 0.35), plot.title = element_text(hjust = 0.35)) +
plot_annotation(theme = theme(plot.margin = margin(l = -60))) +
plot_layout(design = "ABCD
EEEE", heights = c(2, 5)) &
  coord_fixed(ratio = 1.03) & # coord_equal() visually appears stretched due to cell labels
  guides(color = "none") &
  scale_color_manual(breaks = c(TRUE, FALSE), values = c("white", "black")) &
  paletteer::scale_fill_paletteer_c("ggthemes::Blue-Teal", labels = scales::percent, limits = c(0, 1), na.value = "white") &
  scale_y_discrete(limits = rev) &
  theme(panel.grid.major = element_blank())
#ggsave("trees/segmentation/treetops/figures/Figure 08 median confusion matrices.png", height = 12, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/segmentation/treetops/figures/Figure 08 median confusion matrices.svg", height = 12, width = 20, units = "cm", dpi = figureDpi)
#ggsave("trees/segmentation/treetops/figures/Figure 08 median confusion matrices.pdf", height = 12, width = 20, units = "cm", dpi = figureDpi)


## Figure 09: DSM and random forest accuracy distribution by height
radiusDsmAccuracyPowerByHeight = unnest_cv_accuracy_by_height(radiusDsmAccuracyPower)
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
  labs(x = "treetop detection\naccuracy", y = "height above ground, m", color = "local\nmaxima", title = paste("   ", plotLetters[1], "DSM radius")) +
  scale_x_continuous(labels = scales::percent) +
  theme(plot.title = element_text(vjust = 0.66)) +
ggplot() +
  geom_violin(aes(x = treetopAccuracy, y = heightClass, color = after_stat(count), group = heightClass, weight = meanN), randomForestAccuracyByHeight, draw_quantiles = c(0.5), linewidth = 0.3, width = 1.75) +
  geom_segment(aes(x = 0, y = 4.5, xend = 5, yend = 4.5), color = "grey20", linetype = "dashed", linewidth = 0.2) +
  guides(color = "none") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "treetop detection\naccuracy", y = NULL, color = "local\nmaxima", title = paste(plotLetters[2], "DSM forest")) +
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
  coord_cartesian(xlim = c(-0.1, 1)) +
  guides(alpha = "none", fill = guide_colorbar()) +
  scale_fill_gradientn(colors = c("#FF0000", "#FF0000", "#FF0000", "grey80", "#009000", "#009000", "#009000"), values = c(-100, -10, -2, 0, 2, 10, 100), breaks = c(-100, -10, -1, 0, 1, 10, 100), labels = scales::percent(c(-100, -10, -1, 0, 1, 10, 100), scale = 1), limits = c(-100, 100), rescaler = scales::rescale_none, transform = scales::pseudo_log_trans(sigma = 0.25)) +
  scale_x_continuous(breaks = c(-1, -0.1, 0, 0.1, 1), labels = scales::percent, minor_breaks = c(-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), transform = scales::pseudo_log_trans(sigma = 0.02, base = 10)) +
plot_annotation(theme = theme()) +
plot_layout(nrow = 1, widths = c(0.95, 0.95, 0.8, 1), guides = "collect") &
  scale_alpha_manual(breaks = c(TRUE, FALSE), values = c(1, 0.4)) &
  scale_color_gradient(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c(1, 10, 100, 1000, 10000, 100000), limits = c(1, NA), high = "#132B43", low = "#96F1FF", transform = "log10") &
  scale_y_continuous(breaks = seq(0, 90, by = 10), expand = c(0, 1))
#ggsave("trees/segmentation/treetops/figures/Figure 09 DSM accuracy by height class.png", height = 17, width = 20, units = "cm", dpi = figureDpi)


## Figure 10: random forest variable selection and importance
globalImportance = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m9n3 global importance.Rds")
localImportance = readRDS("trees/segmentation/treetops/random forest s4268 458k VSURF Pde m9n3 local importance.Rds")

ggplot() +
  geom_raster(aes(x = "global", y = label, fill = importance), globalImportance) +
  labs(x = NULL, y = NULL, title = paste(plotLetters[1], "global permutation importance")) +
  scale_y_discrete(limits = rev) +
ggplot() +
  geom_raster(aes(x = treetop, y = predictor, fill = importance), localImportance) +
  labs(x = NULL, y = NULL, title = paste(plotLetters[2], "local importance")) +
  scale_x_discrete(labels = c("single top", "merge point", "noise", "processing\nartifact", "other"), limits = c("yes", "merge", "noise", "maybe noise", "no")) +
  scale_y_discrete(labels = NULL, limits = rev) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(nrow = 1, ncol = 2, guides = "collect") &
  coord_equal() &
  labs(x = NULL, y = NULL, fill = "relative\nimportance, %") &
  scale_fill_viridis_c(option = "plasma", limits = c(0, 100 + 1E-14)) &
  theme(axis.text.x = element_text(angle = 90, lineheight = 0.67, hjust = 1, vjust = 0.5), legend.title = element_text(size = 10))
#ggsave("trees/segmentation/treetops/figures/Figure 10 DSM forest importance.png", height = 12, width = 12.5, units = "cm", dpi = figureDpi)

# runtimes
# Get-Dsm: 2,244,000,000 DSM cells from 561 tiles (49364.3 Mpoints) in 06:01: 1739.64 GB at 1.55 tiles/s (88.0 Mpoints/tile, 4.8 GB/s).
# Get-DsmSlopeAndAspect: Found slope and aspect in 561 tiles and generated .vrt in 00:21.
# Get-LocalMaxima: Found 99,598,584 DSM, 20,299,060 CMM, and 96,860,358 CHM maxima within 561 tiles in 19:26 (177,537 DSM maxima/tile).
# DSM optim 14.3329 secs + 1.435918 m 2x25 cross validation 
# clustered random forest: Boruta 36.22h + VSURF + 39.21h + tune 8.823h + cross validation 1.528h + 4.409m fit


## Table 01: DSM dataset content
# Legacy loads as dataset hasn't been updated to the slight (~0.001%) cell elevation differences between the v3 beta and v3 DSMs.
localMaxima46dsm = bind_rows(get_treetop_eligible_maxima("s04200w06810", acceptedTileName = "s04200w06810", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                             get_treetop_eligible_maxima("s04200w06840", acceptedTileName = "s04200w06840", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"),
                             get_treetop_eligible_maxima("s04230w06810", acceptedTileName = "s04230w06810", localMaximaPath = localMaximaPathV3beta, localMaximaLayer = "localMaxima"))
#localMaxima46dsm %>% group_by(tile) %>% summarize(`single top` = sum(treetop == "yes"), `merge point` = sum(treetop == "merge"), `residual noise` = sum(treetop == "noise"), other = sum(treetop == "no"))

stats46dsm = left_join(left_join(localMaxima46dsm %>% group_by(tile, treetop) %>% summarize(maxima = n(), .groups = "drop") %>% 
                                   pivot_wider(id_cols = "tile", names_from = "treetop", values_from = "maxima") %>%
                                   mutate(`total maxima` = no + yes + merge + noise + `maybe noise`,
                                          `processing artifact` = c(0, 2, 4), # s04200w06810, s04230w06810, s04200w06840 
                                          no = no + `maybe noise` - `processing artifact`),
                       acceptedTreetops46dsm %>% group_by(tile) %>% summarize(treetops = n()),
                       by = join_by(tile)),
          missingOrAmbiguous46dsm %>% group_by(tile, notes) %>% summarize(trees = n(), .groups = "drop") %>%
            pivot_wider(id_cols = "tile", names_from = "notes", values_from = "trees") %>%
            mutate(`obscured by noise` = replace_na(`obscured by noise`, 0)),
          by = join_by(tile)) %>%
  bind_rows(summarize(., across(where(is.numeric), sum))) %>%
  rename(`single top` = yes, `merge point` = merge, `residual noise` = noise, other = no) %>%
  mutate(surface = "DSM", tile = replace_na(tile, "total")) %>%
  relocate(surface, tile, `total maxima`, treetops, `single top`, `merge point`, `residual noise`, `processing artifact`, other, `maybe noise`)
stats46dsm
#write_xlsx(stats46dsm, "trees/segmentation/treetops/treetop DSM dataset.xlsx")
 
stats46chm = left_join(treetopDataChm %>% group_by(tile, treetop) %>% summarize(maxima = n(), .groups = "drop") %>%
                         pivot_wider(id_cols = "tile", names_from = "treetop", values_from = "maxima") %>%
                         mutate(`total maxima` = no + yes + merge + noise + `maybe noise`,
                                `processing artifact` = 0,
                                no = no + `maybe noise` - `processing artifact`),
                       acceptedTreetops46chm %>% group_by(tile) %>% summarize(treetops = n()),
                       by = join_by(tile)) %>%
  bind_rows(summarize(., across(where(is.numeric), sum))) %>%
  rename(`single top` = yes, `merge point` = merge, `residual noise` = noise, other = no) %>%
  mutate(surface = "CHM", tile = replace_na(tile, "total")) %>%
  relocate(surface, tile, `total maxima`, treetops, `single top`, `merge point`, `residual noise`, `processing artifact`, other, `maybe noise`)
stats46chm
#write_xlsx(stats46chm, "trees/segmentation/treetops/treetop CHM dataset.xlsx")

stats46cmm = left_join(treetopDataCmm %>% group_by(tile, treetop) %>% summarize(maxima = n(), .groups = "drop") %>%
                         pivot_wider(id_cols = "tile", names_from = "treetop", values_from = "maxima") %>%
                         mutate(`maybe noise` = replace_na(`maybe noise`, 0),
                                `total maxima` = no + yes + merge + noise + `maybe noise`,
                                `processing artifact` = 0,
                                no = no + `maybe noise` - `processing artifact`),
                       acceptedTreetops46cmm %>% group_by(tile) %>% summarize(treetops = n()),
                       by = join_by(tile)) %>%
  bind_rows(summarize(., across(where(is.numeric), sum))) %>%
  rename(`single top` = yes, `merge point` = merge, `residual noise` = noise, other = no) %>%
  mutate(surface = "CMM", tile = replace_na(tile, "total")) %>%
  relocate(surface, tile, `total maxima`, treetops, `single top`, `merge point`, `residual noise`, `processing artifact`, other, `maybe noise`)
stats46cmm
#write_xlsx(stats46cmm, "trees/segmentation/treetops/treetop CMM dataset.xlsx")


## investigatory, prototypes
if (treetopOptions$includeInvestigatory)
{
  # ranges are small so violins are uninformative
  ggplot() +
    geom_violin(aes(x = FALSE, y = trueNegative / total), radiusDsmConfusion, draw_quantiles = c(0.25, 0.5, 0.75)) +
    stat_summary(aes(x = FALSE, y = trueNegative / total), radiusDsmConfusion, fun = mean, geom = "point") +
    labs(x = NULL, y = "true negatives") +
    ggplot() +
    geom_violin(aes(x = TRUE, y = falseNegative / total), radiusDsmConfusion, draw_quantiles = c(0.25, 0.5, 0.75)) +
    stat_summary(aes(x = TRUE, y = falseNegative / total), radiusDsmConfusion, fun = mean, geom = "point") +
    labs(x = NULL, y = "false negatives") +
    ggplot() +
    geom_violin(aes(x = FALSE, y = falsePositive / total), radiusDsmConfusion, draw_quantiles = c(0.25, 0.5, 0.75)) +
    stat_summary(aes(x = FALSE, y = falsePositive / total), radiusDsmConfusion, fun = mean, geom = "point") +
    labs(x = NULL, y = "false positives") +
    ggplot() +
    geom_violin(aes(x = TRUE, y = truePositive / total), radiusDsmConfusion, draw_quantiles = c(0.25, 0.5, 0.75)) +
    stat_summary(aes(x = TRUE, y = truePositive / total), radiusDsmConfusion, fun = mean, geom = "point") +
    labs(x = NULL, y = "true positives") +
    plot_annotation(theme = theme(plot.margin = margin())) +
    plot_layout(nrow = 2, ncol = 2, guides = "collect") &
    coord_cartesian(ylim = c(0, 1)) &
    scale_x_discrete(labels = NULL) &
    scale_y_continuous(labels = scales::percent)
}
