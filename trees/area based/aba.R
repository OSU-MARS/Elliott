library(arrow)
library(dplyr)
library(FNN)
library(ggplot2)
library(furrr)
library(magrittr)
library(patchwork)
library(rsample)
library(terra)
library(tidyr)
library(readxl)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             #legend.key.height = unit(0.85, "line"),
                             legend.spacing.y = unit(0.3, "line"),
                             legend.title = element_text(size = 10),
                             panel.border = element_blank(), 
                             plot.title = element_text(size = 10)))

abaGrid = tibble(originX = 106020, originY = 190480, size = 20) %>% # m EPSG:6556, ABA grid origin and cell size from GIS/Trees/Elliott ABA grid 20 m.gpkg
  mutate(sizeHa = size^2 / 10000)

stands2022 = read_xlsx("GIS/Trees/2015-16 cruise.xlsx") %>% rename(stand = standID2016)
trees2021organon = left_join(read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE) %>% select(-species), # TODO: should 2016 snags be joined? they don't flow through Organon but standing in 2016 won't all have fallen by 2021
                             read_xlsx("trees/Elliott final cruise records 2015-16.xlsx", sheet = "CRUISERECS") %>% rename(stand = StandID, plot = PlotID, tag = TreeID, species = Species) %>% select(stand, plot, tag, species),
                             by = c("stand", "plot", "tag")) %>% # recover undubbed species from original plot measurements
  #mutate(species = case_match(species, 202 ~ "PSME", 351 ~ "ALRU", 263 ~ "TSHE", 312 ~ "ACMA", 242 ~ "THPL", 17 ~ "ABGR", 361 ~ "ARME", 431 ~ "CHCH", 631 ~ "NODE", 920 ~ "Salix", 81 ~ "CADE", 492 ~ "CONU", 815 ~ "QUGA", )) %>%
  filter(year == 2021)

plotHeights = left_join(trees2021organon, stands2022 %>% select(stand, measurePlotsInStand), by = c("stand")) %>% # ~12 seconds, ~9 s due to uncount()
  mutate(treesPerCell = abaGrid$sizeHa * measurePlotsInStand * liveExpansionFactor, # expansion factors from Organon are stand-level; multiplying by the number of measure plots converts the expansion back to plot-level
         treesPerCell = if_else(treesPerCell < 1, 1, round(treesPerCell))) %>% # quantize number of trees per cell for unnest(); since these are 2016 plot measurements grown to 2021 by default at least one tree must be present on the plot to measure and thus at least one tree is present on the grid cell---this does not hold where a tree or snag falls or windsnaps but that data's not available, so assume 100% upright
  uncount(treesPerCell) %>% 
  group_by(plot) %>% arrange(desc(height)) %>%
  summarize(n = n(), species1 = species[1], height1 = height[1], # for now, assume all mortality falls immediately
            species2 = if_else(n() >= 2, species[2], NA_character_), height2 = if_else(n() >= 2, height[2], NA_real_),
            species3 = if_else(n() >= 3, species[3], NA_character_), height3 = if_else(n() >= 3, height[3], NA_real_),
            species4 = if_else(n() >= 4, species[4], NA_character_), height4 = if_else(n() >= 4, height[4], NA_real_),
            species5 = if_else(n() >= 5, species[5], NA_character_), height5 = if_else(n() >= 5, height[5], NA_real_),
            species6 = if_else(n() >= 6, species[6], NA_character_), height6 = if_else(n() >= 6, height[6], NA_real_),
            species7 = if_else(n() >= 7, species[7], NA_character_), height7 = if_else(n() >= 7, height[7], NA_real_),
            species8 = if_else(n() >= 8, species[8], NA_character_), height8 = if_else(n() >= 8, height[8], NA_real_),
            species9 = if_else(n() >= 9, species[9], NA_character_), height9 = if_else(n() >= 9, height[9], NA_real_),
            species10 = if_else(n() >= 10, species[10], NA_character_), height10 = if_else(n() >= 10, height[10], NA_real_),
            species11 = if_else(n() >= 11, species[11], NA_character_), height11 = if_else(n() >= 11, height[11], NA_real_),
            species12 = if_else(n() >= 12, species[12], NA_character_), height12 = if_else(n() >= 12, height[12], NA_real_),
            species13 = if_else(n() >= 13, species[13], NA_character_), height13 = if_else(n() >= 13, height[13], NA_real_),
            species14 = if_else(n() >= 14, species[14], NA_character_), height14 = if_else(n() >= 14, height[14], NA_real_),
            species15 = if_else(n() >= 15, species[15], NA_character_), height15 = if_else(n() >= 15, height[15], NA_real_),
            .groups = "drop") %>%
  mutate(isConifer1 = species1 %in% c("PSME"), isConifer2 = species2 %in% c("PSME"), isConifer3 = species3 %in% c("PSME"),
         isConifer4 = species4 %in% c("PSME"), isConifer5 = species5 %in% c("PSME"), isConifer6 = species6 %in% c("PSME"),
         isConifer7 = species7 %in% c("PSME"), isConifer8 = species8 %in% c("PSME"), isConifer9 = species9 %in% c("PSME"), 
         isConifer10 = species10 %in% c("PSME"), isConifer11 = species11 %in% c("PSME"), isConifer12 = species12 %in% c("PSME"), 
         isConifer13 = species13 %in% c("PSME"), isConifer14 = species14 %in% c("PSME"), isConifer15 = species15 %in% c("PSME"))

plots2016 = read_xlsx("GIS/Trees/2015-16 cruise/CruisePlots_All_20151211.xlsx")
plotMetrics2021 = list() # plot metrics created by metricsJob.R
for (chunkIndex in 1:19)
{
  plotMetrics2021[[chunkIndex]] = vect(file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics", paste0("Elliott stdmetrics plot chunk ", chunkIndex, ".gpkg")))
}
plotMetrics2021 = as_tibble(vect(plotMetrics2021)) %>% # EPSG:6557 but made nonspatial
  filter(PltInteger %in% plots2016$PltInteger, PlotType != "Count") %>% # count plots aren't directly useful to tree list establishment, also exclude dropped plots with known coordinates
  select(-Space, -STAND, -TotalPlots, -Count, -Measure, -Cruiser, -Plot_ID, -PlotType, -stems, -liveTrees, -snags, -primarySpecies, -n, -area) %>% # remove non-LiDAR observable or non-relevant variables
  select(-x, -y) %>% # prevent classification from mapping plots
  select(-zentropy) %>% # always NaN; lidR bug
  rename(plot = PltInteger, zMax = zmax, zMean = zmean, zStdDev = zsd, zSkew = zskew, zKurtosis = zkurt,
         pZaboveZmean = pzabovezmean, pZaboveThreshold = pzabove2, zQ05 = zq5, zQ10 = zq10, zQ15 = zq15, zQ20 = zq20, zQ25 = zq25, zQ30 = zq30, zQ35 = zq35, zQ40 = zq40, zQ45 = zq45, zQ50 = zq50, zQ55 = zq55, zQ60 = zq60, zQ65 = zq65, zQ70 = zq70, zQ75 = zq75, zQ80 = zq80, zQ85 = zq85, zQ90 = zq90, zQ95 = zq95,
         zPcumulative10 = zpcum1, zPcumulative20 = zpcum2, zPcumulative30 = zpcum3, zPcumulative40 = zpcum4, zPcumulative50 = zpcum5, zPcumulative60 = zpcum6, zPcumulative70 = zpcum7, zPcumulative80 = zpcum8, zPcumulative90 = zpcum9,
         intensityTotal = itot, intensityMax = imax, intensityMean = imean, intensityStdDev = isd, intensitySkew = iskew, intensityKurtosis = ikurt, intensityPground = ipground,
         pCumulativeZQ10 = ipcumzq10, pCumulativeZQ30 = ipcumzq30, pCumulativeZQ50 = ipcumzq50, pCumulativeZQ70 = ipcumzq70, pCumulativeZQ90 = ipcumzq90, 
         pFirstReturn = p1th, pSecondReturn = p2th, pThirdRetunr = p3th, pFourthReturn = p4th, pFifthReturn = p5th, pGround = pground)

plotHeights = left_join(plotHeights, plotMetrics2021, by = c("plot")) # join LiDAR metrics to Organon grown cruise data

load("trees/height-diameter/data/trees DSM ring.Rdata") # LiDAR identified treetops from trees.R, a few seconds
trees2021lidar = elliottTreesMod %>% # EPSG:6556
  mutate(abaGridX = floor(1/abaGrid$size * (x - abaGrid$originX)), abaGridY = floor(1/20 * (y - abaGrid$originY))) # ABA grid origin and cell size from GIS/Trees/Elliott ABA grid 20 m.gpkg
rm(elliottTreesMod) # since no way to load a single variable in a .Rdata into a specific variable name


## kNN match of ABA grid cells by heights of detected trees
# 1-64 treetops detected per grid cell, 1-20 trees per plot
if (recalcAbaCellOccupancy)
{
  startTime = Sys.time()
  abaCells = trees2021lidar %>% group_by(abaGridX, abaGridY) %>% arrange(desc(height)) %>% # ~48 minutes; 11.4 M treetops -> million row tibble
    summarize(n = n(), species1 = species[1], height1 = height[1], 
              species2 = if_else(n() >= 2, species[2], NA_character_), height2 = if_else(n() >= 2, height[2], NA_real_),
              species3 = if_else(n() >= 3, species[3], NA_character_), height3 = if_else(n() >= 3, height[3], NA_real_),
              species4 = if_else(n() >= 4, species[4], NA_character_), height4 = if_else(n() >= 4, height[4], NA_real_),
              species5 = if_else(n() >= 5, species[5], NA_character_), height5 = if_else(n() >= 5, height[5], NA_real_),
              species6 = if_else(n() >= 6, species[6], NA_character_), height6 = if_else(n() >= 6, height[6], NA_real_),
              species7 = if_else(n() >= 7, species[7], NA_character_), height7 = if_else(n() >= 7, height[7], NA_real_),
              species8 = if_else(n() >= 8, species[8], NA_character_), height8 = if_else(n() >= 8, height[8], NA_real_),
              species9 = if_else(n() >= 9, species[9], NA_character_), height9 = if_else(n() >= 9, height[9], NA_real_),
              species10 = if_else(n() >= 10, species[10], NA_character_), height10 = if_else(n() >= 10, height[10], NA_real_),
              species11 = if_else(n() >= 11, species[11], NA_character_), height11 = if_else(n() >= 11, height[11], NA_real_),
              species12 = if_else(n() >= 12, species[12], NA_character_), height12 = if_else(n() >= 12, height[12], NA_real_),
              species13 = if_else(n() >= 13, species[13], NA_character_), height13 = if_else(n() >= 13, height[13], NA_real_),
              species14 = if_else(n() >= 14, species[14], NA_character_), height14 = if_else(n() >= 14, height[14], NA_real_),
              species15 = if_else(n() >= 15, species[15], NA_character_), height15 = if_else(n() >= 15, height[15], NA_real_),
              .groups = "drop") %>%
    mutate(isConifer1 = species1 %in% c("PSME"), isConifer2 = species2 %in% c("PSME"), isConifer3 = species3 %in% c("PSME"),
           isConifer4 = species4 %in% c("PSME"), isConifer5 = species5 %in% c("PSME"), isConifer6 = species6 %in% c("PSME"),
           isConifer7 = species7 %in% c("PSME"), isConifer8 = species8 %in% c("PSME"), isConifer9 = species9 %in% c("PSME"), 
           isConifer10 = species10 %in% c("PSME"), isConifer11 = species11 %in% c("PSME"), isConifer12 = species12 %in% c("PSME"), 
           isConifer13 = species13 %in% c("PSME"), isConifer14 = species14 %in% c("PSME"), isConifer15 = species15 %in% c("PSME"))
  Sys.time() - startTime
  save(file = "trees/height-diameter/data/trees by cell.Rdata", abaCells)
} else {
  load("trees/height-diameter/data/trees by cell.Rdata") # abaCells
}

abaMetrics = as.data.frame(project(rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics/grid metrics 20 m.tif"), crs("epsg:6556"), threads = TRUE), xy = TRUE, na.rm = NA) # ~5 s
names(abaMetrics)[47] = "intensityPground" # temporary workaround for C# typo
abaMetrics %<>% mutate(abaGridX = floor(1/abaGrid$size * (x - abaGrid$originX)), abaGridY = floor(1/20 * (y - abaGrid$originY)))

abaCells = left_join(abaCells, abaMetrics %>% rename(nPoints = n), by = c("abaGridX", "abaGridY"))

coniferDeciduousWeight = 5
abaCellsByMatches = abaCells %>% mutate(treesMatched = if_else(n < 8, n, 8),
                                        isConifer1 = coniferDeciduousWeight * isConifer1,
                                        isConifer2 = coniferDeciduousWeight * isConifer2,
                                        isConifer3 = coniferDeciduousWeight * isConifer3,
                                        isConifer4 = coniferDeciduousWeight * isConifer4,
                                        isConifer5 = coniferDeciduousWeight * isConifer5,
                                        isConifer6 = coniferDeciduousWeight * isConifer6,
                                        isConifer7 = coniferDeciduousWeight * isConifer7,
                                        isConifer8 = coniferDeciduousWeight * isConifer8,
                                        isConifer9 = coniferDeciduousWeight * isConifer9,
                                        isConifer10 = coniferDeciduousWeight * isConifer10,
                                        isConifer11 = coniferDeciduousWeight * isConifer11,
                                        isConifer12 = coniferDeciduousWeight * isConifer12,
                                        isConifer13 = coniferDeciduousWeight * isConifer13,
                                        isConifer14 = coniferDeciduousWeight * isConifer14,
                                        isConifer15 = coniferDeciduousWeight * isConifer15,
                                        intensityMean = 0.001 * treesMatched * intensityMean,
                                        pGround = 100 * treesMatched * pGround,
                                        pZaboveThreshold = 100 * treesMatched * pZaboveThreshold,
                                        zMean = treesMatched * zMean,
                                        zQ10 = 10 * treesMatched * zQ10,
                                        zQ20 = 5 * treesMatched * zQ20,
                                        zQ30 = 3.33 * treesMatched * zQ30,
                                        zQ40 = 2.5 * treesMatched * zQ40,
                                        zQ50 = 2 * treesMatched * zQ50,
                                        zQ80 = 1.25 * treesMatched * zQ80) %>%
  filter(is.na(zMean) == FALSE) %>% # 386 cells without LiDAR metrics due to tight south side flight boundary in 2021
  nest(.by = treesMatched) %>% arrange(treesMatched)
plotHeightsScaled = plotHeights %>% mutate(treesMatched = if_else(n < 8, n, 8),
                                           isConifer1 = coniferDeciduousWeight * isConifer1,
                                           isConifer2 = coniferDeciduousWeight * isConifer2,
                                           isConifer3 = coniferDeciduousWeight * isConifer3,
                                           isConifer4 = coniferDeciduousWeight * isConifer4,
                                           isConifer5 = coniferDeciduousWeight * isConifer5,
                                           isConifer6 = coniferDeciduousWeight * isConifer6,
                                           isConifer7 = coniferDeciduousWeight * isConifer7,
                                           isConifer8 = coniferDeciduousWeight * isConifer8,
                                           isConifer9 = coniferDeciduousWeight * isConifer9,
                                           isConifer10 = coniferDeciduousWeight * isConifer10,
                                           isConifer11 = coniferDeciduousWeight * isConifer11,
                                           isConifer12 = coniferDeciduousWeight * isConifer12,
                                           isConifer13 = coniferDeciduousWeight * isConifer13,
                                           isConifer14 = coniferDeciduousWeight * isConifer14,
                                           isConifer15 = coniferDeciduousWeight * isConifer15,
                                           intensityMean = 0.001 * treesMatched * intensityMean,
                                           pGround = 100 * treesMatched * pGround,
                                           pZaboveThreshold = 100 * treesMatched * pZaboveThreshold,
                                           zMean = treesMatched * zMean,
                                           zQ10 = 10 * treesMatched * zQ10,
                                           zQ20 = 5 * treesMatched * zQ20,
                                           zQ30 = 3.33 * treesMatched * zQ30,
                                           zQ40 = 2.5 * treesMatched * zQ40,
                                           zQ50 = 2 * treesMatched * zQ50, 
                                           zQ80 = 1.25 * treesMatched * zQ80) %>%
  filter(is.na(zMean) == FALSE) # 30 plots without LiDAR metrics due to target coordinates not being available

plotHeights1 = plotHeightsScaled %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1) # 9990 plots, plotHeights %>% group_by(n) %>% summarize(n = n())
abaCells1 = abaCellsByMatches$data[[1]] %>% filter(n == 1) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1) # 14547 cells with one treetop, abaCells %>% group_by(n) %>% summarize(n = n())
cellKnnMatches1 = get.knnx(plotHeights1 %>% select(-plot), abaCells1 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights2 = plotHeightsScaled %>% filter(n >= 2) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2) # 9449 plots with two treetops
abaCells2 = abaCellsByMatches$data[[2]] %>% filter(n == 2) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2) # 38063 cells
cellKnnMatches2 = get.knnx(plotHeights2 %>% select(-plot), abaCells2 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights3 = plotHeightsScaled %>% filter(n >= 3) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3) # 8615 plots
abaCells3 = abaCellsByMatches$data[[3]] %>% filter(n == 3) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3) # 72586 cells
cellKnnMatches3 = get.knnx(plotHeights3 %>% select(-plot), abaCells3 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights4 = plotHeightsScaled %>% filter(n >= 4) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4) # 7465 plots
abaCells4 = abaCellsByMatches$data[[4]] %>% filter(n == 4) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4) # 95572 cells 
cellKnnMatches4 = get.knnx(plotHeights4 %>% select(-plot), abaCells4 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights5 = plotHeightsScaled %>% filter(n >= 5) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5) # 5912 plots
abaCells5 = abaCellsByMatches$data[[5]] %>% filter(n == 5) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5) # 98672 cells
cellKnnMatches5 = get.knnx(plotHeights5 %>% select(-plot), abaCells5 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights6 = plotHeightsScaled %>% filter(n >= 6) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6) # 4123 plots
abaCells6 = abaCellsByMatches$data[[6]] %>% filter(n == 6) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6) # 91485 cells 
cellKnnMatches6 = get.knnx(plotHeights6 %>% select(-plot), abaCells6 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights7 = plotHeightsScaled %>% filter(n >= 7) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7) # 2587 plots
abaCells7 = abaCellsByMatches$data[[7]] %>% filter(n == 7) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7) # 82153 cells 
cellKnnMatches7 = get.knnx(plotHeights7 %>% select(-plot), abaCells7 %>% select(-abaGridX, -abaGridY), k = 3)

plotHeights8 = plotHeightsScaled %>% filter(n >= 8) %>% select(plot, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8) # 1669 plots
abaCells8 = abaCellsByMatches$data[[8]] %>% filter(n >= 8) %>% select(abaGridX, abaGridY, pGround, zQ10, zQ20, zQ30, isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8) # 507,011 cells 
cellKnnMatches8 = get.knnx(plotHeights8 %>% select(-plot), abaCells8 %>% select(-abaGridX, -abaGridY), k = 3)

abaCellToPlots = bind_rows(abaCells1 %>% mutate(plot1 = plotHeights1$plot[cellKnnMatches1$nn.index[, 1]], plot2 = plotHeights1$plot[cellKnnMatches1$nn.index[, 2]], plot3 = plotHeights1$plot[cellKnnMatches1$nn.index[, 3]],
                                                distance1 = cellKnnMatches1$nn.dist[, 1], distance2 = cellKnnMatches1$nn.dist[, 2], distance3 = cellKnnMatches1$nn.dist[, 3]),
                           abaCells2 %>% mutate(plot1 = plotHeights2$plot[cellKnnMatches2$nn.index[, 1]], plot2 = plotHeights2$plot[cellKnnMatches2$nn.index[, 2]], plot3 = plotHeights2$plot[cellKnnMatches2$nn.index[, 3]],
                                                distance1 = cellKnnMatches2$nn.dist[, 1], distance2 = cellKnnMatches2$nn.dist[, 2], distance3 = cellKnnMatches2$nn.dist[, 3]),
                           abaCells3 %>% mutate(plot1 = plotHeights3$plot[cellKnnMatches3$nn.index[, 1]], plot2 = plotHeights3$plot[cellKnnMatches3$nn.index[, 2]], plot3 = plotHeights3$plot[cellKnnMatches3$nn.index[, 3]],
                                                distance1 = cellKnnMatches3$nn.dist[, 1], distance2 = cellKnnMatches3$nn.dist[, 2], distance3 = cellKnnMatches3$nn.dist[, 3]),
                           abaCells4 %>% mutate(plot1 = plotHeights4$plot[cellKnnMatches4$nn.index[, 1]], plot2 = plotHeights4$plot[cellKnnMatches4$nn.index[, 2]], plot3 = plotHeights4$plot[cellKnnMatches4$nn.index[, 3]],
                                                distance1 = cellKnnMatches4$nn.dist[, 1], distance2 = cellKnnMatches4$nn.dist[, 2], distance3 = cellKnnMatches4$nn.dist[, 3]),
                           abaCells5 %>% mutate(plot1 = plotHeights5$plot[cellKnnMatches5$nn.index[, 1]], plot2 = plotHeights5$plot[cellKnnMatches5$nn.index[, 2]], plot3 = plotHeights5$plot[cellKnnMatches5$nn.index[, 3]],
                                                distance1 = cellKnnMatches5$nn.dist[, 1], distance2 = cellKnnMatches5$nn.dist[, 2], distance3 = cellKnnMatches5$nn.dist[, 3]),
                           abaCells6 %>% mutate(plot1 = plotHeights6$plot[cellKnnMatches6$nn.index[, 1]], plot2 = plotHeights6$plot[cellKnnMatches6$nn.index[, 2]], plot3 = plotHeights6$plot[cellKnnMatches6$nn.index[, 3]],
                                                distance1 = cellKnnMatches6$nn.dist[, 1], distance2 = cellKnnMatches6$nn.dist[, 2], distance3 = cellKnnMatches6$nn.dist[, 3]),
                           abaCells7 %>% mutate(plot1 = plotHeights7$plot[cellKnnMatches7$nn.index[, 1]], plot2 = plotHeights7$plot[cellKnnMatches7$nn.index[, 2]], plot3 = plotHeights7$plot[cellKnnMatches7$nn.index[, 3]],
                                                distance1 = cellKnnMatches7$nn.dist[, 1], distance2 = cellKnnMatches7$nn.dist[, 2], distance3 = cellKnnMatches7$nn.dist[, 3]),
                           abaCells8 %>% mutate(plot1 = plotHeights8$plot[cellKnnMatches8$nn.index[, 1]], plot2 = plotHeights8$plot[cellKnnMatches8$nn.index[, 2]], plot3 = plotHeights8$plot[cellKnnMatches8$nn.index[, 3]],
                                                distance1 = cellKnnMatches8$nn.dist[, 1], distance2 = cellKnnMatches8$nn.dist[, 2], distance3 = cellKnnMatches8$nn.dist[, 3]))

abaCellPlots = left_join(left_join(abaCells, abaCellToPlots %>% select(abaGridX, abaGridY, plot1, plot2, plot3, distance1, distance2, distance3), by = c("abaGridX", "abaGridY")),
                         plotHeights %>% rename(plot1 = plot, plotN = n, plotSpecies1 = species1, plotSpecies2 = species3, plotSpecies3 = species3, plotSpecies4 = species4, plotSpecies5 = species5, plotSpecies6 = species6, plotSpecies7 = species7, plotSpecies8 = species8, plotHeight1 = height1, plotHeight2 = height2, plotHeight3 = height3, plotHeight4 = height4, plotHeight5 = height5, plotHeight6 = height6, plotHeight7 = height7, plotHeight8 = height8, plotIsConifer1 = isConifer1, plotIsConifer2 = isConifer2, plotIsConifer3 = isConifer3, plotIsConifer4 = isConifer4, plotIsConifer5 = isConifer5, plotIsConifer6 = isConifer6, plotIsConifer7 = isConifer7, plotIsConifer8 = isConifer8) %>% # could use rename_with() but calling a camel casing function with regex column matching isn't simpler than writing out 24 renames
                           select(plot1, plotN, starts_with("plotHeight"), starts_with("plotSpecies"), starts_with("plotIsConifer")), 
                         by = c("plot1")) %>%
  mutate(treesMatched = if_else(n < 8, n, 8),
         speciesError = isConifer1 != plotIsConifer1 + if_else(is.na(isConifer2), 0, isConifer2 != plotIsConifer2) + if_else(is.na(isConifer3), 0, isConifer3 != plotIsConifer3) + if_else(is.na(isConifer4), 0, isConifer4 != plotIsConifer4) + if_else(is.na(isConifer5), 0, isConifer5 != plotIsConifer5) + if_else(is.na(isConifer6), 0, isConifer6 != plotIsConifer6) + if_else(is.na(isConifer7), 0, isConifer7 != plotIsConifer7) + if_else(is.na(isConifer8), 0, isConifer8 != plotIsConifer8),
         heightErrorMae = abs(height1 - plotHeight1) + if_else(is.na(height2), 0, abs(height2 - plotHeight2)) + if_else(is.na(height3), 0, abs(height3 - plotHeight3)) + if_else(is.na(height4), 0, abs(height4 - plotHeight4)) + if_else(is.na(height5), 0, abs(height5 - plotHeight5)) + if_else(is.na(height6), 0, abs(height6 - plotHeight6)) + if_else(is.na(height7), 0, abs(height7 - plotHeight7)) + if_else(is.na(height8), 0, abs(height8 - plotHeight8)))

# 11.4 M treetops detected in matched cells: 40,515 ha in matchable cells, 40,532 ha total ABA grid area
# Elliott (including Hakki): 33727 ha -> 9.5 M LiDAR treetops, 20.7 M imputed treetops
abaCellPlots %>% summarize(cells = nrow(abaCellPlots), totalArea = abaGrid$sizeHa * cells, metricsCells = sum(is.na(zMean) == FALSE), matchedCells = sum(is.na(plotN) == FALSE), matchedArea = matchedCells * abaGrid$sizeHa, treetops = sum(n), matchedTreetops = sum((is.na(plotN) == FALSE) * n), plotTreeEstimate = sum(plotN, na.rm = TRUE)) %>% mutate(treeCountError = matchedTreetops - plotTreeEstimate)


ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 125, yend = 125), color = "grey70", linetype = "longdash", linewidth = 0.3) +
  geom_bin2d(aes(x = n, y = plotN), abaCellPlots, binwidth = c(1, 1)) +
  labs(x = "detected treetops", y = bquote("plot estimate, trees cell"^-1), fill = "plots") +
  scale_fill_viridis_c(labels = scales::label_comma(), trans = "log10")

ggplot() +
  geom_histogram(aes(x = n - plotN, y = after_stat(..count.. / sum(..count..))), abaCellPlots, binwidth = 2) +
  labs(x = "tree count error", y = "probability") +
ggplot() +
  geom_histogram(aes(x = speciesError / treesMatched, y = after_stat(..count.. / sum(..count..))), abaCellPlots, binwidth = 0.01) +
  labs(x = "species error, normalized", y = NULL) +
ggplot() +
  geom_histogram(aes(x = heightErrorMae / treesMatched, y = after_stat(..count.. / sum(..count..))), abaCellPlots, binwidth = 1) +
  labs(x = "tree height MAE, m normalized", y = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  coord_cartesian(ylim = c(0, 0.5)) &
  scale_y_continuous(labels = scales::percent)


# cell and match exploration
tibble(plots = nrow(plotHeights), 
       plot1 = length(unique(abaCellPlots$plot1)), shannon1 = vegan::diversity(abaCellPlots$plot1), # 96% of available plots matched 
       plot2 = length(unique(abaCellPlots$plot2)), shannon2 = vegan::diversity(abaCellPlots$plot2), 
       plot3 = length(unique(abaCellPlots$plot3)), shannon3 = vegan::diversity(abaCellPlots$plot3))

ggplot() +
  geom_histogram(aes(x = distance1 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
  labs(x = "normalized kNN distance\nto most similar plot", y = "probability", fill = "trees\nmatched") +
ggplot() +
  geom_histogram(aes(x = distance2 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
  labs(x = "normalized kNN distance\nto second plot", y = NULL, fill = "trees\nmatched") +
ggplot() +
  geom_histogram(aes(x = distance3 / treesMatched, y = after_stat(..count.. / sum(..count..)), fill = treesMatched, group = treesMatched), abaCellPlots, binwidth = 0.2) +
  labs(x = "normalized kNN distance\nto third plot", y = NULL, fill = "trees\nmatched") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect") &
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.4)) &
  scale_y_continuous(labels = scales::percent)

ggplot() +
  geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..)), fill = "cruise plots,\nunweighted by strata"), plotHeights, binwidth = 1) +
  geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..)), fill = "LiDAR"), abaCells, alpha = 0.5, binwidth = 1) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = bquote("trees 20 x 20 m grid cell"^-1), y = "probability", fill = NULL) +
  scale_y_continuous(labels = scales::percent)

left_join(abaCells %>% mutate(nClass = if_else(n < 16, n, 16)) %>% group_by(nClass) %>% summarize(cells = n()) %>% mutate(pctCells = 100 * cells / sum(cells)),
          plotHeights %>% mutate(nClass = if_else(n < 16, n, 16)) %>% group_by(nClass) %>% summarize(plots = n()) %>% mutate(pctPlots = 100 * plots / sum(plots)),
          by = "nClass")

plots2021 = left_join(trees2021organon, plots2016, by = c("plot")) %>%
  # TODO: adjust fixed radius tree count? 1/200 acre fixed radius plots = 217.6 ft² = 20.215702 m² => expansion to 19.7866 * TreeCount
  select(-PointID, -Dia1, -Ht1, -BHAge, -CrownRatio, -DefectMeasured, -STAND, -elevation, -slope, -aspect, -topographicShelterIndex, -x, -y, -heightDiameterRatio, -imputedBasalArea, -treeBasalAreaPerHectare, -treeBasalAreaPerHectareApprox, -plotsInStand, -standBasalAreaApprox, -measurePlotsInStand, -meanTreesPerBafPlot) %>%
  rename(plot = PlotID) %>%
  arrange(plot, desc(imputedHeight)) %>% 
  group_by(plot) %>%
  mutate(TreeID = row_number(), # renumber trees on plot in descending order of height
         minTreeID = 1, # add columns for cross join on tree ID; TODO: no longer needed, can be removed
         maxTreeID = n()) %>% 
  ungroup()
#plots2021 %>% group_by(speciesGroup) %>% summarize(n = sum(TreeCount), minHt = min(imputedHeight, na.rm = TRUE), maxHt = max(imputedHeight, na.rm = TRUE), isNAdbh = sum(TreeCount * is.na(DBH)), isNAheight = sum(TreeCount * is.na(imputedHeight)))

ggplot() +
  geom_histogram(aes(x = TotalHt, y = after_stat(..count.. / sum(..count..)), weight = TreeCount, fill = "plots"), plots2021, binwidth = 1) +
  geom_histogram(aes(x = height, y = after_stat(..count.. / sum(..count..)), fill = "cells"), trees2021lidar, alpha = 0.5, binwidth = 1) +
  coord_cartesian(xlim = c(0, 90)) +
  labs(x = "tree height, m", y = "probability", fill = NULL) +
  scale_y_continuous(labels = scales::percent)

trees2016 %>% group_by(PlotID) %>% summarize(hasFixedRadiusPlot = any(SamplingMethod == "FIX"), .groups = "drop") %>% group_by(hasFixedRadiusPlot) %>% summarize(plots = n()) # 1924 of 18363 = 10.5% of plots have trees on nested fixed radius 
#trees2016 %>% filter(TreeCount > 1, TotalHt > 4.9) %>% group_by(Species) %>% summarize(trees = sum(TreeCount)) # 973 trees with identical DBH, TotalHt, and TreeCount > 1

ggplot() +
  geom_histogram(aes(x = nNonLidarTrees, y = after_stat(..count.. / sum(..count..))), trees2016 %>% filter(imputedHeight < 4.92) %>% group_by(PlotID) %>% summarize(nNonLidarTrees = sum(if_else(SamplingMethod == "FIX", 19.7866, 1) * TreeCount), .groups = "drop"), binwidth = 20) +
  labs(x = "non-LiDAR trees per 20 m ABA grid cell", y = "probability") +
  scale_y_continuous(labels = scales::percent)

# plot exploration
ggplot() +
  geom_histogram(aes(x = abaGrid$sizeHa * liveExpansionFactor, y = after_stat(..count.. / sum(..count..)), fill = factor(if_else(species %in% c("DF", "RA", "WH", "BM", "OM", "RC"), species, "other"), levels = c("DF", "RA", "WH", "BM", "OM", "RC", "other"))), trees2021organon, binwidth = 0.1) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(x = bquote("individual trees plot"^-1), y = "probability", fill = NULL) +
  scale_fill_manual(breaks = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), limits = c("DF", "RA", "WH", "BM", "OM", "RC", "other"), values = c("forestgreen", "red2", "blue2", "green3", "mediumorchid1", "firebrick", "grey65")) +
  scale_y_continuous(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), labels = scales::percent, trans = scales::pseudo_log_trans(sigma = 0.01))

ggplot() +
  geom_histogram(aes(x = plot1, y = after_stat(..count.. / sum(..count..)), fill = factor(treesMatched), group = treesMatched), abaCellPlots, binwidth = 100) +
  labs(x = "plot ID", y = "probability", fill = "trees\nmatched") +
  scale_y_continuous(labels = scales::percent)


# create plot data from 2016 measurements rather than Organon predictions of 2021 growth
if (use2016trees)
{
  load("trees/height-diameter/data/ACMA3 preferred models.Rdata")
  load("trees/height-diameter/data/ALRU2 preferred models.Rdata")
  load("trees/height-diameter/data/other preferred models.Rdata")
  load("trees/height-diameter/data/PSME preferred models.Rdata")
  load("trees/height-diameter/data/THPL preferred models.Rdata")
  load("trees/height-diameter/data/TSHE preferred models.Rdata")
  load("trees/height-diameter/data/UMCA preferred models.Rdata")
  rm(acmaDiameterFromHeightPreferred, alruDiameterFromHeightPreferred, otherDiameterFromHeightPreferred, psmeDiameterFromHeightPreferred, thplDiameterFromHeightPreferred, tsheDiameterFromHeightPreferred, umcaDiameterFromHeightPreferred)
  
  impute_height = function(treeMeasurements)
  {
    # can't pipe to case_match() outside of a mutate statement (dplyr 1.1.3)
    return(case_match(treeMeasurements$speciesGroup,
                      "DF" ~ predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh, treeMeasurements), # GAM BA+L physio?
                      "RA" ~ predict(alruHeightFromDiameterPreferred$sharmaPartonBalPhysio, treeMeasurements),
                      "WH" ~ predict(tsheHeightFromDiameterPreferred$gamBalRelDbh, treeMeasurements),
                      "BM" ~ predict(acmaHeightFromDiameterPreferred$sharmaParton, treeMeasurements),
                      "OM" ~ predict(umcaHeightFromDiameterPreferred$sharmaPartonPhysio, treeMeasurements),
                      "RC" ~ predict(thplHeightFromDiameterPreferred$sharmaPartonPhysio, treeMeasurements),
                      "other" ~ predict(otherHeightFromDiameterPreferred$gamBal, treeMeasurements)))
  }
  
  plots2021 = trees2016 # from height-diameter/setup.R
    filter(PlotType == "IP", is.na(elevation) == FALSE) %>% # exclude dropped plots (IB, CB), count plots, and measure plots missing coordinates as they also lack ABA statistics for kNN
    mutate(imputedHeight = if_else(is.na(TotalHt) & is.na(Ht2), impute_height(.), if_else(is.na(TotalHt), Ht2, TotalHt))) %>% # several seconds due to GAMs
    # Maybe adjust fixed radius tree count? 1/200 acre fixed radius plots = 217.6 ft² = 20.215702 m² => expansion to 19.7866 * TreeCount
    select(-PointID, -Dia1, -Ht1, -BHAge, -CrownRatio, -DefectMeasured, -STAND, -elevation, -slope, -aspect, -topographicShelterIndex, -x, -y, -heightDiameterRatio, -imputedBasalArea, -treeBasalAreaPerHectare, -treeBasalAreaPerHectareApprox, -plotsInStand, -standBasalAreaApprox, -measurePlotsInStand, -meanTreesPerBafPlot) %>%
    rename(plot = PlotID) %>%
    arrange(plot, desc(imputedHeight)) %>% 
    group_by(plot) %>%
    mutate(TreeID = row_number(), # renumber trees on plot in descending order of height
           minTreeID = 1, # add columns for cross join on tree ID; TODO: no longer needed, can be removed
           maxTreeID = n()) %>% 
    ungroup()
}


splits = vfold_cv(plotMetrics2021, v = 2, repeats = 2)
split = splits[1, ]

getPlotError = function(split)
{
  # available predictors in LiDAR standard metrics
  # elevation
  # slope
  # aspect
  # topographicShelterIndex
  # zmax                    zmean                   zsd                     zskew                   zkurt                  
  # pzabovezmean            pzabove2                zq5                     zq10                    zq15                   
  # zq20                    zq25                    zq30                    zq35                    zq40                   
  # zq45                    zq50                    zq55                    zq60                    zq65                   
  # zq70                    zq75                    zq80                    zq85                    zq90                   
  # zq95                    zpcum1                  zpcum2                  zpcum3                  zpcum4                 
  # zpcum5                  zpcum6                  zpcum7                  zpcum8                  zpcum9                 
  # itot                    imax                    imean                   isd                     iskew                  
  # ikurt                   ipground                ipcumzq10               ipcumzq30               ipcumzq50              
  # ipcumzq70               ipcumzq90               p1th                    p2th                    p3th                   
  # p4th                    p5th                    pground
  trainingPlots = analysis(split$splits[[1]])
  validationPlots = assessment(split$splits[[1]])
  #plotNeighbors = get.knnx(trainingPlots %>% select(-plot), validationPlots %>% select(-plot), k = 10, algorithm = "kd_tree")
  #plotNeighbors = get.knnx(trainingPlots %>% select(starts_with("zq"), -plot), validationPlots %>% select(starts_with("zq"), -plot), k = 10, algorithm = "kd_tree")
  #plotNeighbors = get.knnx(trainingPlots %>% select(zmean), validationPlots %>% select(zmean), k = 10, algorithm = "kd_tree")
  trainingPlotNormalization = trainingPlots %>% summarize(across(starts_with("zq"), .fns = list(mean = mean, sd = sd), .names = "{col}{fn}"))
  trainingPlotsNormalized = trainingPlots %>% select(plot, starts_with("zq")) %>%
    mutate(zq5 = (zq5 - trainingPlotNormalization$zq5mean) / trainingPlotNormalization$zq5sd, # could probably use across() but quality of support seems uncertan, https://community.rstudio.com/t/normalizing-with-group-and-overall-means-using-dplyr/112956/2
           zq10 = (zq10 - trainingPlotNormalization$zq10mean) / trainingPlotNormalization$zq10sd,
           zq15 = (zq15 - trainingPlotNormalization$zq15mean) / trainingPlotNormalization$zq15sd,
           zq20 = (zq20 - trainingPlotNormalization$zq20mean) / trainingPlotNormalization$zq20sd,
           zq25 = (zq25 - trainingPlotNormalization$zq25mean) / trainingPlotNormalization$zq25sd,
           zq30 = (zq30 - trainingPlotNormalization$zq30mean) / trainingPlotNormalization$zq30sd,
           zq35 = (zq35 - trainingPlotNormalization$zq35mean) / trainingPlotNormalization$zq35sd,
           zq40 = (zq40 - trainingPlotNormalization$zq40mean) / trainingPlotNormalization$zq40sd,
           zq45 = (zq45 - trainingPlotNormalization$zq45mean) / trainingPlotNormalization$zq45sd,
           zq50 = (zq50 - trainingPlotNormalization$zq50mean) / trainingPlotNormalization$zq50sd,
           zq55 = (zq55 - trainingPlotNormalization$zq55mean) / trainingPlotNormalization$zq55sd,
           zq60 = (zq60 - trainingPlotNormalization$zq60mean) / trainingPlotNormalization$zq60sd,
           zq65 = (zq65 - trainingPlotNormalization$zq65mean) / trainingPlotNormalization$zq65sd,
           zq70 = (zq70 - trainingPlotNormalization$zq70mean) / trainingPlotNormalization$zq70sd,
           zq75 = (zq75 - trainingPlotNormalization$zq75mean) / trainingPlotNormalization$zq75sd,
           zq80 = (zq80 - trainingPlotNormalization$zq80mean) / trainingPlotNormalization$zq80sd,
           zq85 = (zq85 - trainingPlotNormalization$zq85mean) / trainingPlotNormalization$zq85sd,
           zq90 = (zq90 - trainingPlotNormalization$zq90mean) / trainingPlotNormalization$zq90sd,
           zq95 = (zq95 - trainingPlotNormalization$zq95mean) / trainingPlotNormalization$zq95sd)
  validationPlotsNormalized = validationPlots %>% select(plot, starts_with("zq")) %>%
    mutate(zq5 = (zq5 - trainingPlotNormalization$zq5mean) / trainingPlotNormalization$zq5sd,
           zq10 = (zq10 - trainingPlotNormalization$zq10mean) / trainingPlotNormalization$zq10sd,
           zq15 = (zq15 - trainingPlotNormalization$zq15mean) / trainingPlotNormalization$zq15sd,
           zq20 = (zq20 - trainingPlotNormalization$zq20mean) / trainingPlotNormalization$zq20sd,
           zq25 = (zq25 - trainingPlotNormalization$zq25mean) / trainingPlotNormalization$zq25sd,
           zq30 = (zq30 - trainingPlotNormalization$zq30mean) / trainingPlotNormalization$zq30sd,
           zq35 = (zq35 - trainingPlotNormalization$zq35mean) / trainingPlotNormalization$zq35sd,
           zq40 = (zq40 - trainingPlotNormalization$zq40mean) / trainingPlotNormalization$zq40sd,
           zq45 = (zq45 - trainingPlotNormalization$zq45mean) / trainingPlotNormalization$zq45sd,
           zq50 = (zq50 - trainingPlotNormalization$zq50mean) / trainingPlotNormalization$zq50sd,
           zq55 = (zq55 - trainingPlotNormalization$zq55mean) / trainingPlotNormalization$zq55sd,
           zq60 = (zq60 - trainingPlotNormalization$zq60mean) / trainingPlotNormalization$zq60sd,
           zq65 = (zq65 - trainingPlotNormalization$zq65mean) / trainingPlotNormalization$zq65sd,
           zq70 = (zq70 - trainingPlotNormalization$zq70mean) / trainingPlotNormalization$zq70sd,
           zq75 = (zq75 - trainingPlotNormalization$zq75mean) / trainingPlotNormalization$zq75sd,
           zq80 = (zq80 - trainingPlotNormalization$zq80mean) / trainingPlotNormalization$zq80sd,
           zq85 = (zq85 - trainingPlotNormalization$zq85mean) / trainingPlotNormalization$zq85sd,
           zq90 = (zq90 - trainingPlotNormalization$zq90mean) / trainingPlotNormalization$zq90sd,
           zq95 = (zq95 - trainingPlotNormalization$zq95mean) / trainingPlotNormalization$zq95sd)
  plotNeighbors = get.knnx(trainingPlotsNormalized %>% select(-plot), validationPlotsNormalized %>% select(-plot), k = 10, algorithm = "kd_tree")
  
  validationPlots %<>% mutate(neighbor1 = trainingPlots$plot[plotNeighbors$nn.index[, 1]],
                              neighbor2 = trainingPlots$plot[plotNeighbors$nn.index[, 2]],
                              neighbor3 = trainingPlots$plot[plotNeighbors$nn.index[, 3]])
  validationPlotMeasurements = left_join(plots2021 %>% filter(plot %in% validationPlots$plot),
                                         validationPlots %>% select(plot, starts_with("neighbor")),
                                         by = c("plot"))

  # validation plot: 9
  # nearest training plot: 13044
  #bind_rows(validationPlots %>% filter(plot == 9),
  #          trainingPlots %>% filter(plot == 13044))
  #bind_rows(validationPlotMeasurements %>% filter(plot == 9) %>% select(StandID, plot, TreeID, Species, TotalHt, imputedHeight, DBH),
  #          plotMeasurements %>% filter(plot == validationPlots$neighbor1[1]) %>% select(StandID, plot, TreeID, Species, TotalHt, imputedHeight, DBH))
  plotComparison1 = full_join(validationPlotMeasurements,
                              plots2021, # substantially faster if plotMeasurements isn't grouped
                              by = join_by(x$neighbor1 == y$plot, overlaps(x$minTreeID, x$maxTreeID, y$minTreeID, y$maxTreeID)), suffix = c("", "1")) %>%
    filter((TreeID == TreeID1) | is.na(TreeID1) | ((TreeID == maxTreeID) & (TreeID1 > maxTreeID)) | ((TreeID1 == maxTreeID1) & (TreeID > maxTreeID1))) %>%
    mutate(TreeCount = if_else(TreeID < TreeID1, NA, TreeCount),
           TreeCount1 = if_else(TreeID1 < TreeID, NA, TreeCount1))
  #print(plotComparison1 %>% filter(is.na(TreeID1)) %>% select(StandID, plot, maxTreeID, TreeID, TreeCount, Species, DBH, imputedHeight, StandID1, neighbor1, TreeID1, maxTreeID1, TreeCount1, Species1, DBH1, imputedHeight1), n = 10)
  #length(unique(plotComparison1$plot))
  #plotMeasurements %>% filter(plot == 18276) %>% select(StandID, plot, minTreeID, maxTreeID, TreeID, TreeCount, Species, DBH, imputedHeight)
  plotComparison2 = full_join(validationPlotMeasurements,
                              plots2021,
                              by = join_by(x$neighbor2 == y$plot, overlaps(x$minTreeID, x$maxTreeID, y$minTreeID, y$maxTreeID)), suffix = c("", "2")) %>%
    filter((TreeID == TreeID2) | ((TreeID == maxTreeID) & (TreeID2 > maxTreeID)) | ((TreeID2 == maxTreeID2) & (TreeID > maxTreeID2))) %>%
    mutate(TreeCount = if_else(TreeID < TreeID2, NA, TreeCount),
           TreeCount2 = if_else(TreeID2 < TreeID, NA, TreeCount2))
  #print(plotComparison2 %>% select(StandID, plot, TreeID, TreeCount, Species, DBH, imputedHeight, StandID2, neighbor2, TreeID2, TreeCount2, Species2, DBH2, imputedHeight2), n = 50)
  plotComparison3 = full_join(validationPlotMeasurements,
                              plots2021,
                              by = join_by(x$neighbor3 == y$plot, overlaps(x$minTreeID, x$maxTreeID, y$minTreeID, y$maxTreeID)), suffix = c("", "3")) %>%
    filter((TreeID == TreeID3) | ((TreeID == maxTreeID) & (TreeID3 > maxTreeID)) | ((TreeID3 == maxTreeID3) & (TreeID > maxTreeID3))) %>%
    mutate(TreeCount = if_else(TreeID < TreeID3, NA, TreeCount),
           TreeCount3 = if_else(TreeID3 < TreeID, NA, TreeCount3))
  
  plotDifferences = left_join(left_join(plotComparison1 %>% group_by(plot) %>% summarize(trees = sum(TreeCount, na.rm = TRUE), trees1 = sum(TreeCount1, na.rm = TRUE)),
                                        plotComparison2 %>% group_by(plot) %>% summarize(trees2 = sum(TreeCount2, na.rm = TRUE)),
                                        by = c("plot")),
                              plotComparison3 %>% group_by(plot) %>% summarize(trees3 = sum(TreeCount3, na.rm = TRUE)),
                              by = c("plot"))
  #cor(plotDifferences %>% select(-plot))
  #plotDifferences %>% filter(is.na(trees2))
  plotDifferences %>% mutate(treesDiff1 = trees1 - trees,
                             treesDiff2 = 1/2 * (trees1 + trees2) - trees,
                             treesDiff3 = 1/3 * (trees1 + trees2 + trees3) - trees) %>% 
    summarise(treesMin = min(trees), treesMax = max(trees), treesMean = mean(trees), 
              treesMae = mean(abs(trees - mean(trees))), treesSd = sd(trees), 
              treesMae1 = mean(abs(treesDiff1)), treesRmse1 = sqrt(mean(treesDiff1^2)),
              treesMae2 = mean(abs(treesDiff2)), treesRmse2 = sqrt(mean(treesDiff2^2)),
              treesMae3 = mean(abs(treesDiff3)), treesRmse3 = sqrt(mean(treesDiff3^2)))
  
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 40, yend = 40), color = "grey70", linetype = "longdash") +
    geom_bin2d(aes(x = trees, y = trees1), plotDifferences, binwidth = c(1, 1)) +
    labs(x = "measured tree count", y = "tree count estimate, k = 1", fill = "plots") +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 40, yend = 40), color = "grey70", linetype = "longdash") +
    geom_bin2d(aes(x = trees, y = trees2), plotDifferences, binwidth = c(1, 1)) +
    labs(x = "measured tree count", y = "tree count estimate, k = 2", fill = "plots") +
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 40, yend = 40), color = "grey70", linetype = "longdash") +
    geom_bin2d(aes(x = trees, y = trees3), plotDifferences, binwidth = c(1, 1)) +
    labs(x = "measured tree count", y = "tree count estimate, k = 3", fill = "plots") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_fill_viridis_c(limits = c(1, 200), trans = "log10")
}

split = splits[1, ]
getPlotError(split)

splits$splits[[1]]$id
splitsAndFits = vfold_cv(plotMetrics2021, v = 2, repeats = 1) %>% mutate(fit = future_map(splits, fitFunction))


## correlations between predictors and plot measurements
plotStatistics = plots2021 %>% group_by(plot) %>% 
  summarize(trees = max(TreeID), treeCount = sum(TreeCount), tph = sum(measureTreeTphContribution), 
            h1 = imputedHeight[1], h2 = if_else(n() > 1, imputedHeight[2], NA_real_), h3 = if_else(n() > 2, imputedHeight[3], NA_real_), .groups = "drop")
plotCorrelations = cor(left_join(plotMetrics2021, plotStatistics, by = "plot"), use = "complete.obs")
ggplot() +
  geom_tile(aes(x = var1, y = var2, fill = value), as_tibble(plotCorrelations) %>% mutate(var1 = rownames(plotCorrelations)) %>% pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
              mutate(var1 = as_factor(var1), var2 = as_factor(var2))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_scico(palette = "cork", limits = c(-1, 1)) +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## specialized classification and regression random forests for tree count and tallest tree height
library(caret)
library(ranger)
repeatedCrossValidation = trainControl(method = "repeatedcv", number = 2, repeats = 5, verboseIter = TRUE)

treeCountData = left_join(plotMetrics2021, plotStatistics %>% select(plot, treeCount), by = "plot") %>% select(-plot)
#treeCountForest = ranger(treeCount ~ ., treeCountData, classification = TRUE, num.threads = 12) # 32 s
treeCountForest = train(treeCount ~ ., data = treeCountData %>% sample_n(1000) %>% mutate(treeCount = as_factor(treeCount)), method = "ranger", trControl = repeatedCrossValidation, 
                        tuneGrid = expand.grid(mtry = floor(sqrt(dim(treeCountData)[2])),
                                               splitrule = "gini",
                                               min.node.size = 1))

treeCountData = left_join(plotMetrics2021, plotStatistics %>% select(plot, h1), by = "plot") %>% select(-plot)
h1forest = train(h1 ~ ., data = treeCountData %>% sample_n(1000), method = "ranger", trControl = repeatedCrossValidation, 
                 tuneGrid = expand.grid(mtry = floor(sqrt(dim(treeCountData)[2])),
                                        splitrule = "variance",
                                        min.node.size = 1))


## random forest computational intractability of plot matching
# Initial progress message is slow to appear on large fits and overestimates time by an order of magnitude.
# n       threads   fit time     loopback accuracy, %
#   100    8          0.1 s       0.9 - perfect recall of known plots
#   200    8          0.5 s       1.9
#   500    8          0.5 m       4.8
#  1000    8          2.8 m       9.6
# 10365   12        > 17 h
#library(caret)
#library(ranger)
#fitStart = Sys.time()
#plotForest = ranger(plot ~ ., plotMetrics, classification = TRUE, importance = "impurity", num.threads = 12) # 12 threads keep 14 cores at 100%
#Sys.time() - fitStart
#save(plotForest, file = file.path(getwd(), "trees/area based/plotRandomForestFit.Rdata"))

#sort(importance(plotForest), decreasing = TRUE)
#predictStart = Sys.time()
#plotPredictions = predict(plotForest, plotMetrics)
#Sys.time() - predictStart
#sum(plotPredictions$predictions == plotMetrics$plot)

#plotForest = train(plot ~ ., data = plotMetrics, method = "ranger", trControl = repeatedCrossValidation, 
#                  tuneGrid = expand.grid(mtry = ceiling(sqrt(dim(plotMetrics)[2] - 1)),
#                                         splitrule = 'gini',
#                                         min.node.size = 2))

#fitStart = Sys.time()
#plotKknn = train(plot ~ ., data = plotMetrics %>% sample_n(2000), method = "kknn", tuneGrid = expand.grid(kmax = c(13, 15, 17), distance = c(2, 3, 4), kernel = "optimal")) # 18s @ k = 2, r = 5, n = 2000, 6 element grid
#Sys.time() - fitStart
#plotKnn = train(plot ~ ., data = plotMetrics %>% sample_n(2000), method = "knn", tuneGrid = data.frame(k = c(9, 11, 13, 15, 17, 19))) # 11s @ k = 2, r = 5, n = 2000, 6 k values

## 2009-2021 treetop matching
library(dplyr)
library(FNN)
library(ggplot2)
library(patchwork)
library(sf)

theme_set(theme_bw() + theme(axis.line = element_line(linewidth = 0.3), 
                             legend.title = element_text(size = 10),
                             panel.border = element_blank()))

treetops2009 = st_read(file.path(getwd(), "GIS/DOGAMI/2009 OLC South Coast/treetops DSM ring/treetops clipped.gpkg"), quiet = TRUE) # slow
treetops2021 = st_read(file.path(getwd(), "GIS/DOGAMI/2021 OLC Coos County/treetops DSM ring/treetops clipped.gpkg"), quiet = TRUE) # slow

treetops2021neighbors2009 = get.knnx(st_coordinates(treetops2009)[, c("X", "Y")], st_coordinates(treetops2021)[, c("X", "Y")], k = 3) # ~15 s

treetops2021match = st_drop_geometry(treetops2021) %>%
  mutate(x = st_coordinates(treetops2021)[, "X"], y = st_coordinates(treetops2021)[, "Y"], elevation = st_coordinates(treetops2021)[, "Z"],
         height = 0.3048 * height,
         distance1 = 0.3048 * treetops2021neighbors2009$nn.dist[, 1], height2009_1 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 1]],
         distance2 = 0.3048 * treetops2021neighbors2009$nn.dist[, 2], height2009_2 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 2]],
         distance3 = 0.3048 * treetops2021neighbors2009$nn.dist[, 3], height2009_3 = 0.3048 * treetops2009$height[treetops2021neighbors2009$nn.index[, 3]],
         distance = if_else(height > height2009_1, distance1, if_else(distance2 < 0.3048 * 5, distance2, NA_real_)),
         height2009 = if_else(height > height2009_1, height2009_1, if_else(distance2 < 0.3048 * 5, height2009_2, NA_real_)),
         heightIncrement = height - height2009)

treetops2021match %>% filter(distance < pmin(0.1 * height2009, 0.3048 * 5)) %>%
  mutate(classification = if_else(heightIncrement > 24, "> 2 m/year", if_else(heightIncrement > 12, "1–2 m/year", if_else(heightIncrement < -0.5 * height2009, "replaced", if_else(heightIncrement < -12 * 0.01, "broken", "-0.1–1 m/year"))))) %>% 
  group_by(classification) %>% summarize(n = n()) %>% mutate(pct = 100 * n / sum(n), pctAll = 100 * n / nrow(treetops2021)) %>%
  arrange(desc(n))

ggplot(treetops2021match) +
  geom_histogram(aes(x = heightIncrement, y = after_stat(..count.. / sum(..count..))), binwidth = 1) +
  labs(x = "2009–2021 height growth, m", y = "trees") +
  scale_y_continuous(labels = scales::percent)

ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 100, yend = 100), color = "grey70", linetype = "longdash") +
  geom_bin_2d(aes(x = height2009, y = height), treetops2021match %>% filter(height >= (height2009 - 3.2808 * 5), height <= (height2009 + 3.2808 * 30), distance < pmin(0.1 * height2009, 0.3048 * 5)), binwidth = c(1, 1)) +
  labs(x = "matched tree height in 2009, m", y = "tree height in 2021, m", fill = "trees\nmatched") +
  scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")
ggplot() +
  geom_bin_2d(aes(x = height, y = 1/12 * (height - height2009)), treetops2021match %>% filter(height >= (height2009 - 3.2808 * 5), height <= (height2009 + 3.2808 * 30), distance < pmin(0.1 * height2009, 0.3048 * 5)), binwidth = c(1, 0.1)) +
  labs(x = "matched tree height in 2021, m", y = bquote("mean annual height growth, m year"^-1), fill = "trees\nmatched") +
  scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")


ggplot(treetops2021match %>% filter(distance1 < 30)) +
  geom_histogram(aes(x = 0.3048 * distance1, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
  labs(x = "distance to nearest treetop, m", y = "") +
ggplot(treetops2021match %>% filter(distance2 < 30)) +
  geom_histogram(aes(x = 0.3048 * distance2, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
  labs(x = "distance to second nearest treetop, m", y = NULL) +
ggplot(treetops2021match %>% filter(distance3 < 30)) +
  geom_histogram(aes(x = 0.3048 * distance3, y = after_stat(..count.. / sum(..count..))), binwidth = 0.3048 * 1.5) +
  labs(x = "distance to third nearest treetop, m", y = NULL) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout() &
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 10)) +
  scale_y_continuous(labels = scales::percent)

## Organon growth intervals
#library(readr)
#organon2021csv = read_csv(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.csv"), col_types = cols(stand = "i", plot = "i", species = "c", tag = "i", year = "i", standAge = "i", .default = "d"))
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)

organon2021 = read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE) %>%
  filter(year <= 2021)

ggplot() + 
  geom_bin_2d(aes(x = height2021, y = 1/5 * (height2021 - height2016), weight = liveExpansionFactor2021), organon2021 %>% pivot_wider(id_cols = c(stand, plot, tag, species), names_from = year, names_sep = "", values_from = c(height, liveExpansionFactor)), binwidth = c(1, 0.1)) +
  coord_cartesian(ylim = c(-1.5, 5.75)) +
  labs(x = "predicted tree height in 2021, m", y = bquote("mean annual height growth, m year"^-1), fill = bquote("trees ha"^-1)) +
  scale_fill_viridis_c(breaks = c(1, 10, 100, 1000, 10000), trans = "log10")
