#install.packages(c("arrow", "collapse", "FNN"))
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

abaOptions = tibble(folds = 2, repetitions = 2,
                    coniferDeciduousWeight = 5,
                    includeInvestigatory = FALSE, 
                    recalcAbaCellOccupancy = FALSE)

get_aba_cell_plots = function(validationCells, validationCellsScaled, trainingPlotHeightsScaled, treeMatchBound = 8, lidarMetrics = c("pGround", "zQ10", "zQ20", "zQ30"))
{
  if ((treeMatchBound < 1) | (treeMatchBound > 15))
  {
    stop(paste0("treeMatchBound of ", treeMatchBound, " is not supported."))
  }
  
  # treeMatchBound constrains cells with nCell = 1..treeMatchBound treetops to be matched to plots with nPlot >= nCell trees
  # When forming kNN inputs, the filter on the matchable plot list is therefore always nPlot >= nCell.
  # When forming kNN queries, the cell filters are n == 1, n == 2, ... n = treeMatchBound - 1, n >= treeMatchBound.
  # A set of matched cell slices results, { abaCells1, abaCells2, ... abaCells$(treeMatchBound - 1), abaCells$(treeMatchBound) },
  # which then needs to be row bound back together. Formation of each slice therefore comes from 
  #   cells %>% filter(n == 1..treeMatchBound - 1)
  # and then
  #   cells %>% filter(n >= treeMatchBound)
  # To limit code repetition, the filter clause is evaluated dynamically using rlang, though if { } else { } could also be used.
  cellFilter1 = if (treeMatchBound > 1) { "n == 1" } else { "n >= 1" } # could optimize out n >= 1 check
  plotHeights1 = trainingPlotHeightsScaled %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1)
  abaCells1 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter1)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1)
  cellPlotIndices1 = knnx.index(plotHeights1 %>% select(-plot), abaCells1 %>% select(-abaGridX, -abaGridY), k = 3)
  abaCells1 %<>% mutate(treesMatched = 1, plot1 = plotHeights1$plot[cellPlotIndices1[, 1]], plot2 = plotHeights1$plot[cellPlotIndices1[, 2]], plot3 = plotHeights1$plot[cellPlotIndices1[, 3]])

  if (treeMatchBound >= 2)
  {
    cellFilter2 = if (treeMatchBound > 2) { "n == 2" } else { "n >= 2" }
    plotHeights2 = trainingPlotHeightsScaled %>% filter(n >= 2) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2)
    abaCells2 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter2)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2)
    cellPlotIndices2 = knnx.index(plotHeights2 %>% select(-plot), abaCells2 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells2 %<>% mutate(treesMatched = 2, plot1 = plotHeights2$plot[cellPlotIndices2[, 1]], plot2 = plotHeights2$plot[cellPlotIndices2[, 2]], plot3 = plotHeights2$plot[cellPlotIndices2[, 3]])
  }
  if (treeMatchBound >= 3)
  {
    cellFilter3 = if (treeMatchBound > 3) { "n == 3" } else { "n >= 3" }
    plotHeights3 = trainingPlotHeightsScaled %>% filter(n >= 3) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3)
    abaCells3 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter3)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3)
    cellPlotIndices3 = knnx.index(plotHeights3 %>% select(-plot), abaCells3 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells3 %<>% mutate(treesMatched = 3, plot1 = plotHeights3$plot[cellPlotIndices3[, 1]], plot2 = plotHeights3$plot[cellPlotIndices3[, 2]], plot3 = plotHeights3$plot[cellPlotIndices3[, 3]])
  }
  if (treeMatchBound >= 4)
  {
    cellFilter4 = if (treeMatchBound > 4) { "n == 4" } else { "n >= 4" }
    plotHeights4 = trainingPlotHeightsScaled %>% filter(n >= 4) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4)
    abaCells4 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter4)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4)
    cellPlotIndices4 = knnx.index(plotHeights4 %>% select(-plot), abaCells4 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells4 %<>% mutate(treesMatched = 4, plot1 = plotHeights4$plot[cellPlotIndices4[, 1]], plot2 = plotHeights4$plot[cellPlotIndices4[, 2]], plot3 = plotHeights4$plot[cellPlotIndices4[, 3]])
  }  
  if (treeMatchBound >= 5)
  {
    cellFilter5 = if (treeMatchBound > 5) { "n == 5" } else { "n >= 5" }
    plotHeights5 = trainingPlotHeightsScaled %>% filter(n >= 5) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5)
    abaCells5 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter5)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5)
    cellPlotIndices5 = knnx.index(plotHeights5 %>% select(-plot), abaCells5 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells5 %<>% mutate(treesMatched = 5, plot1 = plotHeights5$plot[cellPlotIndices5[, 1]], plot2 = plotHeights5$plot[cellPlotIndices5[, 2]], plot3 = plotHeights5$plot[cellPlotIndices5[, 3]])
  }
  if (treeMatchBound >= 6)
  {
    cellFilter6 = if (treeMatchBound > 6) { "n == 6" } else { "n >= 6" }
    plotHeights6 = trainingPlotHeightsScaled %>% filter(n >= 6) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6)
    abaCells6 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter6)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6)
    cellPlotIndices6 = knnx.index(plotHeights6 %>% select(-plot), abaCells6 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells6 %<>% mutate(treesMatched = 6, plot1 = plotHeights6$plot[cellPlotIndices6[, 1]], plot2 = plotHeights6$plot[cellPlotIndices6[, 2]], plot3 = plotHeights6$plot[cellPlotIndices6[, 3]])
  }
  if (treeMatchBound >= 7)
  {
    cellFilter7 = if (treeMatchBound > 7) { "n == 7" } else { "n >= 7" }
    plotHeights7 = trainingPlotHeightsScaled %>% filter(n >= 7) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7)
    abaCells7 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter7)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7)
    cellPlotIndices7 = knnx.index(plotHeights7 %>% select(-plot), abaCells7 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells7 %<>% mutate(treesMatched = 7, plot1 = plotHeights7$plot[cellPlotIndices7[, 1]], plot2 = plotHeights7$plot[cellPlotIndices7[, 2]], plot3 = plotHeights7$plot[cellPlotIndices7[, 3]])
  }
  if (treeMatchBound >= 8)
  {
    cellFilter8 = if (treeMatchBound > 8) { "n == 8" } else { "n >= 8" }
    plotHeights8 = trainingPlotHeightsScaled %>% filter(n >= 8) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8)
    abaCells8 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter8)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8)
    cellPlotIndices8 = knnx.index(plotHeights8 %>% select(-plot), abaCells8 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells8 %<>% mutate(treesMatched = 8, plot1 = plotHeights8$plot[cellPlotIndices8[, 1]], plot2 = plotHeights8$plot[cellPlotIndices8[, 2]], plot3 = plotHeights8$plot[cellPlotIndices8[, 3]])
  }
  if (treeMatchBound >= 9)
  {
    cellFilter9 = if (treeMatchBound > 9) { "n == 9" } else { "n >= 9" }
    plotHeights9 = trainingPlotHeightsScaled %>% filter(n >= 9) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9)
    abaCells9 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter9)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9)
    cellPlotIndices9 = knnx.index(plotHeights9 %>% select(-plot), abaCells9 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells9 %<>% mutate(treesMatched = 9, plot1 = plotHeights9$plot[cellPlotIndices9[, 1]], plot2 = plotHeights9$plot[cellPlotIndices9[, 2]], plot3 = plotHeights9$plot[cellPlotIndices9[, 3]])
  }
  if (treeMatchBound >= 10)
  {
    cellFilter10 = if (treeMatchBound > 10) { "n == 10" } else { "n >= 10" }
    plotHeights10 = trainingPlotHeightsScaled %>% filter(n >= 10) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10)
    abaCells10 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter10)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10)
    cellPlotIndices10 = knnx.index(plotHeights10 %>% select(-plot), abaCells10 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells10 %<>% mutate(treesMatched = 10, plot1 = plotHeights10$plot[cellPlotIndices10[, 1]], plot2 = plotHeights10$plot[cellPlotIndices10[, 2]], plot3 = plotHeights10$plot[cellPlotIndices10[, 3]])
  }
  if (treeMatchBound >= 11)
  {
    cellFilter11 = if (treeMatchBound > 11) { "n == 11" } else { "n >= 11" }
    plotHeights11 = trainingPlotHeightsScaled %>% filter(n >= 11) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11)
    abaCells11 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter11)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11)
    cellPlotIndices11 = knnx.index(plotHeights11 %>% select(-plot), abaCells11 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells11 %<>% mutate(treesMatched = 11, plot1 = plotHeights11$plot[cellPlotIndices11[, 1]], plot2 = plotHeights11$plot[cellPlotIndices11[, 2]], plot3 = plotHeights11$plot[cellPlotIndices11[, 3]])
  }
  if (treeMatchBound >= 12)
  {
    cellFilter12 = if (treeMatchBound > 12) { "n == 12" } else { "n >= 12" }
    plotHeights12 = trainingPlotHeightsScaled %>% filter(n >= 12) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12)
    abaCells12 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter12)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12)
    cellPlotIndices12 = knnx.index(plotHeights12 %>% select(-plot), abaCells12 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells12 %<>% mutate(treesMatched = 12, plot1 = plotHeights12$plot[cellPlotIndices12[, 1]], plot2 = plotHeights12$plot[cellPlotIndices12[, 2]], plot3 = plotHeights12$plot[cellPlotIndices12[, 3]])
  }
  if (treeMatchBound >= 13)
  {
    cellFilter13 = if (treeMatchBound > 13) { "n == 13" } else { "n >= 13" }
    plotHeights13 = trainingPlotHeightsScaled %>% filter(n >= 13) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13)
    abaCells13 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter13)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13)
    cellPlotIndices13 = knnx.index(plotHeights13 %>% select(-plot), abaCells13 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells13 %<>% mutate(treesMatched = 13, plot1 = plotHeights13$plot[cellPlotIndices13[, 1]], plot2 = plotHeights13$plot[cellPlotIndices13[, 2]], plot3 = plotHeights13$plot[cellPlotIndices13[, 3]])
  }
  if (treeMatchBound >= 14)
  {
    cellFilter14 = if (treeMatchBound > 14) { "n == 14" } else { "n >= 14" }
    plotHeights14 = trainingPlotHeightsScaled %>% filter(n >= 14) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13, isConifer14, height14)
    abaCells14 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter14)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13, isConifer14, height14)
    cellPlotIndices14 = knnx.index(plotHeights14 %>% select(-plot), abaCells14 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells14 %<>% mutate(treesMatched = 14, plot1 = plotHeights14$plot[cellPlotIndices14[, 1]], plot2 = plotHeights14$plot[cellPlotIndices14[, 2]], plot3 = plotHeights14$plot[cellPlotIndices14[, 3]])
  }
  if (treeMatchBound >= 14)
  {
    cellFilter15 = if (treeMatchBound > 15) { "n == 15" } else { "n >= 15" }
    plotHeights15 = trainingPlotHeightsScaled %>% filter(n >= 15) %>% select(plot, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13, isConifer14, height14, isConifer15, height15)
    abaCells15 = validationCellsScaled %>% filter(!!rlang::parse_expr(cellFilter15)) %>% select(abaGridX, abaGridY, pGround, all_of(lidarMetrics), isConifer1, height1, isConifer2, height2, isConifer3, height3, isConifer4, height4, isConifer5, height5, isConifer6, height6, isConifer7, height7, isConifer8, height8, isConifer9, height9, isConifer10, height10, isConifer11, height11, isConifer12, height12, isConifer13, height13, isConifer14, height14, isConifer15, height15)
    cellPlotIndices15 = knnx.index(plotHeights15 %>% select(-plot), abaCells15 %>% select(-abaGridX, -abaGridY), k = 3)
    abaCells15 %<>% mutate(treesMatched = 15, plot1 = plotHeights15$plot[cellPlotIndices15[, 1]], plot2 = plotHeights15$plot[cellPlotIndices15[, 2]], plot3 = plotHeights15$plot[cellPlotIndices15[, 3]])
  }
  
  # form abaCellToPlots with an if-else stack as R's integer switch is limited and has somewhat cryptic syntax
  if (treeMatchBound == 1) { abaCellToPlots = abaCells1 }
  else if (treeMatchBound == 2) { abaCellToPlots = bind_rows(abaCells1, abaCells2) }
  else if (treeMatchBound == 3) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3) }
  else if (treeMatchBound == 4) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4) }
  else if (treeMatchBound == 5) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5) }
  else if (treeMatchBound == 6) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6) }
  else if (treeMatchBound == 7) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7) }
  else if (treeMatchBound == 8) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8) }
  else if (treeMatchBound == 9) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9) }
  else if (treeMatchBound == 10) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10) }
  else if (treeMatchBound == 11) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10, abaCells11) }
  else if (treeMatchBound == 12) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10, abaCells11, abaCells12) }
  else if (treeMatchBound == 13) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10, abaCells11, abaCells12, abaCells13) }
  else if (treeMatchBound == 14) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10, abaCells11, abaCells12, abaCells13, abaCells14) }
  else if (treeMatchBound == 15) { abaCellToPlots = bind_rows(abaCells1, abaCells2, abaCells3, abaCells4, abaCells5, abaCells6, abaCells7, abaCells8, abaCells9, abaCells10, abaCells11, abaCells12, abaCells13, abaCells14, abaCells15) }
  else { stop(paste0("Unhandled treeMatchBound = ", treeMatchBound, " in cell to plot list formation.")) }
  abaCellPlots = left_join(left_join(validationCells, 
                                     abaCellToPlots %>% select(abaGridX, abaGridY, treesMatched, plot1, plot2, plot3), 
                                     by = join_by(abaGridX, abaGridY)),
                           plotHeights %>% rename(plot1n = n) %>% rename_with(~paste0("plot1", .x), starts_with(c("species", "isConifer", "height"))) %>%
                             select(plot, starts_with("plot1")),
                           by = join_by(x$plot1 == y$plot)) %>%
    # following lines not expressible with across() (dplyr 1.1.4), pivoting longer and back to wider is cumbersome
    mutate(speciesError1matched = isConifer1 != plot1isConifer1 + if_else(treesMatched < 2, 0, isConifer2 != plot1isConifer2) + if_else(treesMatched < 3, 0, isConifer3 != plot1isConifer3) + if_else(treesMatched < 4, 0, isConifer4 != plot1isConifer4) + 
             if_else(treesMatched < 5, 0, isConifer5 != plot1isConifer5) + if_else(treesMatched < 6, 0, isConifer6 != plot1isConifer6) + if_else(treesMatched < 7, 0, isConifer7 != plot1isConifer7) + if_else(treesMatched < 8, 0, isConifer8 != plot1isConifer8) + if_else(treesMatched < 9, 0, isConifer9 != plot1isConifer9) + 
             if_else(treesMatched < 10, 0, isConifer10 != plot1isConifer10) + if_else(treesMatched < 11, 0, isConifer11 != plot1isConifer11) + if_else(treesMatched < 12, 0, isConifer12 != plot1isConifer12) + if_else(treesMatched < 13, 0, isConifer13 != plot1isConifer13) + if_else(treesMatched < 14, 0, isConifer14 != plot1isConifer14) + if_else(treesMatched < 15, 0, isConifer15 != plot1isConifer15),
           speciesError1unmatched = if_else((treesMatched >= 2) | is.na(isConifer2), 0, isConifer2 != plot1isConifer2) + if_else((treesMatched >= 3) | is.na(isConifer3), 0, isConifer3 != plot1isConifer3) + if_else((treesMatched >= 4) | is.na(isConifer4), 0, isConifer4 != plot1isConifer4) + 
             if_else((treesMatched >= 5) | is.na(isConifer5), 0, isConifer5 != plot1isConifer5) + if_else((treesMatched >= 6) | is.na(isConifer6), 0, isConifer6 != plot1isConifer6) + if_else((treesMatched >= 7) | is.na(isConifer7), 0, isConifer7 != plot1isConifer7) + if_else((treesMatched >= 8) | is.na(isConifer8), 0, isConifer8 != plot1isConifer8) + if_else((treesMatched >= 9) | is.na(isConifer9), 0, isConifer9 != plot1isConifer9) + 
             if_else((treesMatched >= 10) | is.na(isConifer10), 0, isConifer10 != plot1isConifer10) + if_else((treesMatched >= 11) | is.na(isConifer11), 0, isConifer11 != plot1isConifer11) + if_else((treesMatched >= 12) | is.na(isConifer12), 0, isConifer12 != plot1isConifer12) + if_else((treesMatched >= 13) | is.na(isConifer13), 0, isConifer13 != plot1isConifer13) + if_else((treesMatched >= 14) | is.na(isConifer14), 0, isConifer14 != plot1isConifer14) + if_else((treesMatched >= 15) | is.na(isConifer15), 0, isConifer15 != plot1isConifer5),
           heightError1matched = abs(height1 - plot1height1) + if_else(treesMatched < 2, 0, abs(height2 - plot1height2)) + if_else(treesMatched < 3, 0, abs(height3 - plot1height3)) + if_else(treesMatched < 4, 0, abs(height4 - plot1height4)) + 
             if_else(treesMatched < 5, 0, abs(height5 - plot1height5)) + if_else(treesMatched < 6, 0, abs(height6 - plot1height6)) + if_else(treesMatched < 7, 0, abs(height7 - plot1height7)) + if_else(treesMatched < 8, 0, abs(height8 - plot1height8)) + if_else(treesMatched < 9, 0, abs(height9 - plot1height9)) +
             if_else(treesMatched < 10, 0, abs(height10 - plot1height10)) + if_else(treesMatched < 11, 0, abs(height11 - plot1height11)) + if_else(treesMatched < 12, 0, abs(height12 - plot1height12)) + if_else(treesMatched < 13, 0, abs(height13 - plot1height13)) + if_else(treesMatched < 14, 0, abs(height14 - plot1height14)) + if_else(treesMatched < 15, 0, abs(height15 - plot1height15)),
           heightError1unmatched = if_else((treesMatched >= 2) | is.na(height2) | is.na(plot1height2), 0, abs(height2 - plot1height2)) + if_else((treesMatched >= 3) | is.na(height3) | is.na(plot1height3), 0, abs(height3 - plot1height3)) + if_else((treesMatched >= 4) | is.na(height4) | is.na(plot1height4), 0, abs(height4 - plot1height4)) + 
             if_else((treesMatched >= 5) | is.na(height5) | is.na(plot1height5), 0, abs(height5 - plot1height5)) + if_else((treesMatched >= 6) | is.na(height6) | is.na(plot1height6), 0, abs(height6 - plot1height6)) + if_else((treesMatched >= 7) | is.na(height7) | is.na(plot1height7), 0, abs(height7 - plot1height7)) + if_else((treesMatched >= 8) | is.na(height8) | is.na(plot1height8), 0, abs(height8 - plot1height8)) + if_else((treesMatched >= 9) | is.na(height9) | is.na(plot1height9), 0, abs(height9 - plot1height9)) + 
             if_else((treesMatched >= 10) | is.na(height10) | is.na(plot1height10), 0, abs(height10 - plot1height10)) + if_else((treesMatched >= 11) | is.na(height11) | is.na(plot1height11), 0, abs(height11 - plot1height11)) + if_else((treesMatched >= 12) | is.na(height12) | is.na(plot1height12), 0, abs(height12 - plot1height12)) + if_else((treesMatched >= 13) | is.na(height13) | is.na(plot1height13), 0, abs(height13 - plot1height13)) + if_else((treesMatched >= 14) | is.na(height14) | is.na(plot1height14), 0, abs(height14 - plot1height14)) + if_else((treesMatched >= 15) | is.na(height15) | is.na(plot1height15), 0, abs(height15 - plot1height15))) %>%
    # if present, exclude grid cells with partial LiDAR coverage: 433 edge cells in the total grid have treetop detections but NA LiDAR metrics
    filter(is.na(nPoints) == FALSE)
  #abaCellPlots %>% summarize(plot1 = sum(is.na(plot1)), plot2 = sum(is.na(plot2)), plot3 = sum(is.na(plot3))) # NA plot IDs indicate some cells were not matched, likely due to partitioning errors above
  #print(abaCellPlots %>% filter(is.na(heightError1matched)), n = 2, width = Inf)
  #abaCellPlots %>% summarize(naSpecies1 = sum(is.na(speciesError1matched)), naSpecies1unmatched = sum(is.na(speciesError1unmatched)), naHeight1 = sum(is.na(heightError1maeMatched)), naHeight1unmatched = sum(is.na(heightError1maeUnmatched)))
  
  return(abaCellPlots)
}

get_tree_lists = function(abaCellPlots, plotTreeCounts, trees2021lidar)
{
  ## missed tree list formation
  # Absent supporting data structures, dplyr requires an overlapped join to list all plot trees cell by cell. Since the 
  # overlap join generates orders of magnitude more rows than necessary, it's both slow and likely intractably memory bound 
  # (>45 s, >128 GB DDR). A workaround is to use pmax()
  # workaround, 
  #startTime = Sys.time()
  abaCellToPlotTreeListJoinFrame = left_join(left_join(left_join(abaCellPlots %>% select(abaGridX, abaGridY, stand1, stand2, plot1, plot2, plot3, n) %>% rename(nCell = n), # 
                                                                 plotTreeCounts %>% select(plot, n) %>% rename(nPlot1 = n),
                                                                 by = join_by(x$plot1 == y$plot), relationship = "many-to-one"), # multiple cells can be matched to any given plot
                                                       plotTreeCounts %>% select(plot, n) %>% rename(nPlot2 = n),
                                                       by = join_by(x$plot2 == y$plot), relationship = "many-to-one"),
                                             plotTreeCounts %>% select(plot, n) %>% rename(nPlot3 = n),
                                             by = join_by(x$plot3 == y$plot), relationship = "many-to-one") %>%
    mutate(n = pmax(nCell, nPlot1, nPlot2, nPlot3)) %>%
    uncount(n) %>% # ~9 s
    group_by(abaGridX, abaGridY) %>% 
    mutate(treeID = row_number()) %>% # uncount() duplicates tag numbers; add a tree ID column which provides a unique ID for each stand-plot-tag tree as expanded to an ABA cell
    ungroup()
  #Sys.time() - startTime
  
  #startTime = Sys.time()
  abaCellTreeLists = left_join(left_join(left_join(left_join(abaCellToPlotTreeListJoinFrame, # ~22 s, ~7 s per join
                                                             trees2021lidar %>% select(abaGridX, abaGridY, standID2016, treeID, species, height) %>% rename(stand = standID2016, speciesL = species, heightL = height), # TODO: include isPlantation, x, y?
                                                             by = join_by(abaGridX, abaGridY, treeID)),
                                                   plotTreeLists %>% select(plot, treeID, species, height) %>% rename(species1 = species, height1 = height),
                                                   by = join_by(x$plot1 == y$plot, treeID)),
                                         plotTreeLists %>% select(plot, treeID, species, height) %>% rename(species2 = species, height2 = height),
                                         by = join_by(x$plot2 == y$plot, treeID)),
                               plotTreeLists %>% select(plot, treeID, species, height) %>% rename(species3 = species, height3 = height),
                               by = join_by(x$plot3 == y$plot, treeID)) %>%
    mutate(stand = if_else(is.na(stand), stand1, stand), 
           isConiferL = if_else(speciesL %in% c("PSME", "TSHE", "THPL", "CHLA", "PISI", "TABR", "ABGR", "conifer", "PICO"), TRUE, FALSE),
           isConifer1 = if_else(species1 %in% c("PSME", "TSHE", "THPL", "CHLA", "PISI", "TABR", "ABGR", "conifer", "PICO"), TRUE, FALSE),
           isConifer2 = if_else(species2 %in% c("PSME", "TSHE", "THPL", "CHLA", "PISI", "TABR", "ABGR", "conifer", "PICO"), TRUE, FALSE),
           isConifer3 = if_else(species3 %in% c("PSME", "TSHE", "THPL", "CHLA", "PISI", "TABR", "ABGR", "conifer", "PICO"), TRUE, FALSE),
           species = if_else(is.na(speciesL) == FALSE, if_else(isConiferL != isConifer1, speciesL, species1), if_else(is.na(species1) == FALSE, species1, if_else(is.na(species2) == FALSE, species2, species3))),
           isConifer = if_else(species %in% c("PSME", "TSHE", "THPL", "CHLA", "PISI", "TABR", "ABGR", "conifer", "PICO"), TRUE, FALSE),
           height = if_else(is.na(heightL) == FALSE, heightL, height1)) %>% 
    filter(is.na(height) == FALSE) # remove rows for candidate trees which weren't accepted
  #Sys.time() - startTime
  # abaCellTreeLists %>% select(-nCell, -nPlot1, -nPlot2, -nPlot3)
  
  return(abaCellTreeLists)
}


## load data
# While Organon works in English units, treelists written by SEEM are in SI units so heights are in m
# Stands are also in SI units.
abaGrid = tibble(originX = 106020, originY = 190480, size = 20) %>% # m EPSG:6556, ABA grid origin and cell size from GIS/Trees/Elliott ABA grid 20 m.gpkg
  mutate(sizeHa = size^2 / 10000)

# stands defined on the Elliott State Research Forest in 2022 plus 
#  1) adjacent stands defined on the Elliott State Forest in 2016
#  2) adjacent stands on other ownerships
stands2022 = left_join(read_xlsx("GIS/Trees/2015-16 cruise.xlsx") %>% rename(stand = standID2016),
                       read_xlsx("GIS/Planning/Elliott Stand Data Feb2022.xlsx") %>% select(StandID, ODSL_VEG_L) %>% rename(vegLabel = ODSL_VEG_L),
                       by = join_by(x$stand == y$StandID)) %>%
  mutate(vegStrata = if_else(vegLabel %in% c("1D1L", "1D2H", "1D2L", "1D3H", "1D4H", "1D5H", "DX1L", "DX2H", "DX2L", "DX34L", "DX3H", "DX4H", "DX5H"), vegLabel, "other")) # cruise strata for cross validation of missing tree imputation from cruise.R
stands2016esf = as_tibble(vect("GIS/Planning/Elliott State Forest + Hakki stands 2016.gpkg")) %>% filter(standID2016 %in% stands2022$stand == FALSE) %>% rename(stand = standID2016)
stands2022 = bind_rows(stands2022, stands2016esf) # picks up areas of bordering stands

trees2021organon = left_join(read_feather(file.path(getwd(), "trees/Organon/Elliott tree lists 2016-2116.feather"), mmap = FALSE) %>% select(-species), # TODO: should 2016 snags be joined? they don't flow through Organon but standing in 2016 won't all have fallen by 2021
                             read_xlsx("trees/Elliott final cruise records 2015-16.xlsx", sheet = "CRUISERECS") %>% rename(stand = StandID, plot = PlotID, tag = TreeID, species = Species) %>% select(stand, plot, tag, species),
                             by = c("stand", "plot", "tag")) %>% # recover undubbed species from original plot measurements
  #mutate(species = case_match(species, 202 ~ "PSME", 351 ~ "ALRU", 263 ~ "TSHE", 312 ~ "ACMA", 242 ~ "THPL", 17 ~ "ABGR", 361 ~ "ARME", 431 ~ "CHCH", 631 ~ "NODE", 920 ~ "Salix", 81 ~ "CADE", 492 ~ "CONU", 815 ~ "QUGA", )) %>%
  filter(year %in% c(2021, 2026)) %>%
  arrange(stand, plot, tag, year) %>%
  group_by(stand, plot, tag) %>%
  mutate(heightGrowth5 = height[2] - height[1], dbhGrowth5 = dbh[2] - dbh[1]) %>% # get 2021-26 five year growth increments predicted by Organon, assume zero spring-summer mortality as moisture stress is low on the Elliott
  ungroup() %>%
  filter(year == 2021) %>%
  mutate(height = height + 0.9 / 5 * heightGrowth5, dbh = dbh + 0.5 / 5 * dbhGrowth5) # approximate sixth growth season: cruise data is from October 2015-February 2016 and the 2021 OLC Coos County flight covered the Elliott on August 30 and 31 2021, so the greater part of six height growth intervals exists between the cruise and LiDAR data (Harrington et al. 2016 https://www.fs.usda.gov/research/treesearch/53128, Gould et al. 2012 https://doi.org/10.1093/treephys/tps106).

stands2021organon = trees2021organon %>% 
  mutate(heightClass = round(height), isConifer = species %in% c("CX", "DF", "GF", "LP", "PC", "PY", "RC", "SS", "WH")) %>% 
  group_by(stand, heightClass) %>%
  summarize(tphLive = sum(liveExpansionFactor), tphConifer = sum(isConifer * liveExpansionFactor), tphHardwood = sum((isConifer == FALSE) * liveExpansionFactor),
            tphDead = sum(deadExpansionFactor), .groups = "drop")

plotTreeLists = left_join(trees2021organon,
                          stands2022 %>% select(stand, vegStrata, isPlantation, measurePlotsInStand), 
                          by = c("stand")) %>%
  mutate(species = forcats::fct_recode(factor(species), "POBA" = "BC", "ACMA" = "BM", "RHPU" = "CA", "Prunus" = "CH", "conifer" = "CX", "PSME" = "DF", "CHCH" = "GC", "ABGR" = "GF", "hardwood" = "HX", "PICO" = "LP", "FRLA" = "OA", "UMCA" = "OM", "CHLA" = "PC", "CONU" = "PD", "ARME" = "PM", "TABR" = "PY", "ALRU" = "RA", "THPL" = "RC", "PISI" = "SS", "NODE" = "TO", "TSHE" = "WH", "Salix" = "WI", "QUGA" = "WO", "other" = "XX"),
         treesPerCell = abaGrid$sizeHa * measurePlotsInStand * liveExpansionFactor, # expansion factors from Organon are stand-level; multiplying by the number of measure plots converts the expansion back to plot-level
         treesPerCell = if_else(treesPerCell < 1, 1, round(treesPerCell))) %>% # quantize number of trees per cell for uncount(); since these are 2016 plot measurements grown to 2021 by default at least one tree must be present on the plot to measure and thus at least one tree is present on the grid cell---this does not hold where a tree or snag falls or windsnaps but that data's not available, so assume 100% upright
  uncount(treesPerCell) %>%
  group_by(stand, plot) %>% 
  arrange(desc(height)) %>%
  mutate(n = n(), treeID = row_number()) %>% # uncount() duplicates tag numbers; add a tree ID column which provides a unique ID for each stand-plot-tag tree as expanded to an ABA cell
  ungroup()

plotTreeCounts = plotTreeLists %>% group_by(stand, plot) %>% summarize(n = n[1], .groups = "drop")

plotHeights = plotTreeLists %>% # ~28 seconds
  group_by(plot) %>% # regroup for summarize(), plotTreeLists is already sorted in descending height order
  summarize(stand = stand[1], n = n[1], vegStrata = vegStrata[1], isPlantation = isPlantation[1],
            species1 = species[1], height1 = height[1], # for now, assume all mortality falls immediately
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
plotMetrics2021 = list() # plot metrics created by metricsJob.R; z values all need to be converted from feet to m
for (chunkIndex in 1:19)
{
  plotMetrics2021[[chunkIndex]] = vect(file.path("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics", paste0("Elliott stdmetrics plot chunk ", chunkIndex, ".gpkg")))
}
plotMetrics2021 = as_tibble(vect(plotMetrics2021)) %>% # EPSG:6557 but geometry dropped, so no need to convert coordinates
  filter(PltInteger %in% plots2016$PltInteger, PlotType != "Count") %>% # count plots aren't directly useful to tree list establishment, also exclude dropped plots with known coordinates
  select(-Space, -STAND, -TotalPlots, -Count, -Measure, -Cruiser, -Plot_ID, -PlotType, -stems, -liveTrees, -snags, -primarySpecies, -n, -area) %>% # remove non-LiDAR observable or non-relevant variables
  select(-x, -y) %>% # prevent classification from mapping plots
  select(-zentropy) %>% # always NaN; lidR bug
  rename(plot = PltInteger, zMax = zmax, zMean = zmean, zStdDev = zsd, zSkew = zskew, zKurtosis = zkurt,
         pZaboveZmean = pzabovezmean, pZaboveThreshold = pzabove2, 
         zQ05 = zq5, zQ10 = zq10, zQ15 = zq15, zQ20 = zq20, zQ25 = zq25, zQ30 = zq30, zQ35 = zq35, zQ40 = zq40, zQ45 = zq45, zQ50 = zq50, zQ55 = zq55, zQ60 = zq60, zQ65 = zq65, zQ70 = zq70, zQ75 = zq75, zQ80 = zq80, zQ85 = zq85, zQ90 = zq90, zQ95 = zq95,
         zPcumulative10 = zpcum1, zPcumulative20 = zpcum2, zPcumulative30 = zpcum3, zPcumulative40 = zpcum4, zPcumulative50 = zpcum5, zPcumulative60 = zpcum6, zPcumulative70 = zpcum7, zPcumulative80 = zpcum8, zPcumulative90 = zpcum9,
         intensityTotal = itot, intensityMax = imax, intensityMean = imean, intensityStdDev = isd, intensitySkew = iskew, intensityKurtosis = ikurt, intensityPground = ipground,
         pCumulativeZQ10 = ipcumzq10, pCumulativeZQ30 = ipcumzq30, pCumulativeZQ50 = ipcumzq50, pCumulativeZQ70 = ipcumzq70, pCumulativeZQ90 = ipcumzq90, 
         pFirstReturn = p1th, pSecondReturn = p2th, pThirdReturn = p3th, pFourthReturn = p4th, pFifthReturn = p5th, pGround = pground) %>%
  mutate(across(starts_with(c("zM", "zQ")), ~0.3048 * .x), # convert mean-median-max and quantiles from feet to meters
         zStdDev = 0.3048 * zStdDev, zSkew = 0.3048^3 * zSkew, zKurtosis = 0.3048^4 * zKurtosis, # remaining English to SI conversions
         pZaboveZmean = 0.01 * pZaboveZmean, pZaboveThreshold = 0.01 * pZaboveThreshold, pGround = 0.01 * pGround, # fix lidR returning probabilities as 0–100 percentages rather than conforming to the 0–1 convention
         across(starts_with(c("pCumulative", "zPcumulative")), ~0.01 * .x), 
         pFirstReturn = 0.01 * pFirstReturn, pSecondReturn = 0.01 * pSecondReturn, pThirdReturn = 0.01 * pThirdReturn, pFourthReturn = 0.01 * pFourthReturn, pFifthReturn = 0.01 * pFifthReturn)

plotHeights = left_join(plotHeights, plotMetrics2021, by = c("plot")) # join LiDAR metrics to Organon grown cruise data

load("trees/height-diameter/data/trees DSM ring.Rdata") # LiDAR identified treetops from Get-Treetops and trees.R, a few seconds
trees2021lidar = elliottTreesMod %>% # ~7.6 s, EPSG:6556; heights and elevations also in m
  mutate(abaGridX = floor(1/abaGrid$size * (x - abaGrid$originX)), # ABA grid origin and cell size from GIS/Trees/Elliott ABA grid 20 m.gpkg
         abaGridY = floor(1/20 * (y - abaGrid$originY))) %>%
  rename(segmentationID = treeID) %>% # because treeID is changed to each detected treetop's height rank below
  group_by(abaGridX, abaGridY) %>%
  arrange(desc(height)) %>%
  mutate(nCell = n(), treeID = row_number()) %>%
  ungroup()
rm(elliottTreesMod) # since no way to load a single variable in a .Rdata into a specific variable name

trees2021lidarByHeightClass = trees2021lidar %>% # ~3 s
  rename(stand = standID2016) %>%
  mutate(heightClass = round(height), isConifer = species %in% c("PSME")) %>% 
  group_by(stand, isConifer, heightClass) %>%
  summarize(standArea = standArea[1], trees = n(), .groups = "drop") %>%
  mutate(tph = trees / standArea)

if (abaOptions$recalcAbaCellOccupancy)
{
  startTime = Sys.time()
  abaCellTreesByStand = trees2021lidar %>% group_by(abaGridX, abaGridY, standID2016) %>% # 1.8 minutes
    summarize(isPlantation = isPlantation[1], nStand = n(), .groups = "drop_last") %>% 
    slice_max(nStand, n = 2, with_ties = FALSE) %>%
    mutate(rowNumber = row_number()) %>%
    ungroup() %>%
    rename(stand = standID2016) %>%
    pivot_wider(id_cols = c("abaGridX", "abaGridY"), names_from = "rowNumber", names_sep = "", values_from = c("stand", "isPlantation", "nStand"))
  Sys.time() - startTime
  #abaCellTreesByStand %>% filter(is.na(stand2) == FALSE)

  # 15,367 cells (1.5%) have zero treetops and are thus silently dropped out of the grid occupancy formed below
  # These are typically canopy gaps with multiple treetops just out of cell.
  startTime = Sys.time()
  abaCellTrees = trees2021lidar %>% # ~49 minutes; 11.4 M treetops -> million row tibble
    summarize(stands = length(unique(standID2016)), n = n(), species1 = species[1], height1 = height[1], 
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
  
  abaCells = left_join(abaCellTreesByStand, abaCellTrees, by = c("abaGridX", "abaGridY")) %>%
    relocate(abaGridX, abaGridY, stands, n)
  save(file = "trees/height-diameter/data/trees by cell.Rdata", abaCells)
} else {
  load("trees/height-diameter/data/trees by cell.Rdata") # abaCells
}

abaMetrics = as.data.frame(project(rast("D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/metrics/grid metrics 20 m.tif"), crs("epsg:6556"), threads = TRUE), xy = TRUE, na.rm = NA) # ~5 s
names(abaMetrics)[47] = "intensityPground" # temporary workaround for C# typo
abaMetrics %<>% mutate(abaGridX = floor(1/abaGrid$size * (x - abaGrid$originX)), 
                       abaGridY = floor(1/20 * (y - abaGrid$originY)),
                       across(starts_with(c("zM", "zQ")), ~0.3048 * .x), # convert mean-median-max and quantiles from feet to meters
                       zStdDev = 0.3048 * zStdDev, zSkew = 0.3048^3 * zSkew, zKurtosis = 0.3048^4 * zKurtosis) # remaining English to SI conversions)

abaCells = left_join(abaCells, abaMetrics %>% rename(nPoints = n), by = c("abaGridX", "abaGridY"))


abaCellsScaled = abaCells %>% mutate(treesMatched = if_else(n < 8, n, 8),
                                     isConifer1 = abaOptions$coniferDeciduousWeight * isConifer1,
                                     isConifer2 = abaOptions$coniferDeciduousWeight * isConifer2,
                                     isConifer3 = abaOptions$coniferDeciduousWeight * isConifer3,
                                     isConifer4 = abaOptions$coniferDeciduousWeight * isConifer4,
                                     isConifer5 = abaOptions$coniferDeciduousWeight * isConifer5,
                                     isConifer6 = abaOptions$coniferDeciduousWeight * isConifer6,
                                     isConifer7 = abaOptions$coniferDeciduousWeight * isConifer7,
                                     isConifer8 = abaOptions$coniferDeciduousWeight * isConifer8,
                                     isConifer9 = abaOptions$coniferDeciduousWeight * isConifer9,
                                     isConifer10 = abaOptions$coniferDeciduousWeight * isConifer10,
                                     isConifer11 = abaOptions$coniferDeciduousWeight * isConifer11,
                                     isConifer12 = abaOptions$coniferDeciduousWeight * isConifer12,
                                     isConifer13 = abaOptions$coniferDeciduousWeight * isConifer13,
                                     isConifer14 = abaOptions$coniferDeciduousWeight * isConifer14,
                                     isConifer15 = abaOptions$coniferDeciduousWeight * isConifer15,
                                     intensityMean = 0.001 * treesMatched * intensityMean,
                                     pGround = 75 * treesMatched * 100 * pGround, # weight for nominal 75 m tree height @ pGround = 1% of all returns
                                     pZaboveThreshold = 75 * treesMatched * pZaboveThreshold,
                                     zMean = treesMatched * zMean,
                                     zQ10 = 4.0 * treesMatched * zQ10, # plot favoring weights from relative height analysis in exploratory section below
                                     zQ20 = 2.7 * treesMatched * zQ20,
                                     zQ30 = 2.0 * treesMatched * zQ30,
                                     zQ40 = 1.8 * treesMatched * zQ40,
                                     zQ50 = 1.6 * treesMatched * zQ50, # 1.6 zQ60, 1.4 zQ70, 1.2 zQ90
                                     zQ80 = 1.3 * treesMatched * zQ80) %>%
  filter(is.na(zMean) == FALSE) # 386 cells without LiDAR metrics due to tight south side flight boundary in 2021
plotHeightsScaled = plotHeights %>% mutate(treesMatched = if_else(n < 8, n, 8),
                                           isConifer1 = abaOptions$coniferDeciduousWeight * isConifer1,
                                           isConifer2 = abaOptions$coniferDeciduousWeight * isConifer2,
                                           isConifer3 = abaOptions$coniferDeciduousWeight * isConifer3,
                                           isConifer4 = abaOptions$coniferDeciduousWeight * isConifer4,
                                           isConifer5 = abaOptions$coniferDeciduousWeight * isConifer5,
                                           isConifer6 = abaOptions$coniferDeciduousWeight * isConifer6,
                                           isConifer7 = abaOptions$coniferDeciduousWeight * isConifer7,
                                           isConifer8 = abaOptions$coniferDeciduousWeight * isConifer8,
                                           isConifer9 = abaOptions$coniferDeciduousWeight * isConifer9,
                                           isConifer10 = abaOptions$coniferDeciduousWeight * isConifer10,
                                           isConifer11 = abaOptions$coniferDeciduousWeight * isConifer11,
                                           isConifer12 = abaOptions$coniferDeciduousWeight * isConifer12,
                                           isConifer13 = abaOptions$coniferDeciduousWeight * isConifer13,
                                           isConifer14 = abaOptions$coniferDeciduousWeight * isConifer14,
                                           isConifer15 = abaOptions$coniferDeciduousWeight * isConifer15,
                                           intensityMean = 0.001 * treesMatched * intensityMean,
                                           pGround = 75 * treesMatched * 100 * pGround, # lidR returns probabilities as 0-100%
                                           pZaboveThreshold = 75 * treesMatched * pZaboveThreshold, # lidR probabilities
                                           zMean = treesMatched * zMean,
                                           zQ10 = 4.0 * treesMatched * zQ10,
                                           zQ20 = 2.7 * treesMatched * zQ20,
                                           zQ30 = 2.0 * treesMatched * zQ30,
                                           zQ40 = 1.8 * treesMatched * zQ40,
                                           zQ50 = 1.6 * treesMatched * zQ50, 
                                           zQ80 = 1.3 * treesMatched * zQ80) %>%
  filter(is.na(zMean) == FALSE) # 30 plots without LiDAR metrics due to target coordinates not being available


## exploratory analysis
if (abaOptions$includeInvestigatory)
{
  # trees of questionable height: simple test for incompletely classified noise
  highTreesByTile = trees2021lidar %>% group_by(tile) %>% reframe(treetops = n(), quantiles = c(0, 10, 25, 50, 75, 90, 100), heightQ = quantile(height, probs = 0.01 * quantiles), heightMean = mean(height), heightOver75 = sum(height > 75)) # 75 m threshold
  highTreesByTileWide = highTreesByTile %>% pivot_wider(id_cols = c("tile", "trees", "heightMean", "heightOver75"), names_prefix = "heightQ", names_from = "quantiles", values_from = "heightQ")
  highTreesByTileWide %>% slice_max(heightOver75, n = 10)
  highTreesByTile %>% filter(tile == "s04080w06720")
  
  ggplot() +
    geom_histogram(aes(x = heightOver75), highTreesByTileWide, binwidth = 5) +
    labs(x = "trees more than 75 m tall", y = "tiles")

  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = 75, yend = 75), color = "grey70", linetype = "longdash", linewidth = 0.3) +
    geom_bin2d(aes(x = heightQ90, y = heightQ100), highTreesByTileWide, binwidth = 1) +
    coord_equal() +
    labs(x = bquote("tile 90"^th*" percentile tree height, m"), y = "tile maximum tree height, m", fill = "tiles") +
    scale_fill_viridis_c()
  
  #writexl::write_xlsx(highTreesByTileWide, "trees/tile treetops.xlsx")

  # distribution of trees per plot after resizing to ABA cells
  # 250 feet = 76.2 m
  ggplot() +
    geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..))), plotHeights, binwidth = 1) +
    labs(x = "plot trees expanded to ABA cell", y = "probability of plot") +
  ggplot() +
    geom_histogram(aes(x = n, y = after_stat(..count.. / sum(..count..))), abaCells, binwidth = 1) +
    labs(x = "treetops detected in ABA cell", y = "probability of cell") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout() &
    coord_trans(x = scales::transform_pseudo_log(sigma = 2), xlim = c(0, 500), ylim = c(0, 0.08)) &
    scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900)) &
    scale_y_continuous(labels = scales::label_percent())
  
  abaCellTreetopEcdf = bind_rows(tibble(n = 0, cells = 15367), abaCells %>% group_by(n) %>% summarize(cells = n())) %>% arrange(n) %>% mutate(ecdf = cumsum(cells) / sum(cells))
  plotTreesEcdf = bind_rows(tibble(n = 0, plots = 0), plotHeights %>% group_by(n) %>% summarize(plots = n())) %>% arrange(n) %>% mutate(ecdf = cumsum(plots) / sum(plots))
  missedTrees = full_join(full_join(tibble(n = seq(0, max(abaCellTreetopEcdf$n))), plotTreesEcdf %>% rename(plotEcdf = ecdf), by = join_by(n)), 
                          abaCellTreetopEcdf %>% rename(cellEcdf = ecdf), by = join_by(n)) %>%
    fill(ends_with("Ecdf"), .direction = "down") %>%
    uncount(weights = 2) %>% # no step option for geom_ribbon() so duplicate each row and then lead n to create ateps
    mutate(n = lead(n, default = max(n) + 1))
  ggplot() +
    geom_ribbon(aes(x = n, ymin = plotEcdf, ymax = if_else(n > 0, cellEcdf, 0), fill = "undetected trees"), missedTrees, alpha = 0.2) +
    geom_step(aes(x = n, y = ecdf, color = "ABA cells"), abaCellTreetopEcdf) +
    geom_step(aes(x = n, y = ecdf, color = "cruise plots"), plotTreesEcdf) +
    coord_trans(x = scales::transform_pseudo_log(sigma = 2), xlim = c(0, 500)) +
    labs(x = "trees per cell", y = "cumulative probability", color = NULL, fill = NULL) +
    scale_color_manual(breaks = c("ABA cells", "cruise plots"), values = c("forestgreen", "blue")) +
    scale_fill_manual(breaks = c("undetected trees"), values = c("green2")) +
    scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900)) +
    scale_y_continuous(labels = scales::label_percent())
  
  # plot and grid cell height quantiles
  # Plots use normalized data, so ground height does not need to be subtracted.
  # Normalized grid cells are also loaded (at time of coding) for consistency with plots.
  abaCellQuantiles = abaMetrics %>% rename(zQ100 = zMax) %>% select(abaGridX, abaGridY, starts_with("zQ")) %>% # ~8.9 seconds
    pivot_longer(cols = starts_with("zQ"), names_pattern = "zQ(.*)", names_to = "quantile", values_to = "zQ") %>% 
    mutate(quantile = as.numeric(quantile)) %>%
    group_by(abaGridX, abaGridY) %>%
    mutate(relativeHeight = zQ / max(zQ)) %>%
    ungroup()
  abaCellCumulativeHeightQuantiles = abaMetrics %>% select(abaGridX, abaGridY, starts_with("zPcumulative")) %>% # ~2.0 s
    pivot_longer(cols = starts_with("zPcumulative"), names_pattern = "zPcumulative(.*)", names_to = "quantile", values_to = "zQ") %>% 
    mutate(quantile = as.numeric(quantile))
  abaCellReturns = abaMetrics %>% rename(pReturn1 = pFirstReturn, pReturn2 = pSecondReturn, pReturn3 = pThirdReturn, pReturn4 = pFourthReturn, pReturn5 = pFifthReturn) %>%
    select(abaGridX, abaGridY, starts_with("pReturn")) %>%
    pivot_longer(cols = starts_with("pReturn"), names_pattern = "pReturn(.*)", names_to = "return", values_to = "probability") %>% 
    mutate(return = as.factor(return))

  plotHeightQuantiles = plotHeights %>% rename(zQ100 = zMax) %>% select(plot, starts_with("zQ")) %>% 
    pivot_longer(cols = starts_with("zQ"), names_pattern = "zQ(.*)", names_to = "quantile", values_to = "zQ") %>% 
    mutate(quantile = as.numeric(quantile)) %>%
    group_by(plot) %>%
    mutate(relativeHeight = zQ / max(zQ)) %>%
    ungroup()
  plotCumulativeHeightQuantiles = plotHeights %>% select(plot, starts_with("zPcumulative")) %>% 
    pivot_longer(cols = starts_with("zPcumulative"), names_pattern = "zPcumulative(.*)", names_to = "quantile", values_to = "zQ") %>% 
    mutate(quantile = as.numeric(quantile))
  plotReturns = plotHeights %>% rename(pReturn1 = pFirstReturn, pReturn2 = pSecondReturn, pReturn3 = pThirdReturn, pReturn4 = pFourthReturn, pReturn5 = pFifthReturn) %>%
    select(plot, starts_with("pReturn")) %>%
    pivot_longer(cols = starts_with("pReturn"), names_pattern = "pReturn(.*)", names_to = "return", values_to = "probability") %>% 
    mutate(return = as.factor(return))
  
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = zQ, group = quantile, fill = quantile), plotHeightQuantiles %>% filter(quantile %in% seq(0, 100, by = 10)), binwidth = 1, na.rm = TRUE) + # drop 5% quantiles for readability
    coord_cartesian(xlim = c(0, 0.04), ylim = c(0, 80)) +
    labs(x = NULL, y = "height on plot, m", fill = "quantile") +
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = relativeHeight, group = quantile, fill = quantile), plotHeightQuantiles %>% filter(quantile %in% seq(0, 95, by = 10)), binwidth = 0.01, na.rm = TRUE) + # drop 5% quantiles for readability
    coord_cartesian(xlim = c(0, 0.025), ylim = c(0, 1)) +
    labs(x = NULL, y = "relative height on plot", fill = "quantile") +
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = zQ, group = quantile, fill = quantile), plotCumulativeHeightQuantiles %>% filter(quantile %in% seq(0, 100, by = 10)), binwidth = 0.01, na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 0.04), ylim = c(0, 1)) +
    labs(x = NULL, y = "cumulative relative\npoint height on plot, m", fill = "quantile") +
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = zQ, group = quantile, fill = quantile), abaCellQuantiles %>% filter(quantile %in% seq(0, 100, by = 10)), binwidth = 1, na.rm = TRUE) + # drop 5% quantiles for readability
    coord_cartesian(xlim = c(0, 0.04), ylim = c(0, 80)) +
    labs(x = "probability", y = "height on cell, m", fill = "quantile") +
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = relativeHeight, group = quantile, fill = quantile), abaCellQuantiles %>% filter(quantile %in% seq(0, 95, by = 10)), binwidth = 0.01, na.rm = TRUE) + # drop 5% quantiles for readability
    coord_cartesian(xlim = c(0, 0.025), ylim = c(0, 1)) +
    labs(x = "probability", y = "relative height on cell", fill = "quantile") +
  ggplot() +
    geom_histogram(aes(x = after_stat(..count.. / sum(..count..)), y = zQ, group = quantile, fill = quantile), abaCellCumulativeHeightQuantiles %>% filter(quantile %in% seq(0, 100, by = 10)), binwidth = 0.01, na.rm = TRUE) +
    coord_cartesian(xlim = c(0, 0.04), ylim = c(0, 1)) +
    labs(x = "probability", y = "cumulative relative\npoint height on cell, m", fill = "quantile") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_fill_viridis_c(limits = c(0, 100)) &
    scale_x_continuous(labels = scales::percent_format())
  #ggsave("trees/area based/figures/Figure S1 cell-plot quantiles.png", height = 10.5, width = 22, units = "cm", dpi = 250)
  
  ggplot() +
    geom_histogram(aes(x = probability, y = after_stat(..count.. / sum(..count..)), group = return, fill = return), plotReturns, binwidth = 0.1 / 100) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.015)) +
    labs(x = "percentage of total returns", y = "probability on plot", fill = "return") +
  ggplot() +
    geom_histogram(aes(x = pGround, y = after_stat(..count.. / sum(..count..))), plotHeights, binwidth = 0.01 / 100) +
    coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.015)) +
    labs(x = "percentage of total returns", y = "ground probability on plot", fill = "return") +
  ggplot() +
    geom_histogram(aes(x = pZaboveThreshold, y = after_stat(..count.. / sum(..count..))), plotHeights, binwidth = 0.01) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.25)) +
    labs(x = bquote("probability z > z"[threshold]), y = "probability on plot", fill = "return") +
  ggplot() +
    geom_histogram(aes(x = probability, y = after_stat(..count.. / sum(..count..)), group = return, fill = return), abaCellReturns, binwidth = 0.1 / 100) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.015)) +
    labs(x = "percentage of total returns", y = "probability on cell", fill = "return") +
  ggplot() +
    geom_histogram(aes(x = pGround, y = after_stat(..count.. / sum(..count..))), abaMetrics, binwidth = 0.01 / 100) +
    coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.015)) +
    labs(x = "percentage of total returns", y = "ground probability on cell", fill = "return") +
  ggplot() +
    geom_histogram(aes(x = pZaboveThreshold, y = after_stat(..count.. / sum(..count..))), abaMetrics, binwidth = 0.01) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.25)) +
    labs(x = bquote("probability z > z"[threshold]), y = "probability on cell", fill = "return") +
  plot_annotation(theme = theme(plot.margin = margin())) +
  plot_layout(guides = "collect") &
    scale_fill_viridis_d(direction = -1) &
    scale_x_continuous(labels = scales::percent_format()) &
    scale_y_continuous(labels = scales::percent_format())
  #ggsave("trees/area based/figures/Figure S2 cell-plot returns and ground.png", height = 10.5, width = 22, units = "cm", dpi = 250)
  
  # guidance for kNN distance scaling of height quantiles
  bind_cols(abaCellQuantiles %>% group_by(quantile) %>% summarise(abaWeight = 1 / mean(relativeHeight, na.rm = TRUE)),
            plotHeightQuantiles %>% group_by(quantile) %>% summarise(abaWeight = 1 / mean(relativeHeight, na.rm = TRUE)))
  
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
}