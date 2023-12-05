# run height-diameter/setup.R to get trees2016 and stands2022

load("trees/height-diameter/data/PSME preferred models.Rdata")
load("trees/height-diameter/data/ALRU2 preferred models.Rdata")
load("trees/height-diameter/data/TSHE preferred models.Rdata")
load("trees/height-diameter/data/ACMA3 preferred models.Rdata")
load("trees/height-diameter/data/UMCA preferred models.Rdata")
load("trees/height-diameter/data/THPL preferred models.Rdata")
load("trees/height-diameter/data/other preferred models.Rdata")


## adjust stands2022 and trees2016 for use with Organon
# TODO: site index accuracy, site indices for species besides Douglas-fir
organonStands = stands2022 %>% 
  mutate(siteIndex = 0.3048 * if_else(is.na(Cruised_Si) == FALSE, Cruised_Si, pmin(ODSL_Site_, ODSL_Physi)), # convert measured SI from ft to m if available, otherwise use lowest modeled value
         siteIndex = replace_na(siteIndex, mean(siteIndex, na.rm = TRUE)), # absent any other approach, fill missing site indices using mean imputation (nearest neighbors or such could be done to clear this, in GIS or elsewhere)
         slopeInPercent = 100 * tan(pi/180 * SlopeMedian), # %
         forwardingRoad = 30, # m
         forwardingUntethered = pmin(15, 2 * RoadDistMean), # m
         forwardingTethered = 2 * RoadDistMean - forwardingUntethered,
         plantingDensityPerHa = 2.47 * 360, # trees per hectare
         yardingFactor = if_else(RoadDistMedian > 0, 0.5 * RoadDistMean / RoadDistMedian, 0.5), # dimensionless
         slopeAbove100PercentFraction = SlopeAbove100PercentFraction) %>% # fraction
  rename(id = StandID, area = standArea, age = standAge2016, allocation = April2021_Allocation) %>%
  select(id, area, siteIndex, age, slopeInPercent, forwardingRoad, forwardingUntethered, forwardingTethered, yardingFactor, plantingDensityPerHa, slopeAbove100PercentFraction)


# Dubbing to proxy species supported by Organon SWO (supported codes are DF, WF, GF, PP, SP, IC, WH, RC, PY, PM, GC, TO, LO?, BM, WO, BO?, RA, WI):
#  BC (Populus trichocarpa) -> RA (Alnus rubra)
#  CA (Frangula purshiana) -> GC (Chrysolepis chrysophylla)
#  CH (Prunus spp, presumably mostly Prunus emarginata) -> WI (Salix)
#  CX (conifer) -> IC (Calocedrus decurrens)
#  HX (hardwood) -> PM (Arbutus menziesii)
#  LP (Pinus contorta) -> WH (Tsuga heterophylla), since Elliott is coastal
#  OA (Fraxinus latifolia) -> BM (Acer macrophyllum)
#  OM (Umbelularia californica) -> BM (Acer macrophyllum)
#  PC (Chamaecyparis lawsoniana) -> RC (Thuja plicata)
#  SS (Picea sitchensis) -> WH (Tsuga heterophylla)
#  XX (unknown, likely predominantly conifer) -> IC (Calocedrus decurrens)
# Assumption here is it's a smaller error to include a tree with the wrong species than to omit the tree.
#
# Some Douglas-fir and red alder are on plots which aren't georeferenced. Fall back to non-physio DBH prediction for these trees.
organonTrees = trees2016 %>% filter(isLive, is.na(DBH) == FALSE) %>%
  mutate(species = fct_collapse(factor(Species, levels = c("BC", "BM", "CA", "CH", "CX", "DF", "GC", "GF", "HX", "LP", "OA", "OM", "PC", "PD", "PM", "PY", "RA", "RC", "SS", "TO", "WH", "WI", "WO", "XX")),
                                RA = c("BC"), GC = c("CA"), WI = c("CH"), IC = c("CX", "XX"), PM = c("HX"), WH = c("LP", "SS"), BM = c("OA", "OM"), RC = c("PC"))) %>%
  group_by(speciesGroup) %>%
  mutate(expansionFactor = meanTreesPerBafPlot / meanTreesPerBafMeasurePlot * measureTreeTphContribution / measurePlotsInStand, # trees per hectare
         imputedHeight = case_when(speciesGroup == "DF" ~ if_else(is.na(elevation) == FALSE, predict(psmeHeightFromDiameterPreferred$sharmaPartonBalPhysioRelDbh, .[cur_group_rows(), ]), predict(psmeHeightFromDiameterPreferred$sharmaPartonBalRelDbh, .[cur_group_rows(), ])),
                                   speciesGroup == "RA" ~ if_else(is.na(elevation) == FALSE, predict(alruHeightFromDiameterPreferred$sharmaPartonBalPhysio, .[cur_group_rows(), ]), predict(alruHeightFromDiameterPreferred$sharmaPartonBal, .[cur_group_rows(), ])),
                                   speciesGroup == "WH" ~ predict(tsheHeightFromDiameterPreferred$gamBalRelDbh, .[cur_group_rows(), ]),
                                   speciesGroup == "BM" ~ predict(acmaHeightFromDiameterPreferred$sharmaParton, .[cur_group_rows(), ]),
                                   speciesGroup == "OM" ~ predict(umcaHeightFromDiameterPreferred$sharmaPartonPhysio, .[cur_group_rows(), ]),
                                   speciesGroup == "RC" ~ predict(thplHeightFromDiameterPreferred$sharmaPartonPhysio, .[cur_group_rows(), ]),
                                   speciesGroup == "other" ~ predict(otherHeightFromDiameterPreferred$gamBal, .[cur_group_rows(), ]))) %>% # a few seconds, likely mainly due to GAM prediction
  rename(stand = StandID, plot = PlotID, tag = TreeID, dbh = DBH, height = imputedHeight, codes = CompCode, heightToBrokenTop = Ht2) %>%
  select(speciesGroup, stand, plot, tag, species, dbh, height, expansionFactor, codes, heightToBrokenTop) %>% # dbh in cm, heights in m, expansion factor in TPH
  ungroup() %>%
  arrange(stand, plot, tag)

# sanity checks
organonStands %>% reframe(quantiles = c(0, 1), age = quantile(age, probs = quantiles, na.rm = TRUE), 
                                               siteIndex = quantile(siteIndex, probs = quantiles), 
                                               slopeInPercent = quantile(slopeInPercent, probs = quantiles), 
                                               forwardingRoad = quantile(forwardingRoad, probs = quantiles),
                                               forwardingUntethered = quantile(forwardingUntethered, probs = quantiles),
                                               forwardingTethered = quantile(forwardingTethered, probs = quantiles),
                                               plantingDensityPerHa = quantile(plantingDensityPerHa, probs = quantiles),
                                               yardingFactor = quantile(yardingFactor, probs = quantiles),
                                               slopeAbove100PercentFraction = quantile(slopeAbove100PercentFraction, probs = quantiles))
organonTrees %>% group_by(speciesGroup) %>% 
  reframe(quantiles = c(0, 0.025, 0.5, 0.975, 1), dbh = quantile(dbh, probs = quantiles), height = quantile(height, probs = quantiles), EF = quantile(expansionFactor, probs = quantiles)) %>%
  pivot_wider(names_from = quantiles, values_from = c(dbh, height, EF))

intensiveStands = left_join(stands2022 %>% filter(April2021_Allocation == "Intensive"),
                            trees2016 %>% filter(isLive) %>% group_by(StandID) %>% summarize(measurePlots = measurePlotsInStand[1], uniqueNonReserveMeasureTrees = sum((is.na(DBH) == FALSE) & (CompCode != "RT"))),
                            by = "StandID") %>%
  filter(is.na(measurePlots) == FALSE, uniqueNonReserveMeasureTrees >= 20) # exclude uncruised stands and stands with too few measure trees to form meaningful prescription guidance

intensiveStands %>% filter(uniqueNonReserveMeasureTrees < 20)
intensiveStands$uniqueMeasureTrees[which(intensiveStands$StandID == 2445)]

intensiveStandIDs = unique(intensiveStands$StandID)
intensiveTrees = organonTrees %>% filter(stand %in% intensiveStandIDs) %>% select(-speciesGroup)
#write_xlsx(list(stands = organonStands, trees = organonTrees %>% select(-speciesGroup), intensiveTrees = intensiveTrees), "trees/Organon/Elliott Organon cruise records 2015-16.xlsx")


## exploratory plots
# distance from road
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1000, yend = 1000), color = "grey70", linetype = "longdash") +
  geom_linerange(aes(x = RoadDistMedian, ymin = RoadDistMin, ymax = RoadDistMax, color = standArea), stands2022, alpha = 0.1, linewidth = 0.3) +
  geom_point(aes(x = RoadDistMedian, y = RoadDistMean, color = standArea), stands2022, alpha = 0.1, size = 0.5) +
  labs(x = "median distance to road, m", y = "distance to road, m", color = "area, ha") +
  scale_color_viridis_c(limits = c(1, 100), trans = "log10") +
  theme(legend.spacing.y = unit(0.3, "line"))

ggplot() +
  geom_histogram(aes(x = RoadDistMin, weight = standArea), stands2022 %>% filter(standArea >= 1), binwidth = 25) +
  coord_cartesian(ylim = c(1, 20000)) +
  labs(x = "minimum distance to road, m", y = "area, ha") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 20000), labels = scales::comma_format())

ggplot() +
  geom_histogram(aes(x = 0.5 * RoadDistMean / RoadDistMedian, weight = standArea), stands2022 %>% filter(standArea >= 1), binwidth = 0.05) +
  coord_cartesian(ylim = c(1, 20000)) +
  labs(x = "crude estimate of yarding factor", y = "area, ha") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 20000), labels = scales::comma_format())

        
# site index: no apparent relationship between cruised and modeled site indices
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 40, yend = 40), color = "grey70", linetype = "longdash") +
  #geom_smooth(aes(x = 0.3048 * Cruised_Si, y = 0.3048 * ODSL_Site_), stands2022, alpha = 0.1, color = "green3", linewidth = 0.6, method = "gam", na.rm = TRUE) +
  geom_point(aes(x = 0.3048 * Cruised_Si, y = 0.3048 * ODSL_Site_, color = standAge2016), stands2022, na.rm = TRUE) +
  labs(x = "cruised site index, m", y = "DSL site index, m", color = "stand age, years") +
ggplot(stands2022) +
  geom_segment(aes(x = 0, y = 0, xend = 40, yend = 40), color = "grey70", linetype = "longdash") +
  #geom_smooth(aes(x = 0.3048 * Cruised_Si, y = 0.3048 * ODSL_Physi), stands2022, alpha = 0.1, color = "green3", linewidth = 0.6, method = "gam", na.rm = TRUE) +
  geom_point(aes(x = 0.3048 * Cruised_Si, y = 0.3048 * ODSL_Physi, color = standAge2016), stands2022, na.rm = TRUE) +
  labs(x = "cruised site index, m", y = "DSL physical site index, m", color = "stand age, years") +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(guides = "collect")

cruisedSImodel = lm(Cruised_Si ~ ODSL_Site_, stands2022 %>% mutate(Cruised_Si = 0.3048 * Cruised_Si, ODSL_Site_ = 0.3048 * ODSL_Site_))
summary(cruisedSImodel)

# measure trees per stand
ggplot() +
  geom_histogram(aes(x = uniqueMeasureTrees), intensiveStands, binwidth = 2) +
  labs(x = "measure trees", y = "intensive stands")


## management allocations
stands2022 %>% filter(StandID %in% unique(trees2016$StandID)) %>% group_by(April2021_Allocation) %>% 
  summarize(n = n(), area = sum(area_ha), age2016min = min(standAge2016), age2016max = max(standAge2016), 
            n20 = sum(standAge2016 <= 20), area20 = sum((standAge2016 <= 20) * area_ha),
            n40 = sum(standAge2016 <= 40), area40 = sum((standAge2016 <= 40) * area_ha))
stands2022 %>% filter(StandID %in% unique(trees2016$StandID), April2021_Allocation == "Intensive")
#stands2022 %>% filter(is.na(Cruised_Si) == FALSE) %>% summarise(ODSL_Physi = cor(Cruised_Si, ODSL_Physi))
#stands2022 %>% filter(is.na(Cruised_Si) == FALSE, is.na(ODSL_Site_) == FALSE) %>% summarise(ODSL_Site = cor(Cruised_Si, ODSL_Site_))

esrfAllocations = read_xlsx("GIS/Planning/ESRF_Allocations.xlsx")
esrfAllocations %>% group_by(April2021_Allocation) %>% summarize(n = n(), areaHa = sum(Acres) / 2.47105)
# approximate management costs (https://www.oregon.gov/dsl/Land/Elliott%20Forest%20Library/ElliottExpensesOverview.pdf)
tibble(managementPerIntensiveHa = 1.7E6 / 5829) # US$/ha-year
