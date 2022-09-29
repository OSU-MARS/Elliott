# Chaney NW, Minasny B, Herman JD, Nauman TP, et al. 2019. POLARIS Soil Properties: 30-m Probabilistic Maps of Soil 
#   Properties Over the Contiguous United States. Water Resources Research 55(4):2916-2938. 
#   https://doi.org/10.1029/2018WR022797
# See Table 1 for parameter definitions. For van Genuchten soil water retention,
#   α values in rasters = log₁₀(α)
#   θs = saturated water content, m³/m³
#   θr = residual water content, m³/m³
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(tidyr)

theme_set(theme_bw() + theme(axis.line = element_line(size = 0.5),
                             legend.background = element_rect(fill = alpha("white", 0.5)),
                             legend.margin = margin(),
                             panel.border = element_blank()))

## multilayer vs single layer water retention
# conversion factor for α: 1 cm H₂O (conventional) = 98.0665 Pa (Thompson and Taylor, Table B.8)
#   Thompson A, Taylor BN. 2008. Guide for the Use of the International System of Units (SI). NIST Special Publication 
#     811, National Institute of Standards and Technology, US Department of Commerce. https://physics.nist.gov/cuu/pdf/sp811.pdf
# Note apropos soil layering literature: Geometric mean is intractable due to large differences in powers, harmonic 
# mean has greater error and bias than arithmetic mean.
getWaterRetention = function(climateCells)
{
  # weighted arithmetic means of alpha and n appears to usually be a slight overestimation of values best
  # matching equipotential multilayer curves: alpha +0.8%, n +0.001%
  rockFraction = 0.05
  waterRetention = crossing(climateCells %>%
                              mutate(alpha = 1 / soilDepthInCMmean * (pmax(0, pmin(5, soilDepthInCMmean)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm0_5mean + # kPa
                                                                      pmax(0, pmin(10, soilDepthInCMmean - 5)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm5_15mean + 
                                                                      pmax(0, pmin(15, soilDepthInCMmean - 15)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm15_30mean + 
                                                                      pmax(0, pmin(30, soilDepthInCMmean - 30)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm30_60mean + 
                                                                      pmax(0, pmin(40, soilDepthInCMmean - 60)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm60_100mean + 
                                                                      pmax(0, pmin(100, soilDepthInCMmean - 100)) * 0.0980665 * 10^vanGenuchtenLog10AlphaInCm100_200mean),
                                     #alphaHarmonic = soilDepthInCMmean / (pmax(0, pmin(5, soilDepthInCMmean)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm0_5mean) +
                                     #                                     pmax(0, pmin(10, soilDepthInCMmean - 5)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm5_15mean) +
                                     #                                     pmax(0, pmin(15, soilDepthInCMmean - 15)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm15_30mean) +
                                     #                                     pmax(0, pmin(30, soilDepthInCMmean - 30)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm30_60mean) + 
                                     #                                     pmax(0, pmin(40, soilDepthInCMmean - 60)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm60_100mean) +
                                     #                                     pmax(0, pmin(100, soilDepthInCMmean - 100)) / (0.0980665 * 10^vanGenuchtenLog10AlphaInCm100_200mean)),
                                     n = 1 / soilDepthInCMmean * (pmax(0, pmin(5, soilDepthInCMmean)) * vanGenuchtenN0_5mean +
                                                                  pmax(0, pmin(10, soilDepthInCMmean - 5)) * vanGenuchtenN5_15mean + 
                                                                  pmax(0, pmin(15, soilDepthInCMmean - 15)) * vanGenuchtenN15_30mean + 
                                                                  pmax(0, pmin(30, soilDepthInCMmean - 30)) * vanGenuchtenN30_60mean + 
                                                                  pmax(0, pmin(40, soilDepthInCMmean - 60)) * vanGenuchtenN60_100mean + 
                                                                  pmax(0, pmin(100, soilDepthInCMmean - 100)) * vanGenuchtenN60_100mean),
                                     #nHarmonic = soilDepthInCMmean / (pmax(0, pmin(5, soilDepthInCMmean)) / vanGenuchtenN0_5mean +
                                     #                                 pmax(0, pmin(10, soilDepthInCMmean - 5)) / vanGenuchtenN5_15mean +
                                     #                                 pmax(0, pmin(15, soilDepthInCMmean - 15)) / vanGenuchtenN15_30mean +
                                     #                                 pmax(0, pmin(30, soilDepthInCMmean - 30)) / vanGenuchtenN30_60mean + 
                                     #                                 pmax(0, pmin(40, soilDepthInCMmean - 60)) / vanGenuchtenN60_100mean +
                                     #                                 pmax(0, pmin(100, soilDepthInCMmean - 100)) / vanGenuchtenN60_100mean),
                                     thetaR = 1 / soilDepthInCMmean * (pmax(0, pmin(5, soilDepthInCMmean)) * vanGenuchtenThetaR0_5mean +
                                                                       pmax(0, pmin(10, soilDepthInCMmean - 5)) * vanGenuchtenThetaR5_15mean + 
                                                                       pmax(0, pmin(15, soilDepthInCMmean - 15)) * vanGenuchtenThetaR15_30mean + 
                                                                       pmax(0, pmin(30, soilDepthInCMmean - 30)) * vanGenuchtenThetaR30_60mean + 
                                                                       pmax(0, pmin(40, soilDepthInCMmean - 60)) * vanGenuchtenThetaR60_100mean + 
                                                                       pmax(0, pmin(100, soilDepthInCMmean - 100)) * vanGenuchtenThetaR60_100mean),
                                     # calculating plant available water as thetaS - thetaR increases error as the difference of the weighted
                                     # sums isn't quite the same as the sum of the weighted differences
                                     #thetaS = 1 / soilDepthInCMmean * (pmax(0, pmin(5, soilDepthInCMmean)) * vanGenuchtenThetaS0_5mean +
                                     #                                  pmax(0, pmin(10, soilDepthInCMmean - 5)) * vanGenuchtenThetaS5_15mean + 
                                     #                                  pmax(0, pmin(15, soilDepthInCMmean - 15)) * vanGenuchtenThetaS15_30mean + 
                                     #                                  pmax(0, pmin(30, soilDepthInCMmean - 30)) * vanGenuchtenThetaS30_60mean + 
                                     #                                  pmax(0, pmin(40, soilDepthInCMmean - 60)) * vanGenuchtenThetaS60_100mean + 
                                     #                                  pmax(0, pmin(100, soilDepthInCMmean - 100)) * vanGenuchtenThetaS60_100mean),
                                     plantAccessibleWater = 1 / soilDepthInCMmean * (pmax(0, pmin(5, soilDepthInCMmean)) * (vanGenuchtenThetaS0_5mean - vanGenuchtenThetaR0_5mean) +
                                                                                     pmax(0, pmin(10, soilDepthInCMmean - 5)) * (vanGenuchtenThetaS5_15mean - vanGenuchtenThetaR5_15mean) +
                                                                                     pmax(0, pmin(15, soilDepthInCMmean - 15)) * (vanGenuchtenThetaS15_30mean - vanGenuchtenThetaR15_30mean) +
                                                                                     pmax(0, pmin(30, soilDepthInCMmean - 30)) * (vanGenuchtenThetaS30_60mean - vanGenuchtenThetaR30_60mean) +
                                                                                     pmax(0, pmin(40, soilDepthInCMmean - 60)) * (vanGenuchtenThetaS60_100mean - vanGenuchtenThetaR60_100mean) +
                                                                                     pmax(0, pmin(100, soilDepthInCMmean - 100)) * (vanGenuchtenThetaS100_200mean - vanGenuchtenThetaR100_200mean))),
                            psiMPa = -10^seq(log10(10), log10(1E-6), length.out = 100)) %>%
    mutate(water = 10 * (1 - rockFraction) * (thetaR * soilDepthInCMmean + pmax(0, pmin(5, soilDepthInCMmean)) * (vanGenuchtenThetaS0_5mean - vanGenuchtenThetaR0_5mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm0_5mean * abs(1000 * psiMPa))^vanGenuchtenN0_5mean)^(1 - 1/vanGenuchtenN0_5mean) +
                                                                           pmax(0, pmin(10, soilDepthInCMmean - 5)) * (vanGenuchtenThetaS5_15mean - vanGenuchtenThetaR5_15mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm5_15mean * abs(1000 * psiMPa))^vanGenuchtenN5_15mean)^(1 - 1/vanGenuchtenN5_15mean) +
                                                                           pmax(0, pmin(15, soilDepthInCMmean - 15)) * (vanGenuchtenThetaS15_30mean - vanGenuchtenThetaR15_30mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm15_30mean * abs(1000 * psiMPa))^vanGenuchtenN15_30mean)^(1 - 1/vanGenuchtenN15_30mean) +
                                                                           pmax(0, pmin(30, soilDepthInCMmean - 30)) * (vanGenuchtenThetaS30_60mean - vanGenuchtenThetaR30_60mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm30_60mean * abs(1000 * psiMPa))^vanGenuchtenN30_60mean)^(1 - 1/vanGenuchtenN30_60mean) +
                                                                           pmax(0, pmin(40, soilDepthInCMmean - 60)) * (vanGenuchtenThetaS60_100mean - vanGenuchtenThetaR60_100mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm60_100mean * abs(1000 * psiMPa))^vanGenuchtenN60_100mean)^(1 - 1/vanGenuchtenN60_100mean) +
                                                                           pmax(0, pmin(100, soilDepthInCMmean - 100)) * (vanGenuchtenThetaS100_200mean - vanGenuchtenThetaR100_200mean) / (1 + (0.0980665 * 10^vanGenuchtenLog10AlphaInCm100_200mean * abs(1000 * psiMPa))^vanGenuchtenN100_200mean)^(1 - 1/vanGenuchtenN100_200mean)), # mm H₂O
           waterSingleLayerLinear = 10 * soilDepthInCMmean * (1 - rockFraction) * (thetaR + plantAccessibleWater / (1 + (alpha * abs(1000 * psiMPa))^n)^(1 - 1/n))) # mm H₂O
  #waterSingleLayerHarmonic = 10 * soilDepthInCMmean * (thetaR + (thetaS - thetaR) / (1 + (alphaHarmonic * abs(psiMPa))^nHarmonic)^(1 - 1/nHarmonic))) # mm H₂O
  return(waterRetention)
}

#climateGridColumnTypes = cols(name = "c", .default = "d")
#climateCells4km = read_csv(file.path(getwd(), 'GIS/iLand/grid 4 km.csv'), col_types = climateGridColumnTypes)
#waterRetention4km = getWaterRetention(climateCells4km)
#climateCells800m = read_csv(file.path(getwd(), 'GIS/iLand/grid 800 m.csv'), col_types = climateGridColumnTypes)
#waterRetention800m = getWaterRetention(climateCells800m)
climateGridColumnTypes = cols(name = "c", name4km = "c", name800m = "c", name400m = "c", name200m = "c", .default = "d")
climateCells100m = read_csv(file.path(getwd(), 'GIS/iLand/grid 100 m.csv'), col_types = climateGridColumnTypes)
waterRetention100m = getWaterRetention(climateCells100m)

ggplot(waterRetention100m) + # slow! around 5 minutes to draw in RStudio, even longer to send to .png if not rendered in RStudio first
  geom_path(aes(x = psiMPa, y = water, group = name), alpha = 0.015, color = "blue2") +
  coord_cartesian(xlim = c(-9, -1.02E-5), ylim = c(10, 1150)) +
  labs(x = bquote(psi*", MPa"), y = bquote("multilayer soil water, mm H"[2]*"O")) +
  scale_x_continuous(breaks = c(-10, -1, -0.1, -0.01, -0.001, -0.0001, -0.00001), 
                     labels = c("-10", "-1", "-0.1", "-0.01", "-0.001", "", "-0.00001"),
                     minor_breaks = c(-5, -2, -0.5, -0.2, -0.05, -0.02, -0.005, -0.002, -0.0005, -0.002, -0.00005, - 0.00002),
                     trans = trans_new(name = "log10negative",
                                       trans = function(x) { return(-log10(-x)) },
                                       inverse = function(x) { return(-10^(-x)) })) +
ggplot(waterRetention100m) +
  geom_path(aes(x = psiMPa, y = waterSingleLayerLinear, group = name), alpha = 0.015, color = "blue2") +
  coord_cartesian(xlim = c(-9, -1.02E-5), ylim = c(10, 1150)) +
  labs(x = bquote(psi*", MPa"), y = bquote("single layer soil water, mm H"[2]*"O")) +
  scale_x_continuous(breaks = c(-10, -1, -0.1, -0.01, -0.001, -0.0001, -0.00001), 
                     labels = c("-10", "-1", "-0.1", "-0.01", "-0.001", "", "-0.00001"),
                     minor_breaks = c(-5, -2, -0.5, -0.2, -0.05, -0.02, -0.005, -0.002, -0.0005, -0.002, -0.00005, - 0.00002),
                     trans = trans_new(name = "log10negative",
                                       trans = function(x) { return(-log10(-x)) },
                                       inverse = function(x) { return(-10^(-x)) })) +
ggplot(waterRetention100m) +
  geom_path(aes(x = psiMPa, y = waterSingleLayerLinear - water, group = name), alpha = 0.015, color = "darkorchid3") +
  coord_cartesian(xlim = c(-9, -1.02E-5)) +
  #coord_cartesian(ylim = c(-4, 4)) +
  labs(x = bquote(psi*", MPa"), y = bquote("single layer error, mm H"[2]*"O")) +
  scale_x_continuous(breaks = c(-10, -1, -0.1, -0.01, -0.001, -0.0001, -0.00001), 
                     labels = c("-10", "-1", "-0.1", "-0.01", "-0.001", "", "-0.00001"),
                     minor_breaks = c(-5, -2, -0.5, -0.2, -0.05, -0.02, -0.005, -0.002, -0.0005, -0.002, -0.00005, - 0.00002),
                     trans = trans_new(name = "log10negative",
                                       trans = function(x) { return(-log10(-x)) },
                                       inverse = function(x) { return(-10^(-x)) }))
ggsave("iLand/gis/resource unit Mualem-van Genuchten curves.png", width = 25, height = 13, units = "cm")

ggplot(climateCells4km) +
  geom_point(aes(x = vanGenuchtenLog10AlphaInCm0_5mean, y = vanGenuchtenN0_5mean, color = "0-5 cm")) +
  geom_point(aes(x = vanGenuchtenLog10AlphaInCm15_30mean, y = vanGenuchtenN15_30mean, color = "15-30 cm")) +
  geom_point(aes(x = vanGenuchtenLog10AlphaInCm30_60mean, y = vanGenuchtenN30_60mean, color = "30-60 cm")) +
  geom_point(aes(x = vanGenuchtenLog10AlphaInCm60_100mean, y = vanGenuchtenN60_100mean, color = "60-100 cm")) +
  geom_point(aes(x = vanGenuchtenLog10AlphaInCm100_200mean, y = vanGenuchtenN100_200mean, color = "100-200 cm")) +
  labs(x = bquote(alpha*", log"[10]*" cm"), y = "n", color = NULL) +
  scale_color_discrete(breaks = c("0-5 cm", "15-30 cm", "30-60 cm", "60-100 cm", "100-200 cm")) +
  theme(legend.justification = c(0, 1), legend.position = c(0.02, 1)) +
  ggplot(climateCells4km) +
  #geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "longdash") + # not visible
  geom_point(aes(x = vanGenuchtenThetaR0_5mean, y = vanGenuchtenThetaS0_5mean, color = "0-5 cm")) +
  geom_point(aes(x = vanGenuchtenThetaR15_30mean, y = vanGenuchtenThetaS15_30mean, color = "15-30 cm")) +
  geom_point(aes(x = vanGenuchtenThetaR30_60mean, y = vanGenuchtenThetaS30_60mean, color = "30-60 cm")) +
  geom_point(aes(x = vanGenuchtenThetaR60_100mean, y = vanGenuchtenThetaS60_100mean, color = "60-100 cm")) +
  geom_point(aes(x = vanGenuchtenThetaR100_200mean, y = vanGenuchtenThetaS100_200mean, color = "100-200 cm")) +
  labs(x = bquote(theta[r]), y = bquote(theta[s]), color = NULL) +
  scale_color_discrete(breaks = c("0-5 cm", "15-30 cm", "30-60 cm", "60-100 cm", "100-200 cm")) +
  theme(legend.position = "none")


## download POLARIS 1° WGS84 soil tiles for Oregon and Washington with up to 1° buffer in United States: ~364 GB
# There is a 14th directory, vrt, with top level virtual rasters that's unlikely to be of interest for regional downloads.
errors = list()
for (polarisVariable in c("theta_s"))
{
  for (polarisStatistic in c("mean", "mode", "p5", "p50", "p95"))
  {
    for (polarisDepth in c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200"))
    {
      for (minLatitude in 41:48) # 41-49 °N: Oregon-California border at 42 °N, Washington-Canada at 49 °N
      {
        for (minLongitude in -125:-117) # 125-116 °W: Pacific Coast to Idaho
        {
          localFilePath = file.path(getwd(), paste0("GIS/POLARIS/", polarisVariable, " ", polarisStatistic, " N", minLatitude, "-", minLatitude + 1, " W", abs(minLongitude), "-", abs(minLongitude + 1), " ", polarisDepth, " cm.tif"))
          if (file.exists(localFilePath))
          {
            next # assume file downloaded successfully and don't reattempt download
          }
          # /POLARIS/PROPERTIES/v1.0/theta_s/p5/5_15/lat4243_lon-120-119.tif since it is 403
          tryCatch(
            {
              polarisFileName = paste0("lat", minLatitude, minLatitude + 1, "_lon-", abs(minLongitude), "-", abs(minLongitude + 1), ".tif")
              download.file(paste0("http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/", polarisVariable, "/", polarisStatistic, "/", str_replace(polarisDepth, "-", "_"), "/", polarisFileName), localFilePath, mode = "wb")
            }, 
            error = function(e)
            {
              errors = append(errors, e)
            })
        }
      }
    }
  }
}
errors