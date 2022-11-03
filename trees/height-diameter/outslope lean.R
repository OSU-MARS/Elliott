library(dplyr)
library(ggplot2)
library(tidyr)

theme_set(theme_bw())

outslope = crossing(height = seq(0, 80, length.out = 100),
                    slope = seq(0, 50, by = 10)) %>% # degrees
              mutate(lean = 3, # degrees
                     slopeDrop = height * sin(pi/180 * lean) * tan(pi/180 * slope))

ggplot(outslope) +
  geom_line(aes(x = height, y = slopeDrop, color = slope, group = slope)) +
  labs(x = "tree height, m", y = "vertical distance of crown projection point below base, m") +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.99))
