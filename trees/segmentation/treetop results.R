# assumes library()s, functions, and local maxima from treetops.R setup

## treetop dataset
treetopData = s4268maxima %>% filter(is.na(treetop) == FALSE) %>%
  mutate(x = 0.3048 * x, y = 0.3048 * y, dsmZ = 0.3048 * dsmZ, height = 0.3048 * height, # convert to metric
         treetop = factor(treetop, levels = c("yes", "merge", "noise", "maybe noise", "no")))

# Figure TBD
localMaximaHistogram = treetopData %>% mutate(heightClass = 0.5 * round(height / 0.5)) %>% group_by(heightClass, treetop) %>%
  summarize(maxima = n(), .groups = "drop_last") %>%
  mutate(maximaInHeightClass = sum(maxima))

ggplot() +
  geom_col(aes(x = maxima, y = heightClass, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 1) +
  labs(x = "local maxima above breast height (1.37 m)", y = "height, m", fill = NULL) +
  scale_x_continuous(labels = scales::comma) +
ggplot() +
  geom_col(aes(x = maxima / maximaInHeightClass, y = heightClass, fill = treetop, group = heightClass), localMaximaHistogram, orientation = "y", width = 1) +
  labs(x = "probability", y = NULL, fill = NULL) +
  scale_x_continuous(labels = scales::percent) +
plot_annotation(theme = theme(plot.margin = margin())) +
plot_layout(widths = c(1, 0.4), guides = "collect") &
  scale_fill_manual(breaks = c("yes", "merge", "noise", "maybe noise", "no"), labels = c("single point treetop", "merge point", "residual noise", "processing artifact", "other"), values = c("purple", "green3", "red", "dodgerblue3", "grey90")) &
  scale_y_continuous(breaks = seq(0, 90, by = 10))

