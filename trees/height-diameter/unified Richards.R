library(dplyr)
library(ggplot2)
library(tidyr)

kU = seq(0.08, 0.12, by = 0.01)
curveFamily = crossing(tibble(x = seq(0, 20, length.out = 100), A = 50, y0 = 0, k = 0.3),
                       d = c(0.85),
                       kU = kU) %>%
   mutate(group = paste0("d = ", d, ", kU = ", kU),
          y = A * (1 + ((y0/A)^(1 - d) - 1) * exp((-kU * x)/d^(d/(1 - d))))^(1/(1 - d)),
          yLogistic = 50 / (1 + exp(-k * (x - 10))))

ggplot(curveFamily) +
  geom_line(aes(x = x, y = yLogistic, linetype = "logistic"), color = "grey25", linetype = "longdash") +
  geom_line(aes(x = x, y = y, color = as.factor(kU), group = group, linetype = "unified Richards")) +
  labs(x = "x", y = "y", color = NULL, linetype = NULL) +
  scale_color_discrete(breaks = kU, labels = sapply(kU, function(kU) { return(bquote(d == 0.85*","~k[U] == .(kU))) })) +
  scale_linetype_manual(breaks = c("unified Richards", "logistic"), values = c("solid", "longdash")) +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.03))
#ggsave("trees/height-diameter/curveFamily.png", units = "cm", width = 9, height = 9, dpi = 90)

dValidity = tibble(d = seq(-4, 4, length.out = 100), divisor = d^(d/(1 - d)))
ggplot(dValidity) +
  geom_line(aes(x = d, y = divisor)) +
  labs(x = "x", y = "unified Richards exponent divisor")
