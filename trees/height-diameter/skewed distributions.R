# skewed distributions
library(PearsonDS)
pearsonPsmeChapmanForm = pearsonFitML(psmeDiameterFromHeightChapmanForm$residuals) # documented type parameter not actually implemented, so tries all seven forms
pearsonPsmeChapmanRichards = pearsonFitML(psmeDiameterFromHeightChapmanRichards$residuals)
pearsonPsmeRuark = pearsonFitML(psmeDiameterFromHeightRuark$residuals)
pearsonPsmeSibbesenForm = pearsonFitML(psmeDiameterFromHeightSibbesenForm$residuals)
tribble(~type, ~m, ~nu, ~location, ~scale,
        pearsonPsmeChapmanForm$type, pearsonPsmeChapmanForm$m, pearsonPsmeChapmanForm$nu, pearsonPsmeChapmanForm$location, pearsonPsmeChapmanForm$scale,
        pearsonPsmeChapmanRichards$type, pearsonPsmeChapmanRichards$m, pearsonPsmeChapmanRichards$nu, pearsonPsmeChapmanRichards$location, pearsonPsmeChapmanRichards$scale,
        pearsonPsmeRuark$type, pearsonPsmeRuark$m, pearsonPsmeRuark$nu, pearsonPsmeRuark$location, pearsonPsmeRuark$scale,
        pearsonPsmeSibbesenForm$type, pearsonPsmeSibbesenForm$m, pearsonPsmeSibbesenForm$nu, pearsonPsmeSibbesenForm$location, pearsonPsmeSibbesenForm$scale)

pearson = crossing(x = seq(-4, 4, length.out = 100), m = 1.72, nu = c(-0.7, -0.65)) %>% # m > 0 narrows width, nu < 0 is right skew
  group_by(m, nu) %>%
  mutate(iv = dpearsonIV(x, m = m[1], nu = nu[1], location = 0, scale = 1))

ggplot(pearson) +
  geom_line(aes(x = x, y = iv, color = as.factor(nu), group = paste(m, nu))) +
  labs(x = "x", y = bquote("pearsonIV(x, m, nu)"), color = bquote(nu)) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))


library(sn)
skewedTpsmeChapmanForm = selm(residuals ~ 1, data = tibble(residuals = psmeDiameterFromHeightChapmanForm$residuals))
skewedTpsmeChapmanRichards = selm(residuals ~ 1, data = tibble(residuals = psmeDiameterFromHeightChapmanRichards$residuals))
skewedTpsmeRuark = selm(residuals ~ 1, data = tibble(residuals = psmeDiameterFromHeightRuark$residuals))
skewedTpsmeSibbesenForm = selm(residuals ~ 1, data = tibble(residuals = psmeDiameterFromHeightSibbesenForm$residuals))
extractSECdistr(skewedTpsmeChapmanForm) # alpha = 2.21, omega = 20.6, xi = -13.0
extractSECdistr(skewedTpsmeChapmanRichards) # alpha = 2.13, omega = 20.6, xi = -12.8
extractSECdistr(skewedTpsmeRuark) # alpha = 2.24, omega = 20.6, xi = -13.0
extractSECdistr(skewedTpsmeSibbesenForm) # alpha = 2.33, omega = 20.8, xi = -13.2
# ggplot() +
#  geom_histogram(aes(x = psmeDiameterFromHeightRuark$residuals, y = ..density..), binwidth = 2.5)

snt = crossing(x = seq(-4, 4, length.out = 100), df = c(5, 6, 7), alpha = c(1.8, 2.0, 2.2, 2.4), omega = 1) %>% # alpha > 0 is right tailed
  group_by(df, alpha, omega) %>%
  mutate(t = dst(x, alpha = alpha[1], nu = df[1], omega = omega[1]))

ggplot(snt) +
  geom_line(aes(x = x, y = t, color = as.factor(alpha), group = paste(df, alpha))) +
  labs(x = "x", y = bquote("t(x, df, "*alpha*")"), color = bquote(alpha)) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))



library(evd)
evd = crossing(x = seq(-4, 4, length.out = 100), scale = c(0.6, 0.7, 0.8, 0.9, 1.0)) %>%
  group_by(scale) %>% # 
  mutate(gumbel = dgumbel(x, loc = 0, scale = scale[1]))

ggplot(evd) +
  geom_line(aes(x = x, y = gumbel, color = as.factor(scale), group = scale)) +
  labs(x = "x", y = bquote("gumbel(x scale)"), color = bquote(scale)) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))


library(emg)
emg = crossing(x = seq(-4, 4, length.out = 100), sigma = 1, lambda = c(2, 3, 5, 10)) %>%
  group_by(sigma, lambda) %>% # 
  mutate(expNormal = demg(x, mu = 0, sigma = sigma[1], lambda = lambda[1]))

ggplot(emg) +
  geom_line(aes(x = x, y = expNormal, color = as.factor(lambda), group = lambda)) +
  labs(x = "x", y = bquote("emg(x "*sigma*", "*lambda*")"), color = bquote(lambda)) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))


library(skewt)
skewt = crossing(x = seq(-4, 4, length.out = 100), df = 5, gamma = c(1.0, 1.1, 1.2, 1.3)) %>% # gamma < 1 is left skew (min 0), gamma > 1 is right skew (max 2)
  group_by(df, gamma) %>% # dskt() returns a matrix if df and gamma are vectors
  mutate(t = dskt(x, df[1], gamma[1]))

ggplot(skewt) +
  geom_line(aes(x = x, y = t, color = as.factor(gamma), group = gamma)) +
  labs(x = "x", y = bquote("t(x df, "*gamma*")"), color = bquote(gamma)) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))
