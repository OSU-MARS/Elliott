library(dplyr)
library(ggplot2)
library(gslnls)
library(mgcv)
library(robustbase)
library(tibble)

theme_set(theme_bw())

## regression setup
data = tibble(x = runif(1000, 0, 100), 
              y = x * (1 + rnorm(1000, mean = 0, sd = 0.1)))
lmFit = lm(y ~ x, data, weight = x^-1)
gamFit = gam(y ~ x, data = data, weights = x^-1)
gslNlsFit = gsl_nls(y ~ m*x + b, data, start = list(m = 1, b = 0), weights = x^-1)
nlrobFit = nlrob(y ~ m*x + b, data, start = list(m = 1, b = 0), weights = x^-1)
nlsFit = nls(y ~ m*x + b, data, start = list(m = 1, b = 0), weights = x^-1)

correctedNlrobFit = nlrobFit
correctedNlrobFit$weights = data$x^-1

residuals = bind_cols(x = data$x, y = data$y, weights = data$x^-1,
                      lmAct = predict(lmFit) - data$y, lmRes = -residuals(lmFit), lmResid = -lmFit$residuals,
                      gamAct = predict(gamFit) - data$y, gamRes = -residuals(gamFit), gamResid = -gamFit$residuals,
                      gslNlsAct = predict(gslNlsFit) - data$y, gslNlsRes = -residuals(gslNlsFit), gslNlsResid = -gslNlsFit$m$resid(),
                      nlrobAct = predict(nlrobFit) - data$y, nlrobRes = -residuals(nlrobFit), nlrobResid = -nlrobFit$m$resid(),
                      nlsAct = predict(nlsFit) - data$y, nlsRes = -residuals(nlsFit), nlsResid = -nlsFit$m$resid())

## AIC
# https://github.com/wch/r-source/blob/trunk/src/library/stats/R/nls.R
# https://github.com/wch/r-source/blob/trunk/src/library/stats/R/logLik.R
#
# Ingdal M, Johnson R, Harrington DA. 2019. The Akaike information criterion in weighted regression of immittance data. 
#  Electrochimica Acta 317:648-653. https://doi.org/10.1016/j.electacta.2019.06.030
# Lindsey C, Sheather S. 2010. Variable Selection in Linear Regression. The Stata Journal 10(4):650-659. 
#  https://doi.org/10.1177/1536867X1101000
#
residuals %>% 
  mutate(gslNlsStdDev = sqrt(1/n() * sum(weights * gslNlsAct^2)) / sqrt(weights)) %>%
  summarize(lmLogLikG = sum(dnorm(lmAct, sd = sqrt(1/n() * sum(weights * lmAct^2)) / sqrt(weights), log = TRUE)),
            gamLogLikG = sum(dnorm(gamAct, sd = sqrt(1/n() * sum(weights * gamAct^2)) / sqrt(weights), log = TRUE)),
            gslNlsLogLikG = sum(dnorm(gslNlsAct, sd = sqrt(1/n() * sum(weights * gslNlsAct^2)) / sqrt(weights), log = TRUE)),
            #gslNlsLogLikG = logLik(gslNlsFit),
            #gslNlsLogLikFormula = -n()/2 * (log(2 * pi) + 1 - log(n()) - sum(log(weights))/n() + log(sum(weights*gslNlsAct^2))),
            #gslNlsLogLikIngdal = -n()/2 * (log(2*pi) + log(sum(weights * gslNlsAct^2)) - log(n()) + 1),
            #gslNlsLogLikLindsey = -n()/2 * (log(sum(weights * gslNlsAct^2)/n()) + 1 + log(2*pi)),
            nlrobLogLikG = sum(dnorm(nlrobAct, sd = sqrt(1/n() * sum(weights * nlrobAct^2)) / sqrt(weights), log = TRUE)),
            #nlrobLogLikG = -n()/2 * (log(2 * pi) + 1 - log(n()) - sum(log(weights))/n() + log(sum(weights*nlrobAct^2))),
            nlsLogLikG = sum(dnorm(nlsAct, sd = sqrt(1/n() * sum(weights * nlsAct^2)) / sqrt(weights), log = TRUE)),
            lmAic = -2*lmLogLikG + 2 * (length(coef(lmFit)) + 1),
            gamAic = -2*gamLogLikG + 2 * (length(coef(gamFit)) + 1),
            gslNlsAic = -2*gslNlsLogLikG + 2 * (length(coef(gslNlsFit)) + 1),
            nlrobAic = -2*nlrobLogLikG + 2 * (length(coef(nlrobFit)) + 1),
            nlsAic = -2*nlsLogLikG + 2 * (length(coef(nlsFit)) + 1),
            gslNlsLogLikT = sum(dt(gslNlsAct / gslNlsStdDev, df = 8, log = TRUE) - log(gslNlsStdDev)),
            gslNlsAicT = -2*gslNlsLogLikT + 2 * (length(coef(gslNlsFit)) + 1))

## BIC
residuals %>% 
  mutate(gslNlsStdDev = sqrt(1/n() * sum(weights * gslNlsAct^2)) / sqrt(weights)) %>%
  summarize(lmLogLikG = sum(dnorm(lmAct, sd = sqrt(1/n() * sum(weights * lmAct^2)) / sqrt(weights), log = TRUE)),
            gamLogLikG = sum(dnorm(gamAct, sd = sqrt(1/n() * sum(weights * gamAct^2)) / sqrt(weights), log = TRUE)),
            gslNlsLogLikG = sum(dnorm(gslNlsAct, sd = sqrt(1/n() * sum(weights * gslNlsAct^2)) / sqrt(weights), log = TRUE)),
            nlrobLogLikG = sum(dnorm(nlrobAct, sd = sqrt(1/n() * sum(weights * nlrobAct^2)) / sqrt(weights), log = TRUE)),
            nlsLogLikG = sum(dnorm(nlsAct, sd = sqrt(1/n() * sum(weights * nlsAct^2)) / sqrt(weights), log = TRUE)),
            lmBic = -2*lmLogLikG + (length(coef(lmFit)) + 1) * log(n()),
            gamBic = -2*gamLogLikG + (length(coef(gamFit)) + 1) * log(n()),
            gslNlsBic = -2*gslNlsLogLikG + (length(coef(gslNlsFit)) + 1) * log(n()),
            nlrobBic = -2*nlrobLogLikG + (length(coef(nlrobFit)) + 1) * log(n()),
            nlsBic = -2*nlsLogLikG + (length(coef(nlsFit)) + 1) * log(n()))

tribble(~fit, ~m, ~b, ~logLik, ~aic, ~bic,
        "lm", lmFit$coefficients["x"], lmFit$coefficients["(Intercept)"], as.numeric(logLik(lmFit)), AIC(lmFit), BIC(lmFit),
        "gam", gamFit$coefficients["x"], gamFit$coefficients["(Intercept)"], as.numeric(logLik(gamFit)), AIC(gamFit), BIC(gamFit),
        "gsl_nls", gslNlsFit$m$getPars()["m"], gslNlsFit$m$getPars()["b"], as.numeric(logLik(gslNlsFit)), AIC(gslNlsFit), BIC(gslNlsFit),
        "nlrob", nlrobFit$m$getPars()["m"], nlrobFit$m$getPars()["b"], as.numeric(logLik(nlrobFit)), AIC(nlrobFit), BIC(nlrobFit),
        "nlrob partially corrected", nlrobFit$m$getPars()["m"], nlrobFit$m$getPars()["b"], as.numeric(logLik(correctedNlrobFit)), AIC(correctedNlrobFit), BIC(correctedNlrobFit),
        "nls", nlsFit$m$getPars()["m"], nlsFit$m$getPars()["b"], as.numeric(logLik(nlsFit)), AIC(nlsFit), BIC(nlsFit))

# check plots
ggplot() +
  geom_point(aes(x = x, y = y), data, alpha = 0.2, color = "grey25", shape = 16) +
  geom_line(aes(x = data$x, y = predict(lmFit), color = "lm")) +
  geom_line(aes(x = data$x, y = predict(gslNlsFit), color = "gsl_nls")) +
  geom_line(aes(x = data$x, y = predict(nlrobFit), color = "nlrob")) +
  geom_line(aes(x = data$x, y = predict(nlsFit), color = "nls")) +
  labs(x = "x", y = "y", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.98))

ggplot() +
  geom_abline(color = "grey70", linetype = "longdash") +
  #geom_point(aes(x = residuals$lmAct, y = residuals$lmRes, color = "lm"), alpha = 0.2, shape = 16) +
  #geom_point(aes(x = residuals$gslNlsAct, y = residuals$gslNlsRes, color = "gsl_nls"), alpha = 0.2, shape = 16) +
  #geom_point(aes(x = residuals$nlrobAct, y = residuals$nlrobRes, color = "nlrob"), alpha = 0.2, shape = 16) +
  #geom_point(aes(x = residuals$nlsAct, y = residuals$nlsRes, color = "nls"), alpha = 0.2, shape = 16) +
  #labs(x = "unweighted residual", y = "residuals(fit)", color = NULL) +
  geom_point(aes(x = residuals$gslNlsAct, y = residuals$gslNlsResid, color = "gsl_nls"), alpha = 0.2, shape = 16) +
  geom_point(aes(x = residuals$nlrobAct, y = residuals$nlrobResid / sqrt(nlrobFit$rweights), color = "nlrob"), alpha = 0.2, shape = 16) +
  geom_point(aes(x = residuals$nlsAct, y = residuals$nlsResid, color = "nls"), alpha = 0.2, shape = 16) +
  labs(x = "weighted residual", y = "residuals(fit)", color = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 0.98))


# https://stats.stackexchange.com/questions/232263/t-distribution-likelihood
#library(dplyr)
#library(ggplot2)
#library(tidyr)
#logLikelihood = crossing(vect = rnorm(1000, sd=0.05), df = c(seq(1, 9) %o% c(1, 10, 100))) %>%
#  group_by(df) %>%
#  summarize(logLikelihood = sum(log(dt(vect / 0.05, df=df)) - log(0.05)))
#ggplot(logLikelihood) +
#  geom_line(aes(x = df, y = logLikelihood)) +
#  labs(x = "t distribution degrees of freedom", y = "log likelihood") +
#  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500), minor_breaks = c(3, 4, 6, 7, 8, 9, 30, 40, 60, 70, 80, 90, 300, 400, 600, 700, 800, 900))
