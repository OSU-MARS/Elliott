## investigation of choices in standard deviation of weighted residuals
# gsl_nls() and nls() both calculate σ
# See https://www.jchau.org/2021/07/12/asymptotic-confidence-intervals-for-nls-regression-in-r/ for background.
# Sibbesen test examples
alruHeightGslNls = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044), weights = pmin(DBH^-2, 1))
alruHeightGslNls10 = gsl_nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044), weights = pmin(10*DBH^-2, 10))
alruHeightNls = nls(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044), weights = pmin(DBH^-2, 1))
alruHeightNlrob = nlrob(TotalHt ~ 1.37 + (a1 + a1p * isPlantation)*DBH^((b1 + b1p * isPlantation)*DBH^(b2 + b2p * isPlantation)), alru2016, start = list(a1 = 0.467, a1p = 0.187, b1 = 1.69, b1p = -0.33, b2 = -0.14, b2p = 0.044), weights = pmin(DBH^-2, 1))

tibble(gsl_nls_sigma = sigma(alruHeightGslNls), gsl_nls_summary = summary(alruHeightGslNls)$sigma, gsl_nls_asymptotic = sqrt(sum(alruHeightGslNls$weights * residuals(alruHeightGslNls)^2) / df.residual(alruHeightGslNls)), 
       gsl_nls10_sigma = sigma(alruHeightGslNls10), gsl_nls10_summary = summary(alruHeightGslNls10)$sigma, gsl_nls10_asymptotic = sqrt(sum(alruHeightGslNls10$weights * residuals(alruHeightGslNls10)^2) / df.residual(alruHeightGslNls10)))

tibble(nls_sigma = sigma(alruHeightNls), nls_asymptotic = sqrt(sum(alruHeightNls$weights * residuals(alruHeightNls)^2) / df.residual(alruHeightNls)))

tibble(nlrob_sigma = sigma(alruHeightNlrob), nlrob_asymptotic = sqrt(sum(alru2016$DBH^-2 * residuals(alruHeightNlrob)^2) / df.residual(alruHeightNlrob))) # sigma(nlrob) doesn't consider weights

summary(alruHeightGslNls)$cov.unscaled
vcov(alruHeightGslNls)
sum(alruHeightGslNls$weights * residuals(alruHeightGslNls)^2)
sum(alruHeightGslNls$weights * residuals(alruHeightGslNls)^2) / df.residual(alruHeightGslNls) * solve(t(alruHeightGslNls$m$gradient()) %*% alruHeightGslNls$m$gradient())

summary(alruHeightGslNls10)$cov.unscaled
vcov(alruHeightGslNls10)
sum(alruHeightGslNls10$weights * residuals(alruHeightGslNls10)^2)
sum(alruHeightGslNls10$weights * residuals(alruHeightGslNls10)^2) / df.residual(alruHeightGslNls10) * solve(t(alruHeightGslNls10$m$gradient()) %*% alruHeightGslNls10$m$gradient())

vcov(alruHeightNls)
sum(residuals(alruHeightNls)^2) / df.residual(alruHeightNls) * solve(t(alruHeightNls$m$gradient()) %*% alruHeightNls$m$gradient())

vcov(alruHeightNlrob)
sum(alruHeightNlrob$weights * alruHeightNlrob$rweights * residuals(alruHeightNlrob)^2) / df.residual(alruHeightNlrob) * solve(t(alruHeightNlrob$m$gradient()) %*% alruHeightNlrob$m$gradient())

bind_cols(name = names(alruHeightGslNls$m$getPars()), gsl_nls = confint2(alruHeightGslNls), gsl_nls10 = confint2(alruHeightGslNls10), nls = confint2(alruHeightNls), nlrob = alruHeightNlrob$m$getPars() + sqrt(diag(sum(alruHeightNlrob$weights * alruHeightNlrob$rweights * residuals(alruHeightNlrob)^2) / df.residual(alruHeightNlrob) * solve(t(alruHeightNlrob$m$gradient()) %*% alruHeightNlrob$m$gradient()))) %o% qt(c(0.025, 0.975), df.residual(alruHeightNlrob))) %>%
  mutate(gsl_nls_min = gsl_nls[, 1], gsl_nls_max = gsl_nls[, 2],
         gsl_nls10_min = gsl_nls10[, 1], gsl_nls10_max = gsl_nls10[, 2],
         nls_min = nls[, 1], nls_max = nls[, 2],
         nlrob_min = nlrob[, 1], nlrob_max = nlrob[, 2]) %>%
  select(-gsl_nls, -gsl_nls10, -nls, -nlrob)


## height residual weighting: all use weighted squared residuals
#   gsl_nls(): wts = weights, .Call(wts), https://github.com/JorisChau/gslnls/blob/master/R/nls.R
#   nls(): w = weights -> w = w^0.5, https://github.com/cran/nls/blob/master/R/nls.R, https://www.gnu.org/software/gsl/doc/html/nls.html#weighted-nonlinear-least-squares
#   nlrob(M estimator): w = psi(resid/Scale) * weights, ._nlrob.w = w, ._nlrob.w = , nls(weights = ._nlrob.w)
alruHeightGslNlsResidualConstPower = gsl_nls(residual ~ a0 + a1*DBH^b1, tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))), start = list(a0 = 0, a1 = 1, b1 = 1))
alruHeightGslNlsResidualPower = gsl_nls(residual ~ a1*DBH^b1, tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))), start = list(a1 = 1, b1 = 1))
alruHeightResidualLm01 = lm(residual ~ 0 + DBH, tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))))
alruHeightResidualLm0.4 = lm(residual ~ I(DBH^0.4), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))))
alruHeightResidualLm0.5 = lm(residual ~ I(sqrt(DBH)), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))))
alruHeightResidualLm1 = lm(residual ~ DBH, tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))))
alruHeightResidualLm2 = lm(residual ~ DBH + I(DBH^2), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls)))) # DBH^2 not significant (p = 0.11)
alruHeightResidualLm3 = lm(residual ~ DBH + I(DBH^2) + I(DBH^3), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightGslNls))))
alruHeightNlrobResidualConstPower = nlrob(residual ~ a0 + a1*DBH^b1, tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightNlrob))), start = list(a0 = 0, a1 = 1, b1 = 1))
alruHeightNlrobResidualPower = nlrob(residual ~ a1*DBH^(b1 + b1p * isPlantation), tibble(DBH = alru2016$DBH, isPlantation = alru2016$isPlantation, residual = abs(residuals(alruHeightNlrob))), start = list(a1 = 1, b1 = 1, b1p = 0)) # a1p, b1p not mutually significant
#confint_nlrob(alruHeightNlrobResidualPower, level = 0.99, weights = rep(1, alruHeightNlrobResidualPower$nobs))

tribble(~name, ~mae, ~aic,
        "gsl_nls const", mean(abs(residuals(alruHeightGslNlsResidualConstPower))), AIC(alruHeightGslNlsResidualConstPower),
        "gsl_nls", mean(abs(residuals(alruHeightGslNlsResidualPower))), AIC(alruHeightGslNlsResidualPower),
        "nlrob const", mean(abs(residuals(alruHeightNlrobResidualConstPower))), AIC(alruHeightNlrobResidualConstPower),
        "nlrob", mean(abs(residuals(alruHeightNlrobResidualPower))), AIC(alruHeightNlrobResidualPower),
        "lm01", mean(abs(residuals(alruHeightResidualLm01))), AIC(alruHeightResidualLm01),
        "lm0.4", mean(abs(residuals(alruHeightResidualLm0.4))), AIC(alruHeightResidualLm0.4),
        "lm0Height.5", mean(abs(residuals(alruHeightResidualLm0.5))), AIC(alruHeightResidualLm0.5),
        "lm1", mean(abs(residuals(alruHeightResidualLm1))), AIC(alruHeightResidualLm1),
        "lm2", mean(abs(residuals(alruHeightResidualLm2))), AIC(alruHeightResidualLm2),
        "lm3", mean(abs(residuals(alruHeightResidualLm3))), AIC(alruHeightResidualLm3)) %>%
  mutate(deltaAic = aic - min(aic)) %>%
  arrange(desc(deltaAic))

tibble(gsl_nls = mean(abs(residuals(alruHeightGslNls))), nlrob = mean(abs(residuals(alruHeightNlrobResidualPower))), lm0.5 = mean(abs(alruHeightResidualLm0.5$residuals)), lm01 = mean(abs(alruHeightResidualLm01$residuals)), lm1 = mean(abs(alruHeightResidualLm1$residuals)), lm2 = mean(abs(alruHeightResidualLm2$residuals)))

ggplot() +
  geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightGslNls)), color = "gsl_nls()"), alpha = 0.1, shape = 16) +
  geom_point(aes(x = alru2016$DBH, y = abs(residuals(alruHeightNlrob)), color = "nlrob()"), alpha = 0.1, shape = 16) +
  geom_smooth(aes(x = alru2016$DBH, y = abs(residuals(alruHeightGslNls)), color = "gsl_nls() GAM", fill = "gsl_nls() GAM"), alpha = 0.1, formula = y ~ s(x, k = 10), method = "gam") +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightGslNlsResidualConstPower), color = "gsl_nls() const"), group = alru2016$isPlantation) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightGslNlsResidualPower), color = "gsl_nls()"), group = alru2016$isPlantation) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightResidualLm0.5), color = "lm0.5"), group = alru2016$isPlantation) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightResidualLm3), color = "lm3"), group = alru2016$isPlantation) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightNlrobResidualConstPower), color = "nlrob() const"), group = alru2016$isPlantation) +
  geom_line(aes(x = alru2016$DBH, y = predict(alruHeightNlrobResidualPower), color = "nlrob()"), group = alru2016$isPlantation) +
  labs(x = "DBH, cm", y = "|height residual|, m", color = NULL, fill = NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.03, 1))
