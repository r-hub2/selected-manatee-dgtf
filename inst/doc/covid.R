## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
dgtf_run_long <- identical(Sys.getenv("DGTF_RUN_EXPENSIVE_VIGNETTES"), "true")

## ----libraries----------------------------------------------------------------
library(dgtf)
library(ggplot2)

## ----eval = FALSE-------------------------------------------------------------
# ytrain <- covid_scz_2020[1:250]
# ytest <- covid_scz_2020[251:264]
# ggplot(data = data.frame(x = 1:264, y = c(ytrain, ytest)), aes(x = x, y = y)) +
#   geom_line() +
#   geom_vline(xintercept = 250, linetype = "dashed") +
#   theme_minimal()

## ----eval = FALSE-------------------------------------------------------------
# set.seed(123)
# mod_DH <- dgtf_model(
#   obs = obs_nbinom(),
#   link = link_identity(),
#   sys = sys_shift(),
#   gain = gain_softplus(),
#   lag = lag_lognormal(
#     meanlog = rnorm(1, 1, 1),
#     sigma2 = 1 / rgamma(1, 10, 1)
#   ),
#   err = err_gaussian(W = 1 / rgamma(1, 10, 1)),
#   seasonality = seas_period(7, init = runif(7, 0, 10))
# )
# 
# mod_DL <- dgtf_distributed_lag(
#   r = 2,
#   kappa = runif(1, 0.1, 0.9),
#   seasonality = seas_weekly(init = runif(7, 0, 10))
# )
# 
# mod_AR <- dgtf_nb_ar(
#   p = 3,
#   W = 1 / rgamma(1, 10, 1),
#   seasonality = seas_weekly(init = runif(7, 0, 10))
# )

## ----eval = FALSE-------------------------------------------------------------
# prior_DH <- dgtf_prior(
#   seasonality = normal(0, 10),
#   W           = inv_gamma(1, 1),
#   rho         = inv_gamma(1, 1),
#   lag         = list(par1 = normal(0, 1), par2 = inv_gamma(1, 1))
# )
# 
# prior_DL <- dgtf_prior(
#   seasonality = normal(0, 10),
#   W           = inv_gamma(1, 1),
#   rho         = inv_gamma(1, 1),
#   lag         = list(par1 = dist_beta(1, 1))
# )
# 
# prior_AR <- dgtf_prior(
#   seasonality = normal(0, 10),
#   W           = inv_gamma(1, 1),
#   rho         = inv_gamma(1, 1)
# )

## ----expensive-vignette-gate, include = FALSE---------------------------------
knitr::opts_chunk$set(eval = dgtf_run_long, purl = dgtf_run_long)

