## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.width = 6,
  fig.height = 4
)
dgtf_run_long <- identical(Sys.getenv("DGTF_RUN_EXPENSIVE_VIGNETTES"), "true")

## ----compile-install-package, eval = FALSE------------------------------------
# # # Compile and install this package locally
# # Rcpp::compileAttributes()
# # devtools::document()
# # devtools::install()

## ----setup--------------------------------------------------------------------
library(dgtf)
library(ggplot2)

## ----build-model--------------------------------------------------------------
mod <- dgtf_model(
  obs         = obs_nbinom(),
  link        = link_identity(),
  sys         = sys_shift(),
  gain        = gain_softplus(),
  lag         = lag_lognormal(meanlog = 1.35, sigma2 = 0.32),
  err         = err_gaussian(W = 0.01),
  seasonality = seas_period(1, init = 2)
)
mod

## ----simulate-----------------------------------------------------------------
sim <- dgtf_simulate_model(mod, ntime = 210, seed = 3269)
ytrain <- c(sim$y)[1:200]
ytest <- c(sim$y)[201:210]
Rt_train <- log1p(exp(c(sim$psi)))[1:200]
Rt_test <- log1p(exp(c(sim$psi)))[201:210]

print(sim)
plot(sim, what = "y", split = 200)
plot(sim, what = "Rt", split = 200)

## ----prior--------------------------------------------------------------------
prior <- dgtf_prior(
  seasonality = normal(0, 10),
  W           = inv_gamma(1, 1),
  rho         = inv_gamma(1, 1),
  lag         = list(par1 = normal(0, 1), par2 = inv_gamma(1, 1))
)

## ----expensive-vignette-gate, include = FALSE---------------------------------
knitr::opts_chunk$set(eval = dgtf_run_long, purl = dgtf_run_long)

