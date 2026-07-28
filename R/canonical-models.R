# canonical-models.R
#
# Canonical model constructors for the dgtf package.
#
# Each function assembles a dgtf_model from existing component constructors
# defined in components.R. The returned object can be passed to dgtf() for
# inference or dgtf_simulate_model() for data generation.
#
# The baseline intensity (intercept) is controlled through the seasonality
# argument. With seas_none() (default) the level is zero; for a nonzero
# baseline use seas_period(1, init = value). See ?dgtf-seasonality.
#
# C++ engine dispatch (from Model::init_model in src/Model.h):
#   sys_eq = "shift"    -> ftrans = "sliding",   lag type = user-specified
#   sys_eq = "nbinom"   -> ftrans = "iterative", lag type = "nbinom" (forced)
#   sys_eq = "identity" -> ftrans = "sliding",   lag type = "uniform" (forced)
#
# Default system structure by constructor:
#   dgtf_hawkes()          -> sys_shift()    (sliding, lognormal lag, truncated by residual_prob)
#   dgtf_distributed_lag() -> sys_nbinom()   (iterative, NB lag, exact Solow form)
#   dgtf_koyck()           -> sys_nbinom()   (iterative, geometric lag, exact)
#   dgtf_nb_ar()           -> sys_identity() (sliding, uniform lag, AR form)
#   dgtf_poisson_ar()      -> sys_identity() (sliding, uniform lag, AR form)
#
# For the distributed-lag and Koyck constructors the user can override
# with sys = "shift" to get the truncated sliding-window form (controlled
# by residual_prob on the lag component). The truncated form is faster
# per step but introduces a small approximation in the lag tail; the
# default exact (iterative) form is mathematically clean and matches
# MCMC's nL = nT fixup.


#' Discretized Hawkes model
#'
#' Constructs a discretized Hawkes model (branching process with
#' immigration). The conditional intensity is
#' \deqn{\lambda_t = \varphi + h(\psi_t) \sum_{l=1}^{L} \phi_l y_{t-l}}
#' where \eqn{h(\cdot)} is the gain function, \eqn{\phi_l} is the
#' discretized serial interval (lognormal), and \eqn{\psi_t} follows a
#' random walk.
#'
#' The discounted variant (M1 in the CSDA paper) is obtained by fitting
#' this same model with \code{use_discount = TRUE} in the inference
#' control, not through a separate constructor.
#'
#' @param obs Character or \code{dgtf_obs} object. Default \code{"nbinom"}.
#' @param gain Character or \code{dgtf_gain} object. Default
#'   \code{"softplus"}.
#' @param meanlog Numeric. Log-scale mean \eqn{\mu} of the serial
#'   interval. Default 1.386.
#' @param sigma2 Numeric. Log-scale \strong{variance} \eqn{\sigma^2} of
#'   the serial interval. Default 0.32.
#' @param W Numeric. Innovation variance. Default 0.01.
#' @param seasonality A \code{\link[=dgtf-seasonality]{dgtf_seasonality}}
#'   object controlling both the seasonal pattern and the baseline level.
#'   Default \code{seas_none()} (zero baseline, no periodicity). For a
#'   nonzero baseline without periodicity use
#'   \code{seas_period(1, init = value)}.
#'
#' @return A \code{dgtf_model} object.
#' @export
#' @examples
#' # Default Hawkes (zero baseline, softplus gain, lognormal lag)
#' m <- dgtf_hawkes()
#'
#' # With nonzero baseline
#' m <- dgtf_hawkes(seasonality = seas_period(1, init = 5))
#'
#' # COVID-19 application with weekly seasonality
#' m <- dgtf_hawkes(meanlog = 1.386, sigma2 = 0.32,
#'                  seasonality = seas_weekly(init = rep(1, 7)))
dgtf_hawkes <- function(obs = "nbinom",
                        gain = "softplus",
                        meanlog = 1.386,
                        sigma2 = 0.32,
                        W = 0.01,
                        seasonality = seas_none()) {

  dgtf_model(
    obs         = as_component(obs, "obs"),
    link        = link_identity(),
    sys         = sys_shift(),
    gain        = as_component(gain, "gain"),
    lag         = lag_lognormal(meanlog = meanlog, sigma2 = sigma2),
    err         = err_gaussian(W = W),
    seasonality = seasonality
  )
}


#' Distributed lag model
#'
#' Constructs a distributed lag model where the lag weights follow a
#' negative binomial distribution parameterized by \eqn{r} (number of
#' successes) and \eqn{\kappa} (failure probability).
#'
#' By default the model uses the exact iterative Solow form
#' (\code{sys = "nbinom"}), which represents the full untruncated NB
#' lag distribution through an \eqn{r+1}-dimensional latent state and
#' the recursion in \code{\link[=dgtf-sys]{sys_nbinom()}}. Pass
#' \code{sys = "shift"} for the truncated sliding-window form, where
#' the lag tail is dropped at the truncation length implied by the
#' lag component's \code{residual_prob} budget. The two forms are
#' mathematically equivalent in the limit of zero truncation budget.
#'
#' When \eqn{r = 1} the lag distribution is geometric and the model
#' reduces to the Koyck model (see \code{\link{dgtf_koyck}}).
#'
#' @param r Integer. Number of successes in the NB lag distribution.
#'   Default 4. Under \code{sys = "nbinom"} this is structural: it
#'   sets the state dimension \eqn{r+1} and must be held fixed at
#'   construction (priors on \code{par2 = r} are refused by HVA).
#' @param kappa Numeric. Failure probability in (0, 1). Default 0.5.
#'   Can be inferred under either system form.
#' @param sys Either \code{"nbinom"} (default, exact iterative form)
#'   or \code{"shift"} (truncated sliding-window form). The exact
#'   form is recommended for fitting and matches MCMC's internal
#'   representation; the truncated form is mainly for users who want
#'   the same residual-probability truncation as
#'   \code{\link{dgtf_hawkes}}.
#' @param obs Character or \code{dgtf_obs} object. Default \code{"nbinom"}.
#' @param gain Character or \code{dgtf_gain} object. Default
#'   \code{"softplus"}.
#' @param W Numeric. Innovation variance. Default 0.01.
#' @param seasonality A \code{dgtf_seasonality} object. Default
#'   \code{seas_none()}.
#'
#' @return A \code{dgtf_model} object.
#' @export
#' @examples
#' m <- dgtf_distributed_lag(r = 4, kappa = 0.5)
#'
#' # DL(2) with weekly seasonality
#' m <- dgtf_distributed_lag(r = 2, seasonality = seas_weekly())
#'
#' # Truncated sliding-window form (matches earlier behaviour)
#' m <- dgtf_distributed_lag(r = 4, kappa = 0.5, sys = "shift")
dgtf_distributed_lag <- function(r = 4L,
                                 kappa = 0.5,
                                 sys = c("nbinom", "shift"),
                                 obs = "nbinom",
                                 gain = "softplus",
                                 W = 0.01,
                                 seasonality = seas_none()) {

  sys <- match.arg(sys)
  sys_comp <- switch(sys,
                     nbinom = sys_nbinom(kappa = kappa),
                     shift  = sys_shift())

  dgtf_model(
    obs         = as_component(obs, "obs"),
    link        = link_identity(),
    sys         = sys_comp,
    gain        = as_component(gain, "gain"),
    lag         = lag_nbinom(r = r, kappa = kappa),
    err         = err_gaussian(W = W),
    seasonality = seasonality
  )
}


#' Koyck distributed lag model
#'
#' Equivalent to \code{dgtf_distributed_lag(r = 1, kappa = kappa)}.
#' The geometric lag weights produce the well-known Koyck recursion
#' \eqn{f_t = \kappa f_{t-1} + (1 - \kappa) h(\psi_t) y_{t-1}},
#' which is exactly the iterative Solow form at \eqn{r = 1}. The
#' default \code{sys = "nbinom"} represents this directly via a
#' two-dimensional state \eqn{(\psi_{t+1}, f_t)}; pass
#' \code{sys = "shift"} for the truncated sliding equivalent.
#'
#' @param kappa Numeric. Memory parameter in (0, 1). Default 0.5.
#' @param sys,obs,gain,W,seasonality See \code{\link{dgtf_distributed_lag}}.
#'
#' @return A \code{dgtf_model} object.
#' @export
#' @examples
#' m <- dgtf_koyck(kappa = 0.7)
dgtf_koyck <- function(kappa = 0.5,
                       sys = c("nbinom", "shift"),
                       obs = "nbinom",
                       gain = "softplus",
                       W = 0.01,
                       seasonality = seas_none()) {

  sys <- match.arg(sys)
  dgtf_distributed_lag(
    r           = 1L,
    kappa       = kappa,
    sys         = sys,
    obs         = obs,
    gain        = gain,
    W           = W,
    seasonality = seasonality
  )
}


#' Negative binomial AR model with time-varying coefficients
#'
#' Constructs a negative binomial AR(p) model. The intensity is
#' \deqn{\lambda_t = \varphi + \sum_{i=1}^{p} h(\psi_{i,t}) y_{t-i}}
#' where each \eqn{\psi_{i,t}} follows an independent random walk.
#'
#' The C++ engine uses \code{sys_eq = "identity"} (\eqn{G = I_p}) and
#' \code{lag_dist = "uniform"} (\eqn{\phi_l \equiv 1}). With
#' \code{full_rank = TRUE}, all p states receive independent noise.
#'
#' @param p Integer. Autoregressive order. Default 3.
#' @param gain Character or \code{dgtf_gain} object. Default
#'   \code{"softplus"}.
#' @param W Numeric. Scalar innovation variance (applied as
#'   \eqn{W \mathbf{I}_p}). Default 0.01.
#' @param seasonality A \code{dgtf_seasonality} object. Default
#'   \code{seas_none()}.
#'
#' @return A \code{dgtf_model} object.
#' @export
#' @examples
#' m <- dgtf_nb_ar(p = 3, seasonality = seas_weekly())
dgtf_nb_ar <- function(p = 3L,
                       gain = "softplus",
                       W = 0.01,
                       seasonality = seas_none()) {

  p <- as.integer(p)
  if (p < 1L) stop("`p` must be a positive integer.", call. = FALSE)

  dgtf_model(
    obs         = obs_nbinom(),
    link        = link_identity(),
    sys         = sys_identity(),
    gain        = as_component(gain, "gain"),
    lag         = lag_uniform(window = p),
    err         = err_gaussian(W = diag(W, nrow = p), full_rank = (p > 1L)),
    seasonality = seasonality
  )
}


#' Poisson AR model with time-varying coefficients
#'
#' Same structure as \code{\link{dgtf_nb_ar}} but with Poisson
#' observations. Denoted M4 in Tang and Prado (2025).
#'
#' @param p Integer. Autoregressive order. Default 5.
#' @param gain,W,seasonality See \code{\link{dgtf_nb_ar}}.
#'
#' @return A \code{dgtf_model} object.
#' @export
#' @examples
#' m <- dgtf_poisson_ar(p = 5)
dgtf_poisson_ar <- function(p = 5L,
                            gain = "softplus",
                            W = 0.01,
                            seasonality = seas_none()) {

  p <- as.integer(p)
  if (p < 1L) stop("`p` must be a positive integer.", call. = FALSE)

  dgtf_model(
    obs         = obs_poisson(),
    link        = link_identity(),
    sys         = sys_identity(),
    gain        = as_component(gain, "gain"),
    lag         = lag_uniform(window = p),
    err         = err_gaussian(W = diag(W, nrow = p), full_rank = (p > 1L)),
    seasonality = seasonality
  )
}
