## Component constructors --------------------------------------------------
##
## Each constructor returns a small list with class
##   c("dgtf_<comp>_<type>", "dgtf_<comp>", "dgtf_component")
##
## where `type` is the C++ engine's string identifier. The `params` slot
## carries any continuous-valued component parameters (lag distribution
## hyperparameters, error variance, etc.). The C++ side is fed via
## `as_settings()` (see `R/model.R`).

new_component <- function(kind, type, params = list(), ...) {
    structure(
        c(list(type = type, params = params), list(...)),
        class = c(paste0("dgtf_", kind, "_", type),
                  paste0("dgtf_", kind),
                  "dgtf_component")
    )
}

#' @export
print.dgtf_component <- function(x, ...) {
    cls <- class(x)[1]
    nm <- sub("^dgtf_", "", cls)
    cat(sprintf("<%s>", nm))
    if (length(x$params)) {
        cat(":")
        for (p in names(x$params))
            cat(sprintf(" %s=%s", p, format(x$params[[p]])))
    }
    cat("\n")
    invisible(x)
}


# ---- Observation distribution -------------------------------------------

#' Observation-distribution components
#'
#' Constructors for the observational level of a `dgtf_model()`.
#' Counts: `obs_poisson()`, `obs_nbinom()` (mean/dispersion form),
#' `obs_nbinom_p()` (number-of-successes / probability form).
#' Continuous: `obs_normal()`.
#'
#' @return A `dgtf_obs` component object.
#' @name dgtf-obs
NULL

#' @rdname dgtf-obs
#' @export
obs_poisson <- function() new_component("obs", "poisson")

#' @rdname dgtf-obs
#' @export
obs_nbinom <- function() new_component("obs", "nbinom")

#' @rdname dgtf-obs
#' @export
obs_nbinom_p <- function() new_component("obs", "nbinomp")

#' @rdname dgtf-obs
#' @export
obs_normal <- function() new_component("obs", "gaussian")


# ---- Link function ------------------------------------------------------

#' Link-function components
#'
#' Maps the conditional mean to a linear predictor.
#'
#' @return A `dgtf_link` component object.
#' @name dgtf-link
NULL

#' @rdname dgtf-link
#' @export
link_identity <- function() new_component("link", "identity")

#' @rdname dgtf-link
#' @export
link_exponential <- function() new_component("link", "exponential")

#' @rdname dgtf-link
#' @export
link_logistic <- function() new_component("link", "logistic")


# ---- Gain function ------------------------------------------------------

#' Gain-function components
#'
#' Differentiable transforms applied to the first state element inside the
#' transfer function.
#'
#' @return A `dgtf_gain` component object.
#' @name dgtf-gain
NULL

#' @rdname dgtf-gain
#' @export
gain_identity <- function() new_component("gain", "identity")

#' @rdname dgtf-gain
#' @export
gain_softplus <- function() new_component("gain", "softplus")

#' @rdname dgtf-gain
#' @export
gain_ramp <- function() new_component("gain", "ramp")

#' @rdname dgtf-gain
#' @export
gain_exponential <- function() new_component("gain", "exponential")

#' @rdname dgtf-gain
#' @export
gain_logistic <- function() new_component("gain", "logistic")


# ---- Lag distribution ---------------------------------------------------

#' Lag-distribution components
#'
#' The discretised PMF that weights past observations in a sliding-window
#' transfer function.
#'
#' \describe{
#'   \item{\code{lag_lognormal(meanlog, sigma2, sdlog, residual_prob)}}{Discretised
#'     log-normal distribution. \strong{Note:} \code{sigma2} is the
#'     \emph{variance} on the log scale, i.e. \eqn{\sigma^2} where
#'     \eqn{\log X \sim N(\mu, \sigma^2)}. This matches the C++ engine
#'     and the CSDA paper. The \code{sdlog} argument is also accepted and
#'     converted internally to \code{sigma2 = sdlog^2}.}
#'   \item{\code{lag_nbinom(r, kappa, residual_prob)}}{Discretised
#'     negative-binomial (Solow form), with \code{r} successes and
#'     \code{kappa} failure probability. Reduces to geometric when
#'     \code{r = 1}.}
#'   \item{\code{lag_uniform(window)}}{Discrete uniform over the window
#'     (used for AR models where each lag gets its own time-varying
#'     coefficient).}
#' }
#'
#' @param meanlog Numeric. Log-normal mean parameter \eqn{\mu}.
#' @param sigma2 Numeric. Log-normal \strong{variance} parameter
#'   \eqn{\sigma^2} on the log scale. Defaults to 0.32 when neither
#'   \code{sigma2} nor \code{sdlog} is supplied.
#' @param sdlog Numeric. Log-normal standard deviation parameter on the
#'   log scale. If supplied, it is squared before passing settings to the
#'   C++ engine. Do not supply both \code{sigma2} and \code{sdlog}.
#' @param r,kappa Negative-binomial hyperparameters.
#' @param residual_prob Numeric in \eqn{(0, 1)}. Tail probability that
#'   the discretised lag distribution is allowed to drop when truncating.
#'   The truncation length \eqn{L} is computed from this quantity as the
#'   smallest integer for which \eqn{1 - \sum_{l=1}^L \phi_l \le}
#'   \code{residual_prob}. Default 0.005.
#' @param window Integer. For \code{lag_uniform()}, sets the AR order
#'   \eqn{p}. Accepted as a deprecated argument by \code{lag_lognormal()}
#'   and \code{lag_nbinom()} for backward compatibility but ignored: the
#'   truncation length is derived from \code{residual_prob}.
#'
#' @return A `dgtf_lag` component object.
#' @name dgtf-lag
NULL

#' @rdname dgtf-lag
#' @export
lag_lognormal <- function(meanlog = 1.386, sigma2 = NULL, sdlog = NULL,
                          residual_prob = 0.005, window = NULL) {
    if (!is.null(sigma2) && !is.null(sdlog))
        stop("Supply only one of `sigma2` or `sdlog`.", call. = FALSE)
    if (!is.null(sdlog)) {
        if (sdlog <= 0) stop("`sdlog` must be positive.", call. = FALSE)
        sigma2 <- sdlog^2
    }
    if (is.null(sigma2))
        sigma2 <- 0.32
    if (sigma2 <= 0) stop("`sigma2` must be positive.", call. = FALSE)
    if (residual_prob <= 0 || residual_prob >= 1)
        stop("`residual_prob` must lie strictly in (0, 1).", call. = FALSE)
    if (!is.null(window))
        warning("`window` is deprecated for `lag_lognormal()`; the ",
                "truncation length is computed from `residual_prob`. ",
                "Ignoring supplied `window`.", call. = FALSE)
    new_component("lag", "lognorm",
                  params = list(par1 = meanlog, par2 = sigma2),
                  residual_prob = as.numeric(residual_prob))
}

#' @rdname dgtf-lag
#' @export
lag_nbinom <- function(r = 1, kappa = 0.5,
                       residual_prob = 0.005, window = NULL) {
    if (r <= 0) stop("`r` must be positive.", call. = FALSE)
    if (kappa <= 0 || kappa >= 1)
        stop("`kappa` must lie strictly in (0, 1).", call. = FALSE)
    if (residual_prob <= 0 || residual_prob >= 1)
        stop("`residual_prob` must lie strictly in (0, 1).", call. = FALSE)
    if (!is.null(window))
        warning("`window` is deprecated for `lag_nbinom()`; the ",
                "truncation length is computed from `residual_prob`. ",
                "Ignoring supplied `window`.", call. = FALSE)
    new_component("lag", "nbinomp",
                  params = list(par1 = kappa, par2 = r),
                  residual_prob = as.numeric(residual_prob))
}

#' @rdname dgtf-lag
#' @export
lag_uniform <- function(window = 1L) {
    if (window < 1L) stop("`window` must be a positive integer.", call. = FALSE)
    new_component("lag", "uniform",
                  params = list(),
                  window = as.integer(window))
}

# Truncation length L that the C++ engine derives from the lag PMF, given
# a residual-probability budget. Mirrors LagDist::get_nlag in the engine,
# capped at `max_lag` and floored at `min_lag` so that the R-side display
# matches what the engine actually uses.
#'
#' @keywords internal
lag_truncation_length <- function(lag, max_lag = 50L, min_lag = 1L) {
    if (is.null(lag)) return(NA_integer_)
    if (identical(lag$type, "uniform")) {
        return(as.integer(lag$window %||% 1L))
    }
    rp <- lag$residual_prob
    if (is.null(rp) || !is.finite(rp) || rp <= 0 || rp >= 1)
        return(NA_integer_)
    prob <- 1 - rp
    nlag <- switch(
        lag$type,
        lognorm  = stats::qlnorm(prob, lag$params$par1,
                                 sqrt(lag$params$par2)),
        nbinomp  = stats::qnbinom(prob, size = lag$params$par2,
                                  prob = 1 - lag$params$par1),
        nbinom   = stats::qnbinom(prob, size = lag$params$par2,
                                  prob = 1 - lag$params$par1),
        NA_real_
    )
    if (!is.finite(nlag)) return(NA_integer_)
    # static_cast<unsigned int> in C++ truncates toward zero; use the same
    # convention here so the R-side display matches the engine.
    as.integer(min(max(floor(nlag), min_lag), max_lag))
}


# ---- System equation ----------------------------------------------------

#' System-equation components
#'
#' Determines the form of the latent state evolution \eqn{\theta_t = g_t(\cdot) + w_t}.
#'
#' - `sys_identity()`: \eqn{G = I}, no state evolution (e.g. Poisson AR).
#' - `sys_shift()`: shift matrix appropriate for windowed Hawkes / DL models.
#' - `sys_nbinom()`: nonlinear evolution for the Solow distributed-lag model
#'   (`kappa` is the geometric / NB decay parameter).
#'
#' @param kappa Decay parameter for the nonlinear `sys_nbinom()` evolution.
#'
#' @return A `dgtf_sys` component object.
#' @name dgtf-sys
NULL

#' @rdname dgtf-sys
#' @export
sys_identity <- function() new_component("sys", "identity")

#' @rdname dgtf-sys
#' @export
sys_shift <- function() new_component("sys", "shift")

#' @rdname dgtf-sys
#' @export
sys_nbinom <- function(kappa = 0.5) {
    if (kappa <= 0 || kappa >= 1)
        stop("`kappa` must lie strictly in (0, 1).", call. = FALSE)
    new_component("sys", "nbinom", params = list(kappa = kappa))
}


# ---- System-error distribution ------------------------------------------

#' System-error distribution components
#'
#' - `err_gaussian(W, full_rank, w0)`: Gaussian system error with variance `W`
#'   (scalar or matrix). `full_rank = TRUE` allows correlated state errors.
#' - `err_constant()`: degenerate (zero) system error.
#'
#' @param W Variance (scalar or square matrix).
#' @param full_rank If `TRUE`, treat `W` as a full covariance matrix.
#' @param w0 Initial state variance offset.
#'
#' @return A `dgtf_err` component object.
#' @name dgtf-err
NULL

#' @rdname dgtf-err
#' @export
err_gaussian <- function(W = 0.01, full_rank = FALSE, w0 = 0) {
    Wm <- if (is.matrix(W)) W else matrix(W)
    if (any(diag(Wm) < 0))
        stop("Diagonal of `W` must be non-negative.", call. = FALSE)
    new_component("err", "gaussian",
                  params = list(var = Wm, w0 = w0, full_rank = isTRUE(full_rank)))
}

#' @rdname dgtf-err
#' @export
err_constant <- function() new_component("err", "constant")


# ---- Coercion of string shortcuts ---------------------------------------

#' @keywords internal
as_component <- function(x, kind) {
    if (inherits(x, paste0("dgtf_", kind))) return(x)
    if (is.character(x) && length(x) == 1L) {
        ctor <- switch(
            kind,
            obs  = list(poisson = obs_poisson, nbinom = obs_nbinom,
                        nbinomp = obs_nbinom_p, normal = obs_normal,
                        gaussian = obs_normal),
            link = list(identity = link_identity, exponential = link_exponential,
                        logistic = link_logistic),
            gain = list(identity = gain_identity, softplus = gain_softplus,
                        ramp = gain_ramp, exponential = gain_exponential,
                        logistic = gain_logistic),
            lag  = list(lognorm = lag_lognormal, lognormal = lag_lognormal,
                        nbinom = lag_nbinom, nbinomp = lag_nbinom,
                        uniform = lag_uniform),
            sys  = list(identity = sys_identity, shift = sys_shift,
                        nbinom = sys_nbinom),
            err  = list(gaussian = err_gaussian, normal = err_gaussian,
                        constant = err_constant)
        )
        f <- ctor[[tolower(x)]]
        if (is.null(f))
            stop(sprintf("Unknown %s component: \"%s\".", kind, x),
                 call. = FALSE)
        return(f())
    }
    stop(sprintf("`%s` must be a dgtf_%s object or a recognised string.",
                 kind, kind), call. = FALSE)
}
