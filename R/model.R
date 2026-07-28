#' Build a DGTF model
#'
#' Assembles the components of a dynamic generalised time-series filter
#' (DGTF) model. The returned object captures all structural choices and
#' fixed-parameter values, and is what `dgtf()` and `dgtf_simulate_model()`
#' operate on.
#'
#' Components may be supplied either as full constructor objects
#' (`obs_nbinom()`, `lag_lognormal(...)`, ...) or as the corresponding
#' string shortcuts (`"nbinom"`, `"lognorm"`, ...). The longhand form is
#' what makes the framework extensible: a custom component can be written
#' by following the protocol in `?dgtf-extending` (forthcoming).
#'
#' Anything supplied as a numeric value (e.g. `intercept = 0`,
#' `err = err_gaussian(W = 0.01)`) is treated as **fixed**. To infer it,
#' supply a prior in [`dgtf_prior()`].
#'
#' @section Intercept and seasonality:
#'
#' The seasonality block already encodes the level of the series, so
#' `intercept` and a non-trivial `seasonality` are mutually exclusive:
#'
#' - `dgtf_model(intercept = a)` (with `seasonality` left `NULL`) is
#'   equivalent to `dgtf_model(seasonality = seas_period(1, init = a))`
#'   — the intercept is the single seasonal effect under period 1.
#' - `dgtf_model(seasonality = seas_weekly(init = c(...)))` includes the
#'   level inside the seasonal pattern; do not also set `intercept`.
#'   Setting both produces a warning and the intercept is ignored.
#'
#' @param obs Observation distribution. See [`dgtf-obs`].
#' @param link Link function. See [`dgtf-link`].
#' @param sys System equation. See [`dgtf-sys`].
#' @param gain Gain function. See [`dgtf-gain`].
#' @param lag Lag distribution. See [`dgtf-lag`].
#' @param err System-error distribution. See [`dgtf-err`].
#' @param intercept Numeric scalar baseline. Used **only** when
#'   `seasonality` is `NULL` / `seas_none()`; otherwise the seasonal
#'   effects already carry the level.
#' @param seasonality A [`dgtf_seasonality`][dgtf-seasonality] object,
#'   `NULL` for none, or a single integer period (uses zero seasonal
#'   effects).
#' @param obs_init Numeric vector of initial observation-distribution
#'   parameters (default `c(0, 30)` matches the historical raw-list API).
#' @param zi Logical; whether the observation model is zero-inflated.
#'
#' @return An object of class `dgtf_model`.
#' @export
#' @examples
#' # Period-1 model: intercept = 5
#' dgtf_model(obs_poisson(), link_identity(), sys_identity(), intercept = 5)
#'
#' # Period-2 deterministic seasonality with explicit pattern
#' dgtf_model(obs_poisson(), link_identity(), sys_identity(),
#'            seasonality = seas_period(period = 2, init = c(1, 2)))
#'
#' # Discretised Hawkes with weekly seasonality
#' dgtf_model(
#'   obs_nbinom(), link_identity(), sys_shift(),
#'   gain        = gain_softplus(),
#'   lag         = lag_lognormal(),
#'   seasonality = seas_weekly(init = c(0.9, 1.0, 1.1, 1.0, 0.9, 0.8, 1.3))
#' )
dgtf_model <- function(obs,
                       link,
                       sys,
                       gain        = gain_identity(),
                       lag         = lag_uniform(window = 1L),
                       err         = err_gaussian(W = 0.01),
                       intercept   = 0,
                       seasonality = NULL,
                       obs_init    = c(0, 30),
                       zi          = FALSE) {

    obs  <- as_component(obs,  "obs")
    link <- as_component(link, "link")
    sys  <- as_component(sys,  "sys")
    gain <- as_component(gain, "gain")
    lag  <- as_component(lag,  "lag")
    err  <- as_component(err,  "err")
    seas <- as_seasonality(seasonality)

    intercept_set <- !isTRUE(all.equal(as.numeric(intercept), 0))
    seas_nontrivial <- !is.null(seas) && !is_trivial_seasonality(seas)

    if (intercept_set && seas_nontrivial)
        warning(
            "Both `intercept` and a non-trivial `seasonality` were supplied. ",
            "The seasonal effects already encode the series level, so ",
            "`intercept` is ignored. Bake the level into `seas_*(init = ...)`.",
            call. = FALSE)

    if (is.null(seas))
        seas <- seas_period(period = 1L, init = as.numeric(intercept))

    structure(
        list(
            obs         = obs,
            link        = link,
            sys         = sys,
            gain        = gain,
            lag         = lag,
            err         = err,
            seasonality = seas,
            obs_init    = as.numeric(obs_init),
            zi          = isTRUE(zi)
        ),
        class = "dgtf_model"
    )
}

#' @export
print.dgtf_model <- function(x, ...) {
    cat("<dgtf_model>\n")
    cat(sprintf("  observation : %s\n", x$obs$type))
    cat(sprintf("  link        : %s\n", x$link$type))
    cat(sprintf("  system eq   : %s%s\n", x$sys$type,
                if (length(x$sys$params))
                    paste0(" (", paste(names(x$sys$params), "=",
                                       format(unlist(x$sys$params)),
                                       collapse = ", "), ")") else ""))
    cat(sprintf("  gain        : %s\n", x$gain$type))
    if (!identical(x$lag$type, "uniform") || isTRUE(x$lag$window > 1L)) {
        lag_info <- if (identical(x$lag$type, "uniform")) {
            sprintf("window=%d", as.integer(x$lag$window %||% 1L))
        } else {
            L_hat <- lag_truncation_length(x$lag)
            sprintf("residual=%g, L=%s",
                    x$lag$residual_prob %||% NA_real_,
                    if (is.na(L_hat)) "?" else as.character(L_hat))
        }
        cat(sprintf("  lag dist    : %s (%s)\n", x$lag$type, lag_info))
        if (length(x$lag$params))
            cat(sprintf("                params: %s\n",
                        paste(names(x$lag$params), "=",
                              format(unlist(x$lag$params)),
                              collapse = ", ")))
    }
    cat(sprintf("  sys error   : %s\n", x$err$type))
    if (x$seasonality$period == 1L) {
        cat(sprintf("  intercept   : %s\n",
                    format(x$seasonality$init[1])))
    } else {
        cat(sprintf("  seasonality : period=%d\n",
                    x$seasonality$period))
        used <- x$seasonality$init[seq_len(x$seasonality$period)]
        cat("                init =", paste(format(used), collapse = ", "), "\n")
    }
    if (x$zi) cat("  zero-inflated: TRUE\n")
    invisible(x)
}

# ---- Lowering to the C++ settings list ----------------------------------

#' Convert a `dgtf_model` to the raw settings list expected by the C++ engine
#'
#' Internal helper. Users normally do not call this directly.
#'
#' @keywords internal
#' @export
as_settings <- function(model) UseMethod("as_settings")

#' @export
as_settings.dgtf_model <- function(model) {
    lag_params <- if (length(model$lag$params))
        unname(unlist(model$lag$params)) else numeric(0)

    # AR models: C++ needs param$lag to be a vector of length p.
    # lag_uniform stores params = list() with window = p separately.
    if (identical(model$sys$type, "identity") &&
        identical(model$lag$type, "uniform") &&
        length(lag_params) == 0L) {
        p <- model$lag$window %||% 1L
        lag_params <- rep(0, p)
    }

    # Residual-probability threshold for discretising lognormal / nbinom
    # lags. The C++ engine takes the complementary `lag_prob_thres`, so the
    # truncation length L is the smallest integer with CDF >= lag_prob_thres.
    param_list <- list(
        obs = model$obs_init,
        lag = lag_params,
        err = list(
            full_rank = isTRUE(model$err$params$full_rank),
            var       = model$err$params$var %||% matrix(0.01),
            w0        = model$err$params$w0  %||% 0
        )
    )
    if (!is.null(model$lag$residual_prob)) {
        param_list$lag_prob_thres <- 1 - as.numeric(model$lag$residual_prob)
    }

    list(
        model = list(
            obs_dist  = model$obs$type,
            link_func = model$link$type,
            sys_eq    = model$sys$type,
            gain_func = model$gain$type,
            lag_dist  = model$lag$type,
            err_dist  = model$err$type
        ),
        param = param_list,
        season = list(
            period   = model$seasonality$period,
            in_state = model$seasonality$in_state,
            init     = model$seasonality$init
        ),
        zero = list(inflated = model$zi)
    )
}

# Tiny null-coalesce
`%||%` <- function(a, b) if (is.null(a)) b else a
