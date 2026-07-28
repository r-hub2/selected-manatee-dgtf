#' Posterior predictive draws
#'
#' Generates posterior predictive samples for an in-sample series under
#' the fitted model. Use this for posterior-predictive checks and
#' goodness-of-fit assessment.
#'
#' @param object A `dgtf_fit`.
#' @param nrep Number of replicate draws **per posterior sample**. The
#'   total number of posterior-predictive draws is
#'   `nrep * <number of posterior samples>` (number of samples is the
#'   second dimension of `object$fit$psi_stored` / `object$fit$Theta`).
#' @param Rt Optional length-`length(y)` vector of "true" intensity
#'   values used as a reference (e.g. `softplus(sim$psi)` in a
#'   simulation study). Pass `NULL` (default) to skip the reference
#'   comparison.
#' @param ... Unused.
#'
#' @param level Credible-interval level for coverage / width / interval-score
#'   reporting (default `0.95`).
#'
#' @return An object of class `dgtf_ppc` (a list). Always present:
#'   * `chi` — averaged chi-square discrepancy of `y`
#'   * `Rt` — `(ntime + 1) x 3` matrix of `R_t` posterior quantiles at
#'     `(alpha/2, 0.5, 1 - alpha/2)` with `alpha = 1 - level`
#'   * `width_Rt` — mean width of the `R_t` credible bands
#'   * `level` — the credible-interval level used
#'
#'   When `nrep > 0`:
#'   * `crps` — averaged continuous ranked probability score for `y`
#'   * `yhat` — `(ntime + 1) x 3` quantile matrix of posterior-predictive `y`
#'   * `coverage_yhat`, `width_yhat`, `interval_score_yhat` — Winkler interval
#'     score and its components for `y`
#'
#'   When a reference `Rt` is supplied (simulation studies):
#'   * `mae_Rt`, `rmse_Rt`, `coverage_Rt`, `interval_score_Rt`
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y, mod, prior, method = "hva", control = vb_control())
#' pp  <- posterior_predict(fit, nrep = 100)
#' }
posterior_predict <- function(object, ...) UseMethod("posterior_predict")

#' @rdname posterior_predict
#' @export
posterior_predict.dgtf_fit <- function(object,
                                       nrep = 100L,
                                       Rt = NULL,
                                       level = 0.95,
                                       ...) {
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1) {
        stop("`level` must be a single number in (0, 1).", call. = FALSE)
    }

    out <- dgtf_posterior_predictive(
        output     = object$fit,
        model_opts = as_settings(object$model),
        y          = as.numeric(object$y),
        nrep       = as.integer(nrep),
        Rt         = if (is.null(Rt)) NULL else as.numeric(Rt),
        level      = level
    )
    class(out) <- c("dgtf_ppc", "list")
    out
}


#' @export
print.dgtf_ppc <- function(x, digits = 3L, ...) {
    lvl <- x$level %||% 0.95
    cat(sprintf("<dgtf posterior predictive check, level = %.0f%%>\n",
                100 * lvl))
    rows <- list(
        c("chi-sq discrepancy",   x$chi),
        c("CRPS (y)",              x$crps),
        c("coverage (y)",          x$coverage_yhat),
        c("width (y)",             x$width_yhat),
        c("interval score (y)",    x$interval_score_yhat),
        c("MAE (R_t)",             x$mae_Rt),
        c("RMSE (R_t)",            x$rmse_Rt),
        c("coverage (R_t)",        x$coverage_Rt),
        c("width (R_t)",           x$width_Rt),
        c("interval score (R_t)",  x$interval_score_Rt)
    )
    rows <- Filter(function(r) !is.null(r[[2]]), rows)
    tab <- data.frame(
        metric = vapply(rows, `[[`, character(1), 1L),
        value  = vapply(rows, function(r)
            formatC(as.numeric(r[[2]]), digits = digits, format = "f"),
            character(1))
    )
    print(tab, row.names = FALSE)
    invisible(x)
}


#' @export
#' @export
summary.dgtf_ppc <- function(object, ...) {
    lvl <- object$level %||% 0.95

    rows <- list(
        list(target = "y",  metric = "chi-square",     value = object$chi,                  level_applies = FALSE),
        list(target = "y",  metric = "CRPS",           value = object$crps,                 level_applies = FALSE),
        list(target = "y",  metric = "coverage",       value = object$coverage_yhat,        level_applies = TRUE),
        list(target = "y",  metric = "width",          value = object$width_yhat,           level_applies = FALSE),
        list(target = "y",  metric = "interval_score", value = object$interval_score_yhat,  level_applies = FALSE),
        list(target = "Rt", metric = "MAE",            value = object$mae_Rt,               level_applies = FALSE),
        list(target = "Rt", metric = "RMSE",           value = object$rmse_Rt,              level_applies = FALSE),
        list(target = "Rt", metric = "coverage",       value = object$coverage_Rt,          level_applies = TRUE),
        list(target = "Rt", metric = "width",          value = object$width_Rt,             level_applies = FALSE),
        list(target = "Rt", metric = "interval_score", value = object$interval_score_Rt,   level_applies = FALSE)
    )
    rows <- Filter(function(r) !is.null(r$value) && is.finite(r$value), rows)
    nominal <- vapply(rows, function(r) if (r$level_applies) lvl else NA_real_, numeric(1))
    metrics <- data.frame(
        target        = vapply(rows, `[[`, character(1), "target"),
        metric        = vapply(rows, `[[`, character(1), "metric"),
        value         = vapply(rows, `[[`, numeric(1),   "value"),
        nominal_level = nominal,
        row.names     = NULL,
        stringsAsFactors = FALSE
    )

    cal_dev <- list()
    if (!is.null(object$coverage_yhat))
        cal_dev$y  <- object$coverage_yhat - lvl
    if (!is.null(object$coverage_Rt))
        cal_dev$Rt <- object$coverage_Rt - lvl

    out <- list(
        level                 = lvl,
        metrics               = metrics,
        calibration_deviation = cal_dev
    )
    class(out) <- c("summary.dgtf_ppc", "list")
    out
}

#' @export
print.summary.dgtf_ppc <- function(x, digits = 3L, ...) {
    cat(sprintf("<dgtf posterior predictive summary, level = %.0f%%>\n",
                100 * x$level))
    cat("\nScalar metrics:\n")
    m <- x$metrics
    m$value         <- formatC(m$value, digits = digits, format = "g")
    m$nominal_level <- ifelse(is.na(m$nominal_level), "",
                              formatC(m$nominal_level, digits = 2, format = "f"))
    print(m, row.names = FALSE)
    # if (length(x$calibration_deviation)) {
    #     cat("\nCalibration deviation (empirical - nominal):\n")
    #     for (nm in names(x$calibration_deviation))
    #         cat(sprintf("  %-3s : %+7.3f\n", nm,
    #                     x$calibration_deviation[[nm]]))
    # }
    invisible(x)
}


#' Out-of-sample forecasts from a DGTF fit
#'
#' \eqn{h}-step-ahead posterior predictive forecasts, optionally
#' scored against held-out future observations.
#'
#' @param object A `dgtf_fit`.
#' @param h Forecast horizon (positive integer).
#' @param nrep Number of replicate draws per posterior sample. The
#'   total number of forecast draws is `nrep * <number of posterior samples>`.
#' @param ypred_true Optional numeric vector of held-out future
#'   observations (length `>= h`). When supplied, scoring metrics
#'   (`coverage`, `chi`, `crps`) are computed against the last `h`
#'   elements (or the last single value when `only_last = TRUE`).
#'   includes per-step coverage and error metrics.
#' @param only_last If `TRUE`, only the `h`-step-ahead prediction is
#'   reported (intermediate horizons are skipped).
#' @param return_samples If `TRUE`, retains the full posterior
#'   predictive forecast draws in `samples`.
#' @param ... Unused.
#'
#' @return An object of class `dgtf_forecast` (a list):
#'   * `quantiles` — numeric matrix with columns `lower`, `median`,
#'     `upper` at the 95% level. Has `h` rows by default, or `1` row
#'     when `only_last = TRUE`.
#'   * `level` — credible-interval level (currently fixed at `0.95`).
#'   * `horizon` — the requested forecast horizon `h`.
#'   * `only_last`, `nsample` — flags / total forecast-draw count.
#'   * `samples` — matrix of posterior predictive draws (only when
#'     `return_samples = TRUE`).
#'   * `ypred_true` — supplied truth vector (when given).
#'   * `coverage`, `chi`, `crps` — scoring metrics (when `ypred_true`
#'     is supplied).
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y[1:180], mod, prior, method = "hva")
#' fc  <- dgtf_forecast_fit(fit, h = 20, nrep = 200,
#'                          ypred_true = sim$y[181:200])
#' fc
#' summary(fc)
#' }
dgtf_forecast_fit <- function(object,
                              h               = 1L,
                              nrep            = 100L,
                              ypred_true      = NULL,
                              only_last       = FALSE,
                              return_samples  = FALSE,
                              ...) {
    if (!inherits(object, "dgtf_fit"))
        stop("`object` must be a `dgtf_fit`.", call. = FALSE)
    if (!is.numeric(h) || length(h) != 1L || h < 1L)
        stop("`h` must be a positive integer.", call. = FALSE)
    if (!is.numeric(nrep) || length(nrep) != 1L || nrep < 1L)
        stop("`nrep` must be a positive integer.", call. = FALSE)
    if (!is.null(ypred_true)) {
        ypred_true <- as.numeric(ypred_true)
        needed <- if (isTRUE(only_last)) 1L else as.integer(h)
        if (length(ypred_true) < needed) {
            stop(sprintf("`ypred_true` must have at least %d element%s.",
                         needed, if (needed == 1L) "" else "s"), call. = FALSE)
        }
    }

    raw <- dgtf_forecast(
        output                 = object$fit,
        model_settings         = as_settings(object$model),
        y                      = as.numeric(object$y),
        nrep                   = as.integer(nrep),
        k                      = as.integer(h),
        ypred_true             = ypred_true,
        only_evaluate_last_one = isTRUE(only_last),
        return_all_samples     = isTRUE(return_samples)
    )

    nsamp_post <- object$control$nsample %||% object$control$n_sample
    nsample <- if (is.null(nsamp_post)) NA_integer_
               else as.integer(nrep) * as.integer(nsamp_post)

    quantiles <- raw$yqt
    colnames(quantiles) <- c("lower", "median", "upper")
    out <- list(
        quantiles = quantiles,
        level     = 0.95,
        horizon   = as.integer(h),
        only_last = isTRUE(only_last),
        nsample   = nsample
    )
    if (isTRUE(return_samples) && !is.null(raw$ypred_samples))
        out$samples <- raw$ypred_samples
    if (!is.null(raw$ypred_true)) {
        out$ypred_true <- as.numeric(raw$ypred_true)
        out$coverage   <- as.numeric(raw$coverage)
        out$chi        <- as.numeric(raw$chi)
        out$crps       <- as.numeric(raw$crps)
    }

    class(out) <- c("dgtf_forecast", "list")
    out
}

#' @export
print.dgtf_forecast <- function(x, digits = 3L, ...) {
    cat(sprintf("<dgtf forecast: h = %d, level = %.0f%%>\n",
                x$horizon, 100 * x$level))
    if (isTRUE(x$only_last))
        cat("(only the final step is reported)\n")
    if (!is.null(x$nsample) && !is.na(x$nsample))
        cat(sprintf("Forecast draws : %s\n",
                    format(x$nsample, big.mark = ",")))

    q <- x$quantiles
    has_truth <- !is.null(x$ypred_true)
    step_ix <- if (isTRUE(x$only_last)) x$horizon else seq_len(nrow(q))

    df <- data.frame(
        step   = step_ix,
        lower  = formatC(q[, "lower"],  digits = digits, format = "g"),
        median = formatC(q[, "median"], digits = digits, format = "g"),
        upper  = formatC(q[, "upper"],  digits = digits, format = "g"),
        stringsAsFactors = FALSE
    )
    if (has_truth) {
        df$truth <- formatC(x$ypred_true, digits = digits, format = "g")
        df$hit   <- ifelse(x$ypred_true >= q[, "lower"] &
                           x$ypred_true <= q[, "upper"], "*", "")
    }

    cat("\nForecast:\n")
    print(df, row.names = FALSE)

    if (has_truth) {
        cat("\nScoring (vs ypred_true):\n")
        rows <- list(
            c("coverage", x$coverage),
            c("chi-sq",   x$chi),
            c("CRPS",     x$crps)
        )
        rows <- Filter(function(r)
            !is.null(r[[2]]) && is.finite(r[[2]]), rows)
        for (r in rows)
            cat(sprintf("  %-9s : %s\n", r[[1]],
                        formatC(as.numeric(r[[2]]), digits = digits,
                                format = "g")))
    }
    invisible(x)
}

#' @export
summary.dgtf_forecast <- function(object, ...) {
    rows <- list()
    if (!is.null(object$coverage))
        rows[[length(rows) + 1L]] <- list(
            metric = "coverage",   value = object$coverage,
            level_applies = TRUE)
    if (!is.null(object$chi))
        rows[[length(rows) + 1L]] <- list(
            metric = "chi-square", value = object$chi,
            level_applies = FALSE)
    if (!is.null(object$crps))
        rows[[length(rows) + 1L]] <- list(
            metric = "CRPS",       value = object$crps,
            level_applies = FALSE)

    metrics <- if (length(rows)) {
        nominal <- vapply(rows, function(r)
            if (r$level_applies) object$level else NA_real_, numeric(1))
        data.frame(
            metric        = vapply(rows, `[[`, character(1), "metric"),
            value         = vapply(rows, `[[`, numeric(1),   "value"),
            nominal_level = nominal,
            row.names     = NULL,
            stringsAsFactors = FALSE
        )
    } else NULL

    out <- list(
        horizon    = object$horizon,
        only_last  = object$only_last,
        level      = object$level,
        nsample    = object$nsample,
        quantiles  = object$quantiles,
        ypred_true = object$ypred_true,
        metrics    = metrics
    )
    class(out) <- c("summary.dgtf_forecast", "list")
    out
}


#' @export
print.summary.dgtf_forecast <- function(x, digits = 3L, ...) {
    cat(sprintf("<dgtf forecast summary: h = %d, level = %.0f%%>\n",
                x$horizon, 100 * x$level))
    if (!is.null(x$nsample) && !is.na(x$nsample))
        cat(sprintf("Forecast draws : %s\n",
                    format(x$nsample, big.mark = ",")))
    if (isTRUE(x$only_last))
        cat("Only the final step is reported.\n")

    if (!is.null(x$metrics)) {
        cat("\nScoring (vs ypred_true):\n")
        m <- x$metrics
        m$value         <- formatC(m$value, digits = digits, format = "g")
        m$nominal_level <- ifelse(is.na(m$nominal_level), "",
                                  formatC(m$nominal_level, digits = 2,
                                          format = "f"))
        print(m, row.names = FALSE)
    } else {
        cat("\nNo `ypred_true` was supplied; no scoring metrics.\n")
    }
    invisible(x)
}