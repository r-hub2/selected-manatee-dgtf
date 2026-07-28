#' Side-by-side comparison of fitted DGTF models, PPC objects, or forecasts
#'
#' Compares an arbitrary named list of \code{dgtf_fit}, \code{dgtf_ppc},
#' or \code{dgtf_forecast} objects on a common set of metrics. The list
#' must be homogeneous: all fits, all PPCs, or all forecasts.
#'
#' For fit lists, the columns are the in-sample posterior-predictive
#' scoring metrics computed via \code{posterior_predict()}, the static
#' parameter count, the HVA log-marginal-likelihood proxy (NA for MCMC),
#' and the runtime breakdown.
#'
#' For PPC lists, the columns are the in-sample scoring metrics on the
#' observed counts and, when the PPC carries a reference $R_t$
#' trajectory, the corresponding $R_t$ recovery metrics.
#'
#' For forecast lists, the columns are the forecast horizon, the number
#' of forecast draws, and the out-of-sample scoring metrics
#' (\code{coverage}, \code{chi}, \code{crps}) when held-out truth was
#' supplied to \code{dgtf_forecast_fit()}.
#'
#' @param fits Named list of objects of homogeneous class.
#' @param ... For fit lists, optional arguments include \code{level},
#'   \code{nrep}, \code{Rt_truth}, \code{truth}, and
#'   \code{force_recompute}; these control on-the-fly posterior predictive
#'   computation and optional truth comparisons. Extra arguments are ignored
#'   for PPC and forecast lists.
#'
#' @return A data frame of class \code{dgtf_comparison}, one row per
#'   input element. The columns depend on the input class.
#'
#' @export
dgtf_compare <- function(fits, ...) UseMethod("dgtf_compare")

#' @export
dgtf_compare.default <- function(fits, ...) {
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits))))
        stop("`fits` must be a named, non-empty list.", call. = FALSE)

    cls <- vapply(fits, function(z) class(z)[1L], character(1L))
    rep_class <- unique(cls)
    if (length(rep_class) > 1L)
        stop("`fits` must be a homogeneous list. Got mixed classes: ",
             paste(rep_class, collapse = ", "), ".", call. = FALSE)

    elt_class <- rep_class[[1L]]
    if (elt_class == "dgtf_fit")
        return(dgtf_compare.dgtf_fit_list(
            structure(fits, class = "dgtf_fit_list"), ...))
    if (elt_class == "dgtf_ppc")
        return(dgtf_compare.dgtf_ppc_list(
            structure(fits, class = "dgtf_ppc_list"), ...))
    if (elt_class == "dgtf_forecast")
        return(dgtf_compare.dgtf_forecast_list(
            structure(fits, class = "dgtf_forecast_list"), ...))
    stop("Elements must be `dgtf_fit`, `dgtf_ppc`, or `dgtf_forecast` ",
         "objects (got: ", elt_class, ").", call. = FALSE)
}

#' @export
dgtf_compare.dgtf_fit_list <- function(fits,
                                       level = 0.95,
                                       nrep = 100L,
                                       Rt_truth = NULL,
                                       truth = NULL,
                                       force_recompute = FALSE,
                                       ...) {
    extras <- list(...)
    if (length(extras))
        warning("Extra arguments to `dgtf_compare` are ignored for fit lists: ",
                paste(names(extras), collapse = ", "), call. = FALSE)

    fits <- unclass(fits)
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1)
        stop("`level` must be a single number in (0, 1).", call. = FALSE)

    rows <- lapply(seq_along(fits), function(i) {
        fit  <- fits[[i]]
        name <- names(fits)[[i]]
        .dgtf_compare_row(fit, name = name, level = level,
                          nrep = nrep, Rt_truth = Rt_truth,
                          truth = truth,
                          force_recompute = force_recompute)
    })

    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    class(df) <- c("dgtf_comparison", "data.frame")
    attr(df, "input_class") <- "dgtf_fit"
    attr(df, "level")       <- level
    attr(df, "has_Rt_truth") <- !is.null(Rt_truth)
    attr(df, "has_truth")   <- !is.null(truth)
    df
}

#' @export
dgtf_compare.dgtf_ppc_list <- function(fits, ...) {
    extras <- list(...)
    if (length(extras))
        warning("Extra arguments to `dgtf_compare` are ignored for PPC lists: ",
                paste(names(extras), collapse = ", "), call. = FALSE)
    fits <- unclass(fits)

    rows <- lapply(seq_along(fits), function(i) {
        ppc  <- fits[[i]]
        name <- names(fits)[[i]]
        out  <- data.frame(
            model            = name,
            crps_y           = .pull(ppc$crps),
            chi_y            = .pull(ppc$chi),
            coverage_y       = .pull(ppc$coverage_yhat),
            width_y          = .pull(ppc$width_yhat),
            interval_score_y = .pull(ppc$interval_score_yhat),
            stringsAsFactors = FALSE
        )
        if (!is.null(ppc$mae_Rt)) {
            out$mae_Rt            <- .pull(ppc$mae_Rt)
            out$rmse_Rt           <- .pull(ppc$rmse_Rt)
            out$coverage_Rt       <- .pull(ppc$coverage_Rt)
            out$width_Rt          <- .pull(ppc$width_Rt)
            out$interval_score_Rt <- .pull(ppc$interval_score_Rt)
        }
        out
    })
    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    class(df) <- c("dgtf_comparison", "data.frame")
    attr(df, "input_class") <- "dgtf_ppc"
    attr(df, "level") <- fits[[1L]]$level %||% 0.95
    attr(df, "has_Rt_truth") <- !is.null(fits[[1L]]$mae_Rt)
    df
}

#' @export
dgtf_compare.dgtf_forecast_list <- function(fits, ...) {
    extras <- list(...)
    if (length(extras))
        warning("Extra arguments to `dgtf_compare` are ignored for forecast lists: ",
                paste(names(extras), collapse = ", "), call. = FALSE)
    fits <- unclass(fits)

    rows <- lapply(seq_along(fits), function(i) {
        fc   <- fits[[i]]
        name <- names(fits)[[i]]
        out  <- data.frame(
            model    = name,
            horizon  = as.integer(fc$horizon %||% NA_integer_),
            n_draws  = as.integer(fc$nsample %||% NA_integer_),
            stringsAsFactors = FALSE
        )
        if (!is.null(fc$coverage)) {
            out$coverage <- .pull(fc$coverage)
            out$chi      <- .pull(fc$chi)
            out$crps     <- .pull(fc$crps)
        }
        out
    })
    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    class(df) <- c("dgtf_comparison", "data.frame")
    attr(df, "input_class") <- "dgtf_forecast"
    attr(df, "level") <- fits[[1L]]$level %||% 0.95
    attr(df, "has_truth_oos") <- !is.null(fits[[1L]]$coverage)
    df
}

# Helper, used by all three methods.
.pull <- function(x, default = NA_real_)
    if (is.null(x)) default else as.numeric(x)

# Internal: extract one row of the comparison table.
.dgtf_compare_row <- function(fit, name, level, nrep,
                              Rt_truth, truth, force_recompute) {

    # 1. Static-parameter count.
    d <- tryCatch(.dgtf_static_draws(fit), error = function(e) NULL)
    n_param <- if (!is.null(d) && !is.null(d$draws_matrix))
                   as.integer(ncol(d$draws_matrix))
               else NA_integer_

    # 2. HVA-only model-evidence proxy.
    log_marglik <- NA_real_
    if (identical(fit$method, "hva") && !is.null(fit$fit$marglik)) {
        cv <- .dgtf_marglik_convergence(fit$fit$marglik)
        if (!is.null(cv)) log_marglik <- cv$tail_median
    }

    # 3. PPC: cached when safe, fresh otherwise. We always recompute
    #    when the caller passes `Rt_truth` (the cached PPC was likely
    #    built without it) or when force_recompute is set.
    use_cache <- !isTRUE(force_recompute) && is.null(Rt_truth)
    ppc <- if (use_cache) attr(fit, "ppc") else NULL
    if (is.null(ppc) || !inherits(ppc, "dgtf_ppc")) {
        ppc <- tryCatch(
            posterior_predict(fit, nrep = nrep, level = level,
                              Rt = Rt_truth),
            error = function(e) {
                warning(sprintf(
                    "posterior_predict() failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
            })
    }

    pull <- function(x, default = NA_real_)
        if (is.null(x)) default else as.numeric(x)
    secs <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)

    out <- data.frame(
        model              = name,
        n_param            = n_param,
        log_marglik        = log_marglik,
        crps_y             = pull(ppc$crps),
        chi_y              = pull(ppc$chi),
        coverage_y         = pull(ppc$coverage_yhat),
        width_y            = pull(ppc$width_yhat),
        interval_score_y   = pull(ppc$interval_score_yhat),
        elapsed_total      = secs(fit$elapsed),
        elapsed_optim      = secs(fit$elapsed_optimization),
        elapsed_sample     = secs(fit$elapsed_sampling),
        stringsAsFactors   = FALSE
    )
    if (!is.null(Rt_truth)) {
        out$mae_Rt            <- pull(ppc$mae_Rt)
        out$rmse_Rt           <- pull(ppc$rmse_Rt)
        out$coverage_Rt       <- pull(ppc$coverage_Rt)
        out$interval_score_Rt <- pull(ppc$interval_score_Rt)
    }

    # Static-parameter recovery (aggregate). The detailed per-parameter
    # view comes from `dgtf_compare_params()`. CRPS and bias are not
    # aggregated here because they mix scales across parameters
    # (W ~ 1e-3 vs rho ~ 30); coverage is a scale-free fraction and
    # reads well as an at-a-glance calibration signal.
    if (!is.null(truth)) {
        rec <- tryCatch(
            param_recovery(fit, truth, level = level),
            error = function(e) {
                warning(sprintf(
                    "param_recovery() failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
                })
        out$n_param_truth <- if (is.null(rec)) NA_integer_ else as.integer(nrow(rec))
        out$param_coverage_rate <- if (is.null(rec) || !nrow(rec)) NA_real_ else mean(rec$coverage)
    }
    out
}

#' @export
print.dgtf_comparison <- function(x, digits = 3L, ...) {
    cls <- attr(x, "input_class") %||% "dgtf_fit"
    cat(sprintf("<dgtf_comparison (%s): %d row%s, level = %.0f%%>\n",
                cls,
                nrow(x), if (nrow(x) == 1L) "" else "s",
                100 * (attr(x, "level") %||% 0.95)))

    pt <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
    int_cols <- c("n_param", "n_param_truth", "horizon", "n_draws")

    for (nm in setdiff(names(pt), "model")) {
        v <- pt[[nm]]
        if (!is.numeric(v)) next
        pt[[nm]] <- if (nm %in% int_cols) {
            ifelse(is.na(v), "", formatC(v, format = "d", big.mark = ","))
        } else if (grepl("^elapsed", nm)) {
            ifelse(is.na(v), "", sprintf("%.1fs", v))
        } else {
            formatC(v, digits = digits, format = "g")
        }
    }
    print(pt, row.names = FALSE)

    if (identical(cls, "dgtf_fit") && any(is.na(x$log_marglik)))
        cat("\nNote: `log_marglik` is the HVA SMC tail estimate of\n",
            "log p(y | gamma_final); MCMC fits have no equivalent stored\n",
            "and show NA. Pointwise log-likelihood is not currently\n",
            "exposed by either inference method.\n", sep = "")

    if (identical(cls, "dgtf_fit") && isTRUE(attr(x, "has_truth")))
        cat("\nUse `dgtf_compare_params(fits, truth)` for the\n",
            "per-parameter recovery table.\n", sep = "")

    invisible(x)
}

#' Per-parameter posterior summaries across multiple fits
#'
#' Long-form table stacking per-parameter posterior summaries for every
#' fit in `fits`, with an added `model` column. One row per
#' (parameter, fit) pair. Use this for the detailed cross-fit view of
#' static-parameter estimation.
#'
#' When `truth` is supplied, the table additionally reports recovery
#' metrics (`bias`, `mae`, `rmse`, `crps`, `coverage`) from
#' [`param_recovery()`]. When `truth` is `NULL`, only the posterior
#' summaries (`mean`, `median`, `sd`, `q_lo`, `q_hi`) are reported,
#' giving a side-by-side view of the marginal posteriors without
#' assuming a known data-generating process.
#'
#' For an at-a-glance summary across fits, use [`dgtf_compare()`]
#' (with `truth` supplied for the recovery columns).
#'
#' @param fits Named list of `dgtf_fit` objects.
#' @param truth Optional named numeric vector or named list of true
#'   static parameter values (same convention as
#'   [`param_recovery()`]). When `NULL` (default), only posterior
#'   summaries are reported.
#' @param level Credible-interval level for the `q_lo`, `q_hi`, and
#'   (when `truth` is supplied) `coverage` columns (default 0.95).
#'
#' @return A data frame of class `dgtf_param_comparison` with columns
#'   `parameter`, `model`, `mean`, `median`, `sd`, `q_lo`, `q_hi`,
#'   and, when `truth` is supplied, `truth`, `bias`, `mae`, `rmse`,
#'   `crps`, `coverage`. Rows are ordered by parameter, then by the
#'   input order of `fits`.
#'
#' @export
dgtf_compare_params <- function(fits, truth = NULL, level = 0.95) {
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits))))
        stop("`fits` must be a named, non-empty list of `dgtf_fit` ",
             "objects.", call. = FALSE)
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1)
        stop("`level` must be a single number in (0, 1).", call. = FALSE)

    rows <- lapply(seq_along(fits), function(i) {
        fit  <- fits[[i]]
        name <- names(fits)[[i]]
        if (!inherits(fit, "dgtf_fit"))
            stop(sprintf("Element `%s` is not a `dgtf_fit`.", name),
                 call. = FALSE)
        rec <- tryCatch(
            .compare_params_row(fit, truth, level),
            error = function(e) {
                warning(sprintf(
                    "compare_params row failed for `%s`: %s",
                    name, conditionMessage(e)), call. = FALSE)
                NULL
            })
        if (is.null(rec) || nrow(rec) == 0L) return(NULL)
        rec$model <- name
        rec
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)

    df <- do.call(rbind, rows)
    other <- setdiff(names(df), c("parameter", "model"))
    df <- df[, c("parameter", "model", other)]
    df$model <- factor(df$model, levels = names(fits))
    df <- df[order(df$parameter, df$model), , drop = FALSE]
    df$model <- as.character(df$model)
    rownames(df) <- NULL
    class(df) <- c("dgtf_param_comparison", "data.frame")
    attr(df, "level")     <- level
    attr(df, "has_truth") <- !is.null(truth)
    df
}

# Internal: produce one fit's worth of summary rows. Dispatches to
# `param_recovery()` when truth is provided, otherwise computes the
# posterior summaries directly from `.dgtf_static_draws()`.
.compare_params_row <- function(fit, truth, level) {
    if (!is.null(truth))
        return(param_recovery(fit, truth, level = level))

    d <- tryCatch(.dgtf_static_draws(fit), error = function(e) NULL)
    if (is.null(d) || is.null(d$draws_matrix) ||
        ncol(d$draws_matrix) == 0L)
        return(NULL)

    dm    <- d$draws_matrix
    alpha <- 1 - level
    qs    <- apply(dm, 2, stats::quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))

    data.frame(
        parameter = colnames(dm),
        mean      = round(apply(dm, 2, mean),         3),
        median    = round(apply(dm, 2, stats::median), 3),
        sd        = round(apply(dm, 2, stats::sd),     3),
        q_lo      = round(qs[1, ],                     3),
        q_hi      = round(qs[2, ],                     3),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

#' @export
print.dgtf_param_comparison <- function(x, digits = 3L, ...) {
    np  <- length(unique(x$parameter))
    nm  <- length(unique(x$model))
    lvl <- attr(x, "level") %||% 0.95
    has_truth <- isTRUE(attr(x, "has_truth"))

    cat(sprintf("<dgtf_param_comparison: %d parameter%s x %d model%s, level = %.0f%%>\n",
                np, if (np == 1L) "" else "s",
                nm, if (nm == 1L) "" else "s",
                100 * lvl))

    pt <- as.data.frame(unclass(x), stringsAsFactors = FALSE)

    # Convert the logical coverage column into a blank-named asterisk
    # column before the numeric formatting loop runs. Only present
    # when truth was supplied.
    if ("coverage" %in% names(pt)) {
        sig <- ifelse(is.na(pt$coverage), " ",
                      ifelse(pt$coverage, "*", " "))
        pt$coverage <- NULL
        pt[[" "]]   <- sig
    }

    for (col in setdiff(names(pt), c("parameter", "model", " "))) {
        pt[[col]] <- formatC(pt[[col]], digits = digits, format = "g")
    }

    print(pt, row.names = FALSE)
    if (has_truth)
        cat(sprintf("---\n'*' indicates the truth lies within the %g%% credible interval\n",
                    100 * lvl))
    invisible(x)
}
