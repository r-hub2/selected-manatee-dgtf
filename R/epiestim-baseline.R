# EpiEstim R_t baseline wrapper for use as a comparison input alongside
# `dgtf_fit` objects in `dgtf_compare_plot()`. EpiEstim is in Suggests,
# not Imports -- this file uses `requireNamespace()` and errors with
# an install hint when the package is missing.

#' EpiEstim R_t baseline for model comparison
#'
#' Thin wrapper around \code{EpiEstim::estimate_R()} and
#' \code{EpiEstim::wallinga_teunis()} that returns posterior R_t bands
#' in a format consumable by \code{dgtf_compare_plot()} -- a
#' `length(y) x 3` numeric matrix with columns `lower`, `median`,
#' `upper`, padded with leading NAs so row `i` corresponds to time
#' `i` in the input series (the same convention as
#' `dgtf_fit$fit$psi`).
#'
#' Intended use: drop a renewal-equation baseline next to `dgtf_fit`
#' objects in the comparison plot.
#'
#' ```
#' fit_hva  <- dgtf(y, mod, prior, method = "hva")
#' fit_mcmc <- dgtf(y, mod, prior, method = "mcmc")
#' epi      <- epiestim_baseline(y, mean_si = 4.7, sd_si = 2.9)
#'
#' dgtf_compare_plot(
#'     list(HVA = fit_hva, MCMC = fit_mcmc, EpiEstim = epi),
#'     what = "Rt")
#' ```
#'
#' @param y Numeric vector of incident counts (the same vector you
#'   pass to \code{dgtf()}). For daily data with weekly seasonality you
#'   typically want to deseasonalise `y` upstream before calling
#'   EpiEstim, since the renewal model assumes a smooth incidence
#'   signal.
#' @param mean_si,sd_si Serial-interval **mean** and **standard
#'   deviation** in the same time units as `y`. Note these are the
#'   mean / sd, not the variance.
#' @param method One of:
#'   * `"parametric_si"` (default) -- known SI mean / sd;
#'   * `"uncertain_si"` -- propagates uncertainty in SI mean / sd
#'     using `mean_sd_si`, `sd_sd_si` as prior standard deviations;
#'   * `"wallinga_teunis"` -- case-reproduction-number method.
#' @param window Sliding-window width in time units (default 7).
#'   Smaller windows give noisier, more responsive estimates.
#' @param level Credible-interval level (default 0.95). EpiEstim
#'   natively reports 0.025 / 0.975 quantiles; for any other level
#'   the bounds are derived from a Normal approximation to the
#'   posterior (mean ± z * Std(R)) and a warning is emitted.
#' @param align Alignment of the EpiEstim window summary in the returned
#'   matrix: `"center"` places each estimate near the middle of its window,
#'   while `"right"` places it at the window end.
#' @param mean_sd_si,sd_sd_si Prior standard deviations on the SI
#'   mean and SI sd respectively. Used only by
#'   `method = "uncertain_si"`.
#' @param ... Additional arguments forwarded to
#'   \code{EpiEstim::make_config()} (e.g. `n1`, `n2`, `n_sim`).
#'
#' @return An object of class `"epiestim_baseline"` -- a
#'   `length(y) x 3` numeric matrix with columns `lower`, `median`,
#'   `upper`. The first `window` rows are `NA` (no estimate yet).
#'   Method / level / SI parameters are stored as attributes.
#'
#' @export
epiestim_baseline <- function(y,
                              mean_si,
                              sd_si,
                              method = c(
                                  "parametric_si",
                                  "uncertain_si",
                                  "wallinga_teunis"
                              ),
                              window = 7L,
                              level = 0.95,
                              align = c("center", "right"),
                              mean_sd_si = 1.8,
                              sd_sd_si = 1.4,
                              ...) {
    method <- match.arg(method)
    align <- match.arg(align)

    if (!requireNamespace("EpiEstim", quietly = TRUE)) {
        stop("Package `EpiEstim` is required for `epiestim_baseline()`. ",
            "Install it with `install.packages(\"EpiEstim\")` ",
            "(it is listed in `Suggests`, not `Imports`).",
            call. = FALSE
        )
    }

    y <- as.numeric(y)
    n <- length(y)
    window <- as.integer(window)

    if (n < window + 2L) {
        stop(sprintf(
            "`y` must have length > window + 1 (got length = %d, window = %d).",
            n, window
        ), call. = FALSE)
    }
    if (!is.numeric(mean_si) || length(mean_si) != 1L || mean_si <= 0) {
        stop("`mean_si` must be a single positive number.", call. = FALSE)
    }
    if (!is.numeric(sd_si) || length(sd_si) != 1L || sd_si <= 0) {
        stop("`sd_si` must be a single positive number.", call. = FALSE)
    }
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1) {
        stop("`level` must be a single number in (0, 1).", call. = FALSE)
    }

    alpha <- 1 - level
    t_start <- 2L:(n - window + 1L)
    t_end <- t_start + (window - 1L)

    cfg_args <- c(list(t_start = t_start, t_end = t_end), list(...))
    if (method == "parametric_si") {
        cfg_args$mean_si <- mean_si
        cfg_args$std_si <- sd_si
    } else if (method == "uncertain_si") {
        cfg_args$mean_si <- mean_si
        cfg_args$std_mean_si <- mean_sd_si
        cfg_args$min_mean_si <- max(0, mean_si - 2 * mean_sd_si)
        cfg_args$max_mean_si <- mean_si + 2 * mean_sd_si
        cfg_args$std_si <- sd_si
        cfg_args$std_std_si <- sd_sd_si
        cfg_args$min_std_si <- max(0, sd_si - 2 * sd_sd_si)
        cfg_args$max_std_si <- sd_si + 2 * sd_sd_si
        cfg_args$n1 <- cfg_args$n1 %||% 100
        cfg_args$n2 <- cfg_args$n2 %||% 100
    } else { # wallinga_teunis
        cfg_args$mean_si <- mean_si
        cfg_args$std_si <- sd_si
        cfg_args$n_sim <- cfg_args$n_sim %||% 100
    }

    config <- do.call(EpiEstim::make_config, cfg_args)
    df_in <- data.frame(I = y)

    out <- if (method == "wallinga_teunis") {
        EpiEstim::wallinga_teunis(df_in,
            method = "parametric_si",
            config = config
        )$R
    } else {
        EpiEstim::estimate_R(df_in,
            method = method,
            config = config
        )$R
    }

    # Pull (lower, central, upper). EpiEstim returns Median(R) for the
    # parametric / uncertain branches and Mean(R) for Wallinga-Teunis,
    # so we accept either as the central column.
    col_lo <- sprintf("Quantile.%g(R)", alpha / 2)
    col_hi <- sprintf("Quantile.%g(R)", 1 - alpha / 2)
    col_md <- if ("Median(R)" %in% names(out)) "Median(R)" else "Mean(R)"

    bands <- if (col_lo %in% names(out) && col_hi %in% names(out)) {
        as.matrix(out[, c(col_lo, col_md, col_hi)])
    } else {
        # Non-default level path: EpiEstim only reports 0.025 / 0.975
        # quantiles, so fall back to a Gaussian approximation around
        # the posterior mean using the reported standard deviation.
        m <- out[[col_md]]
        s <- out[["Std(R)"]] %||%
            ((out[["Quantile.0.975(R)"]] - m) / stats::qnorm(0.975))
        z <- stats::qnorm(1 - alpha / 2)
        warning(sprintf(
            "EpiEstim did not report quantiles at level %.2f; using a Normal approximation (mean +/- %.2f * Std(R)).",
            level, z
        ), call. = FALSE)
        cbind(m - z * s, m, m + z * s)
    }
    colnames(bands) <- c("lower", "median", "upper")

    # Pad leading rows with NA so row i corresponds to time i in `y`.
    # EpiEstim's first estimate is at t_end[1] = window + 1, so we
    # need `window` leading NA rows (covering times 1..window).
    nr <- nrow(bands)
    if (align == "right") {
        pad_lead <- n - nr
        pad_trail <- 0L
    } else {
        # Center: shift left by floor((window - 1) / 2) relative to right-edge.
        half <- as.integer(floor((window - 1L) / 2))
        pad_lead <- n - nr - half
        pad_trail <- half
    }
    if (pad_lead > 0L) {
        bands <- rbind(
            matrix(NA_real_, pad_lead, 3,
                dimnames = list(NULL, colnames(bands))
            ),
            bands
        )
    }
    if (pad_trail > 0L) {
        bands <- rbind(
            bands,
            matrix(NA_real_, pad_trail, 3,
                dimnames = list(NULL, colnames(bands))
            )
        )
    }

    structure(bands,
        class   = c("epiestim_baseline", "matrix", "array"),
        method  = method,
        window  = window,
        align   = align,
        level   = level,
        mean_si = mean_si,
        sd_si   = sd_si
    )
}


#' @export
print.epiestim_baseline <- function(x, ...) {
    cat(sprintf(
        "<epiestim_baseline: method = %s, window = %d, align = %s, level = %.2f>\n",
        attr(x, "method"), attr(x, "window"), attr(x, "align"), attr(x, "level")
    ))
    cat(sprintf(
        "  SI mean = %g, SI sd = %g\n",
        attr(x, "mean_si"), attr(x, "sd_si")
    ))
    cat(sprintf(
        "  %d time points (%d NA, %d estimated)\n",
        nrow(x),
        sum(is.na(x[, "median"])),
        sum(!is.na(x[, "median"]))
    ))
    cat("\nLast 6 rows:\n")
    print(utils::tail(unclass(x), n = 6L))
    invisible(x)
}
