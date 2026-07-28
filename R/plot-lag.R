#' Lag-distribution comparison plot
#'
#' Compares per-lag weights across one or more `dgtf_fit` objects on a
#' single grouped bar chart. Useful for characterizing how different
#' DGTF model variants (Hawkes / distributed lag / Koyck / AR)
#' distribute past-observation influence across lags.
#'
#' Three lag-weight conventions are unified onto a common visual:
#'   * **Hawkes** (`sys = shift`, `lag = lognormal`): the discretized
#'     lognormal PMF at posterior-median (`meanlog`, `sigma2`), via
#'     \code{dlognorm0()}.
#'   * **Distributed lag / Koyck** (`sys = nbinom` or `shift`,
#'     `lag = nbinom`): the discretized NB PMF at posterior-median
#'     (`kappa`, `r`), via \code{dnbinom0()}. Both system forms share the
#'     same NB lag weights; only the underlying state-space
#'     representation differs.
#'   * **AR** (`sys = identity`, `lag = uniform`): the gain-applied
#'     AR coefficient `h(theta_l)` summarized across time and posterior
#'     samples for each lag `l`. These are *not* probabilities — see
#'     `normalize` for how they're placed on a comparable scale.
#'
#' @param fits Either a `dgtf_fit` (treated as a singleton list with
#'   name `"Estimate"`), or a named list of `dgtf_fit` objects.
#' @param window Maximum lag to display. `NULL` (default) uses the
#'   largest lag window across all fits.
#' @param xlim Display-only truncation: maximum lag to show in the
#'   plot. Defaults to `window` (no truncation). Useful when the
#'   underlying lag window is wide (e.g. 30 days) but you only want
#'   to show the first 10 days. Weights are still computed and
#'   normalized over the full `window`, so the visible bars retain
#'   their correct relative mass; the tail past `xlim` is implicit.
#' @param normalize One of:
#'   * `"pmf"` (default): rescale each model's weights to sum to 1.
#'     Compares lag *shapes*; hides the AR / branching-ratio scale.
#'   * `"raw"`: show as-is. Hawkes and DL/Koyck stay near unit mass
#'     (truncation residual aside); AR shows the raw branching
#'     weights `h(theta_l)`. Most honest when magnitude matters.
#'   * `"branching"`: rescale each PMF to integrate to the model's
#'     mean `R_t` (estimated as the time-mean of `h(psi)`); for AR
#'     the natural sum is `sum(h(theta_l))`. Same shape information
#'     as `"pmf"` with absolute magnitudes restored.
#' @param summary Posterior summary applied to each model's
#'   parameter / coefficient draws before computing weights.
#'   `"median"` (default) or `"mean"`.
#' @param alpha Bar transparency (default 0.5).
#' @param palette Optional named character vector of fill colours.
#'   Same convention as \code{dgtf_compare_plot()}.
#' @param xlab,ylab,main Axis / title labels. `ylab = NULL` (default)
#'   picks a label that matches `normalize`.
#' @param theme A ggplot2 theme.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#'
#' @export
#' @examples
#' \dontrun{
#' fit_hawkes <- dgtf(y, dgtf_hawkes(), prior, method = "hva")
#' fit_dl <- dgtf(y, dgtf_distributed_lag(r = 3), prior, method = "hva")
#' fit_ar <- dgtf(y, dgtf_nb_ar(p = 5), prior, method = "hva")
#'
#' plot_dgtf_lag(list(
#'     Hawkes = fit_hawkes,
#'     DL3 = fit_dl,
#'     "AR(5)" = fit_ar
#' ))
#'
#' # Honest magnitudes (AR's raw branching weights look bigger than
#' # the unit-mass shift-system PMFs):
#' plot_dgtf_lag(list(Hawkes = fit_hawkes, "AR(5)" = fit_ar),
#'     normalize = "raw"
#' )
#' }
plot_dgtf_lag <- function(fits,
                          window = NULL,
                          xlim = NULL,
                          normalize = c("pmf", "raw", "branching"),
                          summary = c("median", "mean"),
                          alpha = 0.5,
                          palette = NULL,
                          xlab = "Lag",
                          ylab = NULL,
                          main = NULL,
                          theme = ggplot2::theme_minimal(),
                          ...) {
    normalize <- match.arg(normalize)
    summary <- match.arg(summary)
    summarize <- if (summary == "median") stats::median else mean

    if (inherits(fits, "dgtf_fit")) fits <- list(Estimate = fits)
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits)))) {
        stop("`fits` must be a named, non-empty list of `dgtf_fit` ",
            "objects.",
            call. = FALSE
        )
    }
    if (!all(vapply(fits, inherits, logical(1), "dgtf_fit"))) {
        stop("All elements of `fits` must be `dgtf_fit` objects.",
            call. = FALSE
        )
    }

    # Determine display window. Prefer the model's lag window; for AR
    # this is the AR order p (also stored as lag$window).
    nLs <- vapply(fits, function(f) {
        as.integer(f$model$lag$window %||% 1L)
    }, integer(1))
    if (is.null(window)) window <- max(nLs)
    window <- as.integer(window)

    if (is.null(xlim)) {
        xlim_eff <- window
    } else {
        if (!is.numeric(xlim) || length(xlim) != 1L || xlim < 1L) {
            stop("`xlim` must be a single positive integer.", call. = FALSE)
        }
        xlim_eff <- min(as.integer(xlim), window)
    }

    rows <- lapply(seq_along(fits), function(i) {
        fit <- fits[[i]]
        name <- names(fits)[[i]]
        w <- .dgtf_lag_weights(fit,
            window = window,
            summarize = summarize
        )
        if (normalize == "pmf" && sum(w, na.rm = TRUE) > 0) {
            w <- w / sum(w, na.rm = TRUE)
        }
        if (normalize == "branching") {
            scale_factor <- .dgtf_branching_scale(fit) %||%
                sum(w, na.rm = TRUE)
            if (sum(w, na.rm = TRUE) > 0) {
                w <- (w / sum(w, na.rm = TRUE)) * scale_factor
            }
        }
        data.frame(
            lag = seq_len(window),
            weight = w,
            name = name,
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    if (xlim_eff < window) {
        df <- df[df$lag <= xlim_eff, , drop = FALSE]
    }

    method_names <- unique(df$name)
    cols <- .dgtf_palette(method_names, palette)

    ggplot2::ggplot(df, ggplot2::aes(x = lag, y = weight, fill = name)) +
        ggplot2::geom_col(position = "dodge", alpha = alpha) +
        ggplot2::scale_fill_manual("", values = cols) +
        ggplot2::scale_x_continuous(breaks = seq_len(xlim_eff)) +
        theme +
        ggplot2::labs(
            x = xlab,
            y = ylab %||% switch(normalize,
                pmf       = "Lag PMF",
                raw       = "Lag weight",
                branching = "Branching contribution"
            ),
            title = main
        ) +
        ggplot2::theme(legend.position = "top")
}


# ----------------------------------------------------------------
# Internal helpers (kept in this file to keep the lag plot
# self-contained).
# ----------------------------------------------------------------

# Length-`window` lag weights for one fit. Dispatches on
# (sys, lag) type:
#   shift + lognormal  -> dlognorm0(window, par1, par2)
#   shift + nbinom     -> dnbinom0(window, par1=kappa, par2=r)
#   identity (any lag) -> h(theta_l) summarized across time and draws,
#                         padded / truncated to length `window`.
.dgtf_lag_weights <- function(fit, window, summarize) {
    sys_type <- fit$model$sys$type
    lag_type <- fit$model$lag$type
    L <- as.integer(window)

    pull_param <- function(slot, model_default) {
        v <- fit$fit[[slot]]
        if (!is.null(v) && length(v)) {
            return(summarize(as.numeric(v)))
        }
        as.numeric(model_default)
    }

    if (sys_type == "shift" && lag_type == "lognorm") {
        # par1 = meanlog, par2 = sigma2
        ml <- pull_param("par1", fit$model$lag$params$par1)
        s2 <- pull_param("par2", fit$model$lag$params$par2)
        if (!is.finite(ml) || !is.finite(s2)) {
            stop("Cannot extract lognormal lag parameters.", call. = FALSE)
        }
        return(as.numeric(dlognorm0(L, ml, s2)))
    }
    if ((sys_type == "shift" || sys_type == "nbinom") &&
        (lag_type == "nbinom" || lag_type == "nbinomp")) {
        # par1 = kappa, par2 = r  (see R/components.R:175)
        # Same NB lag PMF whether the engine runs the sliding-truncated
        # or the iterative-exact form; both share the discretized weights.
        kappa <- pull_param("par1", fit$model$lag$params$par1)
        r <- pull_param("par2", fit$model$lag$params$par2)
        if (!is.finite(kappa) || !is.finite(r)) {
            stop("Cannot extract NB-lag parameters.", call. = FALSE)
        }
        return(as.numeric(dnbinom0(L, kappa, r)))
    }
    if (sys_type == "identity") {
        # AR: lag weights = h(theta_l) summarized over time x draws.
        theta <- fit$fit$Theta_stored %||% fit$fit$Theta
        if (is.null(theta)) {
            stop("Cannot extract AR coefficient cube `Theta` from this fit. ",
                "Identity-system lag weights require it (HVA exposes ",
                "this; MCMC AR fits currently do not).",
                call. = FALSE
            )
        }
        if (length(dim(theta)) != 3L) {
            stop("`Theta` must be a 3-D array (nP x ntime x nsample).",
                call. = FALSE
            )
        }
        nP <- dim(theta)[1L]
        gain_name <- fit$model$gain$type %||% "softplus"
        h <- .gain_fun(gain_name)
        weights <- vapply(seq_len(nP), function(l) {
            slice_l <- theta[l, , ] # (ntime+1) x nsample
            summarize(as.numeric(h(as.numeric(slice_l))))
        }, numeric(1))
        if (length(weights) < L) {
            weights <- c(weights, rep(0, L - length(weights)))
        }
        return(weights[seq_len(L)])
    }

    stop(sprintf(
        "Cannot compute lag weights for sys = '%s', lag = '%s'.",
        sys_type, lag_type
    ), call. = FALSE)
}


# Branching-ratio scale: time-mean of the posterior-median `R_t` for
# shift-system fits (psi posterior is in fit$fit$psi as an n x 3
# quantile matrix); for identity-system fits the natural scale comes
# from sum(h(theta_l)) and is computed directly by the caller via
# `sum(w)` (so this function returns NULL for that case and lets the
# fall-through handle it).
.dgtf_branching_scale <- function(fit) {
    sys_type <- fit$model$sys$type
    # Both sliding (sys_shift) and iterative (sys_nbinom) carry psi as
    # the first row of Theta and use R_t = h(psi_t); only the identity
    # (AR) form needs its own scale path.
    if (sys_type == "identity") {
        return(NULL)
    }

    psi <- fit$fit$psi
    if (is.null(psi) || !is.matrix(psi) || ncol(psi) < 2L) {
        return(NULL)
    }

    gain_name <- fit$model$gain$type %||% "softplus"
    h <- .gain_fun(gain_name)
    median_path <- as.numeric(psi[, 2]) # posterior-median psi(t)
    mean(h(median_path), na.rm = TRUE)
}
