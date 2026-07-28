.okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
)
.truth_pattern <- "^(true|truth)$"
.default_ylab <- function(what) {
    switch(what,
        psi  = expression(psi[t]),
        Rt   = expression(R[t]),
        yhat = expression(hat(y)[t])
    )
}


.ts_targets <- c("Rt", "psi", "yhat")
.diag_targets <- c("marglik")

# Soft replacement for match.arg(): time-series targets are an enum, but
# diagnostic targets (e.g. "marglik") have their own dispatch path, and
# any other single string is treated as a static-parameter name (validated
# downstream against `.dgtf_static_draws(fit)$inferred`).
.what_kind <- function(what) {
    if (length(what) != 1L || !is.character(what)) {
        stop("`what` must be a single string.", call. = FALSE)
    }
    if (what %in% .ts_targets) {
        return("ts")
    }
    if (what %in% .diag_targets) {
        return("diag")
    }
    "param"
}


# Mirror of src/GainFunc.h:43 (psi2hpsi). Logistic is registered C++-side
# but falls through to identity in the switch, so we omit it here.
.gain_fun <- function(name) {
    switch(name,
        identity = identity,
        softplus = function(z) log1p(exp(z)),
        exponential = exp,
        ramp = function(z) pmax(z, .Machine$double.eps),
        stop(sprintf("Unknown gain '%s'; extend .gain_fun().", name),
            call. = FALSE
        )
    )
}

# Row-wise quantiles of an n x nsample draws matrix at the central
# (1-level)/2, 0.5, 1-(1-level)/2 probabilities. Returns an n x 3 matrix.
.row_quantiles <- function(draws, level = 0.95) {
    probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
    out <- t(apply(draws, 1L, stats::quantile,
        probs = probs, na.rm = TRUE, names = FALSE
    ))
    colnames(out) <- c("lower", "median", "upper")
    out
}

# Lazy PPC fallback (only triggered when Rt/yhat is requested but no PPC is
# attached and -- for Rt -- psi draws are unavailable).
.require_ppc <- function(fit, level) {
    warning("No PPC attached to fit; running posterior_predict() with ",
        "default settings (nrep = 100, no Rt reference). Attach one ",
        "with `attr(fit, 'ppc') <- posterior_predict(fit, ...)` to ",
        "avoid this on every plot call.",
        call. = FALSE
    )
    posterior_predict(fit, level = level)
}

# Normalize any input to a long-form band data.frame with columns:
#   time (numeric), lower, central, upper, name (character), has_ribbon (logical)
# Dispatches on class(x):
#   dgtf_fit  + what == "psi"             -> x$fit$psi          (T x 3)
#   dgtf_fit  + what %in% c("Rt","yhat")  -> attr(x, "ppc"), or
#                                            posterior_predict(x, level = level)
#                                            on the fly (with a warning)
#   dgtf_ppc  + what %in% c("Rt","yhat")  -> x$Rt or x$yhat     (T+1 x 3)
#   matrix (T x 3)                        -> use as-is (escape hatch for
#                                            external methods like EpiEstim)
#   numeric vector                        -> central only,
#                                            has_ribbon = FALSE  (truth)
.dgtf_band <- function(x, what, level = 0.95, time = NULL, name = NULL) {
    # Returns data.frame:
    #   time       numeric
    #   lower      numeric (NA when has_ribbon = FALSE)
    #   central    numeric (median)
    #   upper      numeric (NA when has_ribbon = FALSE)
    #   name       character (recycled)
    #   has_ribbon logical (recycled)

    # Branch 1: dgtf_fit
    #   psi  -> use the precomputed quantile matrix `fit$fit$psi`.
    #   Rt   -> apply the model's gain function elementwise to the
    #           `n x nsample` posterior draws of psi (via the public
    #           `extract_state_draws()` from R/draws.R:65) and take row-wise
    #           quantiles. No C++ posterior_predict round-trip needed.
    #   yhat -> requires posterior_predict (observation noise sampling), so
    #           use any attached PPC, otherwise warn + run on the fly.
    if (inherits(x, "dgtf_fit")) {
        if (what == "psi") {
            qm <- x$fit$psi # n x 3
        } else if (what == "Rt") {
            psi_draws <- tryCatch(extract_state_draws(x, "psi"),
                error = function(e) NULL
            )
            if (!is.null(psi_draws)) {
                gain_name <- x$model$gain$type %||% "identity"
                h <- .gain_fun(gain_name)
                Rt_draws <- h(psi_draws) # n x nsample
                qm <- .row_quantiles(Rt_draws, level = level)
            } else {
                # Fallback: fits without stored psi draws fall back to PPC.
                ppc <- attr(x, "ppc") %||% .require_ppc(x, level = level)
                qm <- ppc$Rt
                if (is.null(qm)) {
                    stop("PPC has no `Rt` field. Re-run posterior_predict() ",
                        "with a non-NULL `Rt` argument or attach a PPC ",
                        "that contains it.",
                        call. = FALSE
                    )
                }
            }
        } else if (what == "yhat") {
            ppc <- attr(x, "ppc") %||% .require_ppc(x, level = level)
            qm <- ppc$yhat
            if (is.null(qm)) {
                stop("PPC has no `yhat` field. Re-run posterior_predict() ",
                    "with `nrep > 0`.",
                    call. = FALSE
                )
            }
        } else {
            stop(sprintf("Unknown `what = %s` for dgtf_fit.", what),
                call. = FALSE
            )
        }
        return(.dgtf_band_from_matrix(qm,
            time = time, name = name,
            has_ribbon = TRUE
        ))
    }

    # Branch 2: dgtf_ppc
    if (inherits(x, "dgtf_ppc")) {
        if (!what %in% c("Rt", "yhat")) {
            stop(sprintf("Unknown `what = %s` for dgtf_ppc.", what),
                call. = FALSE
            )
        }
        qm <- x[[what]]
        if (is.null(qm)) {
            stop(sprintf("PPC object has no `%s` field.", what),
                call. = FALSE
            )
        }
        return(.dgtf_band_from_matrix(qm,
            time = time, name = name,
            has_ribbon = TRUE
        ))
    }

    # Branch 2b: dgtf_forecast -- only yhat
    if (inherits(x, "dgtf_forecast")) {
        if (!what %in% c("yhat", "y"))
            stop(sprintf(
                "Unknown `what = %s` for dgtf_forecast; only 'yhat' supported.",
                what), call. = FALSE)
        
        qm <- x$quantiles
        if (is.null(qm) || nrow(qm) == 0L)
            stop("Forecast object has no `quantiles`.", call. = FALSE)
        return(.dgtf_band_from_matrix(
            as.matrix(qm), time = time, name = name, has_ribbon = TRUE
        ))
    }

    # Branch 3: T x 3 quantile matrix or data frame (lower/median/upper).
    if ((is.matrix(x) || is.data.frame(x)) && ncol(x) == 3L) {
        return(.dgtf_band_from_matrix(as.matrix(x),
            time = time, name = name,
            has_ribbon = TRUE
        ))
    }

    # Branch 4: numeric vector (or array with at most one non-trivial
    # dimension, e.g. simulate()'s n x 1 psi matrix) -> truth, no ribbon.
    if (is.numeric(x) &&
        (is.null(dim(x)) || sum(dim(x) != 1L) <= 1L)) {
        x_vec <- as.numeric(x)
        n <- length(x_vec)
        if (is.null(time)) {
            time <- seq(0, n - 1)
        } else {
            stopifnot(length(time) == n)
        }
        return(data.frame(
            time = as.numeric(time),
            lower = NA_real_,
            central = x_vec,
            upper = NA_real_,
            name = name %||% "Truth",
            has_ribbon = FALSE,
            stringsAsFactors = FALSE
        ))
    }

    stop("Unsupported input to .dgtf_band(): ",
        paste(class(x), collapse = "/"),
        call. = FALSE
    )
}

# Helper: turn a T x 3 quantile matrix into a long band data.frame.
.dgtf_band_from_matrix <- function(qm, time, name, has_ribbon) {
    stopifnot(is.matrix(qm), ncol(qm) == 3L)
    n <- nrow(qm)
    if (is.null(time)) {
        time <- seq(0, n - 1) # matches plot_ts() convention in
        # R/diagnostics.R:109 (time = 0:(T-1))
    } else {
        # Allow callers to pass `time = band$time` for re-indexing truth onto
        # an existing band's axis. Recycle / truncate to length n if needed.
        if (length(time) != n) {
            if (length(time) > n) {
                time <- time[seq_len(n)]
            } else {
                time <- c(time, seq(max(time) + 1L,
                    length.out = n - length(time)
                ))
            }
        }
    }
    data.frame(
        time = as.numeric(time),
        lower = as.numeric(qm[, 1]),
        central = as.numeric(qm[, 2]),
        upper = as.numeric(qm[, 3]),
        name = name %||% "Estimate",
        has_ribbon = has_ribbon,
        stringsAsFactors = FALSE
    )
}

# Deterministic palette assignment.
# names: character vector of method names (already filtered to drop entries
#        matching .truth_pattern -- truth is rendered separately).
# palette: NULL | character vector | named character vector. Named entries
#        pin specific names to specific hexes; unnamed entries fill from
#        the default generator in alphabetical order of remaining names.
# Returns a named character vector keyed by sort(unique(names)).
.dgtf_palette <- function(names, palette = NULL) {
    names <- unique(names)
    if (!length(names)) {
        return(stats::setNames(character(0), character(0)))
    }
    n <- length(names)
    auto <- if (n == 1L) {
        "black"
    } else if (n <= 8L) {
        .okabe_ito[seq_len(n)]
    } else {
        grDevices::hcl.colors(n, palette = "Set 2")
    }
    # Default assignment: alphabetical order over `names` -> Okabe-Ito sequence.
    out <- stats::setNames(auto, sort(names))[names]

    if (is.null(palette)) {
        return(out)
    }

    if (is.null(names(palette))) {
        # Unnamed vector: pin by position over `names`; if shorter than n,
        # the remaining entries keep the auto colors.
        m <- min(length(palette), n)
        out[seq_len(m)] <- palette[seq_len(m)]
        return(out)
    }

    # Named vector: pin the named entries; auto fills the rest. Unrecognized
    # names in `palette` are silently ignored.
    pinned <- intersect(names(palette), names)
    if (length(pinned)) out[pinned] <- palette[pinned]
    out
}

# Shared ggplot renderer. `bands` is the long data.frame from .dgtf_band
# (possibly rbind'd from many sources).
.plot_dgtf_band <- function(bands, palette = NULL, alpha = 0.5,
                            xlab = "Time", ylab = NULL, main = NULL,
                            theme = ggplot2::theme_minimal(),
                            legend.position = "right",
                            legend.name = "Method",
                            legend.nrow = NULL,
                            hline = NULL) {
    is_truth <- grepl(.truth_pattern, bands$name, ignore.case = TRUE)
    truth <- bands[is_truth, , drop = FALSE]
    methods <- bands[!is_truth, , drop = FALSE]
    method_names <- unique(methods$name)
    cols <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(
        x = time, y = central,
        group = name
    )) +
        (if (length(hline) && all(is.finite(hline))) {
            ggplot2::geom_hline(
                yintercept = as.numeric(hline),
                linetype = "dashed",
                color = "gray50",
                linewidth = 0.4,
                alpha = 0.7
            )
        } else {
            NULL
        }) +
        ggplot2::geom_ribbon(
            data = methods[methods$has_ribbon, , drop = FALSE],
            ggplot2::aes(ymin = lower, ymax = upper, fill = name),
            alpha = alpha, na.rm = TRUE
        ) +
        ggplot2::geom_line(ggplot2::aes(color = name), na.rm = TRUE) +
        ggplot2::scale_color_manual(name = legend.name, values = cols) +
        ggplot2::scale_fill_manual(name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (!is.null(legend.nrow)) {
        p <- p + ggplot2::guides(
            color = ggplot2::guide_legend(nrow = legend.nrow),
            fill  = ggplot2::guide_legend(nrow = legend.nrow)
        )
    }

    if (nrow(truth)) {
        truth_names <- unique(truth$name)
        p <- p + ggplot2::geom_line(
            data = truth,
            ggplot2::aes(
                x = time, y = central,
                linetype = name, group = name
            ),
            color = "black", inherit.aes = FALSE
        ) +
            ggplot2::scale_linetype_manual(
                name = NULL,
                values = stats::setNames(
                    rep("dashed", length(truth_names)),
                    truth_names
                )
            )
    }
    p
}

# Returns a data.frame with columns:
#   name   character  -- method/source label (recycled)
#   index  integer    -- 1..p_dim; identifies which entry of a vector param
#   value  numeric    -- a single draw (or the truth value)
#   kind   character  -- "draws" or "truth"
#
# Dispatches on class(x):
#   dgtf_fit   -> pulls .dgtf_static_draws(x)$draws[[param]] (nsample x p_dim)
#                 and reshapes to long.
#   numeric    -> length 1 OR length p_dim; treated as truth.
#   matrix     -> nsample x p_dim, treated as draws (escape hatch for
#                 externally-fitted methods).
#
# `p_dim_hint` is consulted only for the numeric truth branch (we need to
# know how many indices to broadcast a length-1 truth across, OR validate
# that a length-k truth matches the expected p_dim). It is filled in by the
# compare wrapper after at least one fit has been resolved.
.dgtf_param_long <- function(x, param, name = NULL, p_dim_hint = NULL) {
    if (inherits(x, "dgtf_fit")) {
        sd <- .dgtf_static_draws(x)
        if (!param %in% sd$inferred) {
            stop(sprintf(
                "Parameter '%s' is not inferred for this fit. Inferred: %s",
                param, paste(sd$inferred, collapse = ", ")
            ), call. = FALSE)
        }
        m <- sd$draws[[param]] # nsample x p_dim
        p <- ncol(m)
        return(data.frame(
            name = name %||% "Estimate",
            index = rep(seq_len(p), each = nrow(m)),
            iter = rep(seq_len(nrow(m)), times = p),
            value = as.numeric(m),
            kind = "draws",
            stringsAsFactors = FALSE
        ))
    }
    if (is.matrix(x)) {
        return(data.frame(
            name = name %||% "Estimate",
            index = rep(seq_len(ncol(x)), each = nrow(x)),
            iter = rep(seq_len(nrow(x)), times = ncol(x)),
            value = as.numeric(x),
            kind = "draws",
            stringsAsFactors = FALSE
        ))
    }
    if (is.numeric(x) && is.null(dim(x))) {
        n <- length(x)
        if (!is.null(p_dim_hint) && n != 1L && n != p_dim_hint) {
            stop(sprintf(
                "Truth for '%s' has length %d, expected 1 or %d.",
                param, n, p_dim_hint
            ), call. = FALSE)
        }
        idx <- if (n == 1L && !is.null(p_dim_hint)) {
            seq_len(p_dim_hint)
        } else {
            seq_len(n)
        }
        val <- if (n == 1L && !is.null(p_dim_hint)) rep(x, p_dim_hint) else x
        return(data.frame(
            name = name %||% "Truth",
            index = idx,
            iter = NA_integer_,
            value = as.numeric(val),
            kind = "truth",
            stringsAsFactors = FALSE
        ))
    }
    stop("Unsupported input to .dgtf_param_long(): ",
        paste(class(x), collapse = "/"),
        call. = FALSE
    )
}

.plot_dgtf_param_hist <- function(long, palette = NULL, alpha = 0.5,
                                  bins = 50, xlab = NULL, ylab = "Count",
                                  main = NULL,
                                  theme = ggplot2::theme_minimal(),
                                  legend.position = "right",
                                  legend.name = "Method") {
    is_truth <- long$kind == "truth"
    truth <- long[is_truth, , drop = FALSE]
    methods <- long[!is_truth, , drop = FALSE]
    method_names <- unique(methods$name)
    cols <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(x = value, fill = name)) +
        ggplot2::geom_histogram(
            position = "identity", alpha = alpha,
            bins = bins, na.rm = TRUE
        ) +
        ggplot2::scale_fill_manual(name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (nrow(truth)) {
        p <- p + ggplot2::geom_vline(
            xintercept = truth$value,
            color = "black", linetype = "dashed"
        )
    }
    p
}

.plot_dgtf_param_trace <- function(long, palette = NULL, alpha = 0.5,
                                   xlab = "Iteration", ylab = NULL,
                                   main = NULL,
                                   theme = ggplot2::theme_minimal(),
                                   legend.position = "right",
                                   legend.name = "Method") {
    is_truth <- long$kind == "truth"
    truth <- long[is_truth, , drop = FALSE]
    draws <- long[!is_truth, , drop = FALSE]
    draws$index_f <- factor(draws$index)

    multi_index <- length(unique(draws$index)) > 1L

    if (multi_index) {
        # Per-index discrete palette; method (`name`) collapsed into group.
        p <- ggplot2::ggplot(
            draws,
            ggplot2::aes(
                x = iter, y = value,
                color = index_f,
                group = interaction(name, index_f)
            )
        ) +
            ggplot2::geom_line(alpha = alpha, na.rm = TRUE) +
            ggplot2::scale_color_discrete(name = "Index") +
            theme +
            ggplot2::labs(x = xlab, y = ylab, title = main) +
            ggplot2::theme(legend.position = legend.position)
    } else {
        # Single-index: color by `name` through `.dgtf_palette` so the
        # default for a lone "Estimate" is black + alpha = 0.5 -> the
        # same gray as the histogram. User-passed `color =` flows in
        # via `palette = c(Estimate = color)` from the caller.
        method_names <- unique(draws$name)
        cols <- .dgtf_palette(method_names, palette)
        p <- ggplot2::ggplot(
            draws,
            ggplot2::aes(
                x = iter, y = value,
                color = name, group = name
            )
        ) +
            ggplot2::geom_line(alpha = alpha, na.rm = TRUE) +
            ggplot2::scale_color_manual(name = legend.name, values = cols) +
            theme +
            ggplot2::labs(x = xlab, y = ylab, title = main) +
            ggplot2::theme(legend.position = legend.position) +
            ggplot2::guides(color = "none") # no legend for a single line
    }

    if (nrow(truth)) {
        p <- p + ggplot2::geom_hline(
            yintercept = truth$value,
            color = "black", linetype = "dashed"
        )
    }

    p
}

.plot_dgtf_param_box <- function(long, palette = NULL, alpha = 0.5,
                                 xlab = "Index", ylab = "Value",
                                 main = NULL,
                                 theme = ggplot2::theme_minimal(),
                                 legend.position = "right",
                                 legend.name = "Method") {
    is_truth <- long$kind == "truth"
    truth <- long[is_truth, , drop = FALSE]
    methods <- long[!is_truth, , drop = FALSE]
    methods$index_f <- factor(methods$index)
    method_names <- unique(methods$name)
    cols <- .dgtf_palette(method_names, palette)

    p <- ggplot2::ggplot(methods, ggplot2::aes(
        x = index_f, y = value,
        fill = name
    )) +
        ggplot2::geom_boxplot(alpha = alpha, outlier.size = 0.6) +
        ggplot2::scale_fill_manual(name = legend.name, values = cols) +
        theme +
        ggplot2::labs(x = xlab, y = ylab, title = main) +
        ggplot2::theme(legend.position = legend.position)

    if (nrow(truth)) {
        truth$index_f <- factor(truth$index,
            levels = levels(methods$index_f)
        )
        # Drop method/group from inheritance so the truth marker is
        # rendered once per index regardless of how many methods are
        # plotted. Use shape = 23 (filled diamond) to read clearly on
        # both light and dark themes.
        p <- p + ggplot2::geom_point(
            data = truth,
            ggplot2::aes(x = index_f, y = value),
            inherit.aes = FALSE,
            color = "black", fill = "black",
            shape = 23, size = 3
        )
    }
    p
}

.plot_dgtf_param_single <- function(x, param, truth = NULL,
                                    type = "hist",
                                    color = NULL, alpha = 0.5,
                                    bins = 50, ...) {
    type <- match.arg(type, c("hist", "trace"))
    sd <- .dgtf_static_draws(x)
    if (!param %in% sd$inferred) {
        stop(sprintf(
            "Parameter '%s' is not inferred for this fit. Inferred: %s",
            param, paste(sd$inferred, collapse = ", ")
        ), call. = FALSE)
    }

    if (type == "trace" && !identical(x$method, "mcmc")) {
        stop("Trace plots require an MCMC fit; this fit is '",
            x$method, "'. HVA produces i.i.d. importance samples, ",
            "not chains.",
            call. = FALSE
        )
    }

    p_dim <- ncol(sd$draws[[param]])
    long <- .dgtf_param_long(x, param, name = "Estimate")
    if (!is.null(truth)) {
        long <- rbind(
            long,
            .dgtf_param_long(truth, param,
                name = "Truth",
                p_dim_hint = p_dim
            )
        )
    }

    if (type == "trace") {
        return(.plot_dgtf_param_trace(long,
            alpha = alpha,
            ylab = param, ...
        ))
    }

    palette <- if (!is.null(color)) c(Estimate = color) else NULL
    legend <- if (is.null(truth)) "none" else "right"
    if (p_dim == 1L) {
        .plot_dgtf_param_hist(long,
            palette = palette, alpha = alpha,
            bins = bins, xlab = param,
            legend.position = legend, ...
        )
    } else {
        .plot_dgtf_param_box(long,
            palette = palette, alpha = alpha,
            ylab = param, legend.position = legend, ...
        )
    }
}

.plot_dgtf_param_compare <- function(fits, param,
                                     palette = NULL, alpha = 0.5,
                                     bins = 50,
                                     legend.position = "top", ...) {
    stopifnot(
        is.list(fits), length(fits) >= 1L,
        !is.null(names(fits)), all(nzchar(names(fits)))
    )

    # Resolve p_dim from the first dgtf_fit found; needed to broadcast /
    # validate any numeric-truth entries in the list.
    p_dim_hint <- NULL
    for (e in fits) {
        if (inherits(e, "dgtf_fit")) {
            sd <- .dgtf_static_draws(e)
            if (param %in% sd$inferred) {
                p_dim_hint <- ncol(sd$draws[[param]])
                break
            }
        } else if (is.matrix(e)) {
            p_dim_hint <- ncol(e)
            break
        }
    }
    if (is.null(p_dim_hint)) {
        stop(sprintf(
            "No fit in the input has '%s' inferred (and no draws matrix ",
            "supplied to seed p_dim).", param
        ), call. = FALSE)
    }

    long <- do.call(rbind, Map(
        function(x, nm) {
            .dgtf_param_long(x, param,
                name = nm,
                p_dim_hint = p_dim_hint
            )
        },
        fits, names(fits)
    ))

    if (p_dim_hint == 1L) {
        .plot_dgtf_param_hist(long,
            palette = palette, alpha = alpha,
            bins = bins, xlab = param,
            legend.position = legend.position, ...
        )
    } else {
        .plot_dgtf_param_box(long,
            palette = palette, alpha = alpha,
            ylab = param,
            legend.position = legend.position, ...
        )
    }
}

# Build a marglik convergence plot for a single dgtf_fit.
#
# The HVA marglik trace is heavy-tailed (occasional large outlier
# iterations from degenerate SMC particle weights at moving variational
# parameters). All visual summaries use robust estimators so the
# outliers don't drown the convergence signal.
#
# Layers (bottom -> top):
#   1. Faint raw per-iter trace                                (alpha = alpha_raw)
#   2. +/- 1 rolling-MAD ribbon around the smoothed line       (if show_band)
#   3. Bold rolling-median line (window k)
#   4. Dashed vertical line at the auto-detected stabilization (if any)
#   5. Dotted vertical line at the last large spike            (if any)
#
# `window` overrides the auto window length k = max(50, n/20).
# `y_focus` clips the y-axis to the bulk so a handful of outlier
# iterations don't render the convergence line invisible. Outliers
# remain plotted in the raw layer; they just exit the plot frame.
.plot_dgtf_marglik <- function(x,
                               window = NULL,
                               show_band = TRUE,
                               alpha_raw = 0.25,
                               y_focus = TRUE,
                               y_focus_k = 6,
                               xlab = "Iteration",
                               ylab = NULL,
                               main = NULL,
                               theme = ggplot2::theme_minimal()) {
    ml <- x$fit$marglik
    if (is.null(ml)) {
        stop("No `marglik` trace found in this fit. Convergence plots ",
            "are only available for HVA fits.",
            call. = FALSE
        )
    }
    ml <- as.numeric(ml)
    ml[!is.finite(ml)] <- NA_real_
    n <- length(ml)
    if (n < 5L) {
        stop("`marglik` trace has fewer than 5 points; nothing to plot.",
            call. = FALSE
        )
    }

    # runmed() requires odd window length.
    k <- if (is.null(window)) max(50L, n %/% 20L) else as.integer(window)
    k <- min(k, max(3L, n))
    if (k %% 2L == 0L) k <- k + 1L

    # runmed() requires finite input; impute the global median for the
    # smoother pass only, leaving raw NAs in the displayed line.
    ml_imp <- ml
    if (anyNA(ml_imp)) {
        ml_imp[is.na(ml_imp)] <- stats::median(ml, na.rm = TRUE)
    }

    # Rolling median: robust central-tendency line.
    smth <- as.numeric(stats::runmed(ml_imp, k = k, endrule = "median"))

    # Rolling MAD as a robust noise envelope. Two-pass:
    #   abs_dev   = |x - rolling_median|
    #   mad_band  = 1.4826 * rolling_median(abs_dev)
    # The 1.4826 constant matches stats::mad()'s default and makes the
    # ribbon visually comparable to a +/- 1 SD band under Gaussian
    # noise. Computed once and reused if y_focus needs the same scale.
    sdv <- if (isTRUE(show_band) || isTRUE(y_focus)) {
        abs_dev <- abs(ml_imp - smth)
        1.4826 * as.numeric(stats::runmed(abs_dev,
            k = k,
            endrule = "median"
        ))
    } else {
        NULL
    }

    cv <- .dgtf_marglik_convergence(ml)

    df <- data.frame(
        iter     = seq_len(n),
        marglik  = ml,
        smoothed = smth,
        lower    = if (!is.null(sdv)) smth - sdv else NA_real_,
        upper    = if (!is.null(sdv)) smth + sdv else NA_real_
    )

    plt <- ggplot2::ggplot(df, ggplot2::aes(x = iter)) +
        ggplot2::geom_line(ggplot2::aes(y = marglik),
            alpha = alpha_raw, na.rm = TRUE
        )

    if (isTRUE(show_band)) {
        plt <- plt + ggplot2::geom_ribbon(
            ggplot2::aes(ymin = lower, ymax = upper),
            alpha = 0.15, na.rm = TRUE
        )
    }

    plt <- plt + ggplot2::geom_line(
        ggplot2::aes(y = smoothed),
        linewidth = 0.8, na.rm = TRUE
    )

    if (!is.null(cv) && !is.na(cv$stabilized_iter)) {
        plt <- plt + ggplot2::geom_vline(
            xintercept = cv$stabilized_iter,
            linetype = "dashed", alpha = 0.6
        )
    }
    if (!is.null(cv) && !is.na(cv$last_spike_iter) &&
        # Avoid drawing a second line right on top of the first one.
        (is.na(cv$stabilized_iter) ||
            abs(cv$last_spike_iter - cv$stabilized_iter) > max(5L, cv$roll_win %/% 5L))) {
        plt <- plt + ggplot2::geom_vline(
            xintercept = cv$last_spike_iter,
            linetype = "dotted", alpha = 0.5
        )
    }

    # Y-axis clipping. Use coord_cartesian so points outside the limits
    # are still drawn (lines exit the panel) rather than dropped.
    if (isTRUE(y_focus) && !is.null(sdv)) {
        y_lo <- min(smth - y_focus_k * sdv, na.rm = TRUE)
        y_hi <- max(smth + y_focus_k * sdv, na.rm = TRUE)
        if (is.finite(y_lo) && is.finite(y_hi) && y_lo < y_hi) {
            plt <- plt + ggplot2::coord_cartesian(ylim = c(y_lo, y_hi))
        }
    }

    plt +
        ggplot2::labs(
            x = xlab,
            y = ylab %||%
                expression("log " * hat(p)(y * " | " * gamma) *
                    " (per-iter SMC estimate)"),
            title = main %||% sprintf(
                "HVA marglik trace (rolling median, window = %d)", k
            )
        ) +
        theme
}

#' @export
plot.dgtf_fit <- function(x,
                          what = "Rt",
                          type = c("hist", "trace"),
                          level = 0.95,
                          truth = NULL,
                          time = NULL,
                          color = NULL,
                          alpha = 0.5,
                          xlab = "Time",
                          ylab = NULL,
                          main = NULL,
                          theme = ggplot2::theme_minimal(),
                          hline = NULL,
                          ...) {
    type <- match.arg(type)
    kind <- .what_kind(what)
    if (kind == "param") {
        return(.plot_dgtf_param_single(x, what, truth = truth, type = type, ...))
    }
    if (kind == "diag" && what == "marglik") {
        return(.plot_dgtf_marglik(
            x,
            xlab  = if (identical(xlab, "Time")) "Iteration" else xlab,
            ylab  = ylab,
            main  = main,
            theme = theme,
            ...
        ))
    }

    if (type != "hist") {
        stop("`type = '", type, "'` is only meaningful for static-",
            "parameter plots; got `what = '", what, "'`.",
            call. = FALSE
        )
    }

    if (is.null(hline) && identical(what, "Rt")) hline <- 1
    band <- .dgtf_band(x,
        what = what, level = level, time = time,
        name = "Estimate"
    )
    if (!is.null(truth)) {
        band <- rbind(
            band,
            .dgtf_band(truth,
                what = what,
                time = band$time, name = "Truth"
            )
        )
    }
    .plot_dgtf_band(
        band,
        palette = if (!is.null(color)) c(Estimate = color) else NULL,
        alpha = alpha,
        xlab = xlab,
        ylab = ylab %||% .default_ylab(what),
        main = main,
        theme = theme,
        hline = hline,
        legend.position = if (is.null(truth)) "none" else "right"
    )
}

#' Compare fitted DGTF model plots
#'
#' Overlays posterior bands, forecasts, simulated truth, or parameter summaries
#' from a named list of fitted DGTF-related objects.
#'
#' @param fits Named list of objects accepted by the plotting backend.
#' @param what Quantity to plot, such as `"Rt"`, `"psi"`, `"yhat"`, or a
#'   static parameter name.
#' @param level Credible-interval level.
#' @param time Optional time index.
#' @param palette Optional named color vector.
#' @param alpha Ribbon transparency.
#' @param legend.position,legend.name,legend.nrow Legend controls passed to
#'   ggplot2.
#' @param xlab,ylab,main Axis labels and title.
#' @param theme ggplot2 theme.
#' @param hline Optional horizontal reference line.
#' @param ... Additional arguments for static-parameter plots.
#'
#' @return A ggplot object.
#' @export
dgtf_compare_plot <- function(fits,
                              what = "Rt",
                              level = 0.95,
                              time = NULL,
                              palette = NULL,
                              alpha = 0.5,
                              legend.position = "right",
                              legend.name = "Method",
                              legend.nrow = NULL,
                              xlab = "Time",
                              ylab = NULL,
                              main = NULL,
                              theme = ggplot2::theme_minimal(),
                              hline = NULL,
                              ...) {
    stopifnot(
        is.list(fits), length(fits) >= 1L,
        !is.null(names(fits)), all(nzchar(names(fits)))
    )
    kind <- .what_kind(what)
    if (kind == "param") {
        return(.plot_dgtf_param_compare(fits, what,
            palette = palette, alpha = alpha,
            legend.position = legend.position,
            legend.name = legend.name,
            theme = theme, ...
        ))
    }

    if (is.null(hline) && identical(what, "Rt")) hline <- 1
    bands <- Map(
        function(x, nm) {
            .dgtf_band(x,
                what = what, level = level,
                time = time, name = nm
            )
        },
        fits, names(fits)
    )
    bands <- do.call(rbind, bands)
    .plot_dgtf_band(bands,
        palette         = palette,
        alpha           = alpha,
        xlab            = xlab,
        ylab            = ylab %||% .default_ylab(what),
        main            = main,
        theme           = theme,
        legend.position = legend.position,
        legend.name     = legend.name,
        legend.nrow     = legend.nrow,
        hline           = hline
    )
}

#' @export
plot.dgtf_ppc <- function(x, what = c("Rt", "yhat"),
                          truth = NULL, hline = NULL, ...) {
    what <- match.arg(what)
    band <- .dgtf_band(x, what = what, name = "Estimate")
    if (!is.null(truth)) {
        band <- rbind(
            band,
            .dgtf_band(truth,
                what = what,
                time = band$time, name = "Truth"
            )
        )
    }
    if (is.null(hline) && identical(what, "Rt")) hline <- 1
    .plot_dgtf_band(band, ...,
        legend.position = if (is.null(truth)) "none" else "right",
        ylab = .default_ylab(what),
        hline = hline
    )
}

#' Plot a posterior-predictive forecast
#'
#' Time-series plot of the forecast quantile band (line + ribbon),
#' with an optional truth overlay drawn as a dashed black line. This
#' matches the convention of the other dgtf time-series plotting
#' methods.
#'
#' Single-step forecasts (`only_last = TRUE` in the original
#' [`dgtf_forecast_fit()`] call) are rendered as a point-and-interval
#' instead of a line-and-ribbon, since there's only one (h, y) pair
#' to show; the truth, if present, is drawn as a black "x" at the
#' same time step.
#'
#' @param x A `dgtf_forecast`.
#' @param truth Optional numeric vector of true future values to
#'   overlay (same length convention as `ypred_true` in
#'   [`dgtf_forecast_fit()`]: at least `horizon` entries for the full
#'   trajectory, at least 1 for `only_last`). Defaults to
#'   `x$ypred_true` if it was baked into the forecast object.
#' @param time Optional numeric vector of length `horizon` giving the
#'   x-axis values (e.g. real dates coerced to numeric). Defaults to
#'   `seq_len(horizon)`.
#' @param color Optional color for the forecast band (line + ribbon).
#'   Defaults to black for a single-method forecast plot.
#' @param alpha Ribbon transparency (default 0.5).
#' @param xlab,ylab,main Axis / title labels. NULL falls back to
#'   sensible defaults (`"Forecast step"` and `expression(hat(y))`).
#' @param theme A ggplot2 theme.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#'
#' @export
plot.dgtf_forecast <- function(x,
                               truth = NULL,
                               time = NULL,
                               color = NULL,
                               alpha = 0.5,
                               xlab = NULL,
                               ylab = NULL,
                               main = NULL,
                               theme = ggplot2::theme_minimal(),
                               ...) {
    q <- x$quantiles
    if (is.null(q) || nrow(q) == 0L) {
        stop("Forecast object has no quantiles to plot.", call. = FALSE)
    }

    h <- nrow(q)
    only_last <- isTRUE(x$only_last)

    # x-axis values. only_last collapses to a single time step at
    # x$horizon; otherwise step indices 1..h or user-supplied `time`.
    x_time <- if (only_last) {
        x$horizon
    } else if (is.null(time)) {
        seq_len(h)
    } else {
        if (length(time) != h) {
            stop(sprintf(
                "`time` must have length %d (one entry per forecast step).",
                h
            ), call. = FALSE)
        }
        as.numeric(time)
    }
    xlab <- xlab %||% "Forecast step"
    ylab <- ylab %||% expression(hat(y))

    # Truth: explicit `truth` argument wins, otherwise fall back to
    # the `ypred_true` baked into the forecast object (if any).
    truth_vec <- truth %||% x$ypred_true
    if (!is.null(truth_vec)) {
        truth_vec <- as.numeric(truth_vec)
        truth_vec <- if (only_last) {
            utils::tail(truth_vec, 1L)
        } else if (length(truth_vec) >= h) {
            truth_vec[seq_len(h)]
        } else {
            stop(sprintf(
                "`truth` has length %d but horizon is %d.",
                length(truth_vec), h
            ), call. = FALSE)
        }
    }

    # Single-step path: pointrange + optional truth marker. The
    # multi-step `.plot_dgtf_band()` machinery requires >= 2 points
    # to draw a meaningful line / ribbon.
    if (only_last || h == 1L) {
        df <- data.frame(
            time = as.numeric(x_time),
            lower = as.numeric(q[, "lower"]),
            median = as.numeric(q[, "median"]),
            upper = as.numeric(q[, "upper"]),
            stringsAsFactors = FALSE
        )
        fc_color <- color %||% "black"
        p <- ggplot2::ggplot(df, ggplot2::aes(x = time)) +
            ggplot2::geom_pointrange(
                ggplot2::aes(y = median, ymin = lower, ymax = upper),
                color = fc_color, size = 0.6
            ) +
            theme +
            ggplot2::labs(x = xlab, y = ylab, title = main)
        if (!is.null(truth_vec)) {
            p <- p + ggplot2::geom_point(
                data = data.frame(
                    time = as.numeric(x_time),
                    truth = truth_vec
                ),
                ggplot2::aes(y = truth),
                shape = 4L, size = 3, stroke = 1.1,
                color = "black"
            )
        }
        return(p)
    }

    # Multi-step path: dispatch to the shared band plotter so the
    # forecast looks like every other dgtf time-series plot. Truth
    # rows whose `name` matches `.truth_pattern` ("^(true|truth)$")
    # are auto-styled as a dashed black line by `.plot_dgtf_band()`.
    band <- data.frame(
        time = as.numeric(x_time),
        lower = as.numeric(q[, "lower"]),
        central = as.numeric(q[, "median"]),
        upper = as.numeric(q[, "upper"]),
        name = "Forecast",
        has_ribbon = TRUE,
        stringsAsFactors = FALSE
    )
    if (!is.null(truth_vec)) {
        band <- rbind(band, data.frame(
            time = as.numeric(x_time),
            lower = NA_real_,
            central = truth_vec,
            upper = NA_real_,
            name = "Truth",
            has_ribbon = FALSE,
            stringsAsFactors = FALSE
        ))
    }

    .plot_dgtf_band(
        band,
        palette         = if (!is.null(color)) c(Forecast = color) else NULL,
        alpha           = alpha,
        xlab            = xlab,
        ylab            = ylab,
        main            = main,
        theme           = theme,
        legend.position = if (is.null(truth_vec)) "none" else "right",
        legend.name     = NULL
    )
}

#' In-sample fit + out-of-sample forecast in one plot
#'
#' Plots the posterior-predictive band for `y` over the in-sample
#' window and, optionally, the forecast band for `y` over the
#' out-of-sample window (from a `dgtf_forecast`), distinguished by
#' colour. Truth values for either segment are overlaid as a single
#' dashed black line spanning both windows.
#'
#' The in-sample band can come from either:
#'   * a `dgtf_fit` -- its `attr(., "ppc")` is reused if present,
#'     otherwise [`posterior_predict()`] is run on the fly. The
#'     observed `y` is used as the in-sample truth by default.
#'   * a `dgtf_ppc` -- used directly. No automatic in-sample truth
#'    (the PPC object doesn't carry `y`); pass `truth_in` if you
#'     want one.
#'
#' @param x Either a `dgtf_fit` or a `dgtf_ppc` (the in-sample band
#'   provider).
#' @param forecast Optional `dgtf_forecast` (from
#'   [`dgtf_forecast_fit()`]). When `NULL`, only the in-sample band
#'   is plotted.
#' @param truth_in Optional numeric vector of length `length(fit$y)`
#'   to overlay as the in-sample truth. `NULL` (default) auto-fills
#'   with `fit$y`.
#' @param truth_out Optional numeric vector of length
#'   `forecast$horizon` to overlay as the out-of-sample truth. `NULL`
#'   (default) auto-fills with `forecast$ypred_true` when present.
#' @param show_truth If `FALSE`, suppress all truth overlays
#'   regardless of `truth_in` / `truth_out`. Default `TRUE`.
#' @param time_in,time_out Optional numeric vectors giving the x-axis
#'   values for the two segments (e.g. real dates coerced to
#'   numeric). Defaults: `0:(n-1)` for in-sample and
#'   `n:(n+h-1)` for out-of-sample (with `n = length(fit$y)`,
#'   `h = forecast$horizon`), so the two segments connect.
#' @param level Credible-interval level used when the in-sample PPC
#'   has to be computed on the fly (i.e. `attr(fit, "ppc")` is not
#'   already attached).
#' @param nrep Posterior-predictive replicate count for the
#'   on-the-fly PPC.
#' @param in_color,out_color Hex colours for the two bands. Default
#'   to Okabe-Ito blue / vermillion for visual distinction.
#' @param in_name,out_name Legend labels for the two bands.
#' @param alpha Ribbon transparency (default 0.4).
#' @param xlab,ylab,main Axis / title labels.
#' @param theme A ggplot2 theme.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y[1:180], mod, prior, method = "hva")
#' fc <- dgtf_forecast_fit(fit,
#'     h = 20, nrep = 200,
#'     ypred_true = sim$y[181:200]
#' )
#'
#' # In-sample band + forecast band + dashed truth across the gap.
#' plot_dgtf_y(fit, fc)
#'
#' # In-sample only (forecast = NULL); equivalent to a yhat PPC plot
#' # with the observed series overlaid.
#' plot_dgtf_y(fit)
#'
#' # Real dates on the x-axis.
#' plot_dgtf_y(fit, fc,
#'     time_in  = as.numeric(date_seq[1:180]),
#'     time_out = as.numeric(date_seq[181:200]),
#'     xlab     = "Date"
#' )
#' }
#'
#' @export
plot_dgtf_y <- function(x,
                        forecast = NULL,
                        truth_in = NULL,
                        truth_out = NULL,
                        show_truth = TRUE,
                        time_in = NULL,
                        time_out = NULL,
                        level = 0.95,
                        nrep = 100L,
                        in_color = "#0072B2",
                        out_color = "#D55E00",
                        in_name = "In-sample",
                        out_name = "Forecast",
                        alpha = 0.4,
                        xlab = "Time",
                        ylab = expression(y[t]),
                        main = NULL,
                        theme = ggplot2::theme_minimal(),
                        ...) {
    if (!is.null(forecast) && !inherits(forecast, "dgtf_forecast")) {
        stop("`forecast` must be a `dgtf_forecast` or NULL.",
            call. = FALSE
        )
    }

    # 1. Resolve `x` to (ppc, optional fit). The fit, if available,
    #    is consulted only to auto-fill `truth_in` from `fit$y` --
    #    everything else flows from the PPC. This keeps the bare
    #    `dgtf_ppc` input path fully self-contained.
    fit <- NULL
    if (inherits(x, "dgtf_fit")) {
        fit <- x
        ppc <- attr(fit, "ppc")
        if (is.null(ppc) || !inherits(ppc, "dgtf_ppc") ||
            is.null(ppc$yhat)) {
            ppc <- posterior_predict(fit, nrep = nrep, level = level)
        }
    } else if (inherits(x, "dgtf_ppc")) {
        ppc <- x
    } else {
        stop("`x` must be a `dgtf_fit` or a `dgtf_ppc` (got class: ",
            paste(class(x), collapse = "/"), ").",
            call. = FALSE
        )
    }

    qm_in <- ppc$yhat
    if (is.null(qm_in)) {
        stop("PPC has no `yhat` field. Re-run `posterior_predict()` ",
            "with `nrep > 0` (yhat is only populated when nrep > 0).",
            call. = FALSE
        )
    }

    n_in <- nrow(qm_in)
    t_in <- if (is.null(time_in)) {
        as.numeric(seq_len(n_in) - 1L)
    } else {
        as.numeric(time_in)
    }
    if (length(t_in) != n_in) {
        stop(sprintf(
            "`time_in` must have length %d (the in-sample series length).",
            n_in
        ), call. = FALSE)
    }

    band_in <- data.frame(
        time = t_in,
        lower = as.numeric(qm_in[, 1L]),
        central = as.numeric(qm_in[, 2L]),
        upper = as.numeric(qm_in[, 3L]),
        name = in_name,
        has_ribbon = TRUE,
        stringsAsFactors = FALSE
    )

    # 2. Forecast band (optional). Continue the x-axis from the end
    #    of the in-sample window, using the in-sample step size as
    #    the default increment so dates / integer indices both work.
    band_out <- NULL
    t_out <- numeric(0)
    h_out <- 0L
    if (!is.null(forecast)) {
        qm_out <- forecast$quantiles
        if (is.null(qm_out) || nrow(qm_out) == 0L) {
            stop("Forecast object has no quantiles.", call. = FALSE)
        }
        h_out <- nrow(qm_out)

        if (is.null(time_out)) {
            step <- if (n_in >= 2L) t_in[n_in] - t_in[n_in - 1L] else 1
            t_out <- t_in[n_in] + step * seq_len(h_out)
        } else {
            t_out <- as.numeric(time_out)
            if (length(t_out) != h_out) {
                stop(sprintf(
                    "`time_out` must have length %d (the forecast horizon).",
                    h_out
                ), call. = FALSE)
            }
        }

        band_out <- data.frame(
            time = t_out,
            lower = as.numeric(qm_out[, "lower"]),
            central = as.numeric(qm_out[, "median"]),
            upper = as.numeric(qm_out[, "upper"]),
            name = out_name,
            has_ribbon = TRUE,
            stringsAsFactors = FALSE
        )
    }

    # 3. Truth overlays. Auto-fill defaults; both segments share the
    #    name "Truth" so `.plot_dgtf_band()` draws a single dashed
    #    black line connecting them (visually showing the full
    #    observed-then-future trajectory across the cut).
    if (isTRUE(show_truth)) {
        if (is.null(truth_in) && !is.null(fit)) {
            truth_in <- as.numeric(fit$y)
        }
        if (is.null(truth_out) && !is.null(forecast)) {
            truth_out <- forecast$ypred_true
        }
    } else {
        truth_in <- NULL
        truth_out <- NULL
    }

    truth_band <- NULL
    if (!is.null(truth_in)) {
        truth_in <- as.numeric(truth_in)
        if (length(truth_in) != n_in) {
            stop(sprintf(
                "`truth_in` length %d does not match in-sample length %d.",
                length(truth_in), n_in
            ), call. = FALSE)
        }
        truth_band <- data.frame(
            time = t_in,
            lower = NA_real_,
            central = truth_in,
            upper = NA_real_,
            name = "Truth",
            has_ribbon = FALSE,
            stringsAsFactors = FALSE
        )
    }
    if (!is.null(truth_out) && length(truth_out) > 0L) {
        truth_out <- as.numeric(truth_out)
        if (length(truth_out) < h_out) {
            stop(sprintf(
                "`truth_out` length %d is shorter than horizon %d.",
                length(truth_out), h_out
            ), call. = FALSE)
        }
        truth_out <- truth_out[seq_len(h_out)]
        truth_out_df <- data.frame(
            time = t_out,
            lower = NA_real_,
            central = truth_out,
            upper = NA_real_,
            name = "Truth",
            has_ribbon = FALSE,
            stringsAsFactors = FALSE
        )
        truth_band <- if (is.null(truth_band)) {
            truth_out_df
        } else {
            rbind(truth_band, truth_out_df)
        }
    }

    # 4. Stack all band rows, then dispatch to the shared band
    #    plotter so the visual idiom matches the rest of the
    #    package's time-series plots.
    bands <- band_in
    if (!is.null(band_out)) bands <- rbind(bands, band_out)
    if (!is.null(truth_band)) bands <- rbind(bands, truth_band)

    palette <- stats::setNames(
        c(in_color, out_color),
        c(in_name, out_name)
    )

    .plot_dgtf_band(
        bands,
        palette         = palette,
        alpha           = alpha,
        xlab            = xlab,
        ylab            = ylab,
        main            = main,
        theme           = theme,
        legend.position = "right",
        legend.name     = NULL
    )
}

#' Multi-method in-sample + forecast comparison plot
#'
#' Side-by-side equivalent of [`plot_dgtf_y()`] for an arbitrary number
#' of method/model variants. For each method, an in-sample posterior
#' predictive band is drawn (from the fit's PPC, cached or computed
#' on the fly), and -- when a forecast is supplied -- a forecast band
#' is drawn continuing the time axis past the train/test boundary.
#'
#' Each method gets a single colour applied to *both* its in-sample
#' band and its forecast band; segments are visually distinguished by
#' a vertical dotted line at the train/test split rather than by a
#' second colour channel. This keeps the legend readable across many
#' methods.
#'
#' @param fits Named list. Each entry is a `dgtf_fit` or a `dgtf_ppc`.
#'   `dgtf_fit` entries use their cached PPC if present
#'   (`attr(., "ppc")`), otherwise run [`posterior_predict()`] on the
#'   fly with `nrep` and `level`.
#' @param forecasts Named list of `dgtf_forecast` objects, with names
#'   matching a subset of `names(fits)`. Methods without a forecast
#'   entry are drawn in-sample only. `NULL` (default) draws in-sample
#'   only for every method.
#' @param truth_in Optional numeric vector for the in-sample truth.
#'   `NULL` (default) auto-fills with the first fit's `$y` (assumed
#'   shared across methods).
#' @param truth_out Optional numeric vector for the out-of-sample
#'   truth. `NULL` (default) auto-fills with the first forecast's
#'   `$ypred_true` when present.
#' @param show_truth `FALSE` suppresses all truth overlays.
#' @param time_in,time_out Optional numeric vectors giving the x-axis
#'   values for the two segments. Defaults: `0:(n-1)` for in-sample
#'   and `n:(n+h-1)` for out-of-sample, so the segments connect.
#' @param level,nrep Forwarded to [`posterior_predict()`] for fits
#'   without a cached PPC.
#' @param palette Optional named character vector of colours.
#' @param alpha Ribbon transparency (default 0.4).
#' @param split_line If `TRUE` (default), draw a dotted vertical line
#'   at the train/test boundary. Set `FALSE` to suppress.
#' @param xlab,ylab,main Axis / title labels.
#' @param theme A ggplot2 theme.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#'
#' @export
plot_dgtf_y_compare <- function(fits,
                                forecasts  = NULL,
                                truth_in   = NULL,
                                truth_out  = NULL,
                                show_truth = TRUE,
                                time_in    = NULL,
                                time_out   = NULL,
                                level      = 0.95,
                                nrep       = 100L,
                                palette    = NULL,
                                alpha      = 0.4,
                                split_line = TRUE,
                                xlab       = "Time",
                                ylab       = expression(y[t]),
                                main       = NULL,
                                theme      = ggplot2::theme_minimal(),
                                ...) {
    if (!is.list(fits) || length(fits) == 0L ||
        is.null(names(fits)) || any(!nzchar(names(fits))))
        stop("`fits` must be a named, non-empty list.", call. = FALSE)

    nm <- names(fits)

    if (is.null(forecasts)) forecasts <- list()
    if (!is.list(forecasts))
        stop("`forecasts` must be NULL or a named list.", call. = FALSE)
    extra <- setdiff(names(forecasts), nm)
    if (length(extra))
        stop("`forecasts` has names not in `fits`: ",
             paste(extra, collapse = ", "), call. = FALSE)
    if (length(forecasts) > 0L &&
        !all(vapply(forecasts, inherits, logical(1), "dgtf_forecast")))
        stop("All `forecasts` entries must be `dgtf_forecast` objects.",
             call. = FALSE)

    band_list    <- list()
    in_lengths   <- integer(0)
    t_in_first   <- NULL
    t_out_first  <- NULL
    h_out_first  <- 0L
    fc_seen      <- FALSE

    for (i in seq_along(fits)) {
        name <- nm[i]
        x    <- fits[[i]]

        # Resolve in-sample PPC.
        if (inherits(x, "dgtf_fit")) {
            ppc <- attr(x, "ppc")
            if (is.null(ppc) || !inherits(ppc, "dgtf_ppc") ||
                is.null(ppc$yhat)) {
                ppc <- posterior_predict(x, nrep = nrep, level = level)
            }
        } else if (inherits(x, "dgtf_ppc")) {
            ppc <- x
        } else {
            stop(sprintf(
                "`fits[[\"%s\"]]` must be a `dgtf_fit` or `dgtf_ppc` (got class: %s).",
                name, paste(class(x), collapse = "/")), call. = FALSE)
        }
        qm_in <- ppc$yhat
        if (is.null(qm_in))
            stop(sprintf("`fits[[\"%s\"]]` has no `yhat`. Re-run posterior_predict() with nrep > 0.",
                         name), call. = FALSE)

        n_in <- nrow(qm_in)
        t_in <- if (is.null(time_in)) as.numeric(seq_len(n_in) - 1L)
                else as.numeric(time_in)
        if (length(t_in) != n_in)
            stop(sprintf(
                "`time_in` length %d != in-sample length %d for `%s`.",
                length(t_in), n_in, name), call. = FALSE)
        in_lengths <- c(in_lengths, n_in)
        if (is.null(t_in_first)) t_in_first <- t_in

        band_list[[length(band_list) + 1L]] <- data.frame(
            time       = t_in,
            lower      = as.numeric(qm_in[, 1L]),
            central    = as.numeric(qm_in[, 2L]),
            upper      = as.numeric(qm_in[, 3L]),
            name       = name,         # *same* name as the forecast band
            has_ribbon = TRUE,
            stringsAsFactors = FALSE
        )

        # Optional forecast band.
        if (name %in% names(forecasts)) {
            fc <- forecasts[[name]]
            qm_out <- fc$quantiles
            if (is.null(qm_out) || nrow(qm_out) == 0L) next
            h_out <- nrow(qm_out)
            fc_seen <- TRUE

            t_out <- if (is.null(time_out)) {
                step <- if (n_in >= 2L) t_in[n_in] - t_in[n_in - 1L] else 1
                t_in[n_in] + step * seq_len(h_out)
            } else {
                as.numeric(time_out)
            }
            if (length(t_out) != h_out)
                stop(sprintf(
                    "`time_out` length %d != forecast horizon %d for `%s`.",
                    length(t_out), h_out, name), call. = FALSE)
            if (is.null(t_out_first)) {
                t_out_first <- t_out
                h_out_first <- h_out
            }

            band_list[[length(band_list) + 1L]] <- data.frame(
                time       = t_out,
                lower      = as.numeric(qm_out[, "lower"]),
                central    = as.numeric(qm_out[, "median"]),
                upper      = as.numeric(qm_out[, "upper"]),
                name       = name,         # *same name* -> same color
                has_ribbon = TRUE,
                stringsAsFactors = FALSE
            )
        }
    }

    # Sanity-check in-sample lengths match (otherwise the split line
    # is ambiguous). Warn and use the modal length.
    if (length(unique(in_lengths)) > 1L)
        warning("In-sample lengths differ across fits; using the maximum ",
                "for the split line. Different time axes per method are ",
                "not supported here.", call. = FALSE)
    n_in_eff <- max(in_lengths)
    split_x  <- t_in_first[n_in_eff]

    # Truth overlay. Single "Truth" name keeps the dashed line continuous
    # across the boundary, matching plot_dgtf_y.
    if (isTRUE(show_truth)) {
        if (is.null(truth_in) && inherits(fits[[1L]], "dgtf_fit"))
            truth_in <- as.numeric(fits[[1L]]$y)
        if (is.null(truth_out) && length(forecasts) > 0L)
            truth_out <- forecasts[[1L]]$ypred_true
    } else {
        truth_in <- NULL
        truth_out <- NULL
    }

    truth_rows <- NULL
    if (!is.null(truth_in)) {
        truth_in <- as.numeric(truth_in)
        if (length(truth_in) != n_in_eff)
            stop(sprintf("`truth_in` length %d != in-sample length %d.",
                         length(truth_in), n_in_eff), call. = FALSE)
        truth_rows <- data.frame(
            time = t_in_first, lower = NA_real_,
            central = truth_in, upper = NA_real_,
            name = "Truth", has_ribbon = FALSE,
            stringsAsFactors = FALSE)
    }
    if (!is.null(truth_out) && length(truth_out) > 0L &&
        h_out_first > 0L) {
        truth_out <- as.numeric(truth_out)
        if (length(truth_out) < h_out_first)
            stop(sprintf("`truth_out` length %d < forecast horizon %d.",
                         length(truth_out), h_out_first), call. = FALSE)
        truth_out <- truth_out[seq_len(h_out_first)]
        truth_rows <- rbind(truth_rows, data.frame(
            time = t_out_first, lower = NA_real_,
            central = truth_out, upper = NA_real_,
            name = "Truth", has_ribbon = FALSE,
            stringsAsFactors = FALSE))
    }

    bands <- do.call(rbind, band_list)
    if (!is.null(truth_rows)) bands <- rbind(bands, truth_rows)

    # Build per-method palette. .dgtf_palette handles the named-or-not
    # case; truth gets its own dashed-black styling via .truth_pattern
    # matching inside .plot_dgtf_band.
    method_names <- unique(bands$name[!grepl(.truth_pattern,
                                             bands$name,
                                             ignore.case = TRUE)])
    cols <- .dgtf_palette(method_names, palette)

    p <- .plot_dgtf_band(
        bands,
        palette         = cols,
        alpha           = alpha,
        xlab            = xlab,
        ylab            = ylab,
        main            = main,
        theme           = theme,
        legend.position = "right",
        legend.name     = "Method"
    )

    if (isTRUE(split_line) && fc_seen) {
        p <- p + ggplot2::geom_vline(
            xintercept = split_x,
            linetype   = "dotted",
            color      = "black",
            alpha      = 0.8)
    }
    p
}
