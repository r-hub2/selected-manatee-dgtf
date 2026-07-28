## Method-specific control objects ----------------------------------------
##
## Each `<method>_control()` returns a validated list inheriting from
## `dgtf_control`. The values are merged into the corresponding
## `dgtf_default_algo_settings(method)` C++ defaults at fit time, so the
## user only needs to override the fields they care about.

new_control <- function(method, opts) {
    structure(
        c(list(method = method), opts),
        class = c(paste0("dgtf_", method, "_control"), "dgtf_control")
    )
}

#' Hybrid variational-Bayes control
#'
#' Tuning parameters for the hybrid variational algorithm
#' (`method = "vb"`, alias for the C++ `"hva"` engine).
#'
#' @param iter Number of variational iterations.
#' @param n_sample Number of posterior draws to retain after fitting.
#' @param n_particle Number of particles in the inner SMC step.
#' @param learning_rate Step size for the natural-gradient update.
#' @param eps_step_size Numerical step used for finite-difference checks.
#' @param ... Additional named arguments passed through verbatim to the
#'   C++ engine, useful for advanced experimentation.
#'
#' @return A `dgtf_vb_control` object.
#' @export
#' @examples
#' vb_control(iter = 500, n_sample = 1000)
vb_control <- function(iter          = 1000L,
                       n_sample      = 1000L,
                       n_particle    = 200L,
                       learning_rate = 0.02,
                       eps_step_size = 1e-5,
                       ...) {
    new_control("vb", list(
        niter          = as.integer(iter),
        nsample        = as.integer(n_sample),
        num_particle   = as.integer(n_particle),
        learning_rate  = as.numeric(learning_rate),
        eps_step_size  = as.numeric(eps_step_size),
        extra          = list(...)
    ))
}

#' MCMC control
#'
#' Tuning parameters for the disturbance-MCMC sampler
#' (`method = "mcmc"`).
#'
#' @param iter Number of MCMC iterations (post-warmup retained).
#' @param warmup Burn-in iterations, discarded.
#' @param thin Thinning interval.
#' @param ... Forwarded to the C++ engine.
#'
#' @return A `dgtf_mcmc_control` object.
#' @export
mcmc_control <- function(iter   = 2000L,
                         warmup = 1000L,
                         thin   = 1L,
                         ...) {
    new_control("mcmc", list(
        nsample = as.integer(iter),
        nburnin = as.integer(warmup),
        nthin   = as.integer(thin),
        extra   = list(...)
    ))
}

#' Sequential Monte Carlo control
#'
#' Tuning parameters for SMC-based smoothers
#' (`method = "smc"`, `"ffbs"`, `"tfs"`, `"mcs"`, `"pl"`).
#'
#' @param n_particle Number of particles.
#' @param n_backward Number of backward simulations (smoothing variants).
#' @param resample Resampling strategy: `"systematic"` (default),
#'   `"multinomial"`, or `"stratified"`.
#' @param ess_threshold Effective-sample-size threshold triggering resampling.
#' @param ... Forwarded to the C++ engine.
#'
#' @return A `dgtf_smc_control` object.
#' @export
smc_control <- function(n_particle    = 500L,
                        n_backward    = 100L,
                        resample      = c("systematic", "multinomial",
                                          "stratified"),
                        ess_threshold = 0.5,
                        ...) {
    resample <- match.arg(resample)
    new_control("smc", list(
        num_particle  = as.integer(n_particle),
        num_backward  = as.integer(n_backward),
        resample      = resample,
        ess_threshold = as.numeric(ess_threshold),
        extra         = list(...)
    ))
}

#' Linear-Bayes control
#'
#' Tuning parameters for the linear-Bayes approximation
#' (`method = "linear_bayes"`).
#'
#' @param discount_factor Discount factor in (0, 1].
#' @param ... Forwarded to the C++ engine.
#'
#' @return A `dgtf_lba_control` object.
#' @export
lba_control <- function(discount_factor = 0.95, ...) {
    if (discount_factor <= 0 || discount_factor > 1)
        stop("`discount_factor` must lie in (0, 1].", call. = FALSE)
    new_control("lba", list(
        discount_factor = as.numeric(discount_factor),
        extra           = list(...)
    ))
}

#' @export
print.dgtf_control <- function(x, ...) {
    cat(sprintf("<dgtf_control: %s>\n", x$method))
    for (nm in setdiff(names(x), c("method", "extra"))) {
        v <- x[[nm]]
        if (length(v) > 1L) {
            cat(sprintf("  %s = c(%s)\n", nm,
                        paste(format(v), collapse = ", ")))
        } else {
            cat(sprintf("  %s = %s\n", nm, format(v)))
        }
    }
    if (length(x$extra))
        cat("  extra:", paste(names(x$extra), collapse = ", "), "\n")
    invisible(x)
}

# ---- Lowering control + prior to the C++ method_settings list -----------
#
# Each prior block in the C++ engine has four fields:
#     list(infer = TRUE/FALSE,
#          mh_sd = <numeric>,
#          prior_name = "<dist>",
#          prior_param = c(...))
# We MUST merge into the C++ defaults (not replace), otherwise dropped
# fields like `mh_sd` cause segfaults when the engine reads them.

#' @keywords internal
control_to_method_settings <- function(control, prior, method, model) {
    method_lc <- tolower(method)
 
    # Start from the C++ defaults for this method
    opts <- dgtf_default_algo_settings(method_lc)
 
    # Merge in the values from the control object
    for (nm in setdiff(names(control), c("method", "extra"))) {
        opts[[nm]] <- control[[nm]]
    }
    for (nm in names(control$extra %||% list())) {
        opts[[nm]] <- control$extra[[nm]]
    }
 
    # ---- Guard 1: discount factor + W are mutually exclusive ----------
    #
    # When use_discount = TRUE, W_t is computed adaptively as
    #   W_t = (1/delta - 1) * G_t * C_{t-1} * G_t'
    # There is no fixed W to infer. Inferring W alongside discount
    # produces contradictory updates.
    use_discount <- isTRUE(opts$use_discount)
    if (use_discount && has_prior(prior, "W")) {
        warning(
            "Prior on W is ignored when `use_discount = TRUE`: the ",
            "innovation variance is computed adaptively from the discount ",
            "factor. Remove `W` from the prior specification.",
            call. = FALSE
        )
    }
 
    # ---- Guard 2: AR models have structural lag parameters ------------
    #
    # For AR models (lag = "uniform", sys = "identity"), par1 is the AR
    # order and par2 is unused. Both are integers that set the state
    # dimension nP. Inferring them changes nP mid-run, crashing the SMC.
    lag_is_structural <- identical(model$lag$type, "uniform") ||
                         identical(model$sys$type, "identity")
 
    if (lag_is_structural && !is.null(prior$lag)) {
        lag_requested <- !is.null(prior$lag$par1) || !is.null(prior$lag$par2)
        if (lag_requested) {
            warning(
                "Lag priors are ignored for AR models: the lag parameters ",
                "set the autoregressive order (a structural dimension), ",
                "not continuous distributional parameters. ",
                "Remove `lag` from the prior specification.",
                call. = FALSE
            )
        }
    }
 
    # ---- Count inferred parameters ------------------------------------
    inferred <- character(0)
    if (has_prior(prior, "intercept"))   inferred <- c(inferred, "intercept")
    if (has_prior(prior, "seasonality")) inferred <- c(inferred, "seas")
 
    # W: infer only if prior is set AND discount is not active
    if (has_prior(prior, "W") && !use_discount) inferred <- c(inferred, "W")
 
    if (has_prior(prior, "rho"))         inferred <- c(inferred, "rho")
 
    # Lag parameters: infer only if not structural
    infer_par1 <- FALSE
    infer_par2 <- FALSE
    if (!lag_is_structural && !is.null(prior$lag)) {
        if (!is.null(prior$lag$par1)) {
            inferred <- c(inferred, "par1")
            infer_par1 <- TRUE
        }
        if (!is.null(prior$lag$par2)) {
            inferred <- c(inferred, "par2")
            infer_par2 <- TRUE
        }
    }
 
    # Compute k (rank of variational factor matrix B)
    opts$k <- length(inferred)
 
    # Seasonal period > 1 adds (period - 1) extra dimensions beyond
    # the single "seas" slot already counted.
    if (("seas" %in% inferred) && model$seasonality$period > 1L)
        opts$k <- opts$k + model$seasonality$period - 2L
 
    # When both lag parameters are inferred jointly, they share a
    # factor dimension in the variational family.
    if (infer_par1 && infer_par2)
        opts$k <- opts$k - 1L
 
    opts$k <- max(opts$k, 1L)
 
    # ---- Per-parameter prior blocks -----------------------------------
    #
    # MERGE into the engine defaults (which carry e.g. mh_sd) rather
    # than replacing them entirely.
    update_block <- function(block_name, dist, force_off = FALSE) {
        block <- opts[[block_name]]
        if (!is.list(block)) block <- list()
        if (is.null(dist) || force_off) {
            block$infer <- FALSE
        } else {
            block$infer <- TRUE
            block$prior_name  <- dist$name
            block$prior_param <- dist$params
        }
        opts[[block_name]] <<- block
    }
 
    update_block("seas", prior$seasonality)
    update_block("W",    prior$W,    force_off = use_discount)
    update_block("rho",  prior$rho)
    update_block("par1", prior$lag$par1, force_off = lag_is_structural)
    update_block("par2", prior$lag$par2, force_off = lag_is_structural)
 
    opts
}
