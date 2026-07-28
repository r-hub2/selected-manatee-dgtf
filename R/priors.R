#' Prior distribution constructors
#'
#' Lightweight constructors for the prior distributions accepted by `dgtf()`.
#' Each returns an object of class `dgtf_prior_dist` carrying the
#' distribution name (as the C++ engine recognises it) and a numeric vector
#' of hyperparameters.
#'
#' @details
#' Supported distributions and their (par1, par2) interpretations:
#'
#' - `normal(mean, sd)`
#' - `inv_gamma(shape, scale)`
#' - `dist_gamma(shape, rate)` (renamed to avoid masking `base::gamma`)
#' - `dist_beta(alpha, beta)` (renamed to avoid masking `base::beta`)
#' - `uniform(lower, upper)`
#' - `half_normal(sd)`         (par2 unused; kept for ABI parity)
#' - `half_cauchy(scale)`
#' - `half_t(df, scale)`
#'
#' @param mean,sd,shape,scale,rate,alpha,beta,lower,upper,df Hyperparameters.
#'
#' @return An object of class `dgtf_prior_dist`.
#' @name dgtf-priors
#' @examples
#' normal(0, 1)
#' inv_gamma(2, 1)
NULL

new_prior_dist <- function(name, params) {
    stopifnot(is.character(name), length(name) == 1L)
    stopifnot(is.numeric(params))
    structure(
        list(name = name, params = as.numeric(params)),
        class = c(paste0("dgtf_prior_", name), "dgtf_prior_dist")
    )
}

#' @rdname dgtf-priors
#' @export
normal <- function(mean = 0, sd = 1) {
    if (sd <= 0) stop("`sd` must be positive.", call. = FALSE)
    new_prior_dist("gaussian", c(mean, sd))
}

#' @rdname dgtf-priors
#' @export
inv_gamma <- function(shape = 1, scale = 1) {
    if (shape <= 0 || scale <= 0)
        stop("`shape` and `scale` must be positive.", call. = FALSE)
    new_prior_dist("invgamma", c(shape, scale))
}

#' @rdname dgtf-priors
#' @export
dist_gamma <- function(shape = 1, rate = 1) {
    if (shape <= 0 || rate <= 0)
        stop("`shape` and `rate` must be positive.", call. = FALSE)
    new_prior_dist("gamma", c(shape, rate))
}

#' @rdname dgtf-priors
#' @export
dist_beta <- function(alpha = 1, beta = 1) {
    if (alpha <= 0 || beta <= 0)
        stop("`alpha` and `beta` must be positive.", call. = FALSE)
    new_prior_dist("beta", c(alpha, beta))
}

#' @rdname dgtf-priors
#' @export
uniform <- function(lower = 0, upper = 1) {
    if (lower >= upper)
        stop("`lower` must be strictly less than `upper`.", call. = FALSE)
    new_prior_dist("uniform", c(lower, upper))
}

#' @rdname dgtf-priors
#' @export
half_normal <- function(sd = 1) {
    if (sd <= 0) stop("`sd` must be positive.", call. = FALSE)
    new_prior_dist("halfnormal", c(0, sd))
}

#' @rdname dgtf-priors
#' @export
half_cauchy <- function(scale = 1) {
    if (scale <= 0) stop("`scale` must be positive.", call. = FALSE)
    new_prior_dist("halfcauchy", c(0, scale))
}

#' @rdname dgtf-priors
#' @export
half_t <- function(df = 3, scale = 1) {
    if (df <= 0 || scale <= 0)
        stop("`df` and `scale` must be positive.", call. = FALSE)
    new_prior_dist("halft", c(df, scale))
}

#' @export
print.dgtf_prior_dist <- function(x, ...) {
    cat(sprintf("<dgtf prior: %s(%s)>\n",
                x$name,
                paste(format(x$params, nsmall = 2), collapse = ", ")))
    invisible(x)
}

#' Aggregate priors for inference
#'
#' Collects per-parameter prior specifications. Anything passed here is
#' treated as an unknown to be inferred; anything left out is treated as
#' fixed at the value supplied to `dgtf_model()`.
#'
#' @param intercept Prior on the baseline `a` (regression intercept).
#' @param seasonality Prior on the seasonality effect.
#' @param W Prior on the system-error variance.
#' @param rho Prior on the negative-binomial dispersion.
#' @param lag A named list of priors on the lag-distribution parameters,
#'   e.g. `list(par1 = normal(0, 1), par2 = inv_gamma(1, 1))` (or use
#'   distribution-specific names like `meanlog`/`sdlog`).
#' @param ... Additional named priors on other model parameters (reserved
#'   for future extensions; currently passed through unchanged).
#'
#' @return An object of class `dgtf_prior`.
#' @export
#' @examples
#' dgtf_prior(
#'   seasonality = normal(1, 1),
#'   W           = inv_gamma(1, 1),
#'   rho         = inv_gamma(1, 1),
#'   lag         = list(par1 = normal(0, 1), par2 = inv_gamma(1, 1))
#' )
dgtf_prior <- function(intercept   = NULL,
                       seasonality = NULL,
                       W           = NULL,
                       rho         = NULL,
                       lag         = NULL,
                       ...) {
    extras <- list(...)
    args <- list(
        intercept   = intercept,
        seasonality = seasonality,
        W           = W,
        rho         = rho,
        lag         = lag
    )
    args <- c(args, extras)

    for (nm in names(args)) {
        v <- args[[nm]]
        if (is.null(v)) next
        if (nm == "lag") {
            if (!is.list(v))
                stop("`lag` must be a named list of prior distributions.",
                     call. = FALSE)
            for (sub in names(v)) {
                if (!inherits(v[[sub]], "dgtf_prior_dist"))
                    stop(sprintf("`lag$%s` must be a dgtf prior distribution.",
                                 sub),
                         call. = FALSE)
            }
        } else if (!inherits(v, "dgtf_prior_dist")) {
            stop(sprintf("`%s` must be a dgtf prior distribution (e.g. normal(0, 1)).",
                         nm),
                 call. = FALSE)
        }
    }

    structure(args, class = "dgtf_prior")
}

#' @export
print.dgtf_prior <- function(x, ...) {
    cat("<dgtf_prior>\n")
    nonempty <- Filter(Negate(is.null), x)
    if (length(nonempty) == 0L) {
        cat("  (no priors set; all model parameters fixed)\n")
        return(invisible(x))
    }
    for (nm in names(nonempty)) {
        v <- nonempty[[nm]]
        if (nm == "lag") {
            cat("  lag:\n")
            for (sub in names(v))
                cat(sprintf("    %s = %s(%s)\n",
                            sub, v[[sub]]$name,
                            paste(format(v[[sub]]$params, nsmall = 2),
                                  collapse = ", ")))
        } else {
            cat(sprintf("  %s = %s(%s)\n",
                        nm, v$name,
                        paste(format(v$params, nsmall = 2), collapse = ", ")))
        }
    }
    invisible(x)
}

#' Test whether a parameter has a prior set
#' @keywords internal
has_prior <- function(prior, name) {
    !is.null(prior[[name]])
}
