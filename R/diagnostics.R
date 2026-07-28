## Convergence diagnostics: ESS (single chain) and split-Rhat (>= 1 chain).
## Both are implemented natively in src/Diagnostics.h; this file is the
## R-side dispatch.

#' Effective sample size of an MCMC chain
#'
#' Computes the effective sample size of a 1-D vector of posterior draws
#' using Geyer's initial positive sequence applied to the FFT-based
#' autocorrelation function. Non-finite entries are dropped.
#'
#' @param x A numeric vector of posterior draws (a single chain).
#' @return A scalar in `[1, length(x)]`.
#' @examples
#' set.seed(1); ess(rnorm(2000))
#' x <- as.numeric(stats::arima.sim(list(ar = 0.9), 2000)); ess(x)
#' @export
ess <- function(x) {
    if (!is.numeric(x) || !is.null(dim(x)))
        stop("`x` must be a numeric vector.", call. = FALSE)
    ess_cpp(as.numeric(x))
}

#' Split-Rhat (Gelman-Rubin) convergence diagnostic
#'
#' Each input chain of length `n` is split sequentially into halves of
#' length `floor(n/2)` (the middle sample is dropped when `n` is odd), so
#' a single chain produces 2 sub-chains, two chains produce 4, and so on.
#' Values near 1 indicate the chains have mixed; > 1.01 is a common
#' warning threshold.
#'
#' @param x One of:
#'   * a numeric matrix with rows = iterations, cols = chains;
#'   * a list of equal-length numeric vectors (one per chain), giving
#'     posterior draws of a single parameter;
#'   * a list of one or more `dgtf_fit` objects from independent MCMC
#'     runs. HVA fits are rejected — they produce i.i.d. importance
#'     samples, not chains, so split-Rhat is undefined.
#' @return A scalar for matrix / list-of-vectors input. A named numeric
#'   vector (one entry per scalar static parameter) for list-of-fits input.
#' @export
rhat <- function(x) UseMethod("rhat")

#' @export
rhat.default <- function(x) {
    if (is.numeric(x) && is.null(dim(x)))
        stop("`x` must be a numeric matrix (rows = iterations, ",
             "cols = chains) or a list of chains. Pass ",
             "`matrix(x, ncol = 1)` for a single chain.", call. = FALSE)
    if (!is.matrix(x) || !is.numeric(x))
        stop("`x` must be a numeric matrix.", call. = FALSE)
    split_rhat_cpp(x)
}

#' @export
rhat.list <- function(x) {
    if (length(x) == 0L)
        stop("`x` is empty.", call. = FALSE)
    if (all(vapply(x, inherits, logical(1), "dgtf_fit")))
        return(.rhat_fits(x))
    if (all(vapply(x, function(v) is.numeric(v) && is.null(dim(v)),
                   logical(1)))) {
        ns <- vapply(x, length, integer(1))
        if (length(unique(ns)) != 1L)
            stop("Chains have unequal length: ",
                 paste(ns, collapse = ", "),
                 ". Truncate or rerun to equal length first.",
                 call. = FALSE)
        return(split_rhat_cpp(do.call(cbind, x)))
    }
    stop("`x` must be a list of `dgtf_fit` objects or a list of numeric ",
         "vectors (one per chain).", call. = FALSE)
}

# Multi-fit Rhat: bind chains by parameter, intersect inferred sets,
# error on chain-length mismatch or any non-MCMC fit.
.rhat_fits <- function(fits) {
    methods <- vapply(fits, function(f) f$method %||% NA_character_,
                      character(1L))
    if (any(methods != "mcmc", na.rm = TRUE))
        stop("rhat() requires MCMC fits; HVA produces i.i.d. importance ",
             "samples, not chains.", call. = FALSE)

    ds  <- lapply(fits, .dgtf_static_draws)
    dms <- lapply(ds,   `[[`, "draws_matrix")

    ns <- vapply(dms, nrow, integer(1))
    if (length(unique(ns)) != 1L)
        stop("Chains have unequal length: ",
             paste(ns, collapse = ", "),
             ". Rerun fits with equal `iter`.", call. = FALSE)

    common <- Reduce(intersect, lapply(dms, colnames))
    if (length(common) == 0L)
        stop("No static parameters are inferred in all supplied fits.",
             call. = FALSE)
    dropped <- setdiff(Reduce(union, lapply(dms, colnames)), common)
    if (length(dropped))
        warning("Dropping parameters not present in every fit: ",
                paste(dropped, collapse = ", "), call. = FALSE)

    vapply(common, function(p)
        split_rhat_cpp(do.call(cbind,
            lapply(dms, function(m) m[, p, drop = TRUE]))),
        numeric(1))
}
