## Bridges to the wider Bayesian R ecosystem ------------------------------
##
## These methods let users hand a `dgtf_fit` to packages like `posterior`
## (`as_draws_array`), `loo` (`log_lik`/`loo`), and `bayesplot`. They are
## thin extractors over whatever the C++ engine deposited in `fit$samples`
## and `fit$pointwise_log_lik`.

#' Coerce a fitted DGTF model to a posterior::draws object
#'
#' Returns a `draws_array` (or a `draws_matrix` when only one chain is
#' available) of static-parameter posterior samples. Latent-state
#' trajectories are exposed via [`extract_state_draws()`] for memory
#' reasons.
#'
#' @param x A `dgtf_fit`.
#' @param ... Unused.
#'
#' @return A `posterior::draws_array` if the **posterior** package is
#'   installed and samples are available; otherwise a plain matrix.
#' @export
as_draws <- function(x, ...) UseMethod("as_draws")

#' @rdname as_draws
#' @export
as_draws.dgtf_fit <- function(x, ...) {
    d <- .dgtf_static_draws(x)
    if (is.null(d$draws_matrix) || ncol(d$draws_matrix) == 0L) {
        stop("This fit has no posterior samples for static parameters.",
            call. = FALSE
        )
    }
    if (requireNamespace("posterior", quietly = TRUE)) {
        return(posterior::as_draws_matrix(d$draws_matrix))
    }
    d$draws_matrix
}

#' Pointwise log-likelihood
#'
#' Per-observation log-likelihood matrix with rows = posterior draws,
#' columns = observations. Hand the result to `loo::loo()`.
#'
#' @param object A `dgtf_fit`.
#' @param ... Unused.
#' @export
log_lik <- function(object, ...) UseMethod("log_lik")

#' @rdname log_lik
#' @export
log_lik.dgtf_fit <- function(object, ...) {
    ll <- object$fit$pointwise_log_lik %||% object$fit$log_lik_matrix
    if (is.null(ll))
        stop("`pointwise_log_lik` not stored for this method.", call. = FALSE)
    ll
}

#' Latent-state posterior draws
#'
#' Returns an `n x nsample` matrix for `psi`, or a 3-D array (`nP x (n+1) x nsample`) for `Theta`.
#'
#' @param object A `dgtf_fit`.
#' @param what One of `"psi"` (the scalar `Rt`-driver) or `"Theta"`
#'   (full state matrix).
#' @export
extract_state_draws <- function(object, what = c("psi", "Theta")) {
    what <- match.arg(what)
    out <- switch(what,
        psi   = object$fit$psi_stored, # n x nsample
        Theta = object$fit$Theta_stored %||% object$fit$Theta
    )
    if (is.null(out)) {
        stop(sprintf("State component `%s` not found in fit output.", what),
            call. = FALSE
        )
    }
    out
}