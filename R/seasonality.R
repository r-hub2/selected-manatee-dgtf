## Seasonality constructors -----------------------------------------------
##
## Lower to the C++ `season = list(period, in_state, init)` block. The
## C++ engine reads only the first `period` elements of `init` for the
## seasonal pattern; supplying a longer vector is allowed and the extras
## are ignored (this matches the historical raw-list API where users
## occasionally passed long initialisation vectors for state seeding).

new_seasonality <- function(period, init, in_state = FALSE) {
    structure(
        list(
            period   = as.integer(period),
            init     = as.numeric(init),
            in_state = isTRUE(in_state)
        ),
        class = "dgtf_seasonality"
    )
}

#' Seasonality constructors
#'
#' Specify the period, seasonal effects, and whether seasonality enters
#' via the latent state vector or as a fixed deterministic offset.
#'
#' - `seas_none()`: no seasonality (period 1, init 0). The model has
#'   no intercept either; supply one via `dgtf_model(intercept = ...)`.
#' - `seas_weekly(init, in_state)`: weekly seasonality (period 7).
#' - `seas_period(period, init, in_state)`: arbitrary period.
#'
#' `init` is the seasonal-effect vector. A scalar is recycled to length
#' `period`. A vector of length \eqn{\geq} `period` is accepted as-is;
#' only the first `period` elements form the seasonal pattern, but the
#' extras may be used by the C++ engine to seed the latent state when
#' `in_state = TRUE`.
#'
#' Note that when seasonality has period \eqn{> 1} the seasonal effects
#' already encode the level of the series; there is no separate
#' "intercept". `dgtf_model()` will warn if both `intercept` and a
#' non-trivial `seasonality` are set.
#'
#' @param period Integer seasonal period.
#' @param init Seasonal-effect values. Scalar (recycled), or a vector of
#'   length \eqn{\geq} `period`.
#' @param in_state If `TRUE`, treat the seasonal effects as part of the
#'   latent state evolved by the system equation.
#'
#' @return A `dgtf_seasonality` object.
#' @name dgtf-seasonality
#' @examples
#' seas_none()
#' seas_weekly(init = c(0.9, 1.0, 1.1, 1.0, 0.9, 0.8, 1.3))
#' seas_period(period = 2, init = c(1, 2))
NULL

#' @rdname dgtf-seasonality
#' @export
seas_none <- function() new_seasonality(period = 1L, init = 0)

#' @rdname dgtf-seasonality
#' @export
seas_weekly <- function(init = rep(0, 7L), in_state = FALSE) {
    seas_period(period = 7L, init = init, in_state = in_state)
}

#' @rdname dgtf-seasonality
#' @export
seas_period <- function(period, init = rep(0, period), in_state = FALSE) {
    period <- as.integer(period)
    if (length(period) != 1L || period < 1L)
        stop("`period` must be a positive integer.", call. = FALSE)
    if (length(init) == 1L) init <- rep(init, period)
    if (length(init) < period)
        stop(sprintf("`init` must have at least %d elements (one per period).",
                     period),
             call. = FALSE)
    new_seasonality(period = period, init = init, in_state = in_state)
}

#' @export
print.dgtf_seasonality <- function(x, ...) {
    cat(sprintf("<dgtf seasonality: period=%d>\n", x$period))
    if (x$period > 1L) {
        used <- x$init[seq_len(x$period)]
        cat("  init =", paste(format(used), collapse = ", "))
        if (length(x$init) > x$period)
            cat(sprintf(" (+ %d unused tail values)",
                        length(x$init) - x$period))
        cat("\n")
    }
    invisible(x)
}

#' @keywords internal
as_seasonality <- function(x) {
    if (inherits(x, "dgtf_seasonality")) return(x)
    if (is.null(x)) return(NULL)  # signal "no seasonality specified"
    if (is.numeric(x) && length(x) == 1L) return(seas_period(period = x))
    stop("`seasonality` must be a dgtf_seasonality object, NULL, or a single integer period.",
         call. = FALSE)
}

#' @keywords internal
is_trivial_seasonality <- function(s) {
    is.null(s) || (s$period == 1L && isTRUE(all.equal(s$init[1], 0)))
}
