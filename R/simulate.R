#' Simulate from a DGTF model
#'
#' Draws an observation series and the corresponding latent trajectories
#' from a [`dgtf_model()`].
#'
#' @param model A `dgtf_model` object.
#' @param ntime Length of the simulated series.
#' @param seed Optional RNG seed.
#' @param y0 Initial observation (default 0).
#' @param z Optional binary mask of length `ntime + 1` enforcing
#'   structural zeros at the indicated time points. Meaningful only
#'   for zero-inflated models. Positions with `z = 0` are forced to
#'   `y = 0`, while positions with `z = 1` are simulated from the
#'   observation distribution as usual. Default `NULL` (no mask).
#'
#' @return An object of class `dgtf_sim` (a list) with components
#'   `y`, `lambda`, `psi`, `Theta`, `ft`, `nlag`, plus the originating
#'   `model` and the gain-function name `gain` used to map `psi` to
#'   `R_t`. When the model is zero-inflated, `prob` and `z` are also
#'   present.
#' @export
#' @examples
#' mod <- dgtf_hawkes()
#' sim <- dgtf_simulate_model(mod, ntime = 200, seed = 1)
#' print(sim)
#' plot(sim)
#' plot(sim, what = "Rt")
dgtf_simulate_model <- function(model, ntime, seed = NULL, y0 = 0, z = NULL) {
    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model`.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    settings <- as_settings(model)
    out <- dgtf_simulate(
        settings = settings,
        ntime    = as.integer(ntime),
        y0       = as.numeric(y0),
        z        = if (is.null(z)) NULL else as.numeric(z)
    )
    out$model <- model
    out$gain  <- model$gain$type %||% "identity"
    class(out) <- c("dgtf_sim", "list")
    out
}

#' @export
print.dgtf_sim <- function(x, digits = 3L, ...) {
    y   <- as.numeric(x$y)
    psi <- as.numeric(x$psi)
    h   <- .gain_fun(x$gain %||% "identity")
    Rt  <- h(psi)

    cat("<dgtf simulated series>\n")
    cat(sprintf("Length      : %d  (ntime = %d)\n",
                length(y), length(y) - 1L))
    cat(sprintf("Gain        : %s\n", x$gain %||% "?"))
    if (!is.null(x$nlag))
        cat(sprintf("Lag window  : %d\n", as.integer(x$nlag)))
    cat(sprintf("y           : range [%g, %g], mean %s\n",
                min(y), max(y),
                formatC(mean(y), digits = digits, format = "g")))
    cat(sprintf("R_t         : range [%s, %s], mean %s\n",
                formatC(min(Rt),  digits = digits, format = "g"),
                formatC(max(Rt),  digits = digits, format = "g"),
                formatC(mean(Rt), digits = digits, format = "g")))
    if (!is.null(x$z)) {
        zv <- as.numeric(x$z)
        cat(sprintf("Zero-infl.  : %d / %d states are zero-inflated\n",
                    sum(zv == 0), length(zv)))
    }
    invisible(x)
}


#' Plot a simulated DGTF series
#'
#' Visualises the observed series `y` or the time-varying intensity
#' `R_t = h(psi_t)` from a [`dgtf_simulate_model()`] result.
#'
#' @param x A `dgtf_sim` object.
#' @param what One of `"y"` (default) or `"Rt"`.
#' @param split Optional numeric scalar or vector of time indices at
#'   which to draw a vertical dotted line. Useful for marking a
#'   train/test boundary in simulation studies. Default `NULL`
#'   (no line drawn).
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
plot.dgtf_sim <- function(x, what = c("y", "Rt"), split = NULL, ...) {
    what <- match.arg(what)
    y <- as.numeric(x$y)
    t <- seq_along(y) - 1L

    if (what == "y") {
        df <- data.frame(time = t, value = y)
        p  <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
            ggplot2::geom_line(color = "black") +
            ggplot2::geom_point(color = "black", size = 0.8) +
            ggplot2::labs(x = "t", y = "y")
    } else {
        psi <- as.numeric(x$psi)
        h   <- .gain_fun(x$gain %||% "identity")
        df  <- data.frame(time = t, value = h(psi))
        p   <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
            ggplot2::geom_hline(yintercept = 1, color = "grey50", linetype = "dashed") +
            ggplot2::geom_line(color = "black") +
            ggplot2::labs(x = "t", y = expression(R[t]))
    }
    if (!is.null(split)) {
        p <- p + ggplot2::geom_vline(
            xintercept = as.numeric(split),
            color = "grey50", linetype = "dotted")
    }
    p + ggplot2::theme_minimal()
}

#' Coerce a list of simulated data to a `dgtf_sim` object
#'
#' Wraps externally-generated simulated data (e.g. from another
#' simulator, or a hand-built study) into a `dgtf_sim` so the
#' existing \code{print()} and \code{plot()} methods
#' methods, and any model-check workflows that key off `psi` /
#' `model` as the ground truth, work on it.
#'
#' @param x A list with at least an element `y` (numeric vector of
#'   observations). May also include `psi`, `lambda`, `ft`, `prob`,
#'   `z` (numeric vectors of the same length as `y`), `Theta`
#'   (numeric matrix), and `nlag` (integer scalar).
#' @param model A `dgtf_model` whose **fixed parameter values are
#'   the ground-truth values used to generate `x`**. Its gain
#'   function is also what's used to derive `R_t` from `psi` in
#'   `plot(., what = "Rt")`.
#' @param ... Unused.
#'
#' @return An object of class `dgtf_sim`.
#'
#' @section Important:
#'   The rest of the package treats `psi` and the parameters in
#'   `model` as the **truth** for model-checking purposes (e.g.
#'   `plot(fit, what = "Rt", truth = log1p(exp(sim$psi)))`,
#'   coverage / interval-score calculations in
#'   [`posterior_predict()`]). Pass only values that actually came
#'   from a forward simulation under the same `model`. Arbitrary
#'   numbers will produce meaningless model-checking comparisons.
#'
#' @export
#' @examples
#' \dontrun{
#' # External simulation: any code that draws from a generative
#' # process described by `mod`. Here we just show the shape:
#' mod <- dgtf_hawkes()
#' raw <- list(y = my_y, psi = my_psi)   # at minimum, `y`
#' sim <- as_dgtf_sim(raw, model = mod)
#' plot(sim, what = "Rt")
#' }
as_dgtf_sim <- function(x, model, ...) {
    if (!is.list(x))
        stop("`x` must be a list.", call. = FALSE)
    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model`.", call. = FALSE)
    if (is.null(x$y) || !is.numeric(x$y))
        stop("`x$y` must be a numeric vector.", call. = FALSE)

    n <- length(as.numeric(x$y))

    chk_vec <- function(name) {
        v <- x[[name]]
        if (is.null(v)) return(NULL)
        if (!is.numeric(v))
            stop(sprintf("`x$%s` must be numeric.", name),
                 call. = FALSE)
        nv <- length(as.numeric(v))
        if (nv != n)
            stop(sprintf(
                "`x$%s` has length %d but `x$y` has length %d.",
                name, nv, n), call. = FALSE)
        v
    }

    psi    <- chk_vec("psi")
    lambda <- chk_vec("lambda")
    ft     <- chk_vec("ft")
    prob   <- chk_vec("prob")
    z      <- chk_vec("z")

    Theta <- x$Theta
    if (!is.null(Theta) && !is.numeric(Theta))
        stop("`x$Theta` must be a numeric matrix.", call. = FALSE)

    out <- list(
        y      = x$y,
        nlag   = x$nlag,
        psi    = psi,
        Theta  = Theta,
        ft     = ft,
        lambda = lambda
    )
    if (!is.null(prob)) out$prob <- prob
    if (!is.null(z))    out$z    <- z

    out$model <- model
    out$gain  <- model$gain$type %||% "identity"
    class(out) <- c("dgtf_sim", "list")

    message(
        "Note: `dgtf_sim` is treated as ground truth by the rest of ",
        "the package (`psi` is read as the true latent path; ",
        "parameters in `model` are read as the true static values, ",
        "e.g. for coverage / RMSE in `posterior_predict()`). Use only ",
        "values from a faithful forward simulation under `model` -- ",
        "arbitrary numbers will produce meaningless model-checking ",
        "comparisons.")
    out
}
