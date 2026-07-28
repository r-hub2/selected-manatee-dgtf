#' Posterior recovery of static parameters against truth
#'
#' Compares per-parameter posterior summaries (mean, sd, credible
#' interval) to user-supplied true values. Useful for simulation
#' studies and SBC-style aggregation across replicates.
#'
#' @param object A `dgtf_fit`.
#' @param truth Named numeric vector OR named list of true static
#'   parameter values. Names follow the same convention as the
#'   columns of `as_draws_matrix(fit)`: scalar params use the bare
#'   name (`"W"`, `"rho"`, `"lag_par1"`, `"lag_par2"`,
#'   `"zintercept"`, `"zzcoef"`); multi-element params use bracketed
#'   indices (`"seas[1]"`, `"seas[2]"`, ...). When passing a list,
#'   vector-valued entries are auto-expanded
#'   (`list(seas = c(0.3, 0.5))` becomes `seas[1] = 0.3, seas[2] = 0.5`).
#'   Unknown names are ignored with a warning.
#' @param level Credible-interval level for the `coverage` column
#'   (default `0.95`).
#' @param digits Number of digits to round the numeric columns to
#'   for display (default `3`).
#'
#' @return A data.frame of class `param_recovery` with one row per
#'   matched parameter and columns `parameter`, `truth`, `mean`,
#'   `median`, `sd`, `q_lo`, `q_hi`, `bias` (`mean - truth`), `mae`,
#'   `rmse`, `crps`, `coverage` (logical: truth in `[q_lo, q_hi]`).
#'   The default print method suppresses the `coverage` column and
#'   instead appends a `Signif` column whose value is `"*"` when the
#'   truth lies inside the credible interval, with the convention
#'   explained in a footer (cf. `print.summary.lm`).
#'   Returns `NULL` if there are no posterior draws or no name
#'   matches.
#'
#' @export
#' @examples
#' \dontrun{
#' fit <- dgtf(sim$y, mod, prior, method = "hva")
#' param_recovery(fit,
#'                truth = list(W    = sim$model$err$par1,
#'                             rho  = sim$model$rho,
#'                             seas = sim$model$seasonality$val))
#' }
param_recovery <- function(object, truth, level = 0.95, digits = 3) {
    if (!inherits(object, "dgtf_fit"))
        stop("`object` must be a `dgtf_fit`.", call. = FALSE)
    if (!is.numeric(level) || length(level) != 1L ||
        level <= 0 || level >= 1)
        stop("`level` must be a single number in (0, 1).", call. = FALSE)

    d <- tryCatch(.dgtf_static_draws(object), error = function(e) NULL)
    if (is.null(d) || is.null(d$draws_matrix) || ncol(d$draws_matrix) == 0L)
        return(NULL)

    dm <- d$draws_matrix
    truth_vec <- .normalize_truth(truth, colnames(dm))
    matched   <- intersect(names(truth_vec), colnames(dm))
    if (length(matched) == 0L) return(NULL)

    sub   <- dm[, matched, drop = FALSE]
    alpha <- 1 - level
    qs    <- apply(sub, 2, stats::quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))

    means   <- apply(sub, 2, mean)
    medians <- apply(sub, 2, stats::median)
    sds     <- apply(sub, 2, stats::sd)
    q_lo    <- qs[1, ]
    q_hi    <- qs[2, ]
    tval    <- truth_vec[matched]
    mae <- vapply(matched,
                  function(p) mean(abs(sub[, p] - tval[p])),
                  numeric(1))
    rmse <- vapply(matched,
                   function(p) sqrt(mean((sub[, p] - tval[p])^2)),
                   numeric(1))
    crps <- vapply(matched,
                   function(p) .crps_sample(sub[, p], tval[p]),
                   numeric(1))

    out <- data.frame(
        parameter = matched,
        truth     = unname(tval),
        mean      = round(unname(means), digits),
        median    = round(unname(medians), digits),
        sd        = round(unname(sds), digits),
        q_lo      = round(unname(q_lo), digits),
        q_hi      = round(unname(q_hi), digits),
        bias      = round(unname(means - tval), digits),
        mae       = round(unname(mae), digits),
        rmse      = round(unname(rmse), digits),
        crps      = round(unname(crps), digits),
        coverage  = unname((tval >= q_lo) & (tval <= q_hi)),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    structure(out,
              class = c("param_recovery", "data.frame"),
              level = level)
}

#' @export
#' @export
print.param_recovery <- function(x, ...) {
    lvl <- attr(x, "level") %||% 0.95
    df  <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
    sig <- ifelse(df$coverage, "*", " ")
    df$coverage <- NULL
    df[[" "]]   <- sig
    print(df, row.names = FALSE, ...)
    cat(sprintf("---\n'*' indicates the truth lies within the %g%% credible interval\n",
                100 * lvl))
    invisible(x)
}

# Coerce `truth` (named numeric or named list) into a named numeric
# vector aligned with the column-name convention of
# `.dgtf_static_draws()$draws_matrix`. Vector-valued list entries are
# expanded as `name[1], name[2], ...`. Unknown names produce a single
# warning so that typos are surfaced.
.normalize_truth <- function(truth, params) {
    if (is.null(truth)) return(numeric(0))

    if (is.numeric(truth) && !is.null(names(truth))) {
        out <- stats::setNames(as.numeric(truth), names(truth))
    } else if (is.list(truth)) {
        flat <- list()
        for (nm in names(truth)) {
            v <- as.numeric(truth[[nm]])
            if (length(v) == 1L) {
                flat[[nm]] <- v
            } else if (length(v) > 1L) {
                for (i in seq_along(v))
                    flat[[paste0(nm, "[", i, "]")]] <- v[i]
            }
        }
        out <- unlist(flat)
    } else {
        stop("`truth` must be a named numeric vector or named list.",
             call. = FALSE)
    }

    extra <- setdiff(names(out), params)
    if (length(extra)) {
        warning(sprintf(
            "Truth value(s) for unknown parameter(s) ignored: %s.\nKnown parameters: %s",
            paste(extra, collapse = ", "),
            paste(params, collapse = ", ")), call. = FALSE)
    }
    out
}

.crps_sample <- function(x, y) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0L || !is.finite(y)) return(NA_real_)
    if (n == 1L) return(abs(x - y))
    xs <- sort.int(x)
    term1 <- mean(abs(xs - y))
    term2 <- sum((2 * seq_len(n) - n - 1) * xs) / (n * n)
    term1 - term2
}
