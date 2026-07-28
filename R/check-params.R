valid_static_params <- c(
    "W", "seas", "rho", "par1", "par2", "zintercept", "zzcoef",
    "intercept", "seasonality"
)

#' Validate static parameter names
#'
#' @param params Character vector of parameter names.
#' @param strict If `TRUE`, unknown names are errors; otherwise they produce a
#'   warning and return `FALSE`.
#'
#' @return `TRUE` when all names are known, otherwise `FALSE` in non-strict mode.
#' @keywords internal
check_params <- function(params, strict = TRUE) {
    params <- as.character(params)
    unknown <- setdiff(params, valid_static_params)
    if (!length(unknown))
        return(TRUE)

    msg <- paste("Unknown parameters:", paste(unknown, collapse = ", "))
    if (isTRUE(strict))
        stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    FALSE
}
