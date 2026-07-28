#' Fit a DGTF model
#'
#' The user-facing entry point for inference on a DGTF model. Wraps the
#' low-level C++ inference engines and returns an S3 `dgtf_fit`
#' object with familiar methods (`print`, `summary`, `coef`, `predict`,
#' ...).
#'
#' @section Choosing a method:
#'
#' - `"hva"`: hybrid variational approximation that combines SMC and
#'   variational updates (Tang et al., the engine featured in the
#'   manuscript). Fast; good for real-time applications. Also tuned via
#'   [`vb_control()`].
#' - `"mcmc"`: disturbance Metropolis-Hastings within Gibbs.
#'   Gold-standard uncertainty quantification; slower. Use
#'   [`mcmc_control()`] to tune.
#'
#' @section Inference rule:
#' Anything supplied a prior in `prior` is treated as unknown and
#' inferred. Anything left out of `prior` is treated as fixed at the
#' value supplied to [`dgtf_model()`].
#'
#' @param y Numeric vector of observed counts (or, in future versions,
#'   a matrix for multivariate observations).
#' @param model A [`dgtf_model`] object.
#' @param prior A [`dgtf_prior`] object (default: empty -> all parameters fixed).
#' @param method `"hva"` or `"mcmc"`.
#' @param control A method-specific control object. If `NULL`, sensible
#'   defaults are used.
#' @param seed Optional integer seed.
#' @param verbose If `TRUE`, returns the lowered model and method
#'   settings *without* calling the C++ engine. Useful for debugging
#'   what gets handed to the engines.
#'
#' @return An object of class `dgtf_fit`.
#' @export
#' @examples
#' \dontrun{
#' mod <- dgtf_hawkes()
#' set.seed(1)
#' sim <- dgtf_simulate(mod, ntime = 200)
#' fit <- dgtf(sim$y, mod,
#'             prior = dgtf_prior(W = inv_gamma(1, 1)),
#'             method = "hva",
#'             control = vb_control(iter = 200))
#' }
dgtf <- function(y,
                 model,
                 prior   = dgtf_prior(),
                 method  = c("hva", "mcmc"),
                 control = NULL,
                 seed      = NULL,
                 verbose   = FALSE) {

    if (!inherits(model, "dgtf_model"))
        stop("`model` must be a `dgtf_model` (see `?dgtf_model`).",
             call. = FALSE)
    if (!inherits(prior, "dgtf_prior"))
        stop("`prior` must be a `dgtf_prior` (see `?dgtf_prior`).",
             call. = FALSE)
    method    <- match.arg(method)
    if (!is.numeric(y))
        stop("`y` must be numeric.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)

    if (is.null(control)) {
        control <- switch(
            method,
            vb            = vb_control(),
            mcmc          = mcmc_control(),
            linear_bayes  = lba_control(),
            smc_control()
        )
    }
    if (!inherits(control, "dgtf_control"))
        stop("`control` must be a `dgtf_control` object.", call. = FALSE)

    model_settings  <- as_settings(model)
    method_settings <- control_to_method_settings(control, prior, method, model)

    if (isTRUE(verbose))
        return(list(model_settings = model_settings,
                    method_settings = method_settings,
                    method = method))

    t0 <- Sys.time()
    raw <- dgtf_infer(
        model_settings = model_settings,
        y_in           = as.numeric(y),
        method         = method,
        method_settings = method_settings
    )
    elapsed <- difftime(Sys.time(), t0, units = "secs")

    elapsed_optimization <- NULL
    elapsed_sampling <- NULL
    if (!is.null(raw$fit$elapsed_opt_us)) {
        elapsed_optimization <- as.difftime(
            raw$fit$elapsed_opt_us / 1e6,
            units = "secs"
        )
        raw$fit$elapsed_opt_us <- NULL
    }
    if (!is.null(raw$fit$elapsed_sample_us)) {
        elapsed_sampling <- as.difftime(
            raw$fit$elapsed_sample_us / 1e6,
            units = "secs"
        )
        raw$fit$elapsed_sample_us <- NULL
    }

    structure(
        list(
            fit     = raw$fit,
            y       = as.numeric(y),
            model   = model,
            prior   = prior,
            control = control,
            method  = method,
            elapsed = elapsed,
            elapsed_optimization = elapsed_optimization,
            elapsed_sampling     = elapsed_sampling
        ),
        class = "dgtf_fit"
    )
}
