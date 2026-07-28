#' Default simulated count series with ground-truth model and latents
#'
#' A 200-step series simulated from a Hawkes-style DGTF model. Ships
#' with the originating model and the true latent `psi` path so it
#' can be used as a ground-truth reference for posterior checks.
#'
#' @format A [`dgtf_sim`][dgtf_simulate_model] object: a list
#'   containing `y` (numeric vector of observations), `psi`, `lambda`,
#'   `ft`, `Theta`, `nlag`, plus the originating `model` and `gain`.
#' @source Generated with `dgtf_simulate_model()` using seed 1.
#' @keywords datasets
"sim_hawkes"


#' Daily new COVID-19 cases, Santa Cruz County, California
#'
#' Daily counts of newly reported COVID-19 cases in Santa Cruz
#' County, California from July 1, 2020 through December 1, 2021.
#'
#' @format A numeric vector of length 518.
#' @source California Department of Public Health,
#'   \url{https://lab.data.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state-archived}.
#' @keywords datasets
"covid_scz_2020"
