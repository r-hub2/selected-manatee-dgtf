test_that("component constructors return tagged objects", {
    expect_s3_class(obs_nbinom(),     c("dgtf_obs_nbinom", "dgtf_obs"))
    expect_s3_class(link_identity(),  c("dgtf_link_identity", "dgtf_link"))
    expect_s3_class(gain_softplus(),  c("dgtf_gain_softplus", "dgtf_gain"))
    expect_s3_class(lag_lognormal(),  c("dgtf_lag_lognorm", "dgtf_lag"))
    expect_s3_class(sys_shift(),      c("dgtf_sys_shift", "dgtf_sys"))
    expect_s3_class(err_gaussian(),   c("dgtf_err_gaussian", "dgtf_err"))
})

test_that("string shortcuts coerce to component objects", {
    expect_s3_class(dgtf:::as_component("nbinom",  "obs"), "dgtf_obs_nbinom")
    expect_s3_class(dgtf:::as_component("softplus", "gain"), "dgtf_gain_softplus")
    expect_error(dgtf:::as_component("not_a_link", "link"), "Unknown")
})

test_that("lag constructors validate parameters", {
    expect_equal(
        lag_lognormal(sdlog = 0.5)$params$par2,
        0.25
    )
    expect_equal(
        lag_lognormal(sigma2 = 0.25),
        lag_lognormal(sdlog = 0.5)
    )
    expect_error(
        lag_lognormal(sigma2 = 0.25, sdlog = 0.5),
        "only one"
    )
    expect_error(lag_lognormal(sdlog = -1), "positive")
    expect_error(lag_nbinom(kappa = 1.5), "0, 1")
    expect_error(lag_uniform(window = 0), "positive")
})

test_that("seasonality constructors enforce length matching", {
    expect_s3_class(seas_weekly(),                 "dgtf_seasonality")
    expect_s3_class(seas_period(7),                "dgtf_seasonality")
    expect_s3_class(seas_weekly(init = 5),         "dgtf_seasonality")
    expect_s3_class(seas_weekly(init = 1:7),       "dgtf_seasonality")

    # Vectors longer than the period are accepted (extras are state-init
    # and ignored by the seasonal pattern).
    s <- seas_period(period = 2, init = 1:7)
    expect_s3_class(s, "dgtf_seasonality")
    expect_equal(length(s$init), 7L)
    expect_equal(s$init[1:2], c(1, 2))

    expect_error(seas_period(7, init = 1:5), "at least 7")
    expect_error(seas_period(0),             "positive integer")
})
