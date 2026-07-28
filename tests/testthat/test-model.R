test_that("dgtf_model assembles components and prints", {
    mod <- dgtf_model(
        obs = obs_nbinom(), link = link_identity(),
        sys = sys_shift(), gain = gain_softplus(),
        lag = lag_lognormal(),
        seasonality = seas_weekly()
    )
    expect_s3_class(mod, "dgtf_model")
    out <- capture.output(print(mod))
    expect_true(any(grepl("dgtf_model", out)))
    expect_true(any(grepl("nbinom", out)))
})

test_that("string shortcuts work in dgtf_model", {
    mod <- dgtf_model(obs = "nbinom", link = "identity", sys = "shift",
                      gain = "softplus", lag = "lognorm")
    expect_s3_class(mod$obs, "dgtf_obs_nbinom")
    expect_s3_class(mod$gain, "dgtf_gain_softplus")
})

test_that("intercept lowers to a period-1 seasonality", {
    mod <- dgtf_model(obs_poisson(), link_identity(), sys_identity(),
                      intercept = 5)
    expect_equal(mod$seasonality$period, 1L)
    expect_equal(mod$seasonality$init[1], 5)
})

test_that("intercept and explicit seasonality together warns", {
    expect_warning(
        dgtf_model(obs_poisson(), link_identity(), sys_identity(),
                   intercept   = 5,
                   seasonality = seas_weekly(init = 1:7)),
        "intercept.*ignored"
    )
})

test_that("seas_period accepts longer init vectors (extras ignored)", {
    mod <- dgtf_model(obs_poisson(), link_identity(), sys_identity(),
                      seasonality = seas_period(period = 2, init = 1:7))
    expect_equal(mod$seasonality$period, 2L)
    expect_equal(mod$seasonality$init[1:2], c(1, 2))
})

test_that("as_settings lowers model to the C++ list shape", {
    mod <- dgtf_hawkes()
    s <- as_settings(mod)
    expect_named(s$model, c("obs_dist", "link_func", "sys_eq",
                            "gain_func", "lag_dist", "err_dist"))
    expect_true(is.list(s$param))
    expect_true(is.list(s$season))
    expect_true(is.list(s$zero))
})

test_that("canonical helpers return sensible models", {
    expect_s3_class(dgtf_poisson_ar(p = 3), "dgtf_model")
    expect_s3_class(dgtf_hawkes(), "dgtf_model")
    expect_s3_class(dgtf_koyck(), "dgtf_model")
    expect_s3_class(dgtf_distributed_lag(), "dgtf_model")
})

test_that("dgtf_distributed_lag defaults to the exact iterative Solow form", {
    m <- dgtf_distributed_lag(r = 3, kappa = 0.4)
    expect_identical(m$sys$type, "nbinom")
    expect_identical(m$lag$type, "nbinomp")
})

test_that("dgtf_distributed_lag exposes sys = 'shift' for the truncated form", {
    m_iter  <- dgtf_distributed_lag(r = 3, kappa = 0.4, sys = "nbinom")
    m_shift <- dgtf_distributed_lag(r = 3, kappa = 0.4, sys = "shift")
    expect_identical(m_iter$sys$type, "nbinom")
    expect_identical(m_shift$sys$type, "shift")
    # Both forms keep the NB lag component unchanged.
    expect_identical(m_iter$lag$type,  m_shift$lag$type)
    expect_identical(m_iter$lag$params, m_shift$lag$params)
})

test_that("dgtf_koyck inherits the iterative default and the sys override", {
    m_default <- dgtf_koyck(kappa = 0.5)
    m_shift   <- dgtf_koyck(kappa = 0.5, sys = "shift")
    expect_identical(m_default$sys$type, "nbinom")
    expect_identical(m_shift$sys$type, "shift")
})
