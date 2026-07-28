test_that("posterior_predict generic dispatches correctly", {
    expect_true(is.function(posterior_predict))
    expect_true(inherits(posterior_predict, "function"))
})

test_that("posterior_predict.dgtf_fit errors on non-fit input", {
    expect_error(posterior_predict.dgtf_fit("not a fit"))
})

test_that("dgtf_forecast_fit validates inputs", {
    expect_error(dgtf_forecast_fit("not a fit"),    "must be a `dgtf_fit`")
    fake_fit <- structure(list(fit = list(), model = NULL,
                               y = numeric(0)),
                          class = "dgtf_fit")
    expect_error(dgtf_forecast_fit(fake_fit, h = 0), "positive integer")
})
