test_that("prior constructors produce dgtf_prior_dist objects", {
    expect_s3_class(normal(0, 1), "dgtf_prior_dist")
    expect_s3_class(inv_gamma(2, 1), "dgtf_prior_dist")
    expect_s3_class(dist_gamma(1, 1), "dgtf_prior_dist")
    expect_s3_class(dist_beta(2, 5), "dgtf_prior_dist")
    expect_s3_class(uniform(0, 1), "dgtf_prior_dist")
    expect_s3_class(half_normal(1), "dgtf_prior_dist")
    expect_s3_class(half_cauchy(1), "dgtf_prior_dist")
    expect_s3_class(half_t(3, 1), "dgtf_prior_dist")
})

test_that("prior constructors validate hyperparameters", {
    expect_error(normal(0, -1), "sd")
    expect_error(inv_gamma(-1, 1), "positive")
    expect_error(uniform(1, 0), "less than")
    expect_error(dist_beta(-1, 1), "positive")
})

test_that("dgtf_prior aggregates and validates", {
    p <- dgtf_prior(W = inv_gamma(1, 1), rho = inv_gamma(1, 1))
    expect_s3_class(p, "dgtf_prior")
    expect_true(dgtf:::has_prior(p, "W"))
    expect_false(dgtf:::has_prior(p, "intercept"))

    expect_error(dgtf_prior(W = 1), "prior distribution")
    expect_error(dgtf_prior(lag = list(par1 = 1)), "lag\\$par1")
})
