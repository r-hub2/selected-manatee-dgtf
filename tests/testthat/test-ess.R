test_that("ess returns ~ n for i.i.d. samples", {
    set.seed(1)
    x <- rnorm(5000)
    expect_gt(ess(x), 4500)
    expect_lte(ess(x), length(x))
})

test_that("ess shrinks for autocorrelated AR(1) samples", {
    set.seed(1)
    x <- as.numeric(stats::arima.sim(list(ar = 0.9), 5000))
    # theoretical IID-equivalent: n * (1 - phi) / (1 + phi) ≈ 263 for n=5000
    expect_lt(ess(x), 600)
    expect_gt(ess(x), 1)
})

test_that("ess handles degenerate input", {
    expect_equal(ess(c(1, 1, 1, 1)), 4)        # gamma0 = 0 fallback
    expect_equal(ess(c(1, NA, 2)), ess(c(1, 2)))
    expect_error(ess("a"))
    expect_error(ess(matrix(rnorm(20), 5, 4)))
})
