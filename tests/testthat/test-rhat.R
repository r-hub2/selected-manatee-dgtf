test_that("rhat ~ 1 on i.i.d. samples", {
    set.seed(1)
    M <- matrix(rnorm(4 * 2000), nrow = 2000, ncol = 4)
    expect_lt(abs(rhat(M) - 1), 0.05)
})

test_that("rhat detects between-chain disagreement", {
    set.seed(1)
    M <- cbind(rnorm(2000, mean = 0), rnorm(2000, mean = 1))
    expect_gt(rhat(M), 1.1)
})

test_that("split-Rhat catches within-chain drift on a single chain", {
    set.seed(1)
    x <- c(rnorm(1000), rnorm(1000, mean = 1))           # mid-chain mean shift
    expect_gt(rhat(matrix(x, ncol = 1L)), 1.1)
})

test_that("rhat list-of-vectors form", {
    set.seed(1)
    chains <- list(rnorm(2000), rnorm(2000), rnorm(2000))
    expect_lt(abs(rhat(chains) - 1), 0.05)
    expect_error(rhat(list(rnorm(100), rnorm(50))), "unequal length")
})

test_that("rhat input validation", {
    expect_error(rhat(rnorm(100)), "matrix")
    expect_error(rhat("a"))
    expect_error(rhat(list()), "empty")
})
