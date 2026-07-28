test_that("package loads and namespace is registered", {
    expect_true("dgtf" %in% loadedNamespaces())
})

test_that("check_params validates known parameter names", {
    expect_true(dgtf:::check_params(c("seas", "rho"), strict = FALSE))
    expect_warning(
        out <- dgtf:::check_params(c("seas", "not_a_param"), strict = FALSE),
        "Unknown parameters"
    )
    expect_false(out)
    expect_error(
        dgtf:::check_params(c("seas", "not_a_param"), strict = TRUE),
        "Unknown parameters"
    )
})
