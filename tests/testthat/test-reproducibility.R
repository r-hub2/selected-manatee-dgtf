test_that("set.seed controls the complete HVA calculation", {
    mod <- dgtf_hawkes(
        meanlog = 1.35,
        sigma2 = 0.32,
        W = 0.01,
        seasonality = seas_period(1, init = 2)
    )
    sim <- dgtf_simulate_model(mod, ntime = 25, seed = 7)
    prior <- dgtf_prior(
        W = inv_gamma(1, 1),
        rho = inv_gamma(1, 1)
    )
    control <- vb_control(iter = 3, n_sample = 10, n_particle = 25)

    fit1 <- dgtf(
        sim$y, mod, prior,
        method = "hva", control = control, seed = 123
    )
    fit2 <- dgtf(
        sim$y, mod, prior,
        method = "hva", control = control, seed = 123
    )

    expect_identical(fit1$fit$marglik, fit2$fit$marglik)
    expect_identical(fit1$fit$psi, fit2$fit$psi)
    expect_identical(fit1$fit$W, fit2$fit$W)
})
