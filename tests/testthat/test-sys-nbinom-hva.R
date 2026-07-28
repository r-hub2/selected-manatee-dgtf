test_that("HVA can fit sys_nbinom + lag_nbinom (iterative transfer)", {
    skip_on_cran()

    set.seed(20260526)

    mod <- dgtf_model(
        obs  = obs_nbinom(),
        link = link_identity(),
        sys  = sys_nbinom(),
        gain = gain_softplus(),
        lag  = lag_nbinom(r = 3L, kappa = 0.4),
        err  = err_gaussian(W = 0.01),
        seasonality = seas_period(1, init = 2)
    )

    expect_identical(mod$sys$type, "nbinom")
    expect_identical(mod$lag$type, "nbinomp")

    sim <- dgtf_simulate_model(mod, ntime = 150, seed = 7L)
    y <- as.numeric(sim$y)
    expect_true(all(is.finite(y)))
    expect_true(all(y >= 0))

    prior <- dgtf_prior(
        seasonality = normal(0, 5),
        W           = inv_gamma(1, 1),
        rho         = inv_gamma(1, 1)
    )

    fit <- dgtf(
        y       = y,
        model   = mod,
        prior   = prior,
        method  = "hva",
        control = vb_control(iter = 200, n_sample = 200, n_particle = 200)
    )

    expect_s3_class(fit, "dgtf_fit")
    expect_identical(fit$model$sys$type, "nbinom")

    psi_post <- fit$fit$psi
    expect_true(!is.null(psi_post))
    expect_true(all(is.finite(as.numeric(psi_post))))
})


test_that("HVA recovers kappa under sys_nbinom (par1 prior, iterative grad)", {
    skip_on_cran()

    set.seed(20260601)

    true_kappa <- 0.4

    mod <- dgtf_model(
        obs  = obs_nbinom(),
        link = link_identity(),
        sys  = sys_nbinom(),
        gain = gain_softplus(),
        lag  = lag_nbinom(r = 3L, kappa = true_kappa),
        err  = err_gaussian(W = 0.01),
        seasonality = seas_period(1, init = 2)
    )

    sim <- dgtf_simulate_model(mod, ntime = 200, seed = 11L)
    y <- as.numeric(sim$y)

    prior <- dgtf_prior(
        seasonality = normal(0, 5),
        W           = inv_gamma(1, 1),
        rho         = inv_gamma(1, 1),
        lag         = list(par1 = dist_beta(1, 1))
    )

    fit <- dgtf(
        y       = y,
        model   = mod,
        prior   = prior,
        method  = "hva",
        control = vb_control(iter = 1000, n_sample = 500, n_particle = 500)
    )

    expect_s3_class(fit, "dgtf_fit")

    draws <- .dgtf_static_draws(fit)
    expect_true("par1" %in% draws$inferred)
    par1_draws <- draws$draws[["par1"]]
    expect_true(is.numeric(par1_draws) || is.matrix(par1_draws))
    par1_mean <- mean(as.numeric(par1_draws))
    par1_q <- stats::quantile(as.numeric(par1_draws), c(0.025, 0.975))

    expect_true(par1_mean > 0 && par1_mean < 1)
    expect_lt(abs(par1_mean - true_kappa), 0.20)
    expect_true(par1_q[1] < true_kappa && par1_q[2] > true_kappa)
})


test_that("HVA refuses par2 prior under sys_nbinom (structural dimension)", {
    skip_on_cran()

    mod <- dgtf_model(
        obs  = obs_nbinom(),
        link = link_identity(),
        sys  = sys_nbinom(),
        gain = gain_softplus(),
        lag  = lag_nbinom(r = 3L, kappa = 0.4),
        err  = err_gaussian(W = 0.01),
        seasonality = seas_period(1, init = 2)
    )

    sim <- dgtf_simulate_model(mod, ntime = 80, seed = 17L)
    y <- as.numeric(sim$y)

    prior_bad <- dgtf_prior(
        seasonality = normal(0, 5),
        W           = inv_gamma(1, 1),
        rho         = inv_gamma(1, 1),
        lag         = list(par2 = inv_gamma(2, 4))
    )

    expect_error(
        dgtf(y = y, model = mod, prior = prior_bad, method = "hva",
             control = vb_control(iter = 50, n_sample = 50, n_particle = 100)),
        "par2"
    )
})
