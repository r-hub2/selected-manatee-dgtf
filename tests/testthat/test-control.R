test_that("control constructors return tagged objects", {
    expect_s3_class(vb_control(),    c("dgtf_vb_control", "dgtf_control"))
    expect_s3_class(mcmc_control(),  c("dgtf_mcmc_control", "dgtf_control"))
    expect_s3_class(smc_control(),   c("dgtf_smc_control", "dgtf_control"))
    expect_s3_class(lba_control(),   c("dgtf_lba_control", "dgtf_control"))
})

test_that("lba_control validates discount factor", {
    expect_error(lba_control(0),   "must lie")
    expect_error(lba_control(1.1), "must lie")
})

test_that("smc_control matches resample arg", {
    expect_error(smc_control(resample = "weird"), "should be one of")
})
