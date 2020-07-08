context("Sampler behave well")

test_that("range of parameter draws is somewhat meaningful", {
  lapply(c(res_sv, res_svt, res_svl, res_svl_corrected), function (x) {
    quantile_range <- apply(x$para[, c("mu", "phi", "sigma")], 2, quantile, probs = c(0.1, 0.9))
    expect_gt(quantile_range["10%", "mu"], -30)
    expect_gt(quantile_range["10%", "phi"], -0.5)
    expect_gt(quantile_range["10%", "sigma"], 0.001)
    expect_lt(quantile_range["90%", "mu"], -2)
    expect_lt(quantile_range["90%", "phi"], 0.995)
    expect_lt(quantile_range["90%", "sigma"], 5)
  })
})

test_that("samplers move around", {
  lapply(c(res_sv, res_svt, res_svl, res_svl_corrected), function (x) {
    stdevs <- apply(x$para[, c("mu", "phi", "sigma")], 2, sd)
    expect_gt(stdevs["mu"], 0.1)
    expect_gt(stdevs["phi"], 0.01)
    expect_gt(stdevs["sigma"], 0.01)
  })
})


