context("Sampler behave well")

test_that("range of svsample para draws is somewhat meaningful", {
  lapply(res_sv, function (x) {
    quantile_range <- apply(x$para[, c("mu", "phi", "sigma")], 2, quantile, probs = c(0.1, 0.9))
    expect_gt(quantile_range["10%", "mu"], -30,      label = paste("quantile(mu, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "phi"], -0.5,    label = paste("quantile(phi, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "sigma"], 0.001, label = paste("quantile(sigma, 10%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "mu"], -1,       label = paste("quantile(mu, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "phi"], 0.999,   label = paste("quantile(phi, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "sigma"], 5,     label = paste("quantile(sigma, 90%) with mean model:", x$meanmodel))
  })
})

test_that("range of svtsample para draws is somewhat meaningful", {
  lapply(res_svt, function (x) {
    quantile_range <- apply(x$para[, c("mu", "phi", "sigma", "nu")], 2, quantile, probs = c(0.1, 0.9))
    expect_gt(quantile_range["10%", "mu"], -20,      label = paste("quantile(mu, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "phi"], -0.5,    label = paste("quantile(phi, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "sigma"], 0.001, label = paste("quantile(sigma, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "nu"], 2.01,     label = paste("quantile(nu, 10%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "mu"], 10,       label = paste("quantile(mu, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "phi"], 0.999,   label = paste("quantile(phi, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "sigma"], 5,     label = paste("quantile(sigma, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "nu"], 100,      label = paste("quantile(nu, 90%) with mean model:", x$meanmodel))
  })
})

test_that("range of svlsample para draws is somewhat meaningful", {
  lapply(c(res_svl, res_svl_corrected), function (x) {
    quantile_range <- apply(x$para[, c("mu", "phi", "sigma", "rho")], 2, quantile, probs = c(0.1, 0.9))
    expect_gt(quantile_range["10%", "mu"], -30,      label = paste("quantile(mu, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "phi"], -0.5,    label = paste("quantile(phi, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "sigma"], 0.001, label = paste("quantile(sigma, 10%) with mean model:", x$meanmodel))
    expect_gt(quantile_range["10%", "rho"], -0.9,    label = paste("quantile(rho, 10%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "mu"], -1,       label = paste("quantile(mu, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "phi"], 0.999,   label = paste("quantile(phi, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "sigma"], 5,     label = paste("quantile(sigma, 90%) with mean model:", x$meanmodel))
    expect_lt(quantile_range["90%", "rho"], 0.8,     label = paste("quantile(rho, 90%) with mean model:", x$meanmodel))
  })
})

test_that("samplers move around", {
  lapply(c(res_sv, res_svt, res_svl, res_svl_corrected), function (x) {
  #lapply(c(res_svl, res_svl_corrected), function (x) {
    stdevs <- apply(x$para[, c("mu", "phi", "sigma")], 2, sd)
    expect_gt(stdevs["mu"], 0.1)
    expect_gt(stdevs["phi"], 0.01)
    expect_gt(stdevs["sigma"], 0.01)
  })
})

