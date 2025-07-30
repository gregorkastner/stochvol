test_that("vanilla SV passes Geweke test", {
  set.seed(60)
  result <- .Call(`_stochvol_geweke_fast_cpp`, PACKAGE = "stochvol")
  draws <- result$draws

  for (centered in c(FALSE, TRUE)) {
    priorspec <-
      if (centered) {  # pro-centered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.9),
                       phi = sv_beta(shape1 = 2, shape2 = 1.5),
                       sigma2 = sv_gamma(shape = 0.9, rate = 0.9))
      } else {  # pro-noncentered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                       phi = sv_beta(shape1 = 5, shape2 = 1.5),
                       sigma2 = sv_gamma(shape = 0.9, rate = 9))
      }
    centered_text <- if (centered) "centered" else "noncentered"
    store_para <- t(result[[centered_text]]$para)
    colnames(store_para) <- c("mu", "phi", "sigma")

    thin_skip <- ceiling(draws / min(c(4500, apply(store_para, 2, coda::effectiveSize))))
    thin_index <- seq(1, draws, by = thin_skip)

    expect_gt(shapiro.test((sample(c(-1, 1), draws, replace = TRUE) * store_para[, "sigma"] * sqrt(2 * priorspec$sigma2$rate))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test(((store_para[, "mu"] - priorspec$mu$mean) / priorspec$mu$sd)[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)))[thin_index])$p.value, 1e-5)

    # visual tests for manual checks
    if (FALSE) {
      opar <- par(mfrow = c(4, 1), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      ts.plot(store_para[, "mu"], main = "mu")
      ts.plot(store_para[, "phi"], main = "phi")
      ts.plot(store_para[, "sigma"], main = "sigma")
      abline(h = 0, col = "blue")
      par(opar)
    }
    if (FALSE) {
      opar <- par(mfrow = c(1, 3), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      qqnorm(sample(c(-1, 1), draws, replace = TRUE) * store_para[, "sigma"] * sqrt(2 * priorspec$sigma2$rate)); abline(0, 1, col = "blue")
      qqnorm((store_para[, "mu"] - priorspec$mu$mean) / priorspec$mu$sd); abline(0, 1, col = "blue")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2))); abline(0, 1, col = "blue")
      par(opar)
    }
  }
})

test_that("general SV passes Geweke test", {
  skip_on_cran()
  skip_on_ci()

  set.seed(60)
  result <- .Call(`_stochvol_geweke_general_cpp`, PACKAGE = "stochvol")
  draws <- result$draws / result$thin

  for (centered in c(FALSE, TRUE)) {
    priorspec <-
      if (centered) {  # pro-centered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.9),
                       phi = sv_beta(shape1 = 2, shape2 = 1.5),
                       sigma2 = sv_gamma(shape = 0.9, rate = 0.9),
                       nu = sv_exponential(rate = 0.1),
                       rho = sv_beta(shape1 = 5, shape2 = 5))
      } else {  # pro-noncentered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                       phi = sv_beta(shape1 = 5, shape2 = 1.5),
                       sigma2 = sv_gamma(shape = 0.9, rate = 9),
                       nu = sv_exponential(rate = 0.1),
                       rho = sv_beta(shape1 = 5, shape2 = 5))
      }
    centered_text <- if (centered) "centered" else "noncentered"
    store_para <- t(result[[centered_text]]$para)
    colnames(store_para) <- c("mu", "phi", "sigma", "rho", "nu")


    thin_skip <- ceiling(draws / min(c(4500, apply(store_para, 2, coda::effectiveSize))))
    thin_index <- seq(1, draws, by = thin_skip)

    expect_gt(length(thin_index), 250)

    expect_gt(shapiro.test(qnorm(pgamma(store_para[thin_index, "sigma"]^2, priorspec$sigma2$shape, rate = priorspec$sigma2$rate)))$p.value, 1e-2)
    expect_gt(shapiro.test(qnorm(pnorm(store_para[thin_index, "mu"], priorspec$mu$mean, priorspec$mu$sd)))$p.value, 1e-2)
    expect_gt(shapiro.test(qnorm(pbeta(0.5 * (1 + store_para[thin_index, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)))$p.value, 1e-2)
    expect_gt(shapiro.test(qnorm(pbeta(0.5 * (1 + store_para[thin_index, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2)))$p.value, 1e-2)
    expect_gt(shapiro.test(qnorm(pexp(store_para[thin_index, "nu"] - 2, rate = priorspec$nu$rate)))$p.value, 1e-2)

    # visual tests for manual checks
    if (FALSE) {
      opar <- par(mfrow = c(5, 1), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      ts.plot(store_para[, "mu"], main = "mu")
      ts.plot(store_para[, "phi"], main = "phi")
      ts.plot(store_para[, "sigma"], main = "sigma")
      ts.plot(store_para[, "nu"], main = "nu")
      ts.plot(store_para[, "rho"], main = "rho")
      par(opar)
    }
    if (FALSE) {
      idx <- seq(1, draws)
      opar <- par(mfrow = c(3, 2))
      qqnorm(qnorm(pgamma(store_para[idx, "sigma"]^2, priorspec$sigma2$shape, rate = priorspec$sigma2$rate)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~sigma)); abline(0, 1, col = "red")
      qqnorm(qnorm(pnorm(store_para[idx, "mu"], priorspec$mu$mean, priorspec$mu$sd)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~mu)); abline(0, 1, col = "red")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~phi)); abline(0, 1, col = "red")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~rho)); abline(0, 1, col = "red")
      qqnorm(qnorm(pexp(store_para[idx, "nu"] - 2, rate = priorspec$nu$rate)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~nu)); abline(0, 1, col = "red")
      par(opar)
    }
  }
})

test_that("default fast SV is efficient", {
  skip_on_cran()
  skip_on_ci()

  set.seed(61)
  n <- 150L
  # Pro-centered
  cat("Centered\n")
  priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.9),
                     phi = sv_beta(shape1 = 2, shape2 = 1.5),
                     sigma2 = sv_gamma(shape = 0.9, rate = 0.4))
  ## Simulate data
  sim <- svsim(n, mu = -9, phi = 0.45, sigma = 2)
  samp <- svsample(sim$y, draws = 10000, burnin = 1000, priorspec = priorspec,
                   startpara = list(mu = -9, phi = 0.45, sigma = 2))
  eff_size <- coda::effectiveSize(para(samp)[, sampled_parameters(samp)])
  geweke_test <- 0.5 - abs(0.5 - pnorm(coda::geweke.diag(para(samp)[, sampled_parameters(samp)])$z))

  cat("Effective size\n")
  print(eff_size)

  expect_gt(min(geweke_test), 0.01)
  expect_gt(min(eff_size), 300)
  # Pro-noncentered
  cat("Non-centered\n")
  priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                     phi = sv_beta(shape1 = 20, shape2 = 1.5),
                     sigma2 = sv_gamma(shape = 0.9, rate = 0.9))
  ## Simulate data
  sim <- svsim(n, mu = -9, phi = 0.99, sigma = 1.5)
  samp <- svsample(sim$y, draws = 30000, burnin = 1000, priorspec = priorspec,
                   startpara = list(mu = -9, phi = 0.99, sigma = 1.5))
  eff_size <- coda::effectiveSize(para(samp)[, sampled_parameters(samp)])
  geweke_test <- 0.5 - abs(0.5 - pnorm(coda::geweke.diag(para(samp)[, sampled_parameters(samp)])$z))

  cat("Effective size\n")
  print(eff_size)

  expect_gt(min(geweke_test), 0.01)
  expect_gt(min(eff_size), 300)
})

test_that("default general SV is efficient", {
  skip_on_cran()
  skip_on_ci()

  set.seed(61)
  n <- 150L
  # Pro-centered
  cat("Centered\n")
  priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.9),
                     phi = sv_beta(shape1 = 2, shape2 = 1.5),
                     sigma2 = sv_gamma(shape = 0.9, rate = 0.4),
                     nu = sv_exponential(rate = 0.5),
                     rho = sv_beta(shape1 = 5, shape2 = 5))
  ## Simulate data
  sim <- svsim(n, mu = -9, phi = 0.45, sigma = 2, nu = 5, rho = -0.6)
  samp <- svsample(sim$y, draws = 30000, burnin = 1000, priorspec = priorspec,
                   startpara = list(mu = -9, phi = 0.45, sigma = 2, nu = 12, rho = -0.6))
  eff_size <- coda::effectiveSize(para(samp)[, sampled_parameters(samp)])
  geweke_test <- 0.5 - abs(0.5 - pnorm(coda::geweke.diag(para(samp)[, sampled_parameters(samp)])$z))

  cat("Effective size\n")
  print(eff_size)

  expect_gt(min(geweke_test), 0.01)
  expect_gt(min(eff_size), 120)
  # Pro-noncentered
  cat("Non-centered\n")
  priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                     phi = sv_beta(shape1 = 20, shape2 = 1.5),
                     sigma2 = sv_gamma(shape = 0.9, rate = 0.9),
                     nu = sv_exponential(rate = 0.1),
                     rho = sv_beta(shape1 = 5, shape2 = 5))
  ## Simulate data
  sim <- svsim(n, mu = -9, phi = 0.95, sigma = 1.5, nu = 12, rho = -0.6)
  samp <- svsample(sim$y, draws = 30000, burnin = 1000, priorspec = priorspec,
                   startpara = list(mu = -9, phi = 0.95, sigma = 1.5, nu = 12, rho = -0.6))
  eff_size <- coda::effectiveSize(para(samp)[, sampled_parameters(samp)])
  geweke_test <- 0.5 - abs(0.5 - pnorm(coda::geweke.diag(para(samp)[, sampled_parameters(samp)])$z))

  cat("Effective size\n")
  print(eff_size)

  expect_gt(min(geweke_test), 0.01)
  expect_gt(min(eff_size), 150)
})
