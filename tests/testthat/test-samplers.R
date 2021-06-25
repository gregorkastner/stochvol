context("Samplers are correct")

test_that("vanilla SV passes Geweke test", {
  skip_on_cran()

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
            skip("NYI")
  skip_on_cran()

  for (centered in c(FALSE, TRUE)) {
    priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                     phi = sv_beta(shape1 = 10, shape2 = 5),
                     sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / 1),
                     rho = sv_beta(shape1 = 10, shape2 = 10),
                     nu = sv_exponential(rate = 0.05))
    set.seed(61)
    startpara <- list(mu = mean(priorspec$mu),
                      sigma = sqrt(mean(priorspec$sigma2)),
                      phi = -1 + 2*mean(priorspec$phi),
                      rho = -1 + 2*mean(priorspec$rho),
                      nu = mean(priorspec$nu),
                      #phi = mean(priorspec$phi),
                      #rho = mean(priorspec$rho),
                      #nu = mean(priorspec$nu)+2,
                      beta = mean(priorspec$beta),
                      latent0 = mean(priorspec$mu))

    len <- 30L
    designmatrix <- matrix(NA_real_, 1, 1)
    general_sv <- get_default_general_sv(priorspec)
    general_sv$multi_asis <- 3
    general_sv$starting_parameterization <- if (centered) "centered" else "noncentered"
    general_sv$adaptation_object[[general_sv$starting_parameterization]]$memory <- matrix(NA, 300L, 3L)

    # pre-run to get an efficient proposal for better mixing
    print_settings <- list(quiet = TRUE, chain = 1, n_chains = 1)
    data <- svsim(len, mu = startpara$mu, phi = startpara$phi, sigma = startpara$sigma, rho = startpara$rho, nu = startpara$nu)
    startlatent <- data$latent
    tau <- data$tau
    y <- data$y
    res <- svsample_general_cpp(y, 60000L, 0L, designmatrix, priorspec,
                                1L, 1L, "all",
                                startpara, startlatent, keeptau = FALSE,
                                print_settings = print_settings,
                                correct_model_misspecification = FALSE,
                                interweave = TRUE,
                                myoffset = 0, general_sv = general_sv)

    adaptation_object <- res$general_sv$adaptation_object[[general_sv$starting_parameterization]]
    acceptance_rates <- adaptation_object$memory[, "Acceptance Rate"]
    acceptance_rates <- acceptance_rates[!is.na(acceptance_rates)]
    expect_gt(tail(acceptance_rates, 1), 0.02)

    general_sv$adaptation_object[[general_sv$starting_parameterization]]$memory <- matrix(NA, 0L, 3L)

    general_sv$proposal_diffusion_ken <-
      list(scale = adaptation_object$cached_scale,
           covariance = adaptation_object$cached_covariance)
    general_sv$multi_asis <- 2

    draws <- 20000L
    store_para <- matrix(NA_real_, draws, 6, dimnames = list(NULL, c("mu", "phi", "sigma", "nu", "rho", "h0")))
    store_y <- matrix(NA_real_, draws, len)
    store_h <- matrix(NA_real_, draws, len)
    store_tau <- matrix(NA_real_, draws, len)

    for (tt in seq_len(draws)) {
      for (ttt in seq_len(100)) {
        y[seq_len(len-1)] <- exp(0.5 * head(startlatent, -1)) * (startpara$rho * (tail(startlatent, -1) - startpara$mu - startpara$phi * (head(startlatent, -1) - startpara$mu)) / startpara$sigma + sqrt(1 - startpara$rho^2) * rnorm(len - 1))
        y[len] <- exp(0.5 * tail(startlatent, 1)) * rnorm(1)
        y <- y * sqrt(tau)
        res <- svsample_general_cpp(y, 1L, 0L, designmatrix, priorspec,
                                    1L, 1L, "all",
                                    startpara, startlatent,
                                    keeptau = TRUE,
                                    print_settings = print_settings,
                                    correct_model_misspecification = TRUE,
                                    interweave = FALSE,
                                    myoffset = 0,
                                    general_sv = general_sv)
        param <- res$para[1, ]
        startpara$mu <- param["mu"]
        startpara$phi <- param["phi"]
        startpara$sigma <- param["sigma"]
        startpara$nu <- param["nu"]
        startpara$rho <- param["rho"]
        startpara$latent0 <- res$latent0[1]
        startlatent <- res$latent[1, ]
        tau <- res$tau[1, ]
        general_sv$init_tau <- res$tau[1, ]
      }

      store_para[tt, c("mu", "phi", "sigma", "nu", "rho")] <- res$para[1, c("mu", "phi", "sigma", "nu", "rho")]
      store_para[tt, "h0"] <- res$latent0[1]
      store_y[tt, ] <- y
      store_h[tt, ] <- res$latent[1, ]
      store_tau[tt, ] <- res$tau[1, ]
    }

    thin_skip <- ceiling(draws / min(c(4500, apply(store_para, 2, coda::effectiveSize))))
    thin_index <- seq(1, draws, by = thin_skip)

    expect_gt(shapiro.test((sample(c(-1, 1), draws, replace = TRUE) * store_para[, "sigma"] * sqrt(2 * priorspec$sigma2$rate))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test(((store_para[, "mu"] - priorspec$mu$mean) / priorspec$mu$sd)[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2)))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pexp(store_para[, "nu"] - 2, rate = priorspec$nu$rate)))[thin_index])$p.value, 1e-5)

    # visual tests for manual checks
    if (FALSE) {
      opar <- par(mfrow = c(8, 1), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      ts.plot(store_para[, "h0"], main = "h0")
      ts.plot(store_para[, "mu"], main = "mu")
      ts.plot(store_para[, "phi"], main = "phi")
      ts.plot(store_para[, "sigma"], main = "sigma")
      ts.plot(store_para[, "nu"], main = "nu")
      ts.plot(store_para[, "rho"], main = "rho")
      ts.plot(store_h[, len], main = "h_last")
      ts.plot(store_y[, len], main = "y_last")
      par(opar)
    }
    if (FALSE) {
      idx <- seq(1, draws)
      opar <- par(mfrow = c(3, 2)) #, mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      qqnorm(sample(c(-1, 1), length(idx), replace = TRUE) * store_para[idx, "sigma"] * sqrt(2 * priorspec$sigma2$rate), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~sigma)); abline(0, 1, col = "red")
      qqnorm((store_para[idx, "mu"] - priorspec$mu$mean) / priorspec$mu$sd, ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~mu)); abline(0, 1, col = "red")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~phi)); abline(0, 1, col = "red")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~rho)); abline(0, 1, col = "red")
      qqnorm(qnorm(pexp(store_para[idx, "nu"] - 2, rate = priorspec$nu$rate)), ylab = "Theoretical Quantiles", xlab = "Sample Quantiles", main = bquote("Geweke QQ-Plot for Parameter"~nu)); abline(0, 1, col = "red")
      #qqplot(store_tau[idx, 9], 1/rgamma(2*length(idx), shape = .5*store_para[idx, "nu"], rate = .5*(store_para[idx, "nu"]-2))); abline(0, 1, col = "red")  # keep in mind that this is very heavy-tailed
      par(opar)
    }
  }
})

