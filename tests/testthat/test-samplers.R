context("Samplers are correct")

test_that("vanilla SV passes Geweke test", {
  skip_on_cran()
  skip_if_not_installed(pkg = "magrittr", minimum_version = "1.5")
  library("magrittr")

  # general helper functions
  increment_fun <- function (mu, phi, sigma) {
    function (h_t, eps) { mu + phi * (h_t - mu) + sigma * eps }
  }
  generate_h <- function (len, increment_fun, h0) {
    increment_fun %>%
    Reduce(rnorm(len), h0,
           accumulate = TRUE) %>%
    tail(-1)  # remove h0 from h
  }
  omori_constants <- get_omori_constants()
  p <- omori_constants$prob
  m <- omori_constants$mean
  v <- sqrt(omori_constants$var)

  for (centered in c(FALSE, TRUE)) {
    priorspec <-
      if (centered) {  # pro-centered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                       phi = sv_beta(shape1 = 10, shape2 = 3),
                       sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / 0.1))
      } else {  # pro-noncentered
        specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                       phi = sv_beta(shape1 = 10, shape2 = 1),
                       sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / 0.01))
      }

    set.seed(60)
    startpara <- list(mu = mean(priorspec$mu),
                      phi = -1 + 2*mean(priorspec$phi),
                      sigma = sqrt(mean(priorspec$sigma2)),
                      rho = mean(priorspec$rho),
                      nu = mean(priorspec$nu),
                      beta = mean(priorspec$beta),
                      latent0 = mean(priorspec$mu))

    len <- 30L
    h0 <- startpara$latent0
    #data <- svsim(len, mu = startpara$mu, phi = startpara$phi, sigma = startpara$sigma)
    #h <- 2 * log(data$vol)
    #y <- data$y
    h <- generate_h(len, increment_fun(startpara$mu, startpara$phi, startpara$sigma), h0)
    startlatent <- h
    designmatrix <- matrix(NA_real_, 1, 1)
    fast_sv <- default_fast_sv
    fast_sv$store_indicators <- TRUE
    fast_sv$baseline_parameterization <- if (centered) "centered" else "noncentered"
    fast_sv$init_indicators <- sample.int(10, len, replace = TRUE, prob = p)

    draws <- 20000L
    store_para <- matrix(NA_real_, draws, 4, dimnames = list(NULL, c("mu", "phi", "sigma", "h0")))
    store_y <- matrix(NA_real_, draws, len)
    store_h <- matrix(NA_real_, draws, len)
    store_r <- matrix(NA_real_, draws, len)

    for (tt in seq_len(draws)) {
      z <- fast_sv$init_indicators
      y <- sample(c(-1, 1), len, replace = TRUE) * exp(0.5 * rnorm(len, m[z], v[z])) * exp(0.5 * startlatent)
      res <- svsample_fast_cpp(y, 1L, 30L, designmatrix, priorspec,
                               1L, 1L, "all",
                               startpara, startlatent, keeptau = FALSE,
                               print_settings = list(quiet = TRUE, chain = 1, n_chains = 1),
                               correct_model_misspecification = FALSE,
                               interweave = FALSE,
                               myoffset = 0,
                               fast_sv = fast_sv)
      param <- para(res)[1, ]
      startpara$mu <- param["mu"]
      startpara$phi <- param["phi"]
      startpara$sigma <- param["sigma"]
      if (TRUE) {  # correct way
        startpara$latent0 <- latent0(res)[1]
        fast_sv$init_indicators <- res$indicators[1, ]
        startlatent <- latent(res)[1, ]
      } else {  # wrong way
        startpara$latent0 <- rnorm(1, startpara$mu, startpara$sigma / sqrt(1 - startpara$phi^2))
        fast_sv$init_indicators <- sample.int(10, len, replace = TRUE, prob = p)
        startlatent <- generate_h(len, increment_fun(startpara$mu, startpara$phi, startpara$sigma), startpara$latent0)
      }
      store_para[tt, c("mu", "phi", "sigma")] <- para(res)[1, c("mu", "phi", "sigma")]
      store_para[tt, 4] <- latent0(res)[1]
      store_y[tt, ] <- y
      store_h[tt, ] <- latent(res)[, 1]
      store_r[tt, ] <- res$indicators[, 1]
    }

    thin_skip <- ceiling(draws / min(c(4500, apply(store_para, 2, coda::effectiveSize))))
    thin_index <- seq(1, draws, by = thin_skip)

    expect_gt(shapiro.test((sample(c(-1, 1), draws, replace = TRUE) * store_para[, "sigma"] * sqrt(2 * priorspec$sigma2$rate))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test(((store_para[, "mu"] - priorspec$mu$mean) / priorspec$mu$sd)[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)))[thin_index])$p.value, 1e-5)

    # visual tests for manual checks
    if (FALSE) {
      opar <- par(mfrow = c(5, 1), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      ts.plot(store_para[, "h0"], main = "h0")
      ts.plot(store_para[, "mu"], main = "mu")
      ts.plot(store_para[, "phi"], main = "phi")
      ts.plot(store_para[, "sigma"], main = "sigma")
      #ts.plot(store_h[, len], main = "h_last")
      ts.plot(diff(store_para[, "h0"]), main = "diff(h0)")
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

  for (centered in c(FALSE, TRUE)) {
    priorspec <-
      specify_priors(mu = sv_normal(mean = -9, sd = 0.1),
                     phi = sv_beta(shape1 = 10, shape2 = 5),
                     sigma2 = sv_gamma(shape = 0.5, rate = 0.5 / 1),
                     rho = sv_beta(shape1 = 10, shape2 = 10))
    set.seed(61)
    startpara <- list(mu = mean(priorspec$mu),
                      phi = -1 + 2*mean(priorspec$phi),
                      sigma = sqrt(mean(priorspec$sigma2)),
                      rho = -1 + 2*mean(priorspec$rho),
                      nu = mean(priorspec$nu),
                      beta = mean(priorspec$beta),
                      latent0 = mean(priorspec$mu))

    len <- 30L
    designmatrix <- matrix(NA_real_, 1, 1)
    general_sv <- default_general_sv
    general_sv$multi_asis <- 3
    general_sv$starting_parameterization <- if (centered) "centered" else "noncentered"

    # pre-run to get a good proposal
    print_settings <- list(quiet = TRUE, chain = 1, n_chains = 1)
    data <- svsim(len, mu = startpara$mu, phi = startpara$phi, sigma = startpara$sigma, rho = startpara$rho)
    startlatent <- 2 * log(data$vol)
    y <- data$y
    res <- svsample_general_cpp(y, 60000L, 0L, designmatrix, priorspec,
                                1L, 1L, "all",
                                startpara, startlatent, keeptau = FALSE,
                                print_settings = print_settings,
                                correct_model_misspecification = FALSE,
                                interweave = TRUE,
                                myoffset = 0, general_sv = general_sv)

    expect_gt(tail(res$adaptation[[general_sv$starting_parameterization]]$history[, "Acceptance Rate"], 1), 0.05)

    general_sv$proposal_diffusion_ken <-
      res$adaptation[[general_sv$starting_parameterization]][c("scale", "covariance")]
    general_sv$multi_asis <- 2

    draws <- 20000L
    store_para <- matrix(NA_real_, draws, 5, dimnames = list(NULL, c("mu", "phi", "sigma", "rho", "h0")))
    store_y <- matrix(NA_real_, draws, len)
    store_h <- matrix(NA_real_, draws, len)

    for (tt in seq_len(draws)) {
      for (ttt in seq_len(100)) {
        y[seq_len(len-1)] <- exp(0.5 * head(startlatent, -1)) * (startpara$rho * (tail(startlatent, -1) - startpara$mu - startpara$phi * (head(startlatent, -1) - startpara$mu)) / startpara$sigma + sqrt(1 - startpara$rho^2) * rnorm(len - 1))
        y[len] <- exp(0.5 * tail(startlatent, 1)) * rnorm(1)
        res <- svsample_general_cpp(y, 1L, 0L, designmatrix, priorspec,
                                    1L, 1L, "all",
                                    startpara, startlatent,
                                    keeptau = FALSE,
                                    print_settings = print_settings,
                                    correct_model_misspecification = TRUE,
                                    interweave = FALSE,
                                    myoffset = 0,
                                    general_sv = general_sv)
        param <- para(res)[1, ]
        startpara$mu <- param["mu"]
        startpara$phi <- param["phi"]
        startpara$sigma <- param["sigma"]
        startpara$rho <- param["rho"]
        startpara$latent0 <- latent0(res)[1]
        startlatent <- latent(res)[1, ]
      }

      store_para[tt, c("mu", "phi", "sigma", "rho")] <- para(res)[1, c("mu", "phi", "sigma", "rho")]
      store_para[tt, "h0"] <- latent0(res)[1]
      store_y[tt, ] <- y
      store_h[tt, ] <- latent(res)[1, ]
    }

    thin_skip <- ceiling(draws / min(c(4500, apply(store_para, 2, coda::effectiveSize))))
    thin_index <- seq(1, draws, by = thin_skip)

    expect_gt(shapiro.test((sample(c(-1, 1), draws, replace = TRUE) * store_para[, "sigma"] * sqrt(2 * priorspec$sigma2$rate))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test(((store_para[, "mu"] - priorspec$mu$mean) / priorspec$mu$sd)[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2)))[thin_index])$p.value, 1e-5)
    expect_gt(shapiro.test((qnorm(pbeta(0.5 * (1 + store_para[, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2)))[thin_index])$p.value, 1e-5)

    # visual tests for manual checks
    if (FALSE) {
      opar <- par(mfrow = c(7, 1), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      ts.plot(store_para[, "h0"], main = "h0")
      ts.plot(store_para[, "mu"], main = "mu")
      ts.plot(store_para[, "phi"], main = "phi")
      ts.plot(store_para[, "sigma"], main = "sigma")
      ts.plot(store_para[, "rho"], main = "rho")
      ts.plot(store_h[, len], main = "h_last")
      ts.plot(store_y[, len], main = "y_last")
      par(opar)
    }
    if (FALSE) {
      idx <- seq(1, draws)
      opar <- par(mfrow = c(2, 2), mgp = c(1.6, 0.6, 0), mar = c(1.5, 1.5, 2, 0.5))
      qqnorm(sample(c(-1, 1), length(idx), replace = TRUE) * store_para[idx, "sigma"] * sqrt(2 * priorspec$sigma2$rate)); abline(0, 1, col = "blue")
      qqnorm((store_para[idx, "mu"] - priorspec$mu$mean) / priorspec$mu$sd); abline(0, 1, col = "blue")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "phi"]), priorspec$phi$shape1, priorspec$phi$shape2))); abline(0, 1, col = "blue")
      qqnorm(qnorm(pbeta(0.5 * (1 + store_para[idx, "rho"]), priorspec$rho$shape1, priorspec$rho$shape2))); abline(0, 1, col = "red")
      par(opar)
    }
  }
})

