test_that("zeroes do not interfere with execution", {
  # Create input for the samplers
  y <- c(-0.867537088, 0.424335171, -0.048784030, -0.749483682,
        -0.332278654, 0.503644166, 0.007733118, -0.485776150,
        0.163207939, 0.123654636, 0.535032546, 0.799077152,
        0.306340493, -1.451189223, -1.095406999, 2.765907287,
        -2.196736048, 0.879586318, 0.359074117, -1.700210392)

  draws <- 50
  burnin <- 10

  designmat <- matrix(c(0.18554962, 1.71128430, 2.76005663, 0.27480471, -0.04419060,
                        -1.54191962, 0.03668254, 1.71851215, -1.48924188, 1.45517746,
                        -0.03498587, 0.91532504, 0.61450314, 0.64796734, 0.41568157,
                        -1.79375801, -0.33170476, -0.06896315, -0.83336860, 1.35250078,
                        -2.17504467, 0.19803696, -0.14800168, -0.95317998, 0.92041904,
                        -0.96628316, 0.59918070, 0.81186793, 0.07408481, 0.25670553,
                        -0.56785729, 0.77431129, 1.04092342, -0.52666850, 1.39242409,
                        1.79126519, 0.21054807, 0.19770369, -0.90162245, -0.62361923),
                      length(y), 2)
  ar_values <- c(0, 2)

  designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
  keeptime_values <- c("all", "last")
  thin_values <- c(1, 3)

  pred_steps <- 3
  pred_designmat <- designmat[seq_len(pred_steps), ]

  expect_message(svsample(c(0, y), draws = draws, burnin = burnin, quiet = TRUE), "Argument 'y' *")
})

test_that("svsample executes", {
  # Create input for the samplers
  y <- c(-0.867537088, 0.424335171, -0.048784030, -0.749483682,
        -0.332278654, 0.503644166, 0.007733118, -0.485776150,
        0.163207939, 0.123654636, 0.535032546, 0.799077152,
        0.306340493, -1.451189223, -1.095406999, 2.765907287,
        -2.196736048, 0.879586318, 0.359074117, -1.700210392)

  draws <- 50
  burnin <- 10

  designmat <- matrix(c(0.18554962, 1.71128430, 2.76005663, 0.27480471, -0.04419060,
                        -1.54191962, 0.03668254, 1.71851215, -1.48924188, 1.45517746,
                        -0.03498587, 0.91532504, 0.61450314, 0.64796734, 0.41568157,
                        -1.79375801, -0.33170476, -0.06896315, -0.83336860, 1.35250078,
                        -2.17504467, 0.19803696, -0.14800168, -0.95317998, 0.92041904,
                        -0.96628316, 0.59918070, 0.81186793, 0.07408481, 0.25670553,
                        -0.56785729, 0.77431129, 1.04092342, -0.52666850, 1.39242409,
                        1.79126519, 0.21054807, 0.19770369, -0.90162245, -0.62361923),
                      length(y), 2)
  ar_values <- c(0, 2)

  designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
  keeptime_values <- c("all", "last")
  thin_values <- c(1, 3)

  pred_steps <- 3
  pred_designmat <- designmat[seq_len(pred_steps), ]

  expect_warning(svsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_s3_class("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_s3_class("svdraws")
      }
    }
  }
})

test_that("svsample with nu executes", {
  # Create input for the samplers
  y <- c(-0.867537088, 0.424335171, -0.048784030, -0.749483682,
        -0.332278654, 0.503644166, 0.007733118, -0.485776150,
        0.163207939, 0.123654636, 0.535032546, 0.799077152,
        0.306340493, -1.451189223, -1.095406999, 2.765907287,
        -2.196736048, 0.879586318, 0.359074117, -1.700210392)

  draws <- 50
  burnin <- 10

  designmat <- matrix(c(0.18554962, 1.71128430, 2.76005663, 0.27480471, -0.04419060,
                        -1.54191962, 0.03668254, 1.71851215, -1.48924188, 1.45517746,
                        -0.03498587, 0.91532504, 0.61450314, 0.64796734, 0.41568157,
                        -1.79375801, -0.33170476, -0.06896315, -0.83336860, 1.35250078,
                        -2.17504467, 0.19803696, -0.14800168, -0.95317998, 0.92041904,
                        -0.96628316, 0.59918070, 0.81186793, 0.07408481, 0.25670553,
                        -0.56785729, 0.77431129, 1.04092342, -0.52666850, 1.39242409,
                        1.79126519, 0.21054807, 0.19770369, -0.90162245, -0.62361923),
                      length(y), 2)
  ar_values <- c(0, 2)

  designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
  keeptime_values <- c("all", "last")
  thin_values <- c(1, 3)

  pred_steps <- 3
  pred_designmat <- designmat[seq_len(pred_steps), ]

  expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, quiet = TRUE), NA) %>%
    expect_s3_class("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_s3_class("svdraws")
      }
    }
  }
})

test_that("svsample with rho executes", {
  # Create input for the samplers
  y <- c(-0.867537088, 0.424335171, -0.048784030, -0.749483682,
        -0.332278654, 0.503644166, 0.007733118, -0.485776150,
        0.163207939, 0.123654636, 0.535032546, 0.799077152,
        0.306340493, -1.451189223, -1.095406999, 2.765907287,
        -2.196736048, 0.879586318, 0.359074117, -1.700210392)

  draws <- 50
  burnin <- 10

  designmat <- matrix(c(0.18554962, 1.71128430, 2.76005663, 0.27480471, -0.04419060,
                        -1.54191962, 0.03668254, 1.71851215, -1.48924188, 1.45517746,
                        -0.03498587, 0.91532504, 0.61450314, 0.64796734, 0.41568157,
                        -1.79375801, -0.33170476, -0.06896315, -0.83336860, 1.35250078,
                        -2.17504467, 0.19803696, -0.14800168, -0.95317998, 0.92041904,
                        -0.96628316, 0.59918070, 0.81186793, 0.07408481, 0.25670553,
                        -0.56785729, 0.77431129, 1.04092342, -0.52666850, 1.39242409,
                        1.79126519, 0.21054807, 0.19770369, -0.90162245, -0.62361923),
                      length(y), 2)
  ar_values <- c(0, 2)

  designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
  keeptime_values <- c("all", "last")
  thin_values <- c(1, 3)

  pred_steps <- 3
  pred_designmat <- designmat[seq_len(pred_steps), ]

  expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), quiet = TRUE), NA) %>%
    expect_s3_class("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_s3_class("svdraws")
      }
    }
  }
})

test_that("svsample with constant prior is error-free", {
  sim <- svsim(10, mu = -10, phi = 0.99, sigma = 0.2)
  sv_prior <- specify_priors(mu = sv_constant(-10), phi = sv_constant(0.99), sigma2 = sv_constant(0.2^2))
  expect_no_error(svsample(sim, draws = 50, burnin = 10, priorspec = sv_prior))
})