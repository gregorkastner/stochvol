context("Formula interface")

test_that("formula interface and prediction work", {
  # simulate data
  set.seed(56)
  n <- 50L
  predn <- 10L
  dat <- cbind(x = runif(n, 3, 4),
               z = runif(n, -1, -0.5))
  preddat <- cbind(x = runif(predn, 3, 4),
                   z = runif(predn, -1, -0.5))
  designmatrix <- matrix(c(dat[, "x"], dat[, "x"]^2, log10(dat[, "x"]),
                         dat[, "z"]), ncol = 4)
  predmatrix <- matrix(c(preddat[, "x"], preddat[, "x"]^2, log10(preddat[, "x"]),
                       preddat[, "z"]), ncol = 4)
  betas <- matrix(c(-1, 1, 2, 0), ncol = 1)
  y <- designmatrix %*% betas + svsim(n)$y
  # standard interface
  set.seed(57)
  res1 <- svsample(y, designmatrix = designmatrix, quiet = TRUE)
  pred1 <- suppressWarnings(predict(res1, newdata = predmatrix))
  # formula interface
  data <- as.data.frame(dat)
  data$y <- y
  preddata <- as.data.frame(preddat)
  set.seed(57)
  res2 <- svlm(y ~ 0 + x + I(x^2) + log10(x) + z, data = data, quiet = TRUE)
  pred2 <- suppressWarnings(predict(res2, newdata = preddata))

  expect_equal(para(res1)[, "mu"], para(res2)[, "mu"])
  expect_equal(predy(pred1), predy(pred2))
})
