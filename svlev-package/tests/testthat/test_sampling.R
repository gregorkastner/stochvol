# run this file using
# > test_file('test_sampling.R')

testthat::context("Samplers")

testthat::test_that("fnSampleMixtureState output has correct format", {
  set.seed(1)
  iN <- 30
  dMu <- rnorm(1)
  dSigma <- .1 + rnorm(1)^2
  dRho <- 2*rbeta(1, 2, 1) - 1
  dPhi <- rbeta(1, 2, 10)
  vYStar <- rnorm(iN, sd = 1)
  vH <- rnorm(iN, sd = .1)  # quite sensitive to larger numbers, why??
  vD <- sample(c(-1,1), iN, replace = T)
  
  vMixState <- fnSampleMixtureState(vYStar, vH, vD, dMu, dSigma, dRho, dPhi)
  
  testthat::expect_is(vMixState, "numeric")
  testthat::expect_true(all(vMixState %in% 1:10))
})
