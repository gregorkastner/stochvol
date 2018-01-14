# run this file using
# > test_file('test_posterior_distribution.R')

testthat::context('Posterior distributions')

testthat::test_that('fnMixtureStatePostDist results in the correct format', {
  iN <- 30
  dMu <- 1
  dSigma <- 1
  dRho <- -0.2
  dPhi <- 0.4
  zooDat <- fnGenLogSV(iN, dRho, dSigma, dPhi, dMu, iSeed=1234)
  zooReformDat <- fnReformulateLognSV(zooDat)
  vEpsStar <- as.numeric(zooReformDat[,'eps*'])
  vEta <- as.numeric(zooReformDat[,'eta'][-iN])
  vD <- as.numeric(zooReformDat[,'d'])
  
  mPostDist <- fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, dSigma, dRho)
  # test class
  testthat::expect_is(mPostDist, "matrix")
  # test dimensions
  testthat::expect_equivalent(ncol(mPostDist), 10)
  testthat::expect_equivalent(nrow(mPostDist), iN)
  # test probability
  testthat::expect_true(all(mPostDist >= 0))
  testthat::expect_equivalent(rowSums(mPostDist), rep(1, iN))
})

testthat::test_that('fnMixtureStatePostDist validates the input correctly', {
  iN <- 30
  dMu <- 1
  dSigma <- 1
  dRho <- -0.2
  dPhi <- 0.4
  zooDat <- fnGenLogSV(iN, dRho, dSigma, dPhi, dMu, iSeed=42)
  zooReformDat <- fnReformulateLognSV(zooDat)
  vEpsStar <- as.numeric(zooReformDat[,'eps*'])
  vEta <- as.numeric(zooReformDat[,'eta'][-iN])
  vD <- as.numeric(zooReformDat[,'d'])
  
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, dSigma, dRho), NA)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, -1, dRho), "sigma", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar[-iN], vEta, vD[-iN], dMu, dSigma, dRho), "vEta should have one less element than vEpsStar!", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar[-iN], vEta, vD, dMu, dSigma, dRho), "length", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD[-iN], dMu, dSigma, dRho), "length", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD, c(1, dMu), dSigma, dRho), "single numbers", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, c(1, dSigma), dRho), "single numbers", ignore.case = T)
  testthat::expect_error(fnMixtureStatePostDist(vEpsStar, vEta, vD, dMu, dSigma, c(1, dRho)), "single numbers", ignore.case = T)
})
