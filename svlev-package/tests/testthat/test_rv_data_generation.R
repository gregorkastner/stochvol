# run this file using
# > test_file('test_rv_data_generation.R')

testthat::context('Random number generation')

testthat::test_that('fnMixNormal results in the correct format', {
  n=1000
  dat <- fnMixNormal(n)
  
  testthat::expect_is(dat, "numeric")
  testthat::expect_equal(length(dat), n)
})

testthat::test_that('fnDrawDirichlet results in the correct format', {
  mDist <- matrix(c(.25, .25, .25, .25,
                    .01, .9, .02, .07,
                    0, 0, 1, 0),
                  3, 4, byrow = T)
  
  testthat::expect_is(fnDrawDirichlet(1, mDist), "matrix")
  testthat::expect_is(fnDrawDirichlet(1, mDist[1,]), "matrix")
  testthat::expect_is(fnDrawDirichlet(30, mDist), "matrix")
  testthat::expect_equivalent(dim(fnDrawDirichlet(10, mDist)), c(10, 3))
  testthat::expect_equivalent(dim(fnDrawDirichlet(1, mDist)), c(1, 3))
  testthat::expect_equivalent(dim(fnDrawDirichlet(1, mDist[2,])), c(1, 1))
})

testthat::test_that('an example works for fnDrawDirichlet', {
  n=1000
  mDist <- matrix(c(.25, .25, .25, .25,
                    .01, .9, .02, .07,
                    0, 0, 1, 0),
                  3, 4, byrow = T)
  dat <- fnDrawDirichlet(n, mDist)
  
  testthat::expect_equivalent(dat[,3], rep(3, n))
  testthat::expect_equal(mean(dat[,2] == 2), .9, tolerance = .1)
})

testthat::context('Model generation')

testthat::test_that('fnGenLogSV results in the correct format', {
  n=1000
  rho=-0.2
  sigma=0.5
  phi=0.7
  mu=-0.1
  dat <- fnGenLogSV(iN=n, dRho=rho, dSigma=sigma, dMu=mu, dPhi=phi)
  
  testthat::expect_is(dat, "zoo")
  testthat::expect_equal(nrow(dat), n)
  testthat::expect_equal(ncol(dat), 4)
  testthat::expect_equivalent(colnames(dat), c('y', 'h', 'eps', 'eta'))
})

testthat::test_that('fnReformulateLognSV results in the correct format', {
  n=1000
  rho=-0.2
  sigma=0.5
  phi=0.7
  mu=-0.1
  sv_dat <- fnGenLogSV(iN=n, dRho=rho, dSigma=sigma, dMu=mu, dPhi=phi)
  dat <- fnReformulateLognSV(sv_dat)
  
  testthat::expect_is(dat, "zoo")
  testthat::expect_equal(nrow(dat), n)
  testthat::expect_equal(ncol(dat), 5)
  testthat::expect_equivalent(colnames(dat), c('y*', 'h', 'eps*', 'eta', 'd'))
})

testthat::test_that('fnReformulateLognSV does the correct thing', {
  n=1000
  rho=-0.2
  sigma=0.5
  phi=0.7
  mu=-0.1
  thresh=1e-10
  sv_dat <- fnGenLogSV(iN=n, dRho=rho, dSigma=sigma, dMu=mu, dPhi=phi)
  dat <- fnReformulateLognSV(sv_dat, dThresh=thresh)
  
  logeps2 <- log(sv_dat[,'eps']^2)
  
  testthat::expect_equivalent(dat[,'h'], sv_dat[,'h'])
  testthat::expect_equivalent(dat[,'eta'], sv_dat[,'eta'])
  testthat::expect_true(all((logeps2 >= 2*log(thresh)) | (logeps2 == dat[,'eps*'])))
})

testthat::context('Approximate model generation')

testthat::test_that('fnGenApproxLogSV returns the correct format', {
  iN <- 10
  lDat <- fnGenApproxLogSV(iN)
  
  testthat::expect_equivalent(names(lDat), c('dMu', 'dPhi', 'dSigma', 'dRho',
                                                'vD', 'vS', 'vH', 'vEpsStar', 'vEta',
                                                'vYStar'))
  testthat::expect_length(lDat$dMu, 1)
  testthat::expect_length(lDat$dPhi, 1)
  testthat::expect_length(lDat$dSigma, 1)
  testthat::expect_length(lDat$dRho, 1)
  testthat::expect_length(lDat$vD, iN)
  testthat::expect_length(lDat$vS, iN)
  testthat::expect_length(lDat$vH, iN)
  testthat::expect_length(lDat$vEpsStar, iN)
  testthat::expect_length(lDat$vEta, iN)
  testthat::expect_length(lDat$vYStar, iN)
})

testthat::test_that('fnGenApproxLogSV throws errors when needed', {
  testthat::expect_error(fnGenApproxLogSV())
  testthat::expect_error(fnGenApproxLogSV(iN=10, vD=rep(1,10)))
  testthat::expect_error(fnGenApproxLogSV(vD=rep(0.1,10)))
  testthat::expect_error(fnGenApproxLogSV(10), NA)
  testthat::expect_error(fnGenApproxLogSV(iN=10), NA)
  testthat::expect_error(fnGenApproxLogSV(vD=rep(1,10)), NA)
  testthat::expect_error(fnGenApproxLogSV(vD=c(1,-1,1,1,-1)), NA)
})
