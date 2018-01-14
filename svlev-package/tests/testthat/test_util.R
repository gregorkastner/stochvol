# run this file using
# > test_file('test_util.R')

testthat::context('Utilities')

testthat::test_that('fnAugKalmanFilter returns with the correct format', {
  set.seed(11)
  iT <- 20
  lDat <- fnGenApproxLogSV(iT)
  
  va <- dfModelConstants$a[lDat$vS]
  vb <- dfModelConstants$b[lDat$vS]
  vm <- dfModelConstants$m[lDat$vS]
  vv <- dfModelConstants$v[lDat$vS]
  
  lRes <- fnAugKalmanFilter(0.77, 0.33, -0.9, va, vb, vm, vv, lDat$vD, lDat$vYStar, -1, 40, "centered")
  testthat::expect_equal(names(lRes), c("sigma2", "D", "J1", "L", "f", "F", "hts", "v", "Q", "q", "jt22", "h1var"))
  testthat::expect_length(lRes$D, iT)
  testthat::expect_length(lRes$J1, iT)
  testthat::expect_length(lRes$L, iT)
  testthat::expect_length(lRes$f, iT)
  testthat::expect_length(lRes[["F"]], iT)
  testthat::expect_length(lRes$hts, iT)
  testthat::expect_length(lRes$v, iT)
  testthat::expect_length(lRes$Q, 1)
  testthat::expect_is(lRes$D, "numeric")
  testthat::expect_is(lRes$q, "numeric")
  testthat::expect_is(lRes$Q, "numeric")
  testthat::expect_is(lRes$jt22, "numeric")
  testthat::expect_is(lRes$h1var, "numeric")
  
  lRes <- fnAugKalmanFilter(0.77, 0.33, -0.9, va, vb, vm, vv, lDat$vD, lDat$vYStar, -1, 40, "non-centered")
  testthat::expect_equal(names(lRes), c("sigma2", "D", "J1", "L", "f", "F", "hts", "v", "Q", "q", "jt22", "h1var"))
  testthat::expect_length(lRes$D, iT)
  testthat::expect_length(lRes$J1, iT)
  testthat::expect_length(lRes$L, iT)
  testthat::expect_length(lRes$f, iT)
  testthat::expect_length(lRes[["F"]], iT)
  testthat::expect_length(lRes$hts, iT)
  testthat::expect_length(lRes$v, iT)
  testthat::expect_length(lRes$Q, 1)
  testthat::expect_is(lRes$D, "numeric")
  testthat::expect_is(lRes$q, "numeric")
  testthat::expect_is(lRes$Q, "numeric")
  testthat::expect_is(lRes$jt22, "numeric")
  testthat::expect_is(lRes$h1var, "numeric")
})
