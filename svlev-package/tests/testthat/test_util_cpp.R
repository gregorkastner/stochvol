# run this file using
# > test_file('test_util_cpp.R')

testthat::context('C++ utilities')

testthat::test_that('fnAugKalmanFilterCpp is the same as fnAugKalmanFilter', {
  lDat <- fnGenApproxLogSV(20)
  
  va <- dfModelConstants$a[lDat$vS]
  vb <- dfModelConstants$b[lDat$vS]
  vm <- dfModelConstants$m[lDat$vS]
  vv <- dfModelConstants$v[lDat$vS]
  
  lRes <- fnAugKalmanFilter(0.77, 0.33, -0.9, va, vb, vm, vv,
                            lDat$vD, lDat$vYStar, -1, 40, "centered")
  
  lResCpp <- fnAugKalmanFilterCpp(0.77, 0.33, -0.9, va, vb, vm, vv,
                                  lDat$vD, lDat$vYStar, -1, 40, "centered")
  
  testthat::expect_equal(names(lResCpp), names(lRes))
  testthat::expect_equal(lResCpp$D, lRes$D)
  
  lRes <- fnAugKalmanFilter(0.77, 0.33, -0.9, va, vb, vm, vv,
                            lDat$vD, lDat$vYStar, -1, 40, "non-centered")
  
  lResCpp <- fnAugKalmanFilterCpp(0.77, 0.33, -0.9, va, vb, vm, vv,
                                  lDat$vD, lDat$vYStar, -1, 40, "non-centered")
  
  testthat::expect_equal(names(lResCpp), names(lRes))
  testthat::expect_equal(lResCpp$D, lRes$D)
})
