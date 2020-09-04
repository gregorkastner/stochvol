#include <testthat.h>
#include <Rcpp.h>
#include "densities.h"

bool almost_equal (
    const double a,
    const double b,
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon())) {
  using namespace std;
  return fabs(a - b) < min(fabs(a), fabs(b)) * tol;
}

context("Own density implementations") {
  using namespace stochvol;
  using namespace std;

  test_that("logdnorm works as R::dnorm(..., true)") {
    const double val1 = 0.9,
                 val2 = 0.95;
    const double m = 1,
                 s = 2;
    expect_true(almost_equal(logdnorm(val1, m, s) - logdnorm(val2, m, s), R::dnorm(val1, m, s, true) - R::dnorm(val2, m, s, true)));
    expect_true(almost_equal(logdnorm(val1, m, s + 0.1) - logdnorm(val2, m, s), R::dnorm(val1, m, s + 0.1, true) - R::dnorm(val2, m, s, true)));
    expect_true(almost_equal(logdnorm2(val1, m, s + 0.1, log(s + 0.1)) - logdnorm2(val2, m, s, log(s)), R::dnorm(val1, m, s + 0.1, true) - R::dnorm(val2, m, s, true)));
    expect_true(almost_equal(logdnorm2(val1, m, s) - logdnorm2(val2, m, s), R::dnorm(val1, m, s, true) - R::dnorm(val2, m, s, true)));
    expect_false(almost_equal(logdnorm2(val1, m, s + 0.1) - logdnorm2(val2, m, s), R::dnorm(val1, m, s + 0.1, true) - R::dnorm(val2, m, s, true)));
    expect_false(almost_equal(logdnorm(val1, m, s) - logdnorm(val1, m, s), R::dnorm(val1, m, s, true) - R::dnorm(val2, m, s, true)));
  }

  test_that("logdgamma works as R::dgamma(..., true)") {
    const double val1 = 0.9,
                 val2 = 0.95;
    const double shape = 1,
                 rate = 2;
    expect_true(almost_equal(logdgamma(val1, shape, rate) - logdgamma(val2, shape, rate), R::dgamma(val1, shape, 1/rate, true) - R::dgamma(val2, shape, 1/rate, true)));
    expect_false(almost_equal(logdgamma(val1+0.1, shape, rate) - logdgamma(val2, shape, rate), R::dgamma(val1, shape, 1/rate, true) - R::dgamma(val2, shape, 1/rate, true)));
  }

}

