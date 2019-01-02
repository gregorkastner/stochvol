#include <Rcpp.h>
#include "h-utils.h"
#define _USE_MATH_DEFINES
#include <cmath>

using namespace Rcpp;

double h_log_posterior(const NumericVector h,  // centered
                       const NumericVector y,
                       const double phi,
                       const double rho,
                       const double sigma2,
                       const double mu) {
  const double sigma = sqrt(sigma2);
  const double rho_const = sqrt(1-rho*rho);
  const int n = y.size();
  double result = R::dnorm(h(0), mu, sigma/sqrt(1-phi*phi), true);
  for (int t = 0; t < n-1; t++) {
    result += R::dnorm(h(t+1), mu+phi*(h(t)-mu), sigma, true);
    result += R::dnorm(y(t), exp(.5*h(t))*rho*(h(t+1)-mu-phi*(h(t)-mu))/sigma, exp(.5*h(t))*rho_const, true);
  }
  result += R::dnorm(y(n-1), 0, exp(.5*h(n-1)), true);
  return result;
}

double h_aux_log_posterior(const NumericVector h,  // centered
                           const NumericVector y_star,
                           const NumericVector d,
                           const NumericVector mixing_a,
                           const NumericVector mixing_b,
                           const NumericVector mixing_m,
                           const NumericVector mixing_p,
                           const NumericVector mixing_v,
                           const double phi,
                           const double rho,
                           const double sigma2,
                           const double mu) {
  const int n = y_star.size();
  const int mix_count = mixing_p.size();
  const double sigma = sqrt(sigma2);

  double result = R::dnorm(h(0), mu, sigma/sqrt(1-phi*phi), true);  // log p(h_1 | theta)
  for (int t = 0; t < n; t++) {
    double subresult = 0;
    if (t < n-1) {
      for (int j = 0; j < mix_count; j++) {
        const double C = rho*sigma*exp(mixing_m(j)*.5);  // re-used constant
        const double h_mean = mu+phi*(h(t)-mu)+d(t)*C*mixing_a(j);
        const double h_var = pow(C*mixing_v(j)*mixing_b(j), 2) + sigma2*(1-rho*rho);
        const double yh_cov = d(t)*C*mixing_b(j)*pow(mixing_v(j), 2);
        const double y_mean = h(t)+mixing_m(j) + yh_cov/h_var*(h(t+1)-h_mean);
        const double y_var = (1-pow(yh_cov/mixing_v(j), 2)/h_var) * pow(mixing_v(j), 2);
        subresult += R::dnorm(y_star(t), y_mean, sqrt(y_var), false) * R::dnorm(h(t+1), h_mean, sqrt(h_var), false) * mixing_p(j);
      }
    } else {
      for (int j = 0; j < mix_count; j++) {
        subresult += R::dnorm(y_star(t), h(t)+mixing_m(j), mixing_v(j), false) * mixing_p(j);
      }
    }
    result += log(subresult);
  }

  return result;
}
