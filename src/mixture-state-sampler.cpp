#include <Rcpp.h>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include "mixture-state-sampler.h"
using namespace Rcpp;

NumericMatrix mixture_state_post_dist(const NumericVector eps_star, const NumericVector eta,
                                      const NumericVector d,
                                      const double mu, const double sigma2, const double rho,
                                      const CharacterVector centering,
                                      const DataFrame mixing_constants) {
  
  const NumericVector mixing_a = mixing_constants["a"];
  const NumericVector mixing_b = mixing_constants["b"];
  const NumericVector mixing_m = mixing_constants["m"];
  const NumericVector mixing_p = mixing_constants["p"];
  const NumericVector mixing_v = mixing_constants["v"];
  
  std::string scentering = as<std::string>(centering);
  const int n = eps_star.size();
  const int mix_count = mixing_p.size();
  const double sigma2_used = scentering == "centered" ? sigma2 : 1.0;
  Environment globalenv = Environment::global_env();
  //NumericMatrix result = globalenv["mspd"];
  NumericMatrix result(n, mix_count);
  
  for (int r = 0; r < n; r++) {
    for (int c = 0; c < mix_count; c++) {
      const double a = mixing_a(c);
      const double b = mixing_b(c);
      const double m = mixing_m(c);
      const double v = mixing_v(c);
      const double log_prior = log(mixing_p(c));
      
      double log_eps_star_lik = -0.5 * pow((eps_star[r] - m)/v, 2) - log(v*sqrt(2*M_PI));
      double log_eta_lik;
      if (r < n-1) {
        log_eta_lik = -0.5/(sigma2_used*(1-rho*rho)) * pow(eta[r] - rho*sqrt(sigma2_used)*d[r]*exp(m/2)*(a+b*(eps_star[r]-m)), 2) - 0.5*log(2*M_PI*(1-rho*rho)*sigma2_used);
      } else {
        log_eta_lik = 0.0;
      }
      /*log_*/result(r, c) = log_prior + log_eps_star_lik + log_eta_lik;
    }
    const double max_log_result_row = Rcpp::max(/*log_*/result.row(r));
    /*log_*/result.row(r) = /*log_*/result.row(r) - (max_log_result_row + log(Rcpp::sum(Rcpp::exp(/*log_*/result.row(r)-max_log_result_row))));
    result.row(r) = Rcpp::exp(/*log_*/result.row(r));
    result.row(r) = result.row(r) / Rcpp::sum(result.row(r));
  }
  
  return result;
}

NumericVector draw_s_auxiliary(const NumericVector y_star,
                               const NumericVector d,
                               const NumericVector h,
                               const double phi, const double rho,
                               const double sigma2, const double mu,
                               const CharacterVector centering,
                               const DataFrame mixing_constants) {
  const int n = y_star.size();
  NumericVector eps_star;
  NumericVector eta;
  NumericMatrix post_dist;
  NumericVector unif_vec;
  NumericVector new_states(n);
  
  std::string scentering = as<std::string>(centering);
  if (scentering == "centered") {
    eps_star = y_star - h;
    eta = (tail(h, -1) - mu) - phi*(head(h, -1) - mu);
  } else if (scentering == "non-centered") {
    eps_star = y_star - mu - sqrt(sigma2)*h;
    eta = tail(h, -1) - phi*head(h, -1);
  } else {
    ::Rf_error("Invalid centering");
  }
  post_dist = mixture_state_post_dist(eps_star, eta, d, mu, sigma2, rho, centering, mixing_constants);
  /*  Exact translation  TODO comment out! */
  for (int r = 0; r < n; r++) {
    new_states(r) = sample(10, 1, true, Nullable<NumericVector>(post_dist.row(r)), false)(0);
  }
  
  return new_states;
  /*/
  const int mix_count = post_dist.row(0).size();
  for (int r = 0; r < n; r++) {
    double s = 0;
    for (int c = 0; c < mix_count; c++) {
      s += post_dist(r, c);
      post_dist(r, c) = s;
    }
  }
  
  unif_vec = runif(n);
  for (int r = 0; r < n; r++) {
    new_states(r) = sum(post_dist.row(r) < unif_vec(r));  // C++ indexing
  }
  
  return new_states; /* */
}
