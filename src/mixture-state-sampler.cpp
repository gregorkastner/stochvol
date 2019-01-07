#include "auxmix.h"
#include "mixture-state-sampler.h"
#include <Rcpp.h>
#include "parameterization.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

using namespace Rcpp;

NumericMatrix mixture_state_post_dist(const NumericVector eps_star, const NumericVector eta,
                                      const NumericVector d,
                                      const double mu, const double sigma2, const double rho,
                                      const Parameterization centering) {
  
  const int n = eps_star.size();
  const int mix_count = sizeof(mix_prob)/sizeof(mix_prob[0]);
  const double sigma2_used = centering == Parameterization::CENTERED ? sigma2 : 1.0;
  NumericMatrix result(n, mix_count);
  
  for (int r = 0; r < n; r++) {
    for (int c = 0; c < mix_count; c++) {
      const double a = mix_a[c];
      const double b = mix_b[c];
      const double m = mix_mean[c];
      const double v = mix_var[c];
      const double log_prior = log(mix_prob[c]);
      
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
                               const NumericVector ht,
                               const double phi, const double rho,
                               const double sigma2, const double mu,
                               const Parameterization centering) {

  const int n = y_star.size();
  NumericVector eps_star;
  NumericVector eta;
  NumericMatrix post_dist;
  NumericVector unif_vec;
  NumericVector new_states(n);
  
  switch (centering) {
    case Parameterization::CENTERED:
    eps_star = y_star - h;
    eta = (tail(h, -1) - mu) - phi*(head(h, -1) - mu);
    break;
    case Parameterization::NONCENTERED:
    eps_star = y_star - mu - sqrt(sigma2)*ht;
    eta = tail(ht, -1) - phi*head(ht, -1);
    break;
  }
  post_dist = mixture_state_post_dist(eps_star, eta, d, mu, sigma2, rho, centering);
  /*  Exact translation  TODO comment out! */
  for (int r = 0; r < n; r++) {
    new_states(r) = sample(10, 1, true, Nullable<NumericVector>(post_dist.row(r)), false)(0);
  }
  
  return new_states;

  //const int mix_count = post_dist.row(0).size();
  //for (int r = 0; r < n; r++) {
  //  double s = 0;
  //  for (int c = 0; c < mix_count; c++) {
  //    s += post_dist(r, c);
  //    post_dist(r, c) = s;
  //  }
  //}
  //
  //unif_vec = runif(n);
  //for (int r = 0; r < n; r++) {
  //  new_states(r) = sum(post_dist.row(r) < unif_vec(r));  // C++ indexing
  //}
  //
  //return new_states;
}
