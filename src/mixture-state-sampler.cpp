#include <RcppArmadillo.h>
#include "auxmix.h"
#include "mixture-state-sampler.h"
#include "parameterization.h"
#include <cmath>

arma::uvec draw_s_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Parameterization centering) {
  const int n = y_star.size();
  const double sigma2_used = centering == Parameterization::CENTERED ? sigma2 : 1.0;
  static const int mix_count = mix_a.n_elem;
  arma::vec eps_star;
  arma::vec eta;
  arma::vec unif_vec;
  arma::uvec new_states(n);
  
  switch (centering) {
    case Parameterization::CENTERED:
    eps_star = y_star - h;
    eta = (h.tail(n-1) - mu) - phi*(h.head(n-1) - mu);
    break;
    case Parameterization::NONCENTERED:
    eps_star = y_star - mu - sqrt(sigma2)*ht;
    eta = ht.tail(n-1) - phi*ht.head(n-1);
    break;
  }

  //Rcpp::Rcout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;
  
  static const arma::vec::fixed<10> mix_log_prob = arma::log(mix_prob);
  static const arma::vec::fixed<10> likelihood_normalizer = 0.5 * arma::log(2 * arma::datum::pi * mix_var);
  static const arma::vec::fixed<10> help_eta_mean = rho * std::sqrt(sigma2_used) * arma::exp(0.5 * mix_mean);
  const double log_eta_coefficient = -0.5 / (sigma2_used * (1 - rho * rho));
  const double log_eta_constant = -0.5 * std::log(2 * arma::datum::pi * sigma2_used * (1 - rho * rho));
  for (int r = 0; r < n; r++) {
    arma::vec::fixed<mix_count> post_dist;
    for (int c = 0; c < mix_count; c++) {
      const double a = mix_a[c];
      const double b = mix_b[c];
      const double m = mix_mean[c];
      const double v2 = mix_var[c];
      const double log_prior = mix_log_prob[c];
      
      double log_eps_star_lik = -0.5 * std::pow((eps_star[r] - m), 2) / v2 - likelihood_normalizer[c];
      double log_eta_lik;
      if (r < n - 1) {
        log_eta_lik = log_eta_coefficient * std::pow(eta[r] - d[r] * help_eta_mean[c] * (a + b * (eps_star[r] - m)), 2) + log_eta_constant;
      } else {
        log_eta_lik = 0.0;
      }
      /*log_*/post_dist[c] = log_prior + log_eps_star_lik + log_eta_lik;
    }
    const double max_log_post_dist = arma::max(/*log_*/post_dist);
    post_dist = arma::cumsum(arma::exp(/*log_*/post_dist - max_log_post_dist));
    post_dist = post_dist / post_dist[mix_count-1];

    auto binary_search_result = std::lower_bound(post_dist.cbegin(), post_dist.cend(), R::unif_rand());
    new_states[r] = std::distance(post_dist.cbegin(), binary_search_result);
  }
  
  return new_states;
}

