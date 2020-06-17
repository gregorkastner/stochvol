#include <RcppArmadillo.h>
#include "sampling_latent_states.h"
#include "h-utils.h"
#include "auxmix.h"
#include "parameterization.h"

using namespace Rcpp;

arma::vec draw_latent(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const bool correct) {
  // Draw h from AUX
  const arma::uvec s = draw_s_auxiliary(y_star, d, h, ht, phi, rho, sigma2, mu, Parameterization::CENTERED);
  const arma::vec proposed = draw_h_auxiliary(y_star, d, s, phi, rho, sigma2, mu, priormu_mu, priormu_sigma, Parameterization::CENTERED);

  if (correct) {
    return correct_latent_auxiliaryMH(y, y_star, d, h, ht, proposed, phi, rho, sigma2, mu);
  } else {
    return proposed;
  }
}

arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::uvec& z,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering) {
  // arrays for
  //   A_z_{+; -}, X_z.beta, W_z_{+; -}.beta
  //   m_t, Lambda_t_inv, temp_t
  //   alpha_t (the result)
  // only double for
  //   A_ij_{t-1; t}
  //   W_{t-1; t}.beta, X_{t-1; t}.beta
  //   Omega_{t-1,t; tt; t,t+1}

  // Compile time constants
  static const arma::vec::fixed<10> exp_m2 = arma::exp(.5 * mix_mean);
  static const arma::vec::fixed<10> b_exp_m2 = mix_b % exp_m2;
  static const arma::vec::fixed<10> b2_exp_m = arma::square(b_exp_m2);

  // Runtime constants
  const int n = y_star.n_elem;
  const double sigma = std::sqrt(sigma2),
               phi2 = std::pow(phi, 2),
               rho_sigma = rho * sigma,
               help_mu_phi = mu * (1 - phi),
               a_1 = mu,
               P_1_inv = (1 - phi2) / sigma2;

  // A as the function of z
  const double A_z_22 = 1 / (sigma2 * (1 - std::pow(rho, 2)));
  const arma::vec::fixed<10> A_z_11 = mix_varinv - rho_sigma * b2_exp_m * A_z_22,
                             A_z_12 = rho_sigma * b_exp_m2 * A_z_22;  // must be multiplied by d_t
  //const arma::vec::fixed<10>& A_z_21 = A_z_12;

  // X.beta and W.beta as the function of z
  const arma::vec::fixed<10>& X_z__beta = mix_mean,
                             help_W_z = rho_sigma * mix_a % exp_m2;

  // Partial calculations
  double A_tm1_11, A_t_11,
         A_tm1_12, A_t_12,
         //A_tm1_21, A_t_21,
         A_tm1_22, A_t_22;
  double W_tm1__beta, W_t__beta;
  double X_tm1__beta, X_t__beta;
  double Omega_tm1_t, Omega_tt, Omega_t_tp1;
  double c_t;

  // Stored partial calculations
  arma::vec Lambda_inv(n),
            m(n),
            help_Omega_Lambda(n - 1);

  // Result object
  arma::vec alpha(n);

  // McCausland pre-calculations
  {
    const int t = 0;
    const int z_t = z[t];
    const int d_t = d[t];

    A_t_11 = A_z_11[z_t];
    A_t_12 = d_t * A_z_12[z_t];
    //A_t_21 = A_t_12;
    A_t_22 = A_z_22;

    X_t__beta = X_z__beta[z_t];
    W_t__beta = help_mu_phi + d_t * help_W_z[z_t];

    Omega_tt = A_t_11 + 2 * phi * A_t_12 + phi2 * A_t_22 + P_1_inv;
    Omega_t_tp1 = -A_t_12 - phi * A_t_22;
    c_t = (A_t_11 + phi * A_t_12) * (y_star[t] - X_t__beta) + Omega_t_tp1 * W_t__beta + a_1 * P_1_inv;

    Lambda_inv[t] = 1 / std::sqrt(Omega_tt);
    help_Omega_Lambda[t] = Omega_t_tp1 * Lambda_inv[t];
    m[t] = c_t * std::pow(Lambda_inv[t], 2);
  }

  for (int t = 1; t < n - 1; t++) {
    const int z_t = z[t];
    const int d_t = d[t];

    // Running partial results
    A_tm1_11 = A_t_11;
    A_tm1_12 = A_t_12;
    //A_tm1_21 = A_t_21;
    A_tm1_22 = A_t_22;
    W_tm1__beta = W_t__beta;
    X_tm1__beta = X_t__beta;
    Omega_tm1_t = Omega_t_tp1;

    // New partial results
    A_t_11 = A_z_11[z_t];
    A_t_12 = d_t * A_z_12[z_t];
    //A_t_21 = A_t_12;
    A_t_22 = A_z_22;

    X_t__beta = X_z__beta[z_t];
    W_t__beta = help_mu_phi + d_t * help_W_z[z_t];

    Omega_tt = A_t_11 + 2 * phi * A_t_12 + phi2 * A_t_22 + A_tm1_22;
    Omega_t_tp1 = -A_t_12 - phi * A_t_22;
    c_t = (A_t_11 + phi * A_t_12) * (y_star[t] - X_t__beta) + Omega_t_tp1 * W_t__beta -
      A_tm1_12 * (y_star[t - 1] - X_tm1__beta) + A_tm1_22 * W_tm1__beta;

    Lambda_inv[t] = 1 / std::sqrt(Omega_tt - std::pow(help_Omega_Lambda[t - 1], 2));
    help_Omega_Lambda[t] = Omega_t_tp1 * Lambda_inv[t];
    m[t] = (c_t - Omega_tm1_t * m[t - 1]) * std::pow(Lambda_inv[t], 2);
  }

  {
    const int t = n - 1;
    const int z_t = z[t];
    const int d_t = d[t];

    // Running partial results
    A_tm1_11 = A_t_11;
    A_tm1_12 = A_t_12;
    //A_tm1_21 = A_t_21;
    A_tm1_22 = A_t_22;
    W_tm1__beta = W_t__beta;
    X_tm1__beta = X_t__beta;
    Omega_tm1_t = Omega_t_tp1;

    // New partial results
    A_t_11 = A_z_11[z_t];
    A_t_12 = d_t * A_z_12[z_t];
    //A_t_21 = A_t_12;
    A_t_22 = A_z_22;

    X_t__beta = X_z__beta[z_t];
    W_t__beta = help_mu_phi + d_t * help_W_z[z_t];

    Omega_tt = mix_varinv[z_t] + A_tm1_22;
    c_t = mix_varinv[z_t] * (y_star[t] - X_t__beta) - A_tm1_12 * (y_star[t - 1] - X_tm1__beta) + A_tm1_22 * W_tm1__beta;

    Lambda_inv[t] = 1 / std::sqrt(Omega_tt - std::pow(help_Omega_Lambda[t - 1], 2));
    m[t] = (c_t - Omega_tm1_t * m[t - 1]) * std::pow(Lambda_inv[t], 2);
  }

  // McCausland smoothing
  {
    const int t = n - 1;
    alpha[t] = m[t] + R::norm_rand() * Lambda_inv[t];
  }

  for (int t = n - 2; t >= 0; t--) {
    alpha[t] = m[t] + (R::norm_rand() - help_Omega_Lambda[t] * alpha[t + 1]) * Lambda_inv[t];
  }

  return alpha;
}

arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& proposed,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
    //const CharacterVector centering,

  // Calculate MH acceptance ratio
  const double hlp1 = h_log_posterior(proposed, y, phi, rho, sigma2, mu);
  const double hlp2 = h_log_posterior(h, y, phi, rho, sigma2, mu);
  const double halp1 = h_aux_log_posterior(proposed, y_star, d, phi, rho, sigma2, mu);
  const double halp2 = h_aux_log_posterior(h, y_star, d, phi, rho, sigma2, mu);
  const double log_acceptance = hlp1-hlp2-(halp1-halp2);
  arma::vec result;
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = proposed;
  } else {
    result = h;
  }

  return result;
}

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

