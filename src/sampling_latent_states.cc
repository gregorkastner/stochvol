/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2018
 *     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

/*
 * sampling_latent_states.cc
 * 
 * Definitions of the functions declared in sampling_latent_states.h.
 * Documentation can also be found in sampling_latent_states.h.
 */

#include <RcppArmadillo.h>
#include <expert.hpp>
#include "densities.h"
#include "sampling_latent_states.h"
#include "utils.h"
#include "utils_latent_states.h"

using namespace Rcpp;

namespace stochvol {

namespace fast_sv {

// Main sampling function for the latent states.
// See documentation above the declaration
LatentVector draw_latent(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  const arma::vec& y = data;  // rename
  const unsigned int T = y.n_elem;

  double omega_offdiag;  // contains off-diag element of precision matrix (const)
  arma::vec omega_diag(T+1),  // contains diagonal elements of precision matrix
    covector(T+1);  // holds covector (see McCausland et al. 2011)
  const double sigma2 = std::pow(sigma, 2),
               sigma2inv = 1. / sigma2,
               Bh0inv = determine_Bh0inv(phi, prior_spec);

  switch (expert.baseline) {
    case Parameterization::CENTERED:  // fill precision matrix omega and covector c for CENTERED para:
      {
        const double phi2 = std::pow(phi, 2);

        omega_diag[0] = (Bh0inv + phi2) * sigma2inv;
        covector[0] = mu * (Bh0inv - phi*(1-phi)) * sigma2inv;

        for (unsigned int j = 1; j < T; j++) {
          omega_diag[j] = mix_varinv[r[j-1]] + (1+phi2)*sigma2inv; 
          covector[j] = (data[j-1] - mix_mean[r[j-1]])*mix_varinv[r[j-1]]
            + mu*(1-phi)*(1-phi)*sigma2inv;
        }
        omega_diag[T] = mix_varinv[r[T-1]] + sigma2inv;
        covector[T] = (data[T-1] - mix_mean[r[T-1]])*mix_varinv[r[T-1]] + mu*(1-phi)*sigma2inv;
        omega_offdiag = -phi*sigma2inv;
      }
      break;
    case Parameterization::NONCENTERED:  // fill precision matrix omega and covector c for NONCENTERED para:
      {
        const double phi2 = std::pow(phi, 2);

        omega_diag[0] = Bh0inv + phi2;
        covector[0] = 0.;

        for (unsigned int j = 1; j < T; j++) {
          omega_diag[j] = mix_varinv[r[j-1]] * sigma2 + 1 + phi2;
          covector[j] = mix_varinv[r[j-1]] * sigma * (data[j-1] - mix_mean[r[j-1]] - mu);
        }
        omega_diag[T] = mix_varinv[r[T-1]] * sigma2 + 1;
        covector[T] = mix_varinv[r[T-1]] * sigma * (data[T-1] - mix_mean[r[T-1]] - mu);
        omega_offdiag = -phi;
      }
      break;
    default:
      ::Rf_error("draw_latent: This part of the code should never be reached.");
  } 

  // Cholesky decomposition
  const auto cholesky_matrix = cholesky_tridiagonal(omega_diag, omega_offdiag);
  const arma::vec& chol_diag = cholesky_matrix.chol_diag;
  const arma::vec& chol_offdiag = cholesky_matrix.chol_offdiag;

  // Solution of Chol*x = covector ("forward algorithm")
  arma::vec htmp = forward_algorithm(chol_diag, chol_offdiag, covector);
  htmp.transform( [](const double h_elem) -> double { return h_elem + R::norm_rand(); });

  // Solution of (Chol')*x = htmp ("backward algorithm")
  const arma::vec hnew = backward_algorithm(chol_diag, chol_offdiag, htmp);

  return {hnew[0], hnew.tail(T)};
}

// Main sampling function for the mixture indicators.
// See documentation above the declaration
arma::uvec draw_mixture_indicators(
    const arma::vec& data,
    const double mu,
    const double phi,
    const double sigma,
    const arma::vec& h) {
  const int T = data.n_elem;

  // calculate non-normalized CDF of posterior indicator probs
  const arma::vec mixprob = find_mixture_indicator_cdf(data-h);

  // find correct indicators
  const arma::uvec r = inverse_transform_sampling(mixprob, T);

  return r;
}

}  // END namespace fast_sv

namespace general_sv {

// Sample from the proposal distribution of the latent
// states. It is the conditional posterior distribution
// for the latent states in the conditionally Gaussian
// state space model by Omori et al (2007).
// Parameter 'centering' sets the parameterization.
arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::uvec& z,
    const double h0,
    const Parameterization centering);  // currently only CENTERED

// Compute the Metropolis-Hasting acceptance-rejection
// ratio and accept or reject 'proposed' against the current
// draw 'h'.
arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const double h0,
    const arma::vec& h,
    const arma::vec& proposed,
    const double s_integral);

// Output struct for 'draw_s_auxiliary'
struct Draw_s_auxiliary_out {
  arma::uvec s;
  double s_integral;
};

// Sample the vector of mixture indicators. For more
// information see Omori et al (2007).
Draw_s_auxiliary_out draw_s_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::vec& h,
    const arma::vec& ht,
    const Parameterization centering,
    const bool do_compute_s_integral);

// Main sampling function for the latent states.
// See documentation above the declaration
LatentVector draw_latent(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::vec& h,
    const arma::vec& ht,
    const PriorSpec& prior_spec,
    const ExpertSpec_GeneralSV& expert) {
  // Draw h0 | h1, mu, phi, sigma
  const double phi2 = std::pow(phi, 2),
               sigma2 = std::pow(sigma, 2),
               B02 = sigma2 / determine_Bh0inv(phi, prior_spec),
               h1 = h[0],
               denominator = B02 * phi2 + sigma2,
               mean = (mu * sigma2 + phi * B02 * (h1 - mu * (1 - phi))) / denominator,
               var = sigma2 * B02 / denominator;
  const double h0 = R::rnorm(mean, std::sqrt(var));

  // Draw h from AUX
  const auto draw_s_out = draw_s_auxiliary(y_star, d, mu, phi, sigma, rho, h, ht, Parameterization::CENTERED, expert.correct_latent_draws);
  const arma::vec proposed = draw_h_auxiliary(y_star, d, mu, phi, sigma, rho, draw_s_out.s, h0, Parameterization::CENTERED);

  return {h0,
    expert.correct_latent_draws ?
      correct_latent_auxiliaryMH(y, y_star, d, mu, phi, sigma, rho, h0, h, proposed, draw_s_out.s_integral) :
      proposed};
}

// See documentation above the declaration
arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::uvec& z,
    const double h0,
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
  static const arma::vec::fixed<10> exp_m_half = arma::exp(.5 * mix_mean);
  static const arma::vec::fixed<10> b_exp_m2 = mix_b % exp_m_half;
  static const arma::vec::fixed<10> b2_exp_m = arma::square(b_exp_m2);

  // Runtime constants
  const int n = y_star.n_elem;
  const double sigma2 = std::pow(sigma, 2),
               phi2 = std::pow(phi, 2),
               rho2 = std::pow(rho, 2),
               rho_const = 1 / (1 - rho2),
               help_mu_phi = mu * (1 - phi),
               a_1 = help_mu_phi + phi * h0,
               P_1_inv = 1 / sigma2;

  // A as the function of z
  const double A_z_22 = rho_const / sigma2;
  const arma::vec::fixed<10> A_z_11 = mix_varinv + rho2 * b2_exp_m * rho_const,
                             A_z_12 = -rho * b_exp_m2 * rho_const / sigma;  // must be multiplied by d_t
  //const arma::vec::fixed<10>& A_z_21 = A_z_12;

  // X.beta and W.beta as the function of z
  const arma::vec::fixed<10>& X_z__beta = mix_mean,
                              help_W_z = rho * sigma * mix_a % exp_m_half;

  // Partial calculations
  double /*A_tm1_11,*/ A_t_11,
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
  // TODO Omega and c can be done in parallel
  {
    const int t = 0;
    const unsigned int z_t = z[t];
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
    const unsigned int z_t = z[t];
    const int d_t = d[t];

    // Running partial results
    //A_tm1_11 = A_t_11;
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
    const int unsigned z_t = z[t];
    const int d_t = d[t];

    // Running partial results
    //A_tm1_11 = A_t_11;
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

// See documentation above the declaration
arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const double h0,
    const arma::vec& h,
    const arma::vec& proposed,
    const double s_integral) {
  //const CharacterVector centering,

  // Calculate MH acceptance ratio
  const double hlp1 = h_log_posterior(proposed, y, phi, rho, sigma, mu, h0);
  const double hlp2 = h_log_posterior(h, y, phi, rho, sigma, mu, h0);
  const double halp1 = h_aux_log_posterior(proposed, y_star, d, phi, rho, sigma, mu, h0);
  const double halp2 = s_integral + logdnorm2(h[0], mu + phi * (h0 - mu), sigma);  //h_aux_log_posterior(h, y_star, d, phi, rho, sigma, mu, h0);
  const double log_acceptance = hlp1-hlp2-(halp1-halp2);
  arma::vec result;
  if (log_acceptance > 0 || std::exp(log_acceptance) > R::unif_rand()) {
    result = proposed;
  } else {
    result = h;
  }

  return result;
}

// See documentation above the declaration
Draw_s_auxiliary_out draw_s_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const double mu,
    const double phi,
    const double sigma,
    const double rho,
    const arma::vec& h,
    const arma::vec& ht,
    const Parameterization centering,
    const bool do_compute_s_integral) {
  const int n = y_star.size();
  const double sigma2_used = centering == Parameterization::CENTERED ? std::pow(sigma, 2) : 1.0;
  static const int mix_count = mix_a.n_elem;
  arma::vec eps_star;
  arma::vec eta;
  arma::uvec new_states(n);
  double s_integral = 0;

  switch (centering) {
    case Parameterization::CENTERED:
      eps_star = y_star - h;
      eta = (h.tail(n-1) - mu) - phi*(h.head(n-1) - mu);
      break;
    case Parameterization::NONCENTERED:
      eps_star = y_star - mu - sigma*ht;
      eta = ht.tail(n-1) - phi*ht.head(n-1);
      break;
  }

  static const arma::vec::fixed<10> mix_log_prob = arma::log(mix_prob);
  static const arma::vec::fixed<10> likelihood_normalizer = 0.5 * arma::log(2 * arma::datum::pi * mix_var);
  static const arma::vec::fixed<10> exp_m_half = arma::exp(0.5 * mix_mean);

  const arma::vec::fixed<10> help_eta_mean = rho * std::sqrt(sigma2_used) * exp_m_half;
  const double log_eta_coefficient = -0.5 / (sigma2_used * (1 - rho * rho));
  const double log_eta_constant = -0.5 * std::log(2 * arma::datum::pi * sigma2_used * (1 - rho * rho));
  for (int t = 0; t < n; t++) {
    arma::vec::fixed<mix_count> post_dist;
    for (int c = 0; c < mix_count; c++) {
      const double a = mix_a[c];
      const double b = mix_b[c];
      const double m = mix_mean[c];
      const double v2 = mix_var[c];
      const double log_prior = mix_log_prob[c];

      const double log_eps_star_lik = -0.5 * std::pow((eps_star[t] - m), 2) / v2 - likelihood_normalizer[c];
      const double log_eta_lik = t == n - 1 ?
        0 :
        log_eta_coefficient * std::pow(eta[t] - d[t] * help_eta_mean[c] * (a + b * (eps_star[t] - m)), 2) + log_eta_constant;
      /*log_*/post_dist[c] = log_prior + log_eps_star_lik + log_eta_lik;
    }
    s_integral += do_compute_s_integral ? std::log(arma::sum(arma::exp(/*log_*/post_dist))) : 0;
    const double max_log_post_dist = arma::max(/*log_*/post_dist);
    post_dist = arma::cumsum(arma::exp(/*log_*/post_dist - max_log_post_dist));
    post_dist = post_dist / post_dist[mix_count-1];

    const auto binary_search_result = std::lower_bound(post_dist.cbegin(), post_dist.cend(), R::unif_rand());
    new_states[t] = std::distance(post_dist.cbegin(), binary_search_result);
  }

  return {new_states, s_integral};
}

}  // END namespace general_sv

}

