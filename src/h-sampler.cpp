#include <RcppArmadillo.h>
#include "h-sampler.h"
#include "h-utils.h"
#include "mixture-state-sampler.h"
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
  // TODO optimize draw_s_auxiliary
  const arma::uvec s = draw_s_auxiliary(y_star, d, h, ht, phi, rho, sigma2, mu, Parameterization::CENTERED);
  Rcpp::Rcout << "run draw_latent" << std::endl;
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
  const arma::vec::fixed<10> A_z_11 = 2 * mix_2varinv - rho_sigma * b2_exp_m * A_z_22,
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

    Rcpp::Rcout << Omega_tt << ' ';
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

    Rcpp::Rcout << Omega_tt << ' ';
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

    Omega_tt = 2 * mix_2varinv[z_t] + A_tm1_22;
    c_t = 2 * mix_2varinv[z_t] * (y_star[t] - X_t__beta) - A_tm1_12 * (y_star[t - 1] - X_tm1__beta) + A_tm1_22 * W_tm1__beta;

    Lambda_inv[t] = 1 / std::sqrt(Omega_tt - std::pow(help_Omega_Lambda[t - 1], 2));
    m[t] = (c_t - Omega_tm1_t * m[t - 1]) * std::pow(Lambda_inv[t], 2);

    Rcpp::Rcout << Omega_tt << std::endl;
  }

  Rcpp::Rcout << Lambda_inv.t() << std::endl <<
    m.t() << std::endl << std::endl;

  // McCausland smoothing
  {
    const int t = n - 1;
    alpha[t] = m[t] + R::norm_rand() * Lambda_inv[t];
  }

  for (int t = n - 2; t >= 0; t--) {
    alpha[t] = m[t] + (R::norm_rand() - help_Omega_Lambda[t] * alpha[t + 1]) * Lambda_inv[t];
  }

  Rcpp::Rcout << alpha.t() << std::endl << std::endl << std::endl;
  return alpha;
}
 
//arma::vec draw_h_auxiliary_old(
//    const arma::vec& y_star,
//    const arma::ivec& d,
//    const arma::vec& s,
//    const double phi,
//    const double rho,
//    const double sigma2,
//    const double mu,
//    const double priormu_mu,
//    const double priormu_sigma,
//    const Parameterization centering) {
//  arma::vec mixing_a(s.size()); std::transform(s.cbegin(), s.cend(), mixing_a.begin(), [](const int selem) -> double {return mix_a[selem];});
//  arma::vec mixing_b(s.size()); std::transform(s.cbegin(), s.cend(), mixing_b.begin(), [](const int selem) -> double {return mix_b[selem];});
//  arma::vec mixing_m(s.size()); std::transform(s.cbegin(), s.cend(), mixing_m.begin(), [](const int selem) -> double {return mix_mean[selem];});
//  arma::vec mixing_v(s.size()); std::transform(s.cbegin(), s.cend(), mixing_v.begin(), [](const int selem) -> double {return sqrt(mix_var[selem]);});
//  
//  const List filter_result = aug_kalman_filter(phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, priormu_mu, pow(priormu_sigma, 2), centering);
//  
//  const List smoothing_result = simulation_smoother(mu, filter_result, centering);
//  const arma::vec eta = smoothing_result["eta"];
//  const double eta0 = as<NumericVector>(smoothing_result["eta0"])[0];
//
//  const int n = as<NumericVector>(filter_result["D"]).size();
//  arma::vec h(n, arma::fill::zeros);
//  arma::vec dt;
//  switch (centering) {
//    case Parameterization::CENTERED:
//    h[0] = mu + eta0;
//    dt = mu*(1-phi) + rho*sqrt(sigma2)*(d%mixing_a%exp(mixing_m/2));
//    break;
//    case Parameterization::NONCENTERED:
//    h[0] = eta0;
//    dt = rho*(d%mixing_a%exp(mixing_m/2));
//    break;
//  }
//
//  for (int i = 0; i < n-1; i++) {
//    h[i+1] = dt[i] + phi*h[i] + eta[i];
//  }
//
//  return h;
//}

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

