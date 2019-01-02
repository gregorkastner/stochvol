#include "theta-sampler.h"
#include "h-sampler.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix svlsample_cpp (
    const int n,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector d,
    const int thin,
    const int burnin,
    const double phi_init,
    const double rho_init,
    const double sigma2_init,
    const double mu_init,
    const Rcpp::NumericVector& h_init,
    const double prior_phi_a,
    const double prior_phi_b,
    const double prior_rho_a,
    const double prior_rho_b,
    const double prior_sigma2_shape,
    const double prior_sigma2_rate,
    const double prior_mu_mu,
    const double prior_mu_sigma2,
    const double stdev,
    const Rcpp::CharacterVector& strategy,
    const Rcpp::DataFrame& mixing_constants) {

  NumericVector h = h_init, ht = (h_init-mu_init)/sqrt(sigma2_init);
  double phi = phi_init, rho = rho_init, sigma2 = sigma2_init, mu = mu_init;
  NumericVector theta = {phi, rho, sigma2, mu};

  NumericMatrix params(n/thin, 4);
  NumericMatrix vols(n/thin, y.length());

  for (int i = -burnin+1; i < n+1; i++) {
    const bool thinning_round = (thin > 1) && (i % thin != 0);  // is this a thinning round?
    const int i_plus_burnin = i+burnin;
    const int n_plus_burnin = n+burnin;
    if (i_plus_burnin % (n_plus_burnin/10) == 0) {
      if (i_plus_burnin <= burnin) {
        Rcout << "Burnin:   ";
      } else {
        Rcout << "Sampling: ";
      }
      Rcout << i_plus_burnin/n_plus_burnin * 100 << "%" << std::endl;
    }

    // only centered
    h = draw_latent_auxiliaryMH(y, y_star, d, h, phi, rho, sigma2, mu, mixing_constants);
    ht = (h-mu)/sqrt(sigma2);

    for (int ind_strategy = 0; ind_strategy < strategy.length(); ind_strategy++) {
      if (as<std::string>(strategy(ind_strategy)) == "centered") {
        theta = draw_theta_rwMH(phi, rho, sigma2, mu, y, h,
            NumericVector::create(prior_phi_a, prior_phi_b),
            NumericVector::create(prior_rho_a, prior_rho_b),
            NumericVector::create(prior_sigma2_shape, prior_sigma2_rate),
            NumericVector::create(prior_mu_mu, prior_mu_sigma2),
            wrap("centered"),
            stdev);
      } else if (as<std::string>(strategy(ind_strategy)) == "non-centered") {
        theta = draw_theta_rwMH(phi, rho, sigma2, mu, y, ht,
            NumericVector::create(prior_phi_a, prior_phi_b),
            NumericVector::create(prior_rho_a, prior_rho_b),
            NumericVector::create(prior_sigma2_shape, prior_sigma2_rate),
            NumericVector::create(prior_mu_mu, prior_mu_sigma2),
            wrap("non-centered"),
            stdev);
      } else {
        Rf_error("problem1");
      }

      phi = theta(0);
      rho = theta(1);
      sigma2 = theta(2);
      mu = theta(3);

      if (as<std::string>(strategy(ind_strategy)) == "centered") {
        ht = (h-mu)/sqrt(sigma2);
      } else if (as<std::string>(strategy(ind_strategy)) == "non-centered") {
        h = sqrt(sigma2)*ht+mu;
      } else {
        Rf_error("problem2");
      }
    }

    if ((i >= 1) && !thinning_round) {
      params.row(i/thin-1) = theta;
      vols.row(i/thin-1) = Rcpp::exp(h/2);
    }
  }

  return cbind(params, vols);
}

