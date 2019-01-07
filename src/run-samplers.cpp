#include <algorithm>
#include "progutils.h"
#include "run-samplers.h"
#include "theta-sampler.h"
#include "h-sampler.h"
#include "parameterization.hpp"

using namespace Rcpp;

Rcpp::List svlsample_cpp (
    const int draws,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    const int burnin,
    const int thinpara,
    const int thinlatent,
    const int thintime,
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
    const double prior_mu_sigma,
    const bool verbose,
    const double stdev,
    const bool gammaprior,
    const Rcpp::CharacterVector& strategy_rcpp) {

  const int N = burnin + draws;

  NumericVector h = h_init, ht = (h_init-mu_init)/sqrt(sigma2_init);
  double phi = phi_init, rho = rho_init, sigma2 = sigma2_init, mu = mu_init;
  NumericVector theta = {phi, rho, sigma2, mu};

  NumericMatrix params(draws/thinpara, 4);
  NumericMatrix latent(draws/thinlatent, y.length()/thintime);

  // don't use strings or RcppCharacterVector
  std::vector<Parameterization> strategy(strategy_rcpp.length());
  std::transform(strategy_rcpp.cbegin(), strategy_rcpp.cend(), strategy.begin(),
      [](const SEXP& par) -> Parameterization {
        if (as<std::string>(par) == "centered") return Parameterization::CENTERED;
        else if (as<std::string>(par) == "non-centered") return Parameterization::NONCENTERED;
        else Rf_error("Illegal parameterization");
  });

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  const int show = verbose ? progressbar_init(N) : 0;

  for (int i = -burnin+1; i < draws+1; i++) {
    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    // only centered
    h = draw_latent_auxiliaryMH(y, y_star, d, h, ht, phi, rho, sigma2, mu, prior_mu_mu, prior_mu_sigma);
    ht = (h-mu)/sqrt(sigma2);

    for (auto par : strategy) {
      switch (par) {
        case Parameterization::CENTERED:
        theta = draw_theta_rwMH(phi, rho, sigma2, mu, y, h,
            NumericVector::create(prior_phi_a, prior_phi_b),
            NumericVector::create(prior_rho_a, prior_rho_b),
            NumericVector::create(prior_sigma2_shape, prior_sigma2_rate),
            NumericVector::create(prior_mu_mu, prior_mu_sigma),
            par,
            stdev,
            gammaprior);
        break;
        case Parameterization::NONCENTERED:
        theta = draw_theta_rwMH(phi, rho, sigma2, mu, y, ht,
            NumericVector::create(prior_phi_a, prior_phi_b),
            NumericVector::create(prior_rho_a, prior_rho_b),
            NumericVector::create(prior_sigma2_shape, prior_sigma2_rate),
            NumericVector::create(prior_mu_mu, prior_mu_sigma),
            par,
            stdev,
            gammaprior);
        break;
      }

      phi = theta(0);
      rho = theta(1);
      sigma2 = theta(2);
      mu = theta(3);

      switch (par) {
        case Parameterization::CENTERED:
        ht = (h-mu)/sqrt(sigma2);
        break;
        case Parameterization::NONCENTERED:
        h = sqrt(sigma2)*ht+mu;
        break;
      }
    }

    if ((i >= 1) && !thinpara_round) {
      params.row(i/thinpara-1) = theta;
    }
    if ((i >= 1) && !thinlatent_round) {
      for (int volind = 0, thincol = thintime-1; thincol < h.length(); volind++, thincol += thintime) {
        latent(i/thinlatent-1, volind) = h(thincol);
      }
    }
  }

  if (verbose) progressbar_finish(N);  // finalize progress bar

  return Rcpp::List::create(
      Rcpp::_["para"] = params,
      Rcpp::_["latent"] = latent);
}

