#include <RcppArmadillo.h>
#include <algorithm>
#include "progutils.h"
#include "sampler-leverage.h"
#include "theta-sampler.h"
#include "h-sampler.h"
#include "parameterization.hpp"

using namespace Rcpp;

Rcpp::List svlsample_cpp (
    const Rcpp::NumericVector& y_in,
    const int draws,
    const int burnin,
    const Rcpp::NumericMatrix& X,
    const int thinpara,
    const int thinlatent,
    const int thintime,
    const Rcpp::List& theta_init,
    const Rcpp::NumericVector& h_init,
    const double prior_phi_a,
    const double prior_phi_b,
    const double prior_rho_a,
    const double prior_rho_b,
    const double prior_sigma2_shape,
    const double prior_sigma2_rate,
    const double prior_mu_mu,
    const double prior_mu_sigma,
    const double prior_beta_mu,
    const double prior_beta_sigma,
    const bool verbose,
    const double offset,
    const double stdev,
    const bool gammaprior,
    const Rcpp::CharacterVector& strategy_rcpp) {

  const int N = burnin + draws;
  const bool regression = !ISNA(X.at(0,0));
  const int T = y_in.length();
  const int p = X.ncol();

  NumericVector y = y_in;
  NumericVector y_star = Rcpp::log(y*y + offset);
  NumericVector d(T); std::transform(y_in.cbegin(), y_in.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });

  double phi = theta_init["phi"];
  double rho = theta_init["rho"];
  double sigma2 = pow(theta_init["sigma"], 2);
  double mu = theta_init["mu"];
  NumericVector h = h_init, ht = (h_init-mu)/sqrt(sigma2);
  arma::vec beta(p); beta.fill(0.0);

  arma::mat betas(regression * draws/thinpara, p, arma::fill::zeros);
  NumericMatrix params(draws/thinpara, 4);
  NumericMatrix latent(draws/thinlatent, T/thintime);

  // priors in objects
  const NumericVector prior_phi = {prior_phi_a, prior_phi_b};
  const NumericVector prior_rho = {prior_rho_a, prior_rho_b};
  const NumericVector prior_sigma2 = {prior_sigma2_shape, prior_sigma2_rate};
  const NumericVector prior_mu = {prior_mu_mu, prior_mu_sigma};

  // don't use strings or RcppCharacterVector
  Rcpp::IntegerVector strategy(strategy_rcpp.length());
  std::transform(strategy_rcpp.cbegin(), strategy_rcpp.cend(), strategy.begin(),
      [](const SEXP& par) -> int {
        if (as<std::string>(par) == "centered") return int(Parameterization::CENTERED);
        else if (as<std::string>(par) == "non-centered") return int(Parameterization::NONCENTERED);
        else Rf_error("Illegal parameterization");
  });

  // some stuff for the regression part
  // prior mean and precision matrix for the regression part (currently fixed)
  const arma::vec y_in_arma(y_in.begin(), T);
  const arma::vec priorbetamean = arma::ones(p) * prior_beta_mu;
  const arma::mat priorbetaprec = arma::eye(p, p) / pow(prior_beta_sigma, 2);
  arma::vec normalizer(T);
  arma::mat X_reg(T, p);
  arma::vec y_reg(T);
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);
  const arma::vec h_arma(h.begin(), h.length(), false);  // create view
  const arma::vec ht_arma(ht.begin(), ht.length(), false);  // create view

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  const int show = verbose ? progressbar_init(N) : 0;

  for (int i = -burnin+1; i < draws+1; i++) {
    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {  // slightly circumstantial due to the combined use of Rcpp and arma
      std::copy(X.cbegin(), X.cend(), X_reg.begin());  // important!
      y = y_in_arma - X_reg*beta;
      y_star = Rcpp::log(y*y + offset);
      std::transform(y.cbegin(), y.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });
    }

    // update theta and h
    update_svl (y, y_star, d,
      phi, rho, sigma2, mu,
      h, ht,
      prior_phi, prior_rho,
      prior_sigma2, prior_mu,
      stdev, gammaprior, strategy);

    // update beta
    if (regression) {
      y_reg = y_in_arma;
      y_reg.head(T-1) -= rho * (arma::exp(h_arma.head(T-1)/2) % (ht_arma.tail(T-1) - phi*ht_arma.head(T-1)));

      normalizer = arma::exp(-h_arma/2);
      normalizer.head(T-1) /= sqrt(1-pow(rho, 2));
      // X has already been copied to X_reg
      X_reg.each_col() %= normalizer;
      y_reg %= normalizer;

      // cholesky factor of posterior precision matrix
      postprecchol = arma::chol(X_reg.t() * X_reg + priorbetaprec);

      // inverse cholesky factor of posterior precision matrix 
      postpreccholinv = arma::inv(arma::trimatu(postprecchol));

      // posterior covariance matrix and posterior mean vector
      postcov = postpreccholinv * postpreccholinv.t();
      postmean = postcov * (X_reg.t() * y_reg + priorbetaprec * priorbetamean);

      armadraw.imbue([]() -> double {return R::rnorm(0, 1);});  // equivalent to armadraw = Rcpp::rnorm(p); but I don't know if rnorm creates a vector

      // posterior betas
      beta = postmean + postpreccholinv * armadraw;
    }

    // store draws
    if ((i >= 1) && !thinpara_round) {
      params.at(i/thinpara-1, 0) = mu;
      params.at(i/thinpara-1, 1) = phi;
      params.at(i/thinpara-1, 2) = sqrt(sigma2);
      params.at(i/thinpara-1, 3) = rho;
      if (regression) {
        betas.row(i/thinpara-1) = beta.t();
      }
    }
    if ((i >= 1) && !thinlatent_round) {
      for (int volind = 0, thincol = thintime-1; thincol < h.length(); volind++, thincol += thintime) {
        latent.at(i/thinlatent-1, volind) = h[thincol];
      }
    }
  }

  if (verbose) progressbar_finish(N);  // finalize progress bar

  return Rcpp::List::create(
      Rcpp::_["para"] = params,
      Rcpp::_["latent"] = latent,
      Rcpp::_["beta"] = betas);
}

void update_svl (
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    Rcpp::NumericVector& h,
    Rcpp::NumericVector& ht,
    const Rcpp::NumericVector& prior_phi,
    const Rcpp::NumericVector& prior_rho,
    const Rcpp::NumericVector& prior_sigma2,
    const Rcpp::NumericVector& prior_mu,
    const double stdev,
    const bool gammaprior,
    const Rcpp::IntegerVector& strategy) {

  // only centered
  h = draw_latent_auxiliaryMH(y, y_star, d, h, ht, phi, rho, sigma2, mu, prior_mu[0], prior_mu[1]);
  ht = (h-mu)/sqrt(sigma2);

  for (int ipar : strategy) {
    const Parameterization par = Parameterization(ipar);
    switch (par) {
      case Parameterization::CENTERED:
        draw_theta_rwMH(
            phi, rho, sigma2, mu, y, h, ht,
            prior_phi,
            prior_rho,
            prior_sigma2,
            prior_mu,
            par,
            stdev,
            gammaprior);
        break;
      case Parameterization::NONCENTERED:
        draw_theta_rwMH(
            phi, rho, sigma2, mu, y, h, ht,
            prior_phi,
            prior_rho,
            prior_sigma2,
            prior_mu,
            par,
            stdev,
            gammaprior);
        break;
    }

    switch (par) {
      case Parameterization::CENTERED:
        ht = (h-mu)/sqrt(sigma2);
        break;
      case Parameterization::NONCENTERED:
        h = sqrt(sigma2)*ht+mu;
        break;
    }
  }
  return;
}

