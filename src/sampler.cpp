#include <RcppArmadillo.h>
#include "sampler.h"
#include "update_functions.h"
#include "auxmix.h"
#include "progutils.h"
#include "regression.h"
#include "parameterization.hpp"

using namespace Rcpp;

Rcpp::List svsample_cpp(
    const Rcpp::NumericVector& y_in,
    const int draws,
    const int burnin,
    const Rcpp::NumericMatrix& X_in,
    const double bmu,
    const double Bmu,
    const double a0,
    const double b0,
    const double Bsigma,
    const int thin,
    const int timethin,
    const Rcpp::List& startpara_in,
    const Rcpp::NumericVector& startvol_in,
    const bool keeptau,
    const bool quiet,
    const int para,
    const int MHsteps,
    const double B011,
    const double B022,
    const double mhcontrol,
    const bool gammaprior,
    const bool truncnormal,
    const double offset,
    const bool dontupdatemu,
    const Rcpp::NumericVector& priordf_in,
    const Rcpp::NumericVector& priorbeta_in,
    const double priorlatent0) {

  NumericVector y(y_in), startvol(startvol_in);
  arma::vec armay(y.begin(), y.length(), false);

  NumericMatrix X(X_in);

  int T = y.length();
  int p = X.ncol();
  arma::mat armaX(X.begin(), T, p, false);

  // should we model the mean as well?
  bool regression;
  if (ISNA(X.at(0,0))) regression = false; else regression = true;
  NumericVector priorbeta(priorbeta_in);

  // prior mean and precision matrix for the regression part (currently fixed)
  arma::vec priorbetamean(p); priorbetamean.fill(priorbeta[0]);
  arma::mat priorbetaprec(p, p, arma::fill::zeros);
  priorbetaprec.diag() += 1./(priorbeta[1]*priorbeta[1]);

  List startpara(startpara_in);

  // number of MCMC draws
  int N 	    = burnin + draws;

  // verbosity control
  bool verbose = !quiet;

  // "expert" settings:
  double B011inv         = 1/B011;
  double B022inv         = 1/B022;
  bool Gammaprior        = gammaprior;
  double MHcontrol       = mhcontrol;
  int parameterization   = para;
  bool centered_baseline = parameterization % 2; // 0 for C, 1 for NC baseline

  // t-errors:
  bool terr;
  NumericVector priordf(priordf_in);

  if (ISNA(priordf[0])) terr = false; else terr = true;

  // moment-matched IG-prior
  double c0 = 2.5;
  double C0 = 1.5*Bsigma;

  // pre-calculation of a posterior parameter
  double cT = -10000; 
  if (Gammaprior) {  
    if (MHsteps == 2 || MHsteps == 3) cT = T/2.0; // we want IG(-.5,0) as proposal
    else if (MHsteps == 1) cT = (T-1)/2.0; // we want IG(-.5,0) as proposal
  }
  else {
    if (MHsteps == 2) cT = c0 + (T+1)/2.0;  // pre-calculation outside the loop
    else Rf_error("This setup is not yet implemented");
  }

  if (dontupdatemu == true && MHsteps == 1) { // not implemented (would be easy, though)
    Rf_error("This setup is not yet implemented");
  }

  // initialize the variables:
  NumericVector sigma2inv(N+1, pow(as<double>(startpara["sigma"]), -2));
  NumericVector phi(N+1, as<double>(startpara["phi"]));
  NumericVector mu(N+1, as<double>(startpara["mu"]));

  NumericVector h(T);  // contains h1 to hT, but not h0!
  for (int i = 0; i < T; i++) h[i] = startvol[i];  // need to make a manual copy (!)
  if (!centered_baseline) h = (h-mu[0])*sqrt(sigma2inv[0]);

  double h0;

  IntegerVector r(T);  // mixture indicators

  int hstorelength = (T-1)/timethin+1;
  NumericMatrix hstore(hstorelength, draws/thin);
  NumericVector h0store(draws/thin);

  NumericMatrix mixprob(10, T);  // mixture probabilities

  NumericVector data = log(y*y + offset);  // commonly used transformation
  arma::vec armadata(data.begin(), data.length(), false);

  NumericVector datastand = clone(data);  // standardized "data" (different for t-errors)
  arma::vec armadatastand(datastand.begin(), datastand.length(), false);

  NumericVector curpara(3);  // holds mu, phi, sigma in every iteration
  curpara[0] = mu[0];
  curpara[1] = phi[0];
  curpara[2] = 1/sqrt(sigma2inv[0]);

  // some stuff for the t-errors
  double nu;
  if (terr) nu = as<double>(startpara["nu"]);
  NumericVector tau(T, 1.);

  NumericVector nustore(terr * (N+1), nu);
  NumericMatrix taustore(hstorelength * keeptau, draws/thin);

  // some stuff for the regression part
  NumericVector curbeta(p, .012345);
  arma::vec armabeta(curbeta.begin(), curbeta.length(), false);
  arma::vec normalizer;
  arma::mat Xnew = as<arma::mat>(clone(X));
  arma::vec ynew = as<arma::vec>(clone(y));
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);
  NumericMatrix betastore(regression * (N+1), p);
  arma::mat armabetastore(betastore.begin(), betastore.nrow(), betastore.ncol(), false);

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) show = progressbar_init(N);

  for (int i = 0; i < N; i++) {  // BEGIN main MCMC loop

    // print a progress sign every "show" iterations
    if (verbose) if (!(i % show)) progressbar_print();

    if (regression) {
      armadatastand = armadata = log(square(armay - armaX * armabeta));
    }

    if (terr) {
      if (centered_baseline) update_terr(data - h, tau, nu, priordf[0], priordf[1]);
      else update_terr(data - curpara[0] - curpara[2] * h, tau, nu, priordf[0], priordf[1]);

      datastand = data - log(tau);
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_sv(datastand, curpara, h, h0, mixprob, r, centered_baseline, C0, cT,
        Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior,
        truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);

    if (regression) { // update betas (regression)
      normalizer = exp(-h/2);
      Xnew = armaX;
      Xnew.each_col() %= normalizer;
      ynew = as<arma::vec>(y) % normalizer;

      // cholesky factor of posterior precision matrix
      postprecchol = arma::chol(Xnew.t() * Xnew + priorbetaprec);  // TODO handle exceptions the R way

      // inverse cholesky factor of posterior precision matrix 
      postpreccholinv = arma::solve(arma::trimatu(postprecchol), arma::eye<arma::mat>(p, p));

      // posterior covariance matrix and posterior mean vector
      postcov = postpreccholinv * postpreccholinv.t();
      postmean = postcov * (Xnew.t() * ynew + priorbetaprec * priorbetamean);

      armadraw = rnorm(p);

      // posterior betas
      armabeta = postmean + postpreccholinv * armadraw;
    }

    // storage:
    if (!((i+1) % thin)) if (i >= burnin) {  // this means we should store h
      if (centered_baseline) {
        for (int j = 0; j < hstorelength; j++) hstore.at(j, (i-burnin)/thin) = h[timethin*j];
        h0store[(i-burnin)/thin] = h0;
      } else {
        for (int j = 0; j < hstorelength; j++) hstore.at(j, (i-burnin)/thin) = curpara[0] + curpara[2]*h[timethin*j];
        h0store[(i-burnin)/thin] = curpara[0] + curpara[2]*h0;
      }
      if (keeptau && terr) {
        for (int j = 0; j < hstorelength; j++) taustore.at(j, (i-burnin)/thin) = tau[timethin*j];
      }
    }
    mu[i+1] = curpara[0];
    phi[i+1] = curpara[1];
    sigma2inv[i+1] = 1/(curpara[2]*curpara[2]);
    if (terr) nustore[i+1] = nu;
    if (regression) armabetastore.row(i+1) = armabeta.t();
  }  // END main MCMC loop

  if (verbose) progressbar_finish(N);  // finalize progress bar

  // Prepare return value and return
  return cleanUp(mu, phi, sqrt(1/sigma2inv), hstore, h0store, nustore, taustore, betastore);
}

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

