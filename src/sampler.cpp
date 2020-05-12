#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "sampler.h"
#include "update_functions.h"
#include "auxmix.h"
#include "progutils.h"
#include "regression.h"
#include "parameterization.h"

using namespace Rcpp;

List svsample_cpp(
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X_in,
    const double bmu,
    const double Bmu,
    const double a0,
    const double b0,
    const double Bsigma,
    const int thin,
    const int thintime,
    const List& startpara,
    const arma::vec& startvol,
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
    const arma::vec& priordf_in,
    const arma::vec& priorbeta_in,
    const double priorlatent0) {

  arma::vec y(y_in);

  arma::mat X(X_in);

  int T = y.size();
  int p = X.n_cols;

  // should we model the mean as well?
  bool regression;
  if (ISNA(X.at(0,0))) regression = false; else regression = true;
  arma::vec priorbeta(priorbeta_in);

  // prior mean and precision matrix for the regression part (currently fixed)
  arma::vec priorbetamean(p); priorbetamean.fill(priorbeta[0]);
  arma::mat priorbetaprec(p, p, arma::fill::zeros);
  priorbetaprec.diag() += 1./(priorbeta[1]*priorbeta[1]);

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
  arma::vec priordf(priordf_in);

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
  arma::vec sigma2inv = arma::ones(N+1) * std::pow(as<double>(startpara["sigma"]), -2);
  arma::vec phi = arma::ones(N+1) * as<double>(startpara["phi"]);
  arma::vec mu = arma::ones(N+1) * as<double>(startpara["mu"]);

  arma::vec h = startvol;  // contains h1 to hT, but not h0!
  if (!centered_baseline) h = (h-mu[0])*sqrt(sigma2inv[0]);

  double h0;

  arma::ivec r(T);  // mixture indicators

  int hstorelength = T/thintime;  // thintime must either be 1 or T
  arma::mat hstore(hstorelength, draws/thin);
  arma::vec h0store(draws/thin);

  arma::mat mixprob(10, T);  // mixture probabilities
  arma::vec mixprob_vec(mixprob.begin(), mixprob.n_elem, false);

  arma::vec data = log(y%y + offset);  // commonly used transformation

  arma::vec datastand = data;  // standardized "data" (different for t-errors)

  arma::vec curpara(3);  // holds mu, phi, sigma in every iteration
  curpara[0] = mu[0];
  curpara[1] = phi[0];
  curpara[2] = 1/sqrt(sigma2inv[0]);

  // some stuff for the t-errors
  double nu = -1;
  if (terr) nu = as<double>(startpara["nu"]);
  arma::vec tau = arma::ones(T);

  arma::vec nustore = arma::ones(terr * (N+1)) * nu;
  arma::mat taustore(hstorelength * keeptau, draws/thin);

  // some stuff for the regression part
  arma::vec curbeta(p, arma::fill::ones);
  arma::vec normalizer;
  arma::mat Xnew = X;
  arma::vec ynew = y;
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);
  arma::mat betastore(regression * (N+1), p);

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) show = progressbar_init(N);

  for (int i = 0; i < N; i++) {  // BEGIN main MCMC loop
    R_CheckUserInterrupt();

    // print a progress sign every "show" iterations
    if (verbose) if (!(i % show)) progressbar_print();

    if (regression) {
      datastand = data = log(square(y - X*curbeta));
    }

    if (terr) {
      if (centered_baseline) update_terr(data - h, tau, nu, priordf[0], priordf[1]);
      else update_terr(data - curpara[0] - curpara[2] * h, tau, nu, priordf[0], priordf[1]);

      datastand = data - log(tau);
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_sv(datastand, curpara, h, h0, mixprob_vec, r, centered_baseline, C0, cT,
        Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior,
        truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);

    if (regression) { // update betas (regression)
      normalizer = exp(-h/2);
      Xnew = X;
      Xnew.each_col() %= normalizer;
      ynew = y % normalizer;

      // cholesky factor of posterior precision matrix
      postprecchol = arma::chol(Xnew.t() * Xnew + priorbetaprec);  // TODO handle exceptions the R way

      // inverse cholesky factor of posterior precision matrix 
      postpreccholinv = arma::solve(arma::trimatu(postprecchol), arma::eye<arma::mat>(p, p));

      // posterior covariance matrix and posterior mean vector
      postcov = postpreccholinv * postpreccholinv.t();
      postmean = postcov * (Xnew.t() * ynew + priorbetaprec * priorbetamean);

      armadraw = rnorm(p);

      // posterior betas
      curbeta = postmean + postpreccholinv * armadraw;
    }

    // storage:
    if (!((i+1) % thin)) if (i >= burnin) {  // this means we should store h
      if (centered_baseline) {
        for (int j = 0; j < hstorelength; j++) hstore.at(j, (i-burnin)/thin) = h[thintime * (j + 1) - 1];
        h0store[(i-burnin)/thin] = h0;
      } else {
        for (int j = 0; j < hstorelength; j++) hstore.at(j, (i-burnin)/thin) = curpara[0] + curpara[2]*h[thintime * (j + 1) - 1];
        h0store[(i-burnin)/thin] = curpara[0] + curpara[2]*h0;
      }
      if (keeptau && terr) {
        for (int j = 0; j < hstorelength; j++) taustore.at(j, (i-burnin)/thin) = tau[thintime * (j + 1) - 1];
      }
    }
    mu[i+1] = curpara[0];
    phi[i+1] = curpara[1];
    sigma2inv[i+1] = 1/(curpara[2]*curpara[2]);
    if (terr) nustore[i+1] = nu;
    if (regression) betastore.row(i+1) = curbeta.t();
  }  // END main MCMC loop

  if (verbose) progressbar_finish(N);  // finalize progress bar

  // Prepare return value and return
  return cleanUp(mu, phi, sqrt(1/sigma2inv), hstore, h0store, nustore, taustore, betastore);
}

List svlsample_cpp (
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X,
    const int thinpara,
    const int thinlatent,
    const int thintime,
    const List& theta_init,
    const arma::vec& h_init,
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
    const bool use_mala,
    const bool gammaprior,
    const bool correct,
    const CharacterVector& strategy_rcpp,
    const bool dontupdatemu) {

  const int N = burnin + draws;
  const bool regression = !ISNA(X.at(0,0));
  const int T = y_in.size();
  const int p = X.n_cols;

  arma::vec y = y_in;
  arma::vec y_star = log(y%y + offset);
  arma::ivec d(T); std::transform(y_in.cbegin(), y_in.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });

  double phi = theta_init["phi"];
  double rho = theta_init["rho"];
  double sigma2 = std::pow(as<double>(theta_init["sigma"]), 2);
  double mu = theta_init["mu"];
  arma::vec h = h_init, ht = (h_init-mu)/sqrt(sigma2);
  arma::vec beta(p); beta.fill(0.0);

  arma::mat betas(regression * draws/thinpara, p, arma::fill::zeros);
  arma::mat params(draws/thinpara, 4);

  const int hstorelength = T/thintime;  // thintime must either be 1 or T
  arma::mat latent(draws/thinlatent, hstorelength);

  // priors in objects
  const arma::vec prior_phi = {prior_phi_a, prior_phi_b};
  const arma::vec prior_rho = {prior_rho_a, prior_rho_b};
  const arma::vec prior_sigma2 = {prior_sigma2_shape, prior_sigma2_rate};
  const arma::vec prior_mu = {prior_mu_mu, prior_mu_sigma};

  // don't use strings or RcppCharacterVector
  arma::ivec strategy(strategy_rcpp.length());
  std::transform(strategy_rcpp.cbegin(), strategy_rcpp.cend(), strategy.begin(),
      [](const SEXP& par) -> int {
        if (as<std::string>(par) == "centered") return int(Parameterization::CENTERED);
        else if (as<std::string>(par) == "noncentered") return int(Parameterization::NONCENTERED);
        else Rf_error("Illegal parameterization");
  });

  // some stuff for the regression part
  // prior mean and precision matrix for the regression part (currently fixed)
  const arma::vec y_in_arma(y_in.begin(), T);
  const arma::vec priorbetamean = arma::ones(p) * prior_beta_mu;
  const arma::mat priorbetaprec = arma::eye(p, p) / std::pow(prior_beta_sigma, 2);
  arma::vec normalizer(T);
  arma::mat X_reg(T, p);
  arma::vec y_reg(T);
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);

  // adaptive MH
  stochvol::Adaptation adaptation(
      4,
      draws + burnin,
      use_mala ? 100 : 100,
      use_mala ? 0.574 : 0.234,
      0.05,
      use_mala ? 0.001 : 0.1);

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  const int show = verbose ? progressbar_init(N) : 0;

  for (int i = -burnin+1; i < draws+1; i++) {
    R_CheckUserInterrupt();

    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {
      y = y_in_arma - X*beta;
      y_star = arma::log(arma::square(y));
      std::transform(y.cbegin(), y.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });
    }

    // update theta and h
    update_svl(
        y, y_star, d,
        phi, rho, sigma2, mu,
        h, ht,
        adaptation,
        prior_phi, prior_rho,
        prior_sigma2, prior_mu,
        use_mala,
        gammaprior,
        correct,
        strategy,
        dontupdatemu);

    // update beta
    if (regression) {
      y_reg = y_in_arma;
      y_reg.head(T-1) -= rho * (arma::exp(h.head(T-1)/2) % (ht.tail(T-1) - phi*ht.head(T-1)));

      normalizer = arma::exp(-h/2);
      normalizer.head(T-1) /= std::sqrt(1 - std::pow(rho, 2));
      // X has already been copied to X_reg
      std::copy(X.cbegin(), X.cend(), X_reg.begin());  // important!
      X_reg.each_col() %= normalizer;
      y_reg %= normalizer;

      // cholesky factor of posterior precision matrix
      postprecchol = arma::chol(X_reg.t() * X_reg + priorbetaprec);

      // inverse cholesky factor of posterior precision matrix 
      postpreccholinv = arma::inv(arma::trimatu(postprecchol));

      // posterior covariance matrix and posterior mean vector
      postcov = postpreccholinv * postpreccholinv.t();
      postmean = postcov * (X_reg.t() * y_reg + priorbetaprec * priorbetamean);

      armadraw.imbue([]() -> double {return R::rnorm(0, 1);});  // equivalent to armadraw = rnorm(p); but I don't know if rnorm creates a vector

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
      for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
        latent.at(i/thinlatent-1, volind) = h[thintime * (thincol + 1) - 1];
      }
    }
  }

  if (verbose) progressbar_finish(N);  // finalize progress bar

  return List::create(
      _["para"] = params,
      _["adaptation"] = Rcpp::wrap(adaptation.get_storage()),
      _["latent"] = latent,
      _["beta"] = betas);
}

