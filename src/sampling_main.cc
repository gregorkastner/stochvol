#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "sampling_main.h"
#include "single_update.h"
#include "type_definitions.h"
#include "utils_main.h"

using namespace Rcpp;

namespace stochvol {

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
    const int thinpara,
    const int thinlatent,
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
    const double priordf,
    const arma::vec& priorbeta_in,
    const double priorlatent0) {

  arma::vec y(y_in);

  arma::mat X(X_in);

  const int T = y.size();
  const int p = X.n_cols;

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
  const bool terr = priordf > 0;

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
  arma::vec sigma2inv_store(draws / thinpara); sigma2inv_store[0] = std::pow(as<double>(startpara["sigma"]), -2);
  arma::vec phi_store(draws / thinpara); phi_store[0] = as<double>(startpara["phi"]);
  arma::vec mu_store(draws / thinpara); mu_store[0] = as<double>(startpara["mu"]);

  arma::vec h = startvol;  // contains h1 to hT, but not h0!
  if (!centered_baseline) h = (h-mu_store[0])*sqrt(sigma2inv_store[0]);

  double h0;

  arma::ivec r(T);  // mixture indicators

  int hstorelength = T/thintime;  // thintime must be either 1 or T
  arma::mat hstore(hstorelength, draws/thinlatent);
  arma::vec h0store(draws/thinlatent);

  arma::mat mixprob(10, T);  // mixture probabilities
  arma::vec mixprob_vec(mixprob.begin(), mixprob.n_elem, false);

  arma::vec data = log(y%y + offset);  // commonly used transformation

  arma::vec datastand = data;  // standardized "data" (different for t-errors)

  arma::vec curpara(3);  // holds mu, phi, sigma in every iteration
  curpara[0] = mu_store[0];
  curpara[1] = phi_store[0];
  curpara[2] = 1 / std::sqrt(sigma2inv_store[0]);

  // some stuff for the t-errors
  double nu = -1;
  if (terr) nu = as<double>(startpara["nu"]);
  arma::vec tau(T, arma::fill::ones);

  arma::vec nustore;
  if (terr) {
    nustore.resize(draws / thinpara);
    nustore[0] = nu;
  }
  arma::mat taustore;
  if (keeptau) {
    taustore.resize(hstorelength, draws / thinlatent);
  }

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
  arma::mat betastore;
  if (regression) {
    betastore.resize(draws / thinpara, p);
  }

  // Prior specification
  const PriorSpec prior_spec {
    (priorlatent0 <= 0) ? PriorSpec::Latent0() : PriorSpec::Latent0(PriorSpec::Constant(priorlatent0)),
    dontupdatemu ? PriorSpec::Mu(PriorSpec::Constant(0)) : PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),
    PriorSpec::Phi(PriorSpec::Beta(a0, b0)),
    Gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma)) : PriorSpec::Sigma2(PriorSpec::InverseGamma(2.5, C0)),
    priordf > 0 ? PriorSpec::Nu(PriorSpec::Exponential(priordf)) : PriorSpec::Nu(PriorSpec::Infinity())
  };
  // Expert (sampler specific) settings
  const ExpertSpec_VanillaSV expert {
    parameterization > 2,  // interweave
    parameterization % 2 ? Parameterization::CENTERED : Parameterization::NONCENTERED,  // centered_baseline
    B011inv,
    B022inv,
    MHsteps,
    MHcontrol < 0 ? ExpertSpec_VanillaSV::ProposalSigma2::INDEPENDENT : ExpertSpec_VanillaSV::ProposalSigma2::LOG_RANDOM_WALK,
    MHcontrol,
    truncnormal ? ExpertSpec_VanillaSV::ProposalPhi::TRUNCATED_NORMAL : ExpertSpec_VanillaSV::ProposalPhi::ACCEPT_REJECT_NORMAL
  };

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) show = progressbar_init(N);

  for (int i = -burnin + 1; i < draws + 1; i++) {  // BEGIN main MCMC loop
    R_CheckUserInterrupt();

    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {
      datastand = data = log(square(y - X*curbeta));
    }

    if (terr) {
      if (centered_baseline) {
        update_df_svt(data - h, tau, nu, prior_spec);
      } else {
        update_df_svt(data - curpara[0] - curpara[2] * h, tau, nu, prior_spec);
      }

      datastand = data - arma::log(tau);
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    double mu = curpara[0],
           phi = curpara[1],
           sigma2 = std::pow(curpara[2], 2);
    update_vanilla_sv(datastand, mu, phi, sigma2, h0, h, r, prior_spec, expert);
    curpara = {mu, phi, std::sqrt(sigma2)};

    if (regression) { // update betas (regression)
      normalizer = exp(-h/2);
      Xnew = X;
      Xnew.each_col() %= normalizer;
      ynew = y % normalizer;

      bool success = true;
      // cholesky factor of posterior precision matrix
      success = success && arma::chol(postprecchol, Xnew.t() * Xnew + priorbetaprec);
      // inverse cholesky factor of posterior precision matrix 
      success = success && arma::inv(postpreccholinv, arma::trimatu(postprecchol));
      if (!success) {
        Rcpp::stop("Cholesky or its inverse failed");
      }

      // posterior covariance matrix and posterior mean vector
      postcov = postpreccholinv * postpreccholinv.t();
      postmean = postcov * (Xnew.t() * ynew + priorbetaprec * priorbetamean);

      armadraw = rnorm(p);

      // posterior betas
      curbeta = postmean + postpreccholinv * armadraw;
    }

    // store draws
    if ((i >= 1) && !thinpara_round) {
      const int index = i / thinpara - 1;
      mu_store(index) = curpara[0];
      phi_store(index) = curpara[1];
      sigma2inv_store(index) = std::pow(curpara[2], -2);
      if (terr) {
        nustore[index] = nu;
      }
      if (regression) {
        betastore.row(index) = curbeta.t();
      }
    }
    if ((i >= 1) && !thinlatent_round) {
      const int index = i / thinlatent - 1;
      if (centered_baseline) {
        h0store[index] = h0;
        for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
          hstore.at(volind, index) = h[thintime * (thincol + 1) - 1];
        }
      } else {
        h0store[index] = curpara[0] + curpara[2]*h0;
        for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
          hstore.at(volind, index) = curpara[0] + curpara[2]*h[thintime * (thincol + 1) - 1];
        }
      }
      if (keeptau && terr) {
        for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
          taustore.at(volind, index) = tau[thintime * (thincol + 1) - 1];
        }
      }
    }
  }  // END main MCMC loop

  if (verbose) progressbar_finish(N);  // finalize progress bar

  // Prepare return value and return
  return cleanup(mu_store, phi_store, sqrt(1/sigma2inv_store), hstore, h0store, nustore, taustore, betastore);
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
  arma::vec y_star = arma::log(y%y + offset);
  arma::ivec d(T); std::transform(y_in.cbegin(), y_in.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });

  double phi = theta_init["phi"];
  double rho = theta_init["rho"];
  double sigma2 = std::pow(Rcpp::as<double>(theta_init["sigma"]), 2);
  double mu = theta_init["mu"];
  double h0 = std::numeric_limits<double>::quiet_NaN();
  arma::vec h = h_init, ht = (h_init-mu)/sqrt(sigma2);
  arma::vec beta(p); beta.fill(0.0);

  arma::mat betas(regression * draws/thinpara, p, arma::fill::zeros);
  arma::mat params(draws/thinpara, 4);

  const int hstorelength = T/thintime;  // thintime must be either 1 or T
  arma::vec latent0(draws/thinlatent);
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

  // Prior specification
  const PriorSpec prior_spec {
    PriorSpec::Latent0(),
    dontupdatemu ? PriorSpec::Mu(PriorSpec::Constant(0)) : PriorSpec::Mu(PriorSpec::Normal(prior_mu[0], prior_mu[1])),
    PriorSpec::Phi(PriorSpec::Beta(prior_phi[0], prior_phi[1])),
    gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(prior_sigma2[0], prior_sigma2[1])) : PriorSpec::Sigma2(PriorSpec::InverseGamma(prior_sigma2[0] + 2, prior_sigma2[1] / (prior_sigma2[0] * (prior_sigma2[0] + 1)))),  // moment matched inverse gamma
    PriorSpec::Nu(PriorSpec::Infinity()),
    PriorSpec::Rho(PriorSpec::Beta(prior_rho[0], prior_rho[1]))
  };
  // Expert (sampler specific) settings
  std::vector<Parameterization> strategy_vector(strategy.n_elem);
  std::transform(strategy.cbegin(), strategy.cend(), strategy_vector.begin(), [](const int ipar) -> Parameterization { return Parameterization(ipar); });
  const ExpertSpec_GeneralSV expert {
    strategy_vector,
    correct,
    use_mala
  };

  // adaptive MH
  AdaptationCollection adaptation_collection(
      4,
      (draws + burnin) / 200 + 1,
      200,
      use_mala ? 0.35 : 0.16,  //0.574 : 0.234,
      0.10,  // between 0 and 1: the larger the value the stronger and longer the adaptation
      0.001);

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
    update_general_sv(y, y_star, d, mu, phi, sigma2, rho, h0, h, adaptation_collection, prior_spec, expert);
    ht = (h - mu) / std::sqrt(sigma2);

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

      bool success = true;
      // cholesky factor of posterior precision matrix
      success = success && arma::chol(postprecchol, X_reg.t() * X_reg + priorbetaprec);
      // inverse cholesky factor of posterior precision matrix 
      success = success && arma::inv(postpreccholinv, arma::trimatu(postprecchol));
      if (!success) {
        Rcpp::stop("Cholesky or its inverse failed");
      }

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
      latent0[i/thinlatent-1] = h0;
      for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
        latent.at(i/thinlatent-1, volind) = h[thintime * (thincol + 1) - 1];
      }
    }
  }

  if (verbose) progressbar_finish(N);  // finalize progress bar

  return List::create(
      _["para"] = params,
      _["adaptation_centered"] = Rcpp::wrap(adaptation_collection.centered.get_storage()),
      _["adaptation_noncentered"] = Rcpp::wrap(adaptation_collection.noncentered.get_storage()),
      _["latent"] = latent,
      _["latent0"] = latent0,
      _["beta"] = betas);
}

}

