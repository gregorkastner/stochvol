#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "sampling_main.h"
#include "single_update.h"
#include "type_definitions.h"
#include "utils_main.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

List svsample_cpp(
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keep_tau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {

  arma::vec y(y_in);

  const int T = y.size();
  const int p = X.n_cols;

  // should we model the mean as well?
  const bool regression = !::ISNA(X.at(0,0));
  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_VanillaSV expert = list_to_vanilla_sv(expert_in, interweave);

  if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT && expert.mh_blocking_steps == 1) { // not implemented (would be easy, though)
    ::Rf_error("Single block update leaving mu constant is not yet implemented");
  }

  // shortcuts / precomputations
  const int thintime = ([&keeptime_in, T] () -> int {
      const std::string keeptime = as<std::string>(keeptime_in);
      if (keeptime == "all") {
        return 1;
      } else if (keeptime == "last") {
        return T;
      } else {
        ::Rf_error("Unknown value for 'keeptime'; got \"%s\"", keeptime.c_str());
      }
  })();

  // prior mean and precision matrix for the regression part (currently fixed)
  const arma::vec priorbetamean = arma::ones(p) * prior_spec.beta.normal.mean;
  const arma::mat priorbetaprec = arma::eye(p, p) / std::pow(prior_spec.beta.normal.sd, 2);

  // number of MCMC draws
  const int N = burnin + draws;

  // verbosity control
  const bool verbose = !quiet;

  // t-errors:
  const bool terr = prior_spec.nu.distribution == PriorSpec::Nu::EXPONENTIAL;

  // initialize the variables:
  double mu = startpara["mu"],
         phi = startpara["phi"],
         sigma2 = std::pow(as<double>(startpara["sigma"]), 2);
  arma::vec sigma2inv_store(draws / thinpara),
            phi_store(draws / thinpara),
            mu_store(draws / thinpara);
  sigma2inv_store[0] = 1. / sigma2;
  phi_store[0] = phi;
  mu_store[0] = mu;

  arma::vec h = startlatent;  // contains h1 to hT, but not h0!
  double h0 = startpara["latent0"];

  const bool keep_r = expert_in["store_indicators"];
  arma::ivec r = expert_in["init_indicators"];  // mixture indicators
  if (r.n_elem == 1) {
    const double r_elem = r[0];
    r.set_size(T);
    r.fill(r_elem);
  } else if (r.n_elem != T) {
    ::Rf_error("Bad initialization for the vector of mixture indicators. Should have length %d, received length %d, first element %f", T, r.n_elem, r[0]);
  }
  arma::imat r_store;
  if (keep_r) {
    r_store.set_size(T / thintime, draws / thinlatent);
  }

  const int hstorelength = T / thintime;  // thintime must be either 1 or T
  arma::mat h_store(hstorelength, draws / thinlatent);
  arma::vec h0_store(draws / thinlatent);

  arma::vec data = arma::log(y%y + offset);  // commonly used transformation

  arma::vec datastand = data;  // standardized "data" (different for t-errors)

  // some stuff for the t-errors
  double nu = -1;
  if (terr) nu = startpara["nu"];
  arma::vec tau(T, arma::fill::ones);

  arma::vec nu_store;
  if (terr) {
    nu_store.set_size(draws / thinpara);
    nu_store[0] = nu;
  }
  arma::mat tau_store;
  if (keep_tau) {
    tau_store.set_size(hstorelength, draws / thinlatent);
  }

  // some stuff for the regression part
  arma::vec beta = startpara["beta"];
  arma::vec normalizer;
  arma::mat Xnew = X;
  arma::vec ynew = y;
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);
  arma::mat beta_store;
  if (regression) {
    beta_store.resize(draws / thinpara, p);
  }

  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  int show = 0;
  if (verbose) show = progressbar_init(N);

  for (int i = -burnin + 1; i < draws + 1; i++) {  // BEGIN main MCMC loop
    //Rcpp::Rcout << "round " << i << std::endl;
    R_CheckUserInterrupt();

    const bool thinpara_round = (thinpara > 1) && (i % thinpara != 0);  // is this a parameter thinning round?
    const bool thinlatent_round = (thinlatent > 1) && (i % thinlatent != 0);  // is this a latent thinning round?

    // print a progress sign every "show" iterations
    if (verbose && (i % show == 0)) progressbar_print();

    if (regression) {
      datastand = data = arma::log(arma::square(y - X*beta));
    }

    if (terr) {
      update_df_svt(data - h, tau, nu, prior_spec);
      datastand = data - arma::log(tau);
    }

    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_vanilla_sv(datastand, mu, phi, sigma2, h0, h, r, prior_spec, expert);

    if (regression) { // update betas (regression)
      normalizer = arma::exp(-h/2);
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
      beta = postmean + postpreccholinv * armadraw;
    }

    // store draws
    if ((i >= 1) && !thinpara_round) {
      const int index = i / thinpara - 1;
      mu_store(index) = mu;
      phi_store(index) = phi;
      sigma2inv_store(index) = 1. / sigma2;
      if (terr) {
        nu_store[index] = nu;
      }
      if (regression) {
        beta_store.row(index) = beta.t();
      }
    }
    if ((i >= 1) && !thinlatent_round) {
      const int index = i / thinlatent - 1;
      h0_store[index] = h0;
      for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
        h_store.at(volind, index) = h[thintime * (thincol + 1) - 1];
      }
      if (keep_tau && terr) {
        for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
          tau_store.at(volind, index) = tau[thintime * (thincol + 1) - 1];
        }
      }
      if (keep_r) {
        for (int volind = 0, thincol = 0; thincol < hstorelength; volind++, thincol++) {
          r_store.at(volind, index) = r[thintime * (thincol + 1) - 1];
        }
      }
    }
  }  // END main MCMC loop

  if (verbose) progressbar_finish(N);  // finalize progress bar

  // Prepare return value and return
  return cleanup(mu_store, phi_store, sqrt(1/sigma2inv_store), h_store, h0_store, nu_store, tau_store, beta_store, r_store);
}

List svlsample_cpp (
    const arma::vec& y_in,
    const int draws,
    const int burnin,
    const arma::mat& X,
    const Rcpp::List& priorspec_in,
    const int thinpara,
    const int thinlatent,
    const Rcpp::CharacterVector& keeptime_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const bool keeptau,
    const bool quiet,
    const bool correct_model_specification,
    const bool interweave,
    const double offset,
    const Rcpp::List& expert_in) {

  const int N = burnin + draws;
  const bool regression = !ISNA(X.at(0,0));
  const int T = y_in.size();
  const int p = X.n_cols;

  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  const ExpertSpec_GeneralSV expert = list_to_general_sv(expert_in, correct_model_specification, interweave);

  const bool verbose = !quiet;
  const int thintime = ([&keeptime_in, T] () -> int {
      const std::string keeptime = as<std::string>(keeptime_in);
      if (keeptime == "all") {
        return 1;
      } else if (keeptime == "last") {
        return T;
      } else {
        Rf_error("Unknown value for 'keeptime'; got \"%s\"", keeptime.c_str());
      }
  })();

  arma::vec y = y_in;
  arma::vec y_star = arma::log(y%y + offset);
  arma::ivec d(T); std::transform(y_in.cbegin(), y_in.cend(), d.begin(), [](const double y_elem) -> int { return y_elem > 0 ? 1 : -1; });

  double mu = startpara["mu"];
  double phi = startpara["phi"];
  double sigma2 = std::pow(as<double>(startpara["sigma"]), 2);
  double rho = startpara["rho"];
  double h0 = startpara["latent0"];
  arma::vec h = startlatent, ht = centered_to_noncentered(mu, std::sqrt(sigma2), h);
  arma::vec beta = startpara["beta"];

  arma::mat betas(regression * draws/thinpara, p, arma::fill::zeros);
  Rcpp::NumericMatrix para(draws/thinpara, 4);
  arma::mat params(para.begin(), para.nrow(), para.ncol(), false, false);

  const int hstorelength = T/thintime;  // thintime must be either 1 or T
  arma::vec latent0(draws/thinlatent);
  arma::mat latent(draws/thinlatent, hstorelength);

  // some stuff for the regression part
  // prior mean and precision matrix for the regression part (currently fixed)
  const arma::vec y_in_arma(y_in.begin(), T);
  const arma::vec priorbetamean = arma::ones(p) * prior_spec.beta.normal.mean;
  const arma::mat priorbetaprec = arma::eye(p, p) / std::pow(prior_spec.beta.normal.sd, 2);
  arma::vec normalizer(T);
  arma::mat X_reg(T, p);
  arma::vec y_reg(T);
  arma::mat postprecchol(p, p);
  arma::mat postpreccholinv(p, p);
  arma::mat postcov(p, p);
  arma::vec postmean(p);
  arma::vec armadraw(p);

  // adaptive MH
  const int batch_size = 200,
            memory_size = expert.adapt ? expert.strategy.size() * (draws + burnin) / batch_size + 1 : 1;
  const double target_acceptance = expert.proposal_para == ExpertSpec_GeneralSV::ProposalPara::METROPOLIS_ADJUSTED_LANGEVIN_ALGORITHM ? 0.574 : 0.234,  //0.35 : 0.16,
               lambda = 0.1,
               init_scale = 0.001;
  AdaptationCollection adaptation_collection(
      4,
      memory_size,
      batch_size,
      target_acceptance,
      lambda,  // between 0 and 1: the larger the value the stronger and longer the adaptation
      init_scale);

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
    ht = centered_to_noncentered(mu, std::sqrt(sigma2), h);

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
        Rcpp::stop("Cholesky or its inverse failed during the sampling of the betas");
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

  Rcpp::CharacterVector coln(para.ncol());
  coln[0] = "mu"; coln[1] = "phi"; coln[2] = "sigma"; coln[3] = "rho";
  colnames(para) = coln;

  return List::create(
      _["para"] = para,
      _["adaptation"] = List::create(
        _["centered"] = List::create(
          _["history"] = adaptation_collection.centered.get_storage(),
          _["scale"] = wrap(adaptation_collection.centered.get_proposal().get_scale()),
          _["covariance"] = wrap(adaptation_collection.centered.get_proposal().get_covariance())),
        _["noncentered"] = List::create(
          _["history"] = adaptation_collection.noncentered.get_storage(),
          _["scale"] = wrap(adaptation_collection.noncentered.get_proposal().get_scale()),
          _["covariance"] = wrap(adaptation_collection.noncentered.get_proposal().get_covariance()))),
      _["latent"] = latent,
      _["latent0"] = latent0,
      _["beta"] = betas);
}

}

