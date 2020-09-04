#include <RcppArmadillo.h>
#include "utils_main.h"
#include "densities.h"
#include <type_definitions.h>

using namespace Rcpp;

namespace stochvol {

List cleanup(
    const arma::vec& mu,
    const arma::vec& phi,
    const arma::vec& sigma,
    const arma::mat& h_store,
    const arma::vec& h0_store,
    const arma::vec& nu_store,
    const arma::mat& tau_store,
    const arma::mat& beta_store,
    const arma::imat& r_store) {
  int paracols;
  if (nu_store.size() > 0) paracols = 4; else paracols = 3;

  CharacterVector coln(paracols);
  NumericMatrix res(mu.n_elem, paracols);
  arma::mat res_arma(res.begin(), mu.n_elem, paracols, false, false); 
  res_arma.col(0) = mu; coln.at(0) = "mu";
  res_arma.col(1) = phi; coln.at(1) = "phi";
  res_arma.col(2) = sigma; coln.at(2) = "sigma";
  if (nu_store.size() > 0) {
    res_arma.col(3) = nu_store; coln.at(3) = "nu";
  }
  colnames(res) = coln;

  List val = List::create(
      _["para"] = res,
      _["latent"] = h_store,
      _["latent0"] = h0_store,
      _["beta"] = beta_store,
      _["tau"] = tau_store,
      _["indicators"] = r_store);

  return val;
}

int progressbar_init(
    const int N) {
  int show;
  REprintf("\n      ");
  if (N >= 2500) {
    for (int i = 0; i < 50+1; i++) REprintf(" ");
    show = N/50;
  }
  else {
    for (int i = 0; i < (N-1)/50+1; i++) REprintf(" ");
    show = 50;
  }
  REprintf("] 100%%\r  0%% [");
  R_FlushConsole();
  return show;
}

void progressbar_finish(
    const int N) {
  if (!(N % 50) && N >= 2500) REprintf("+");
  ::REprintf("] 100%%\n\n");
  ::R_FlushConsole();
}

PriorSpec list_to_priorspec(
    const Rcpp::List& list) {
  PriorSpec priorspec;
  const SEXP priorlatent0_sexp = list["latent0"];
  const List priormu = list["mu"],
             priorphi = list["phi"],
             priorsigma2 = list["sigma2"],
             priornu = list["nu"],
             priorrho = list["rho"],
             priorbeta = list["beta"];
  // latent0
  if (::Rf_isString(priorlatent0_sexp)) {
    const CharacterVector priorlatent0_rcpp = priorlatent0_sexp;
    const std::string priorlatent0 = as<std::string>(priorlatent0_rcpp);
    if (priorlatent0 == "stationary") {
      priorspec.latent0.variance = PriorSpec::Latent0::STATIONARY;
    } else {
      ::Rf_error("The prior specification of latent0 should be either the string \"stationary\" or an sv_constant object; got string \"%s\". See function specify_priors", priorlatent0.c_str());
    }
  } else if (::Rf_isList(priorlatent0_sexp)) {
    const List priorlatent0 = priorlatent0_sexp;
    if (priorlatent0.inherits("sv_constant")) {
      priorspec.latent0.variance = PriorSpec::Latent0::CONSTANT;
      priorspec.latent0.constant.value = as<double>(priorlatent0["value"]);
    } else {
      const CharacterVector classes_rcpp = priorlatent0.attr("class");
      const std::string classes = as<std::string>(classes_rcpp.at(0));
      ::Rf_error("The prior specification of latent0 should be either the string \"stationary\" or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
    }
  } else {
    ::Rf_error("The prior specification of latent0 should be either the string \"stationary\" or an sv_constant object; got type number %d. See function specify_priors", TYPEOF(priorlatent0_sexp));
  }
  // mu
  if (priormu.inherits("sv_normal")) {
    priorspec.mu.distribution = PriorSpec::Mu::NORMAL;
    priorspec.mu.normal.mean = as<double>(priormu["mean"]);
    priorspec.mu.normal.sd = as<double>(priormu["sd"]);
  } else if (priormu.inherits("sv_constant")) {
    priorspec.mu.distribution = PriorSpec::Mu::CONSTANT;
    priorspec.mu.constant.value = as<double>(priormu["value"]);
  } else {
    const CharacterVector classes_rcpp = priormu.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of mu should be either an sv_normal object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // phi
  if (priorphi.inherits("sv_beta")) {
    priorspec.phi.distribution = PriorSpec::Phi::BETA;
    priorspec.phi.beta.alpha = as<double>(priorphi["alpha"]);
    priorspec.phi.beta.beta = as<double>(priorphi["beta"]);
  } else if (priorphi.inherits("sv_constant")) {
    priorspec.phi.distribution = PriorSpec::Phi::CONSTANT;
    priorspec.phi.constant.value = as<double>(priorphi["value"]);
  } else {
    const CharacterVector classes_rcpp = priorphi.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of phi should be either an sv_beta object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // sigma2
  if (priorsigma2.inherits("sv_gamma")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::GAMMA;
    priorspec.sigma2.gamma.shape = as<double>(priorsigma2["shape"]);
    priorspec.sigma2.gamma.rate = as<double>(priorsigma2["rate"]);
  } else if (priorsigma2.inherits("sv_inverse_gamma")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::INVERSE_GAMMA;
    priorspec.sigma2.inverse_gamma.shape = as<double>(priorsigma2["shape"]);
    priorspec.sigma2.inverse_gamma.scale = as<double>(priorsigma2["scale"]);
  } else if (priorsigma2.inherits("sv_constant")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::CONSTANT;
    priorspec.sigma2.constant.value = as<double>(priorsigma2["value"]);
  } else {
    const CharacterVector classes_rcpp = priorsigma2.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of sigma2 should be an object of one of sv_gamma, sv_inverse_gamma, or sv_constant classes; got list with class %s. See function specify_priors", classes.c_str());
  }
  // nu
  if (priornu.inherits("sv_exponential")) {
    priorspec.nu.distribution = PriorSpec::Nu::EXPONENTIAL;
    priorspec.nu.exponential.rate = as<double>(priornu["rate"]);
  } else if (priornu.inherits("sv_constant")) {
    priorspec.nu.distribution = PriorSpec::Nu::CONSTANT;
    priorspec.nu.constant.value = as<double>(priornu["value"]);
  } else if (priornu.inherits("sv_infinity")) {
    priorspec.nu.distribution = PriorSpec::Nu::INIFINITY;
  } else {
    const CharacterVector classes_rcpp = priornu.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of nu should be an object of one of sv_exponential, sv_infinity, or sv_constant classes; got list with class %s. See function specify_priors", classes.c_str());
  }
  // rho
  if (priorrho.inherits("sv_beta")) {
    priorspec.rho.distribution = PriorSpec::Rho::BETA;
    priorspec.rho.beta.alpha = as<double>(priorrho["alpha"]);
    priorspec.rho.beta.beta = as<double>(priorrho["beta"]);
  } else if (priorrho.inherits("sv_constant")) {
    priorspec.rho.distribution = PriorSpec::Rho::CONSTANT;
    priorspec.rho.constant.value = as<double>(priorrho["value"]);
  } else {
    const CharacterVector classes_rcpp = priorrho.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of rho should be either an sv_beta object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // beta
  if (priorbeta.inherits("sv_normal")) {
    priorspec.beta.distribution = PriorSpec::Covariates::NORMAL;
    priorspec.beta.normal.mean = as<double>(priorbeta["mean"]);
    priorspec.beta.normal.sd = as<double>(priorbeta["sd"]);
  } else if (priorbeta.inherits("sv_constant")) {
    priorspec.beta.distribution = PriorSpec::Covariates::CONSTANT;
    priorspec.beta.constant.value = as<double>(priorbeta["value"]);
  } else {
    const CharacterVector classes_rcpp = priorbeta.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of beta should be either an sv_normal object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }

  return priorspec;
}

ExpertSpec_VanillaSV list_to_vanilla_sv(
    const Rcpp::List& list,
    const bool interweave) {
  const std::string baseline_parameterization_str = as<std::string>(list["baseline_parameterization"]),
                    proposal_phi_str = as<std::string>(list["proposal_phi"]),
                    proposal_sigma2_str = as<std::string>(list["proposal_sigma2"]);
  const double proposal_intercept_var = as<double>(list["proposal_intercept_var"]),
               proposal_phi_var = as<double>(list["proposal_phi_var"]),
               proposal_sigma2_rw_scale = as<double>(list["proposal_sigma2_rw_scale"]);
  const int mh_blocking_steps = as<int>(list["mh_blocking_steps"]);

  Parameterization baseline_parameterization;
  if (baseline_parameterization_str == "centered") {
    baseline_parameterization = Parameterization::CENTERED;
  } else if (baseline_parameterization_str == "noncentered") {
    baseline_parameterization = Parameterization::NONCENTERED;
  } else {
    ::Rf_error("Unknown value of baseline_parameterization in expert$fast_sv == \"%s\"; should be either \"centered\" or \"noncentered\"", baseline_parameterization_str.c_str());
  }

  ExpertSpec_VanillaSV::ProposalPhi proposal_phi;
  if (proposal_phi_str == "immediate acceptance-rejection") {
    proposal_phi = ExpertSpec_VanillaSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL;
  } else if (proposal_phi_str == "repeated acceptance-rejection") {
    proposal_phi = ExpertSpec_VanillaSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL;
  } else {
    ::Rf_error("Unknown value of proposal_phi in expert$fast_sv == \"%s\"; should be either \"immediate acceptance-rejection\" or \"repeated acceptance-rejection\"", proposal_phi_str.c_str());
  }

  ExpertSpec_VanillaSV::ProposalSigma2 proposal_sigma2;
  if (proposal_sigma2_str == "independence") {
    proposal_sigma2 = ExpertSpec_VanillaSV::ProposalSigma2::INDEPENDENCE;
  } else if (proposal_sigma2_str == "log random walk") {
    proposal_sigma2 = ExpertSpec_VanillaSV::ProposalSigma2::LOG_RANDOM_WALK;
  } else {
    ::Rf_error("Unknown value of proposal_sigma2 in expert$fast_sv == \"%s\"; should be either \"independence\" or \"log random walk\"", proposal_sigma2_str.c_str());
  }

  return {
    interweave,
    baseline_parameterization,
    1 / proposal_intercept_var,
    1 / proposal_phi_var,
    mh_blocking_steps,
    proposal_sigma2,
    proposal_sigma2_rw_scale,
    proposal_phi
  };
}

ExpertSpec_GeneralSV list_to_general_sv(
    const Rcpp::List& list,
    const bool correct_latent_draws,
    const bool interweave) {
  const int multi_asis = as<int>(list["multi_asis"]);
  const std::string starting_parameterization_str = as<std::string>(list["starting_parameterization"]);
  const std::string proposal_para_str = as<std::string>(list["proposal_para"]);
  const SEXP proposal_diffusion_ken_sexp = list["proposal_diffusion_ken"];

  // starting parameterization
  Parameterization starting_parameterization;
  if (starting_parameterization_str == "centered") {
    starting_parameterization = Parameterization::CENTERED;
  } else if (starting_parameterization_str == "noncentered") {
    starting_parameterization = Parameterization::NONCENTERED;
  } else {
    ::Rf_error("Unknown parameterization setting in expert$general_sv$starting_parameterization == \"%s\"; should be \"centered\" or \"noncentered\"", starting_parameterization_str.c_str());
  }
  const Parameterization other_parameterization = starting_parameterization == Parameterization::CENTERED ? Parameterization::NONCENTERED : Parameterization::CENTERED;

  // proposal strategy for the parameters
  ExpertSpec_GeneralSV::ProposalPara proposal_para;
  if (proposal_para_str == "random walk") {
    proposal_para = ExpertSpec_GeneralSV::ProposalPara::RANDOM_WALK;
  } else if (proposal_para_str == "metropolis-adjusted langevin algorithm") {
    proposal_para = ExpertSpec_GeneralSV::ProposalPara::METROPOLIS_ADJUSTED_LANGEVIN_ALGORITHM;
  } else {
    ::Rf_error("Unknown proposal setting in expert$general_sv$proposal_para == \"%s\"; should be \"random walk\" or \"metropolis-adjusted langevin algorithm\"", proposal_para_str.c_str());
  }

  // parameterization strategy
  std::vector<Parameterization> strategy;
  int strategy_size;
  if (interweave) {
    strategy_size = multi_asis * 2;
  } else {
    strategy_size = multi_asis;
  }
  strategy.reserve(strategy_size);
  for (int i = 0; i < multi_asis; i++) {
    strategy.push_back(starting_parameterization);
    if (interweave) {
      strategy.push_back(other_parameterization);
    }
  }

  // adaptation or fix proposal diffusion
  const bool adapt = ::Rf_isLogical(proposal_diffusion_ken_sexp);
  ProposalDiffusionKen proposal_diffusion_ken;
  if (!adapt) {
    const List proposal_diffusion_ken_rcpp = proposal_diffusion_ken_sexp;
    proposal_diffusion_ken.set(
        as<double>(proposal_diffusion_ken_rcpp["scale"]),
        as<arma::mat>(proposal_diffusion_ken_rcpp["covariance"]));
  }

  return {
    strategy,
    correct_latent_draws,
    proposal_para,
    adapt,
    proposal_diffusion_ken
  };
}

}

