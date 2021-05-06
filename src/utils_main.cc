/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
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
 * utils_main.cc
 *
 * Definitions of the functions declared in utils_main.h.
 * Documentation can also be found in utils_main.h.
 */

#include <RcppArmadillo.h>
#include <expert.hpp>
#include "utils_main.h"
#include "utils_latent_states.h"
#include "densities.h"

using namespace Rcpp;

namespace stochvol {

namespace fast_sv {

double compute_correction_weight(
    const arma::vec& data,
    const arma::vec& log_data2,
    const arma::vec& h,
    const arma::vec& exp_h_half) {
  static const arma::vec::fixed<10> log_mix_sd {arma::log(mix_sd)};
  const unsigned int T = data.n_elem;
  double exact_log_lik = 0,
         aux_log_lik = 0;
  for (unsigned int i = 0; i < T; i++) {
    exact_log_lik += logdnorm2(data[i], 0, exp_h_half[i], .5 * h[i]);
    double aux_lik_i = 0;
    for (unsigned int j = 0; j < mix_prob.n_elem; j++) {
      aux_lik_i += std::exp(logdnorm2(log_data2[i], h[i] + mix_mean[j], mix_sd[j], log_mix_sd[j])) * mix_prob[j];
    }
    aux_log_lik += std::log(aux_lik_i);
  }
  return exact_log_lik - aux_log_lik;
}

}  // END namespace fast_sv

void transpose_and_rename(
    const int T,
    NumericMatrix& para,
    NumericMatrix& latent0,
    NumericMatrix& latent,
    NumericMatrix& tau,
    NumericMatrix& betas) {
  para = Rcpp::transpose(para);
  latent = Rcpp::transpose(latent);
  tau = Rcpp::transpose(tau);
  betas = Rcpp::transpose(betas);

  {  // colnames in para
    const Rcpp::CharacterVector col_names {"mu", "phi", "sigma", "nu", "rho"};
    colnames(para) = col_names;
  }
  {  // colnames in latent0
    colnames(latent0) = Rcpp::CharacterVector({"h_0"});
  }
  {  // colnames in latent
    const unsigned int ncol = latent.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 1; c <= ncol; c++) {
      std::string name = "h_";
      name += std::to_string(T-ncol+c);
      col_names[c-1] = name;
    }
    colnames(latent) = col_names;
  }
  {  // colnames in betas
    const unsigned int ncol = betas.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 0; c < ncol; c++) {
      std::string name = "beta_";
      name += std::to_string(c);
      col_names[c] = name;
    }
    colnames(betas) = col_names;
  }
  {  // colnames in tau
    const unsigned int ncol = tau.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 1; c <= ncol; c++) {
      std::string name = "tau_";
      name += std::to_string(T-ncol+c);
      col_names[c-1] = name;
    }
    colnames(tau) = col_names;
  }
}

List cleanup(
    const int T,
    NumericMatrix& para,
    NumericMatrix& latent0,
    NumericMatrix& latent,
    NumericMatrix& tau,
    NumericMatrix& betas,
    IntegerMatrix& mixture_indicators,
    NumericVector& correction_weight_para,
    NumericVector& correction_weight_latent) {
  transpose_and_rename(T, para, latent0, latent, tau, betas);

  mixture_indicators = Rcpp::transpose(mixture_indicators);

  {  // colnames in mixture_indicators
    const unsigned int ncol = mixture_indicators.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 1; c <= ncol; c++) {
      std::string name = "r_";
      name += std::to_string(T-ncol+c);
      col_names[c-1] = name;
    }
    colnames(mixture_indicators) = col_names;
  }

  if (correction_weight_para.size() > 0) {
    correction_weight_para = Rcpp::exp(correction_weight_para - Rcpp::max(correction_weight_para));
    correction_weight_para = correction_weight_para / Rcpp::sum(correction_weight_para);
  }
  if (correction_weight_latent.size() > 0) {
    correction_weight_latent = Rcpp::exp(correction_weight_latent - Rcpp::max(correction_weight_latent));
    correction_weight_latent = correction_weight_latent / Rcpp::sum(correction_weight_latent);
  }

  List val = List::create(
      _["para"] = para,
      _["latent"] = latent,
      _["latent0"] = latent0,
      _["beta"] = betas,
      _["tau"] = tau,
      _["indicators"] = mixture_indicators + 1u,
      _["correction_weight_para"] = correction_weight_para,
      _["correction_weight_latent"] = correction_weight_latent);

  return val;
}

List cleanup(
    const int T,
    NumericMatrix& para,
    NumericMatrix& latent0,
    NumericMatrix& latent,
    NumericMatrix& tau,
    NumericMatrix& betas,
    AdaptationCollection& adaptation_collection) {
  transpose_and_rename(T, para, latent0, latent, tau, betas);

  return List::create(
      _["para"] = para,
      _["adaptation"] = adaptation_collection.serialize(),
      _["latent"] = latent,
      _["latent0"] = latent0,
      _["tau"] = tau,
      _["beta"] = betas);
}

int progressbar_init(
    const int N) {
  int show;
  ::REprintf("\n      ");
  if (N >= 2500) {
    for (int i = 0; i < 50+1; i++) ::REprintf(" ");
    show = N/50;
  }
  else {
    for (int i = 0; i < (N-1)/50+1; i++) ::REprintf(" ");
    show = 50;
  }
  ::REprintf("] 100%%\r  0%% [");
  ::R_FlushConsole();
  return show;
}

void progressbar_finish(
    const int N) {
  if (!(N % 50) && N >= 2500) ::REprintf("+");
  ::REprintf("] 100%%\n\n");
  ::R_FlushConsole();
}

PriorSpec list_to_priorspec(
    const Rcpp::List& list) {
  PriorSpec priorspec;
  const SEXP priorlatent0_sexp = list["latent0_variance"];
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
      ::Rf_error("The prior specification for the variance of latent0 should be either the string \"stationary\" or an sv_constant object; got string \"%s\". See function specify_priors", priorlatent0.c_str());
    }
  } else if (::Rf_isVectorList(priorlatent0_sexp)) {
    const List priorlatent0 = priorlatent0_sexp;
    if (priorlatent0.inherits("sv_constant")) {
      priorspec.latent0.variance = PriorSpec::Latent0::CONSTANT;
      priorspec.latent0.constant.value = priorlatent0["value"];
    } else {
      const CharacterVector classes_rcpp = priorlatent0.attr("class");
      const std::string classes = as<std::string>(classes_rcpp.at(0));
      ::Rf_error("The prior specification for the variance of latent0 should be either the string \"stationary\" or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
    }
  } else {
    ::Rf_error("The prior specification for the variance of latent0 should be either the string \"stationary\" or an sv_constant object; got type number %d. See function specify_priors", TYPEOF(priorlatent0_sexp));
  }
  // mu
  if (priormu.inherits("sv_normal")) {
    priorspec.mu.distribution = PriorSpec::Mu::NORMAL;
    priorspec.mu.normal.mean = priormu["mean"];
    priorspec.mu.normal.sd = priormu["sd"];
  } else if (priormu.inherits("sv_constant")) {
    priorspec.mu.distribution = PriorSpec::Mu::CONSTANT;
    priorspec.mu.constant.value = priormu["value"];
  } else {
    const CharacterVector classes_rcpp = priormu.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of mu should be either an sv_normal object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // phi
  if (priorphi.inherits("sv_beta")) {
    priorspec.phi.distribution = PriorSpec::Phi::BETA;
    priorspec.phi.beta.alpha = priorphi["shape1"];
    priorspec.phi.beta.beta = priorphi["shape2"];
  } else if (priorphi.inherits("sv_normal")) {
    priorspec.phi.distribution = PriorSpec::Phi::NORMAL;
    priorspec.phi.normal.mean = priorphi["mean"];
    priorspec.phi.normal.sd = priorphi["sd"];
  } else if (priorphi.inherits("sv_constant")) {
    priorspec.phi.distribution = PriorSpec::Phi::CONSTANT;
    priorspec.phi.constant.value = priorphi["value"];
  } else {
    const CharacterVector classes_rcpp = priorphi.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of phi should be either an sv_beta object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // sigma2
  if (priorsigma2.inherits("sv_gamma")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::GAMMA;
    priorspec.sigma2.gamma.shape = priorsigma2["shape"];
    priorspec.sigma2.gamma.rate = priorsigma2["rate"];
  } else if (priorsigma2.inherits("sv_inverse_gamma")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::INVERSE_GAMMA;
    priorspec.sigma2.inverse_gamma.shape = priorsigma2["shape"];
    priorspec.sigma2.inverse_gamma.scale = priorsigma2["scale"];
  } else if (priorsigma2.inherits("sv_constant")) {
    priorspec.sigma2.distribution = PriorSpec::Sigma2::CONSTANT;
    priorspec.sigma2.constant.value = priorsigma2["value"];
  } else {
    const CharacterVector classes_rcpp = priorsigma2.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of sigma2 should be an object of one of sv_gamma, sv_inverse_gamma, or sv_constant classes; got list with class %s. See function specify_priors", classes.c_str());
  }
  // nu
  if (priornu.inherits("sv_exponential")) {
    priorspec.nu.distribution = PriorSpec::Nu::EXPONENTIAL;
    priorspec.nu.exponential.rate = priornu["rate"];
  } else if (priornu.inherits("sv_constant")) {
    priorspec.nu.distribution = PriorSpec::Nu::CONSTANT;
    priorspec.nu.constant.value = priornu["value"];
  } else if (priornu.inherits("sv_infinity")) {
    priorspec.nu.distribution = PriorSpec::Nu::INFINITE;
  } else {
    const CharacterVector classes_rcpp = priornu.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of nu should be an object of one of sv_exponential, sv_infinity, or sv_constant classes; got list with class %s. See function specify_priors", classes.c_str());
  }
  // rho
  if (priorrho.inherits("sv_beta")) {
    priorspec.rho.distribution = PriorSpec::Rho::BETA;
    priorspec.rho.beta.alpha = priorrho["shape1"];
    priorspec.rho.beta.beta = priorrho["shape2"];
  } else if (priorrho.inherits("sv_constant")) {
    priorspec.rho.distribution = PriorSpec::Rho::CONSTANT;
    priorspec.rho.constant.value = priorrho["value"];
  } else {
    const CharacterVector classes_rcpp = priorrho.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of rho should be either an sv_beta object or an sv_constant object; got list with class %s. See function specify_priors", classes.c_str());
  }
  // beta
  if (priorbeta.inherits("sv_multinormal")) {
    try {
      priorspec.beta.multivariate_normal.mean = as<arma::vec>(priorbeta["mean"]);
      priorspec.beta.multivariate_normal.precision = as<arma::mat>(priorbeta["precision"]);
    } catch (...) {
      //Rcout << "Received prior specification for beta:" << std::endl << priorbeta << std::endl;
      ::Rf_error("Unable to convert priorspec$priorbeta to a mean vector and a precision matrix");
    }
    if (!priorspec.beta.multivariate_normal.precision.is_sympd()) {
      //Rcout << "Received precision matrix as the prior specification for beta:" << std::endl << priorspec.beta.multivariate_normal.precision << std::endl;
      ::Rf_error("The precision matrix of the prior specification for beta is not symmetric and/or positive definite.");
    }
  } else {
    const CharacterVector classes_rcpp = priorbeta.attr("class");
    const std::string classes = as<std::string>(classes_rcpp.at(0));
    ::Rf_error("The prior specification of beta should be an sv_multinormal object; got list with class %s. See function specify_priors", classes.c_str());
  }

  return priorspec;
}

ExpertSpec_FastSV list_to_fast_sv(
    const Rcpp::List& list,
    const bool interweave) {
  const std::string baseline_parameterization_str = as<std::string>(list["baseline_parameterization"]),
                    proposal_phi_str = as<std::string>(list["proposal_phi"]),
                    proposal_sigma2_str = as<std::string>(list["proposal_sigma2"]);
  const double proposal_intercept_var = list["proposal_intercept_var"],
               proposal_phi_var = list["proposal_phi_var"],
               proposal_sigma2_rw_scale = list["proposal_sigma2_rw_scale"];
  const int mh_blocking_steps = as<int>(list["mh_blocking_steps"]);
  const Rcpp::List update_list = list["update"];

  Parameterization baseline_parameterization;
  if (baseline_parameterization_str == "centered") {
    baseline_parameterization = Parameterization::CENTERED;
  } else if (baseline_parameterization_str == "noncentered") {
    baseline_parameterization = Parameterization::NONCENTERED;
  } else {
    ::Rf_error("Unknown value of baseline_parameterization in expert$fast_sv == \"%s\"; should be either \"centered\" or \"noncentered\"", baseline_parameterization_str.c_str());
  }

  ExpertSpec_FastSV::ProposalPhi proposal_phi;
  if (proposal_phi_str == "immediate acceptance-rejection") {
    proposal_phi = ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL;
  } else if (proposal_phi_str == "repeated acceptance-rejection") {
    proposal_phi = ExpertSpec_FastSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL;
  } else {
    ::Rf_error("Unknown value of proposal_phi in expert$fast_sv == \"%s\"; should be either \"immediate acceptance-rejection\" or \"repeated acceptance-rejection\"", proposal_phi_str.c_str());
  }

  ExpertSpec_FastSV::ProposalSigma2 proposal_sigma2;
  if (proposal_sigma2_str == "independence") {
    proposal_sigma2 = ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE;
  } else if (proposal_sigma2_str == "log random walk") {
    proposal_sigma2 = ExpertSpec_FastSV::ProposalSigma2::LOG_RANDOM_WALK;
  } else {
    ::Rf_error("Unknown value of proposal_sigma2 in expert$fast_sv == \"%s\"; should be either \"independence\" or \"log random walk\"", proposal_sigma2_str.c_str());
  }

  ExpertSpec_FastSV::Update update;
  update.latent_vector = update_list["latent_vector"];
  update.mixture_indicators = update_list["mixture_indicators"];
  update.parameters = update_list["parameters"];

  return {
    interweave,
    baseline_parameterization,
    1 / proposal_intercept_var,
    1 / proposal_phi_var,
    mh_blocking_steps,
    proposal_sigma2,
    proposal_sigma2_rw_scale,
    proposal_phi,
    update
  };
}

ExpertSpec_GeneralSV list_to_general_sv(
    const Rcpp::List& list,
    const bool correct_latent_draws,
    const bool interweave) {
  const int multi_asis = as<int>(list["multi_asis"]);
  const std::string starting_parameterization_str = as<std::string>(list["starting_parameterization"]);
  const SEXP proposal_diffusion_ken_sexp = list["proposal_diffusion_ken"];
  const Rcpp::List update_list = list["update"];

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
  ExpertSpec_GeneralSV::ProposalPara proposal_para = ExpertSpec_GeneralSV::ProposalPara::RANDOM_WALK;

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
        proposal_diffusion_ken_rcpp["scale"],
        as<arma::mat>(proposal_diffusion_ken_rcpp["covariance"]));
  }

  ExpertSpec_GeneralSV::Update update;
  update.latent_vector = update_list["latent_vector"];
  update.parameters = update_list["parameters"];

  return {
    strategy,
    correct_latent_draws,
    proposal_para,
    adapt,
    proposal_diffusion_ken,
    update
  };
}

/*
ProposalDiffusionKen list_to_proposal_ken (
    const List& list) {
  return ProposalDiffusionKen(list["scale"], list["covariance"]);
}

List proposal_ken_to_list (
    const ProposalDiffusionKen& ken) {
  return List::create(
      _["scale"] = wrap(ken.get_scale()),
      _["covariance"] = wrap(ken.get_covariance()));
}
*/

AdaptationCollection list_to_adaptationcollection (
    const List& list) {
  return {
    list_to_adaptation(list["centered"]),
    list_to_adaptation(list["noncentered"])};
}

Adaptation list_to_adaptation (
    const List& list) {
  using Memory = std::vector<Adaptation::Storage>;
  const NumericMatrix memory_rcpp = list["memory"];
  Memory memory;
  memory.reserve(memory_rcpp.nrow());
  for (int i = 0; i < memory_rcpp.nrow() and not std::isnan(memory_rcpp(i, 0)); i++) {
    memory.push_back({memory_rcpp(i, 0), memory_rcpp(i, 1), memory_rcpp(i, 2)});
  }

  const NumericVector mu_rcpp (as<NumericVector>(list["mu"]));
  const NumericMatrix Sigma_rcpp (as<NumericMatrix>(list["Sigma"]));
  const NumericMatrix draws_batch_rcpp (as<NumericMatrix>(list["draws_batch"]));
  const NumericMatrix cached_covariance_rcpp (as<NumericMatrix>(list["cached_covariance"]));

  const arma::vec mu (mu_rcpp.cbegin(), mu_rcpp.length());
  const arma::mat Sigma (Sigma_rcpp.cbegin(), Sigma_rcpp.nrow(), Sigma_rcpp.ncol());
  const arma::mat draws_batch (draws_batch_rcpp.cbegin(), draws_batch_rcpp.nrow(), draws_batch_rcpp.ncol());
  const arma::mat cached_covariance (cached_covariance_rcpp.cbegin(), cached_covariance_rcpp.nrow(), cached_covariance_rcpp.ncol());

  return {
    as<int>(list["dim"]),
    memory,  //std::move(memory),
    as<int>(list["batch_size"]),
    as<double>(list["target_acceptance"]),
    as<double>(list["lambda"]),
    as<double>(list["scale"]),
    as<double>(list["C"]),
    as<double>(list["alpha"]),
    as<double>(list["gamma"]),
    as<int>(list["count_acceptance"]),
    as<int>(list["i_batch"]),
    mu,  //std::move(mu),
    Sigma,  //std::move(Sigma),
    draws_batch,  //std::move(draws_batch),
    as<bool>(list["updated_proposal"]),
    as<double>(list["cached_scale"]),
    cached_covariance};
}

}

