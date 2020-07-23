#ifndef _SINGLE_UPDATES_H_
#define _SINGLE_UPDATES_H_

#include <RcppArmadillo.h>
#include <adaptation.hpp>
#include "type_definitions.h"

namespace stochvol {

// single MCMC update approx. SV
void update_vanilla_sv(
    const arma::vec& log_data2,
    double& mu,
    double& phi,
    double& sigma2,
    double& h0,
    arma::vec& h,
    arma::ivec& r,
    const PriorSpec& prior_spec,  // C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, truncnormal, dontupdatemu, priorlatent0 feed into this
    const ExpertSpec_VanillaSV& expert);  // parameterization, centered_baseline, B011inv, B022inv, Mhcontrol, MHsteps feed into this

// single MCMC update degrees of freedom SVt
void update_df_svt(
    const arma::vec& log_data2_minus_h,
    arma::vec& tau,
    double& nu,
    const PriorSpec& prior_spec);

// single MCMC update SVl
void update_general_sv(
    const arma::vec& data,
    const arma::vec& log_data2,
    const arma::ivec& sign_data,
    double& mu,
    double& phi,
    double& sigma2,
    double& rho,
    double& h0,
    arma::vec& h,
    AdaptationCollection& adaptation,
    const PriorSpec& prior_spec,  // prior_mu, prior_phi, prior_sigma2, prior_rho, gammaprior, dontupdatemu feed into this (plus priorlatent0, truncnormal nyi)
    const ExpertSpec_GeneralSV& expert);  // strategy, correct, use_mala feed into this

// OLD FUNCTIONS

inline
void update_sv [[gnu::deprecated]] (
    const arma::vec& data,
    arma::vec& curpara,
    arma::vec& h,
    double& h0,
    arma::vec& mixprob,
    arma::ivec& r,
    const bool centered_baseline,
    const double C0,
    const double cT,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const double B011inv,
    const double B022inv,
    const bool Gammaprior,
    const bool truncnormal,
    const double MHcontrol,
    const int MHsteps,
    const int parameterization,
    const bool dontupdatemu,
    const double priorlatent0) {
  double mu = curpara[0],
         phi = curpara[1],
         sigma2 = std::pow(curpara[2], 2);
  const PriorSpec prior_spec {
    (priorlatent0 <= 0) ? PriorSpec::Latent0() : PriorSpec::Latent0(PriorSpec::Constant(priorlatent0)),
    dontupdatemu ? PriorSpec::Mu(PriorSpec::Constant(0)) : PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),
    PriorSpec::Phi(PriorSpec::Beta(a0, b0)),
    Gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma)) : PriorSpec::Sigma2(PriorSpec::InverseGamma(2.5, C0))
  };
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
  update_vanilla_sv(data, mu, phi, sigma2, h0, h, r, prior_spec, expert);
  curpara = {mu, phi, std::sqrt(sigma2)};
}

inline
void update_terr [[gnu::deprecated]] (
    const arma::vec& data,
    arma::vec& tau,
    double& nu,
    const double lambda) {
  update_df_svt(data, tau, nu, {{}, PriorSpec::Constant(0), PriorSpec::Constant(0), PriorSpec::Constant(0), PriorSpec::Exponential(lambda), PriorSpec::Constant(0)});
}

inline
void update_svl [[gnu::deprecated]] (
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    double& h0,
    arma::vec& h,
    arma::vec& ht,
    AdaptationCollection& adaptation,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const bool use_mala,
    const bool gammaprior,
    const bool correct,
    const arma::ivec& strategy,
    const bool dontupdatemu) {
  const PriorSpec prior_spec {
    PriorSpec::Latent0(),
    dontupdatemu ? PriorSpec::Mu(PriorSpec::Constant(0)) : PriorSpec::Mu(PriorSpec::Normal(prior_mu[0], prior_mu[1])),
    PriorSpec::Phi(PriorSpec::Beta(prior_phi[0], prior_phi[1])),
    gammaprior ? PriorSpec::Sigma2(PriorSpec::Gamma(prior_sigma2[0], prior_sigma2[1])) : PriorSpec::Sigma2(PriorSpec::InverseGamma(prior_sigma2[0] + 2, prior_sigma2[1] / (prior_sigma2[0] * (prior_sigma2[0] + 1)))),  // moment matched inverse gamma
    PriorSpec::Nu(PriorSpec::Infinity()),
    PriorSpec::Rho(PriorSpec::Beta(prior_rho[0], prior_rho[1]))
  };
  std::vector<Parameterization> strategy_vector(strategy.n_elem);
  std::transform(strategy.cbegin(), strategy.cend(), strategy_vector.begin(), [](const int ipar) -> Parameterization { return Parameterization(ipar); });
  const ExpertSpec_GeneralSV expert {
    strategy_vector,
    correct,
    use_mala
  };
  update_general_sv(y, y_star, d, mu, phi, sigma2, rho, h0, h, adaptation, prior_spec, expert);
  ht = (h - mu) / std::sqrt(sigma2);
}

}

#endif  // _SINGLE_UPDATES_H_
