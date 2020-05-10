#include <RcppArmadillo.h>
#include <cmath>
#include "theta-sampler.h"
#include "theta-utils.h"
#include "parameterization.h"

bool draw_theta(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const Parameterization centering,
    const stochvol::Adaptation<4>::Result& adaptation_proposal,
    const bool gammaprior,
    const Proposal sampler) {
  arma::vec proposed;
  switch (sampler) {
    case Proposal::RWMH:
      proposed = theta_propose_rwmh(phi, rho, sigma2, mu, y, h, ht, adaptation_proposal);
      break;
    case Proposal::MALA:
      proposed = theta_propose_mala(phi, rho, sigma2, mu, y, h, prior_phi, prior_rho, prior_sigma2, prior_mu, adaptation_proposal);
      break;
  }
  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
    mu_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];
  if (centering == Parameterization::NONCENTERED) {
    exp_h_half_proposal_nc = arma::exp(.5 * (std::sqrt(sigma2_prop) * ht + mu_prop));
  }
  const arma::vec& exp_h_half_proposal = centering == Parameterization::CENTERED ? exp_h_half : exp_h_half_proposal_nc;

  const double log_acceptance =
    (theta_log_prior(phi_prop, rho_prop, sigma2_prop, mu_prop, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
     theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu_prop, y, h, ht, exp_h_half_proposal, centering)) -
    (theta_log_prior(phi, rho, sigma2, mu, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
     theta_log_likelihood(phi, rho, sigma2, mu, y, h, ht, exp_h_half, centering)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 || std::exp(log_acceptance) > R::runif(0, 1);
  if (accepted) {
    phi = phi_prop;
    rho = rho_prop;
    sigma2 = sigma2_prop;
    mu = mu_prop;
  }

  return accepted;
}

// TODO reimplement
//void draw_thetamu_rwMH(
//    double& phi,
//    double& rho,
//    double& sigma2,
//    const double mu,
//    const arma::vec& y,
//    const arma::vec& h,
//    const arma::vec& ht,
//    const arma::vec& prior_phi,
//    const arma::vec& prior_rho,
//    const arma::vec& prior_sigma2,
//    const Parameterization centering,
//    const arma::mat& proposal_chol,
//    const arma::mat& proposal_chol_inv,
//    const bool gammaprior) {
//  const arma::vec proposed = thetamu_propose(phi, rho, sigma2, y, h, ht, proposal_chol, proposal_chol_inv);
//  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
//    prop_old_logdens = proposed[3], prop_new_logdens = proposed[4];
//  const double log_acceptance = (thetamu_log_prior(phi_prop, rho_prop, sigma2_prop, prior_phi, prior_rho, prior_sigma2, gammaprior) +
//    theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu, y, h, ht, centering)) -
//    (thetamu_log_prior(phi, rho, sigma2, prior_phi, prior_rho, prior_sigma2, gammaprior) +
//    theta_log_likelihood(phi, rho, sigma2, mu, y, h, ht, centering)) -
//    (prop_new_logdens - prop_old_logdens);
//  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
//    phi = phi_prop;
//    rho = rho_prop;
//    sigma2 = sigma2_prop;
//  }
//}

