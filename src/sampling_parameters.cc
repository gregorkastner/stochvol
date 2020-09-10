#include <RcppArmadillo.h>
#include <cmath>
#include "sampling_parameters.h"
#include "utils_parameters.h"
#include "utils_latent_states.h"
#include "type_definitions.h"
#include "densities.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

ReturnRegression regression_centered(
    const double h0,
    const arma::vec& h,
    double mu,
    double phi,
    double sigma,
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
    const bool dontupdatemu,
    const double priorlatent0) {

  const double Bh0inv = (priorlatent0 < 0) ? 1-phi*phi : 1./priorlatent0;

  if (dontupdatemu) mu = 0;  // TODO change

  const int T = h.size();
  double sigma_prop;

  // first calculate bT and BT:
  double sum1 = h[0];
  double sum3 = h0*h[0];
  double sum4 = h[0]*h[0];
  for (int j = 1; j < T-1; j++) {
    sum1 += h[j];
    sum3 += h[j-1]*h[j];
    sum4 += h[j]*h[j];
  }
  double sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
  sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
  sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
  sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2

  const double tmp1 = 1/(((sum4 + B011inv)*(T+B022inv)-sum1*sum1)),
               BT11 = (T + B022inv)*tmp1,
               BT12 = -sum1*tmp1,
               BT22 = (sum4+B011inv)*tmp1,
               bT1 = BT11*sum3 + BT12*sum2,
               bT2 = BT12*sum3 + BT22*sum2;

  // draw sigma^2 
  if (MHsteps == 2 || MHsteps == 3 || dontupdatemu) { // draw sigma^2 from full conditional
    double z = std::pow(((h[0]-mu)-phi*(h0-mu)), 2);  // TODO: more efficiently via sum1, sum2, etc.
    for (int j = 0; j < (T-1); j++) {
      z += std::pow((h[j+1]-mu)-phi*(h[j]-mu), 2);
    }
    z += std::pow(h0-mu, 2) * Bh0inv;
    if (MHcontrol > 0) {  // let's do a log normal random walk
      const double sigma2_prop = std::exp(R::rnorm(2 * std::log(sigma), MHcontrol)),
                   log_ar_prob = logacceptrateRW(sigma2_prop, std::pow(sigma, 2), Bsigma, T, z);

      if (log_ar_prob >= 0 || R::unif_rand() < std::exp(log_ar_prob)) {
        sigma = std::sqrt(sigma2_prop);
      }
    } else {  // either IG(-.5,0)-proposal or IG(1.5,1.5*Bsigma)-prior
      if (Gammaprior) {
        const double CT = .5*z,
                     sigma2_prop = 1/R::rgamma(cT, 1/CT);
        if (std::log(R::unif_rand()) <
            logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma)) {
          sigma = std::sqrt(sigma2_prop);
        }
      } else {
        const double CT = C0+.5*z;
        sigma = std::sqrt(1/R::rgamma(cT, 1/CT));
      }
    }
    sigma_prop = sigma;
  } else if (Gammaprior && MHsteps == 1) {  // draw sigma^2 marginalized over gamma, phi
    const double CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
    sigma_prop = std::sqrt(1/R::rgamma(cT, 1/CT));
  } else {
    ::Rf_error("The 1-block sampler is only implemented with a gamma prior for sigma^2");
  }


  // sampling of "betas" (i.e. phi and gamma)
  if (MHsteps == 3 || dontupdatemu) {

    {
      // sampling of phi from full conditional:
      const double gamma = (1-phi)*mu,  // = 0 if mu = 0
                   BTsqrt = sigma/std::sqrt(sum4+B011inv),
                   bT = (sum3-gamma*sum1)/(sum4+B011inv),
                   phi_prop = R::rnorm(bT, BTsqrt);

      double ar_prob = 0;
      if (priorlatent0 < 0.) { // needed only if prior of h0 depends on phi
        ar_prob += logdnorm(h0, mu, sigma/std::sqrt(1-phi_prop*phi_prop));
        ar_prob -= logdnorm(h0, mu, sigma/std::sqrt(1-phi*phi));
      } 
      ar_prob += logdbeta((phi_prop+1)/2, a0, b0);
      ar_prob -= logdbeta((phi+1)/2, a0, b0);
      ar_prob += logdnorm(phi, 0, sigma/std::sqrt(B011inv));
      ar_prob -= logdnorm(phi_prop, 0, sigma/std::sqrt(B011inv));

      if (std::log(R::runif(0, 1)) < ar_prob) {
        phi = phi_prop;
      }
    }

    if (!dontupdatemu) {
      // sampling of gamma from full conditional:
      const double gamma = (1-phi)*mu,
                   BTsqrt = sigma/std::sqrt(T+B022inv),
                   bT = (sum2-phi*sum1)/(T+B022inv),
                   gamma_prop = R::rnorm(bT, BTsqrt);

      double ar_prob = 0;
      ar_prob = logdnorm(h0, gamma_prop/(1-phi), sigma/std::sqrt(Bh0inv));
      ar_prob -= logdnorm(h0, gamma/(1-phi), sigma/std::sqrt(Bh0inv));
      ar_prob += logdnorm(gamma_prop, bmu*(1-phi), std::sqrt(Bmu)*(1-phi));
      ar_prob -= logdnorm(gamma, bmu*(1-phi), std::sqrt(Bmu)*(1-phi));
      ar_prob += logdnorm(gamma, 0, sigma/std::sqrt(B022inv));
      ar_prob -= logdnorm(gamma_prop, 0, sigma/std::sqrt(B022inv));

      if (std::log(R::unif_rand()) < ar_prob) {
        mu = gamma_prop/(1-phi);
      }
    }
  } else {
    double chol11 = std::sqrt(BT11),
           chol12 = (BT12/chol11),
           chol22 = std::sqrt(BT22-chol12*chol12);
    chol11 *= sigma;
    chol12 *= sigma;
    chol22 *= sigma;

    double phi_prop,
           gamma_prop;
    if (truncnormal) { // draw from truncated normal via inversion method
      const double quant_low = R::pnorm(-1, bT1, chol11, true, false),
                   quant_high = R::pnorm(1, bT1, chol11, true, false);
      phi_prop = R::qnorm((quant_low + R::unif_rand()*(quant_high-quant_low)),
          bT1, chol11, true, false);
      gamma_prop = R::rnorm(bT2 + chol12*((phi_prop-bT1)/chol11),
          chol22);
    } else { // draw from normal and reject (return) if outside
      const double innov = R::norm_rand();
      phi_prop = bT1 + chol11*innov;
      if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
        return {mu, phi, sigma};
      }
      gamma_prop = bT2 + chol12*innov + chol22*R::norm_rand();
    }

    // acceptance probability exp(R) calculated on a log scale
    const double phi_prop_const = 1-phi_prop,  // some temps used for acceptance probability
                 phi_const = 1-phi;

    double ar_prob;
    if (MHsteps == 2) {
      ar_prob = 0;
    } else {  // if (MHsteps == 1)
      ar_prob = logacceptrateGamma(std::pow(sigma_prop, 2), sigma*sigma, Bsigma);  // initialize R
    }

    if (priorlatent0 < 0.) {
      ar_prob += logdnorm(h0, gamma_prop/phi_prop_const, sigma_prop/std::sqrt(1-phi_prop*phi_prop));
      ar_prob -= logdnorm(h0, mu, sigma/std::sqrt(1-phi*phi));
    } else {
      ar_prob += logdnorm(h0, gamma_prop/phi_prop_const, std::sqrt(priorlatent0)*sigma_prop);
      ar_prob -= logdnorm(h0, mu, std::sqrt(priorlatent0)*sigma);
    }
    ar_prob += logdnorm(gamma_prop, bmu*phi_prop_const, std::sqrt(Bmu)*phi_prop_const);
    ar_prob -= logdnorm(mu*phi_const, bmu*phi_const, std::sqrt(Bmu)*phi_const);
    ar_prob += logdbeta((phi_prop+1)/2, a0, b0);
    ar_prob -= logdbeta((phi+1)/2, a0, b0);
    ar_prob += logdnorm(phi, 0, sigma/std::sqrt(B011inv));
    ar_prob -= logdnorm(phi_prop, 0, sigma_prop/std::sqrt(B011inv));
    ar_prob += logdnorm(mu*phi_const, 0, sigma/std::sqrt(B011inv));
    ar_prob -= logdnorm(gamma_prop, 0, sigma_prop/std::sqrt(B011inv));

    // accept/reject
    if (std::log(R::unif_rand()) < ar_prob) {
      mu = gamma_prop/(1-phi_prop);
      phi = phi_prop;
      if (MHsteps == 1) sigma = sigma_prop;
    }
  }

  return {mu, phi, sigma};
}

ReturnRegression regression_noncentered(
    const arma::vec& data,
    const double ht0,
    const arma::vec& ht,
    const arma::ivec& r,
    double mu,
    double phi,
    double sigma,
    const double Bsigma,
    const double a0,
    const double b0,
    const double bmu,
    const double Bmu,
    const bool truncnormal,
    const int MHsteps,
    const bool dontupdatemu,
    const double priorlatent0) {
  if (dontupdatemu) mu = 0;

  const int T = ht.size();

  if (MHsteps == 3 || dontupdatemu) {  // Gibbs-sample mu|sigma,... and sigma|mu,...

    // first, draw sigma from the full conditional posterior:
    double tmp1 = 0,
           tmp2 = 0;
    for (int j = 0; j < T; j++) {
      tmp1 += ht[j]*ht[j]*mix_varinv[r[j]];
      tmp2 += ht[j]*(data[j]-mix_mean[r[j]]-mu)*mix_varinv[r[j]];
    }
    const double BT11 = 1/(tmp1+1/Bsigma),
                 bT1 = BT11*tmp2;
    sigma = R::rnorm(bT1, std::sqrt(BT11));

    // TODO: check w.r.t. sign of sigma (low priority, 3 block is
    // practically useless anyway if dontupdatemu==false)
    if (!dontupdatemu) {
      // second, draw mu from the full conditional posterior:
      double tmp1 = 0,
             tmp2 = 0;
      for (int j = 0; j < T; j++) {
        tmp1 += mix_varinv[r[j]];
        tmp2 += (data[j]-mix_mean[r[j]]-sigma*ht[j])*mix_varinv[r[j]];
      }
      const int BT22 = 1/(tmp1+1/Bmu);
      const int bT2 = BT22*(tmp2 + bmu/Bmu);
      mu = R::rnorm(bT2, std::sqrt(BT22));
    }

  } else {  // Gibbs-sample mu and sigma jointly (regression) 
    double BT11 = 1/Bmu,
           BT12 = 0,
           BT22 = 1/Bsigma,
           bT1 = 0,
           bT2 = bmu/Bmu;

    for (int j = 0; j < T; j++) {
      const double tmp1 = mix_varinv[r[j]],
                   tmp2 = (data[j]-mix_mean[r[j]])*tmp1,
                   tmp3 = ht[j]*tmp1;
      BT11 += tmp1;
      BT12 -= tmp3;
      BT22 += ht[j]*tmp3;
      bT1 += ht[j]*tmp2;
      bT2 += tmp2;
    }

    {
      const double tmp = BT11*BT22-BT12*BT12;
      BT11 /= tmp;
      BT12 /= tmp;
      BT22 /= tmp;
    }

    {
      const double tmp = bT1;
      bT1 = BT11*tmp + BT12*bT2;
      bT2 = BT12*tmp + BT22*bT2;
    }

    const double chol11 = std::sqrt(BT11),
                 chol12 = (BT12/chol11),
                 chol22 = std::sqrt(BT22-chol12*chol12);

    arma::vec2 innov;
    innov.imbue([]() -> double { return R::norm_rand(); });
    sigma = bT1 + chol11*innov[0];
    mu = bT2 + chol12*innov[0] + chol22*innov[1];
  }


  // Sampling phi: find posterior mean muT and posterior variance sigma2T

  double sum_covar = ht0*ht[0],
         sum_var = ht0*ht0;
  double phi_prop;
  for (int j = 0; j < T-1; j++) {
    sum_covar += ht[j]*ht[j+1];
    sum_var += ht[j]*ht[j];
  }
  const double mean = sum_covar/sum_var,
               sd = 1/std::sqrt(sum_var);

  // actual sampling
  if (truncnormal) {  // draw from truncated normal via inversion method
    const double quant_low = R::pnorm(-1, mean, sd, true, false),
                 quant_high = R::pnorm(1, mean, sd, true, false);
    phi_prop = R::qnorm((quant_low + R::unif_rand()*(quant_high-quant_low)),
        mean, sd, true, false);
  } else {  // draw from normal and reject (return) if outside
    phi_prop = R::rnorm(mean, sd);
    if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
      return {mu, phi, std::fabs(sigma)};
    }
  }

  // now for the MH step, acceptance prob ar_prob
  double ar_prob;
  if (priorlatent0 < 0) { // only needed if prior of ht0 depends on phi
    ar_prob = std::exp(logdnorm(ht0, 0, 1/std::sqrt(1-phi_prop*phi_prop))
        - logdnorm(ht0, 0, 1/std::sqrt(1-phi*phi)));
  } else {
    ar_prob = 1;
  }
  ar_prob *= propBeta((phi_prop+1)/2, (phi+1)/2, a0, b0);
  //         ^^note that factor 1/2 from transformation of densities cancels

  // accept/reject
  if (R::unif_rand() < ar_prob) {
    phi = phi_prop;
  }

  return {mu, phi, std::fabs(sigma)};
}

bool draw_theta(
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    const arma::vec& y,
    const double h0,
    const double ht0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const Parameterization centering,
    const ProposalDiffusionKen& adaptation_proposal,
    const bool gammaprior,
    const Proposal sampler) {
  arma::vec proposed;
  switch (sampler) {
    case Proposal::RWMH:
      proposed = theta_propose_rwmh(phi, rho, sigma2, mu, adaptation_proposal);
      break;
    case Proposal::MALA:
      proposed = theta_propose_mala(phi, rho, sigma2, mu, y, h0, h, prior_phi, prior_rho, prior_sigma2, prior_mu, adaptation_proposal);
      break;
  }
  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
    mu_prop = proposed[3], prop_old_logdens = proposed[4], prop_new_logdens = proposed[5];
  if (centering == Parameterization::NONCENTERED) {
    exp_h_half_proposal_nc = arma::exp(.5 * noncentered_to_centered(mu_prop, std::sqrt(sigma2_prop), ht));
  }
  const arma::vec& exp_h_half_proposal = centering == Parameterization::CENTERED ? exp_h_half : exp_h_half_proposal_nc;

  const double log_acceptance =
    (theta_log_prior(phi_prop, rho_prop, sigma2_prop, mu_prop, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
     theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu_prop, y, h0, ht0, h, ht, exp_h_half_proposal, centering)) -
    (theta_log_prior(phi, rho, sigma2, mu, prior_phi, prior_rho, prior_sigma2, prior_mu, gammaprior) +
     theta_log_likelihood(phi, rho, sigma2, mu, y, h0, ht0, h, ht, exp_h_half, centering)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 || std::exp(log_acceptance) > R::unif_rand();
  if (accepted) {
    phi = phi_prop;
    rho = rho_prop;
    sigma2 = sigma2_prop;
    mu = mu_prop;
  }

  return accepted;
}

bool draw_thetamu_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    const double mu,
    const arma::vec& y,
    const double h0,
    const double ht0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const Parameterization centering,
    const ProposalDiffusionKen& adaptation_proposal,
    const bool gammaprior) {
  const arma::vec5 proposed = thetamu_propose(phi, rho, sigma2, adaptation_proposal);
  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
    prop_old_logdens = proposed[3], prop_new_logdens = proposed[4];
  if (centering == Parameterization::NONCENTERED) {
    exp_h_half_proposal_nc = arma::exp(.5 * (std::sqrt(sigma2_prop) * ht + mu));
  }
  const arma::vec& exp_h_half_proposal = centering == Parameterization::CENTERED ? exp_h_half : exp_h_half_proposal_nc;

  const double log_acceptance =
    (thetamu_log_prior(phi_prop, rho_prop, sigma2_prop, prior_phi, prior_rho, prior_sigma2, gammaprior) +
     theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu, y, h0, ht0, h, ht, exp_h_half_proposal, centering)) -
    (thetamu_log_prior(phi, rho, sigma2, prior_phi, prior_rho, prior_sigma2, gammaprior) +
     theta_log_likelihood(phi, rho, sigma2, mu, y, h0, ht0, h, ht, exp_h_half, centering)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 || std::exp(log_acceptance) > R::unif_rand();
  if (accepted) {
    phi = phi_prop;
    rho = rho_prop;
    sigma2 = sigma2_prop;
  }

  return accepted;
}

}

