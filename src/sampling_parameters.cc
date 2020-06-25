#include <RcppArmadillo.h>
#include <cmath>
#include "sampling_parameters.h"
#include "utils_parameters.h"
#include "utils_latent_states.h"  // TODO sample rather from true SV
#include "type_definitions.h"
#include "densities.h"

using namespace Rcpp;

arma::vec regression_centered(
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

  double Bh0inv = 1./priorlatent0;
  if (priorlatent0 < 0) Bh0inv = 1-phi*phi;

  if (dontupdatemu) mu = 0;

  const int T = h.size();
  double z, CT, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, tmp1,
         BT11, BT12, BT22, bT1 = 0, bT2 = 0, chol11, chol12, chol22, phi_prop,
         gamma_prop, tmpR, tmpR2, logR;
  double R = -10000.;
  double sigma2_prop = -10000.;
  double sigma_prop = -10000.;
  arma::vec innov(1);
  arma::vec quant(2);

  // first calculate bT and BT:
  sum1 = h[0];
  sum3 = h0*h[0];
  sum4 = h[0]*h[0];
  for (int j = 1; j < T-1; j++) {
    sum1 += h[j];
    sum3 += h[j-1]*h[j];
    sum4 += h[j]*h[j];
  }
  sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
  sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
  sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
  sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2

  tmp1 = 1/(((sum4 + B011inv)*(T+B022inv)-sum1*sum1));
  BT11 = (T + B022inv)*tmp1;
  BT12 = -sum1*tmp1;
  BT22 = (sum4+B011inv)*tmp1;

  bT1 = BT11*sum3 + BT12*sum2;
  bT2 = BT12*sum3 + BT22*sum2;

  // draw sigma^2 
  if (MHsteps == 2 || MHsteps == 3 || dontupdatemu == true) { // draw sigma^2 from full conditional
    z = pow(((h[0]-mu)-phi*(h0-mu)),2);  // TODO: more efficiently via sum1, sum2, etc.
    for (int j = 0; j < (T-1); j++) {
      z += pow((h[j+1]-mu)-phi*(h[j]-mu),2);
    }
    z += (h0-mu)*(h0-mu)*Bh0inv;
    if (MHcontrol > 0) {  // let's do a log normal random walk
      sigma2_prop = exp(R::rnorm(log(sigma*sigma), MHcontrol));
      logR = logacceptrateRW(sigma2_prop, sigma*sigma, Bsigma, T, z);

      if (log(R::runif(0, 1)) < logR) sigma = sqrt(sigma2_prop);
    } else {  // either IG(-.5,0)-proposal or IG(1.5,1.5*Bsigma)-prior
      if (Gammaprior) {
        CT = .5*z;
        sigma2_prop = 1/R::rgamma(cT, 1/CT);
        if (log(R::runif(0, 1)) <
            logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma)) {
          sigma = sqrt(sigma2_prop);
        }
      } else {
        CT = C0+.5*z;
        sigma = sqrt(1/R::rgamma(cT, 1/CT));
      }
    }
  } else if (MHsteps == 1) {  // draw sigma^2 marginalized over gamma, phi
    if (Gammaprior) {
      CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
      sigma2_prop = 1/R::rgamma(cT, 1/CT);
    }
  }


  // sampling of "betas" (i.e. phi and gamma)
  if (MHsteps == 3 || dontupdatemu == true) {

    // sampling of phi from full conditional:
    double gamma = (1-phi)*mu;  // = 0 if mu = 0
    double BTsqrt = sigma/sqrt(sum4+B011inv);
    double bT = (sum3-gamma*sum1)/(sum4+B011inv);
    phi_prop = R::rnorm(bT, BTsqrt);

    R = 0;
    if (priorlatent0 < 0.) { // needed only if prior of h0 depends on phi
      R += logdnorm(h0, mu, sigma/sqrt(1-phi_prop*phi_prop));
      R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
    } 
    R += logdbeta((phi_prop+1)/2, a0, b0);
    R -= logdbeta((phi+1)/2, a0, b0);
    R += logdnorm(phi, 0, sigma/sqrt(B011inv));
    R -= logdnorm(phi_prop, 0, sigma/sqrt(B011inv));

    if (log(R::runif(0, 1)) < R) {
      phi = phi_prop;
    }

    if (dontupdatemu == false) {
      // sampling of gamma from full conditional:
      gamma = (1-phi)*mu;
      BTsqrt = sigma/sqrt(T+B022inv);
      bT = (sum2-phi*sum1)/(T+B022inv);
      gamma_prop = R::rnorm(bT, BTsqrt);

      R = logdnorm(h0, gamma_prop/(1-phi), sigma/sqrt(Bh0inv));
      R -= logdnorm(h0, gamma/(1-phi), sigma/sqrt(Bh0inv));
      R += logdnorm(gamma_prop, bmu*(1-phi), sqrt(Bmu)*(1-phi));
      R -= logdnorm(gamma, bmu*(1-phi), sqrt(Bmu)*(1-phi));
      R += logdnorm(gamma, 0, sigma/sqrt(B022inv));
      R -= logdnorm(gamma_prop, 0, sigma/sqrt(B022inv));

      if (log(R::runif(0, 1)) < R) {
        mu = gamma_prop/(1-phi);
      }
    }
  } else { 
    // Some more temps needed for sampling the betas
    chol11 = sqrt(BT11);
    chol12 = (BT12/chol11);
    chol22 = sqrt(BT22-chol12*chol12);
    chol11 *= sigma;
    chol12 *= sigma;
    chol22 *= sigma;

    if (truncnormal) { // draw from truncated normal via inversion method
      quant[0] = R::pnorm(-1, bT1, chol11, true, false);
      quant[1] = R::pnorm(1, bT1, chol11, true, false);
      phi_prop = R::qnorm((quant[0] + R::runif(0, 1)*(quant[1]-quant[0])),
          bT1, chol11, true, false);
      gamma_prop = R::rnorm(bT2 + chol12*((phi_prop-bT1)/chol11),
          chol22);
    }
    else { // draw from normal and reject (return) if outside
      innov[0] = R::rnorm(0, 1);
      phi_prop = bT1 + chol11*innov[0];
      if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
        return {mu, phi, sigma};
      }
      else gamma_prop = bT2 + chol12*innov[0] + chol22*R::rnorm(0, 1);
    }

    // acceptance probability exp(R) calculated on a log scale
    tmpR = 1-phi_prop;  // some temps used for acceptance probability
    tmpR2 = 1-phi;

    R = 0.;  // initialize R
    if (MHsteps == 2) {
      sigma_prop = sigma;  // sigma was accepted/rejected independently
    } else if (MHsteps == 1) {
      sigma_prop = sqrt(sigma2_prop);  // accept sigma jointly with "betas"
      R = logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma);  // initialize R
    }

    if (priorlatent0 < 0.) {
      R += logdnorm(h0, gamma_prop/tmpR, sigma_prop/sqrt(1-phi_prop*phi_prop));
      R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
    } else {
      R += logdnorm(h0, gamma_prop/tmpR, sqrt(priorlatent0)*sigma_prop);
      R -= logdnorm(h0, mu, sqrt(priorlatent0)*sigma);
    }
    R += logdnorm(gamma_prop, bmu*tmpR, sqrt(Bmu)*tmpR);
    R -= logdnorm(mu*tmpR2, bmu*tmpR2, sqrt(Bmu)*tmpR2);
    R += logdbeta((phi_prop+1)/2, a0, b0);
    R -= logdbeta((phi+1)/2, a0, b0);
    R += logdnorm(phi, 0, sigma/sqrt(B011inv));
    R -= logdnorm(phi_prop, 0, sigma_prop/sqrt(B011inv));
    R += logdnorm(mu*tmpR2, 0, sigma/sqrt(B011inv));
    R -= logdnorm(gamma_prop, 0, sigma_prop/sqrt(B011inv));

    // accept/reject
    if (log(R::runif(0, 1)) < R) {
      mu = gamma_prop/(1-phi_prop);
      phi = phi_prop;
      if (MHsteps == 1) sigma = sigma_prop;
    }
  }

  return {mu, phi, sigma};
}

arma::vec regression_noncentered(
    const arma::vec& data,
    const double h0,
    const arma::vec& h,
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

  const int T = h.size();
  double sumtmp1, sumtmp2, expR, phi_prop, BT11, BT12, BT22, bT1, bT2, tmp1,
         tmp2, tmp3, tmp, chol11, chol12, chol22, tmpmean, tmpsd;
  arma::vec innov(2);
  arma::vec quant(2);

  if (MHsteps == 3 || dontupdatemu) {  // Gibbs-sample mu|sigma,... and sigma|mu,...

    // first, draw sigma from the full conditional posterior:
    tmp1 = 0; 
    tmp2 = 0; 
    for (int j = 0; j < T; j++) {
      tmp1 += h[j]*h[j]*mix_varinv[r[j]];
      tmp2 += h[j]*(data[j]-mix_mean[r[j]]-mu)*mix_varinv[r[j]];
    }
    BT11 = 1/(tmp1+1/Bsigma);
    bT1 = BT11*tmp2; 
    //  REprintf("old: %f, new: mean %f and sd %f\n", sigma, bT1, sqrt(BT11));
    sigma = R::rnorm(bT1, sqrt(BT11));

    // TODO: check w.r.t. sign of sigma (low priority, 3 block is
    // practically useless anyway if dontupdatemu==false)
    if (!dontupdatemu) {
      // second, draw mu from the full conditional posterior:
      tmp1 = 0; 
      tmp2 = 0; 
      for (int j = 0; j < T; j++) {
        tmp1 += mix_varinv[r[j]];
        tmp2 += (data[j]-mix_mean[r[j]]-sigma*h[j])*mix_varinv[r[j]];
      }
      BT22 = 1/(tmp1+1/Bmu);
      bT2 = BT22*(tmp2 + bmu/Bmu);
      //  REprintf("old: %f, new: mean %f and sd %f\n\n", mu, bT2, sqrt(BT22));
      mu = R::rnorm(bT2, sqrt(BT22));
    }

  } else {  // Gibbs-sample mu and sigma jointly (regression) 
    BT11 = 1/Bmu;
    BT12 = 0;
    BT22 = 1/Bsigma;
    bT1 = 0;
    bT2 = bmu/Bmu;

    for (int j = 0; j < T; j++) {
      tmp1 = mix_varinv[r[j]];
      tmp2 = (data[j]-mix_mean[r[j]])*tmp1;
      tmp3 = h[j]*tmp1;
      BT11 += tmp1;
      BT12 -= tmp3;
      BT22 += h[j]*tmp3;
      bT1 += h[j]*tmp2;
      bT2 += tmp2;
    }

    tmp = BT11*BT22-BT12*BT12;
    BT11 /= tmp;
    BT12 /= tmp;
    BT22 /= tmp;

    tmp = bT1;
    bT1 = BT11*tmp + BT12*bT2;
    bT2 = BT12*tmp + BT22*bT2;

    chol11 = sqrt(BT11);
    chol12 = (BT12/chol11);
    chol22 = sqrt(BT22-chol12*chol12);

    innov = rnorm(2);
    sigma = bT1 + chol11*innov[0];
    mu = bT2 + chol12*innov[0] + chol22*innov[1];
  }


  // Sampling phi: find posterior mean muT and posterior variance sigma2T

  sumtmp1 = h0*h[0];
  sumtmp2 = h0*h0;
  for (int j = 0; j < T-1; j++) {
    sumtmp1 += h[j]*h[j+1];
    sumtmp2 += h[j]*h[j];
  }
  tmpmean = sumtmp1/sumtmp2;
  tmpsd = 1/sqrt(sumtmp2);

  // actual sampling
  if (truncnormal) {  // draw from truncated normal via inversion method
    quant[0] = R::pnorm(-1, tmpmean, tmpsd, true, false);
    quant[1] = R::pnorm(1, tmpmean, tmpsd, true, false);
    phi_prop = R::qnorm((quant[0] + R::runif(0, 1)*(quant[1]-quant[0])),
        tmpmean, tmpsd, true, false);
  }
  else {  // draw from normal and reject (return) if outside
    phi_prop = R::rnorm(tmpmean, tmpsd); 
    if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
      return {mu, phi, fabs(sigma)};
    }
  }

  // now for the MH step, acceptance prob expR
  if (priorlatent0 < .0) { // only needed if prior of ho depends on phi
    expR  = exp(logdnorm(h0, 0, 1/sqrt(1-phi_prop*phi_prop))
        - logdnorm(h0, 0, 1/sqrt(1-phi*phi)));
  } else expR = 1;
  expR *= propBeta((phi_prop+1)/2, (phi+1)/2, a0, b0);
  // ^^note that factor 1/2 from transformation of densities cancels

  // accept/reject
  if (R::runif(0, 1) < expR) phi = phi_prop;

  return {mu, phi, fabs(sigma)};
}

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
    const stochvol::Adaptation::Result& adaptation_proposal,
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

bool draw_thetamu_rwMH(
    double& phi,
    double& rho,
    double& sigma2,
    const double mu,
    const arma::vec& y,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& exp_h_half,
    arma::vec& exp_h_half_proposal_nc,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const Parameterization centering,
    const stochvol::Adaptation::Result& adaptation_proposal,
    const bool gammaprior) {
  const arma::vec5 proposed = thetamu_propose(phi, rho, sigma2, y, h, ht, adaptation_proposal);
  const double phi_prop = proposed[0], rho_prop = proposed[1], sigma2_prop = proposed[2],
    prop_old_logdens = proposed[3], prop_new_logdens = proposed[4];
  if (centering == Parameterization::NONCENTERED) {
    exp_h_half_proposal_nc = arma::exp(.5 * (std::sqrt(sigma2_prop) * ht + mu));
  }
  const arma::vec& exp_h_half_proposal = centering == Parameterization::CENTERED ? exp_h_half : exp_h_half_proposal_nc;

  const double log_acceptance =
    (thetamu_log_prior(phi_prop, rho_prop, sigma2_prop, prior_phi, prior_rho, prior_sigma2, gammaprior) +
     theta_log_likelihood(phi_prop, rho_prop, sigma2_prop, mu, y, h, ht, exp_h_half_proposal, centering)) -
    (thetamu_log_prior(phi, rho, sigma2, prior_phi, prior_rho, prior_sigma2, gammaprior) +
     theta_log_likelihood(phi, rho, sigma2, mu, y, h, ht, exp_h_half, centering)) -
    (prop_new_logdens - prop_old_logdens);

  const bool accepted = log_acceptance > 0 || std::exp(log_acceptance) > R::runif(0, 1);
  if (accepted) {
    phi = phi_prop;
    rho = rho_prop;
    sigma2 = sigma2_prop;
  }

  return accepted;
}

