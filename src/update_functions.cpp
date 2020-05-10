#include <RcppArmadillo.h>
#include "update_functions.h"
#include "auxmix.h"
#include "progutils.h"
#include "densities.h"
#include "regression.h"
#include "h-sampler.h"
#include "theta-sampler.h"
#include "theta-utils.h"

using namespace Rcpp;

// update_terr performs an update of latent tau and df parameter nu
void update_terr(
    const arma::vec& data, 
    arma::vec& tau,
    double &nu,
    const double lower,
    const double upper) {

  int T = data.size();

  // **** STEP 1: Update tau ****

  double sumtau = 0.;

  for (int i = 0; i < T; i++) {
    // Watch out, R::rgamma(shape, scale), not Rf_rgamma(shape, rate)
    tau[i] = 1./R::rgamma((nu + 1.) / 2., 2. / (nu + exp(data[i])));
    sumtau += log(tau[i]) + 1/tau[i];
  }


  // **** STEP 2: Update nu ****

  double numean = newtonRaphson(nu, sumtau, T, lower, upper);
  double auxsd = sqrt(-1/ddlogdnu(numean, T)); 
  double nuprop = R::rnorm(numean, auxsd);
  double logR = logdnu(nuprop, sumtau, T) - logdnu(nu, sumtau, T) +
    logdnorm(nu, numean, auxsd) - logdnorm(nuprop, numean, auxsd);

  if (log(R::runif(0.,1.)) < logR && nuprop < upper && nuprop > lower) nu = nuprop;
}

// update performs one MCMC sampling step (normal errors):
void update_sv(
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
  int T = data.size();

  if (dontupdatemu) curpara[0] = 0; // just to be sure

  arma::vec omega_diag(T+1);  // contains diagonal elements of precision matrix
  double omega_offdiag;  // contains off-diag element of precision matrix (const)
  arma::vec chol_offdiag(T), chol_diag(T+1);  // Cholesky-factor of Omega
  arma::vec covector(T+1);  // holds covector (see McCausland et al. 2011)
  arma::vec htmp(T+1);  // intermediate vector for sampling h
  arma::vec hnew(T+1);  // intermediate vector for sampling h

  const double mu = curpara[0];
  const double phi = curpara[1];
  const double sigma2inv = pow(curpara[2], -2);

  double Bh0inv = 1./priorlatent0;
  if (priorlatent0 < 0) Bh0inv = 1-phi*phi;

  /*
   * Step (c): sample indicators
   */

  // calculate non-normalized CDF of posterior indicator probs

  if (centered_baseline) findMixCDF(mixprob, data-h);
  else findMixCDF(mixprob, data-mu-curpara[2]*h); 

  // find correct indicators (currently by inversion method)
  invTransformSampling(mixprob, r, T);

  /*
   * Step (a): sample the latent volatilities h:
   */

  if (centered_baseline) { // fill precision matrix omega and covector c for CENTERED para:

    omega_diag[0] = (Bh0inv + phi*phi) * sigma2inv;
    covector[0] = mu * (Bh0inv - phi*(1-phi)) * sigma2inv;

    for (int j = 1; j < T; j++) {
      omega_diag[j] = mix_varinv[r[j-1]] + (1+phi*phi)*sigma2inv; 
      covector[j] = (data[j-1] - mix_mean[r[j-1]])*mix_varinv[r[j-1]]
        + mu*(1-phi)*(1-phi)*sigma2inv;
    }
    omega_diag[T] = mix_varinv[r[T-1]] + sigma2inv;
    covector[T] = (data[T-1] - mix_mean[r[T-1]])*mix_varinv[r[T-1]]
      + mu*(1-phi)*sigma2inv;
    omega_offdiag = -phi*sigma2inv;  // omega_offdiag is constant

  } else { // fill precision matrix omega and covector c for NONCENTERED para:

    const double sigmainvtmp = sqrt(sigma2inv);
    const double phi2tmp = phi*phi;

    omega_diag[0] = phi2tmp + Bh0inv;
    covector[0] = 0.;

    for (int j = 1; j < T; j++) {
      omega_diag[j] = mix_varinv[r[j-1]]/sigma2inv + 1 + phi2tmp; 
      covector[j] = mix_varinv[r[j-1]]/sigmainvtmp*(data[j-1] - mix_mean[r[j-1]] - mu);
    }
    omega_diag[T] = mix_varinv[r[T-1]]/sigma2inv + 1;
    covector[T] = mix_varinv[r[T-1]]/sigmainvtmp*(data[T-1] - mix_mean[r[T-1]] - mu);
    omega_offdiag = -phi;  // omega_offdiag is constant
  } 

  // Cholesky decomposition
  cholTridiag(omega_diag, omega_offdiag, chol_diag, chol_offdiag);

  // Solution of Chol*x = covector ("forward algorithm")
  forwardAlg(chol_diag, chol_offdiag, covector, htmp);

  htmp += as<arma::vec>(rnorm(T+1));

  // Solution of (Chol')*x = htmp ("backward algorithm")
  backwardAlg(chol_diag, chol_offdiag, htmp, hnew);

  h = hnew.tail(T);
  h0 = hnew[0];

  /*
   * Step (b): sample mu, phi, sigma
   */

  if (centered_baseline) {  // this means we have C as base
    curpara = regressionCentered(h0, h, mu, phi, curpara[2],
        C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv,
        B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps, dontupdatemu, priorlatent0);

    if (parameterization == 3) {  // this means we should interweave
      double h0_alter;
      htmp = (h-curpara[0])/curpara[2];
      h0_alter = (h0-curpara[0])/curpara[2];
      curpara = regressionNoncentered(data, h0_alter, htmp, r,
          curpara[0], curpara[1], curpara[2], Bsigma, a0, b0, bmu, Bmu,
          truncnormal, MHsteps, dontupdatemu, priorlatent0);
      h = curpara[0] + curpara[2]*htmp;
      h0 = curpara[0] + curpara[2]*h0_alter;
    }


  } else {  // NC as base

    curpara = regressionNoncentered(data, h0, h, r, mu, phi, curpara[2],
        Bsigma, a0, b0, bmu, Bmu, truncnormal, MHsteps,
        dontupdatemu, priorlatent0);

    if (parameterization == 4) {  // this means we should interweave
      double h0_alter;
      htmp = curpara[0] + curpara[2]*h;
      h0_alter = curpara[0] + curpara[2]*h0;
      curpara = regressionCentered(h0_alter, htmp, curpara[0], curpara[1], curpara[2],
          C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv,
          Gammaprior, truncnormal, MHcontrol, MHsteps,
          dontupdatemu, priorlatent0);
      h = (htmp-curpara[0])/curpara[2];
      h0 = (h0_alter-curpara[0])/curpara[2];
    }
  }
}

void update_svl (
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    double& phi,
    double& rho,
    double& sigma2,
    double& mu,
    arma::vec& h,
    arma::vec& ht,
    stochvol::Adaptation<4>& adaptation,
    const arma::vec& prior_phi,
    const arma::vec& prior_rho,
    const arma::vec& prior_sigma2,
    const arma::vec& prior_mu,
    const bool use_mala,
    const bool gammaprior,
    const bool correct,
    const arma::ivec& strategy) {
  // only centered
  h = draw_latent(y, y_star, d, h, ht, phi, rho, sigma2, mu, prior_mu[0], prior_mu[1], correct);
  ht = (h - mu) / std::sqrt(sigma2);
  arma::vec exp_h_half = arma::exp(.5 * h);  // cache exp() calculations
  arma::vec exp_h_half_tilde = arma::exp(.5 * (std::sqrt(sigma2) * ht + mu));

  const Proposal proposal = use_mala ? Proposal::MALA : Proposal::RWMH;
  const auto adapted_proposal = adaptation.get_proposal();
  for (int ipar : strategy) {
    const Parameterization par = Parameterization(ipar);
  //  if (dontupdatemu) {
  //      draw_thetamu_rwMH(
  //          phi, rho, sigma2, mu, y, h, ht,
  //          prior_phi,
  //          prior_rho,
  //          prior_sigma2,
  //          par,
  //          proposal_chol.submat(0, 0, 2, 2),
  //          proposal_chol_inv.submat(0, 0, 2, 2),
  //          gammaprior);
  //  } else {
    const bool theta_updated = draw_theta(
          phi, rho, sigma2, mu,
          y, h, ht, exp_h_half, exp_h_half_tilde,
          prior_phi,
          prior_rho,
          prior_sigma2,
          prior_mu,
          par,
          adapted_proposal,
          gammaprior,
          proposal);
//  }

    if (theta_updated) {
      switch (par) {
        case Parameterization::CENTERED:
          ht = (h - mu)/ std::sqrt(sigma2);
          exp_h_half_tilde = arma::exp(.5 * (std::sqrt(sigma2) * ht + mu));
          break;
        case Parameterization::NONCENTERED:
          h = std::sqrt(sigma2) * ht + mu;
          exp_h_half = arma::exp(.5 * h);
          break;
      }
    }
  }
  adaptation.register_sample(theta_transform_inv(phi, rho, sigma2, mu));

  return;
}

