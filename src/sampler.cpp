// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "auxmix.h"
#include "progutils.h"
#include "densities.h"

using namespace Rcpp;

// Step (b): sample mu, phi, sigma - __CENTERED__ version:
NumericVector regressionCentered(
    double h0, const Rcpp::NumericVector &h,
    double mu, double phi, double sigma,
    double C0, double cT, double Bsigma,
    double a0, double b0,
    double bmu, double Bmu,
    double B011inv, double B022inv,
    bool Gammaprior, bool truncnormal, double MHcontrol, int MHsteps,
    const bool dontupdatemu, const double priorlatent0) {
  
  double Bh0inv = 1./priorlatent0;
  if (priorlatent0 < 0) Bh0inv = 1-phi*phi;
  
  if (dontupdatemu) mu = 0;
  
  int T = h.length();
  double z, CT, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, tmp1,
    BT11, BT12, BT22, bT1 = 0, bT2 = 0, chol11, chol12, chol22, phi_prop,
    gamma_prop, tmpR, tmpR2, logR;
  double R = -10000.;
  double sigma2_prop = -10000.;
  double sigma_prop = -10000.;
  Rcpp::NumericVector innov(2);
  Rcpp::NumericVector quant(2);
  
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
      sigma2_prop = exp(Rcpp::rnorm(1,log(sigma*sigma), MHcontrol)[0]);
      logR = logacceptrateRW(sigma2_prop, sigma*sigma, Bsigma, T, z);
      
      if (log(Rcpp::runif(1)[0]) < logR) sigma = sqrt(sigma2_prop);
    } else {  // either IG(-.5,0)-proposal or IG(1.5,1.5*Bsigma)-prior
      if (Gammaprior) {
        CT = .5*z;
        sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
        if (log(Rcpp::as<double>(Rcpp::runif(1))) <
          logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma)) {
          sigma = sqrt(sigma2_prop);
        }
      } else {
        CT = C0+.5*z;
        sigma = sqrt(1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT)));
      }
    }
  } else if (MHsteps == 1) {  // draw sigma^2 marginalized over gamma, phi
    if (Gammaprior) {
      CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
      sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
    }
  }
  
  
  // sampling of "betas" (i.e. phi and gamma)
  if (MHsteps == 3 || dontupdatemu == true) {
    
    // sampling of phi from full conditional:
    double gamma = (1-phi)*mu;  // = 0 if mu = 0
    double BTsqrt = sigma/sqrt(sum4+B011inv);
    double bT = (sum3-gamma*sum1)/(sum4+B011inv);
    phi_prop = as<double>(Rcpp::rnorm(1, bT, BTsqrt));
    
    R = 0;
    if (priorlatent0 < 0.) { // needed only if prior of h0 depends on phi
      R += logdnorm(h0, mu, sigma/sqrt(1-phi_prop*phi_prop));
      R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
    } 
    R += logdbeta((phi_prop+1)/2, a0, b0);
    R -= logdbeta((phi+1)/2, a0, b0);
    R += logdnorm(phi, 0, sigma/sqrt(B011inv));
    R -= logdnorm(phi_prop, 0, sigma/sqrt(B011inv));
    
    if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
      phi = phi_prop;
    }
    
    if (dontupdatemu == false) {
      // sampling of gamma from full conditional:
      gamma = (1-phi)*mu;
      BTsqrt = sigma/sqrt(T+B022inv);
      bT = (sum2-phi*sum1)/(T+B022inv);
      gamma_prop = as<double>(Rcpp::rnorm(1, bT, BTsqrt));
      
      R = logdnorm(h0, gamma_prop/(1-phi), sigma/sqrt(Bh0inv));
      R -= logdnorm(h0, gamma/(1-phi), sigma/sqrt(Bh0inv));
      R += logdnorm(gamma_prop, bmu*(1-phi), sqrt(Bmu)*(1-phi));
      R -= logdnorm(gamma, bmu*(1-phi), sqrt(Bmu)*(1-phi));
      R += logdnorm(gamma, 0, sigma/sqrt(B022inv));
      R -= logdnorm(gamma_prop, 0, sigma/sqrt(B022inv));
      
      if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
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
      quant = Rcpp::pnorm(Rcpp::NumericVector::create(-1,1), bT1, chol11);
      phi_prop = (Rcpp::qnorm((quant[0] + Rcpp::runif(1)*(quant[1]-quant[0])),
                              bT1, chol11))[0];
      gamma_prop = (Rcpp::rnorm(1, bT2 + chol12*((phi_prop-bT1)/chol11),
                                chol22))[0];
    }
    else { // draw from normal and reject (return) if outside
      innov = Rcpp::rnorm(1);
      phi_prop = bT1 + chol11*innov[0];
      if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
        Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
        return ret;
      }
      else gamma_prop = bT2 + chol12*innov[0] + chol22*Rcpp::rnorm(1)[0];
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
    if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
      mu = gamma_prop/(1-phi_prop);
      phi = phi_prop;
      if (MHsteps == 1) sigma = sigma_prop;
    }
  }
  
  Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
  return ret;
}

// Step (b): sample mu, phi, sigma - __NONCENTERED__ version:
Rcpp::NumericVector regressionNoncentered(
    const Rcpp::NumericVector &data,
    double h0, const Rcpp::NumericVector &h,
    const int * const r,
    double mu, double phi, double sigma,
    double Bsigma, double a0, double b0,
    double bmu, double Bmu,
    bool truncnormal, int MHsteps,
    const bool dontupdatemu, const double priorlatent0) {
  if (dontupdatemu) mu = 0;
  
  int T = h.length();
  double sumtmp1, sumtmp2, expR, phi_prop, BT11, BT12, BT22, bT1, bT2, tmp1,
  tmp2, tmp3, tmp, chol11, chol12, chol22, tmpmean, tmpsd;
  Rcpp::NumericVector innov(2);
  Rcpp::NumericVector quant(2);
  
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
    sigma = as<double>(rnorm(1, bT1, sqrt(BT11)));
    
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
      mu = as<double>(Rcpp::rnorm(1, bT2, sqrt(BT22)));
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
    quant = Rcpp::pnorm(Rcpp::NumericVector::create(-1,1), tmpmean, tmpsd);
    phi_prop = (Rcpp::qnorm((quant[0] + Rcpp::runif(1)*(quant[1]-quant[0])),
                            tmpmean, tmpsd))[0];
  }
  else {  // draw from normal and reject (return) if outside
    phi_prop = (Rcpp::rnorm(1, tmpmean, tmpsd))[0]; 
    if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
      Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, fabs(sigma));
      return ret;
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
  if (Rcpp::as<double>(Rcpp::runif(1)) < expR) phi = phi_prop;
  
  Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, fabs(sigma));
  return ret;
}

// update_terr performs an update of latent tau and df parameter nu
void update_terr(const NumericVector &data, 
                 double *tau, double &nu,
                 const double lower, const double upper) {
  
  int T = data.length();
  
  // **** STEP 1: Update tau ****
  
  double sumtau = 0.;
  
  for (int i = 0; i < T; i++) {
    // Watch out, Rf_rgamma(shape, scale), not Rf_rgamma(shape, rate)
    tau[i] = 1./::Rf_rgamma((nu + 1.) / 2., 2. / (nu + exp(data(i))));
    sumtau += log(tau[i]) + 1/tau[i];
  }
  
  
  // **** STEP 2: Update nu ****
  
  double numean = newtonRaphson(nu, sumtau, T, lower, upper);
  double auxsd = sqrt(-1/ddlogdnu(numean, T)); 
  double nuprop = ::Rf_rnorm(numean, auxsd);
  double logR = logdnu(nuprop, sumtau, T) - logdnu(nu, sumtau, T) +
    logdnorm(nu, numean, auxsd) - logdnorm(nuprop, numean, auxsd);
  
  if (log(::Rf_runif(0.,1.)) < logR && nuprop < upper && nuprop > lower) nu = nuprop;
}

// update performs one MCMC sampling step (normal errors):
// [[Rcpp::interfaces(r, cpp)]]
void update(const NumericVector &data, double *curpara_in, double *h_in,
            double &h0, double *mixprob, int *r,
            const bool centered_baseline, const double C0, const double cT,
            const double Bsigma, const double a0, const double b0,
            const double bmu, const double Bmu, const double B011inv,
            const double B022inv, const bool Gammaprior, const bool truncnormal,
            const double MHcontrol, const int MHsteps, const int parameterization,
            const bool dontupdatemu, const double priorlatent0) {
  
  int T = data.length();
  
  NumericVector h(T);  // h needs to be NumericVector in current implementation
  for (int j = 0; j < T; j++) h(j) = h_in[j];  // maybe revisit to avoid copying
  
  NumericVector curpara(3);  // curpara needs to be NumericVector in current implementation
  for (int j = 0; j < 3; j++) curpara(j) = curpara_in[j];  // maybe revisit to avoid copying
  
  if (dontupdatemu) curpara(0) = 0; // just to be sure
  
  NumericVector omega_diag(T+1);  // contains diagonal elements of precision matrix
  double omega_offdiag;  // contains off-diag element of precision matrix (const)
  NumericVector chol_offdiag(T), chol_diag(T+1);  // Cholesky-factor of Omega
  NumericVector covector(T+1);  // holds covector (see McCausland et al. 2011)
  NumericVector htmp(T+1);  // intermediate vector for sampling h
  NumericVector hnew(T+1);  // intermediate vector for sampling h
  
  
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
  cholTridiag(omega_diag, omega_offdiag, &chol_diag(0), &chol_offdiag(0));
  
  // Solution of Chol*x = covector ("forward algorithm")
  forwardAlg(chol_diag, chol_offdiag, covector, &htmp(0));
  
  htmp = htmp + rnorm(T+1);
  
  // Solution of (Chol')*x = htmp ("backward algorithm")
  backwardAlg(chol_diag, chol_offdiag, htmp, &hnew(0));
  
  for (int j = 0; j < T; j++) h(j) = hnew(j+1);  // TODO: REVISIT!!
  h0 = hnew(0);
  
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
  
  for (int j = 0; j < T; j++) h_in[j] = h(j);  // maybe revisit to avoid copying
  for (int j = 0; j < 3; j++) curpara_in[j] = curpara(j);  // maybe revisit to avoid copying
  }

// [[Rcpp::export]]
SEXP sampler(const SEXP y_in,
             const SEXP draws_in,
             const SEXP burnin_in,
             const SEXP X_in,
             const SEXP bmu_in,
             const SEXP Bmu_in,
             const SEXP a0_in,
             const SEXP b0_in,
             const SEXP Bsigma_in,
             const SEXP thin_in,
             const SEXP timethin_in,
             const SEXP startpara_in,
             const SEXP startvol_in,
             const SEXP keeptau_in,
             const SEXP quiet_in,
             const SEXP para_in,
             const SEXP MHsteps_in,
             const SEXP B011_in,
             const SEXP B022_in,
             const SEXP mhcontrol_in,
             const SEXP gammaprior_in,
             const SEXP truncnormal_in,
             const SEXP offset_in,
             const SEXP dontupdatemu_in,
             const SEXP priordf_in,
             const SEXP priorbeta_in,
             const SEXP priorlatent0_in) {

 //RNGScope scope;       // just in case no seed has been set at R level
 GetRNGstate(); // "by hand" because RNGScope isn't safe if return
                // variables are declared afterwards

 // convert SEXP into Rcpp-structures (no copy at this point)
 NumericVector y(y_in), startvol(startvol_in);
 arma::vec armay(y.begin(), y.length(), false);

 NumericMatrix X(X_in);

 int T = y.length();
 int p = X.ncol();
 arma::mat armaX(X.begin(), T, p, false);
 
 // should we model the mean as well?
 bool regression;
 if (ISNA(X(0,0))) regression = false; else regression = true;
 NumericVector priorbeta(priorbeta_in);

 // prior mean and precision matrix for the regression part (currently fixed)
 arma::vec priorbetamean(p); priorbetamean.fill(priorbeta(0));
 arma::mat priorbetaprec(p, p, arma::fill::zeros);
 priorbetaprec.diag() += 1./(priorbeta(1)*priorbeta(1));
 
 List startpara(startpara_in);

 double priorlatent0 = as<double>(priorlatent0_in);

 // number of MCMC draws
 int burnin = as<int>(burnin_in);
 int draws  = as<int>(draws_in);
 int N 	    = burnin + draws;
 
 // prior parameters
 double bmu    = as<double>(bmu_in);
 double Bmu    = as<double>(Bmu_in);
 double a0     = as<double>(a0_in);
 double b0     = as<double>(b0_in);
 double Bsigma = as<double>(Bsigma_in);
 
 // thinning parameters
 int timethin = as<int>(timethin_in);
 int thin     = as<int>(thin_in);

 bool keeptau = as<bool>(keeptau_in);
 
 // verbosity control
 bool verbose = !as<bool>(quiet_in);

 // "expert" settings:
 double B011inv         = 1/as<double>(B011_in);
 double B022inv         = 1/as<double>(B022_in);
 bool Gammaprior        = as<bool>(gammaprior_in);
 bool truncnormal       = as<bool>(truncnormal_in);
 double MHcontrol       = as<double>(mhcontrol_in);
 int MHsteps            = as<int>(MHsteps_in);
 int parameterization   = as<int>(para_in);
 bool centered_baseline = as<int>(para_in) % 2; // 0 for C, 1 for NC baseline

 // offset:
 double offset		= as<double>(offset_in);

 // t-errors:
 bool terr;
 NumericVector priordf(priordf_in);

 if (ISNA(priordf(0))) terr = false; else terr = true;

 // indicator telling us whether assume mu = 0 fixed
 bool dontupdatemu      = as<bool>(dontupdatemu_in);

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
  else return LogicalVector::create(NA_LOGICAL);  // not implemented!
 }

 if (dontupdatemu == true && MHsteps == 1) { // not implemented (would be easy, though)
  return LogicalVector::create(NA_LOGICAL);
 }
 
 // initialize the variables:
 NumericVector sigma2inv(N+1, pow(as<double>(startpara["sigma"]), -2));
 NumericVector phi(N+1, as<double>(startpara["phi"]));
 NumericVector mu(N+1, as<double>(startpara["mu"]));
 
 NumericVector h(T);  // contains h1 to hT, but not h0!
 for (int i = 0; i < T; i++) h[i] = startvol[i];  // need to make a manual copy (!)
 if (!centered_baseline) h = (h-mu[0])*sqrt(sigma2inv[0]);
 
 double h0;
 
 IntegerVector r(T);  // mixture indicators
 
 int hstorelength = (T-1)/timethin+1;
 NumericMatrix hstore(hstorelength, draws/thin);
 NumericVector h0store(draws/thin);

 NumericMatrix mixprob(10, T);  // mixture probabilities
 
 NumericVector data = log(y*y + offset);  // commonly used transformation
 arma::vec armadata(data.begin(), data.length(), false);
 
 NumericVector datastand = clone(data);  // standardized "data" (different for t-errors)
 arma::vec armadatastand(datastand.begin(), datastand.length(), false);
 
 NumericVector curpara(3);  // holds mu, phi, sigma in every iteration
 curpara[0] = mu[0];
 curpara[1] = phi[0];
 curpara[2] = 1/sqrt(sigma2inv[0]);
 
 // some stuff for the t-errors
 double nu;
 if (terr) nu = as<double>(startpara["nu"]);
 NumericVector tau(T, 1.);

 NumericVector nustore(terr * (N+1), nu);
 NumericMatrix taustore(hstorelength * keeptau, draws/thin);

 // some stuff for the regression part
 NumericVector curbeta(p, .012345);
 arma::vec armabeta(curbeta.begin(), curbeta.length(), false);
 arma::vec normalizer;
 arma::mat Xnew = as<arma::mat>(clone(X));
 arma::vec ynew = as<arma::vec>(clone(y));
 arma::mat postprecchol(p, p);
 arma::mat postpreccholinv(p, p);
 arma::mat postcov(p, p);
 arma::vec postmean(p);
 arma::vec armadraw(p);
 NumericMatrix betastore(regression * (N+1), p);
 arma::mat armabetastore(betastore.begin(), betastore.nrow(), betastore.ncol(), false);

 // initializes the progress bar
 // "show" holds the number of iterations per progress sign
 int show = 0;
 if (verbose) show = progressbar_init(N);
 
 for (int i = 0; i < N; i++) {  // BEGIN main MCMC loop

  // print a progress sign every "show" iterations
  if (verbose) if (!(i % show)) progressbar_print();

  if (regression) {
   armadatastand = armadata = log(square(armay - armaX * armabeta));
  }
  
  if (terr) {
   if (centered_baseline) update_terr(data - h, &tau(0), nu, priordf(0), priordf(1));
   else update_terr(data - curpara(0) - curpara(2) * h, &tau(0), nu, priordf(0), priordf(1));
   
   datastand = data - log(tau);
  }

  // a single MCMC update: update indicators, latent volatilities,
  // and parameters ONCE
  update(datastand, &curpara(0), &h(0), h0, &mixprob(0,0), &r(0), centered_baseline, C0, cT,
         Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior,
	 truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);

  if (regression) { // update betas (regression)
   normalizer = exp(-h/2);
   Xnew = armaX;
   Xnew.each_col() %= normalizer;
   ynew = as<arma::vec>(y) % normalizer;
   
   // cholesky factor of posterior precision matrix
   postprecchol = arma::chol(Xnew.t() * Xnew + priorbetaprec);
   
   // inverse cholesky factor of posterior precision matrix 
   postpreccholinv = arma::solve(arma::trimatu(postprecchol), arma::eye<arma::mat>(p, p));

   // posterior covariance matrix and posterior mean vector
   postcov = postpreccholinv * postpreccholinv.t();
   postmean = postcov * (Xnew.t() * ynew + priorbetaprec * priorbetamean);

   armadraw = rnorm(p);

   // posterior betas
   armabeta = postmean + postpreccholinv * armadraw;
  }
  
  // storage:
  if (!((i+1) % thin)) if (i >= burnin) {  // this means we should store h
   store_h(&h(0), &hstore(0, (i-burnin)/thin), timethin, hstorelength,
           h0, &h0store((i-burnin)/thin), curpara, centered_baseline);
   if (keeptau && terr) {
    store_tau(&tau(0), &taustore(0, (i-burnin)/thin), timethin, hstorelength);
   }
  }
  mu[i+1] = curpara[0];
  phi[i+1] = curpara[1];
  sigma2inv[i+1] = 1/(curpara[2]*curpara[2]);
  if (terr) nustore[i+1] = nu;
  if (regression) armabetastore.row(i+1) = armabeta.t();
 }  // END main MCMC loop

 if (verbose) progressbar_finish(N);  // finalize progress bar

 // Prepare return value and return
 PutRNGstate();
 return cleanUp(mu, phi, sqrt(1/sigma2inv), hstore, h0store, nustore, taustore, betastore);
}

