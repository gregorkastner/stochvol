#ifndef _SAMPLING_LATENT_STATES_H_
#define _SAMPLING_LATENT_STATES_H_

// Functions related to the sampling of the latent states.
// This includes the all-without-a-loop (AWOL, McCausland et al., 2011)
// sampling of the latent vector in the auxiliary model from Omori et al. (2007)
// and the corresponding correction step thereafter

#include <RcppArmadillo.h>
#include "type_definitions.h"

namespace stochvol {

struct LatentVector {
  double h0;
  arma::vec h;
};

LatentVector draw_latent(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const bool correct);

arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::uvec& z,
    const double h0,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Parameterization centering);

arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double h0,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& proposed,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu);

arma::uvec draw_s_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const Parameterization centering);

}

#endif  // _SAMPLING_LATENT_STATES_H_

