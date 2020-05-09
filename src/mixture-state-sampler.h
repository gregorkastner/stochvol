#ifndef MIXTURE_STATE_SAMPLER_H
#define MIXTURE_STATE_SAMPLER_H

#include <RcppArmadillo.h>
#include "parameterization.h"

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

#endif  // MIXTURE_STATE_SAMPLER_H

