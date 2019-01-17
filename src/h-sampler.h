#ifndef H_SAMPLER_H
#define H_SAMPLER_H

#include <Rcpp.h>
#include "parameterization.hpp"

void draw_h_auxiliary(
    Rcpp::NumericVector& h,
    Rcpp::List& cache,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    const Rcpp::NumericVector& s,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering);

void draw_latent_auxiliaryMH(
    Rcpp::NumericVector& h,
    Rcpp::List& cache,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& y_star,
    const Rcpp::NumericVector& d,
    const Rcpp::NumericVector& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma);

#endif  // H_SAMPLER_H
