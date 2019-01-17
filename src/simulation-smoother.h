#ifndef SIMULATION_SMOOTHER_H
#define SIMULATION_SMOOTHER_H

#include <Rcpp.h>
#include "parameterization.hpp"

void simulation_smoother(
    double& eta0,
    Rcpp::List& cache,
    const double mu,
    const Parameterization centering,
    const double sigma2);

void simulation_smoother_c(
    Rcpp::NumericVector& eta,
    double& eta0,
    const double mu,
    const Rcpp::List& filter_results);

void simulation_smoother_nc(
    Rcpp::NumericVector& eta,
    double& eta0,
    const double mu,
    const Rcpp::List& filter_results,
    const double sigma2);

#endif  // SIMULATION_SMOOTHER_H
