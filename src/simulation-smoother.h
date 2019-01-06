#ifndef SIMULATION_SMOOTHER_H
#define SIMULATION_SMOOTHER_H

#include <Rcpp.h>

Rcpp::List simulation_smoother(const double mu,
                               const Rcpp::List filter_results,
                               const Rcpp::CharacterVector centering);

Rcpp::List simulation_smoother_c(const double mu,
                                 const Rcpp::List filter_results);
  
Rcpp::List simulation_smoother_nc(const double mu,
                                  const Rcpp::List filter_results);

#endif  // SIMULATION_SMOOTHER_H
