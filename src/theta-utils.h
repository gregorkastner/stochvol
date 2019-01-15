#ifndef THETA_UTILS_H
#define THETA_UTILS_H

#include <Rcpp.h>
#include "parameterization.hpp"

double theta_log_likelihood(const double phi, const double rho,
                            const double sigma2, const double mu,
                            const Rcpp::NumericVector y,
                            const Rcpp::NumericVector h,
                            const Rcpp::NumericVector ht,
                            const Parameterization centering);

double theta_log_prior(const double phi, const double rho,
                       const double sigma2, const double mu,
                       const Rcpp::NumericVector prior_phi,
                       const Rcpp::NumericVector prior_rho,
                       const Rcpp::NumericVector prior_sigma2,
                       const Rcpp::NumericVector prior_mu,
                       const bool gammaprior);

Rcpp::NumericVector theta_transform(const double f, const double r,
                                    const double s, const double m);

Rcpp::NumericVector theta_transform_inv(const double phi, const double rho,
                                        const double sigma2, const double mu);

double theta_transform_log_det_jac(const double f, const double r,
                                   const double s, const double m);

double theta_transform_inv_log_det_jac(const double phi, const double rho,
                                       const double sigma2, const double mu);

Rcpp::NumericVector theta_proposal_stdev(const double phi, const double rho,
                                         const double sigma2, const double mu,
                                         const Rcpp::NumericVector y,
                                         const Rcpp::NumericVector h,
                                         const Rcpp::NumericVector ht,
                                         const double stdev);

Rcpp::NumericVector theta_propose(const double phi, const double rho,
                                  const double sigma2, const double mu,
                                  const Rcpp::NumericVector y,
                                  const Rcpp::NumericVector h,
                                  const Rcpp::NumericVector ht,
                                  const double stdev);

double theta_log_likelihood_c(const double phi, const double rho,
                              const double sigma2, const double mu,
                              const Rcpp::NumericVector y, const Rcpp::NumericVector h);

double theta_log_likelihood_nc(const double phi, const double rho,
                               const double sigma2, const double mu,
                               const Rcpp::NumericVector y, const Rcpp::NumericVector h);

#endif  // THETA_UTILS_H
