#ifndef THETA_UTILS_H
#define THETA_UTILS_H

#include <Rcpp.h>

// [[Rcpp::export]]
double theta_log_likelihood(const double phi, const double rho,
                            const double sigma2, const double mu,
                            const Rcpp::NumericVector y, const Rcpp::NumericVector h,
                            const Rcpp::CharacterVector centering);

double auxtheta_log_prior(const double phi, const double rho,
                          const double sigma2,
                          const Rcpp::NumericVector prior_phi,
                          const Rcpp::NumericVector prior_rho,
                          const Rcpp::NumericVector prior_sigma2);

// [[Rcpp::export]]
double theta_log_prior(const double phi, const double rho,
                       const double sigma2, const double mu,
                       const Rcpp::NumericVector prior_phi,
                       const Rcpp::NumericVector prior_rho,
                       const Rcpp::NumericVector prior_sigma2,
                       const Rcpp::NumericVector prior_mu);

// [[Rcpp::export]]
Rcpp::NumericVector theta_transform(const double f, const double r,
                                    const double s, const double m);

// [[Rcpp::export]]
Rcpp::NumericVector theta_transform_inv(const double phi, const double rho,
                                        const double sigma2, const double mu);

// [[Rcpp::export]]
double theta_transform_log_det_jac(const double f, const double r,
                                   const double s, const double m);

// [[Rcpp::export]]
double theta_transform_inv_log_det_jac(const double phi, const double rho,
                                       const double sigma2, const double mu);

// [[Rcpp::export]]
Rcpp::NumericVector theta_proposal_stdev(const double phi, const double rho,
                                         const double sigma2, const double mu,
                                         const Rcpp::NumericVector y, const Rcpp::NumericVector h,
                                         const double stdev = .1);

// [[Rcpp::export]]
Rcpp::NumericVector theta_propose(const double phi, const double rho,
                                  const double sigma2, const double mu,
                                  const Rcpp::NumericVector y, const Rcpp::NumericVector h,
                                  const double stdev = .1);

double theta_log_likelihood_c(const double phi, const double rho,
                              const double sigma2, const double mu,
                              const Rcpp::NumericVector y, const Rcpp::NumericVector h);

double theta_log_likelihood_nc(const double phi, const double rho,
                               const double sigma2, const double mu,
                               const Rcpp::NumericVector y, const Rcpp::NumericVector h);

// [[Rcpp::export]]
double auxtheta_log_posterior(const Rcpp::NumericVector auxtheta,
                              const Rcpp::NumericVector a,
                              const Rcpp::NumericVector b,
                              const Rcpp::NumericVector mm,
                              const Rcpp::NumericVector v,
                              const Rcpp::NumericVector d,
                              const Rcpp::NumericVector y_star,
                              const Rcpp::NumericVector prior_phi,
                              const Rcpp::NumericVector prior_rho,
                              const Rcpp::NumericVector prior_sigma2,
                              const Rcpp::NumericVector prior_mu);

#endif  // THETA_UTILS_H
