#ifndef _UTILS_MAIN_H_
#define _UTILS_MAIN_H_

// Helper functions for the main sampler

#include <RcppArmadillo.h>

// Sums up results and prepares return value
Rcpp::List cleanup(
    const arma::vec& mu,
    const arma::vec& phi,
    const arma::vec& sigma,
    const arma::mat& hstore,
    const arma::vec& h0store,
    const arma::vec& nustore,
    const arma::mat& taustore,
    const arma::mat& betastore);

// sets up the progress bar
int progressbar_init(
    const int N);

// adds one '+' to progress bar
inline
void progressbar_print() {
  ::REprintf("+");
  ::R_FlushConsole();
}

// finalizes progress bar
void progressbar_finish(
    const int N);

#endif  // _UTILS_MAIN_H_

