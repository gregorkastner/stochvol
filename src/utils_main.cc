#include <RcppArmadillo.h>
#include "utils_main.h"
#include "densities.h"

using namespace Rcpp;

namespace stochvol {

List cleanup(
    const arma::vec& mu,
    const arma::vec& phi,
    const arma::vec& sigma,
    const arma::mat& hstore,
    const arma::vec& h0store,
    const arma::vec& nustore,
    const arma::mat& taustore,
    const arma::mat& betastore) {
  int paracols;
  if (nustore.size() > 0) paracols = 4; else paracols = 3;

  CharacterVector coln(paracols);
  NumericMatrix res(mu.size(), paracols);
  arma::mat res_arma(res.begin(), mu.size(), paracols, false, false); 
  res_arma.col(0) = mu; coln[0] = "mu";
  res_arma.col(1) = phi; coln[1] = "phi";
  res_arma.col(2) = sigma; coln[2] = "sigma";
  if (nustore.size() > 0) {
    res_arma.col(3) = nustore; coln[3] = "nu";
  }
  colnames(res) = coln;

  List val = List::create(
      _["para"] = res,
      _["latent"] = hstore,
      _["latent0"] = h0store,
      _["beta"] = betastore,
      _["tau"] = taustore);

  return val;
}

int progressbar_init(
    const int N) {
  int show;
  REprintf("\n      ");
  if (N >= 2500) {
    for (int i = 0; i < 50+1; i++) REprintf(" ");
    show = N/50;
  }
  else {
    for (int i = 0; i < (N-1)/50+1; i++) REprintf(" ");
    show = 50;
  }
  REprintf("] 100%%\r  0%% [");
  R_FlushConsole();
  return show;
}

void progressbar_finish(
    const int N) {
  if (!(N % 50) && N >= 2500) REprintf("+");
  REprintf("] 100%%\n\n");
  R_FlushConsole();
}

}

