#include <RcppArmadillo.h>
#include "progutils.h"
#include "densities.h"

using namespace Rcpp;

/* Contains the following code modules:

   a) some helper functions such as progress bar tools and return value
      prettifier

   b) some functions related to the Cholesky decomposition used for 
      sampling AWOL and efficiently solving the systems of linear
      equations

   c) function for inverse transform sampling

   d) a very basic Newton-Raphson algorithm for finding the root
      of dlogdnu (defined in densities.h)
*/

// a)
// Sums up results and prepares return value
List cleanUp(
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
 
 arma::mat res(mu.size(), paracols); 
 res.col(0) = mu;
 res.col(1) = phi;
 res.col(2) = sigma;
 if (nustore.size() > 0) res.col(3) = nustore;

 List val = List::create(
   _["para"] = res,
   _["latent"] = hstore,
   _["latent0"] = h0store,
   _["beta"] = betastore,
   _["tau"] = taustore);

 return val;
}

// sets up the progress bar
int progressbar_init(int N) {
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

// finalizes progress bar
void progressbar_finish(int N) {
 if (!(N % 50) && N >= 2500) REprintf("+");
 REprintf("] 100%%\n\n");
 R_FlushConsole();
}

// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(
    const arma::vec& omega_diag,
    double omega_offdiag,
    arma::vec& chol_diag,
    arma::vec& chol_offdiag) {
 chol_diag[0] = sqrt(omega_diag[0]);  // maybe speed up via iterators?
 for (int j = 1; j < int(omega_diag.size()); j++) {
  chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
  chol_diag[j] = sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
 }
}

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector,
    arma::vec& htmp) {
 htmp[0] = covector[0]/chol_diag[0];
 for (int j = 1; j < int(chol_diag.size()); j++) {
  htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
 }
}

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp,
    arma::vec& h) {
 int T = chol_diag.size();
 h[T-1] = htmp[T-1]/chol_diag[T-1];
 for (int j = T-2; j >= 0; j--) {
  h[j] = (htmp[j] - chol_offdiag[j]*h[j+1])/chol_diag[j];
 }
}

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(
    const arma::vec& mixprob,
    arma::ivec& r,
    int T) {
 int index;
 arma::vec innov = runif(T); 
 double temp;
 bool larger, smaller;
 for (int j = 0; j < T; j++) {
  index = (10-1)/2;  // start searching in the middle
  temp = innov[j]*mixprob[9 + 10*j];  // current (non-normalized) value
  larger = false;  // indicates that we already went up
  smaller = false; // indicates that we already went down
  while(true) {
   if (temp > mixprob[index +  10*j]) {
    if (smaller == true) {
     index++;
     break;
    }
    else {
    index++;
    larger = true;
    }
   }
   else {
    if (larger == true) {
     break;
    }
    else {
     if (index == 0) {
      break;
     }
     else {
      index--;
      smaller = true;
     }
    } 
   }
  }
 r[j] = index;
 }
}

// d)
// find the root of a function (Newton-Raphson)
double newtonRaphson(
    double startval,
    double sumtau,
    int n,
    double lower,
    double upper,
    double tol,
    int maxiter) {
 double x = startval;
 double error = R_PosInf;
 double xnew;
 bool converged = false;
 for (int i = 0; i < maxiter; i++) {
  xnew = x - dlogdnu(x, sumtau, n)/ddlogdnu(x, n);
  if (xnew > upper) xnew = upper; else if (xnew < lower) xnew = lower;
  error = fabs(xnew - x);
  x = xnew;
  if (error < tol) {
   converged = true;
   break;
  }
 }
 if (!converged) x = NA_REAL;
 return x;
}

