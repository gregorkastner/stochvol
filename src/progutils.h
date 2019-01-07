#ifndef _PROGUTILS_H_
#define _PROGUTILS_H_

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

#include <Rcpp.h>

// a)
// Sums up results and prepares return value
Rcpp::List cleanUp(const Rcpp::NumericVector & mu,
                   const Rcpp::NumericVector & phi,
		   const Rcpp::NumericVector & sigma,
		   const Rcpp::NumericMatrix & hstore,
		   const Rcpp::NumericVector & h0store,
		   const Rcpp::NumericVector & nustore,
		   const Rcpp::NumericMatrix & taustore,
		   const Rcpp::NumericMatrix & betastore);

// sets up the progress bar
int progressbar_init(int N);

// adds one '+' to progress bar
inline void progressbar_print() {
 REprintf("+");
 R_FlushConsole();
}

// finalizes progress bar
void progressbar_finish(int N);

// to store (some) values of h
inline void store_h(double * h, double * hstore, int timethin,
                    int hstorelength,
		    double h0, double * h0store,
		    Rcpp::NumericVector curpara,
		    bool centered_baseline) {

 if (centered_baseline) {
  for (int j = 0; j < hstorelength; j++) hstore[j] = h[timethin*j];
  *h0store = h0;
 } else {
  for (int j = 0; j < hstorelength; j++) hstore[j] = curpara[0] + curpara[2]*h[timethin*j];
  *h0store = curpara[0] + curpara[2]*h0;
 }
}

// to store (some) values of tau
inline void store_tau(double * tau, double * taustore, int timethin, int hstorelength) {
 for (int j = 0; j < hstorelength; j++) taustore[j] = tau[timethin*j];
}

// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(const Rcpp::NumericVector & omega_diag, double omega_offdiag,
                 double * chol_diag, double * chol_offdiag);

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag,
                const Rcpp::NumericVector & covector, double * htmp);

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag,
                 const Rcpp::NumericVector & htmp, double * h);

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(const double * const, int * r, int T);

// d)
// find the root of a function (Newton-Raphson)
double newtonRaphson(double startval, double sumtau, int n,
                     double lower = R_NegInf, double upper = R_PosInf,
		     double tol = 1e-03, int maxiter = 50);
#endif
