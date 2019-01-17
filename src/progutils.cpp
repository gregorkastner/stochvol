#include "progutils.h"
#include "densities.h"

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
Rcpp::List cleanUp(const Rcpp::NumericVector & mu,
                   const Rcpp::NumericVector & phi,
		   const Rcpp::NumericVector & sigma,
		   const Rcpp::NumericMatrix & hstore,
		   const Rcpp::NumericVector & h0store,
		   const Rcpp::NumericVector & nustore,
		   const Rcpp::NumericMatrix & taustore,
		   const Rcpp::NumericMatrix & betastore) {
 int paracols;
 if (nustore.size() > 0) paracols = 4; else paracols = 3;
 
 Rcpp::NumericMatrix res(mu.length(), paracols); 
 res(Rcpp::_,0) = mu;
 res(Rcpp::_,1) = phi;
 res(Rcpp::_,2) = sigma;
 if (nustore.size() > 0) res(Rcpp::_,3) = nustore;

/* res.attr("dimnames") = Rcpp::List::create(
   R_NilValue, 
   Rcpp::CharacterVector::create("mu", "phi", "sigma")); */

 Rcpp::List val = Rcpp::List::create(
   Rcpp::_["para"] = res,
   Rcpp::_["latent"] = hstore,
   Rcpp::_["latent0"] = h0store,
   Rcpp::_["beta"] = betastore,
   Rcpp::_["tau"] = taustore);

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
void cholTridiag(const Rcpp::NumericVector & omega_diag, double omega_offdiag, Rcpp::NumericVector& chol_diag, Rcpp::NumericVector& chol_offdiag)
{
 chol_diag[0] = sqrt(omega_diag[0]);  // maybe speed up via iterators?
 for (int j = 1; j < omega_diag.length(); j++) {
  chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
  chol_diag[j] = sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
 }
}

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag, const Rcpp::NumericVector & covector, Rcpp::NumericVector& htmp)
{
 htmp[0] = covector[0]/chol_diag[0];
 for (int j = 1; j < chol_diag.length(); j++) {
  htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
 }
}

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag, const Rcpp::NumericVector & htmp, Rcpp::NumericVector& h)
{
 int T = chol_diag.length();
 h[T-1] = htmp[T-1]/chol_diag[T-1];
 for (int j = T-2; j >= 0; j--) {
  h[j] = (htmp[j] - chol_offdiag[j]*h[j+1])/chol_diag[j];
 }
}

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(const Rcpp::NumericVector& mixprob, Rcpp::IntegerVector& r, int T) {
 int index;
 Rcpp::NumericVector innov = Rcpp::runif(T); 
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
double newtonRaphson(double startval, double sumtau, int n,
                     double lower, double upper, double tol,
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

