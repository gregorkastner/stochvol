// Functions related to sampling the indicators
// Constants relating to the approximation of log(chisq) through
// normal mixture (Omori et al., 2007) can be found in auxmix.h

#include <Rcpp.h>
#include "auxmix.h"

// Non-normalized posterior probabilities
void findMixprobs(Rcpp::NumericVector& mixprob, const Rcpp::NumericVector & datanorm)  {
 int T = datanorm.length();
 int tmp; 
 for (int c = 0; c < T; c++) {  // SLOW (10*T calls to exp)!
  tmp = 10*c;
  for (int r = 0; r < 10; r++) {
   mixprob(tmp+r) = exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
  }
 }
}

// Cumulative sum over columns of a matrix
void colCumsums(Rcpp::NumericVector& x, int const nrow, int const ncol) {
 int tmp;
 for (int c = 0; c < ncol; c++) {
  tmp = c*nrow;
  for (int r = 1; r < nrow; r++) {
   x(tmp+r) = x(tmp+r-1) + x(tmp+r);
  }
 }
}

// Combines findMixprobs() and colCumsums() (see above) into one function
void findMixCDF(Rcpp::NumericVector& mixprob, const Rcpp::NumericVector & datanorm)  {
 int T = datanorm.length();
 int tmp; 
 for (int c = 0; c < T; c++) {  // SLOW (10*T calls to exp)!
  tmp = 10*c;
  mixprob(tmp) = exp(mix_pre[0]-(datanorm[c]-mix_mean[0])*(datanorm[c]-mix_mean[0])*mix_2varinv[0]);
  for (int r = 1; r < 10; r++) {
   mixprob(tmp+r) = mixprob(tmp+r-1) + exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
  }
 }
}
