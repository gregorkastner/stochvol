#include <Rcpp.h>
#include "auxmix.h"
#include "init_cache.h"

using namespace Rcpp;

void init_cache(
    List& l,
    const NumericVector& y) {
  const int n = y.length();
  if (!(l.containsElementNamed("initialized") && as<int>(l["length"]) == n)) {
    Rcout << "a";
    const int mix_count = sizeof(mix_prob)/sizeof(mix_prob[0]);

    l["post_dist"] = NumericMatrix(n, mix_count);
    l["eta"] = NumericVector(n);
    l["mixing_a"] = NumericVector(n);
    l["mixing_b"] = NumericVector(n);
    l["mixing_m"] = NumericVector(n);
    l["mixing_v"] = NumericVector(n);
    l["s"] = NumericVector(n);
    l["proposed"] = NumericVector(n);
    l["D"] = NumericVector(n);
    l["J1"] = NumericVector(n);
    l["L"] = NumericVector(n);
    l["f"] = NumericVector(n);
    l["F"] = NumericVector(n);
    l["hts"] = NumericVector(n);

    l["initialized"] = true;
    l["length"] = n;
  }
}
