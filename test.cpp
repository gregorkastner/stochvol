#include <RcppArmadillo.h>
#include <RInside.h>
#include "src/sampler.h"
#include "src/sampler-leverage.h"

using namespace std;
using namespace Rcpp;

int main(int argc, char** argv) {
  RInside R(argc, argv);
  const string evalstr = R"(
    mydevlib <- "~/R/under_development"
    library(coda)
    library(stochvol, lib.loc = mydevlib)

    set.seed(421)
    len <- 100
    dmat <- cbind(rep_len(1, len), rgamma(len, 0.5, 0.25))
    datreg <- svsim(len, phi=0.9, mu=-10, sigma=0.1, rho=-0.4)
    datreg$y <- datreg$y + as.numeric(dmat %*% rbind(-1, 2))
    y <- datreg$y
    priormu <- c(0, 100); priorphi = c(5, 1.5); priorsigma = 1
    priorrho <- c(3, 5); priorbeta = c(0, 10000)
    priornu <- NA_real_
    thinpara <- 1; thinlatent <- 1; thintime <- 1
    quiet <- FALSE
    startpara <- list(phi=0.95, mu=-10, rho=-0.4, sigma=0.1)
    startlatent <- rep_len(-10, len)
    strategies <- c("centered", "non-centered")
    expert <- list(parameterization = rep(strategies, 5),
                   mhcontrol = 0.1,
                   gammaprior = TRUE)
    strategy <- expert$parameterization
    )";
  R.parseEvalQ(evalstr);

  NumericMatrix designmatrix = R["dmat"];
  NumericVector y = R["y"];
  List startpara = R["startpara"];
  NumericVector startvol = R["startlatent"];
  NumericVector priorbeta = R["priorbeta"];
  NumericVector priordf = R["priornu"];
  CharacterVector strategy = R["strategy"];
  List svdraws = svsample_cpp(
      y, 10000, 1000, designmatrix,
      0, 100*100, 5, 1.5, 1, 1, 1,
      startpara, startvol,
      false, true, 3, 2,
      100000000, 1000000000000,
      -1, true, false, 1e-9, false,
      priordf, priorbeta, -1);
  List svldraws = svlsample_cpp(
      y, 10000, 1000, designmatrix,
      1, 1, 1, startpara, startvol,
      5, 1.5, 3, 5, 0.5, 0.5,
      0, 100, priorbeta(0), priorbeta(1),
      false, 1e-9, 0.1, true, strategy);
  return 0;
}
