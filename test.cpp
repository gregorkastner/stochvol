#include <RcppArmadillo.h>
#include <RInside.h>

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
    draws = 10000
    burnin = 1000
    designmatrix = NA
    priormu = c(0, 100); priorphi = c(5, 1.5); priorsigma = 1
    priorrho = c(3, 5); priorbeta = c(0, 10000)
    thinpara = 1; thinlatent = 1; thintime = 1
    quiet = FALSE
    startpara = list(phi=0.95, mu=-10, rho=-0.4, sigma=0.1)
    startlatent=rep_len(-10, len)
    strategies <- c("centered", "non-centered")
    expert <- list(parameterization = rep(strategies, 5),
                   mhcontrol = 0.1,
                   gammaprior = TRUE)
    )";
  R.parseEvalQ(evalstr);

  Rcpp::NumericMatrix designmatrix = R["dmat"];
  cout << designmatrix(0, 0) << std::endl;
  return 0;
}
