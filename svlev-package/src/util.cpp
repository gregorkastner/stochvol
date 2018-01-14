#include <Rcpp.h>
#include <string>
using namespace Rcpp;

List fnAugKalmanFilterCenteredCpp(const double phi, const double sigma2, const double rho,
                          const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                          const NumericVector d, const NumericVector y_star,
                          const double mu_mu, const double sigma2_mu) {
  const R_xlen_t n = y_star.size();
  
  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * sigma * exp(m/2);
  const NumericVector h_ts = b * v * gamma_ts;
  const double j_t22 = sigma2 * (1-rho*rho);
  const double h1_var = sigma2/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = -1;
  double P = h1_var;
  
  // Init returned values
  NumericVector D(n);
  NumericVector J(n);
  NumericVector L(n);
  NumericVector f(n);
  NumericVector F(n);
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < int(n); i++) {
    D[i] = P + v[i]*v[i];
    K = (phi*P + h_ts[i]*v[i])/D[i];
    L[i] = phi - K;
    J[i] = h_ts[i] - K*v[i];
    P = phi*P*L[i] + h_ts[i]*J[i] + j_t22;
    f[i] = y_star[i] - m[i] - a_star;
    F[i] = -A_star;
    a_star = a[i]*gamma_ts[i] + phi*a_star + K*f[i];
    A_star = phi*(A_star+1) - 1 + K*F[i];
    
    q += F[i]*f[i]/D[i];
    Q += F[i]*F[i]/D[i];
  }
  
  return List::create(
    _["sigma2"] = sigma2,
    _["D"] = D,
    _["J1"] = J,
    _["L"] = L,
    _["f"] = f,
    _["F"] = F,
    _["hts"] = h_ts,
    _["v"] = v,
    _["Q"] = Q,
    _["q"] = q,
    _["jt22"] = j_t22,
    _["h1var"] = h1_var
  );
}

List fnAugKalmanFilterNonCenteredCpp(const double phi, const double sigma2, const double rho,
                                  const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                                  const NumericVector d, const NumericVector y_star,
                                  const double mu_mu, const double sigma2_mu) {
  const R_xlen_t n = y_star.size();
  
  const double sigma = sqrt(sigma2);
  const NumericVector gamma_ts = d * rho * exp(m/2);
  const NumericVector h_ts = b * v * gamma_ts;
  const double j_t22 = (1-rho*rho);
  const double h1_var = 1/(phi*phi > 1-1e-8 ? 1e-8 : (1-phi*phi));
  
  double a_star = 0;
  double A_star = 0;
  double P = h1_var;
  
  // Init returned values
  NumericVector D(n);
  NumericVector J(n);
  NumericVector L(n);
  NumericVector f(n);
  NumericVector F(n);
  double Q = 1/sigma2_mu;
  double q = mu_mu * Q;
  
  // Init for the loop
  double K;
  
  // Main loop
  for (int i = 0; i < int(n); i++) {
    D[i] = sigma2*P + v[i]*v[i];
    K = (sigma*phi*P + h_ts[i]*v[i])/D[i];
    L[i] = phi - sigma*K;
    J[i] = h_ts[i] - K*v[i];
    P = phi*P*L[i] + h_ts[i]*J[i] + j_t22;
    f[i] = y_star[i] - m[i] - sigma*a_star;
    F[i] = 1.0 - sigma*A_star;
    a_star = a[i]*gamma_ts[i] + phi*a_star + K*f[i];
    A_star = phi*A_star + K*F[i];
    
    q += F[i]*f[i]/D[i];
    Q += F[i]*F[i]/D[i];
  }
  
  return List::create(
    _["sigma2"] = sigma2,
    _["D"] = D,
    _["J1"] = J,
    _["L"] = L,
    _["f"] = f,
    _["F"] = F,
    _["hts"] = h_ts,
    _["v"] = v,
    _["Q"] = Q,
    _["q"] = q,
    _["jt22"] = j_t22,
    _["h1var"] = h1_var
  );
}

//' Augmented Kalman filter
//' 
//' Implementation of the augmented Kalman filter for the special
//' case of the current model. The augmented Kalman filter is a
//' modified version of the standard Kalman filter that
//' considers the the error terms' correlation and the constant
//' \eqn{\mu} in the state equation.
//' Used to calculate the likelihood of \eqn{\theta=(\phi,\sigma,\rho)},
//' and to calculate the posterior of \eqn{\mu}.
//' For a specification see Appendix B of Nakajima, 2009.
//' @author hdarjus \email{hdarjus@gmail.com}
//' @param dPhi parameter \eqn{\phi}
//' @param dSigma2 square of parameter \eqn{\sigma}
//' @param dRho parameter \eqn{\rho}
//' @param va numeric vector of mixture states constants
//' @param vb numeric vector of mixture states constants
//' @param vm numeric vector of mixture states constants
//' @param vv numeric vector of mixture states constants
//' @param vD vector of signs
//' @param vYStar vector of y*
//' @param dMuMu prior mean of \eqn{\mu}
//' @param dSigma2Mu prior variance of \eqn{\mu}
//' @return List with elements \code{D, J1, L, f, F, hts, v, Q, q, jt22, h1var}, from which
//'   \code{D, J1, L, f, F, h1var, v} are numeric vectors of length T,
//'   and \code{Q, q, jt22, h1var} are numbers. All of them are just partial results helping
//'   later calculations.
//' @references Nakajima, Jouchi, and Yasuhiro Omori.
//'   "Leverage, heavy-tails and correlated jumps in stochastic volatility models."
//'   Computational Statistics & Data Analysis 53.6 (2009): 2335-2353.
// [[Rcpp::export]]
List fnAugKalmanFilterCpp(const double phi, const double sigma2, const double rho,
                          const NumericVector a, const NumericVector b, const NumericVector m, const NumericVector v,
                          const NumericVector d, const NumericVector y_star,
                          const double mu_mu, const double sigma2_mu, const StringVector centering) {
  List result;
  std::string scentering = as<std::string>(centering);
  if (scentering == "centered") {
    result = fnAugKalmanFilterCenteredCpp(phi, sigma2, rho, a, b, m, v, d, y_star, mu_mu, sigma2_mu);
  } else if (scentering == "non-centered") {
    result = fnAugKalmanFilterNonCenteredCpp(phi, sigma2, rho, a, b, m, v, d, y_star, mu_mu, sigma2_mu);
  } else {
    ::Rf_error("c++: unknown type of centering");
  }
  return result;
}
