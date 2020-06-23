#include <RcppArmadillo.h>
#include "utils_latent_states.h"
#include <cmath>

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

void invTransformSampling(
    const arma::vec& mixprob,
    arma::ivec& r,
    int T) {
  int index;
  arma::vec innov = Rcpp::runif(T); 
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


void findMixprobs(
    arma::vec& mixprob,
    const arma::vec& datanorm)  {
  int T = datanorm.size();
  int tmp; 
  for (int c = 0; c < T; c++) {  // TODO slow (10*T calls to exp)!
    tmp = 10*c;
    for (int r = 0; r < 10; r++) {
      mixprob[tmp+r] = exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
    }
  }
}

void colCumsums(
    arma::vec& x,
    int const nrow,
    int const ncol) {
  int tmp;
  for (int c = 0; c < ncol; c++) {
    tmp = c*nrow;
    for (int r = 1; r < nrow; r++) {
      x[tmp+r] = x[tmp+r-1] + x[tmp+r];
    }
  }
}

void findMixCDF(
    arma::vec& mixprob,
    const arma::vec& datanorm)  {
  int T = datanorm.size();
  int tmp; 
  for (int c = 0; c < T; c++) {  // TODO slow (10*T calls to exp)!
    tmp = 10*c;
    mixprob[tmp] = exp(mix_pre[0]-(datanorm[c]-mix_mean[0])*(datanorm[c]-mix_mean[0])*mix_2varinv[0]);
    for (int r = 1; r < 10; r++) {
      mixprob[tmp+r] = mixprob[tmp+r-1] + exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
    }
  }
}

double h_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  const double sigma = std::sqrt(sigma2);
  const double sigma2_inv = 1 / sigma2;
  const double rho_const = std::sqrt(1 - rho * rho);
  const int n = y.size();
  const arma::vec exp_h_half = arma::exp(0.5 * h);
  double result = -.5 * std::pow(h[0] - mu, 2) * sigma2_inv * (1 - phi * phi);  // log p(h_1 | theta)
  for (int t = 0; t < n - 1; t++) {
    result += -.5 * std::pow(h[t + 1] - (mu + phi * (h[t] - mu)), 2) * sigma2_inv;
    result += -.5 * std::pow((y[t] - (exp_h_half[t] * rho * (h[t + 1] - mu - phi * (h[t] - mu)) / sigma)) / (exp_h_half[t] * rho_const), 2) - .5 * h[t];
  }
  result += -.5 * std::pow(y[n - 1] / exp_h_half[n - 1], 2) - .5 * h[n - 1];
  return result;
}

double h_aux_log_posterior(
    const arma::vec& h,  // centered
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
  const int n = y_star.size();
  static const int mix_count = mix_a.n_elem;
  const double sigma = std::sqrt(sigma2);
  static const arma::vec::fixed<10> exp_m_half = arma::exp(mix_mean * .5);
  const arma::vec C = rho * sigma * exp_m_half;  // re-used constant

  double result = -.5 * std::pow(h[0] - mu, 2) / sigma2 * (1 - phi * phi);  // log p(h_1 | theta)
  for (int t = 0; t < n; t++) {
    double subresult = 0;
    if (t < n - 1) {
      for (int j = 0; j < mix_count; j++) {
        const double h_mean = mu + phi * (h[t] - mu) + d[t] * C[j] * mix_a[j];
        const double h_var = mix_var[j] * std::pow(C[j] * mix_b[j], 2) + sigma2 * (1 - rho * rho);
        const double yh_cov = d[t] * C[j] * mix_b[j] * mix_var[j];
        const double y_mean = h[t] + mix_mean[j] + yh_cov / h_var * (h[t + 1] - h_mean);
        const double y_var = (1 - std::pow(yh_cov, 2) / (mix_var[j] * h_var)) * mix_var[j];
        subresult += std::exp(-.5 * (std::pow(y_star[t] - y_mean, 2) / y_var + std::pow(h[t + 1] - h_mean, 2) / h_var)) / std::sqrt(y_var * h_var) * mix_prob[j];
      }
    } else {
      for (int j = 0; j < mix_count; j++) {
        subresult += std::exp(-.5 * std::pow(y_star[t] - (h[t] + mix_mean[j]), 2) / mix_var[j]) / std::sqrt(mix_var[j]) * mix_prob[j];
      }
    }
    result += std::log(subresult);
  }

  return result;
}

