#ifndef ADAPTATION_H
#define ADAPTATION_H

#include <RcppArmadillo.h>

// Encapsulate adaptation logic
// Adaptation happens after every batch_size draws.
//
// Notation from DOI 10.1007/s11222-008-9110-y
template<int dim, int batch_size>
class Adaptation {
public:
  struct Result {
    const double scale;
    const arma::mat::fixed<dim, dim> Covariance;
    const arma::mat::fixed<dim, dim> Precision;  // Covariance_inv
    const arma::mat::fixed<dim, dim> Covariance_chol;
    const arma::mat::fixed<dim, dim> Covariance_chol_inv;
  };

  Adaptation (const double _lambda)
    : lambda{_lambda}, alpha{calculate_alpha(lambda)} { }

  Adaptation ()
    : Adaptation(0.1) { }

  inline
  void step_gamma () {
    const long i_previous = std::lround(std::pow(gamma / C, alpha));
    gamma = C * std::pow(i_previous + 1, -alpha);
  }

  inline
  void register_sample (const arma::vec::fixed<dim>& sample) {
    const int i = state.get_i_batch();
    draws_batch.col(i) = sample;
    if (i == batch_size - 1) {
      step_gamma();
      state.update_covariance(draws_batch, gamma);
      updated_proposal = true;
    }
  }

  inline
  const Result get_proposal () {
    // TODO cache all of these
    // TODO update scale
    updated_proposal = false;
    const auto covariance = state.get_covariance();
    const auto covariance_chol = arma::chol(covariance, "lower");
    return {scale, covariance, arma::inv_sympd(covariance), covariance_chol, arma::inv(arma::trimatl(covariance))};
  }

private:
  class State {
  public:
    State () {
      mu.fill(0);
      Sigma = arma::eye(arma::size(Sigma));
    }

    inline
    int get_i_batch () {
      const int i_batch_current = i_batch;
      i_batch = (i_batch + 1) % batch_size;
      return i_batch_current;
    }

    // overwrites draws
    inline
    void update_covariance(arma::mat::fixed<dim, batch_size> draws, const double gamma) {
      draws.each_col() -= mu;
      mu += gamma * (arma::sum(draws, 1) / batch_size - mu);
      Sigma += gamma * (draws * draws.t() / (batch_size - 1) - Sigma);
    }

    inline
    const arma::mat::fixed<dim, dim> get_covariance () const {
      return Sigma;
    }

  private:
    int i_batch = 0;  // use get_i_batch() to access i_batch
    arma::vec::fixed<dim> mu;
    arma::mat::fixed<dim, dim> Sigma;
  };

  State state;
  arma::mat::fixed<dim, batch_size> draws_batch;
  bool updated_proposal = false;

  const double lambda;  // controls the speed of adaptation
  const double alpha;  // log-speed of adaptation, derived from lambda
  const double C = 20;  // constant factor in the speed of adaptation
  double gamma = C;  // initialize gamma
  double scale = 0.2;

  inline
  static double calculate_alpha (const double l) {
    const double lambda_reciprocal_value = 1 / (1 + l);
    return lambda_reciprocal_value + (1 - lambda_reciprocal_value) / 64;
  }
};

#endif

