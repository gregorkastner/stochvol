#ifndef ADAPTATION_H
#define ADAPTATION_H

#include <RcppArmadillo.h>
#include <vector>

namespace stochvol {
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

    struct Storage {
      const double gamma;
      const double scale;
      const double rate_acceptance;
    };

    Adaptation (const int _memory_size, const double _target_acceptance = 0.234, const double _lambda = 0.1)
      : lambda{_lambda}, target_acceptance{_target_acceptance}, alpha{calculate_alpha(_lambda)} {
      if (target_acceptance <= 0.1 || target_acceptance >= 0.8) {
        Rcpp::warning("Target acceptance rate should be between 10\% and 80\%");
      }
      memory.reserve(_memory_size);
    }

    inline
    double step_gamma (const double gamma_previous) const {
      const long i_previous = std::lround(std::pow(C / gamma_previous, 1.0 / alpha));
      return C * std::pow(i_previous + 1, -alpha);
    }

    inline
    void register_sample (const arma::vec::fixed<dim>& sample) {
      const int i = state.get_i_batch();
      draws_batch.col(i) = sample;
      if (i >= 1) {
        count_acceptance += 1e-8 < arma::sum(arma::square((draws_batch.col(i - 1) - sample) / (arma::abs(draws_batch.col(i - 1)) + 1e-8)));
      }
      if (i == batch_size - 1) {
        if (count_acceptance > 1) {  // 1 acceptance seems still too few
          // TODO P(scale changes) -> 0
          const double rate_acceptance = count_acceptance / (batch_size - 1.);
          if (rate_acceptance < 0.05) {
            // no update_covariance
            scale *= .1;
          } else if (rate_acceptance > 0.95) {
            // no update_covariance
            scale *= 10.;
          } else {
            updated_proposal = state.update_covariance(draws_batch, gamma);
            if (rate_acceptance < 0.1) {  // [0.05, 0.1]
              scale *= .75;
            } else if (rate_acceptance < target_acceptance) {  // [0.1, target_acceptance]
              scale *= 0.95;
            } else if (rate_acceptance < 0.8) {  // [target_acceptance, 0.8]
              scale *= 1.05;
            } else {  // [0.8, 0.95]
              scale *= 1.5;
            }
          }
        } else if (gamma > C / 1024.) {  // no such updates after a point
          scale *= .01;
        }
        store_statistics();
        gamma = step_gamma(gamma);
        count_acceptance = 0;
      }
    }

    inline
    arma::mat get_storage () const {
      arma::mat storage(memory.size(), 3);
      for (int i = 0; i < storage.n_rows; i++) {
        storage.row(i) = arma::rowvec({memory[i].gamma, memory[i].scale, memory[i].rate_acceptance});
      }
      return storage;
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
      bool update_covariance(arma::mat::fixed<dim, batch_size> draws, const double gamma) {
        draws.each_col() -= mu;
        mu += gamma * (arma::sum(draws, 1) / batch_size - mu);
        Sigma += gamma * (draws * draws.t() / (batch_size - 1) - Sigma);
        return true;
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

    std::vector<Storage> memory;

    const double target_acceptance;
    const double lambda;  // controls the speed of adaptation
    const double alpha;  // log-speed of adaptation, derived from lambda
    const double C = 0.99;  // constant factor in the speed of adaptation
    double gamma = C;  // initialize gamma
    double scale = 1.;
    int count_acceptance = 0;

    inline
    void store_statistics () {
      memory.push_back({gamma, scale, count_acceptance / (batch_size - 1.)});
    }

    inline
    static double calculate_alpha (const double l) {
      const double lambda_reciprocal_value = 1 / (1 + l);
      return lambda_reciprocal_value + (1 - lambda_reciprocal_value) / 64;
    }
  };
}

#endif

