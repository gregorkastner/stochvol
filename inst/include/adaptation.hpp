#ifndef ADAPTATION_H
#define ADAPTATION_H

#include <RcppArmadillo.h>
#include <vector>

namespace stochvol {
  // Encapsulate adaptation logic
  // Adaptation happens after every batch_size draws.
  //
  // Notation from DOI 10.1007/s11222-008-9110-y
  class Adaptation {
  public:
    // Useful structs
    struct Result {
      Result (
          const double _scale,
          const arma::mat& _covariance) {
        set_result(_scale, _covariance);
      }

      inline
      void set_result (
          const double _scale,
          const arma::mat& _covariance) {
        scale = _scale;
        covariance = _covariance;
        const bool success = arma::inv_sympd(precision, _covariance) &&
          arma::chol(covariance_chol, _covariance, "lower") &&
          arma::inv(covariance_chol_inv, arma::trimatl(covariance_chol));
        if (!success) {
          Rcpp::stop("Failed to take Cholesky or to take inverse");
        }
      }
      
      double scale;
      arma::mat covariance;
      arma::mat precision;  // Covariance_inv
      arma::mat covariance_chol;
      arma::mat covariance_chol_inv;
    };

    struct Storage {
      const double gamma;
      const double scale;
      const double rate_acceptance;
    };

    // Constructor
    Adaptation (
        const int _dim,
        const int _memory_size,
        const int _batch_size = 100,
        const double _target_acceptance = 0.234,
        const double _lambda = 0.1,
        const double _scale = 0.1,
        const double _C = 0.99)
      : dim{_dim},
        target_acceptance{_target_acceptance},
        lambda{_lambda},
        alpha{calculate_alpha(_lambda)},
        C{_C},
        scale{_scale},
        state(_dim, _batch_size),
        draws_batch(dim, _batch_size),
        cache_result(scale, arma::mat(_dim, _dim, arma::fill::eye)) {
      if (target_acceptance <= 0.1 || target_acceptance >= 0.8) {
        Rcpp::warning("Target acceptance rate should be between 10% and 80%");
      }
      memory.reserve(_memory_size);
    }
    Adaptation (const Adaptation& other)
      : dim{other.dim},
        target_acceptance{other.target_acceptance},
        lambda{other.lambda},
        alpha{other.alpha},
        C{other.C},
        scale{other.scale},
        state(other.state.dim, other.state.batch_size),
        draws_batch(arma::size(other.draws_batch)),
        cache_result(other.scale, arma::mat(arma::size(other.cache_result.covariance), arma::fill::eye)) {
      memory.reserve(other.memory.capacity());
    }

    // Calculate next gamma based on the current one
    inline
    double step_gamma (const double gamma_previous) const {
      const long i_previous = std::lround(std::pow(C / gamma_previous, 1.0 / alpha));
      return C * std::pow(i_previous + 1, -alpha);
    }

    // Store new samples in the batch state
    // Do all kinds of adaptation logic
    inline
    void register_sample (const arma::vec& sample) {
      const int i = state.get_i_batch();
      draws_batch.col(i) = sample;
      if (i >= 1) {
        count_acceptance += 1e-20 < arma::sum(arma::abs(draws_batch.col(i - 1) - sample));  // is this sophisticated enough?
      }
      if (i == state.batch_size - 1) {
        store_statistics();
        if (count_acceptance > 1) {  // 1 acceptance seems still too few
          const auto randomize_1_or = [] (const double value, const double probability_of_value) -> double {
            return (probability_of_value >= 1. || R::unif_rand() < probability_of_value) ? value : 1.;
          };
          const double probability_of_change = 100. * gamma / C;
          const double rate_acceptance = count_acceptance / (state.batch_size - 1.);
          if (rate_acceptance < 0.05) {
            // no update_covariance
            scale *= randomize_1_or(.1, probability_of_change);
          } else if (rate_acceptance > 0.95) {
            // no update_covariance
            scale *= randomize_1_or(10., probability_of_change);
          } else {
            updated_proposal = state.update_covariance(draws_batch, gamma);
            if (rate_acceptance < 0.1) {  // [0.05, 0.1]
              scale *= randomize_1_or(.75, probability_of_change);
            } else if (rate_acceptance < target_acceptance) {  // [0.1, target_acceptance]
              scale *= randomize_1_or(0.95, probability_of_change);
            } else if (rate_acceptance < 0.8) {  // [target_acceptance, 0.8]
              scale *= randomize_1_or(1.05, probability_of_change);
            } else {  // [0.8, 0.95]
              scale *= randomize_1_or(1.5, probability_of_change);
            }
          }
        } else if (gamma > C * 0.001) {  // no such updates after a point
          scale *= .01;
        }
        gamma = step_gamma(gamma);
        count_acceptance = 0;
      }
    }

    inline
    Rcpp::NumericMatrix get_storage () const {
      const Rcpp::CharacterVector coln({"Gamma", "Scale", "Acceptance Rate"});
      Rcpp::NumericMatrix storage(memory.size(), 3);
      for (int i = 0; i < storage.nrow(); i++) {
        storage(i, Rcpp::_) = Rcpp::NumericVector({memory[i].gamma, memory[i].scale, memory[i].rate_acceptance});
      }
      colnames(storage) = coln;
      return storage;
    }

    inline
    const Result& get_proposal () {
      if (updated_proposal) {
        updated_proposal = false;
        cache_result.set_result(scale, state.get_covariance());
      } else {
        cache_result.scale = scale;
      }
      return cache_result;
    }

  private:
    class State {
    public:
      const int batch_size;
      const int dim;

      State (const int _dim, const int _batch_size)
        : batch_size{_batch_size}, dim{_dim} {
        mu.zeros(dim);
        Sigma.zeros(dim, dim);
      }

      inline
      int get_i_batch () {
        const int i_batch_current = i_batch;
        i_batch = (i_batch + 1) % batch_size;
        return i_batch_current;
      }

      // overwrites draws
      inline
      bool update_covariance(arma::mat draws, const double gamma) {  // TODO maybe prevent copy of 'draws'
        draws.each_col() -= mu;
        mu += gamma * (arma::sum(draws, 1) / batch_size - mu);
        Sigma += gamma * (draws * draws.t() / (batch_size - 1) - Sigma);
        return true;
      }

      inline
      const arma::mat& get_covariance () const {
        return Sigma;
      }

    private:
      int i_batch = 0;  // use get_i_batch() to access i_batch
      arma::vec mu;
      arma::mat Sigma;
    };

    const int dim;
    const double target_acceptance;
    const double lambda;  // controls the speed of adaptation
    const double alpha;  // log-speed of adaptation, derived from lambda
    const double C;  // constant factor in the speed of adaptation
    double gamma = C;  // initialize gamma
    double scale;
    int count_acceptance = 0;

    State state;
    arma::mat draws_batch;
    bool updated_proposal = false;

    std::vector<Storage> memory;
    Result cache_result;

    inline
    void store_statistics () {
      memory.push_back({gamma, scale, count_acceptance / (state.batch_size - 1.)});
    }

    inline
    static double calculate_alpha (const double l) {
      const double lambda_reciprocal_value = 1 / (1 + l);
      return lambda_reciprocal_value + (1 - lambda_reciprocal_value) / 64;
    }
  };

  // Convenience class
  struct AdaptationCollection {
    Adaptation centered;
    Adaptation noncentered;

    AdaptationCollection (
        const int _dim,
        const int _memory_size,
        const int _batch_size = 100,
        const double _target_acceptance = 0.234,
        const double _lambda = 0.1,
        const double _scale = 0.1,
        const double _C = 0.99)
      : centered(_dim, _memory_size, _batch_size, _target_acceptance, _lambda, _scale, _C),
        noncentered(centered) { }
  };
}

#endif

