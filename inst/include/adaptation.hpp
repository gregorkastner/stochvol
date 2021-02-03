/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2020
 *     Darjus Hosszejni Copyright (C) 2019-2020
 *  
 *  This file is part of the R package stochvol: Efficient Bayesian
 *  Inference for Stochastic Volatility Models.
 *  
 *  The R package stochvol is free software: you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 2 or
 *  any later version of the License.
 *  
 *  The R package stochvol is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with the R package stochvol. If that is not the case, please
 *  refer to <http://www.gnu.org/licenses/>.
 */

#ifndef _STOCHVOL_ADAPTATION_H_
#define _STOCHVOL_ADAPTATION_H_

#include <RcppArmadillo.h>
#include <vector>
#include "type_definitions.hpp"

namespace stochvol {
  // Encapsulate adaptation logic
  // Adaptation happens after every batch_size draws.
  //
  // Notation from DOI 10.1007/s11222-008-9110-y
  class Adaptation {
  public:
    // Useful structs
    struct Storage {
      double gamma;
      double scale;
      double rate_acceptance;
    };

    // Constructor / initializes from zero knowledge
    Adaptation (
        const int _dim,
        const int _memory_size,
        const int _batch_size = 100,
        const double _target_acceptance = 0.234,
        const double _lambda = 0.1,
        const double _scale = 0.1,
        const double _C_in = 0.99)  // _C is reserved name in the C++ standard
      : target_acceptance{_target_acceptance},
        lambda{_lambda},
        alpha{calculate_alpha(_lambda)},
        C{_C_in},
        scale{_scale},
        state(_dim, _batch_size),
        draws_batch(_dim, _batch_size),
        cache_result(scale, arma::mat(_dim, _dim, arma::fill::eye)) {
      if (target_acceptance <= 0.1 || target_acceptance >= 0.8) {
        Rcpp::warning("Target acceptance rate should be between 10% and 80%");
      }
      memory.reserve(_memory_size);
    }

    // Empty object
    Adaptation ()
      : Adaptation(0, 0) {}

    // "Exact" constructor sets all fields
    // use move semantics as much as possible
    Adaptation (
        const int _dim,
        std::vector<Storage>&& _memory,
        const int _batch_size,
        const double _target_acceptance,
        const double _lambda,
        const double _scale,
        const double _C_in,
        const double _alpha,
        const double _gamma,
        const int _count_acceptance,
        const int _i_batch,
        arma::vec&& _mu,
        arma::mat&& _Sigma,
        arma::mat&& _draws_batch,
        const bool _updated_proposal,
        const double _cached_scale,
        const arma::mat& _cached_covariance)
      : target_acceptance{_target_acceptance}, lambda{_lambda},
        alpha{_alpha}, C{_C_in}, gamma{_gamma}, scale{_scale},
        count_acceptance{_count_acceptance},
        state(_dim, _batch_size, _i_batch, std::move(_mu), std::move(_Sigma)),
        draws_batch(std::move(_draws_batch)), updated_proposal{_updated_proposal},
        memory(std::move(_memory)),
        cache_result(_cached_scale, _cached_covariance) {}

    // "Exact" constructor sets all fields
    // but uses Armadillo memory re-usage instead of move
    Adaptation (
        const int _dim,
        std::vector<Storage>&& _memory,
        const int _batch_size,
        const double _target_acceptance,
        const double _lambda,
        const double _scale,
        const double _C_in,
        const double _alpha,
        const double _gamma,
        const int _count_acceptance,
        const int _i_batch,
        arma::vec& _mu,
        arma::mat& _Sigma,
        arma::mat& _draws_batch,
        const bool _updated_proposal,
        const double _cached_scale,
        const arma::mat& _cached_covariance)
      : target_acceptance{_target_acceptance}, lambda{_lambda},
        alpha{_alpha}, C{_C_in}, gamma{_gamma}, scale{_scale},
        count_acceptance{_count_acceptance},
        state(_dim, _batch_size, _i_batch, _mu, _Sigma),
        draws_batch(_draws_batch.begin(), _draws_batch.n_rows, _draws_batch.n_cols, false),
        updated_proposal{_updated_proposal},
        memory(std::move(_memory)),
        cache_result(_cached_scale, _cached_covariance) {}

    // "Exact" constructor sets all fields
    // but copies everything
    Adaptation (
        const int _dim,
        const std::vector<Storage>& _memory,
        const int _batch_size,
        const double _target_acceptance,
        const double _lambda,
        const double _scale,
        const double _C_in,
        const double _alpha,
        const double _gamma,
        const int _count_acceptance,
        const int _i_batch,
        const arma::vec& _mu,
        const arma::mat& _Sigma,
        const arma::mat& _draws_batch,
        const bool _updated_proposal,
        const double _cached_scale,
        const arma::mat& _cached_covariance)
      : target_acceptance{_target_acceptance}, lambda{_lambda},
        alpha{_alpha}, C{_C_in}, gamma{_gamma}, scale{_scale},
        count_acceptance{_count_acceptance},
        state(_dim, _batch_size, _i_batch, _mu, _Sigma),
        draws_batch{_draws_batch},
        updated_proposal{_updated_proposal},
        memory(_memory),
        cache_result(_cached_scale, _cached_covariance) {
      memory.reserve(_memory.capacity());
    }

    // copy constructor
    Adaptation (const Adaptation& other)
      : Adaptation(
          other.draws_batch.n_rows,
          other.memory,
          other.state.get_batch_size(),
          other.target_acceptance,
          other.lambda,
          other.scale,
          other.C,
          other.alpha,
          other.gamma,
          other.count_acceptance,
          other.state.get_i_batch(),
          other.state.get_mean(),
          other.state.get_covariance(),
          other.draws_batch,
          other.updated_proposal,
          other.cache_result.get_scale(),
          other.cache_result.get_covariance()) {}

    /*
    // move constructor
    Adaptation (Adaptation&& other)
      : Adaptation(
          other.draws_batch.n_rows,
          std::move(other.memory),
          other.state.get_batch_size(),
          other.target_acceptance,
          other.lambda,
          other.scale,
          other.C,
          other.alpha,
          other.gamma,
          other.count_acceptance,
          other.state.get_i_batch(),
          std::move(other.state.get_mean()),
          std::move(other.state.get_covariance()),
          std::move(other.draws_batch),
          other.updated_proposal,
          other.cache_result.get_scale(),
          other.cache_result.get_covariance()) {}  // this can't be moved as of now
          */

    // Calculate next gamma based on the current one
    inline
    double step_gamma (const double gamma_previous) const {
      const long i_previous = std::lround(std::pow(C / gamma_previous, 1.0 / alpha));
      return C * std::pow(i_previous + 1, -alpha);
    }

    // Store new samples in the batch state
    // Do all kinds of adaptation logic
    inline
    void register_sample (const bool accepted, const arma::vec& sample) {
      const int i = state.get_and_increment_i_batch();
      draws_batch.col(i) = sample;
      count_acceptance += accepted;
      if (i == state.get_batch_size() - 1) {
        store_statistics();
        if (count_acceptance > 1) {  // 1 acceptance seems still too few
          const auto randomize_1_or = [] (const double value, const double probability_of_value) -> double {
            return (probability_of_value >= 1. or R::unif_rand() < probability_of_value) ? value : 1.;
          };
          const double probability_of_change = 100. * gamma / C;
          const double relative_acceptance = count_acceptance / (target_acceptance * state.get_batch_size());
          if (relative_acceptance < 0.2) {  // less than 20% of target_acceptance
            // no update_covariance
            scale *= randomize_1_or(.1, probability_of_change);
          } else if (relative_acceptance > 3) {  // more than 3x target_acceptance
            // no update_covariance
            scale *= randomize_1_or(10., probability_of_change);
          } else {
            updated_proposal = state.update_covariance(draws_batch, gamma);
            if (relative_acceptance < 0.5) {  // [20%, 50%] of target_acceptance
              scale *= randomize_1_or(.75, probability_of_change);
            } else if (relative_acceptance < 1) {  // [50%, 100%] of target_acceptance
              scale *= randomize_1_or(0.95, probability_of_change);
            } else if (relative_acceptance < 0.8) {  // [100%, 180%] of target_acceptance
              scale *= randomize_1_or(1.05, probability_of_change);
            } else {  // [180%, 300%] of target_acceptance
              scale *= randomize_1_or(1.5, probability_of_change);
            }
          }
        } else if (gamma > C * 0.001) {  // if acceptance is too low (but no such updates after a point)
          scale *= .01;
        }
        gamma = step_gamma(gamma);
        count_acceptance = 0;
      }
    }

    inline
    const ProposalDiffusionKen& get_proposal () {
      if (updated_proposal) {
        updated_proposal = false;
        cache_result.set(scale, state.get_covariance());
      } else {
        cache_result.set_scale(scale);
      }
      return cache_result;
    }

    // Functions for communication with R
    inline
    Rcpp::NumericMatrix get_storage () const {
      const Rcpp::CharacterVector coln({"Gamma", "Scale", "Acceptance Rate"});
      Rcpp::NumericMatrix storage(memory.capacity(), 3);
      storage.fill(::R_NaReal);
      for (unsigned int i = 0; i < memory.size(); i++) {
        storage(i, 0) = memory[i].gamma;
        storage(i, 1) = memory[i].scale;
        storage(i, 2) = memory[i].rate_acceptance;
      }
      colnames(storage) = coln;
      return storage;
    }
    
    inline
    Rcpp::List serialize () const {
      return Rcpp::List::create(
        Rcpp::_["dim"] = Rcpp::wrap(draws_batch.n_rows),
        Rcpp::_["memory"] = get_storage(),
        Rcpp::_["batch_size"] = Rcpp::wrap(state.get_batch_size()),
        Rcpp::_["target_acceptance"] = Rcpp::wrap(target_acceptance),
        Rcpp::_["lambda"] = Rcpp::wrap(lambda),
        Rcpp::_["scale"] = Rcpp::wrap(scale),
        Rcpp::_["C"] = Rcpp::wrap(C),
        Rcpp::_["alpha"] = Rcpp::wrap(alpha),
        Rcpp::_["gamma"] = Rcpp::wrap(gamma),
        Rcpp::_["count_acceptance"] = Rcpp::wrap(count_acceptance),
        Rcpp::_["i_batch"] = Rcpp::wrap(state.get_i_batch()),
        Rcpp::_["mu"] = Rcpp::wrap(state.get_mean()),
        Rcpp::_["Sigma"] = Rcpp::wrap(state.get_covariance()),
        Rcpp::_["draws_batch"] = Rcpp::wrap(draws_batch),
        Rcpp::_["updated_proposal"] = Rcpp::wrap(updated_proposal),
        Rcpp::_["cached_scale"] = Rcpp::wrap(cache_result.get_scale()),
        Rcpp::_["cached_covariance"] = Rcpp::wrap(cache_result.get_covariance()));
    }

  private:
    class State {
    public:
      State (const int _dim, const int _batch_size)
        : batch_size{_batch_size},
          mu(_dim, arma::fill::zeros),
          Sigma(_dim, _dim, arma::fill::eye) {
      }

      // constructor not meant to be used others than stochvol developers
      // move semantics
      State (const int _dim, const int _batch_size,
          const int _i_batch, arma::vec&& _mu,
          arma::mat&& _Sigma)
        : batch_size{_batch_size},
          i_batch{_i_batch}, mu(std::move(_mu)),
          Sigma(std::move(_Sigma)) {}

      // constructor not meant to be used others than stochvol developers
      // Armadillo memory re-usage
      State (const int _dim, const int _batch_size,
          const int _i_batch, arma::vec& _mu,
          arma::mat& _Sigma)
        : batch_size{_batch_size},
          i_batch{_i_batch},
          mu(_mu.begin(), _mu.n_elem, false),
          Sigma(_Sigma.begin(), _Sigma.n_rows, _Sigma.n_cols, false) {}

      // constructor not meant to be used others than stochvol developers
      // copy everything
      State (const int _dim, const int _batch_size,
          const int _i_batch, const arma::vec& _mu,
          const arma::mat& _Sigma)
        : batch_size{_batch_size},
          i_batch{_i_batch},
          mu{_mu},
          Sigma{_Sigma} {}

      inline
      int get_and_increment_i_batch () {
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
      int get_i_batch () const {
        return i_batch;
      }

      inline
      int get_batch_size () const {
        return batch_size;
      }

      inline
      const arma::vec& get_mean () const {
        return mu;
      }

      inline
      const arma::mat& get_covariance () const {
        return Sigma;
      }

      inline
      arma::vec&& get_mean () {
        return std::move(mu);
      }

      inline
      arma::mat&& get_covariance () {
        return std::move(Sigma);
      }

    private:
      int batch_size;
      int i_batch = 0;  // use get_and_increment_i_batch() or get_i_batch() to access i_batch
      arma::vec mu;
      arma::mat Sigma;
    };

    double target_acceptance;  // target acceptance rate
    const double lambda;  // controls the speed of adaptation
    double alpha;  // log-speed of adaptation, derived from lambda
    const double C;  // constant factor in the speed of adaptation
    double gamma = C;  // initialize gamma
    double scale;  // current scaling for the covariance matrix
    int count_acceptance = 0;  // count the number of accepted proposals in this batch

    State state;  // current mean vector and covariance matrix for the random walk
    arma::mat draws_batch;  // storage for all the draws in the current batch
    bool updated_proposal = false;  // is the current cache up to date?

    std::vector<Storage> memory;  // store gamma, scale, and count_acceptance
    ProposalDiffusionKen cache_result;  // cache the current scale and covariance matrix (and its inverse, cholesky decomposition)

    // store gamma, scale, and the acceptance rate (computed from count_acceptance)
    inline
    void store_statistics () {
      if (memory.capacity() > memory.size()) {
        memory.push_back({gamma, scale, count_acceptance / double(state.get_batch_size())});
      }
    }

    inline
    static double calculate_alpha (const double l) {  // this code is duplicated in R/wrappers.R:get_default_adaptation(priorspec)
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
        const double _C_in = 0.99)  // _C is reserved name in the C++ standard
      : centered(_dim, _memory_size, _batch_size, _target_acceptance, _lambda, _scale, _C_in),
        noncentered(centered) { }  // this is the copy (and not the move) constructor because 'centered' is an lvalue

    AdaptationCollection ()
      : centered(), noncentered() {}

    AdaptationCollection (
        const Adaptation& _centered,
        const Adaptation& _noncentered)
      : centered(_centered), noncentered(_noncentered) {}
    /*
    AdaptationCollection (
        Adaptation&& _centered,
        Adaptation&& _noncentered)
      : centered(std::move(_centered)), noncentered(std::move(_noncentered)) {}
      */

    inline
    Adaptation& operator[] (const Parameterization par) {
      switch (par) {
        case Parameterization::CENTERED:
          return centered;
        case Parameterization::NONCENTERED:
          return noncentered;
        default:
          ::Rf_error("Adaptation operator[]: Mistake in the switch-case");
      }
    }

    inline
    const Adaptation& operator[] (const Parameterization par) const {
      switch (par) {
        case Parameterization::CENTERED:
          return centered;
        case Parameterization::NONCENTERED:
          return noncentered;
        default:
          ::Rf_error("Adaptation operator[]: Mistake in the switch-case");
      }
    }

    inline
    Rcpp::List serialize () const {
      return Rcpp::List::create(
          Rcpp::_["centered"] = centered.serialize(),
          Rcpp::_["noncentered"] = noncentered.serialize());
    }
  };
}

#endif

