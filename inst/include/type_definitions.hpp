/*
 * R package stochvol by
 *     Gregor Kastner Copyright (C) 2013-2021
 *     Darjus Hosszejni Copyright (C) 2019-2021
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

#ifndef _TYPE_DEFINITIONS_H_
#define _TYPE_DEFINITIONS_H_

#include <RcppArmadillo.h>

namespace stochvol {

  enum class Parameterization {CENTERED, NONCENTERED};

  enum class Proposal {RWMH, MALA};

  // Class specifying the prior distributions
  // Implementation stays at the "end user functions"
  struct PriorSpec {
    // Recognized distributions
    struct Constant {
      double value;
      Constant (const double _v) : value{_v} {}
    };
    struct Normal {
      double mean, sd;
      Normal (const double _m, const double _s) : mean{_m}, sd{_s} {}
    };
    struct MultivariateNormal {
      arma::vec mean;
      arma::mat precision;
    };
    struct Gamma {
      double shape, rate;
      Gamma (const double _s, const double _r) : shape{_s}, rate{_r} {}
    };
    struct InverseGamma {
      double shape, scale;
      InverseGamma (const double _sh, const double _sc) : shape{_sh}, scale{_sc} {}
    };
    struct Beta {
      double alpha, beta;
      Beta (const double _a, const double _b) : alpha{_a}, beta{_b} {}
    };
    struct Exponential {
      double rate;
      Exponential (const double _r) : rate{_r} {}
    };
    struct Infinity {};

    // Parameter classes
    struct Mu {
      enum {CONSTANT, NORMAL} distribution;
      union {
        Constant constant;
        Normal normal;
      };

      Mu (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Mu (const Normal& _n) : distribution{NORMAL}, normal{_n} {}
    };

    struct Phi {
      enum {CONSTANT, BETA, NORMAL} distribution;
      union {
        Constant constant;
        Beta beta;
        Normal normal;
      };

      Phi (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Phi (const Beta& _b) : distribution{BETA}, beta{_b} {}
      Phi (const Normal& _n) : distribution{NORMAL}, normal{_n} {}
    };

    struct Sigma2 {
      enum {CONSTANT, GAMMA, INVERSE_GAMMA} distribution;
      union {
        Constant constant;
        Gamma gamma;
        InverseGamma inverse_gamma;
      };

      Sigma2 (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Sigma2 (const Gamma& _g) : distribution{GAMMA}, gamma{_g} {}
      Sigma2 (const InverseGamma& _i) : distribution{INVERSE_GAMMA}, inverse_gamma{_i} {}
    };

    struct Nu {
      enum {CONSTANT, EXPONENTIAL, INFINITE} distribution;
      union {
        Constant constant;
        Exponential exponential;
      };

      Nu (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Nu (const Exponential& _e) : distribution{EXPONENTIAL}, exponential{_e} {}
      Nu (const Infinity& _i) : distribution{INFINITE} {}
    };

    struct Rho {
      enum {CONSTANT, BETA} distribution;
      union {
        Constant constant;
        Beta beta;
      };

      Rho (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Rho (const Beta& _b) : distribution{BETA}, beta{_b} {}
    };

    struct Latent0 {
      enum {CONSTANT, STATIONARY} variance;
      union {
        Constant constant;
      };

      Latent0 (const Constant& _c) : variance{CONSTANT}, constant{_c} {}
      Latent0 () : variance{STATIONARY} {}  // default: stationary
    };

    struct Covariates {
      MultivariateNormal multivariate_normal;

      Covariates (const MultivariateNormal& _n) : multivariate_normal {_n.mean, _n.precision} { }
    };

    // Members
    Latent0 latent0;
    Mu mu;
    Phi phi;
    Sigma2 sigma2;
    Nu nu;
    Rho rho;
    Covariates beta;

    // This constructor is implicitly well defined with initializer lists from C++14 on
    PriorSpec(
        const Latent0& _l = Latent0 {},
        const Mu& _m = Mu {Normal(0, 100)},
        const Phi& _p = Phi {Beta(15, 1.5)},
        const Sigma2& _s = Sigma2 {Gamma(0.5, 0.5)},
        const Nu& _n = Nu {Infinity()},
        const Rho& _r = Rho {Constant(0)},
        const Covariates& _b = Covariates{MultivariateNormal{arma::zeros(1), arma::eye(1, 1)}})
      : latent0 {_l},
        mu {_m},
        phi {_p},
        sigma2 {_s},
        nu {_n},
        rho {_r},
        beta {_b} {}
  };

  // Proposal scale and covariance matrix used in the random walk and the Metropolis-adjusted Langevin algorithm
  // The English word 'ken' means range of sight.
  struct ProposalDiffusionKen {
    ProposalDiffusionKen () : ProposalDiffusionKen(0, {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}) {}
    ProposalDiffusionKen (
        const double _scale,
        const arma::mat& _covariance) {
      set(_scale, _covariance);
    }

    inline
    void set (
        const double _scale,
        const arma::mat& _covariance) {
      set_scale(_scale);
      covariance = _covariance;
      bool success = arma::inv_sympd(precision, _covariance) and
        arma::chol(covariance_chol, _covariance, "lower");
      success = success and arma::inv(covariance_chol_inv, arma::trimatl(covariance_chol));
      if (!success) {
        Rcpp::stop("Failed to take Cholesky or to take inverse");
      }
    }

    inline
    void set_scale (const double _scale) {
      scale = _scale;
    }

    inline
    double get_scale () const {
      return scale;
    }

    inline
    const arma::mat& get_covariance () const {
      return covariance;
    }

    inline
    const arma::mat& get_precision () const {
      return precision;
    }

    inline
    const arma::mat& get_covariance_chol () const {
      return covariance_chol;
    }

    inline
    const arma::mat& get_covariance_chol_inv () const {
      return covariance_chol_inv;
    }

    private:
    double scale;
    arma::mat covariance;
    arma::mat precision;  // Covariance_inv
    arma::mat covariance_chol;
    arma::mat covariance_chol_inv;
  };

}

#endif  // _TYPE_DEFINITIONS_H_

