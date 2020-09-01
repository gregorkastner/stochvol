#ifndef _TYPE_DEFINITIONS_H_
#define _TYPE_DEFINITIONS_H_

#include <RcppArmadillo.h>
#include <vector>

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
      double mean, stdev;
      Normal (const double _m, const double _s) : mean{_m}, stdev{_s} {}
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
      enum {CONSTANT, EXPONENTIAL, INIFINITY} distribution;
      union {
        Constant constant;
        Exponential exponential;
      };

      Nu (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Nu (const Exponential& _e) : distribution{EXPONENTIAL}, exponential{_e} {}
      Nu (const Infinity& _i) : distribution{INIFINITY} {}
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
      enum {CONSTANT, NORMAL} distribution;
      union {
        Constant constant;
        Normal normal;
      };

      Covariates (const Constant& _c) : distribution{CONSTANT}, constant{_c} {}
      Covariates (const Normal& _n) : distribution{NORMAL}, normal{_n} {}
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
        const Covariates& _b = Covariates{Normal(0, 100)})
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
      const bool success = arma::inv_sympd(precision, _covariance) &&
        arma::chol(covariance_chol, _covariance, "lower") &&
        arma::inv(covariance_chol_inv, arma::trimatl(covariance_chol));
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

  struct ExpertSpec_VanillaSV {
    enum class ProposalPhi {IMMEDIATE_ACCEPT_REJECT_NORMAL, REPEATED_ACCEPT_REJECT_NORMAL};  // scoped enums would not be needed in C++14
    enum class ProposalSigma2 {INDEPENDENCE, LOG_RANDOM_WALK};  // scoped enums would not be needed in C++14

    bool interweave;
    Parameterization baseline;
    double proposal_intercept_varinv,  // B011inv
           proposal_phi_varinv;  // B022inv
    int mh_blocking_steps;  // MHsteps
    ProposalSigma2 proposal_sigma2;  // MHcontrol
    double proposal_sigma2_rw_scale;  // MHcontrol
    ProposalPhi proposal_phi;  // truncnormal

    // This constructor is implicitly well defined with initializer lists from C++14 on
    ExpertSpec_VanillaSV(
        const bool _interweave = true,
        const Parameterization _baseline = Parameterization::CENTERED,
        const double _proposal_intercept_varinv = 1e-12,
        const double _proposal_phi_varinv = 1e-8,
        const int _mh_blocking_steps = 2,
        const ProposalSigma2 _proposal_sigma2 = ProposalSigma2::INDEPENDENCE,
        const double _proposal_sigma2_rw_scale = 0.1,
        const ProposalPhi _proposal_phi = ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL)
      : interweave {_interweave},
        baseline {_baseline},
        proposal_intercept_varinv {_proposal_intercept_varinv},
        proposal_phi_varinv {_proposal_phi_varinv},
        mh_blocking_steps {_mh_blocking_steps},
        proposal_sigma2 {_proposal_sigma2},
        proposal_sigma2_rw_scale {_proposal_sigma2_rw_scale},
        proposal_phi {_proposal_phi} {}
  };

  struct ExpertSpec_GeneralSV {
    using StrategyVector = std::vector<Parameterization>;

    enum class ProposalPara {RANDOM_WALK, METROPOLIS_ADJUSTED_LANGEVIN_ALGORITHM};

    StrategyVector strategy;
    bool correct_latent_draws;
    ProposalPara proposal_para;
    bool adapt;
    ProposalDiffusionKen proposal_diffusion_ken;

    ExpertSpec_GeneralSV(
        const StrategyVector& _strategy = {Parameterization::CENTERED, Parameterization::NONCENTERED},
        const bool _correct_latent_draws = false,
        const ProposalPara _proposal_para = ProposalPara::RANDOM_WALK,
        const bool _adapt = true,
        const ProposalDiffusionKen& _proposal_diffusion_ken = {0, {}})
      : strategy {_strategy},
        correct_latent_draws {_correct_latent_draws},
        proposal_para {_proposal_para},
        adapt {_adapt},
        proposal_diffusion_ken {_proposal_diffusion_ken} {}
  };

}

#endif  // _TYPE_DEFINITIONS_H_

