#ifndef _TYPE_DEFINITIONS_H_
#define _TYPE_DEFINITIONS_H_

#include <vector>

namespace stochvol {

  enum class Parameterization {CENTERED, NONCENTERED};

  enum class Proposal {RWMH, MALA};

  // Class specifying the prior distributions
  // Implementation stays at the "end user functions"
  struct PriorSpec {
    // Recognized distributions
    struct Constant {
      const double value;
      Constant (const double _v) : value{_v} {}
    };
    struct Normal {
      const double mean, stdev;
      Normal (const double _m, const double _s) : mean{_m}, stdev{_s} {}
    };
    struct Gamma {
      const double shape, rate;
      Gamma (const double _s, const double _r) : shape{_s}, rate{_r} {}
    };
    struct InverseGamma {
      const double shape, scale;
      InverseGamma (const double _sh, const double _sc) : shape{_sh}, scale{_sc} {}
    };
    struct Beta {
      const double alpha, beta;
      Beta (const double _a, const double _b) : alpha{_a}, beta{_b} {}
    };
    struct Exponential {
      const double rate;
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

    // Members
    Latent0 latent0;
    Mu mu;
    Phi phi;
    Sigma2 sigma2;
    Nu nu;
    Rho rho;

    // This constructor is implicitly well defined with initializer lists from C++14 on
    PriorSpec(
        const Latent0& _l,
        const Mu& _m,
        const Phi& _p,
        const Sigma2& _s,
        const Nu& _n = Nu {Infinity()},
        const Rho& _r = Rho {Constant(0)})
      : latent0 {_l},
        mu {_m},
        phi {_p},
        sigma2 {_s},
        nu {_n},
        rho {_r} {}
  };

  struct ExpertSpec_VanillaSV {
    enum class ProposalPhi {ACCEPT_REJECT_NORMAL, TRUNCATED_NORMAL};  // scoped enums would not be needed in C++14
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
        const ProposalPhi _proposal_phi = ProposalPhi::ACCEPT_REJECT_NORMAL)
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

    StrategyVector strategy;
    bool correct_latent_draws,
         proposal_mala;

    ExpertSpec_GeneralSV(
        const StrategyVector& _strategy = {Parameterization::CENTERED, Parameterization::NONCENTERED},
        const bool _correct_latent_draws = false,
        const bool _proposal_mala = false)
      : strategy {_strategy},
        correct_latent_draws {_correct_latent_draws},
        proposal_mala {_proposal_mala} {}
  };

}

#endif  // _TYPE_DEFINITIONS_H_

