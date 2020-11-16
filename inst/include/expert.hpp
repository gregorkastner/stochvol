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

#ifndef _STOCHVOL_EXPERT_H_
#define _STOCHVOL_EXPERT_H_

#include <RcppArmadillo.h>
#include <vector>
#include "type_definitions.hpp"

namespace stochvol {

  struct ExpertSpec_FastSV {
    struct Update {
      bool mixture_indicators,
           latent_vector,
           parameters;
    };

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
    Update update;

    // This constructor is implicitly well defined with initializer lists from C++14 on
    ExpertSpec_FastSV(
        const bool _interweave = true,
        const Parameterization _baseline = Parameterization::CENTERED,
        const double _proposal_intercept_varinv = 1e-12,
        const double _proposal_phi_varinv = 1e-8,
        const int _mh_blocking_steps = 2,
        const ProposalSigma2 _proposal_sigma2 = ProposalSigma2::INDEPENDENCE,
        const double _proposal_sigma2_rw_scale = 0.1,
        const ProposalPhi _proposal_phi = ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL,
        const Update& _update = {true, true, true})
      : interweave {_interweave},
        baseline {_baseline},
        proposal_intercept_varinv {_proposal_intercept_varinv},
        proposal_phi_varinv {_proposal_phi_varinv},
        mh_blocking_steps {_mh_blocking_steps},
        proposal_sigma2 {_proposal_sigma2},
        proposal_sigma2_rw_scale {_proposal_sigma2_rw_scale},
        proposal_phi {_proposal_phi},
        update {_update} {}
  };

  struct ExpertSpec_GeneralSV {
    using StrategyVector = std::vector<Parameterization>;

    struct Update {
      bool latent_vector,
           parameters;
    };

    enum class ProposalPara {RANDOM_WALK};  //, METROPOLIS_ADJUSTED_LANGEVIN_ALGORITHM};

    StrategyVector strategy;
    bool correct_latent_draws;
    ProposalPara proposal_para;
    bool adapt;
    ProposalDiffusionKen proposal_diffusion_ken;
    Update update;

    ExpertSpec_GeneralSV(
        const StrategyVector& _strategy = {Parameterization::CENTERED, Parameterization::NONCENTERED},
        const bool _correct_latent_draws = false,
        const ProposalPara _proposal_para = ProposalPara::RANDOM_WALK,
        const bool _adapt = true,
        const ProposalDiffusionKen& _proposal_diffusion_ken = {0, {}},
        const Update& _update = {true, true})
      : strategy {_strategy},
        correct_latent_draws {_correct_latent_draws},
        proposal_para {_proposal_para},
        adapt {_adapt},
        proposal_diffusion_ken {_proposal_diffusion_ken},
        update {_update} {}
  };

}

#endif

