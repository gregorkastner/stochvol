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

#ifndef stochvol_H_
#define stochvol_H_

#include <RcppArmadillo.h>
#include "type_definitions.h"
#include "adaptation.hpp"
#include "expert.hpp"

namespace stochvol {

    inline
    void update_fast_sv (const arma::vec& log_data2, double& mu, double& phi, double& sigma, double& h0, arma::vec& h, arma::uvec& r, const PriorSpec& prior_spec, const ExpertSpec_FastSV& expert) {
        typedef void(*CppFunction)(const arma::vec&, double&, double&, double&, double&, arma::vec&, arma::uvec&, const PriorSpec&, const ExpertSpec_FastSV&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_fast_sv");
        }
        {
            cpp_function(log_data2, mu, phi, sigma, h0, h, r, prior_spec, expert);
        }
    }

    inline
    void update_general_sv (const arma::vec& data, const arma::vec& log_data2, const arma::ivec& sign_data, double& mu, double& phi, double& sigma, double& rho, double& h0, arma::vec& h, AdaptationCollection& adaptation_collection, const PriorSpec& prior_spec, const ExpertSpec_GeneralSV& expert) {
        typedef void(*CppFunction)(const arma::vec&, const arma::vec&, const arma::ivec&, double&, double&, double&, double&, double&, arma::vec&, AdaptationCollection&, const PriorSpec&, const ExpertSpec_GeneralSV&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_general_sv");
        }
        {
            cpp_function(data, log_data2, sign_data, mu, phi, sigma, rho, h0, h, adaptation_collection, prior_spec, expert);
        }
    }

    inline
    void update_t_error (const arma::vec& homosked_data, arma::vec& tau, const arma::vec& mean, const arma::vec& sd, double& nu, const PriorSpec& prior_spec, const bool do_tau_acceptance_rejection = true) {
        typedef void(*CppFunction)(const arma::vec&, arma::vec&, const arma::vec&, const arma::vec&, double&, const PriorSpec&, const bool);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_t_error");
        }
        {
            cpp_function(homosked_data, tau, mean, sd, nu, prior_spec, do_tau_acceptance_rejection);
        }
    }

    inline
    void update_regressors (const arma::vec& dependent_variable, const arma::mat& independent_variables, arma::vec& beta, const PriorSpec& prior_spec) {
        typedef void(*CppFunction)(const arma::vec&, const arma::mat&, arma::vec&, const PriorSpec&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_regressors");
        }
        {
            cpp_function(dependent_variable, independent_variables, beta, prior_spec);
        }
    }

    inline
    PriorSpec list_to_priorspec (const Rcpp::List& list) {
        typedef PriorSpec(*CppFunction)(const Rcpp::List&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "list_to_priorspec");
        }
        {
            return cpp_function(list);
        }
    }

    inline
    ExpertSpec_GeneralSV list_to_general_sv (const Rcpp::List& list, const bool correct_model_specification, const bool interweave) {
        typedef ExpertSpec_GeneralSV(*CppFunction)(const Rcpp::List&, const bool, const bool);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "list_to_general_sv");
        }
        {
            return cpp_function(list, correct_model_specification, interweave);
        }
    }

    inline
    ExpertSpec_FastSV list_to_fast_sv (const Rcpp::List& list, const bool interweave) {
        typedef ExpertSpec_FastSV(*CppFunction)(const Rcpp::List&, const bool);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "list_to_fast_sv");
        }
        {
            return cpp_function(list, interweave);
        }
    }


    inline
    void update_sv //[[gnu::deprecated]]
    (const arma::vec& data, arma::vec& curpara, arma::vec& h, double& h0, arma::vec& mixprob, arma::ivec& r, const bool centered_baseline, const double C0, const double cT, const double Bsigma, const double a0, const double b0, const double bmu, const double Bmu, const double B011inv, const double B022inv, const bool Gammaprior, const bool truncnormal, const double MHcontrol, const int MHsteps, const int parameterization, const bool dontupdatemu, const double priorlatent0) {
        typedef void(*CppFunction)(const arma::vec&, arma::vec&, arma::vec&, double&, arma::vec&, arma::ivec&, const bool, const double, const double, const double, const double, const double, const double, const double, const double, const double, const bool, const bool, const double, const int, const int, const bool, const double);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_sv");
        }
        {
            cpp_function(data, curpara, h, h0, mixprob, r, centered_baseline, C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);
        }
    }

}

#endif // stochvol_H_
