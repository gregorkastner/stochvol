#ifndef stochvol_H_
#define stochvol_H_

#include <RcppArmadillo.h>
#include "adaptation.hpp"
#include "type_definitions.h"

namespace stochvol {

    inline
    void update_fast_sv (const arma::vec& log_data2, double& mu, double& phi, double& sigma2, double& h0, arma::vec& h, arma::ivec& r, const PriorSpec& prior_spec, const ExpertSpec_FastSV& expert) {
        typedef void(*CppFunction)(const arma::vec&, double&, double&, double&, double&, arma::vec&, arma::ivec&, const PriorSpec&, const ExpertSpec_FastSV&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_fast_sv");
        }
        {
            cpp_function(log_data2, mu, phi, sigma2, h0, h, r, prior_spec, expert);
        }
    }

    inline
    void update_general_sv (const arma::vec& data, const arma::vec& log_data2, const arma::ivec& sign_data, double& mu, double& phi, double& sigma2, double& rho, double& h0, arma::vec& h, AdaptationCollection& adaptation, const PriorSpec& prior_spec, const ExpertSpec_GeneralSV& expert) {
        typedef void(*CppFunction)(const arma::vec&, const arma::vec&, const arma::ivec&, double&, double&, double&, double&, double&, arma::vec&, AdaptationCollection&, const PriorSpec&, const ExpertSpec_GeneralSV&);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_general_sv");
        }
        {
            cpp_function(data, log_data2, sign_data, mu, phi, sigma2, rho, h0, h, adaptation, prior_spec, expert);
        }
    }


    inline
    void update_sv [[gnu::deprecated]] (const arma::vec& data, arma::vec& curpara, arma::vec& h, double& h0, arma::vec& mixprob, arma::ivec& r, const bool centered_baseline, const double C0, const double cT, const double Bsigma, const double a0, const double b0, const double bmu, const double Bmu, const double B011inv, const double B022inv, const bool Gammaprior, const bool truncnormal, const double MHcontrol, const int MHsteps, const int parameterization, const bool dontupdatemu, const double priorlatent0) {
        typedef void(*CppFunction)(const arma::vec&, arma::vec&, arma::vec&, double&, arma::vec&, arma::ivec&, const bool, const double, const double, const double, const double, const double, const double, const double, const double, const double, const bool, const bool, const double, const int, const int, const bool, const double);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_sv");
        }
        {
            cpp_function(data, curpara, h, h0, mixprob, r, centered_baseline, C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);
        }
    }

    inline
    void update_svl [[gnu::deprecated]] (const arma::vec& y, const arma::vec& y_star, const arma::ivec& d, double& phi, double& rho, double& sigma2, double& mu, double& h0, arma::vec& h, arma::vec& ht, Adaptation& adaptation_proposal, const arma::vec& prior_phi, const arma::vec& prior_rho, const arma::vec& prior_sigma2, const arma::vec& prior_mu, const bool use_mala, const bool gammaprior, const bool correct, const arma::ivec& strategy, const bool dontupdatemu) {
        typedef void(*CppFunction)(const arma::vec&, const arma::vec&, const arma::ivec&, double&, double&, double&, double&, double&, arma::vec&, arma::vec&, Adaptation&, const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&, const bool, const bool, const bool, const arma::ivec&, const bool);
        static CppFunction cpp_function = NULL;
        if (cpp_function == NULL) {
            cpp_function = (CppFunction)R_GetCCallable("stochvol", "update_svl");
        }
        {
            cpp_function(y, y_star, d, phi, rho, sigma2, mu, h0, h, ht, adaptation_proposal, prior_phi, prior_rho, prior_sigma2, prior_mu, use_mala, gammaprior, correct, strategy, dontupdatemu);
        }
    }

}

#endif // stochvol_H_
