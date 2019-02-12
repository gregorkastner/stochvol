#ifndef stochvol_H_
#define stochvol_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace stochvol {

    using namespace Rcpp;

    inline void update_sv(const arma::vec& data, arma::vec& curpara, arma::vec& h, double& h0, arma::vec& mixprob, arma::ivec& r, const bool centered_baseline, const double C0, const double cT, const double Bsigma, const double a0, const double b0, const double bmu, const double Bmu, const double B011inv, const double B022inv, const bool Gammaprior, const bool truncnormal, const double MHcontrol, const int MHsteps, const int parameterization, const bool dontupdatemu, const double priorlatent0) {
        typedef void(*Update_sv)(const arma::vec&, arma::vec&, arma::vec&, double&, arma::vec&, arma::ivec&, const bool, const double, const double, const double, const double, const double, const double, const double, const double, const double, const bool, const bool, const double, const int, const int, const bool, const double);
        static Update_sv p_update_sv = NULL;
        if (p_update_sv == NULL) {
            //validateSignature("void(*update_sv)(const arma::vec&,arma::vec&,arma::vec&,double&,arma::vec&,arma::ivec&,const bool,const double,const double,const double,const double,const double,const double,const double,const double,const double,const bool,const bool,const double,const int,const int,const bool,const double)");
            p_update_sv = (Update_sv)R_GetCCallable("stochvol", "update_sv");
        }
        {
            p_update_sv(data, curpara, h, h0, mixprob, r, centered_baseline, C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);
        }
    }

    inline void update_svl(const arma::vec& y, const arma::vec& y_star, const arma::ivec& d, double& phi, double& rho, double& sigma2, double& mu, arma::vec& h, arma::vec& ht, const arma::vec& prior_phi, const arma::vec& prior_rho, const arma::vec& prior_sigma2, const arma::vec& prior_mu, const arma::mat& proposal_chol, const arma::mat& proposal_chol_inv, const bool gammaprior, const bool correct, const arma::ivec& strategy) {
        typedef void(*Update_svl)(const arma::vec&, const arma::vec&, const arma::ivec&, double&, double&, double&, double&, arma::vec&, arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&, const arma::vec&, const arma::mat&, const arma::mat&, const bool, const bool, const arma::ivec&);
        static Update_svl p_update_svl = NULL;
        if (p_update_svl == NULL) {
            //validateSignature("void(*update_svl)(const arma::vec&,const arma::vec&,const arma::ivec&,double&,double&,double&,double&,arma::vec&,arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,const arma::vec&,const arma::mat&, const arma::mat&,const bool,const bool,const arma::ivec&)");
            p_update_svl = (Update_svl)R_GetCCallable("stochvol", "update_svl");
        }
        {
            p_update_svl(y, y_star, d, phi, rho, sigma2, mu, h, ht, prior_phi, prior_rho, prior_sigma2, prior_mu, proposal_chol, proposal_chol_inv, gammaprior, correct, strategy);
        }
    }

}

#endif // stochvol_H_
