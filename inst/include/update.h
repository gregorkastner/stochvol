#include <R_ext/Rdynload.h>
#include <Rcpp.h>

#ifndef _STOCHVOL_UPDATE_H_
#define _STOCHVOL_UPDATE_H_

inline void update(const Rcpp::NumericVector &data, double *curpara, double *h_in,
            double &h0, double *mixprob, int *r, const bool centered_baseline, const double C0,
	    const double cT, const double Bsigma, const double a0, const double b0,
	    const double bmu, const double Bmu, const double B011inv, const double B022inv,
	    const bool Gammaprior, const bool truncnormal, const double MHcontrol,
	    const int MHsteps, const int parameterization, const bool dontupdatemu,
	    const double priorlatent0) {
 
 typedef void(*Update)(const Rcpp::NumericVector &, double *, double *,
            double &, double *, int *, const bool, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const bool, const bool, const double, const int, const int,
	    const bool, const double);
 
 static Update update = (Update)R_GetCCallable("stochvol", "update");

 update(data, curpara, h_in, h0, mixprob, r, centered_baseline, C0,
             cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv, Gammaprior,
             truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu,
	     priorlatent0);
}

#endif  // _STOCHVOL_UPDATE_H_

