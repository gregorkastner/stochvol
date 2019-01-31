svsample_cpp <- function(y_in, draws, burnin, X_in, bmu, Bmu, a0, b0, Bsigma, thin, timethin, startpara_in, startvol_in, keeptau, quiet, para, MHsteps, B011, B022, mhcontrol, gammaprior, truncnormal, offset, dontupdatemu, priordf_in, priorbeta_in, priorlatent0) {
    .Call(`_stochvol_svsample_cpp`, y_in, draws, burnin, X_in, bmu, Bmu, a0, b0, Bsigma, thin, timethin, startpara_in, startvol_in, keeptau, quiet, para, MHsteps, B011, B022, mhcontrol, gammaprior, truncnormal, offset, dontupdatemu, priordf_in, priorbeta_in, priorlatent0, PACKAGE = "stochvol")
}

svlsample_cpp <- function(y, draws, burnin, X, thinpara, thinlatent, thintime, theta_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, prior_beta_mu, prior_beta_sigma, verbose, offset, stdev, gammaprior, correct, strategy) {
    .Call(`_stochvol_svlsample_cpp`, y, draws, burnin, X, thinpara, thinlatent, thintime, theta_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, prior_beta_mu, prior_beta_sigma, verbose, offset, stdev, gammaprior, correct, strategy, PACKAGE = "stochvol")
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_stochvol_RcppExport_registerCCallable', PACKAGE = 'stochvol')
})
