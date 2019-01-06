svlsample_cpp <- function(draws, y, y_star, d, burnin, thinpara, thinlatent, thintime, phi_init, rho_init, sigma2_init, mu_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, verbose, stdev, gammaprior, strategy, mixing_constants) {
    .Call("_stochvol_svlsample_cpp", draws, y, y_star, d, burnin, thinpara, thinlatent, thintime, phi_init, rho_init, sigma2_init, mu_init, h_init, prior_phi_a, prior_phi_b, prior_rho_a, prior_rho_b, prior_sigma2_shape, prior_sigma2_rate, prior_mu_mu, prior_mu_sigma, verbose, stdev, gammaprior, strategy, mixing_constants, PACKAGE = "stochvol")
}

sampler <- function(y, draws, burnin, designmatrix,
          priormu1, priormu2, priorphi1, priorphi2, priorsigma, 
          thinlatent, thintime, startpara, startlatent, keeptau, myquiet, para,
          mhsteps, B011, B022, mhcontrol, gammaprior, truncnormal,
          myoffset, dontupdatemu, priornu, priorbeta, priorlatent0) {
  .Call("_stochvol_sampler", y, draws, burnin, designmatrix,
          priormu1, priormu2, priorphi1, priorphi2, priorsigma, 
          thinlatent, thintime, startpara, startlatent, keeptau, myquiet, para,
          mhsteps, B011, B022, mhcontrol, gammaprior, truncnormal,
          myoffset, dontupdatemu, priornu, priorbeta, priorlatent0, PACKAGE = "stochvol")
}


