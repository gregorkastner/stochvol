#include <RcppArmadillo.h>
#include "single_update.h"
#include "densities.h"
#include "type_definitions.h"
#include "sampling_latent_states.h"
#include "sampling_parameters.h"
#include "utils_latent_states.h"
#include "utils_parameters.h"
#include "utils.h"

using namespace Rcpp;

namespace stochvol {

void update_vanilla_sv(
    const arma::vec& log_data2,
    double& mu,
    double& phi,
    double& sigma2,
    double& h0,
    arma::vec& h,
    arma::ivec& r,
    const PriorSpec& prior_spec,  // C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, dontupdatemu, priorlatent0 feed into this
    const ExpertSpec_VanillaSV& expert) {  // parameterization, centered_baseline, B011inv, B022inv, truncnormal, Mhcontrol, MHsteps feed into this
  // TODO setup validation (put in a "debug" environment in the end)
  // phi is beta
  // sigma2 is either inverse_gamma with shape == 2.5 or gamma with shape == 0.5
  // inverse_gamma prior and mh_blocking_steps != 2 is nyi
  // constant mu implies mh_blocking_steps == 3 (mh_blocking_steps == 1 is nyi and 2 does not make sense because mu is drawn jointly with either phi or sigma2 in the 2-block)

  const arma::vec& data = log_data2;
  const int T = data.size();
  double sigma = std::sqrt(sigma2);
  double ht0 = centered_to_noncentered(mu, sigma, h0);
  arma::vec ht = centered_to_noncentered(mu, sigma, h);

  // PriorSpec (this part will be removed at later stages of the refactorization)
  const double a0 = prior_spec.phi.beta.alpha,
               b0 = prior_spec.phi.beta.beta;
  const double bmu = prior_spec.mu.distribution == PriorSpec::Mu::NORMAL ? prior_spec.mu.normal.mean : std::numeric_limits<double>::lowest(),
               Bmu = prior_spec.mu.distribution == PriorSpec::Mu::NORMAL ? std::pow(prior_spec.mu.normal.sd, 2) : std::numeric_limits<double>::lowest();
  const bool Gammaprior = prior_spec.sigma2.distribution == PriorSpec::Sigma2::GAMMA,
             dontupdatemu = prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT;
  const double C0 = [&prior_spec]() -> double {
    switch (prior_spec.sigma2.distribution) {
      case PriorSpec::Sigma2::INVERSE_GAMMA:
        return prior_spec.sigma2.inverse_gamma.scale;
      default:
        return std::numeric_limits<double>::lowest();
    }
  }();
  const double cT = [&prior_spec, &expert, T]() -> double {
    switch (prior_spec.sigma2.distribution) {
      case PriorSpec::Sigma2::GAMMA:
        if (expert.mh_blocking_steps == 1) {
          return 0.5 * (T - 1);
        } else {
          return 0.5 * T;
        }
      case PriorSpec::Sigma2::INVERSE_GAMMA:
        if (expert.mh_blocking_steps == 2) {
          return prior_spec.sigma2.inverse_gamma.shape + 0.5 * (T + 1);
        } else {
          return std::numeric_limits<double>::lowest();
        }
      default:
        return std::numeric_limits<double>::lowest();
    }
  }();
  const double Bsigma = [&prior_spec]() -> double {
    switch (prior_spec.sigma2.distribution) {
      case PriorSpec::Sigma2::GAMMA:
        return 0.5 / prior_spec.sigma2.gamma.rate;
      default:
        return std::numeric_limits<double>::lowest();
    }
  }();
  const double priorlatent0 = [&prior_spec]() -> double {
    switch (prior_spec.latent0.variance) {
      case PriorSpec::Latent0::STATIONARY:
        return -1;
      case PriorSpec::Latent0::CONSTANT:
        return prior_spec.latent0.constant.value;
    }
  }();
  const double Bh0inv = [&prior_spec, phi]() -> double {
    switch (prior_spec.latent0.variance) {
      case PriorSpec::Latent0::STATIONARY:
        return 1. - std::pow(phi, 2);
      case PriorSpec::Latent0::CONSTANT:
        return 1. / prior_spec.latent0.constant.value;
    }
  }();

  // ExpertSpec (this part will be removed at later stages of the refactorization)
  const double B011inv = expert.proposal_intercept_varinv,
               B022inv = expert.proposal_phi_varinv;
  const bool truncnormal = expert.proposal_phi == ExpertSpec_VanillaSV::ProposalPhi::REPEATED_ACCEPT_REJECT_NORMAL;
  const int MHsteps = expert.mh_blocking_steps;
  const double MHcontrol = [&expert]() -> double {
    switch (expert.proposal_sigma2) {
      case ExpertSpec_VanillaSV::ProposalSigma2::INDEPENDENCE:
        return -1;
      case ExpertSpec_VanillaSV::ProposalSigma2::LOG_RANDOM_WALK:
        return expert.proposal_sigma2_rw_scale;
    }
  }();


  if (dontupdatemu) {
    mu = 0;
  }

  /*
   * Step (c): sample indicators
   */
  {
    // calculate non-normalized CDF of posterior indicator probs
    arma::vec mixprob(10 * T);

    find_mixture_indicator_cdf(mixprob, data-h);

    if (false) {  // TODO remove
      switch (expert.baseline) {
        case Parameterization::CENTERED:
          find_mixture_indicator_cdf(mixprob, data-h);
          break;
        case Parameterization::NONCENTERED:
          find_mixture_indicator_cdf(mixprob, data-mu-sigma*h); 
          break;
      }
    }

    // find correct indicators
    inverse_transform_sampling(mixprob, r, T);
  }

  /*
   * Step (a): sample the latent volatilities h:
   */
  {
    double omega_offdiag;  // contains off-diag element of precision matrix (const)
    arma::vec omega_diag(T+1),  // contains diagonal elements of precision matrix
              chol_offdiag(T), chol_diag(T+1),  // Cholesky-factor of Omega
              covector(T+1),  // holds covector (see McCausland et al. 2011)
              htmp(T+1),  // intermediate vector for sampling h
              hnew(T+1);  // intermediate vector for sampling h
    const double sigma2inv = 1. / sigma2;

    switch (expert.baseline) {
      case Parameterization::CENTERED:  // fill precision matrix omega and covector c for CENTERED para:
        {
          const double phi2 = std::pow(phi, 2);

          omega_diag[0] = (Bh0inv + phi2) * sigma2inv;
          covector[0] = mu * (Bh0inv - phi*(1-phi)) * sigma2inv;

          for (int j = 1; j < T; j++) {
            omega_diag[j] = mix_varinv[r[j-1]] + (1+phi2)*sigma2inv; 
            covector[j] = (data[j-1] - mix_mean[r[j-1]])*mix_varinv[r[j-1]]
              + mu*(1-phi)*(1-phi)*sigma2inv;
          }
          omega_diag[T] = mix_varinv[r[T-1]] + sigma2inv;
          covector[T] = (data[T-1] - mix_mean[r[T-1]])*mix_varinv[r[T-1]] + mu*(1-phi)*sigma2inv;
          omega_offdiag = -phi*sigma2inv;
        }
        break;
      case Parameterization::NONCENTERED:  // fill precision matrix omega and covector c for NONCENTERED para:
        {
          const double phi2 = std::pow(phi, 2);

          omega_diag[0] = Bh0inv + phi2;
          covector[0] = 0.;

          for (int j = 1; j < T; j++) {
            omega_diag[j] = mix_varinv[r[j-1]] * sigma2 + 1 + phi2;
            covector[j] = mix_varinv[r[j-1]] * sigma * (data[j-1] - mix_mean[r[j-1]] - mu);
          }
          omega_diag[T] = mix_varinv[r[T-1]] * sigma2 + 1;
          covector[T] = mix_varinv[r[T-1]] * sigma * (data[T-1] - mix_mean[r[T-1]] - mu);
          omega_offdiag = -phi;
        }
        break;
    } 

    // Cholesky decomposition
    cholesky_tridiagonal(omega_diag, omega_offdiag, chol_diag, chol_offdiag);

    // Solution of Chol*x = covector ("forward algorithm")
    forward_algorithm(chol_diag, chol_offdiag, covector, htmp);

    htmp += as<arma::vec>(rnorm(T+1));
    //htmp.transform( [](double h_elem) -> double { return h_elem + R::norm_rand(); });

    // Solution of (Chol')*x = htmp ("backward algorithm")
    backward_algorithm(chol_diag, chol_offdiag, htmp, hnew);

    switch (expert.baseline) {
      case Parameterization::CENTERED:
        h = hnew.tail(T);
        h0 = hnew[0];

        ht = centered_to_noncentered(mu, sigma, h);
        ht0 = centered_to_noncentered(mu, sigma, h0);
        break;
      case Parameterization::NONCENTERED:
        ht = hnew.tail(T);
        ht0 = hnew[0];

        h = noncentered_to_centered(mu, sigma, ht);
        h0 = noncentered_to_centered(mu, sigma, ht0);
        break;
    }
  }

  /*
   * Step (b): sample mu, phi, sigma
   */
  {
    switch (expert.baseline) {
      case Parameterization::CENTERED:
        {
          {
            const auto parameter_draw = regression_centered(h0, h, mu, phi, sigma,
                C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv,
                B022inv, Gammaprior, truncnormal, MHcontrol, MHsteps, dontupdatemu, priorlatent0);
            mu = parameter_draw.mu;
            phi = parameter_draw.phi;
            sigma = parameter_draw.sigma;
            sigma2 = std::pow(sigma, 2);

            ht = centered_to_noncentered(mu, sigma, h);
            ht0 = centered_to_noncentered(mu, sigma, h0);
          }

          if (expert.interweave) {
            const auto parameter_draw = regression_noncentered(data, ht0, ht, r,
                mu, phi, sigma, Bsigma, a0, b0, bmu, Bmu,
                truncnormal, MHsteps, dontupdatemu, priorlatent0);
            mu = parameter_draw.mu;
            phi = parameter_draw.phi;
            sigma = parameter_draw.sigma;
            sigma2 = std::pow(sigma, 2);
            h = noncentered_to_centered(mu, sigma, ht);
            h0 = noncentered_to_centered(mu, sigma, ht0);
          }
        }
        break;
      case Parameterization::NONCENTERED:
        {
          {
            const auto parameter_draw = regression_noncentered(data, ht0, ht, r, mu, phi, sigma,
                Bsigma, a0, b0, bmu, Bmu, truncnormal, MHsteps,
                dontupdatemu, priorlatent0);
            mu = parameter_draw.mu;
            phi = parameter_draw.phi;
            sigma = parameter_draw.sigma;
            sigma2 = std::pow(sigma, 2);

            h = noncentered_to_centered(mu, sigma, ht);
            h0 = noncentered_to_centered(mu, sigma, ht0);
          }

          if (expert.interweave) {
            const auto parameter_draw = regression_centered(h0, h, mu, phi, sigma,
                C0, cT, Bsigma, a0, b0, bmu, Bmu, B011inv, B022inv,
                Gammaprior, truncnormal, MHcontrol, MHsteps,
                dontupdatemu, priorlatent0);
            mu = parameter_draw.mu;
            phi = parameter_draw.phi;
            sigma = parameter_draw.sigma;
            sigma2 = std::pow(sigma, 2);
            // updating ht0 and ht is unnecessary, we don't return them
            //ht = centered_to_noncentered(mu, sigma, h);
            //ht0 = centered_to_noncentered(mu, sigma, h0);
          }
        }
    }
  }
}

void update_df_svt(
    const arma::vec& log_data2_minus_h,
    arma::vec& tau,
    double& nu,
    const PriorSpec& prior_spec) {

  R_assert(prior_spec.nu.distribution == prior_spec.nu.EXPONENTIAL, "Call to update_df_svt: Non-matching model specification. Prior for nu should be exponential.");

  const arma::vec& data = log_data2_minus_h;
  const double lambda = prior_spec.nu.exponential.rate;
  const int T = data.size();

  // **** STEP 1: Update tau ****

  // Watch out, R::rgamma(shape, scale), not Rf_rgamma(shape, rate)
  std::transform(
      data.cbegin(), data.cend(),
      tau.begin(),
      [nu](const double data_i) -> double { return 1./R::rgamma((nu + 1.) / 2., 2. / (nu + exp(data_i))); });
  // TODO since C++17 this can be done with transform_reduce (in parallel)
  const double sumtau = std::accumulate(
      tau.cbegin(), tau.cend(),
      0.,
      [](const double partial_sum, const double tau_i) -> double { return partial_sum + std::log(tau_i) + 1. / tau_i;});

  // **** STEP 2: Update nu ****

  const double numean = newton_raphson(nu, sumtau, T, lambda);
  const double auxsd = sqrt(-1/ddlogdnu(numean, T)); 
  const double nuprop = R::rnorm(numean, auxsd);
  const double logR =
    logdnu(nuprop, sumtau, lambda, T) - logdnorm(nuprop, numean, auxsd) -
    (logdnu(nu, sumtau, lambda, T) - logdnorm(nu, numean, auxsd));

  if (logR >= 1. || std::log(R::unif_rand()) < logR) {
    nu = nuprop;
  }
}

void update_general_sv(
    const arma::vec& data,
    const arma::vec& log_data2,
    const arma::ivec& sign_data,
    double& mu,
    double& phi,
    double& sigma2,
    double& rho,
    double& h0,
    arma::vec& h,
    AdaptationCollection& adaptation_collection,
    const PriorSpec& prior_spec,  // prior_mu, prior_phi, prior_sigma2, prior_rho, gammaprior, dontupdatemu feed into this (plus priorlatent0, truncnormal nyi)
    const ExpertSpec_GeneralSV& expert) {  // strategy, correct, use_mala feed into this

  const arma::vec& y = data,
                 & y_star = log_data2;
  const arma::ivec& d = sign_data;
  arma::vec ht = centered_to_noncentered(mu, std::sqrt(sigma2), h);

  // PriorSpec
  const arma::vec prior_mu {prior_spec.mu.normal.mean, prior_spec.mu.normal.sd},  // invalid if prior_spec.mu.distribution != NORMAL
                  prior_phi {prior_spec.phi.beta.alpha, prior_spec.phi.beta.beta},  // invalid if prior_spec.phi.distribution != BETA
                  prior_sigma2 = [&prior_spec]() -> arma::vec2 {
                    switch (prior_spec.sigma2.distribution) {
                      case PriorSpec::Sigma2::GAMMA:
                        return {prior_spec.sigma2.gamma.shape, prior_spec.sigma2.gamma.rate};
                      case PriorSpec::Sigma2::INVERSE_GAMMA:
                        return {prior_spec.sigma2.inverse_gamma.shape, prior_spec.sigma2.inverse_gamma.scale};
                      default:
                        return {};
                    }
                  }(),
                  prior_rho {prior_spec.rho.beta.alpha, prior_spec.rho.beta.beta};  // invalid if prior_spec.rho.distribution != BETA
  const bool dontupdatemu = prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT,
             gammaprior = prior_spec.sigma2.distribution == PriorSpec::Sigma2::GAMMA;
  
  // ExpertSpec
  const bool correct = expert.correct_latent_draws;
  const std::vector<Parameterization>& strategy = expert.strategy;
  const bool use_mala = expert.proposal_para == ExpertSpec_GeneralSV::ProposalPara::METROPOLIS_ADJUSTED_LANGEVIN_ALGORITHM;  // TODO remove, only 'proposal' is needed
  const bool adapt = expert.adapt;

  // only centered
  {
    const LatentVector h_full = draw_latent(y, y_star, d, h, ht, phi, rho, sigma2, mu, correct);
    h0 = h_full.h0;
    h = h_full.h;  // std::move(h_full.h); ?
  }
  ht = centered_to_noncentered(mu, std::sqrt(sigma2), h);
  arma::vec exp_h_half = arma::exp(.5 * h);  // cache exp() calculations
  arma::vec exp_h_half_proposal_nc;

  const Proposal proposal = use_mala ? Proposal::MALA : Proposal::RWMH;
  for (const Parameterization par : strategy) {
    Adaptation& adaptation = par == Parameterization::CENTERED ? adaptation_collection.centered : adaptation_collection.noncentered;
    const ProposalDiffusionKen& adapted_proposal = adapt ? adaptation.get_proposal() : expert.proposal_diffusion_ken;
    bool theta_updated = false;
    if (dontupdatemu) {
      theta_updated = draw_thetamu_rwMH(
          phi, rho, sigma2, mu,
          y, h0, h, ht,
          exp_h_half, exp_h_half_proposal_nc,
          prior_phi,
          prior_rho,
          prior_sigma2,
          par,
          adapted_proposal,
          gammaprior);
      if (adapt) {
        adaptation.register_sample(theta_transform_inv(phi, rho, sigma2, mu).head(3));
      }
    } else {
      theta_updated = draw_theta(
          phi, rho, sigma2, mu,
          y, h0, h, ht,
          exp_h_half, exp_h_half_proposal_nc,
          prior_phi,
          prior_rho,
          prior_sigma2,
          prior_mu,
          par,
          adapted_proposal,
          gammaprior,
          proposal);
      if (adapt) {
        adaptation.register_sample(theta_transform_inv(phi, rho, sigma2, mu));
      }
    }

    if (theta_updated) {
      switch (par) {
        case Parameterization::CENTERED:
          ht = centered_to_noncentered(mu, std::sqrt(sigma2), h);
          break;
        case Parameterization::NONCENTERED:
          h = noncentered_to_centered(mu, std::sqrt(sigma2), ht);
          exp_h_half = std::move(exp_h_half_proposal_nc);
          break;
      }
    }
  }

  return;
}

}

