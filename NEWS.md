# stochvol 3.0.10

- NEW FUNCTION 'svlm', which has a formula interface;
  it is a wrapper around 'svsample'; many thanks
  to Peter Knaus for his help
- Turn on printing on Windows
- Implement Geweke test in C++; it is feasible to
  execute it as a CRAN test now
- Small change in the behavior of 'predict.svdraws':
  when 'newdata' is given then 'steps' is ignored;
  A warning is shown if this is relevant
- Updated examples for 'svsample' and 'predict' so
  that they use the extractors as intended
- Simplified vignette: cache some results to
  reduce dependencies

# stochvol 3.0.6

- Correct wrongly submitted vignettes

# stochvol 3.0.5

- Re-added first vignette with the necessary updates
- Bugfix in 'residuals.svdraws' and 'plot.svresid'; thanks
  to David Zoltan Szabo

# stochvol 3.0.4

- Updated vignette
- Prevent some compilation warnings
- Updated CITATION file

# stochvol 3.0.3

- Bugfix in 'svsample' when the input data contains 0s

# stochvol 3.0.2

- Bugfix in the fast access functions 'svsample_fast_cpp' and
  'svsample_general_cpp': dimensions were incorrect
- Bring out the adaptation object for 'svsample_general_cpp' to
  the R level; used to be accessible from C++ only
- Update to the adaptation object in C++

# stochvol 3.0.1

- Bugfix for calls to 'svsample' with draws = 1
- Fix some #include directives
- Avoid compilation problems on Solaris

# stochvol 3.0.0

- New model: heavy-tailed SV with leverage and its sampling
  function 'svtlsample'
- Change in the heavy-tailed models: now they operate with a
  normalized Student's t-distribution, i.e.\ its standard
  deviation is 1 for all degrees of freedom
- New, optional API for specifying prior distributions
- It is possible to set any or all of 'mu', 'phi', 'sigma',
  'nu', 'rho' to a constant value, i.e.\ set a Dirac prior
- Replaced the prior distribution for the degrees of freedom
  parameter 'nu' with an exponential distribution; it used to be
  a uniform distribution
- New set of unit tests including a Geweke test
- Entirely refactored C++ backend
- New, rethought set of exported C++ functions: 'update_fast_sv',
  'update_general_sv', 'update_regression', and 'update_t_error';
  all of them are documented
- Vastly improved the computational efficiency of the former
  'svlsample' code
- New feature: run independent MCMC chains on the same data set
- Integration of the 'parallel' package: option to run independent
  MCMC chains on a 'SNOW' cluster or using the 'multicore' strategy
- Modified backend for 'svdraws' and 'svpredict' objects to
  to incorporate independent chains: they contain 'mcmc.list'
  objects instead of plain 'mcmc' objects.
  WARNING! This may break code that exploited the backend!
- Additional extractor functions for 'svdraws' objects: 'vola',
  'sv_beta', 'sv_tau', 'sampled_parameters', index chains via [],
  'as.array'
- New extractors for 'svpredict' objects: 'predy', 'predlatent',
  'predvola', index chains via []
- CITATION file updated
- New vignette
- SV with leverage also includes 'latent0' from now on
- New feature: rolling window estimation via the functions
  'svsample_roll', 'svtsample_roll', 'svlsample_roll', and
  'svtlsample_roll'
- Unified plotting between models: removed the 'scaling' plot
  from heavy-tailed model outputs
- Many new, small examples
- Small bugfixes
- New "fast-access" functions to circumvent input validation:
  'svsample_fast_cpp' and 'svsample_general_cpp'; their arguments
  are slightly different from the 'svsample' family of functions,
  they come with documentation
- 'svsample2' is set to deprecated

# stochvol 2.0.5

- Bugfix in predict when a constant mean model is used for prediction
- Test suite added
- Updated CITATION file

# stochvol 2.0.4

- Bugfix in 'svlpredict' when an inverse gamma prior is used.
- Function 'svtpredict' introduced for convenience when estimating
  SV with heavy tails.
- New generic: logret; old logret function is now logret.default
  (this change shouldn't bother the end user).
- New generic: paratraceplot; old paratraceplot function is now
  paratraceplot.svdraws (this change shouldn't bother the end user).

# stochvol 2.0.3

- The exported C++ function 'update_svl' can draw posterior values
  with a fixed 'mu'.
- Controlling the amount of latent variable draws to be stored can now be
  controlled through the 'keeptime' argument to 'svsample', 'svsample2',
  'svlsample', and 'svlsample2'.

# stochvol 2.0.2

- Bugfix in 'predict.svdraws'

# stochvol 2.0.1

- Introduced a minimalistic 'plot.svpredict' / 'plot.svlpredict'
- More fine-grained control over the covariance matrix in 'svlsample',
  the default employs an approximate covariance matrix coming from 'svsample'
- Bugfix in 'predict.svdraws' and minor changes in documentation thereof

# stochvol 2.0.0

- New collaborator: Darjus Hosszejni
- New functionality:
  - leverage effect through 'svlsample'
  - simulation of asymmetric returns with 'svsim'
  - prediction using designmatrices and heavy-tails
- Updated functionality:
  - the exported C++ function 'update' uses RcppArmadillo objects and
    it has been renamed to 'update_sv'
  - new examples
- Bugfix:
  - correct plotting after sampling with 'gammaprior = FALSE'
  - fixed naming of latent states in printing and plotting in the case of time
    series with conditionally heavy-tailed innovations
  - in 'svsample': input check for length of burnin is now fixed (thanks
    to Nikolas Kuschnig)
  - other small fixes
- Deprecated:
  - functions 'arpredict' and '.svsample'
  - parameter 'thintime' in function 'svsample'

# stochvol 1.3.4

- DESCRIPTION file updated.

# stochvol 1.3.3

- Registered native routines.
- Thinned graphs in vignette a bit to stay below the required 5MB mark.
- Workaround for typo appearing in conversion of help files to LaTeX document.
  Thanks to Evelyn Mitchell for pointing this out.

# stochvol 1.3.2

- Some more changes in the AWOL sampler because version 1.3.1 appeared to be
  unstable under Solaris: Sampling of h (including h0) is now done AWOL.
- Turned progress indicator off again (not sure whether this causes issues
  with Solaris).

# stochvol 1.3.1

- Parameter starting values are now set to their respective prior means
  (if not specified by the user).
- Sampling h[-0] is now done conditionally on h0. This change should not
  make a difference for the standard use cases of stochvol but is necessary
  to yield correct results for the new prior for h0 which was introduced in
  version 1.3.0.
- Turned on progress indicator also for Windows (need %% instead of \045
  or \% to escape a percent symbol in Rprintf).

# stochvol 1.3.0

- Implemented new prior for initial log-volatility h0. Can be used to
  stabilize the level of the volatilities, especially when stochvol
  is used within the context of a more general MCMC algorithm.
- Shortened the heavytails-vignette (a bit) to stay below the 5MB mark.
- Fixed error in help files of svsample and svsample2 concerning startpara.
  Thanks to Dominic Cervicek for pointing this out.

# stochvol 1.2.4

- Fixed typo which appeared in pdf package documentation. Thanks to Stefan
  Voigt for pointing this out.

# stochvol 1.2.3

- Updated CITATION to cater for newly published article in the Journal of
  Statistical Software 69(5), 1-30, doi: 10.18637/jss.v069.i05.
- The function logret can now handle xts objects. Thanks to Niko Hauzenberger
  for pointing this out.

# stochvol 1.2.2

- Fixed bug in initialization of svsample2 that was introduced in 1.2.1.
  Thanks to John Kerpel for pointing this out.

# stochvol 1.2.1

- Added option 'keeptau' to the store the "variance inflation factors"
  used for the sampler with conditional t innovations. Thanks to Sergey Egiev
  for pointing out that this may be useful to at what point(s) in time the
  normal disturbance had to be "upscaled" by a mixture factor and when the
  series behaved "normally".
- Fixed residuals.svdraws to cater for a potential 'designmatrix'.

# stochvol 1.2.0

- Allows for incorporation of a simple mean (regression-type) model.
  For details, please see ?svsample and ?arpredict. And, hopefully
  soon, the corresponding vignette.
- Improved input-checking (somewhat).
- Argument 'thintime' in svsample may now also be 'firstlast',
  meaning that latent volatility draws are only kept for the first
  and the last point in time.

# stochvol 1.1.4

- Included non-base default packages that are imported in NAMESPACE
  into DESCRIPTION. Thanks to Kurt Hornik for pointing this out.

# stochvol 1.1.3

- Fixed a typo in Makefile for building vignettes.

# stochvol 1.1.2

- Fixed missing imports from non-base default packages.

# stochvol 1.1.1

- Fixed "uninitialized variable" in function newtonRaphson
  (progutils.cpp). Warning appeared when building for Windows.
- Minor changes in main vignette (article.Rnw) so that the included
  plots have smaller file size and the package does not exceed the 5MB
  limit (was the case on r-patched-solaris-sparc).

# stochvol 1.1.0

- Introduced sampling and simulating conditional t-innvoations.
  (Hopefully) all affected functions
  - svsample
  - svsample2
  - volplot
  - paratraceplot
  - paradensplot
  - plot.svdraws
  - plot.svresid
  - svsim
  - updatesummary
  - predict.svdraws
  - residuals.svdraws
  - print.svsim
  - print.summary.svsim
  - print.summary.svdraws
  and the corresponding help files have been adapted.
  See package vignette "heavytails" for details.
  Still a "beta-feature", please use with care.
  Bug reports and/or comments are warmly welcome!
- The plot functions now take an optional svsim object.
- predict.svdraws is now in R/utilities_svdraws.R (was in R/plotting.R)

# stochvol 1.0.2

- Fixed some more typos.
- Plotting functions are now in R/plotting.R

# stochvol 1.0.1

- Updated keywords in vignette.
- Changed VignetteKeyword entries to adhere to CRAN style guide.

# stochvol 1.0.0

- First stable release.
- Minor stylistic changes in package vignette.
- Updated CITATION file.
- Version numbering changed from X.Y-Z to X.Y.Z.

# stochvol 0.9-2

- Fixed some typographical errors in help files. Thanks to Angela Bitto
  for pointing these out.

# stochvol 0.9-1

- Fixed typographical errors in vignette.
- Added additional results in Table 1 of vignette.
- Minor changes in DESCRIPTION to adhere to CRAN style guide.

# stochvol 0.9-0

- Requires now Rcpp >= 0.11.0 and consequently R >= 3.0.0.
- Replaced RNGScope in sampler.cpp by "manual" GetRNGstate() and
  PutRNGstate() statements because the former is not safe if return
  variables are declared afterwards, cf.
  https://www.mail-archive.com/rcpp-devel@lists.r-forge.r-project.org/msg07519.html
  and follow-ups for further information.
- Minimal overhead sampling function .svsample was renamed to
  svsample2. The former name will be faded out; please use svsample2
  from now on.
- Substantial rewrite of vignette, especially Section 5; includes now
  a comparison of SV and GARCH.
- Updates and corrections in .Rd files.
- Added some info when package is attached.
- Sample size is now denoted as "n" instead of "T", in order to avoid
  confusion with "TRUE" and/or the transpose symbol.

# stochvol 0.8-4

- Fixed a bug introduced in version 0.8-2, where I forgot to let
  .svsample know about the changes in update(...). Thanks to Keiran
  Thompson for pointing this out.
- Updated CITATION file.

# stochvol 0.8-3

- Minor changes in vignette.

# stochvol 0.8-2

- update(...) now takes an additional boolean argument preventing
  an update of mu (but instead hold mu = 0 fixed). This feature is
  needed e.g. for factor stochastic volatility models.
- C++ function "update" rewritten to use ordinary pointers instead of
  Rcpp objects. This became necessary because Rcpp objects cannot be
  constructed re-using memory due to the SEXPREC structure. See
  https://www.mail-archive.com/rcpp-devel@lists.r-forge.r-project.org/msg06811.html
  and follow-ups for a more detailed explanation.

# stochvol 0.8-1

- Added a -I directive for compiler in Makevars.win so that winbuilder
  finds headers in inst/include/.
- Re-ran all code used in vignette.

# stochvol 0.8-0

- Re-structured main sampler code in sampler.cpp: new function "update"
  performs a single MCMC iteration. Additionally, this function is also
  made available to be called by C/C++ code in another package. To use
  this function within C/C++ code, simply #include <update.h> (defined
  in inst/include) in the C/C++ code of your package, which itself
  Imports:/Depends: stochvol (>= 0.8). See package factorstochvol for
  an example using this mechanism.
- Minor stylistic changes in helper functions and subroutines at
  C++-level).
- Fixed a minor bug in sampler.cpp where regressionNoncentered
  returned sigma instead of fabs(sigma) when phi is drawn outside of
  the unit sphere.
- Fixed DESCRIPTION to Import: Rcpp (instead of Depend:).
- Fixed NAMESPACE: Now has "importFrom(Rcpp, evalCpp)" as required by
  Rcpp 0.11.0.

# stochvol 0.7-2

- Updated CITATION file.
- Updated vignette (mainly stylistic changes).

# stochvol 0.7-1

- Manual is now provided as a vignette.

# stochvol 0.7-0

- Added a manual.
- Changed the default prior for mu from c(-10, 3) to c(0, 100) to
  avoid disaster when percentage log-returns (instead of log-returns)
  are used but the prior is not adapted accordingly.
- Included the date in the exrates dataset.
- Changed svsample to also accept vectors with zero-returns; it throws
  a warning now (instead of an error) and adds an offset constant.
- Fixed plot.svdraws to reset par to old values.
- Added a "residual plot" feature. Thanks to Hedibert Lopes for
  pointing this out.
- Added the convenience function "logret" for taking log-returns of a
  given series, with the possibility of de-meaning.
- Cleaned up summary.svdraws, which now returns a "summary.svdraws"
  object. For actual printing, print.summary.svdraws is used.
- Cleaned up summary.svsim, which now returns a "summary.svsim" object.
  For actual printing, print.summary.svsim is used.

# stochvol 0.6-1

- Fixed typo in wrapper.R which previously disallowed changing the
  "expert" argument "proposalvar4sigmaphi": Replaced
  "proposalvar2sigmaphi" by "proposalvar4sigmaphi" on line 76.
- Included EUR exchange rates from the European Bank's Statistical Data
  Warehouse. Use with "data(exrates)".
- Added CITATION file.
- Replaced Rprintf by REprintf and cat(...) by cat(..., file=stderr())
  for status reports. Thanks to Kurt Hornik for pointing this out.

# stochvol 0.6-0

- Introduced ".svsample" for minimal overhead sampling. Intended to
  be used as a plug-in into other MCMC samplers. No input checking, use
  with proper care!
- Disabled progress bar in non-Unix-like systems due to problems with
  console flushing.
- Some bug fixes for solaris. Seems to build fine now. Thanks to Brian
  Ripley for reporting the bugs.

# stochvol 0.5-1

- Replaced all paste0(...) calls by paste(..., sep='') for
  compatibility reasons.

# stochvol 0.5-0

- First CRAN release version.

# Initial TODO

- Code updatesummary in C (apply w/ quantiles is slow).
- Make AWOL sampler available for fixed parameters. Thanks to Hedibert
  Lopes for pointing this out.
