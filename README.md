# stochvol
Partial re-write of the R package stochvol to allow for asymmetry (leverage).

## TODOs

### API
* -plot,svldraws
* -plot,svlresid
* -plot,svlsim
* -predict,svldraws
* -print,summary.svldraws
* -print,summary.svlsim
* -print,svldraws
* -print,svlsim
* -residuals,svldraws
* -summary,svldraws
* -summary,svlsim
* arpredict
* -logret
* -paradensplot
* -paratraceplot
* -svlsample
* svlsample2
* -svlsim
* -updatesummary
* -volplot
* -mytraceplot
* -mydensplot
* -svlsample_cpp

### Features
* regression, betas
* 'update' function for svlsample
* 'update' in Rcpp
* factorstochvol and update in Rcpp

### Code cleanup
* -remove unnecessary svlsamplr code

### Optimize
* 'sampler' should do smarter thinning
* svlsamplr should not copy that much

### Document
* roxygen

### Other
* Depends: methods ("Rscript --vanilla" doesn't load it)
* -prior_mu in draw_h_auxiliary
* ggplot2
* maybe c('svldraws', 'svdraws') should be the new class structure
* G: change svpredict's format to include prediction about 'y' as well. or simply discard the y's predicted by svlpredict?
* create joint file structure
* G: inverse gamma prior parameters!
* G: plot inverse gamma prior

