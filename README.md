# stochvol
Partial re-write of the R package stochvol to allow for asymmetry (leverage).

## To discuss

### svpredict
`svlpredict` makes more sense with predicted `y` values included. For coherence, both `svpredict` and `svlpredict` could include predicted `y` values. Additionally, there's no really need for `arpredict` then, the usual `predict` could do it, and a designmatrix could be added for prediction as well.

### update function
It could use an `Rcpp` interface. It copies now, it wouldn't copy then. However, then `factorstochvol` would need an update as well.

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
* arpredict?
* -logret
* -paradensplot
* -paratraceplot
* -svlsample
* -svlsample2
* -svlsim
* -updatesummary
* -volplot
* -mytraceplot
* -mydensplot
* -svlsample_cpp

### Features
* -regression, betas
* -`update` function for svlsample

### Code cleanup
* -remove unnecessary `svlsamplr` code

### Optimize
* `sampler` should do smarter thinning
* `svlsamplr` should not copy that much

### Document
* roxygen

### Other
* -Depends: methods (`Rscript --vanilla` doesn't load it)
* -prior_mu in `draw_h_auxiliary`
* `ggplot2`
* -maybe c('svldraws', 'svdraws') should be the new class structure
* create joint file structure
* check TODOs in the code
* -plot correct prior if gammaprior=FALSE
