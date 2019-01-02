# stochvol
Partial re-write of the R package stochvol to allow for asymmetry (leverage).

## TODOs

### API
* plot,svldraws
* plot,svlresid
* plot,svlsim
* predict,svldraws
* print,summary.svldraws
* print,summary.svlsim
* print,svldraws
* print,svlsim
* residuals,svldraws
* summary,svldraws
* summary,svlsim
* arpredict
* logret
* paradensplot
* paratraceplot
* svlsample
* svlsample2
* svlsim
* updatesummary
* volplot

### Features
* regression
* 'update' function

### Code cleanup
* remove unnecessary svlsamplr code

### Optimize
* 'sampler' should do smarter thinning
* svlsamplr should not copy that much

### Document
* roxygen

### Other
* Depends: methods ("Rscript --vanilla" doesn't load it)
* ggplot2

