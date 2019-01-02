# Make sure that the unwanted version of stochvol is not loaded

generate <- FALSE  # generate test file, or do comparison?
if (generate) {
  library(stochvol)
} else {
  mydevlib <- "~/R/under_development"
  devtools::install(".", args = paste0("--library=", mydevlib), repos = NULL, dependencies = FALSE)
  library(RcppArmadillo)
  library("methods")
  library(stochvol, lib.loc = mydevlib)
}

set.seed(1)
foo <- svsample(rnorm(50), 10, 100)
filename <- "test-with-old-stochvol.RDS"

if (generate) {
  saveRDS(foo, filename)
} else {
  foo.old <- readRDS(filename)
  foo$runtime <- foo.old$runtime <- NULL
  if (all.equal(foo$para, foo.old$para) &&
      all.equal(foo$latent, foo.old$latent)) {
    message("Test passed!")
  } else {
    stop("Test failed!")
  }
}
