# stochvol <img src="man/figures/logo.svg" align="right" padding-left="3px" />

This is the development repository of the [`R`](https://www.r-project.org/) package [`stochvol`](https://cran.r-project.org/package=stochvol).

## Features

The package provides methods to estimate the stochastic volatility model, potentially with conditionally heavy tails and/or with leverage.
Using functions `svsample`, `svtsample`, `svlsample`, and `svtlsample`, one can conduct Bayesian inference on all parameters, including the time-varying volatilities (the states in the state space).

Additional features:

* Prediction, plotting, residual extraction work with the usual functions in `R` (`predict`, `plot`, and `residuals`)
* Choose from a range of prior distrubutions; see `help("specify_priors", package="stochvol")`
* Built-in support for linear regression and autoregressive processes with stochastic volatility errors; look for function argument `designmatrix`
* Easy interfacing with [`bayesplot`](https://cran.r-project.org/web/packages/bayesplot/) functions through the `as.array()` specialization
* Rolling or expanding window estimation can be used for backtesting; see `help("svsample_roll", package="stochvol")`
* Run independent Markov chains using `R`'s cross-platform parallelization; look for function arguments `n_chains`, `parallel`, `n_cpus`, and `cl` (for "cluster")
* For plug&play Bayesian modeling, when stochastic volatility is part of a larger model, fast-access functions can speed up execution in `R`; see `help("svsample_fast_cpp", package="stochvol")`
* For advanced users, there is a `C++` interface; see e.g. `help("update_fast_sv", package="stochvol")`
* For teaching purposes, you can fix any parameter to a known value using `sv_constant` as the prior specification

## Install CRAN Version

Type into your `R` session:

```r
install.packages("stochvol")
```

For more information, please visit the [CRAN page](https://cran.r-project.org/package=stochvol) of the package.

## Install Latest Development Version

Type into your `R` session:

```r
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github(
  repo = "https://github.com/gregorkastner/stochvol")
```

## Documentation

The best introduction is the combination of our vignettes:

* [Dealing with Stochastic Volatility in Time Series Using the R Package stochvol](https://www.jstatsoft.org/index.php/jss/article/view/v069i05)
* [Modeling Univariate and Multivariate Stochastic Volatility in R with stochvol and factorstochvol](https://arxiv.org/abs/1906.12123)

For individual functions, please refer to the help pages after installing the package.
For instance, for `svsample`, execute

```r
help("svsample", package = "stochvol")
```

For more information, please visit stochvol's [CRAN page](https://cran.r-project.org/package=stochvol).
