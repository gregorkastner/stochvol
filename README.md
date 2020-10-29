# stochvol
This is the development repository of the [`R`](https://www.r-project.org/) package [`stochvol`](https://cran.r-project.org/package=stochvol).

# Install CRAN Version
The code chunk below should work.
Type into your `R` session:
```r
install.packages("stochvol")
```
For more information, please visit the [CRAN page](https://cran.r-project.org/package=stochvol) of the package.

# Install Latest Stable Version
Type into your `R` session:
```r
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("https://github.com/gregorkastner/stochvol")
```

# Install Latest Development Version
Type into your `R` session:
```r
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github(
  repo = "https://github.com/gregorkastner/stochvol",
  ref = "dev")
```

# Usage
The best introduction is the combination of our vignettes:

* [Dealing with Stochastic Volatility in Time Series Using the R Package stochvol](https://www.jstatsoft.org/index.php/jss/article/view/v069i05/v69i05.pdf)
* [Modeling Univariate and Multivariate Stochastic Volatility in R with stochvol and factorstochvol](https://cran.r-project.org/web/packages/stochvol/vignettes/article2.pdf)

For more information, please visit stochvol's [CRAN page](https://cran.r-project.org/package=stochvol).
