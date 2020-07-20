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

* [Dealing with Stochastic Volatility in Time Series Using the R Package stochvol](https://cran.r-project.org/package=stochvol/vignettes/article.pdf)
* [Heavy-Tailed Innovations in the R Package stochvol](https://cran.r-project.org/package=stochvol/vignettes/heavytails.pdf).

For more information, please visit the [CRAN page](https://cran.r-project.org/package=stochvol) of the package.
