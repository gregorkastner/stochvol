# Stochastic volatility with leverage
An implementation of the stochastic volatility model with leverage [1, 2].

## Install
To install the package from GitHub, you should first install `devtools`, then the rest:
```r
install.packages("devtools")
library("devtools")
install_github("hdarjus/stoch-vol-with-leverage@0.1")
```

# Development details
The variable names use the Hungarian notation (https://en.wikipedia.org/wiki/Hungarian_notation).

# References
[1] Nakajima, Jouchi, and Yasuhiro Omori. "Leverage, heavy-tails and correlated jumps in stochastic volatility models." Computational Statistics & Data Analysis 53.6 (2009): 2335-2353.

[2] Omori, Yasuhiro, et al. "Stochastic volatility with leverage: Fast and efficient likelihood inference." Journal of Econometrics 140.2 (2007): 425-449.
