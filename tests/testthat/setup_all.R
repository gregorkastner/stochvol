# Create input for the samplers
y <- c(-0.867537088, 0.424335171, -0.048784030, -0.749483682,
       -0.332278654, 0.503644166, 0.007733118, -0.485776150,
       0.163207939, 0.123654636, 0.535032546, 0.799077152,
       0.306340493, -1.451189223, -1.095406999, 2.765907287,
       -2.196736048, 0.879586318, 0.359074117, -1.700210392)

draws <- 500
burnin <- 100

designmat <- matrix(c(0.18554962, 1.71128430, 2.76005663, 0.27480471, -0.04419060,
                      -1.54191962, 0.03668254, 1.71851215, -1.48924188, 1.45517746,
                      -0.03498587, 0.91532504, 0.61450314, 0.64796734, 0.41568157,
                      -1.79375801, -0.33170476, -0.06896315, -0.83336860, 1.35250078,
                      -2.17504467, 0.19803696, -0.14800168, -0.95317998, 0.92041904,
                      -0.96628316, 0.59918070, 0.81186793, 0.07408481, 0.25670553,
                      -0.56785729, 0.77431129, 1.04092342, -0.52666850, 1.39242409,
                      1.79126519, 0.21054807, 0.19770369, -0.90162245, -0.62361923),
                    length(y), 2)
ar_values <- c(0, 2)

designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
keeptime_values <- c("all", "last")
thin_values <- c(1, 3)

# Run samplers (repeated code not nice but readable)
res_sv <- list(svsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svt <- list(svtsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svl <- list(svlsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svl_corrected <- list(svlsample(y, draws = draws, burnin = burnin, expert = list(correct.latent.draws = TRUE), quiet = TRUE))  # default values
res_svl_mala <- list(svlsample(y, draws = draws, burnin = burnin, expert = list(mhcontrol = list(use.mala = TRUE), correct.latent.draws = FALSE), quiet = TRUE))  # default values
for (dm in designmatrix_values) {
  for (kt in keeptime_values) {
    for (th in thin_values) {
      res_sv <- c(res_sv, list(  # non-default values
        svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)
      ))
      res_svt <- c(res_svt, list(  # non-default values
        svtsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)
      ))
      res_svl <- c(res_svl, list(  # non-default values
        svlsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)
      ))
      res_svl_corrected <- c(res_svl_corrected, list(  # non-default values
        svlsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, expert = list(correct.latent.draws = TRUE), quiet = TRUE)
      ))
      res_svl_mala <- c(res_svl_mala, list(  # non-default values
        svlsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, expert = list(mhcontrol = list(use.mala = TRUE), correct.latent.draws = TRUE), quiet = TRUE)
      ))
    }
  }
}

draws <- 30
burnin <- 10

