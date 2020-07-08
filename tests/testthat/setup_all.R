# Create input for the samplers
y <- c(0.01, -0.01, 0.1, -0.09, -0.01, 0.001, -0.001, 0.01, -0.009, -0.001)

draws <- 300
burnin <- 10

designmat <- matrix(c(0.3485322, 0.536128, 0.2978863, -2.045419, 1.390987, -1.76708, 0.2594924, -1.576583, -1.303889, -0.6729478, -0.6123462, 0.6823703, 1.824069, 0.7192591, -1.059588, -1.303618, 0.3878489, -0.4601865, 1.912537, -0.6549144), 10, 2)
ar_values <- c(0, 2)

designmatrix_values <- c(list(NA), as.list(paste0("ar", ar_values)), list(designmat))
keeptime_values <- c("all", "last")
thin_values <- c(1, 3)

# Run samplers (repeated code not nice but readable)
res_sv <- list(svsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svt <- list(svtsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svl <- list(svlsample(y, draws = draws, burnin = burnin, quiet = TRUE))  # default values
res_svl_corrected <- list(svlsample(y, draws = draws, burnin = burnin, expert = list(correct.latent.draws = TRUE), quiet = TRUE))  # default values
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
    }
  }
}

