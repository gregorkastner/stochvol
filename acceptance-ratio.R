library("tidyverse")

devtools::load_all()

if (!file.exists("apple.RDS")) {
  dat_q <- quantmod::getSymbols("AAPL", env = NULL, from = "2015-01-01")
  dat <- data.frame(Date = zoo::index(dat), AAPL = zoo::coredata(dat_q$AAPL.Adjusted))
  saveRDS(dat, "apple.RDS")
} else if (!exists("dat")) {
  dat <- readRDS("apple.RDS")
}

y <- dat$AAPL
y <- y[is.finite(y)]
y <- diff(log(y))

results <- tibble(T = integer(),
                  Runtime = numeric(),
                  Mean.Means = numeric(),
                  Mean.Stdevs = numeric(),
                  Stdev.Means = numeric(),
                  Stdev.Stdevs = numeric(),
                  Rho = numeric(),
                  Phi = numeric(),
                  Sigma = numeric(),
                  Mu = numeric())
i <- 0
for (n in c(10, 100, 1000)) {
  cat(n, "\n")
  for (phi in c(0.95, 0.99)) {
    cat(phi, "\n")
    for (rho in c(0, -0.2, -0.4)) {
      cat(rho, "\n")
      for (sigma in c(0.1, 0.3)) {
        cat(sigma, "\n")
        for (mu in c(-7, -10)) {
          cat(mu, "\n")
          i <- i+1
          cat("i/n = ", i, "/72\n\n")
          yy <- head(y, n)
          yy <- yy - mean(yy)
          startpara <- list(phi = phi, sigma = sigma, mu = mu, rho = rho)
          res_sv <- svsample(yy, quiet = TRUE)
          startlatent <- apply(res_sv$latent, 2, median)
          runt <- system.time({
            res <- svlsample2(yy, draws = 300, burnin = 50,
                              priormu = c(0, 100), priorphi = c(5, 1.5), priorsigma = 1, priorrho = c(4, 4),
                              thinpara = 1, thinlatent = 1, thintime = NULL, keeptime = "all",
                              quiet = TRUE, startpara, startlatent)
          })
          res_df <- as.data.frame(res)
          res_df$sd <- sqrt(res_df$var)
          colMeans(res_df)
          apply(res_df, 2, sd)
          results <- results %>%
            add_row(T = n,
                    Runtime = summary(runt)["user"],
                    Mean.Means = mean(res_df$mean),
                    Mean.Stdevs = mean(res_df$sd),
                    Stdev.Means = sd(res_df$mean),
                    Stdev.Stdevs = sd(res_df$sd),
                    Rho = rho,
                    Phi = phi,
                    Sigma = sigma,
                    Mu = mu)
        }
      }
    }
  }
}

saveRDS(results, "results.RDS")

