myplot <- function (X) {
  par(mfrow = c(4, 1))
  for (i in 1:4) {
    plot(X[, i], type = "l")
  }
}

myplot2 <- function (X) {
  #par(mfrow = c(4, 1))
  for (i in 1:4) {
    traceplot(X[, i])
  }
}

acceptance.rate <- function (X) {
  mean(tail(X[, 1], -1) != head(X[, 1], -1))
}

rangetransform <- function (x, r, funcinv, func) {
  tx <- funcinv(x)
  txl <- tx-r
  txu <- tx+r
  func(txu)-func(txl)
}

funcinv.phi <- function (x) {
  .5*log(2/(1-x)-1)
}
func.phi <- function (tx) {
  1-2/(exp(2*tx)+1)
}

funcinv.sigma2 <- log
func.sigma2 <- exp

rangetransform.phi <- function (x, r) {
  rangetransform(x, r, funcinv.phi, func.phi)
}
rangetransform.sigma2 <- function (x, r) {
  rangetransform(x, r, funcinv.sigma2, func.sigma2)
}

asisprint <- function (x, censtring) {
  if (length(x) %% 2 == 0) {
    toshorten <- rep(as.character(censtring), length(x)/2)
    if (identical(toshorten, as.character(x)))
      return(sprintf("ASISx%d", length(x)/2))
  }
  paste0("(", paste(x, collapse=", "), ")")
}
