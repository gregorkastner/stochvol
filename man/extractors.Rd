\name{Extractors}
\alias{para}
\alias{latent}
\alias{latent0}
\alias{priors}
\alias{thinning}
\alias{runtime}
\title{
Common Extractors for 'svdraws' Objects
}
\description{
Some simple extractors returning the corresponding element of an \code{svdraws} object.
}
\usage{
para(x)
latent(x)
latent0(x)
priors(x)
thinning(x)
runtime(x)
}
\arguments{
  \item{x}{\code{svdraws} object.}
}
\value{
 The return value depends on the actual funtion:
 \item{para(x)}{extracts the parameter draws and returns them as an \code{mcmc} object.}
 \item{latent(x)}{extracts the latent contemporaneous log-volatility draws and returns them as an \code{mcmc} object.}
 \item{latent0(x)}{extracts the latent initial log-volatility draws and returns as an \code{mcmc} object.}
 \item{priors(x)}{extracts the prior parameters used and returns them in a \code{list}.}
 \item{thinning(x)}{extracts the thinning parameters used and returns them in a \code{list}.}
 \item{runtime(x)}{extracts the runtime and returns it as a \code{proc_time} object.}
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}

\keyword{utilities}
