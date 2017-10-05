\name{exrates}
\alias{exrates}
\docType{data}
\title{
  Euro exchange rate data
}
\description{
The data set contains the daily bilateral prices of one Euro in 23 currencies from January 3, 2000, until April 4, 2012. Conversions to New Turkish Lira and Fourth Romanian Leu have been incorporated.
}

\usage{
data(exrates)
}

\source{
ECB Statistical Data Warehouse (\url{http://sdw.ecb.europa.eu})
}

\seealso{
 \code{\link{svsample}}
}

\examples{
\dontrun{
data(exrates)
dat <- logret(exrates$USD, demean = TRUE)  ## de-meaned log-returns
res <- svsample(dat)                       ## run MCMC sampler
plot(res, forecast = 100)                  ## display results
}
}
\keyword{datasets}
