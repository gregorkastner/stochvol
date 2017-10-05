\name{arpredict}
\alias{arpredict}
\title{
Dynamic prediction for the AR-SV model
}
\description{
Simulates draws from the posterior predictive density of a fitted AR-SV model.
}
\usage{
arpredict(object, volpred)
}
\arguments{
  \item{object}{
   \code{svdraws} object as returned from \code{\link{svsample}}.
}
  \item{volpred}{
   \code{svpredict} object as returned from \code{\link{predict.svdraws}}.
}
}
\value{
 Returns an object of class \code{c("distpredict", "mcmc")} containing simulations from the posterior predictive density of \code{y_(n+1),...,y_(n+steps)}.
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}
\note{
You can use the usual \code{coda} methods for \code{mcmc} objects to print, plot, or summarize the predictions.
}

\seealso{
\code{\link{predict.svdraws}}.
}

\examples{\dontrun{
data(exrates)
y <- exrates$USD

## Fit AR(1)-SV model to EUR-USD exchange rates
res <- svsample(y, designmatrix = "ar1")

## Use predict.svdraws to obtain predictive volatilities
ahead <- 100
predvol <- predict(res, steps = ahead)

## Use arpredict to obtain draws from the posterior predictive
preddraws <- arpredict(res, predvol)

## Calculate predictive quantiles
predquants <- apply(preddraws, 2, quantile, c(.1, .5, .9))

## Visualize
ts.plot(y, xlim = c(length(y) - ahead, length(y) + ahead),
	ylim = range(predquants))
for (i in 1:3) {
 lines((length(y) + 1):(length(y) + ahead), predquants[i,],
       col = 3, lty = c(2, 1, 2)[i])
}
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
