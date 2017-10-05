\name{predict.svdraws}
\alias{predict.svdraws}
\title{
Prediction of Future Log-Volatilities
}
\description{
Simulates draws from the predictive density of the latent log-volatility process.
}
\usage{
\method{predict}{svdraws}(object, steps = 1, ...)
}
\arguments{
  \item{object}{
   \code{svdraws} object.
}
  \item{steps}{
   single number, coercible to integer. Denotes the number of steps to forecast.
}
  \item{...}{
   currently ignored.
}
}
\value{
 Returns an object of class \code{c("svpredict", "mcmc")} containing simulations from the predictive density of \code{h_(n+1),...,h_(n+steps)}.
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}
\note{
You can use the usual \code{coda} methods for \code{mcmc} objects to print, plot, or summarize the predictions, or use them within \code{\link{volplot}} or \code{\link{plot.svdraws}}.
}

\seealso{
\code{\link{plot.svdraws}}, \code{\link{volplot}}.
}

\examples{
## Simulate a short and highly persistent SV process 
sim <- svsim(100, mu = -10, phi = 0.99, sigma = 0.2)

## Obtain 5000 draws from the sampler (that's not a lot)
draws <- svsample(sim$y, draws = 5000, burnin = 100,
		  priormu = c(-10, 1), priorphi = c(20, 1.5), priorsigma = 0.2)

## Predict 10 days ahead
fore <- predict(draws, 10)

## Check out the results
summary(fore)
plot(draws, forecast = fore)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
