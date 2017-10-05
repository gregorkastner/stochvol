\name{paradensplot}
\alias{paradensplot}
\title{
Probability Density Function Plot for the Parameter Posteriors
}
\description{
Displays a plot of the density estimate for the posterior distribution of the parameters \code{mu}, \code{phi}, \code{sigma} (and potentially \code{nu}), computed by the \code{\link[stats]{density}} function.
}
\usage{
paradensplot(x, showobs = TRUE, showprior = TRUE, showxlab = TRUE,
	     mar = c(1.9, 1.9, 1.9, 0.5), mgp = c(2, 0.6, 0),
	     simobj = NULL, ...)
}
\arguments{
  \item{x}{
   \code{svdraws} object.
}
  \item{showobs}{
   logical value, indicating whether the observations should be displayed along the x-axis. If many draws have been obtained, the default (\code{TRUE}) can render the plotting to be quite slow, and you might want to try setting \code{showobs} to \code{FALSE}.
}
  \item{showprior}{
   logical value, indicating whether the prior distribution should be displayed. The default value is \code{TRUE}.
}
  \item{showxlab}{
   logical value, indicating whether the x-axis should be labelled with the number of iterations and the bandwith obtained from \code{\link[stats]{density}}. The default value is \code{TRUE}.
}
  \item{mar}{
   numerical vector of length 4, indicating the plot margins. See \code{\link[graphics]{par}} for details. The default value is \code{c(1.9, 1.9, 1.9, 0.5)}, which is slightly smaller than the R-defaults.
}
  \item{mgp}{
   numerical vector of length 3, indicating the axis and label positions. See \code{\link[graphics]{par}} for details. The default value is \code{c(2, 0.6, 0)}, which is slightly smaller than the R-defaults.
}
  \item{simobj}{object of class \code{svsim} as returned by the SV simulation function \code{\link{svsim}}. If provided, ``true'' data generating values will be added to the plots.
}
  \item{\dots}{
   further arguments are passed on to the invoked \code{plot} function.
}
}
\details{
\code{paradensplot} is modeled after \code{\link[coda]{densplot}} in the \code{coda} package, with some modifications for parameters that have (half-)bounded support.
}

\value{
 Called for its side effects. Returns argument \code{x} invisibly.
}

\note{
 You can call this function directly, but it is more commonly called by the \code{\link{plot.svdraws}} method.
}

\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}

\seealso{
\code{\link{paratraceplot}}, \code{\link{volplot}}, \code{\link{plot.svdraws}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
