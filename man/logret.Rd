\name{logret}
\alias{logret}
\title{
Computes (de-meaned) log returns.
}
\description{
Small utlity function returning either \code{diff(log(x))} in case the argument \code{demean} is set to \code{FALSE}, or \code{diff(log(x)) - mean(diff(log(x)))} in case that \code{demean} is \code{TRUE}.
}
\usage{
logret(x, demean = FALSE)
}
\arguments{
  \item{x}{Real-valued vector.}
  \item{demean}{A single logical value indicating whether the returns should be de-meaned. Defaults to \code{FALSE}.}
}
\value{
 A vector of length \code{length(x) - 1}, containing (de-meaned) returns.
}
\author{
Gregor Kastner \email{gregor.kastner@wu.ac.at}
}

\keyword{utilities}
