\name{rsd}
\alias{rsd}
\title{Relative standard deviation using paired measurements}
\usage{
  rsd(x, y = NULL)
}
\arguments{
  \item{x}{dataframe, matrix or vector}

  \item{y}{dataframe, matrix or vector compatible with x}
}
\value{
  Returns a matrix of rsd values.
}
\description{
  Calculates an unbiased estimate of the the relative
  squared difference of x.  Although most RSD calculations
  use the following formula: mean(abs(x_i-y_i)/(x_i+y_i))
  we here use: mean(sqrt(pi)*abs(x_i-y_i)/(x_i+y_i)) The
  latter is an unbiased estimate of the relative standard
  deviation (sigma/mu) when data are normal distributed.
}
\examples{
# TBA
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

