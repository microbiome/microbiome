\name{outlierPvalue}
\alias{outlierPvalue}
\title{Statistical test for outliers in a vector}
\usage{
  outlierPvalue(x, avgvar)
}
\arguments{
  \item{x}{numerical vector}

  \item{avgvar}{the estimated variance of the distribution
  Returns:}
}
\value{
  Returns the p-value of the outlier test
}
\description{
  Function for determining the significances (p-values) of
  outliers in vector. Uses chisq.out.test from the package
  'outliers'.
}
\details{
  Arguments:
}
\examples{
# TBA
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

