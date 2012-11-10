\name{calculate.hclust}
\alias{calculate.hclust}
\title{calculate.hclust}
\usage{
  calculate.hclust(dat, method = "complete",
    metric = "correlation")
}
\arguments{
  \item{dat}{data matrix}

  \item{method}{hierarchical clustering method (see
  ?hclust)}

  \item{metric}{clustering metrics (euclidean /
  correlation)

  Returns:}
}
\value{
  hclust object for log10 and for absolute scale data
}
\description{
  Description: Calculate hierarchical clustering for
  standard selections in profiling script
}
\details{
  Arguments:
}
\examples{
# hc <- calculate.hclust(dat, "complete", "correlation")
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

