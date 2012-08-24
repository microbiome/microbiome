\name{ScaleProfile}
\alias{ScaleProfile}
\title{Description: Between-arrays normalization}
\usage{
  ScaleProfile(r.feature, method = "minmax",
    bg.adjust = NULL, minmax.quantiles = c(0.005, 0.995))
}
\arguments{
  \item{r.feature}{data matrix in logarithmic scale}

  \item{method}{normalization method}

  \item{bg.adjust}{background adjustment}

  \item{minmax.quantiles}{quantiles for minmax Returns:}
}
\value{
  Normalized data matrix
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

