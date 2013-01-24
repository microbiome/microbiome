\name{ScaleProfile}
\alias{ScaleProfile}
\title{Description: Between-arrays normalization}
\usage{
  ScaleProfile(dat, method = "minmax", bg.adjust = NULL,
    minmax.quantiles = c(0.005, 0.995),
    minmax.points = NULL)
}
\arguments{
  \item{dat}{data matrix in original absolute scale}

  \item{method}{normalization method}

  \item{bg.adjust}{background adjustment}

  \item{minmax.quantiles}{quantiles for minmax}

  \item{minmax.points}{minmax end points

  Returns:}
}
\value{
  Normalized data matrix in absolute scale
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

