\name{scaling.minmax}
\alias{scaling.minmax}
\title{Description: Minmax scaling.}
\usage{
  scaling.minmax(dat, quantile.points = NULL,
    minmax.points = NULL, robust = FALSE)
}
\arguments{
  \item{dat}{data matrix in original absolute scale}

  \item{quantile.points}{quantiles for minmax}

  \item{minmax.points}{pre-calculated quantiles for minmax
  in log10 scale; overrides quantile.points argument}

  \item{robust}{Select minmax version.

  Returns:}
}
\value{
  normalized data matrix in absolute scale
}
\description{
  Arguments:
}
\note{
  With robust = FALSE, the standard minmax is carried out.
  This shifts and scales each array such that their min and
  max values are identical across arrays. The robust = TRUE
  will perform the scaling such that the upper quantiles
  and minimum values of the data match (instead of maximum
  values).
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

