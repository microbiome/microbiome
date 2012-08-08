\name{scaling.minmax}
\alias{scaling.minmax}
\title{Description: Minmax scaling.}
\usage{
  scaling.minmax(r, quantile.points, robust = FALSE)
}
\arguments{
  \item{r}{data matrix in log10 scale}

  \item{quantile.points}{quantiles for minmax}

  \item{robust}{Select minmax version.

  Returns:}
}
\value{
  normalized data matrix
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
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

