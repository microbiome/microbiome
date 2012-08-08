\name{estimate.diversity}
\alias{estimate.diversity}
\title{Description: Estimate diversities for each sample (column)
Aliases: get.diversity.estimates. Also replaces the function ST: diversity.indices}
\usage{
  estimate.diversity(dat, diversity.index = "shannon",
    det.th = NULL)
}
\arguments{
  \item{dat}{data matrix (phylotypes x samples) in original
  (non-log) scale}

  \item{diversity.index}{diversity index}

  \item{det.th}{detection threshold. Used for richness and
  evenness estimation. Not used in diversity estimation.

  Returns:}
}
\value{
  diversity indicators
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

