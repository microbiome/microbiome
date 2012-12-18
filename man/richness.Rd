\name{richness}
\alias{richness}
\title{richness}
\usage{
  richness(dat, det.th = NULL)
}
\arguments{
  \item{dat}{data matrix (phylotypes x samples) in original
  (non-log) scale}

  \item{det.th}{detection threshold. Used for richness and
  evenness estimation. Not used in diversity estimation.

  Returns:}
}
\value{
  Vector containing richness estimate for each sample
}
\description{
  Description: Estimate richness for each sample with a
  given threshold
}
\details{
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

