\name{diversity}
\alias{diversity}
\title{diversity}
\usage{
diversity(dat, diversity.index = "shannon", det.th = 0)
}
\arguments{
  \item{dat}{data matrix (phylotypes x samples) in original
  (non-log) scale}

  \item{diversity.index}{diversity index (shannon or
  invsimpson)}

  \item{det.th}{detection threshold. Used for richness and
  evenness estimation. Not used in diversity estimation.

  Returns:}
}
\value{
Vector containing diversity estimate for each sample
}
\description{
Description: Estimate diversity for each sample with a
given threshold
}
\details{
Arguments:
}
\examples{
data(peerj32); div <- diversity(10^t(peerj32$microbes), det.th = 0)
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

