\name{evenness}
\alias{evenness}
\title{evenness}
\usage{
evenness(dat, det.th = NULL)
}
\arguments{
  \item{dat}{data matrix (phylotypes x samples) in original
  (non-log) scale}

  \item{det.th}{detection threshold. Used for richness and
  evenness estimation. Not used in diversity estimation.

  Returns:}
}
\value{
Vector containing evenness estimate for each sample
}
\description{
Description: Estimate evenness for each sample with a given
threshold
}
\details{
Arguments:
}
\examples{
data(peerj32); eve <- evenness(10^t(peerj32$microbes), det.th = 100)
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

