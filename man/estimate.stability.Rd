\name{estimate.stability}
\alias{estimate.stability}
\title{Description: Stability analysis. Calculates average Pearson '
 correlation between samples in the input data and picks the lower '
 triangular matrix to avoid duplicating the correlations. Returns
 correlations and stability estimate (average of the correlations).
 Can also be used to calculate stability between two data sets.
 Then provide two data sets as inputs.}
\usage{
<<<<<<< HEAD:man/estimate.stability.Rd
  estimate.stability(dat1, dat2 = NULL, method = "pearson")
=======
calculate.stability(dat1, dat2 = NULL, method = "pearson")
>>>>>>> 405f396f9007cea42de70fa994a387fffaea0aeb:man/calculate.stability.Rd
}
\arguments{
  \item{dat1}{data matrix phylotypes vs. samples (in log10
  scale)}

  \item{dat2}{Optional. Second data matrix phylotypes vs.
  samples.  Provide this to calculate stability between two
  (paired) data sets.}

  \item{method}{Correlation method (see ?cor)

  Returns:}
}
\value{
List with correlations and astability estimate
}
\description{
Arguments:
}
\examples{
data(peerj32); s <- estimate.stability(t(peerj32$microbes)[, 1:5])
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

