\name{calculate.stability}
\alias{calculate.stability}
\title{Description: Stability analysis. Calculates average Pearson '
 correlation between samples in the input data and picks the lower '
 triangular matrix to avoid duplicating the correlations. Returns
 correlations and stability estimate (average of the correlations).
 Can also be used to calculate stability between two data sets.
 Then provide two data sets as inputs.}
\usage{
  calculate.stability(dat1, dat2 = NULL,
    method = "pearson")
}
\arguments{
  \item{dat1}{data matrix phylotypes vs. samples}

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
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

