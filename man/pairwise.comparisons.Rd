\name{pairwise.comparisons}
\alias{pairwise.comparisons}
\title{Description: Perform pairwise comparison between factor levels}
\usage{
  pairwise.comparisons(x, y, qth = 0.05, resdir = NULL)
}
\arguments{
  \item{x}{data matrix: samples x features}

  \item{y}{factor levels for the samples}

  \item{qth}{q-value threshold for included features}

  \item{resdir}{Optional. If given, output files will be
  written in the result directory

  Returns:}
}
\value{
  List with top findings from pairwise comparisons and
  their q-values
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

