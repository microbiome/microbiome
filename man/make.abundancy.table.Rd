\name{make.abundancy.table}
\alias{make.abundancy.table}
\title{Description: Calculate abundancies
Discretize Hitchip matrix to form abundancy table
of form j, nj where j is number of counts and nj is number
of phylotypes with the corresponding counts
this format is often required by richness estimation}
\usage{
  make.abundancy.table(dat, det.th,
    discretization.resolution = 1)
}
\arguments{
  \item{dat}{data matrix}

  \item{det.th}{detection threshold}

  \item{discretization.resolution}{discretization
  resolution

  Returns:}
}
\value{
  abundancy table
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

