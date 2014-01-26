\name{levelmap}
\alias{levelmap}
\title{levelmap}
\usage{
  levelmap(phylotypes = NULL, level.from, level.to,
    phylogeny.info)
}
\arguments{
  \item{phylotypes}{phylotypes to convert; if NULL then
  considering all phylotypes in the phylogeny.info}

  \item{level.from}{level.from}

  \item{level.to}{level.to}

  \item{phylogeny.info}{phylogeny.info

  Returns:}
}
\value{
  mappings
}
\description{
  Description: map phylotypes between hierarchy levels
}
\details{
  Arguments:
}
\examples{
# NOT RUN: phylogeny.info <- GetPhylogeny(); levelmap(phylotypes = "Akkermansia", level.from = "L2", level.to = "L1", phylogeny.info)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

