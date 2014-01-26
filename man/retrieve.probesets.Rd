\name{retrieve.probesets}
\alias{retrieve.probesets}
\title{retrieve.probesets}
\usage{
  retrieve.probesets(phylogeny.info, level = "species",
    name = NULL)
}
\arguments{
  \item{phylogeny.info}{data.frame with oligo - phylotype
  mapping information}

  \item{level}{phylotype level for probesets}

  \item{name}{specify phylotypes to check (optional)

  Returns:}
}
\value{
  A list. Probes for each phylotype.
}
\description{
  Description: List probes for each probeset
}
\details{
  Arguments:
}
\examples{
# NOT RUN phylogeny.info <- GetPhylogeny(); sets <- retrieve.probesets(phylogeny.info, level = "species", name = "Weissella confusa")
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

