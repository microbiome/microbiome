\name{diversity.table}
\alias{diversity.table}
\title{diversity.table}
\usage{
diversity.table(dat, level.from, level.to, phylogeny.info = NULL,
  diversity.index = "shannon", det.th = 0, min.probes = 0)
}
\arguments{
  \item{dat}{data matrix (phylotypes x samples) in original
  (non-log) scale; this should present the highest level
  phylogeny (typically probe-level) used in the analysis}

  \item{phylogeny.info}{mapping table between taxonomic
  levels}

  \item{level.from}{higher-level taxonomic groups for which
  the diversity will be calculated}

  \item{level.to}{lower taxonomic level used for the
  diversity calculations. Must correspond to the input data
  matrix argument 'dat'}

  \item{diversity.index}{diversity index (shannon or
  invsimpson)}

  \item{det.th}{Optional detection threshold.}

  \item{min.probes}{Minimum number of probes for a
  phylotype

  Returns:}
}
\value{
Table with various richness, evenness, and diversity
indicators
}
\description{
Description: Estimate diversity within each taxonomic group
across the samples
}
\details{
Arguments:
}
\examples{
# NOT RUN divtab <- diversity.table(dat, level.from = "L2", level.to = "oligo")
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

