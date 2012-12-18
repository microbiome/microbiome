\name{summarize.probesets}
\alias{summarize.probesets}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets(phylogeny.info, oligo.data, method,
    level, verbose = TRUE, rm.phylotypes = NULL)
}
\arguments{
  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{level}{summarization level}

  \item{verbose}{print intermediate messages}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, L1, L2) Returns:}
}
\value{
  summarized data matrix
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

