\name{summarize.probesets.species}
\alias{summarize.probesets.species}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets.species(phylogeny.info, oligo.data,
    method, verbose = TRUE,
    rm.species = c("Victivallis vadensis"))
}
\arguments{
  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{verbose}{print intermediate messages}

  \item{rm.species}{Species to exclude Returns:}
}
\value{
  summarized data matrix in log10 scale
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

