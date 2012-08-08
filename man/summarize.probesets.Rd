\name{summarize.probesets}
\alias{summarize.probesets}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets(oligo.map, oligo.data, method, level,
    verbose = TRUE, rm.phylotypes = NULL, rm.oligos = NULL)
}
\arguments{
  \item{oligo.map}{oligo - phylotype matching data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{level}{summarization level}

  \item{verbose}{print intermediate messages}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, level 1, level 2)}

  \item{rm.oligos}{Oligos to remove Returns:}
}
\value{
  summarized data matrix
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

