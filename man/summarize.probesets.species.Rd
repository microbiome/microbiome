\name{summarize.probesets.species}
\alias{summarize.probesets.species}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets.species(oligo.map, oligo.data,
    method, verbose = TRUE,
    rm.species = c("Victivallis vadensis"))
}
\arguments{
  \item{oligo.map}{oligo - phylotype matching data.frame}

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
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

