\name{ngp}
\alias{ngp}
\title{Description: ngp: Non-negative probabilistic solution to probe mixing problem.}
\usage{
  ngp(oligo.data, phylogeny.info, level, lambda = 0.001,
    alpha = 1, beta = 1)
}
\arguments{
  \item{oligo.data}{oligo.data in absolute domain}

  \item{phylogeny.info}{oligo - phylotype mapping data
  frame}

  \item{level}{taxonomic level}

  \item{lambda}{- stregth of gamma prior. Default: 0.001.}

  \item{alpha}{- alpha parameter of gamma prior. Default:
  1.}

  \item{beta}{- beta parameter of gamma prior. Default: 1.

  Returns:}
}
\value{
  list
}
\description{
  Uses gamma prior to solve known cross-hybridisation
  problems. Perfect probe targets are given in
  phylogeny.info data frame.
}
\details{
  Arguments:
}
\author{
  Contact: Jarkko Salojarvi
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

