\name{deconvolution.nonneg}
\alias{deconvolution.nonneg}
\title{Description: Deconvolution}
\usage{
  deconvolution.nonneg(oligo.data, phylogeny.info, level,
    block.solution = T, verbose = F, ...)
}
\arguments{
  \item{oligo.data}{oligo.data in absolute domain}

  \item{phylogeny.info}{oligo - phylotype mapping data
  frame}

  \item{level}{taxonomic level}

  \item{block.solution}{block.solution}

  \item{verbose}{verbose}

  \item{...}{parameters to be passed to ngp function

  Returns:}
}
\value{
  list
}
\description{
  !OBSOLETE! Used nmf package (currently removed from CRAN)
  for solving cross-hyb problem. Calls function ngp with
  standard arguments (SEE: ngp).
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

