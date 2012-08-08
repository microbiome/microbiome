\name{get.phylogeny}
\alias{get.phylogeny}
\title{Description: Get phylogeny}
\usage{
  get.phylogeny(phylogeny, rmoligos = NULL, dbuser, dbpwd,
    dbname, verbose = TRUE,
    remove.nonspecific.oligos = FALSE, chip = "HITChip")
}
\arguments{
  \item{phylogeny}{phylogeny}

  \item{rmoligos}{oligos to exclude}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{verbose}{verbose}

  \item{remove.nonspecific.oligos}{Logical. Remove oligos
  with multiple targets.}

  \item{chip}{chip type Returns:}
}
\value{
  oligo.map
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

