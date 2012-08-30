\name{get.oligomap}
\alias{get.oligomap}
\title{get.oligomap}
\usage{
  get.oligomap(phylogeny = "16S", rmoligos = NULL, dbuser,
    dbpwd, dbname, verbose = TRUE,
    remove.nonspecific.oligos = FALSE, chip = "HITChip")
}
\arguments{
  \item{phylogeny}{phylogeny (default: 16S)}

  \item{rmoligos}{oligos to exclude}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{verbose}{verbose}

  \item{remove.nonspecific.oligos}{Logical. Remove oligos
  with multiple targets.}

  \item{chip}{chip type

  Returns:}
}
\value{
  oligomap
}
\description{
  Description: Get phylogeny
}
\details{
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

