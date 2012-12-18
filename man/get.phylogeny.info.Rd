\name{get.phylogeny.info}
\alias{get.phylogeny.info}
\title{get.phylogeny.info}
\usage{
  get.phylogeny.info(phylogeny = "16S", rmoligos = NULL,
    dbuser, dbpwd, dbname, verbose = TRUE,
    remove.nonspecific.oligos = FALSE, chip = "HITChip",
    host = NULL, port = NULL)
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

  \item{chip}{chip type}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections Returns:}
}
\value{
  phylogeny.info
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

