\name{get.probedata}
\alias{get.probedata}
\title{Description: Get probedata}
\usage{
  get.probedata(hybridization.ids, rmoligos, dbuser, dbpwd,
    dbname, mc.cores = 1)
}
\arguments{
  \item{hybridization.ids}{Specify the hybridizations to
  retrieve}

  \item{rmoligos}{oligos to exclude}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{mc.cores}{Number of cores for multicore computing
  Returns:}
}
\value{
  list with data (features x hybridizations matrix) and
  info (features x info) fields
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

