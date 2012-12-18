\name{get.probedata}
\alias{get.probedata}
\title{Description: Get probedata}
\usage{
  get.probedata(hybridization.ids, rmoligos, dbuser, dbpwd,
    dbname, mc.cores = 1, host = NULL, port = NULL)
}
\arguments{
  \item{hybridization.ids}{Specify the hybridizations to
  retrieve}

  \item{rmoligos}{oligos to exclude}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{mc.cores}{Number of cores for multicore computing}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections Returns:}
}
\value{
  list with data (features x hybridizations matrix) and
  info (features x info) fields
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

