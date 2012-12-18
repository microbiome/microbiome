\name{preprocess.chipdata}
\alias{preprocess.chipdata}
\title{Description: Profiling preprocessing script}
\usage{
  preprocess.chipdata(dbuser, dbpwd, dbname, mc.cores = 1,
    verbose = TRUE, host = NULL, port = NULL)
}
\arguments{
  \item{dbuser}{MySQL username}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySQL database name}

  \item{mc.cores}{Optional. Number of cores if
  parallelization is used.}

  \item{verbose}{monitor processing through intermediate
  messages}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections

  Returns:}
}
\value{
  Preprocessed data and parameters
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

