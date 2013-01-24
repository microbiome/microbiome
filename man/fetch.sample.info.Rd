\name{fetch.sample.info}
\alias{fetch.sample.info}
\title{fetch.sample.info}
\usage{
  fetch.sample.info(allowed.projects, chiptype = NULL,
    dbuser, dbpwd, dbname, selected.samples = NULL,
    host = NULL, port = NULL)
}
\arguments{
  \item{allowed.projects}{list projects for which to fetch
  the data}

  \item{chiptype}{chiptype (eg. new.chip)}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{selected.samples}{Sample to investigate. By default
  all.}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections Returns:}
}
\value{
  project.info data.frame
}
\description{
  Description: Fetch sample information from HITChip atlas
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

