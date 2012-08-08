\name{fetch.sample.info}
\alias{fetch.sample.info}
\title{Description: Fetch sample information from HITChip atlas}
\usage{
  fetch.sample.info(allowed.projects, chiptype = NULL,
    dbuser, dbpwd, dbname, selected.samples = NULL)
}
\arguments{
  \item{allowed.projects}{list projects for which to fetch
  the data}

  \item{chiptype}{chiptype (eg. new.chip)}

  \item{dbuser}{MySQL user}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySqL database name}

  \item{selected.samples}{Sample to investigate. By default
  all. Returns:}
}
\value{
  project.info data.frame
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

