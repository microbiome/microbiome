\name{importProjectFromPhyloarrayToR}
\alias{importProjectFromPhyloarrayToR}
\title{Description: Import Project From Phyloarray To R. This function imports a
project's ln(data) from probeLevel from phyloarray mysql-database to R.
Modified from avgSignalPerTaxon.sql Query 2 by JN
NOTE!!: Must be modified to average over the identical probes!!
The desired projectName is given as an argument
as well as the algorithm versions.}
\usage{
  importProjectFromPhyloarrayToR(projectName = "ProvasI",
    spatNormAlgv = "1.1", normAlgv = "1.1", dbuser = NULL,
    dbpwd = NULL, dbname = NULL)
}
\arguments{
  \item{projectName}{Project to import}

  \item{spatNormAlgv}{Spatial normalization algorithm
  version}

  \item{normAlgv}{Normalization algorithm version}

  \item{dbuser}{mysql database username}

  \item{dbpwd}{mysql database password}

  \item{dbname}{mysql database name

  Returns:}
}
\value{
  data
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

