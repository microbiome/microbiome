\name{fetch.measurements}
\alias{fetch.measurements}
\title{Fetch measurements associated with one or more extractions from a phyloarray database}
\usage{
  fetch.measurements(con, extrList,
    which.data = c("raw", "rawBGcor", "spatnorm", "qnorm"),
    transformation = "none")
}
\arguments{
  \item{con}{MySQL connection}

  \item{extrList}{An extractionList object obtained from
  fetch.extractions}

  \item{which.data}{One or more of raw, rawBGcor, spatnorm
  and qnorm}

  \item{transformation}{One of none, avg, log or avglog}
}
\value{
  A dataframe with the selected measurements
}
\description{
  The function fetches records from the
  \emph{featureMeasurements} table, including some fields
  from the array, hybridisation and featureExtraction
  tables. The \code{which.data} argument indicates which
  data to fetch from the database, where \code{raw} will
  fetch the raw data (\code{FGsignal}), \code{rawBGcor}
  will fetch the raw background corrected data
  (\code{FGsignal-BGsignal}), \code{spatnorm} will fetch
  the spatially normalized data and \code{qnorm} will fetch
  the sample-quantile normalized data. The
  \code{transformation} argument indicates which
  transformation to apply to the data indicated by the
  \code{which.data} argument. Only one type of
  transformation is allowed, where "\code{none}" will not
  apply a transformation, "\code{avg}" will take the values
  averaged over oligo's and featureExtractions including
  standard deviation, "\code{log}" will take the 10-base
  logarithm of te data and "\code{avglog}" will return the
  oligo and featureExtraction averaged values of the
  logarithm of the data, including standard deviation.
}
\examples{
# TBA
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

