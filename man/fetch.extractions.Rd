\name{fetch.extractions}
\alias{fetch.extractions}
\title{fetch.extractions
Fetch extractions from the phyloarray database.}
\usage{
  fetch.extractions(con, condition = NULL)
}
\arguments{
  \item{con}{MySQL connection}

  \item{condition}{list of lists with field-value pairs}
}
\value{
  An \code{extractionList} object consisting of a list with
  one entry called "extractions" which is a dataframe
  consisting of the records obtained from the database.
}
\description{
  The function fetches complete records from the
  \emph{extractions} table, including some fields from the
  array and hybridisation tables. If \code{condition=NULL}
  it will fetch all records and if \code{condition} is
  defined it will fetch only those records that comply with
  the condition.  Condition must be a list of lists, each
  having at least the fields \emph{field} and \emph{value}.
  The \emph{field} field must be a character vector of
  length 1 and \emph{value} must be a vector.  Each of the
  values will be evaluated as a optional value for the
  \emph{field} field.  This is an example:
  "\code{list(list(field='projectName',value=c(A,B)))}"
  that will be evaluated to the SQL condition "\code{WHERE
  (projectName='A' OR projectName='B')}".
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

