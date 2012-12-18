\name{fetch.projects}
\alias{fetch.projects}
\title{fetch.projects}
\usage{
  fetch.projects(con, condition = NULL)
}
\arguments{
  \item{con}{MySQL connection}

  \item{condition}{list of lists with field-value pairs

  Fetch complete records from the \emph{projects} table. If
  \code{condition=NULL}; will fetch all records and if
  \code{condition} is defined it will fetch only those
  records that comply with the condition.  Condition must
  be a list of lists, each having at least the fields
  \emph{field} and \emph{value}.  The \emph{field} field
  must be a character vector of length 1 and \emph{value}
  must be a vector.  Each of the values will be evaluated
  as a optional value for the \emph{field} field.  This is
  an example:
  "\code{list(list(field='projectName',value=c(A,B)))}"
  that will be evaluated to the SQL condition "\code{WHERE
  (projectName='A' OR projectName='B')}".}
}
\value{
  A dataframe with the selected records from the projects
  table
}
\description{
  Fetch projects from the phyloarray database
}
\examples{
# TBA
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

