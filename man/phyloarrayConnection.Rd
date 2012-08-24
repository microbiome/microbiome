\name{phyloarrayConnection}
\alias{phyloarrayConnection}
\title{Tests whether the database connection is a phyloarray connection.
Expands one element (one "field", "value" pair list) from a list
of "field", "value" pair lists}
\usage{
  phyloarrayConnection(con)
}
\arguments{
  \item{con}{a MySQL database connection.}
}
\value{
  TRUE when the test succeeds. Otherwise a program halt.
}
\description{
  Tests whether the database connection is a phyloarray
  connection. Expands one element (one "field", "value"
  pair list) from a list of "field", "value" pair lists
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

