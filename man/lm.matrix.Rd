\name{lm.matrix}
\alias{lm.matrix}
\title{Description:  Computes the given lm model for each variable (row) in the given matrix/data frame
 		 Designed for HITchip matrices, but is applicable to any other matrix also.
 		 NOTE: does not take log automatically!
FIXME: merge with lm.matrix2}
\usage{
  lm.matrix(formula, d.matrix, d.info, type = "coef")
}
\arguments{
  \item{formula}{formula}

  \item{d.matrix}{data matrix}

  \item{d.info}{information}

  \item{type}{information to return

  Returns:}
}
\value{
  matrix
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

