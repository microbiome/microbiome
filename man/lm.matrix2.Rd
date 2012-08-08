\name{lm.matrix2}
\alias{lm.matrix2}
\title{Description:  Computes the given lm model for each variable (row) in the given matrix/data frame
 		 Designed for HITchip matrices, but is applicable to any other matrix also.
 		 NOTE: does not take log automatically!}
\usage{
  lm.matrix2(formula, d.matrix, d.info, contr = NULL,
    Cvar = NULL, fit.contrast = NULL)
}
\arguments{
  \item{formula}{formula}

  \item{d.matrix}{data matrix}

  \item{d.info}{information}

  \item{contr}{TBA}

  \item{Cvar}{TBA}

  \item{fit.contrast}{Returns:}
}
\value{
  lmMatrixRes
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

