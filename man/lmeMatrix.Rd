\name{lmeMatrix}
\alias{lmeMatrix}
\title{Description: Computes the given lme model for each variable (row) in the given data frame
     Designed for HITchip matrices, but is applicable to any other matrix also.
     NOTE: does not take log automatically!}
\usage{
  lmeMatrix(fixed, vars, d.matrix, d.info, random,
    correlation = NULL, weights = NULL, subset,
    method = c("REML", "ML"), na.action = na.fail,
    control = list(), contrasts = NULL, keep.data = TRUE)
}
\arguments{
  \item{fixed}{TBA}

  \item{vars}{TBA}

  \item{d.matrix}{TBA}

  \item{d.info}{TBA}

  \item{random}{TBA}

  \item{correlation}{TBA}

  \item{weights}{TBA}

  \item{subset}{TBA}

  \item{method}{REML or ML}

  \item{na.action}{TBA}

  \item{control}{TBA}

  \item{contrasts}{TBA}

  \item{keep.data}{TBA

  Returns:}
}
\value{
  lmeMatrixRes
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

