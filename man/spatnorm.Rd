\name{spatnorm}
\alias{spatnorm}
\title{2D Spatial normalisation of microarrays}
\usage{
  spatnorm(data, x = "x", y = "y", signal = "signal",
    method = c("loess"), span = 0.03, degree = 2,
    family = "symmetric", subset = 1 == 1)
}
\arguments{
  \item{data}{dataframe with positional and signal
  information}

  \item{x}{name of the x column in the dataframe}

  \item{y}{name of the y column in the dataframe}

  \item{signal}{name of the signal column in the dataframe}

  \item{method}{method name, currently only "loess" is
  supported}

  \item{span}{span}

  \item{degree}{degree}

  \item{family}{family}

  \item{subset}{logical expression indicating rows to use
  for normalisation}
}
\value{
  TBA
}
\description{
  2D Spatial normalisation of microarrays
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

