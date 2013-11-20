\name{densityplot}
\alias{densityplot}
\title{DensityPlot}
\usage{
  densityplot(mat, title = NULL, x.ticks = 10,
    rounding = 0, add.points = TRUE, point.color = "red",
    adjust = 1)
}
\arguments{
  \item{mat}{Data matrix to plot. The first two columns
  will be visualized as a cross-plot.}

  \item{title}{title text}

  \item{x.ticks}{Number of ticks on the X axis}

  \item{rounding}{Rounding for X axis tick values}

  \item{add.points}{Plot the data points as well}

  \item{point.color}{Color of the data points}

  \item{adjust}{Kernel width adjustment

  Returns:}
}
\value{
  ggplot2 object
}
\description{
  Description: Plots densities of data points in addition
  to cross-plot points.
}
\details{
  Arguments:
}
\examples{
p <- densityplot(cbind(rnorm(100), rnorm(100)))
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

