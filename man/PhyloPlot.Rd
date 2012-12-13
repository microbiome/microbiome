\name{PhyloPlot}
\alias{PhyloPlot}
\title{Description:
Visualize the histogram, density estimate, and estimated Gaussian latent classes for a given vector}
\usage{
  PhyloPlot(vec, max.responses = 8, mixture.method = "bic",
    bic.threshold = 0, title.text = NULL, xlab.text = NULL,
    ylab.text = "Frequency", binwidth = 0.05)
}
\arguments{
  \item{vec}{numeric data vector}

  \item{max.responses}{maximum number of latent classes to
  consider}

  \item{mixture.method}{information criterion for
  estimating the number of latent classes}

  \item{bic.threshold}{threshold for adding new components
  (latent classes) in the mixture}

  \item{title.text}{title.text}

  \item{xlab.text}{xlab.text}

  \item{ylab.text}{ylab.text

  Returns:}
}
\value{
  List: ggplot2 object and NetResponse model object
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

