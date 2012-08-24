\name{bootstrap.microbes}
\alias{bootstrap.microbes}
\title{bootstrap.microbes}
\usage{
  bootstrap.microbes(D, Nsample = NULL, minPrev = 2,
    Nboot = 100, I.thr = 1.8, ncore = 1)
}
\arguments{
  \item{D}{data}

  \item{Nsample}{sample size}

  \item{minPrev}{minimum prevalence}

  \item{Nboot}{bootstrap sample size}

  \item{I.thr}{threshold}

  \item{ncore}{number of nodes for parallelization

  Returns:}
}
\value{
  TBA
}
\description{
  Description: bootstrap.microbes
}
\details{
  Arguments:
}
\examples{
# data("toydata-hitchip"); bootstrap.microbes(phylum.matrix)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

