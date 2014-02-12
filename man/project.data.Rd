\name{project.data}
\alias{project.data}
\title{Description: Project high-dimensional data on two-dimensional plane by various methods}
\usage{
project.data(amat, type = "PCA")
}
\arguments{
  \item{amat}{data matrix (samples x features)}

  \item{type}{projection type (options: PCA, MDS.classical,
  MDS.nonmetric, Sammon) Returns:}
}
\value{
projected data matrix
}
\description{
Arguments:
}
\examples{
data(peerj32); xy <- project.data(peerj32$microbes[,1:3])
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

