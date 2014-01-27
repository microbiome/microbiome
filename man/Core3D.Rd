\name{Core3D}
\alias{Core3D}
\title{Core3D}
\usage{
Core3D(coreMat, title = "Core microbiota", cex.axis = 0.7)
}
\arguments{
  \item{coreMat}{core matrix}

  \item{title}{title}

  \item{cex.axis}{axis text size

  Returns:}
}
\value{
Used for its side effects
}
\description{
Description: Core visualization 3D
}
\details{
Arguments:
}
\examples{
data(peerj32); c3d <- Core3D(createCore(t(peerj32$microbes)))
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

