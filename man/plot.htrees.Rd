\name{plot.htrees}
\alias{plot.htrees}
\title{plot.htrees
Description: Plot hierarchical clustering for the input data in absolute
and log10 scale using euclidean and pearson correlation similarities.
Intended for internal use in the run.profiling.script function.}
\usage{
  plot.htrees(dat, method = "complete")
}
\arguments{
  \item{dat}{oligoprofile data in original (non-log)
  domain}

  \item{method}{hierarchical clustering method Returns:}
}
\value{
  NULL
}
\description{
  Arguments:
}
\examples{
# dat <- read.profiling(params$wdir, "species", "rpa"); tmp <- plot.htrees(dat)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

