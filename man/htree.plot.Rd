\name{htree.plot}
\alias{htree.plot}
\title{htree.plot
Description: Plot hierarchical clustering for the input data in absolute
and log10 scale using euclidean and pearson correlation similarities.
Intended for internal use in the run.profiling.script function.}
\usage{
  htree.plot(dat, method = "complete", metric = "spearman")
}
\arguments{
  \item{dat}{oligoprofile data in original (non-log)
  domain}

  \item{method}{hierarchical clustering method}

  \item{metric}{clustering similarity measure Returns:}
}
\value{
  NULL
}
\description{
  Arguments:
}
\examples{
# dat <- read.profiling(params$wdir, "species", "rpa"); tmp <- htree.plot(dat)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

