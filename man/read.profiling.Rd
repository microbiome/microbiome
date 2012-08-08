\name{read.profiling}
\alias{read.profiling}
\title{Description: Function to read run.profiling.script output into R}
\usage{
  read.profiling(level, method, data.dir, log10 = FALSE)
}
\arguments{
  \item{level}{phylogenetic level ("oligo" / "species" /
  "level1" / "level2" / "level0") or "phylogeny"}

  \item{method}{("rpa" / "sum")}

  \item{data.dir}{Profiling script output directory for
  reading the data. If not given, GUI will ask to specify
  the file and overruns the possible level / method
  arguments in the function call.}

  \item{log10}{Logical. Logarithmize the data TRUE/FALSE.
  By default, the data is in original non-log scale.

  Returns:}
}
\value{
  data matrix (phylo x samples)
}
\description{
  Arguments:
}
\examples{
# params <- run.profiling.script(...); dat <- read.profiling.data(params$wdir, "level1", "rpa")
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

