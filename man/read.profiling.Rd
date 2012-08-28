\name{read.profiling}
\alias{read.profiling}
\title{read.profiling}
\usage{
  read.profiling(level, method = "sum", data.dir = NULL,
    log10 = TRUE)
}
\arguments{
  \item{level}{phylogenetic level ("oligo" / "species" /
  "L1" / "L2" / "L0") or "oligomap"}

  \item{method}{("rpa" / "sum" / "ave" / "nmf")}

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
  Description: read run.profiling.script output into R
}
\details{
  Arguments:
}
\examples{
# params <- run.profiling.script(...); dat <- read.profiling.data(params$wdir, "L1", "rpa")
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

