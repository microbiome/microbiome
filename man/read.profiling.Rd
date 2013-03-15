\name{read.profiling}
\alias{read.profiling}
\title{read.profiling}
\usage{
  read.profiling(level = NULL, method = "frpa", data.dir,
    log10 = TRUE, impute = TRUE)
}
\arguments{
  \item{level}{phylogenetic level ("oligo" / "species" /
  "L1" / "L2" / "L0") or "phylogeny.info"}

  \item{method}{("frpa" / "rpa" / "sum" / "ave" / "nmf")}

  \item{data.dir}{Profiling script output directory for
  reading the data. If not given, GUI will ask to specify
  the file and overruns the possible level / method
  arguments in the function call.}

  \item{log10}{Logical. Logarithmize the data TRUE/FALSE.
  By default, the data is in original non-log scale.}

  \item{impute}{impute missing oligo signals

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
# params <- run.profiling.script(...); dat <- read.profiling(params$wdir, "L1", "rpa")
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

