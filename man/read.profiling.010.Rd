\name{read.profiling.010}
\alias{read.profiling.010}
\title{read.profiling.010}
\usage{
  read.profiling.010(level = NULL, method = "rpa",
    data.dir = NULL, log10 = TRUE, impute = TRUE)
}
\arguments{
  \item{level}{phylogenetic level ("oligo" / "species" /
  "L1" / "L2" / "L0") or "phylogeny.info"}

  \item{method}{("rpa" / "sum" / "ave" / "nmf")}

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
  Description: read profiling script output from profiling
  script v. 010 into R
}
\details{
  Arguments:
}
\examples{
# params <- read.profiling.010(...);
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

