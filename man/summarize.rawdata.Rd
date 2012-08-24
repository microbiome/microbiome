\name{summarize.rawdata}
\alias{summarize.rawdata}
\title{Description: Compress normalized raw data matrix into final probe-level matrix:
             summarize oligos into probes and hybridisations into samples}
\usage{
  summarize.rawdata(fdat, fdat.hybinfo, fdat.oligoinfo,
    oligo.ids)
}
\arguments{
  \item{fdat}{normalized raw data matrix oligos x
  hybridisations}

  \item{fdat.hybinfo}{hybridization info table}

  \item{fdat.oligoinfo}{oligo info table}

  \item{oligo.ids}{oligo.ids Returns:}
}
\value{
  probes x samples matrix
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

