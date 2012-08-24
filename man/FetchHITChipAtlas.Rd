\name{FetchHITChipAtlas}
\alias{FetchHITChipAtlas}
\title{FetchHITChipAtlas}
\usage{
  FetchHITChipAtlas(allowed.projects, dbuser, dbpwd,
    dbname, result.path, my.scaling = "minmax",
    minmax.quantiles = c(0.005, 0.995), bgc.method = NULL,
    mc.cores = 3, remove.nonspecific.oligos = FALSE)
}
\arguments{
  \item{allowed.projects}{Vector of projectNames to
  include}

  \item{dbuser}{MySQL username}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySQL database name}

  \item{result.path}{Directory for storing the results.}

  \item{my.scaling}{Normalization method. Default
  "minmax".}

  \item{minmax.quantiles}{Quantiles for minmax
  normalization. Default c(0.005, 0.995).}

  \item{bgc.method}{Optional. Background correction method.
  By default no background correction.}

  \item{mc.cores}{Number of cores for parallel computation}

  \item{remove.nonspecific.oligos}{Optional logical set to
  remove nonspecific oligos (TRUE). By default keeping the
  non-specific oligos (FALSE) Returns:}
}
\value{
  List with the following elements: training.data: training
  data matrices (including alternative preprocessing
  methods) test.data: test data matrices full.data: full
  data matrices atlas.sampleinfo: sample information matrix
  parameters: input parameters and R session information
}
\description{
  Description: Complete preprocessing pipeline for the
  HITChip atlas with default parameters
}
\details{
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

