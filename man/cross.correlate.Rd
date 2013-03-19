\name{cross.correlate}
\alias{cross.correlate}
\title{Description: Cross-correlate columns of the input matrices}
\usage{
  cross.correlate(x, y, method = "pearson", qth = NULL,
    cth = NULL, order = FALSE, n.signif = 0,
    mode = "table", qvalues = TRUE)
}
\arguments{
  \item{x}{matrix (samples x features if annotation
  matrix)}

  \item{y}{matrix (samples x features if cross-correlated
  with annotations)}

  \item{method}{association method (pearson, spearman for
  continuous; categorical for discrete)}

  \item{qth}{q-value threshold to include features}

  \item{cth}{correlation threshold to include features}

  \item{order}{order the results}

  \item{n.signif}{mininum number of significant
  correlations for each element}

  \item{mode}{Specify output format ("table" or "matrix")}

  \item{qvalues}{calculate qvalues

  Returns:}
}
\value{
  List with cor, pval, qval
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

