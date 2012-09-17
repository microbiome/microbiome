\name{cross.correlate}
\alias{cross.correlate}
\title{Description: Cross-correlate input variables}
\usage{
  cross.correlate(annot, dat, method = "pearson",
    qth = NULL, cth = NULL, order = FALSE, n.signif = 0,
    mode = "table")
}
\arguments{
  \item{annot}{annotation matrix: samples x features}

  \item{dat}{numeric data matrix: samples x features}

  \item{method}{association method (pearson, spearman for
  continuous; categorical for discrete}

  \item{qth}{q-value threshold for included features}

  \item{cth}{correlation threshold}

  \item{order}{order the results}

  \item{n.signif}{mininum number of significant
  correlations for each element}

  \item{mode}{Specify the output format ("table" or
  "matrix")

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

