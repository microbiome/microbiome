\name{cross.correlate}
\alias{cross.correlate}
\title{Description: Cross-correlate input variables}
\usage{
  cross.correlate(annot, dat, method = "pearson",
    qth = NULL, cth = NULL, order = FALSE, n.signif = 0,
    verbose = TRUE, mode = "list")
}
\arguments{
  \item{annot}{annotation matrix: samples x features}

  \item{dat}{data matrix: samples x features}

  \item{method}{association method (pearson, spearman for
  continuous; categorical for discrete}

  \item{qth}{q-value threshold for included features}

  \item{cth}{correlation threshold}

  \item{order}{order the results}

  \item{n.signif}{mininum number of significant
  correlations for each element}

  \item{verbose}{verbose

  Returns:}
}
\value{
  List with cor, pval, qval, N
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

