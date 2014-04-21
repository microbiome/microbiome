\name{cross.correlate}
\alias{cross.correlate}
\title{Description: Cross-correlate columns of the input matrices}
\usage{
cross.correlate(x, y = NULL, method = "pearson", p.adj.threshold = Inf,
  cth = NULL, order = FALSE, n.signif = 0, mode = "table",
  p.adj.method = "fdr", verbose = F, filter.self.correlations = F)
}
\arguments{
  \item{x}{matrix (samples x features if annotation
  matrix)}

  \item{y}{matrix (samples x features if cross-correlated
  with annotations)}

  \item{method}{association method (pearson, spearman for
  continuous; categorical for discrete)}

  \item{p.adj.threshold}{q-value threshold to include
  features}

  \item{cth}{correlation threshold to include features}

  \item{order}{order the results}

  \item{n.signif}{mininum number of significant
  correlations for each element}

  \item{mode}{Specify output format ("table" or "matrix")}

  \item{p.adj.method}{p-value multiple testing correction
  method. Either "qvalue" or one of the methods in p.adjust
  function ("BH" and others; see help(p.adjust)). Default:
  "fdr"}

  \item{verbose}{verbose}

  \item{filter.self.correlations}{Filter out correlations
  between identical items.

  Returns:}
}
\value{
List with cor, pval, pval.adjusted
}
\description{
Arguments:
}
\details{
As the method=categorical (discrete) association measure
for nominal (no order for levels) variables ' we using
Goodman and Kruskal tau based on
http://www.r-bloggers.com/measuring-associations-between-non-numeric-variables/
}
\examples{
data(peerj32); cc <- cross.correlate(peerj32$microbes[1:20, 1:10], peerj32$lipids[1:20,1:10])
}
\author{
Contact: Leo Lahti
\email{microbiome-admin@googlegroups.com}
}
\references{
See citation("microbiome")
}
\keyword{utilities}

