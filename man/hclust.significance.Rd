\name{hclust.significance}
\alias{hclust.significance}
\title{Description: asses significance of hierarchical clustering
clusters based on multiscale bootstrap resamping}
\usage{
  hclust.significance(dat, my.groups = NULL, R,
    sample.sizes = 1, min.size, corr.th, replace = TRUE,
    metric = "pearson", verbose = TRUE,
    pvalue.threshold = 0.05, remove.nested.clusters = TRUE)
}
\arguments{
  \item{dat}{data matrix samples to be clustered on the
  rows}

  \item{my.groups}{groups from hierarchical clustering}

  \item{R}{bootstrap iterations}

  \item{sample.sizes}{sample sizes for multiscale bootstrap
  as fractions of the total data (0...1..etc.). Use vector
  of different sizes to perform multiscale bootstrap
  resampling.}

  \item{min.size}{minimum cluster size}

  \item{corr.th}{correlation threshold for determining
  hclust groups}

  \item{replace}{sampling with replacement TRUE/FALSE}

  \item{metric}{correlation metrix}

  \item{verbose}{verbose}

  \item{pvalue.threshold}{pvalue threshold: filter out
  non-significant groups}

  \item{remove.nested.clusters}{Remove nested clusters, or
  keep all? Default: TRUE}
}
\value{
  List of significant clusters (corresponding to the given
  p-value)
}
\description{
  Description: asses significance of hierarchical
  clustering clusters based on multiscale bootstrap
  resamping
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{methods}

