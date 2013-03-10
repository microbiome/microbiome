\name{hclust.significance}
\alias{hclust.significance}
\title{Description: asses significance of hierarchical clustering
clusters based on multiscale bootstrap resamping}
\usage{
  hclust.significance(dat, my.groups, R,
    sample.sizes = c(0.2, 0.5, 1, 2), min.size, corr.th,
    replace = TRUE)
}
\arguments{
  \item{dat}{data matrix samples to be clustered on the
  rows}

  \item{my.groups}{groups from hierarchical clustering}

  \item{R}{bootstrap iterations}

  \item{sample.sizes}{sample sizes for multiscale bootstrap
  as fractions of the total data (0...1)}

  \item{min.size}{minimum cluster size}

  \item{corr.th}{correlation threshold for determining
  hclust groups}

  \item{replace}{sampling with replacement TRUE/FALSE}
}
\value{
  p-value corresponding to each input cluster in my.groups
  argument
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

