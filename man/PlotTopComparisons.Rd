\name{PlotTopComparisons}
\alias{PlotTopComparisons}
\title{Description: Visualization of top findings from pairwise.comparisons function}
\usage{
  PlotTopComparisons(x, y, oligo.map, color.level,
    bar.level, top.findings, top.findings.qvals,
    qth = 0.05, qth.star = 0.01, mode = "barplot")
}
\arguments{
  \item{x}{data matrix}

  \item{y}{annotaction factor}

  \item{oligo.map}{mapping between features}

  \item{color.level}{feature level to color}

  \item{bar.level}{feature level for bars}

  \item{top.findings}{output from pairwise.comparisons
  function}

  \item{top.findings.qvals}{output from
  pairwise.comparisons function}

  \item{qth}{q-value threshold}

  \item{qth.star}{q-value threshold for stars}

  \item{mode}{barplot / heatmap Returns:}
}
\value{
  List with top findings from pairwise comparisons and
  their q-values
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

