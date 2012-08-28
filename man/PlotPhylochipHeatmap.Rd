\name{PlotPhylochipHeatmap}
\alias{PlotPhylochipHeatmap}
\title{PlotPhylochipHeatmap}
\usage{
  PlotPhylochipHeatmap(data, oligomap, metric = "pearson",
    level = "L2", include.tree = TRUE,
    palette = "white/blue", fontsize = 12,
    figureratio = 12, hclust.method = "complete")
}
\arguments{
  \item{data}{oligoprofile data in original (non-log)
  domain}

  \item{oligomap}{oligo-phylotype mappings}

  \item{metric}{clustering metric}

  \item{level}{taxonomic level to show (L0 / L1 / L2 /
  species)}

  \item{include.tree}{include.tree}

  \item{palette}{color palette}

  \item{fontsize}{font size}

  \item{figureratio}{figure ratio}

  \item{hclust.method}{hierarchical clustering method

  Returns:}
}
\value{
  parameters
}
\description{
  Description: Plots heatmap of the oligo profiles together
  with phylotype grouping and sample clustering
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

