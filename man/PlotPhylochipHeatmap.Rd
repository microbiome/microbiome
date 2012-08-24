\name{PlotPhylochipHeatmap}
\alias{PlotPhylochipHeatmap}
\title{PlotPhylochipHeatmap}
\usage{
  PlotPhylochipHeatmap(data, phylogeny,
    metric = "euclidian", tax.level = "level 2",
    include.tree = FALSE, palette = "black/yellow/white",
    fontsize = 12, figureratio = 12,
    hclust.method = "ward")
}
\arguments{
  \item{data}{oligoprofile data in original (non-log)
  domain}

  \item{phylogeny}{oligo-phylotype mappings}

  \item{metric}{clustering metric}

  \item{tax.level}{taxonomic level to show}

  \item{include.tree}{include.tree}

  \item{palette}{color palette}

  \item{fontsize}{font size}

  \item{figureratio}{figure ratio}

  \item{hclust.method}{hierarchical clustering method

  Returns:}
}
\value{
  plots the profile
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

