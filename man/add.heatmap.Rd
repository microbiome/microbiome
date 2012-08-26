\name{add.heatmap}
\alias{add.heatmap}
\title{add.heatmap
Description: Add oligprofile heatmap into output directory}
\usage{
  add.heatmap(dat, output.dir, output.file = NULL,
    phylogeny, ppcm = 150, hclust.method = "complete",
    palette = "white/blue", level = "level 2",
    metric = "pearson", figureratio = 12, fontsize = 12,
    include.tree = TRUE)
}
\arguments{
  \item{dat}{oligoprofile data in original (non-log)
  domain}

  \item{output.dir}{output data directory}

  \item{output.dir}{output file name}

  \item{phylogeny}{oligo-phylotype mappings}

  \item{ppcm}{figure size}

  \item{hclust.method}{hierarchical clustering method}

  \item{palette}{color palette}

  \item{level}{taxonomic level to show}

  \item{metric}{clustering metric}

  \item{figureratio}{figure ratio}

  \item{fontsize}{font size}

  \item{include.tree}{include.tree

  Returns:}
}
\value{
  Plotting parameters
}
\description{
  Arguments:
}
\examples{
# dat <- read.profiling(params$wdir, "species", "rpa"); hc <- add.heatmap(dat)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

