\name{PlotCrosshyb}
\alias{PlotCrosshyb}
\title{Description: Cross-hybridization between multimodal taxa as percentages of shared probes.
The number indicates how many percent of oligos for the row taxon are also hybridizing
the corresponding column taxon.}
\usage{
  PlotCrosshyb(tax.level = "L1", chip = "HITChip",
    selected.taxa = NULL, show.plot = TRUE,
    order.rows = TRUE, order.cols = TRUE,
    keep.empty = FALSE, rounding = 1,
    phylogeny.info = NULL)
}
\arguments{
  \item{tax.level}{Taxonomic level to investigate}

  \item{chip}{Chip type (e.g. "HITChip")}

  \item{selected.taxa}{Restrict cross-hyb analysis to the
  selected groups.}

  \item{show.plot}{Produce the plot}

  \item{order.rows}{Order table rows}

  \item{order.cols}{Order table columns}

  \item{keep.empty}{Keep taxa that do not show any
  cross-hybridization}

  \item{rounding}{Rounding of the cell contents}

  \item{phylogeny.info}{phylogeny.info

  Returns:}
}
\value{
  A list containing cross-hybridization table and plot
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

