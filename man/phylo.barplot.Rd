\name{phylo.barplot}
\alias{phylo.barplot}
\title{phylo.barplot}
\usage{
  phylo.barplot(x, color.level = "L1",
    phylogeny.info = NULL, title = NULL, plot = TRUE)
}
\arguments{
  \item{x}{Data vector across taxa (each element should be
  named by taxon)}

  \item{color.level}{Higher-order phylogenetic level to
  indicate by colors}

  \item{phylogeny.info}{oligo-phylotype mappings}

  \item{title}{title}

  \item{plot}{draw plot TRUE/FALSE

  Returns:}
}
\value{
  ggplot2 object
}
\description{
  Description: Barplot for *ITChip sample (across taxa)
  with higher-level taxonomic groups indicated by colours.
}
\details{
  Arguments:
}
\examples{
# p <- phylo.barplot(signal, color.level = "L1", phylogeny.info = NULL)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

