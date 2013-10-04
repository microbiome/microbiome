\name{Phylogeneticenrichments}
\alias{Phylogeneticenrichments}
\title{Phylogeneticenrichments}
\usage{
  Phylogeneticenrichments(x, phylogeny.info,
    origlevel = colnames(phylogeny.info)[3],
    maplevel = "L1", onlyEnriched = T, p.th = 0.05)
}
\arguments{
  \item{x}{x}

  \item{phylogeny.info}{phylogeny.info}

  \item{origlevel}{origlevel}

  \item{maplevel}{maplevel}

  \item{onlyEnriched}{onlyEnriched}

  \item{p.th}{p-value threshold

  Returns:}
}
\value{
  enrichments
}
\description{
  Description: Computes the dependence of the higher
  phylogenetic level groups of the given list x, which can
  be a vector/data frame of either oligos, species, or
  level2 groups. Computes only the enrichments if
  onlyEnrich=T (neglects the unexpected disappearences of
  the groups in the given list). Outputs only the
  enrichments below the given p-value threshold p.th.
  Requires phylogeny.info=the first 4 columns from
  phylogenyprofile, origlevel=the level from which the list
  is, and maplevel= the level to which the mapping should
  be made and the enrichments of which computed. Outputs a
  four component list: p-values from Fisher's tests, input
  tables for tests, the actual tests outputs, and the full
  phylogenetic mapping information for the given list.
  Changes: Version 1 computed two-tailed p-values also for
  single occurrences on origlevel. v2: two-tailed p-values,
  but not for single occurrences v3: one-tailed p-values
  (enrichments), if wanted, outputting only enrichments
  under the given p-value
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

