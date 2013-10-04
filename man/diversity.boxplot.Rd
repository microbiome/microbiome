\name{diversity.boxplot}
\alias{diversity.boxplot}
\title{diversity.boxplot}
\usage{
  diversity.boxplot(dat, sample.groups,
    diversity.index = "shannon", det.th = NULL,
    title = NULL, col.list = NULL, ylim = NULL)
}
\arguments{
  \item{dat}{data matrix in the original (non-log) scale:
  features x samples.}

  \item{sample.groups}{specify the distinct sample groups}

  \item{diversity.index}{Options: "shannon" / "invsimpson"
  / "richness" / "evenness"}

  \item{det.th}{Detection threshold for richness and
  evenness calculations}

  \item{title}{figure title}

  \item{col.list}{Optional: colors for the sample groups}

  \item{ylim}{y-axis limits

  Returns:}
}
\value{
  Sample group list corresponding to the boxplot groups.
}
\description{
  Description: diversity.boxplot
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

