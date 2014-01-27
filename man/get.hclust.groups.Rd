\name{get.hclust.groups}
\alias{get.hclust.groups}
\title{Description: Get hierarchical clustering groups, including the nested ones if recursive = TRUE}
\usage{
get.hclust.groups(dat, corr.th, recursive = FALSE, min.size = 2,
  metric = "pearson")
}
\arguments{
  \item{dat}{data matrix}

  \item{corr.th}{correlation threshold for hclust}

  \item{recursive}{get also nested groups TRUE/FALSE}

  \item{min.size}{minimum cluster size}

  \item{metric}{similarity measure (spearman / pearson)}
}
\value{
groups
}
\description{
Description: Get hierarchical clustering groups, including
the nested ones if recursive = TRUE
}
\examples{
data(peerj32); cl <- get.hclust.groups(peerj32$microbes, corr.th = 0.8)
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
See citation("microbiome")
}
\keyword{methods}

