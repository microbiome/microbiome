\name{get.hclust.groups}
\alias{get.hclust.groups}
\title{Description: Get hierarchical clustering groups, including the nested ones if recursive = TRUE}
\usage{
  get.hclust.groups(dat, corr.th, recursive = FALSE,
    min.size = 2)
}
\arguments{
  \item{dat}{data matrix}

  \item{corr.th}{correlation threshold for hclust}

  \item{recursive}{get also nested groups TRUE/FALSE}

  \item{min.size}{minimum cluster size}
}
\value{
  groups
}
\description{
  Description: Get hierarchical clustering groups,
  including the nested ones if recursive = TRUE
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

