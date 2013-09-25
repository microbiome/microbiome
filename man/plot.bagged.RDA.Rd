\name{plot.bagged.RDA}
\alias{plot.bagged.RDA}
\title{Description: Bagged RDA visualization}
\usage{
  plot.bagged.RDA(Bag.res, Y,
    which.bac = 1:nrow(Bag.res$loadings), ptype = "spider",
    comp = 1:2, cex.bac = 0.5, plot.names = T,
    group.cols = as.numeric(unique(Y)), ...)
}
\arguments{
  \item{Bag.res}{Output from Bagged.RDA.Feature.Selection
  function}

  \item{Y}{vector with names(Y)=rownames(X), for example}

  \item{which.bac}{TBA}

  \item{ptype}{TBA}

  \item{comp}{TBA}

  \item{cex.bac}{TBA}

  \item{plot.names}{Plot names}

  \item{group.controls}{TBA

  Returns:}
}
\value{
  TBA
}
\description{
  Arguments:
}
\examples{
#
  #library(microbiome)
  #data(peerj32)
  #x <- t(peerj32$microbes)
  #y <- factor(peerj32$meta$time); names(y) <- rownames(peerj32$meta)
  #Bag.res <- Bagged.RDA.Feature.Selection(x, y, sig.thresh=0.05, nboot=100)
  #plot.bagged.RDA(Bag.res, y)
}
\author{
  Contact: Jarkko Salojarvi
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

