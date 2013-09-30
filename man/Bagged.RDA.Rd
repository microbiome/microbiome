\name{Bagged.RDA}
\alias{Bagged.RDA}
\title{Description: BAGGED RDA Bootstrap solutions that follows the
Jack-knife estimation of PLS by Martens and Martens, 2000.  Solves
rotational invariance of latent space by orthogonal procrustes
rotations}
\usage{
  Bagged.RDA(X, Y, boot = 1000)
}
\arguments{
  \item{X}{a matrix, samples on columns, variables
  (bacteria) on rows.}

  \item{Y}{vector with names(Y)=rownames(X), for example}

  \item{boot}{Number of bootstrap iterations

  Returns:}
}
\value{
  List with items: loadings: bagged loadings scores: bagged
  scores significance: significances of X variables etc.
  TBA.
}
\description{
  Arguments:
}
\author{
  Contact: Jarkko Salojarvi
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

