\name{InstallMarginal}
\alias{InstallMarginal}
\title{Install marginal dependencies on-demand from other functions.}
\usage{
  InstallMarginal(package, repo = "CRAN", ...)
}
\arguments{
  \item{package}{String name for the name to be installed}

  \item{repo}{CRAN or BioC}

  \item{...}{further arguments passed on to
  install.packages}
}
\value{
  Invisible NULL
}
\description{
  Function is always supposed to be called from a parent
  function and if the marginal depedency package is not
  installed, the function will report the name of the
  parent function requiring the package. Note that the
  whole call stack is not reported, only the immediate
  parent.
}
\note{
  meant for package internal use only
}
\author{
  Contact: \email{microbiome-admin@googlegroups.com}
}
\keyword{utilities}

