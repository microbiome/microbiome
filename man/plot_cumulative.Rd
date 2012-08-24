\name{plot_cumulative}
\alias{plot_cumulative}
\title{plot_cunulative}
\usage{
  plot_cumulative(d.sub, writedir, fname, i.set = NULL,
    type = "cumulative", ylim = NULL, phylogeny)
}
\arguments{
  \item{d.sub}{d.sub}

  \item{writedir}{output directory}

  \item{fname}{output file name}

  \item{i.set}{i.set}

  \item{type}{plot type}

  \item{ylim}{y axis limits}

  \item{phylogeny}{phylogeny matrix

  Returns:}
}
\value{
  TBA
}
\description{
  Description: plot_cumulative
}
\details{
  Arguments:
}
\examples{
# data("toydata-hitchip"); bs <- bootstrap.microbes(phylum.matrix); plot_cumulative(bs, writedir = "./", "tmp-", phylogeny = phylogeny)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

