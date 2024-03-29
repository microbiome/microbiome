% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inequality.R
\name{inequality}
\alias{inequality}
\title{Gini Index}
\usage{
inequality(x)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object}
}
\value{
A vector of Gini indices
}
\description{
Calculate Gini indices for a phyloseq object.
}
\details{
Gini index is a common measure for relative inequality in
economical income, but can also be used as a community diversity
measure. Gini index is between [0,1], and increasing gini index implies
increasing inequality.
}
\examples{
data(dietswap)
d <- inequality(dietswap)
}
\references{
Relative Distribution Methods in the Social Sciences. Mark S.
Handcock and Martina Morris, Springer-Verlag, Inc., New York,
1999. ISBN 0387987789.
}
\seealso{
diversity, reldist::gini (inspired by that implementation but
independently written here to avoid external dependencies)
}
\author{
Contact: Leo Lahti \email{microbiome-admin@googlegroups.com}
}
\keyword{utilities}
