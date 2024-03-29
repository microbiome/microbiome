% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rare_members.R
\name{rare_members}
\alias{rare_members}
\title{Rare Taxa}
\usage{
rare_members(x, detection = 1/100, prevalence = 50/100, include.lowest = FALSE)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object}

\item{detection}{Detection threshold for absence/presence
(strictly greater by default).}

\item{prevalence}{Prevalence threshold (in [0, 1]). The
required prevalence is strictly greater by default. To include the
limit, set include.lowest to TRUE.}

\item{include.lowest}{Include the lower boundary of the detection and
prevalence cutoffs. FALSE by default.}
}
\value{
Vector of rare taxa
}
\description{
Determine members of the rare microbiota with given abundance
and prevalence threshold.
}
\details{
For phyloseq object, lists taxa that are less prevalent than the
given prevalence threshold. Optionally, never exceeds the given abundance
threshold (by default, all abundaneces accepted). For matrix, lists
columns that satisfy these criteria.
}
\examples{
data(dietswap)
# Detection threshold: the taxa never exceed the given detection threshold
# Prevalence threshold 20 percent (strictly greater by default)
a <- rare_members(dietswap, detection=100/100, prevalence=20/100)
}
\references{
To cite the microbiome R package, see citation('microbiome')
}
\seealso{
core_members
}
\author{
Leo Lahti \email{microbiome-admin@googlegroups.com}
}
\keyword{utilities}
