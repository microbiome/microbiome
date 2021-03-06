% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_bmi.R
\name{group_bmi}
\alias{group_bmi}
\title{Body-Mass Index (BMI) Classes}
\usage{
group_bmi(
  x,
  breaks = "standard",
  n = 10,
  labels = NULL,
  include.lowest = TRUE,
  right = FALSE,
  dig.lab = 3,
  ordered_result = FALSE
)
}
\arguments{
\item{x}{Numeric vector (BMI)}

\item{breaks}{Class break points. Either a vector of breakpoints, or one of
the predefined options ("standard", "standard_truncated", "even").}

\item{n}{Number of groups for the breaks = "even" option.}

\item{labels}{labels for the levels of the resulting category.  By default,
    labels are constructed using \code{"(a,b]"} interval notation.  If
    \code{labels = FALSE}, simple integer codes are returned instead of
    a factor.}

\item{include.lowest}{logical, indicating if an \sQuote{x[i]} equal to
    the lowest (or highest, for \code{right = FALSE}) \sQuote{breaks}
    value should be included.}

\item{right}{logical, indicating if the intervals should be closed on
    the right (and open on the left) or vice versa.}

\item{dig.lab}{integer which is used when labels are not given.  It
    determines the number of digits used in formatting the break numbers.}

\item{ordered_result}{logical: should the result be an ordered factor?}
}
\value{
Factor of BMI groups.
}
\description{
Cut BMI information to standard discrete factors.
}
\details{
Regarding the breaks arguments, the "even" option aims to
cut the samples in groups with approximately the same size (by
quantiles). The "standard" option corresponds to standard obesity
categories defined by the cutoffs <18.5 (underweight); <25 (lean);
<30 (obese); <35 (severe obese); <40 (morbid obese); <45 (super
obese). The standard_truncated combines the severe, morbid and super
obese into a single group.
}
\examples{
bmi.numeric <- range(rnorm(100, mean = 25, sd = 3))
bmi.factor <- group_bmi(bmi.numeric)
}
\references{
See citation('microbiome')
}
\seealso{
base::cut
}
\author{
Contact: Leo Lahti \email{microbiome-admin@googlegroups.com}
}
\keyword{utilities}
