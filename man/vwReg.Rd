\name{vwReg}
\alias{vwReg}
\title{Description: Draw regression curve with smoothed error bars
based on the Visually-Weighted Regression by Solomon M. Hsiang; see
http://www.fight-entropy.com/2012/07/visually-weighted-regression.html
The R implementation is based on Felix Schönbrodt's code from
http://www.nicebread.de/visually-weighted-regression-in-r-a-la-solomon-hsiang/}
\usage{
  vwReg(formula, data, B = 1000, shade = TRUE,
    spag = FALSE, mweight = TRUE, show.lm = FALSE,
    show.CI = FALSE)
}
\arguments{
  \item{formula}{formula}

  \item{data}{data}

  \item{B}{number bootstrapped smoothers}

  \item{shade:}{plot the shaded confidence region?}

  \item{spag:}{plot spaghetti lines?}

  \item{mweight:}{should the median smoother be visually
  weighted?}

  \item{show.lm:}{should the linear regresison line be
  plotted?}

  \item{show.CI:}{should the 95% CI limits be plotted?

  Returns:}
}
\value{
  ggplot2 object
}
\description{
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

