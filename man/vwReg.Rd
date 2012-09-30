\name{vwReg}
\alias{vwReg}
\title{Description: Draw regression curve with smoothed error bars
based on the Visually-Weighted Regression by Solomon M. Hsiang; see
http://www.fight-entropy.com/2012/07/visually-weighted-regression.html
The R implementation is based on Felix Schönbrodt's code (under MIT/FreeBSD license) from http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/}
\usage{
  vwReg(formula, data, B = 1000, shade = TRUE,
    spag = FALSE, mweight = TRUE, show.lm = FALSE,
    show.CI = FALSE)

  vwReg(formula, data, B = 1000, shade = TRUE,
    spag = FALSE, mweight = TRUE, show.lm = FALSE,
    show.CI = FALSE)
}
\arguments{
  \item{formula}{formula}

  \item{data}{data}

  \item{title}{title}

  \item{B}{number bootstrapped smoothers}

  \item{shade}{plot the shaded confidence region?}

  \item{shade.alpha}{should the CI shading fade out at the
  edges? (by reducing alpha; 0 = no alpha decrease, 0.1 =
  medium alpha decrease, 0.5 = strong alpha decrease)}

  \item{spag}{plot spaghetti lines?}

  \item{mweight}{should the median smoother be visually
  weighted?}

  \item{show.lm}{should the linear regresison line be
  plotted?}

  \item{show.median}{show median smoother}

  \item{median.col}{median color}

  \item{show.CI}{should the 95\% CI limits be plotted?}

  \item{method}{the fitting function for the spaghettis;
  default: loess}

  \item{bw}{define a default b/w-palette (TRUE/FALSE)}

  \item{slices}{number of slices in x and y direction for
  the shaded region. Higher numbers make a smoother plot,
  but takes longer to draw. I wouldn'T go beyond 500}

  \item{palette}{provide a custom color palette for the
  watercolors}

  \item{ylim}{restrict range of the watercoloring}

  \item{quantize}{either "continuous", or "SD". In the
  latter case, we get three color regions for 1, 2, and 3
  SD (an idea of John Mashey)}

  \item{...}{further parameters passed to the fitting
  function, in the case of loess, for example, "span = .9",
  or "family = 'symmetric'"

  Returns:}

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

  ggplot2 object
}
\description{
  Arguments:

  Arguments:
}
\author{
  Based on the original version from Felix Schönbrodt.
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}

  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")

  See citation("microbiome")
}
\keyword{utilities}

