\name{correlation.heatmap}
\alias{correlation.heatmap}
\title{correlation.heatmap}
\usage{
  correlation.heatmap(df, Xvar, Yvar, fill, star = "p.adj",
    p.adjusted.threshold = 0.05, correlation.threshold = 0,
    step = 0.2,
    colours = c("darkblue", "blue", "white", "red", "darkred"),
    limits = c(-1, 1), legend.text = fill,
    order.rows = TRUE, order.cols = TRUE, text.size = 10,
    filter.significant = TRUE)
}
\arguments{
  \item{df}{Data frame. Each row corresponds to a pair of
  correlated variables. The columns give variable names,
  correlations and significance estimates.}

  \item{Xvar}{X axis variable column name. For instance
  "X".}

  \item{Yvar}{Y axis variable column name. For instance
  "Y".}

  \item{fill}{Column to be used for heatmap coloring. For
  instance "correlation".}

  \item{star}{Column to be used for cell highlighting. For
  instance "p.adj".}

  \item{p.adjusted.threshold}{Significance threshold for
  the stars.}

  \item{correlation.threshold}{Include only elements that
  have absolute correlation higher than this value}

  \item{step}{color interval}

  \item{colours}{heatmap colours}

  \item{limits}{colour scale limits}

  \item{legend.text}{legend text}

  \item{order.rows}{Order rows to enhance visualization
  interpretability}

  \item{order.cols}{Order columns to enhance visualization
  interpretability}

  \item{text.size}{Adjust text size}

  \item{filter.significant}{Keep only the elements with at
  least one significant entry

  Returns:}
}
\value{
  ggplot2 object
}
\description{
  Description: Visualizes n x m correlation table as
  heatmap. See examples for details.
}
\details{
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

