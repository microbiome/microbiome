\name{check.wilcoxon}
\alias{check.wilcoxon}
\title{Description: Calculate Wilcoxon test (unpaired; BH correction) for the specified sample groups.
Either provide the input data as matrix, file path, or select the file through GUI.}
\usage{
  check.wilcoxon(dat = NULL, fnam = NULL,
    p.adjust.method = "BH", sort = FALSE)
}
\arguments{
  \item{dat}{data matrix (features x samples)}

  \item{file}{data file (if data matrix not provided)}

  \item{p.adjust.method}{p-value correction method for
  p.adjust function (default "BH")}

  \item{sort}{sort the results Returns:}
}
\value{
  (Corrected) p-values for two-group comparison.
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

