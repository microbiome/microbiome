\name{check.wilcoxon}
\alias{check.wilcoxon}
\title{Description: Calculate Wilcoxon test (unpaired; BH correction) for the specified sample groups.
Either provide the input data as matrix, file path, or select the file through GUI.}
\usage{
  check.wilcoxon(dat = NULL, fnam = NULL, G1, G2,
    p.adjust.method = "BH", sort = FALSE)
}
\arguments{
  \item{dat}{data matrix (features x samples)}

  \item{fnam}{data file (if data matrix not provided)}

  \item{G1}{Sample group 1 (for comparison)}

  \item{G2}{Sample group 2 (for comparison)}

  \item{p.adjust.method}{p-value correction method for
  p.adjust function (default "BH"). If NULL, no correction
  will be performed.}

  \item{sort}{sort the results

  Returns:}
}
\value{
  (Corrected) p-values for two-group comparison.
}
\description{
  Arguments:
}
\examples{
data(peerj32); pval <- check.wilcoxon(t(peerj32$microbes), G1 = 1:22, G2 = 23:44)
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

