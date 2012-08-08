\name{sampleOutliers}
\alias{sampleOutliers}
\title{Outlier detection based on replicate measurements of a sample}
\usage{
  sampleOutliers(con, sampleID, significance = 0.001)
}
\arguments{
  \item{con}{MySQL connection}

  \item{sampleID}{sampleID}

  \item{significance}{significance}
}
\value{
  A sampleOutliers object which is a list containing a
  'spots' dataframe, a 'dbcon' connection to the phyloarray
  database from which the sample data were taken and a
  'sampleID' length 1 character vector containing the
  sampleID.
}
\description{
  This function takes all arrays on which a sample was
  analysed and calculates, based on replicates of oligo's
  within and between arrays which spots on which arrays are
  likely to be outliers.  The isOutlier flag of the spot is
  set in the database if the spot is likely to be an
  outlier.  The routine needs normalized data. If these are
  unavailable for an array on which the sample was
  hybridised then the calculation is skipped and a warning
  is issued in the R console.
}
\details{
  Note: sampleOutliers only considers probes in the class
  'ssRNA', and does not take into account hybridisations of
  which the 'isDiscarded' flag is set, or
  featureExtractions of which the 'noReproCheck' is set.
}
\examples{
# TBA
}
\author{
  Douwe Molenaar. Maintainer: Leo Lahti
  \email{leo.lahti@iki.fi}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

