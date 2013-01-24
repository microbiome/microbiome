\name{pick.training.samples}
\alias{pick.training.samples}
\title{pick.training.samples}
\usage{
  pick.training.samples(sample.info,
    training.fraction = 0.8, rseed = 1463)
}
\arguments{
  \item{sample.info}{Sample metadata data frame. Include
  the fields 'sampleID' and 'projectName', specifying the
  sample identifier and project name for each sample.}

  \item{training.fraction}{Specify the fraction of the data
  to go into training set}

  \item{rseed}{random seed Returns:}
}
\value{
  List with two elements: training and test, listing the
  corresponding sample sets.
}
\description{
  Description: Split data set into training and test
  samples
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

