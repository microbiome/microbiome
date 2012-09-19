# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Generate toydata for the package
#'
#' Arguments:
#'   @param oligo.matrix.nolog.simulated
#'   @param phylogeny.info phylogeny.info
#'   @param output.file output file name
#'
#' Returns:
#'   @return output file name
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GenerateSimulatedData <- function (oligo.matrix.nolog.simulated, phylogeny.info, output.file = "toydata-hitchip.rda") {

  set.seed(344)
  N <- ncol(oligo.matrix.nolog.simulated)
  #s <- sample(ncol(oligo.matrix.nolog.simulated), N)
  s <- colnames(oligo.matrix.nolog.simulated)
  colnames(oligo.matrix.nolog.simulated) <- paste("Sample-", 1:N, sep = "")

  #species.matrix.log10.simulated <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), "sum", "species", verbose = TRUE, phylotype.rm.list("HITChip"))

  genus.matrix.log10.simulated <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), "sum", "L2", verbose = TRUE, phylotype.rm.list("HITChip"))

  #phylum.matrix.log10.simulated <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), "sum", "L1", verbose = TRUE, phylotype.rm.list("HITChip"))

  metadata.simulated <- data.frame(list(
              sampleID = colnames(oligo.matrix.nolog.simulated),
	      time = rep(1:4, 5),
              age = runif(N, 0, 100),
              bmi = runif(N, 20, 40),
	      subjectID = paste("subjectID", rep(1:4, 5), sep = "-"),
	      group = sample(paste("group", rep(1:4, 5))),
              gender = sample(c("M", "F"), N, replace = TRUE),
              diet = sample(c("Apricots", "Beverages", "Carrots"), N, replace = TRUE)))


  message(paste("Saving the output in", output.file))
  #save(metadata.simulated, oligo.matrix.nolog.simulated, species.matrix.log10.simulated, genus.matrix.log10.simulated, phylum.matrix.log10.simulated, file = output.file, compress = "xz")

  save(metadata.simulated, oligo.matrix.nolog.simulated, genus.matrix.log10.simulated, file = output.file, compress = "xz")

  output.file

}


#' toydata-hitchip
#'
#' The toydata-hitchip contains the following simulated data matrices: 
#'
#' metadata.simulated sample metadata
#' oligo.matrix.nolog.simulated: oligo-level data matrix (oligos vs. samples)
#' genus.matrix.log10.simulated: genus (L2) level data matrix (phylotypes vs. samples)
#'
#' Generate with GenerateSimulatedData function
#'
#' @name toydata-hitchip
#' @docType data
#' @author Leo Lahti \email{microbiome-devel@@googlegroups.com}
#' @references \url{microbiome.github.com}
#' @keywords data
NULL

#' phylogeny
#'
#' The phylogeny contains the oligo-species-L2-L1 mappings for HITChip
#'
#' @name phylogeny
#' @docType data
#' @author Leo Lahti \email{microbiome-devel@@googlegroups.com}
#' @references \url{microbiome.github.com}
#' @keywords data
NULL
