# Copyright (C) 2011-2013 Leo Lahti and Jarkko Salojarvi 
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
#'   @param output.dir output directory name
#'
#' Returns:
#'   @return output file name
#'
#' @export
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GenerateSimulatedData <- function (output.dir) {

  data.directory <- system.file("extdata/", package = "microbiome")

  phylogeny.info <- read.profiling(level = "phylogeny.info", data.dir = data.directory)[, 1:6]

  oligo.matrix.nolog.simulated <- read.profiling(level = "oligo", data.dir = data.directory, log10 = FALSE)
  N <- ncol(oligo.matrix.nolog.simulated)
  colnames(oligo.matrix.nolog.simulated) <- paste("Sample.", 1:N, sep = "")

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.matrix.nolog.simulated
  levels <- c("species", "L2", "L1")
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa", "nmf")) {

        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), method = method, level = level, rm.phylotypes = phylotype.rm.list("HITChip"))
      			       	          
        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  # Write summary matrices into the output directory
  outd <- HITChipDB::WriteChipData(finaldata, output.dir, phylogeny.info)

  set.seed(344)
  metadata.simulated <- data.frame(list(
              sampleID = I(colnames(oligo.matrix.nolog.simulated)),
	      time = rep(1:4, 5),
              age = runif(N, 0, 100),
              bmi = runif(N, 20, 40),
	      subjectID = I(paste("subjectID", rep(1:4, 5), sep = "-")),
	      group = I(sample(paste("group", rep(1:4, 5), sep = "-"))),
              gender = I(sample(c("M", "F"), N, replace = TRUE)),
              diet = I(sample(c("Apricots", "Beverages", "Carrots"), N, replace = TRUE))))
  
  write.table(metadata.simulated, file = paste(output.dir, "/metadata.tab", sep = ""), quote = FALSE, sep = "\t")

  output.dir

}


