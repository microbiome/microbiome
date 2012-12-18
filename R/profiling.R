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

#' run.profiling.script
#' 
#' Description: Profiling main script
#'
#' Arguments:
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name (HITChip: "Phyloarray"; MITChip: "Phyloarray_MIT"; PITChip old: "Phyloarray_PIT"; PITChip new: "pitchipdb")
#'   @param verbose verbose
#'
#' Returns:
#'   @return Profiling parameters. Also writes output to the user-specified directory.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

run.profiling.script <- function (dbuser, dbpwd, dbname, verbose = TRUE) {

  # Fetch and preprocess the data		     
  chipdata  <- preprocess.chipdata(dbuser, dbpwd, dbname)

  finaldata <- chipdata$data
  params    <- chipdata$params
  phylogeny.info <- chipdata$phylogeny.info

  ## Write preprocessed data in tab delimited file
  outd <- WriteChipData(finaldata, params$wdir, phylogeny.info, verbose = verbose)

  # Add oligo heatmap into output directory
  # Provide oligodata in the _original (non-log) domain_
  hc.params <- add.heatmap(finaldata[["oligo"]], output.dir = params$wdir, phylogeny.info = phylogeny.info)

  # Plot hierachical clustering trees into the output directory
  if (ncol(finaldata[["oligo"]])>2) { 

    method <- "complete"
    dat <- finaldata[["oligo"]]
	print(head(dat))
    png(paste(params$wdir, "hclust_oligo_euclidean_log10.png", sep = ""), width = 480, height = 480 * ncol(dat)/20)
    hc <- hclust(dist(t(log10(dat + 1))), method = method)
    plot(hc, hang = -1, main = "hclust/euclid/oligo/log10", xlab = "Samples")
    dev.off()

    png(paste(params$wdir, "hclust_oligo_euclidean_absolute.png", sep = ""), width = 480, height = 480 * ncol(dat)/20)
    hc <- hclust(dist(t(dat)), method = method)
    plot(hc, hang = -1, main = "hclust/euclid/oligo/raw", xlab = "Samples")
    dev.off()

    png(paste(params$wdir, "hclust_oligo_pearson_log10.png", sep = ""), width = 480, height = 480 * ncol(dat)/20)
    hc <- hclust(as.dist(1 - cor(log10(dat + 1), use = "complete.obs")), method = method)
    plot(hc, hang = -1, main = "hclust/pearson/oligo/log10", xlab = "Samples")
    dev.off()

    png(paste(params$wdir, "hclust_oligo_pearson_absolute.png", sep = ""), width = 480, height = 480 * ncol(dat)/20)
    hc <- hclust(as.dist(1 - cor(dat, use = "complete.obs")), method = method)
    plot(hc, hang = -1, main = "hclust/pearson/oligo/raw", xlab = "Samples")
    dev.off()

  }

  # Plot hclust trees on screen
  tmp <- htree.plot(finaldata[["oligo"]])

  # Write parameters into log file
  tmp <- WriteLog(chipdata$naHybs, params)
  params$logfilename <- tmp$log.file
  params$paramfilename <- tmp$parameter.file

  ## featurelevel data: fdat.orig, fdat.hybinfo, fdat.oligoinfo, 
  # save(finaldata, phylogeny.info, params, file = paste(params$wdir, "/sourcefiles.RData", sep = ""), compress = "xz")

  params

}


#' calculate.hclust
#' 
#' Description: Calculate hierarchical clustering for standard selections in profiling script
#'
#' Arguments:
#'   @param dat data matrix 
#'   @param method hierarchical clustering method (see ?hclust)
#'   @param metric clustering metrics (euclidean / correlation)
#'
#' Returns:
#'   @return hclust object for log10 and for absolute scale data
#'
#' @export
#' @examples # hc <- calculate.hclust(dat, "complete", "correlation") 
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

calculate.hclust <- function (dat, method = "complete", metric = "correlation") {

  if (metric == 'euclidean') {
    hc <- hclust(dist(t(dat)), method = method)
  } else if (metric == 'correlation') {
    hc <- hclust(as.dist(1 - cor(dat, use = "complete.obs")), method = method)
  } else {  
    stop("Provide proper metric for calculate.hclust!")
  }

  hc

}


#' Description: get probeset data matrix
#' 
#' Arguments:
#'   @param name name
#'   @param level taxonomic level
#'   @param phylogeny.info phylogeny.info
#'   @param oligo.matrix oligos vs. samples preprocessed data matrix; absolute scale
#'   @param log10 Logical. Log or no log?
#'
#' Returns:
#'   @return probeset data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.probeset <- function (name, level, phylogeny.info, oligo.matrix, log10 = TRUE) {

  # Pick probes for this entity
  probes <- retrieve.probesets(phylogeny.info, level, name)

  sets <- vector(length = length(probes), mode = "list")
  names(sets) <- names(probes)

  for (nam in names(probes)) {

    # Pick expression for particular probes (absolute scale)
    p <- intersect(probes[[nam]], rownames(oligo.matrix))
    dat <- NULL
    if (length(p) > 0) {
      dat <- oligo.matrix[p, ]
  
      dat <- matrix(dat, nrow = length(probes[[nam]]))
      rownames(dat) <- probes[[nam]]
      colnames(dat) <- colnames(oligo.matrix)

      # Logarithmize probeset?
      if ( log10 ) { dat <- log10(dat) } 
    } 
    sets[[nam]] <- dat

  }
  
  if (length(sets) == 1) {sets <- sets[[1]]}

  # Return
  sets

}

