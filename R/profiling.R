# Copyright (C) 2006-2012 Leo Lahti, and Jarkko Salojarvi, Janne
# Nikkila, and Douwe Molenaar. All rights reserved.
# Contact: <leo.lahti@iki.fi>

# This file is a part of the microbiome R package
#
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
#'   @param dbname MySQL database name
#'   @param verbose verbose
#'
#' Returns:
#'   @return Profiling parameters. Also writes output to the user-specified directory.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

run.profiling.script <- function (dbuser, dbpwd, dbname, verbose = TRUE) {

  # Fetch and preprocess the data		     
  chipdata  <- preprocess.chipdata(dbuser, dbpwd, dbname)
  finaldata <- chipdata$data
  params    <- chipdata$params
  phylogeny <- chipdata$phylogeny

  ## Write preprocessed data in tab delimited file
  outd <- WriteChipData(finaldata, params$wdir, phylogeny, verbose = verbose)

  # Make basic plots
  # feed in here oligodata in _original (non-log) domain_
  plot.params <- add.hclust.plots(finaldata[["oligo"]], data.dir = params$wdir) 

  # Write log file
  tmp <- WriteLog(chipdata$naHybs, params)
  params$logfilename <- tmp$log.file
  params$paramfilename <- tmp$parameter.file

  ## featurelevel data: fdat.orig, fdat.hybinfo, fdat.oligoinfo, 
  # save(finaldata, phylogeny, params, file = paste(params$wdir, "/sourcefiles.RData", sep = ""), compress = "xz")

  params

}


#' calculate.hclust
#' 
#' Description: Calculate hierarchical clustering for standard selections in profiling script
#'
#' Arguments:
#'   @param dat data matrix for clustering in log10 domain
#'   @param hclust.method clustering method
#'   @param metric clustering metrics
#'
#' Returns:
#'   @return hclust object for log10 and for absolute scale data
#'
#' @export
#' @examples # TBA
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

calculate.hclust <- function (dat, hclust.method = "ward", metric = "correlation") {

  if (metric == 'euclidian') {

    hclog <- hclust(dist(t(dat)),          method = hclust.method)
    hcraw <- hclust(dist(t(10^(dat) - 1)), method = hclust.method)

  } else if (metric == 'correlation') {

    hclog <- hclust(as.dist(1 - cor(dat, use = "complete.obs")),          method = hclust.method)
    hcraw <- hclust(as.dist(1 - cor(10^(dat) - 1, use = "complete.obs")), method = hclust.method)

  } else {
  
    stop("Provide proper metric for calculate.hclust!")

  }

  list(log10 = hclog, raw = hcraw)

}


#' Description: get probeset data matrix
#' 
#' Arguments:
#'   @param name name
#'   @param level level
#'   @param phylo phylo
#'   @param data data
#'   @param log10 Logical. Log or no log?
#'
#' Returns:
#'   @return probeset data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.probeset <- function (name, level, phylo, data, log10 = TRUE) {

  # Pick probes for this entity
  probes <- list.probes(name, level, phylo)
  
  # Pick expression for particular probes
  dat <- data$probes[probes,]

  # Log?
  if ( log10 ) { dat <- log10(dat) } 
  
  # Return
  dat
}




