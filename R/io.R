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


#' read.profiling
#' 
#' Description: read run.profiling.script output into R
#'
#' Arguments:
#'   @param level phylogenetic level ("oligo" / "species" / "L1" / "L2" / "L0") or "phylogeny"
#'   @param method ("rpa" / "sum" / "ave" / "nmf")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @examples # params <- run.profiling.script(...); dat <- read.profiling.data(params$wdir, "level1", "rpa")
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling <- function(level, method = "sum", data.dir = NULL, log10 = TRUE){

  svDialogsT <- require(svDialogs)
  if(!svDialogsT) { install.packages("svDialogs") }

  ##  Select file
  if (is.null(data.dir)) {
    f <- choose.files(multi=F)
  } else {
    if (level %in% c("level0", "level1", "level2", "species")) {
      f <- paste(data.dir, "/", level, "_log10_", method, ".tab", sep = "")
    } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    } else if (level == "phylogeny") {
      f <- paste(data.dir, "/phylogenyinfo.tab", sep = "")
    }
  }

  # Recognize level and method from the file name
  level <- NULL; method <- NULL
  if (grep("oligo", f)) { level <- "oligo"}
  if (grep("species", f)) { level <- "species"}
  if (grep("level0", f)) { level <- "level0"}
  if (grep("level1", f)) { level <- "level1"}
  if (grep("level2", f)) { level <- "level2"}
  if (grep("phylogeny", f)) { level <- "phylogeny"}
  if (grep("rpa", f)) { method <- "rpa"}
  if (grep("sum", f)) { method <- "sum"}
  if (grep("ave", f)) { method <- "ave"}

  if (level %in% c("level0", "level1", "level2", "species")) {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1)

  } else if (level == "oligo") {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1)

  } else if (level == "phylogeny") {

    tab <- read.csv(f, header = TRUE, sep = "\t")
    colnames(tab)[[which(colnames(tab) == "level.1")]] <- "level 1"
    colnames(tab)[[which(colnames(tab) == "level.2")]] <- "level 2"

  }

  if (log10 && (level %in% c("oligo", "species", "level0", "level1", "level2"))) {
    message("Logarithmizing the data")
    tab <- log10(tab)        
  }

  tab    

}
