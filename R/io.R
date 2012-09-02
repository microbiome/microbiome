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
#'   @param level phylogenetic level ("oligo" / "species" / "L1" / "L2" / "L0") or "phylogeny.info"
#'   @param method ("rpa" / "sum" / "ave" / "nmf")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @examples # params <- run.profiling.script(...); dat <- read.profiling.data(params$wdir, "L1", "rpa")
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling <- function(level, method = "sum", data.dir = NULL, log10 = TRUE, impute = TRUE){

  # level <- "oligo"; method = "sum"; data.dir = "test/"; log10 = TRUE

  svDialogsT <- require(svDialogs)
  if(!svDialogsT) { install.packages("svDialogs") }

  ##  Select file
  if (is.null(data.dir)) {

    f <- tk_choose.files(multi = F)

    # Recognize level and method from the file name 
    level <- NULL; method <- NULL
    if (!length(grep("oligo", f))==0) { level <- "oligo"}
    if (!length(grep("species", f)) == 0) { level <- "species"}
    if (!length(grep("L0", f)) == 0) { level <- "L0"}
    if (!length(grep("L1", f)) == 0) { level <- "L1"}
    if (!length(grep("L2", f)) == 0) { level <- "L2"}
    if (!length(grep("phylogeny", f)) == 0) { level <- "phylogeny.info"}
    if (!length(grep("rpa", f)) == 0) { level <- "rpa"}
    if (!length(grep("sum", f)) == 0) { level <- "sum"}
    if (!length(grep("ave", f)) == 0) { level <- "ave"}
    if (!length(grep("nmf", f)) == 0) { level <- "nmf"}

    tclServiceMode(FALSE)

  } else {
    if (level %in% c("L0", "L1", "L2", "species")) {
      f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")
    } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    } else if (level == "phylogeny.info") {
      f <- paste(data.dir, "/phylogeny.info.tab", sep = "")
    }
  }

  if (level %in% c("L0", "L1", "L2", "species")) {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1)

  } else if (level == "oligo") {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1)

  } else if (level == "phylogeny.info") {

    tab <- read.csv(f, header = TRUE, sep = "\t")
    tab <- polish.phylogeny.info(tab)

  }

  # Convert to numeric
  if (level %in% c("oligo", "species", "L0", "L1", "L2")) {
  
    rnams <- rownames(tab)
    cnams <- colnames(tab)
    tab <- apply(tab, 2, as.numeric)
    rownames(tab) <- rnams
    colnames(tab) <- cnams

  }

  if (impute && any(is.na(tab))) {
    warning(paste("The matrix has ", sum(is.na(tab)), " missing values - imputing."))
    tab <- 10^t(impute(t(log10(tab))))
  }

  if (log10 && (level %in% c("oligo", "species", "L0", "L1", "L2"))) {
    message("Logarithmizing the data")
    tab <- log10(tab)        
  }

  tab    

}
