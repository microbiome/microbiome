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


#' read.profiling
#' 
#' Description: read run.profiling.script output into R
#'
#' Arguments:
#'   @param level phylogenetic level ("oligo" / "species" / "L1" / "L2" / "L0") or "phylogeny.full", "phylogeny.filtered"
#'   @param method ("frpa" / "rpa" / "sum" / "ave" / "nmf")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'   @param impute impute missing oligo signals
#' 
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @examples # params <- run.profiling.script(...); dat <- read.profiling(params$wdir, "L1", "rpa")
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling <- function(level = NULL, method = "frpa", data.dir, log10 = TRUE, impute = TRUE){

  # level <- "oligo"; method = "sum"; data.dir = "test/"; log10 = TRUE
  if (level %in% c("L0", "L1", "L2", "species")) {
      f <- paste(data.dir, "/", level, "-", method, ".tab", sep = "")
    } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile.tab", sep = "")
    } else if (level == "phylogeny.full") {
      f <- paste(data.dir, "/phylogeny.full.tab", sep = "")
    } else if (level %in% c("phylogeny.filtered")) {
      f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")
    } else if (level %in% c("phylogeny.info", "phylogeny.full", "phylogeny")) {
      f <- paste(data.dir, "/phylogeny.full.tab", sep = "")
    }
  
  message(paste("Reading", f))

  if (level %in% c("L0", "L1", "L2", "species")) {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)

  } else if (level == "oligo") {

    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)

  } else if (length(grep("phylogeny", level)) > 0) {

    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
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

get.file.method <- function (f) {

  method <- NULL

    if (f == "rpa")  { method <- "rpa"}
    if (f == "frpa") { method <- "frpa"}
    if (!length(grep("sum", f)) == 0) { method <- "sum"}
    if (!length(grep("ave", f)) == 0) { method <- "ave"}
    if (!length(grep("nmf", f)) == 0) { method <- "nmf"}

  method

}
get.file.level <- function (f) {

    level <- NULL

    if (!length(grep("oligo", f))==0) { level <- "oligo"}
    if (!length(grep("species", f)) == 0) { level <- "species"}
    if (!length(grep("L0", f)) == 0) { level <- "L0"}
    if (!length(grep("L1", f)) == 0) { level <- "L1"}
    if (!length(grep("L2", f)) == 0) { level <- "L2"}
    if (!length(grep("phylogeny", f)) == 0) { level <- "phylogeny.info"}

  level

}

#' read.profiling.010
#' 
#' Description: read profiling script output from profiling script v. 010 into R 
#'
#' Arguments:
#'   @param level phylogenetic level ("oligo" / "species" / "L1" / "L2" / "L0") or "phylogeny.info"
#'   @param method ("rpa" / "sum" / "ave" / "nmf")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'   @param impute impute missing oligo signals
#'
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @examples # params <- read.profiling.010(...); 
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling.010 <- function(level = NULL, method = "rpa", data.dir = NULL, log10 = TRUE, impute = TRUE){

  if (level == "L0") { level <- "level0"}
  if (level == "L1") { level <- "level1"}
  if (level == "L2") { level <- "level2"}
  if (method == "rpa") {method <- "RPA"}
  if (method == "sum") {method <- "Sum"}
  if (method == "ave") {method <- "log10Ave"}
  if (method == "nmf") {method <- "NMF"}

  InstallMarginal("svDialogs")

  ##  Select file
  if (is.null(data.dir)) {

    f <- svDialogs::tk_choose.files(multi = F)

    # Recognize level and method from the file name 
    level <- NULL; method <- NULL
    if (!length(grep("oligo", f)) == 0) { level <- "oligo"}
    if (!length(grep("species", f)) == 0) { level <- "species"}
    if (!length(grep("level0", f)) == 0) { level <- "level0"}
    if (!length(grep("level1", f)) == 0) { level <- "level1"}
    if (!length(grep("level2", f)) == 0) { level <- "level2"}
    if (!length(grep("phylogeny", f)) == 0) { level <- "phylogeny.info"}
    if (!length(grep("RPA", f)) == 0) { method <- "RPA"}
    if (!length(grep("Sum", f)) == 0) { method <- "Sum"}
    if (!length(grep("Ave", f)) == 0) { method <- "log10Ave"}
    if (!length(grep("NMF", f)) == 0) { method <- "NMF"}

    svDialogs::tclServiceMode(FALSE)

  } else {
    if (level %in% c("level0", "level1", "level2", "species")) {
      f <- paste(data.dir, "/", level, "_", method, "_010.tab", sep = "")
    } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile_010.tab", sep = "")
    } else if (level == "phylogeny.info") {
      f <- paste(data.dir, "/phylogenyinfo_010.tab", sep = "")
    }
  }

  # level2_Sum_010.tab  
  # oligoprofile_010.tab  
  # species_log10Ave_010.tab

  # Read the data
  if (level == "phylogeny.info") {
    tab <- read.table(f, sep = "\t", header = T, row.names = NULL)
  } else {
    tab <- read.table(f, sep = "\t", header = T, row.names = 1)
  }

  # Check that the data is logarithmized as required in the arguments
  if (level == "phylogeny.info") {
    tab <- tab
  } else if (!length(grep("log10", f)) == 0 && !log10) { 
    tab <- 10^tab
  } else if (length(grep("log10", f)) == 0 && log10) {
    tab <- log10(tab)
  } else {
    tab <- tab
  }

  # Always impute by rows and for log10 data
  if (impute && any(is.na(tab))) {
    warning(paste("The matrix has ", sum(is.na(tab)), " missing values - imputing."))
    if (!log10) {
      tab <- 10^t(impute(t(log10(tab))))
    } else {
      tab <- t(impute(t(tab)))
    }
  }

  as.matrix(tab)

}

