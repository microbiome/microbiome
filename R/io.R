# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <leo.lahti@iki.fi>. All rights reserved.

# This file is a part of the microbiome R package

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Function to read run.profiling.script output into R
#'
#' Arguments:
#'   @param level phylogenetic level ("oligo" / "species" / "level1" / "level2" / "level0") or "phylogeny"
#'   @param method ("rpa" / "sum")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @examples # params <- run.profiling.script(...); dat <- read.profiling.data(params$wdir, "level1", "rpa")
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
read.profiling <- function(level, method, data.dir, log10 = FALSE){

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

#' Description: Write log file
#'
#' Arguments:
#'   @param naHybs hybridisation that were removed due to NAs  
#'   @param params parameters
#'
#' Returns:
#'   @return List of scaling methods
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

WriteLog <- function (naHybs, params) {

  scaling <- list.scaling.methods()
  scriptVersion <- sessionInfo()$otherPkgs$microbiome$Version # microbiome package number

  ## Write log of parameters used in profiling in to the file
  tmpTime <- strsplit(as.character(Sys.time()), split=" ")[[1]]
  tmpDate <- tmpTime[1]
  tmpTime <- paste(strsplit(tmpTime[2], split=":")[[1]], collapse=".")
  profTime <- paste(tmpDate,tmpTime,sep="_")
  logfilename <- paste(params$wdir,"/",profTime,"_profiling_log.txt", sep="")

  cat("Log of profiling script\n", "\n", file=logfilename)
  cat("profiling date: ",profTime, "\n", file=logfilename, append=T)
  cat("script version: ", scriptVersion,  "\n",file=logfilename, append=T)
  cat("data retrieved from db: ",params$useDB,  "\n", file=logfilename, append=T)
  cat("project IDs: ",params$prj$projectID,  "\n", file=logfilename, append=T)
  cat("sample IDs: ",params$samples$sampleID,  "\n", file=logfilename, append=T)
  cat("excluded oligos: ",params$rm.phylotypes$oligos,  "\n", file=logfilename, append=T)
  cat("excluded species: ",params$rm.phylotypes[["species"]], "\n", file=logfilename, append=T)
  cat("excluded level 1: ",params$rm.phylotypes[["level 1"]], "\n", file=logfilename, append=T)
  cat("excluded level 2: ",params$rm.phylotypes[["level 2"]], "\n", file=logfilename, append=T)
  cat("excluded hybridisations: ",naHybs,  "\n", file=logfilename, append=T)
  cat("remove non-specific oligos: ",params$remove.nonspecific.oligos, "\n",file=logfilename, append=T)
  cat("phylogeny: ",params$phylogeny,  "\n", file=logfilename, append=T)
  cat("scaling: ",params$scal,  "\n", file=logfilename, append=T)
  cat("data in directory: ",params$wdir, "\n",file=logfilename, append=T)

  # Now graphics are outside of this function
  # cat("clustering tree in: ",params$clusterGraphFile,  "\n", file=logfilename, append=T)
  # cat("tree ratio: ",params$figureratio, "\n",file=logfilename, append=T)
  # cat("clustering metric: ",params$clmet, "\n",file=logfilename, append=T)
  # cat("phylogeny level in figure: ",params$lev, "\n",file=logfilename, append=T)
  # cat("figure coloring: ", params$pal, "\n",file=logfilename, append=T)
  # cat("figure fontsize: ", params$fontsize, "\n",file=logfilename, append=T)
  # cat("data saved: ", save.data, "\n",file=logfilename, append=T)

  ## Save profiling parameters 
  paramfilename <- paste(params$wdir,"/",profTime,"_profiling_params.Rdata", sep="")
  save(logfilename, profTime, scriptVersion, params, naHybs, file = paramfilename)  

  list(log.file = logfilename, parameter.file = paramfilename)
  
}




#' Description: Write preprocessed data into the output directory
#'
#' Arguments:
#'   @param finaldata preprocessed data matrices in absolute scale (from the preprocess.chipdata function)
#'   @param output.dir output directory
#'   @param phylogeny phylogeny
#'   @param verbose verbose
#'
#' Returns:
#'   @return Preprocessed data in absolute scale, phylogeny, and parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

WriteChipData <- function (finaldata, output.dir, phylogeny, verbose = TRUE) {

  ## featurelevel data
  ## WriteMatrix(cbind(fdat.oligoinfo, d.scaled), paste("featureprofile_", scriptVersion, ".tab", sep = ""), params)

  ## Write oligoprofile in original (non-log) domain
  fname <- paste(output.dir, "/oligoprofile.tab", sep = "")
  mydat <- finaldata[["oligo"]]
  WriteMatrix(cbind(rownames(mydat), mydat), fname, verbose)
    
  ## Write the other levels in log domain
  for (level in setdiff(names(finaldata), "oligo")) {
    for (method in names(finaldata[[level]])) {
      
      fname <- paste(output.dir, "/", level, "-", method, ".tab", sep = "")
      mydat <- finaldata[[level]][[method]]
      WriteMatrix(cbind(rownames(mydat), mydat), fname, verbose)

    }
  }

  ## Write oligo specificity at species level (the number of species for the oligo targets)
  fname <- paste(output.dir, "/phylogenyinfo.tab", sep = "")
  nSpeciesPerOligo <- sapply(split(phylogeny, phylogeny$oligoID), function(x) length(unique(x$species)))
  WriteMatrix(cbind(phylogeny, nSpeciesPerOligo = nSpeciesPerOligo[phylogeny$oligoID]), fname, verbose)

  # Return path to the output directory 
  output.dir
 
}


#' Description: Write matrix in tab file
#'
#' Arguments:
#'   @param dat data matrix
#'   @param filename output file
#'   @param verbose verbose
#' Returns:
#'   @return output file location
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

WriteMatrix <- function (dat, filename, verbose = FALSE) { 

  if (verbose) { message(paste("Writing output in ", filename)) }
  write.table(dat, file = filename, quote = FALSE, sep = "\t", row.names = FALSE)
  filename

}
