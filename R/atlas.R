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

#' FetchHITChipAtlas
#'
#' Description: Complete preprocessing pipeline for the HITChip atlas 
#' with default parameters
#'
#' Arguments:
#'   @param allowed.projects Vector of projectNames to include
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name
#'   @param result.path Directory for storing the results.
#'   @param my.scaling Normalization method. Default "minmax".
#'   @param mc.cores Number of cores for parallel computation
#'   @param remove.nonspecific.oligos Optional logical set to remove nonspecific oligos (TRUE). By default keeping the non-specific oligos (FALSE)
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#' Returns:
#'   @return  List with the following elements: 
#'     	      training.data: training data matrices (including alternative preprocessing methods)
#'	      test.data: test data matrices 
#'	      full.data: full data matrices	      
#'            atlas.sampleinfo: sample information matrix
#'	      parameters: input parameters and R session information
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

FetchHITChipAtlas <- function (allowed.projects, dbuser, dbpwd, dbname, 
		     	       result.path,
		     	       my.scaling = "minmax", 
			       mc.cores = 3, 
			       remove.nonspecific.oligos = FALSE, host = NULL, port = NULL) 
{

  # source("~/scripts/R/HITchip/atlas.R"); library(microbiome); fs <- list.files("~/Rpackages/microbiome/microbiome/R/", full.names = T); for (f in fs) {source(f)}; allowed.projects <- ListAtlasProjects(); dbuser = 'lmlahti'; dbpwd = 'passu'; dbname = 'Phyloarray'; result.path <- "~/tmp/"; my.scaling = "minmax"; mc.cores = 3; remove.nonspecific.oligos = FALSE

  # Install new MySQL dump of the database with: 
  # mysql -u"dbuser" -p"dbpwd" dbname < dump.sql
  InstallMarginal("parallel")

  # Probes and species to exclude
  rm.phylotypes <- phylotype.rm.list("HITChip")

  # Define the chip type to use
  chiptype   <- "Agilent-016089" 

  full.data.file <- paste(result.path, "/atlas.full.RData", sep = "")
  training.data.file <- paste(result.path, "/atlas.train.RData", sep = "")
  test.data.file <- paste(result.path, "/atlas.test.RData", sep = "")
  parameter.data.file <- paste(result.path, "/atlas.parameters.RData", sep = "")
  a <- try(save(chiptype, file = parameter.data.file))
  if (!is.null(a)) { stop("Create result directory in advance!") }

  message("Extract sample information from the HITChip database")
  project.info <- HITChipDB::fetch.sample.info(allowed.projects, chiptype, dbuser, dbpwd, dbname, host = host, port = port)

  phylogeny.info <- HITChipDB::get.phylogeny.info(phylogeny = "16S", rm.phylotypes$oligos, dbuser, dbpwd, dbname, remove.nonspecific.oligos = remove.nonspecific.oligos, host = host, port = port)

  message("Get probe-level data for the selected hybridisations")
  tmp <- HITChipDB::get.probedata(unique(project.info[["hybridisationID"]]), 
      	 		rm.phylotypes$oligos, dbuser, dbpwd, dbname, mc.cores = mc.cores, host = host, port = port)  

  fdat.orig <- tmp$data       # features x hybs, original non-log scale
  fdat.oligoinfo <- tmp$info      # oligoinfo

  # Annotations for selected hybridisations
  fdat.hybinfo <- project.info[match(colnames(fdat.orig), project.info$hybridisationID), ]  	       	                    
  rownames(fdat.hybinfo) <- colnames(fdat.orig)

  ## Discard the hybs that contain only NAs
  onlyNA <- colMeans(is.na(fdat.orig)) == 1
  naHybs <- colnames(fdat.orig)[onlyNA]
  if(sum(onlyNA) > 0){
    message("Removing the following hybs, because they contain only NAs:\n")
    message(naHybs,"\n\n")
    fdat.orig <- fdat.orig[, !onlyNA]
    fdat.hybinfo <- fdat.hybinfo[, !onlyNA]
  }

  # calculate quantile points in original scale 
  # hard-code to unify all analyses; these values were calculated manually from the HITChip atlas with 3200 samples, 
  # and rounded to 3 significant digits
  minmax.points <- c(30, 133000) 

  # Normalize (input required in log-scale)
  #fdat <- ScaleProfile(fdat.orig, method = my.scaling, minmax.quantiles = c(0.005, 0.995))
  fdat <- HITChipDB::ScaleProfile(fdat.orig, method = my.scaling, minmax.points = minmax.points)

  # Summarize probes into oligos and hybridisations into samples
  oligo.log10 <- HITChipDB::summarize.rawdata(log10(fdat), fdat.hybinfo, fdat.oligoinfo, 
  	     		oligo.ids = sort(unique(phylogeny.info$oligoID)))
	 
  # First produce full preprocessed data matrices
  data.matrices.full <- list(oligo = oligo.log10)
  for (level in c("species", "L1", "L2")) {
    for (method in c("sum", "rpa")) { 
      message(paste(level, method))
      data.matrices.full[[level]][[method]] <- summarize.probesets(phylogeny.info, oligo.log10, method, level, rm.phylotypes = rm.phylotypes)
    }
  }

  # Sample annotation matrix
  sample.info.full <- fdat.hybinfo[match(colnames(oligo.log10), 
  		                         fdat.hybinfo$sampleID),]
  rownames(sample.info.full) <- as.character(sample.info.full$sampleID)

  # -----------------------------------------------------------

  # Split the data into training and test sets
  set.seed(3652)
  splitted <- pick.training.samples(sample.info.full, training.fraction = 0.80)

  data.matrices <- list()
  sample.info <- list()

  for (sample.set in c("training", "test")) {

    samples <- splitted[[sample.set]]
    sample.info[[sample.set]] <- sample.info.full[samples,]

    data.matrices[[sample.set]] <- list(oligo = oligo.log10[, samples])
    for (level in names(data.matrices.full)) {
      for (method in names(data.matrices.full[[level]])) {
        if (!(level == "species" && method == "nmf")) {
          data.matrices[[sample.set]][[level]][[method]] <- data.matrices.full[[level]][[method]][, samples]
        } else {
          data.matrices[[sample.set]][[level]][[method]] <- NULL
        }
      }
    }
  }

  message(paste("Saving training data in ", training.data.file), sep = "")
  atlas <- data.matrices[["training"]]
  atlas.sampleinfo <- sample.info[["training"]]
  save(atlas, atlas.sampleinfo, file = training.data.file, compress = "xz")

  message(paste("Saving test data in ", test.data.file), sep = "")
  atlas <- data.matrices[["test"]]
  atlas.sampleinfo <- sample.info[["test"]]
  save(atlas, atlas.sampleinfo, file = test.data.file, compress = "xz") 

  message(paste("Saving full data matrix in ", full.data.file), sep = "")
  atlas <- data.matrices.full
  atlas.sampleinfo <- sample.info.full 
  save(atlas, atlas.sampleinfo, file = full.data.file, compress = "xz")

  # Save parameters
  session.info <- sessionInfo()
  params <- list(dbuser = dbuser, dbpwd = NA, dbname = dbname, host = host, port = port,
  	         my.scaling = my.scaling, # minmax.quantiles = minmax.quantiles, 
		 minmax.points = minmax.points, 
		 result.path = result.path, 
		 allowed.projects = allowed.projects, 
		 rm.oligos = rm.phylotypes$oligos, rm.phylotypes = rm.phylotypes, 
		 files = list(full.data.file, training.data.file, test.data.file, parameter.data.file),
		 session.info = session.info, minmax.points = minmax.points, date = date())

  # Save version info
  save(phylogeny.info, params, file = parameter.data.file, compress = "xz")

  list(training.data = data.matrices[["training"]], 
       test.data = data.matrices[["test"]], 
       full.data = data.matrices.full, 
       sampleinfo = sample.info.full, 
       phylogeny.info = phylogeny.info, 
       parameters = params)

}






#' pick.training.samples
#'
#' Description: Split data set into training and test samples
#' 
#' Arguments:
#'   @param sample.info Sample metadata data frame. Include the fields 'sampleID' and 'projectName', specifying the sample identifier and project name for each sample. 
#'   @param training.fraction Specify the fraction of the data to go into training set
#'   @param rseed random seed
#' Returns:
#'   @return List with two elements: training and test, listing the corresponding sample sets.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

pick.training.samples <- function (sample.info, training.fraction = 0.80, rseed = 1463) {

  if (!"sampleID" %in% colnames(sample.info)) {
    warning("sampleID field missing - using row names as identifiers")
    sample.info$sampleID <- rownames(sample.info)
  }

  if (!"projectName" %in% colnames(sample.info)) {
    warning("projectName field missing - all samples classified in MyProject")
    sample.info$projectName <- rep("MyProject", nrow(sample.info))
  }

  # Split data into training and test sets
  if (sum(duplicated(sample.info$sampleID)) == 0) {
    set.seed(rseed)
    project.names <- unique(sample.info$projectName)
    training.inds <- sapply(table(sample.info$projectName), function (n) {sample(n, round(training.fraction*n))})
    training.samples <- lapply(project.names, function (pn) {sample.info[sample.info$projectName == pn, "sampleID"][training.inds[,pn]]})
  } else {
    stop("Duplicate sampleIDs. Handle before splitting into training and test sets.")
  }

  training.set <- unlist(training.samples)
  test.set <- setdiff(sample.info$sampleID, training.set)
 
  list(training = training.set, test = test.set)

}

