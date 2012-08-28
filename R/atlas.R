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
#'   @param minmax.quantiles Quantiles for minmax normalization. Default c(0.005, 0.995).
#'   @param bgc.method Optional. Background correction method. By default no background correction.
#'   @param mc.cores Number of cores for parallel computation
#'   @param remove.nonspecific.oligos Optional logical set to remove nonspecific oligos (TRUE). By default keeping the non-specific oligos (FALSE)
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
			       minmax.quantiles = c(0.005, 0.995), 
			       bgc.method = NULL, 
			       mc.cores = 3, 
			       remove.nonspecific.oligos = FALSE) 
{

  # library(microbiome); source("~/scripts/R/HITchip/atlas.R"); allowed.projects <- ListAtlasProjects()[1:2]; dbuser = 'lmlahti'; dbpwd = 'passu'; dbname = 'Phyloarray'; result.path <- "~/data/HITChip/Atlas/20120828/"; my.scaling = "minmax"; minmax.quantiles = c(0.005, 0.995); bgc.method = NULL; mc.cores = 3; remove.nonspecific.oligos = FALSE

  # Install new MySQL dump of the database with: 
  # mysql -u"dbuser" -p"dbpwd" dbname < dump.sql
  require(parallel)

  # Probes and species to exclude
  rm.phylotypes <- phylotype.rm.list("HITChip")

  # Define the chip type to use
  chiptype   <- "Agilent-016089" 

  full.data.file <- paste(result.path, "/atlas.full.RData", sep = "")
  training.data.file <- paste(result.path, "/atlas.train.RData", sep = "")
  test.data.file <- paste(result.path, "/atlas.test.RData", sep = "")
  parameter.data.file <- paste(result.path, "/atlas.parameters.RData", sep = "")
  a <- try(save(chiptype, file = parameter.data.file))
  if (!is.null(a)) {stop("Create result directory in advance!")}

  message("Extract sample information from the HITChip database")
  project.info <- fetch.sample.info(allowed.projects, chiptype, dbuser, dbpwd, dbname)

  oligomap <- get.oligomap(phylogeny = "16S", rm.phylotypes$oligos, dbuser, dbpwd, dbname, remove.nonspecific.oligos = remove.nonspecific.oligos)

  message("Get probe-level data for the selected hybridisations")
  tmp <- get.probedata(unique(project.info[["hybridisationID"]]), 
      	 		rm.phylotypes$oligos, dbuser, dbpwd, dbname, mc.cores = mc.cores)  

  fdat.orig <- 1 + tmp$data       # features x hybs, original non-log scale; ensure smallest value is 1
  fdat.oligoinfo <- tmp$info      # oligoinfo

  # Annotations for selected hybridisations
  fdat.hybinfo <- project.info[match(colnames(fdat.orig), project.info$hybridisationID), ]  	       	                    
  rownames(fdat.hybinfo) <- colnames(fdat.orig)

  # calculate quantile points in original scale 
  maxabs <- mean(apply(fdat.orig, 2, quantile, max(minmax.quantiles), na.rm = TRUE))
  minabs <- mean(apply(fdat.orig, 2, quantile, min(minmax.quantiles), na.rm = TRUE))  
  minmax.points <- c(minabs, maxabs)

  # Normalize (input required in log-scale)
  fdat <- ScaleProfile(fdat.orig, method = my.scaling, minmax.quantiles = c(0.005, 0.995))

  # Summarize probes into oligos and hybridisations into samples
  oligo.data <- summarize.rawdata(fdat, fdat.hybinfo, fdat.oligoinfo, 
  	     		oligo.ids = sort(unique(oligomap$oligoID)))
	 			
  # Background correction
  if (!is.null(bgc.method)) { 
    oligo.data <- oligo.bg.correction(oligo.data, bgc.method)
  }

  # First produce full preprocessed data matrices
  data.matrices.full <- list(oligo = oligo.data)
  for (level in c("species", "L1", "L2")) {
    for (method in c("ave", "sum", "rpa", "nmf")) { 
      data.matrices.full[[level]][[method]] <- summarize.probesets(oligomap, oligo.data, method, level, rm.phylotypes = rm.phylotypes)
    }
  }

  # Sample annotation matrix
  sample.info.full <- fdat.hybinfo[match(colnames(oligo.data), 
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

    data.matrices[[sample.set]] <- list(oligo = oligo.data[, samples])
    for (level in names(data.matrices.full)) {
      for (method in names(data.matrices.full[[level]])) {
        data.matrices[[sample.set]][[level]][[method]] <- data.matrices.full[[level]][[method]][, samples]
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
  params <- list(dbuser = dbuser, dbpwd = NA, dbname = dbname, 
  	         my.scaling = my.scaling, minmax.quantiles = minmax.quantiles, 
		 bgc.method = bgc.method, 
		 result.path = result.path, 
		 allowed.projects = allowed.projects, 
		 rm.oligos = rm.phylotypes$oligos, rm.phylotypes = rm.phylotypes, 
		 files = list(full.data.file, training.data.file, test.data.file, parameter.data.file),
		 session.info = session.info, minmax.points = minmax.points, date = date())

  # Save version info
  save(oligomap, params, file = parameter.data.file, compress = "xz")

  list(training.data = data.matrices[["training"]], 
       test.data = data.matrices[["test"]], 
       full.data = data.matrices.full, 
       sampleinfo = sample.info.full, 
       oligomap = oligomap, 
       parameters = params)

}




#' fetch.sample.info
#'
#' Description: Fetch sample information from HITChip atlas
#'
#' Arguments:
#'   @param allowed.projects list projects for which to fetch the data
#'   @param chiptype chiptype (eg. new.chip)
#'   @param dbuser MySQL user
#'   @param dbpwd MySQL password
#'   @param dbname MySqL database name
#'   @param selected.samples Sample to investigate. By default all.
#' Returns:
#'   @return project.info data.frame
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

fetch.sample.info <- function (allowed.projects, chiptype = NULL, 
		  dbuser, dbpwd, dbname, selected.samples = NULL) { 

  require(RMySQL)
  drv <- dbDriver("MySQL")
  con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)

  # Fetch all data from the database
   rs <- dbSendQuery(con, paste("SELECT p.projectName,p.projectID,s.subjectID,s.sampleID,s.samplingDate,s.normAlgVersion,h.hybridisationID,h.dye,a.arrayID,a.barcode,sl.designID,s.reproducibility,s.normalisationFinished,s.imageID,fe.extractionID,fe.extractionName,fe.noSampleNormalisation,h.isDiscarded,fe.hasReproCheck
     FROM sample s               
     JOIN hybridisation h USING (sampleID) JOIN featureextraction fe USING (hybridisationID)
     JOIN project p USING (projectID)
     JOIN array a USING (arrayID)
     JOIN slide sl USING (barcode)
     ORDER BY s.projectID, s.sampleID, h.hybridisationID, fe.extractionID"))

  message("Fetch selected projects and samples")
  project.info.all <- fetch(rs, n = -1)

  # If no chiptype specified, use all
  if (is.null(chiptype)) {chiptype <- unique(project.info.all$designID)}
  if (is.null(selected.samples)) {selected.samples <- unique(project.info.all$sampleID)}

  # Close MySQL connection
  dbDisconnect(con) 

  # Filter out samples based on predefined criteria
  rkeep <- project.info.all$projectName %in% allowed.projects &
           !as.logical(project.info.all$isDiscarded) &
           !as.logical(project.info.all$noSampleNormalisation) &
           as.logical(project.info.all$normalisationFinished) &
           project.info.all$hasReproCheck & 
           project.info.all$designID %in% chiptype &
           project.info.all$normAlgVersion == 1.1 &
    	   project.info.all$sampleID %in% selected.samples
 
  # Remove annotations which are identical for all samples
  ckeep <- sapply(project.info.all, function (x) {!length(unique(x)) == 1})

  message("Filter the data")
  project.info <- project.info.all[rkeep, ckeep]

  project.info                  
     
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
    training.samples <- lapply(project.names, function (pn) {sample.info[sample.info$projectName == pn, "sampleID"][training.inds[, pn]]})
  } else {
    stop("Duplicate sampleIDs. Handle before splitting into training and test sets.")
  }

  training.set <- unlist(training.samples)
  test.set <- setdiff(sample.info$sampleID, training.set)
 
  list(training = training.set, test = test.set)

}

