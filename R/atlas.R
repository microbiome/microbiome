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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

FetchHITChipAtlas <- function (allowed.projects, dbuser, dbpwd, dbname, 
		     	       result.path,
		     	       my.scaling = "minmax", 
			       minmax.quantiles = c(0.005, 0.995), 
			       bgc.method = NULL, 
			       mc.cores = 3, 
			       remove.nonspecific.oligos = FALSE) 
{

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

  oligo.map <- get.phylogeny(phylogeny = "16S", rm.phylotypes$oligos, dbuser, dbpwd, dbname, remove.nonspecific.oligos = remove.nonspecific.oligos)

  message("Get probe-level data for the selected hybridisations")
  tmp <- get.probedata(unique(project.info[["hybridisationID"]]), 
      	 		rm.phylotypes$oligos, dbuser, dbpwd, dbname, mc.cores = mc.cores)  
  fdat.orig <- tmp$data       # features x hybs, original non-log scale
  fdat.oligoinfo <- tmp$info  # oligoinfo

  # Annotations for selected hybridisations
  fdat.hybinfo <- project.info[match(colnames(fdat.orig), 
  	       	                     project.info$hybridisationID), ]
  rownames(fdat.hybinfo) <- colnames(fdat.orig)

  # Normalize (input required in log-scale)
  fdat <- ScaleProfile(log10(fdat.orig), method = my.scaling, 
       	  				 minmax.quantiles = minmax.quantiles)
 			
  # Summarize probes into oligos and hybridisations into samples
  oligo.data <- summarize.rawdata(fdat, fdat.hybinfo, fdat.oligoinfo, 
  	     			oligo.ids = sort(unique(oligo.map$oligoID)))
	 			
  # Background correction
  if (!is.null(bgc.method)) { 
    oligo.data <- oligo.bg.correction(oligo.data, bgc.method)
  }

  # First produce full preprocessed data matrices
  data.matrices.full <- list(oligo = oligo.data)
  for (level in c("species", "level 1", "level 2")) {
    for (method in c("ave", "sum", "rpa")) {
      data.matrices.full[[level]][[method]] <- summarize.probesets(oligo.map, oligo.data, method, level, rm.phylotypes = rm.phylotypes, rm.oligos = rm.phylotypes$oligos)
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
		 session.info = session.info)

  # Save version info
  save(oligo.map, params, file = parameter.data.file, compress = "xz")

  list(training.data = data.matrices[["training"]], 
       test.data = data.matrices[["test"]], 
       full.data = data.matrices.full, 
       sampleinfo = sample.info.full, 
       oligomap = oligo.map, 
       parameters = params)

}


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
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
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



#' Description: List projects in MySQL database
#'
#' Arguments:
#'   @param dbuser MySQL user
#'   @param dbpwd MySQL password
#'   @param dbname MySqL database name
#' Returns:
#'   @return project names vector
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.mysql.projects <- function (dbuser, dbpwd, dbname) { 

  require(RMySQL)
  drv <- dbDriver("MySQL")
  con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)

  # Fetch all data from the database
   rs <- dbSendQuery(con, paste("SELECT p.projectName FROM project p"))

  project.info <- fetch(rs, n = -1) 

  unique(project.info$projectName)

}



#' Description: Get phylogeny
#' 
#' Arguments:
#'   @param hybridization.ids Specify the hybridizations to retrieve
#'   @param rmoligos oligos to exclude
#'   @param dbuser MySQL user
#'   @param dbpwd  MySQL password
#'   @param dbname MySqL database name
#'   @param mc.cores Number of cores for multicore computing
#' Returns:
#'   @return list with data (features x hybridizations matrix) and info (features x info) fields 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.probedata <- function (hybridization.ids, rmoligos, dbuser, dbpwd, dbname, mc.cores = 1) {

  # List unique hybridisations for the selected samples
  hids <- mysql.format(hybridization.ids)
                      
  require(RMySQL)
  drv <- dbDriver("MySQL")
  con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)

  rs <- dbSendQuery(con, statement = paste("SELECT featureID,extractionID,fe.hybridisationID,spatNormSignal,isOutlier
      		FROM featuremeasurement 
		JOIN featureextraction fe USING (extractionID)
		JOIN hybridisation h USING (hybridisationID)
                JOIN arrayfeature af USE INDEX (PRIMARY) USING (featureID)
		WHERE fe.hybridisationID IN", hids))
  rawdata <- fetch(rs, n = -1)

  ## Check if there is any data
  rawdataDim <- dim(rawdata)
  if(rawdataDim[1]==0) {
    stop("No data found for these samples (perhaps they are not normalized yet?).\n\n")
  }

  message("Remove outliers")
  rawdata$spatNormSignal[as.logical(rawdata$isOutlier)] <- NA

  message("Split data into arrays")
  rawdata.esplit <- split(rawdata, rawdata$hybridisationID)

  message("Remove NAs")
  na.inds <- sapply(rawdata.esplit, function (x) all(is.na(x$spatNormSignal)))
  rawdata.esplit <- rawdata.esplit[!na.inds]

  # Get probeID - featureID - oligoID mappings
  rs <- dbSendQuery(con, "SELECT fe.featureID,p.probeID,p.oligoID,fe.arrayCol,fe.arrayRow FROM arrayfeature fe JOIN probe p USING (probeID)")
  probes <- fetch(rs, n = -1) 
  
  # Remove specified oligos
  probes <- probes[!probes$oligoID %in% rmoligos,]

  ftab.info <- data.frame(list(featureID = unique(rawdata$featureID)))
  ftab.info[["probeID"]] <- probes$probeID[match(ftab.info$featureID, probes$featureID)]
  ftab.info[["oligoID"]] <- probes$oligoID[match(ftab.info$probeID, probes$probeID)]

  message("Remove NA oligos")
  keep <- !is.na(ftab.info$oligoID)
  ftab.info <- ftab.info[keep, ]
  rownames(ftab.info) <- ftab.info$featureID

  # LL 4.4.2012. With HITChip atlas we encountered some cases where the arrays had different number of entries
  # due to duplicates on some arrays. Now added automated handling here to avoid problems with other array types
  # that may have different natural number of elements on the array.
  if (length(table(sapply(rawdata.esplit, nrow))) == 2) {
    ntmp <- max(sapply(rawdata.esplit, nrow))
    message(paste("Remove elements containing duplicated entries (", round(100*mean(!sapply(rawdata.esplit, nrow) == ntmp), 2), "%)", sep = ""))
    

    # ntmp == !10799 encountered with HITChip atlas, not yet elsewhere
    rawdata.esplit <- rawdata.esplit[!sapply(rawdata.esplit, nrow) == ntmp]
  } else if (length(table(sapply(rawdata.esplit, nrow))) > 2) {
    stop("Error 10799. Arrays are not comparable. Contact R package admins.")
  }

  # Form features x hybridizations matrix 
  inds <- match(rownames(ftab.info), rawdata.esplit[[1]]$featureID)
  ftab <- matrix(NA, nrow = nrow(ftab.info), ncol = length(rawdata.esplit))
  rownames(ftab) <- rownames(ftab.info)
  colnames(ftab) <- names(rawdata.esplit)
  for (hid in names(rawdata.esplit)) { ftab[, hid] <- I(rawdata.esplit[[hid]][inds, "spatNormSignal"]) }

  # Close MySQL connection
  dbDisconnect(con)

  # Clean up memory
  gc()

  list(data = ftab, info = ftab.info)
}



#' Description: Split data set into training and test samples
#' 
#' Arguments:
#'   @param sample.info Sample information table
#'   @param training.fraction Specify the fraction of the data to go into training set
#'   @param rseed random seed
#' Returns:
#'   @return List with two elements: training and test, listing the corresponding sample sets.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

pick.training.samples <- function (sample.info, training.fraction = 0.80, rseed = 1463) {

  # Split data into training and test sets
  if (sum(duplicated(sample.info$sampleID)) == 0) {
    set.seed(rseed)
    project.names <- unique(sample.info$projectName)
    training.inds <- sapply(table(sample.info$projectName), function (n) {sample(n, round(training.fraction*n))})
    training.samples <- lapply(project.names, function (pn) {sample.info[sample.info$projectName == pn, "sampleID"][training.inds[[pn]]]})
  } else {
    stop("Duplicate sampleIDs. Handle before splitting into training and test sets.")
  }

  training.set <- unlist(training.samples)
  test.set <- setdiff(sample.info$sampleID, training.set)
 
  list(training = training.set, test = test.set)

}


#' Description: Get phylogeny
#' 
#' Arguments:
#'   @param phylogeny phylogeny
#'   @param rmoligos oligos to exclude
#'   @param dbuser MySQL user
#'   @param dbpwd MySQL password
#'   @param dbname MySqL database name
#'   @param verbose verbose
#'   @param remove.nonspecific.oligos Logical. Remove oligos with multiple targets.
#'   @param chip chip type
#' Returns:
#'   @return oligo.map
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.phylogeny <- function (phylogeny, rmoligos = NULL, dbuser, dbpwd, dbname, verbose = TRUE, remove.nonspecific.oligos = FALSE, chip = "HITChip") {   

  if (verbose) { message("Load phylogeny info") }

  require(RMySQL)
  drv <- dbDriver("MySQL")
  con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)

  ## Collect the full phylogenetic information for oligos
  message("Collect the full phylogeny")
  excloligos <- ifelse(length(rmoligos>0),
                       paste('AND NOT (', paste("o.oligoID='",rmoligos,"'",sep="",collapse=" OR "), ') ', sep=''), '')
                       

  if (chip == "MITChip") {

    # Also get level0

    full16Squery <- paste("SELECT 		   		 
   		   		 l0.name AS 'level 0', 
   		                 l1.name AS 'level 1', ", "                        	 
   		                 l2.name AS 'level 2', ", "                        	 
				 species.name AS 'species', 
				 specimen.name AS 'specimen', 
				 ot.oligoID AS 'oligoID', ","
                        	 o.pmTm, 
				 ot.Tm, 
				 ot.mismatch, 
				 ot.complement ",

                        "FROM phylogeny ph ",

                        "JOIN taxon l0 USING (phylogenyID) ",
                        "JOIN taxtotax tt0 ON (tt0.parentID=l0.taxonID) ",

                        "JOIN taxon l1 ON (tt0.childID=l1.taxonID) ",
                        "JOIN taxtotax tt1 ON (tt1.parentID=l1.taxonID) ",

                        "JOIN taxon l2 ON (tt1.childID=l2.taxonID) ",
                        "JOIN taxtotax tt2 ON (tt2.parentID=l2.taxonID) ",

                        "JOIN taxon species ON (tt2.childID=species.taxonID) ",
                        "JOIN taxtotax tt3 ON (tt3.parentID=species.taxonID) ",
                        "JOIN taxon specimen ON (tt3.childID=specimen.taxonID) ",
                        "JOIN oligotargetpair ot ON (ot.targetID=specimen.targetID) ",
                        "JOIN oligo o ON (ot.oligoID=o.oligoID) ",

                        "AND l0.taxonLevel='level 0' ",

                        "AND tt0.nodeDistance=1 ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        "AND tt3.nodeDistance=1 ",

                        excloligos,
                        "ORDER BY l0.name, l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")
  } else {

    # get levels 1-3
    full16Squery <- paste("SELECT l1.name AS 'level 1', l2.name AS 'level 2', ", 
                        "species.name AS 'species', specimen.name AS 'specimen', ot.oligoID AS 'oligoID', ",
                        "o.pmTm, ot.Tm, ot.mismatch, ot.complement ",
                        "FROM phylogeny ph ",
                        "JOIN taxon l1 USING (phylogenyID) ",
                        "JOIN taxtotax tt1 ON (tt1.parentID=l1.taxonID) ",
                        "JOIN taxon l2 ON (tt1.childID=l2.taxonID) ",
                        "JOIN taxtotax tt2 ON (tt2.parentID=l2.taxonID) ",
                        "JOIN taxon species ON (tt2.childID=species.taxonID) ",
                        "JOIN taxtotax tt3 ON (tt3.parentID=species.taxonID) ",
                        "JOIN taxon specimen ON (tt3.childID=specimen.taxonID) ",
                        "JOIN oligotargetpair ot ON (ot.targetID=specimen.targetID) ",
                        "JOIN oligo o ON (ot.oligoID=o.oligoID) ",
                        "AND l1.taxonLevel='level 1' ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        excloligos,
                        "ORDER BY l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")

  }

  rs <- dbSendQuery(con, full16Squery)
  full16S <- fetch(rs, n = -1)
  
  # Close MySQL connection
  dbDisconnect(con)

  # Apply the standard filters
  if (verbose) {
    message("Hybridisation temperature: Tm >= pmTm - 2.5\n")
    message("Number of mismatches: 0\n")
    message("Must be a complement sequence\n")
    message("No requirement for a full-length hybridisation\n\n")
  }

  pruned16S <- prune16S(full16S, pmTm.margin = 2.5, complement = 1, mismatch = 0)

  rmoligos2 <- rmoligos
  if (remove.nonspecific.oligos) {
    if (verbose) {message("Removing oligos with multiple targets")}
    nPhylotypesPerOligo <- n.phylotypes.per.oligo(pruned16S, "level 2") 
    nonspecific.oligos <- setdiff(pruned16S$oligoID, names(which(nPhylotypesPerOligo == 1)))
    rmoligos2 <- c(rmoligos, nonspecific.oligos)
  } 

  pruned16S <- pruned16S[!pruned16S$oligoID %in% rmoligos2, ]

  pruned16S        

}


#' Description: filter 16S data
#' 
#' Arguments:
#'   @param full16S full16S
#'   @param pmTm.margin default 2.5
#'   @param complement logical
#'   @param mismatch logical
#'
#' Returns:
#'   @return filtered 16S data
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

prune16S <- function (full16S, pmTm.margin = 2.5, complement = 1, mismatch = 0) {
 
  keep <- full16S$Tm >= full16S$pmTm-pmTm.margin &
          full16S$complement == complement &
          full16S$mismatch == mismatch

  pruned16S <- full16S[keep, ]

  pruned16S
}


#' Description: Format string vector to mysql query format
#' 
#' Arguments:
#'   @param s string vector
#'
#' Returns:
#'   @return mysql query version
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities


mysql.format <- function (s) {
  paste("('",paste(as.character(s),collapse="','",sep="'"),"')",sep="")
}
  

#' Description: Check number of matching phylotypes for each probe
#' 
#' Arguments:
#'   @param oligo.map oligo - phylotype matching data.frame
#'   @param level phylotype level
#'
#' Returns:
#'   @return number of matching phylotypes for each probe
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

n.phylotypes.per.oligo <- function (oligo.map, level) {
  sapply(split(oligo.map[, c("oligoID", level)], oligo.map$oligoID), function(x) length(unique(x[[level]])))
}
  

#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param oligo.map oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param level summarization level
#'   @param verbose print intermediate messages
#'   @param rm.phylotypes Phylotypes to exclude (a list with fields species, level 1, level 2)
#'   @param rm.oligos Oligos to remove
#' Returns:
#'   @return summarized data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

summarize.probesets <- function (oligo.map, oligo.data, method, level, verbose = TRUE, rm.phylotypes = NULL, rm.oligos = NULL) {

  #level <- gsub(".", " ", level) # "level.1" -> "level 1"

  if ("level.0" %in% colnames(oligo.map)) {
    level0 <- "level.0"
  } else if ("level 0" %in% colnames(oligo.map)) {
    level0 <- "level 0"
  } else {
    level0 <- NULL
  }

  if ("level.1" %in% colnames(oligo.map)) {
    level1 <- "level.1"
  } else if ("level 1" %in% colnames(oligo.map)) {
    level1 <- "level 1"
  } else {
    level1 <- NULL
  }

  if ("level.2" %in% colnames(oligo.map)) {
    level2 <- "level.2"
  } else if ("level 2" %in% colnames(oligo.map)) {
    level2 <- "level 2"
  } else {
    level2 <- NULL
  }

  # Start by summarizing into species level
  rm.species <- unique(c(unique(oligo.map[oligo.map[[level1]] %in% rm.phylotypes[[level1]], "species"]), 
  	     		 unique(oligo.map[oligo.map[[level2]] %in% rm.phylotypes[[level2]], "species"]),
			 rm.phylotypes$species))
			 
  # Ensure that all L2 groups below specified L1 are removed as well
  rm.phylotypes[[level2]] <- unique(c(unique(oligo.map[oligo.map[[level1]] %in% rm.phylotypes[[level1]], level2]), 
			 rm.phylotypes[[level2]]))
			 	
  # Remove specified oligos
  if (!is.null(rm.oligos)) { oligo.data <- oligo.data[setdiff(rownames(oligo.data), rm.oligos), ]}
  oligo.map <- oligo.map[!oligo.map$oligoID %in% rm.oligos, ]

  # Get species matrix in original scale
  species.matrix <- 10^summarize.probesets.species(oligo.map, oligo.data, method, verbose, rm.species)
   
  if (level == "species") {

    summarized.matrix <- species.matrix

  } else if (level %in% c(level0, level1, level2)) {

    phylogroups <- level.species(level, oligo.map)
    
    # Remove the specified phylogroups
    phylogroups <- phylogroups[setdiff(names(phylogroups), rm.phylotypes[[level]])]

    summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(oligo.data))
    rownames(summarized.matrix) <- names(phylogroups)
    colnames(summarized.matrix) <- colnames(oligo.data)

    for (pg in names(phylogroups)) {
      specs <- unique(phylogroups[[pg]])
      mat <- matrix(species.matrix[specs,], nrow = length(specs))
      if (method == "ave") { vec <- colMeans(mat) }
      if (method == "sum") { vec <- colSums(mat)  } 
      if (method == "rpa") { vec <- colSums(mat)  } # For RPA, use the sum for L1/L2
      summarized.matrix[pg, ] <- vec
    }
  } else {

    message(level)
    message(nchar(level))
    message(colnames(oligo.map))
    stop("Provide proper level!")

  }

  # Return in the original log10 domain    
  log10(summarized.matrix)

}



#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param oligo.map oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose print intermediate messages
#'   @param rm.species Species to exclude
#' Returns:
#'   @return summarized data matrix in log10 scale
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

summarize.probesets.species <- function (oligo.map, oligo.data, method, verbose = TRUE, rm.species = c("Victivallis vadensis")) {

  level <- "species"			    

  probesets <- retrieve.probesets(oligo.map, level = level)
  probesets <- probesets[setdiff(names(probesets), rm.species)]
  
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(oligo.map, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets), ncol(oligo.data)), 
  		    	      dimnames = list(names(probesets), colnames(oligo.data))) 

  for (set in names(probesets)) {

    if (verbose) { message(set) }

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "rpa") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      vec <- rpa.fit(dat, sigma2.method = "robust", alpha = 1 + 0.1*ncol(oligo.data)/2, beta = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)$mu

    } else if (method == "ave") {

      vec <- log10(colMeans((10^dat), na.rm = T))

    } else if (method == "sum") {

      # Weight each probe by the inverse of the number of matching phylotypes
      # Then calculate sum -> less specific probes are downweighted
      # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      dat2 <- (10^dat) / nPhylotypesPerOligo[rownames(dat)]
      dat2[dat2 < 1] <- 1
      vec <- log10(colSums(dat2, na.rm = T))

    }
    
    summarized.matrix[set, ] <- vec 

  }

  summarized.matrix
  
}


#' Description: RPA for HITChip
#' 
#' Arguments:
#'   @param level level
#'   @param phylo phylo
#'   @param oligo.data oligo.data
#'
#' Returns:
#'   @return RPA preprocessed data
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

calculate.rpa <- function (level, phylo, oligo.data) {

  # List entities (e.g. species)
  phylo.list <- split(phylo, phylo[[level]])
  entities <- names(phylo.list)

  # initialize
  summarized.matrix <- array(NA, dim = c(length(entities), ncol(oligo.data)), dimnames = list(entities, colnames(oligo.data)))
  noise.list <- list() 

  for (entity in names(phylo.list)) {
    message(entity)

    # Pick expression for particular probes
    probes <- unique(phylo.list[[entity]][, "oligoID"])

    # oligo.data is already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 

    # dat: probes x samples
    if (nrow(dat) < 2) {
      vec <- as.vector(dat) # NOTE: circumvent RPA if there are no replicates 
      noise <- NA
    } else {
      res <- rpa.fit(dat)
      vec <- res$mu
      noise <- sqrt(res$sigma2)
      names(noise) <- probes
    }

    noise.list[[entity]] <- noise
    summarized.matrix[entity, ] <- vec #, epsilon, alpha, beta, sigma2.method, d.method)

  }

  list(emat = summarized.matrix, noise = noise.list)
  
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



#' Description: list corresponding probes
#' 
#' Arguments:
#'   @param name name
#'   @param level level
#'   @param phylo phylo
#'
#' Returns:
#'   @return probes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.probes <- function (name, level, phylo) {

  as.character(unique(phylo[which(phylo[[level]] == name), "oligoID"]))
  
}

#' Description: list corresponding entities
#' 
#' Arguments:
#'   @param level level
#'   @param phylo phylo
#'
#' Returns:
#'   @return entities
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.entities <- function (level, phylo) {

  as.character(unique(phylo[[level]]))

}

