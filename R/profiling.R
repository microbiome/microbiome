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



#' Description: Profiling preprocessing script
#'
#' Arguments:
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name
#'   @param mc.cores Optional. Number of cores if parallelization is used.
#'   @param verbose monitor processing through intermediate messages
#'
#' Returns:
#'   @return Preprocessed data and parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

preprocess.chipdata <- function (dbuser, dbpwd, dbname, mc.cores = 1, verbose = TRUE) {

  # for Phyloarray database version 0.8 or 0.9

  ## ask parameters or read from R-file
  con <- dbConnect(dbDriver("MySQL"), username = dbuser, password = dbpwd, dbname = dbname)

  params <- ReadParameters(con)  
  params$chip <- detect.chip(dbname)
  params$rm.phylotypes <- phylotype.rm.list(params$chip) # List oligos and phylotypes to remove by default

  # Get sample information matrix for the selected projects	
  project.info <- fetch.sample.info(params$prj$projectName, chiptype = NULL, 
  	       	  		    dbuser, dbpwd, dbname, 
  	       	  		    selected.samples = params$samples$sampleID)

  message("Get probe-level data for the selected hybridisations")
  tmp <- get.probedata(unique(project.info[["hybridisationID"]]), params$rm.phylotypes$oligos, dbuser, dbpwd, dbname, mc.cores = mc.cores)
  fdat.orig <- tmp$data       # features x hybs, original non-log scale
  fdat.oligoinfo <- tmp$info  # oligoinfo

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
  
  ##################################
  ## GET OLIGO-PHYLOTYPE MAPPINGS
  ##################################

  # This handles also pmTm, complement and mismatch filtering
  pruned16S <- get.phylogeny(params$phylogeny, 
    	       		     rmoligos = params$rm.phylotypes$oligos, 
	    		     dbuser, dbpwd, dbname, verbose = verbose, 
			     remove.nonspecific.oligos = params$remove.nonspecific.oligos, 
			     chip = params$chip)

  ##################
  ## COMPUTE SCALING
  ##################

  # selected scaling for featurelevel data
  # Background correction after this step, if any. 
  # Order of normalization / bg correction was validated empirically.
  # bg.adjust intentionally set to NULL here. 
  # bg correction done _after_ oligo summarization, if any (see next steps)
  d.scaled <- ScaleProfile(log10(fdat.orig), params$normalization, bg.adjust = NULL) 

  ####################
  ## COMPUTE SUMMARIES
  ####################

  # Summarize probes into oligos and hybridisations into samples
  d.oligo2 <- summarize.rawdata(d.scaled, 
  	      			fdat.hybinfo, 
				fdat.oligoinfo = fdat.oligoinfo, 
				oligo.ids = sort(unique(pruned16S$oligoID)))

  # Then apply background correction if required:
  # if (!is.null(params$bgc.method)) { d.oligo2 <- oligo.bg.correction(d.oligo2, bgc.method = params$bgc.method) }
  # d.oligo2 <- oligo.bg.correction(d.oligo2, bgc.method = NULL)

  oligo.log10 <- d.oligo2
  # Return to the original scale
  oligo.abs <- matrix(10^d.oligo2, nrow = nrow(d.oligo2)) # - 1  
  rownames( oligo.abs ) <- rownames( d.oligo2 )
  colnames( oligo.abs ) <- colnames( d.oligo2 )

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.abs
  levels <- c("species", "level 2", "level 1")
  if (params$chip == "MITChip") {levels <- c(levels, "level 0")}
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa", "ave")) {

    	summarized.log10 <- summarize.probesets(pruned16S, oligo.log10, 
      			       	          method = method, level = level, 	
					  rm.phylotypes = params$rm.phylotypes,
					  rm.oligos = params$rm.phylotypes$oligos)

        # Store the data in absolute scale					  
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  list(data = finaldata, phylogeny = pruned16S, naHybs = naHybs, params = params)

}



#' Description: Detect the chip type (H/M/PITChip) from database name
#'
#' Arguments:
#'   @param dbname MySQL database name
#' 
#' Returns:
#'   @return chip name
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

detect.chip <- function (dbname) {

  if (dbname %in% c("PhyloArray_MIT", "phyloarray_mit", "Phyloarray_MIT")) {
    chip <- "MITChip"
  } else if (dbname %in% c("PhyloArray_PIT", "phyloarray_pit", "Phyloarray_PIT", "pitchipdb")) {
    chip <- "PITChip"
  } else if (dbname %in% c("PhyloArray_HIT", "phyloarray_hit", "Phyloarray_HIT", "Phyloarray")) {
    chip <- "HITChip"
  } else {
    stop("Check database name (dbname)!")
  }

  chip

}



#' Description: Define parameters in select box, or bring from earlier session
#'
#' Arguments:
#'   @param con Output from dbConnect(dbDriver("MySQL"), username = dbuser, password = dbpwd, dbname = 'PhyloArray')
#' 
#' Returns:
#'   @return list with defined parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

ReadParameters <- function (con) {

  scaling <- list.scaling.methods()	       

  ## Determine the working directory
  wdir <- tclvalue(tkchooseDirectory(title = "Save output files into directory:")) 
        
  ## Choose samples to display
  prj <- microbiome::choose.projects(con, multi = TRUE, condition = NULL)

  if(nrow(prj) < 1) { stop("Choose at least 1 project") }
  samples <- choose.samples(con, multi=TRUE, title='Select samples', 
  	       			   condition=list(list(field='projectID', value=prj$projectID)))
  if(nrow(samples) < 2) { stop("Choose at least 2 samples") }

  defaults <- list(phylogeny = "16S", remove.nonspecific.oligos = FALSE, normalization = "minmax")
  s <- NULL; for (nam in names(defaults)) {s <- paste(s, paste(nam, ":", defaults[[nam]], sep = ""), "; ", sep = "")}

  use.default.parameters <- tk_select.list(c(paste("Yes, use the defaults:", s), "No, proceed to parameter selection"), preselect = paste("Yes, use the defaults:", s), multiple = FALSE, title = paste('Use default parameters?'))

  rs <- dbSendQuery(con, "SELECT phylogenyID, name FROM phylogeny WHERE NOT name='ROOT'")
  phylogenies <- fetch(rs, n = -1)
  phylogenies <- unique(phylogenies$name)
  defaults$phylogeny <- phylogenies[grep(defaults$phylogeny, phylogenies)][[1]] 

  if (substr(use.default.parameters, 1, 2) == "No") {    

    ## Choose the phylogeny and the (lowest) summary taxonomic level to use
    if (length(phylogenies) > 1) {
      phylogeny <- tk_select.list(phylogenies, preselect = defaults$phylogeny, multiple = FALSE, title = 'Select phylogeny for profiling')
    } else {
      phylogeny <- phylogenies[[1]]
    }

    # Exclude non-specific oligos?
    remove.nonspecific.oligos <- tk_select.list(c("Yes", "No"), multiple = FALSE, preselect = defaults$remove.nonspecific.oligos, title = "Remove non-specific oligos?")
    if (remove.nonspecific.oligos == "Yes") {remove.nonspecific.oligos <- TRUE}
    if (remove.nonspecific.oligos == "No") {remove.nonspecific.oligos <- FALSE}
   
    ## Normalization method
    scal <- tk_select.list(names(scaling), preselect = defaults$normalization, multiple = FALSE, title = "Select normalization method")

  } else {

    phylogeny <- defaults$phylogeny
    remove.nonspecific.oligos <- defaults$remove.nonspecific.oligos
    scal <- defaults$normalization

  }

  # LL and JS decided to remove default BG correction 29.3.2012
  # based on various validation tests.
  # bgc.method <- select.list(c("2*sd bkg intensity", "6*sd bkg intensity"), multiple = FALSE, preselect = "6*sd bkg intensity", title = "Select background correction method:")
  bgc.method <- NULL # Intentional

  list(wdir = wdir, prj = prj, samples = samples, phylogeny = phylogeny, normalization = scal, bgc.method = bgc.method, remove.nonspecific.oligos = remove.nonspecific.oligos)

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



#' Description: List scaling methods
#'
#' Arguments:
#'
#' Returns:
#'   @return List of scaling methods
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.scaling.methods <- function () {

  list('none'='none',
                #'minimum/median'='minmed',
                'minimum/maximum'='minmax',
                'minmax'='minmax',
                #'median'='med',
                'quantile'='quant'
                #'normExp+MedianFC'='normExpMedianFC',
                #'normExp+quantile'='normExpQuant'
   )

}


#' Description: List clustering metrics
#'
#' Arguments:
#'
#' Returns:
#'   @return list of clustering metrics
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.clustering.metrics <- function () {

  list('Pearsons correlation coefficient'='correlation',
                 'euclidian'='euclidian')
}



#' Description: List color scales
#'
#' Arguments:
#'
#' Returns:
#'   @return list of color scales
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

list.color.scales <- function () {
  ## Different colour scales
  list('white/blue'=colorRampPalette(c("white","darkblue"),interpolate='linear')(100),
       'white/black'=colorRampPalette(c("white","black"),interpolate='linear')(100),
       'black/yellow/white'=colorRampPalette(c("black","yellow","white"),bias=0.5,interpolate='linear')(100))

}



#' Description: Fetch data from the database
#'
#' Arguments:
#'   @param params params 
#'   @param con con
#'   @param scriptVersion scriptVersion
#'   @param save.data save.data
#'   @param scaling scaling
#'   @param cmetrics cmetrics
#'
#' Returns:
#'   @return data 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

FetchData <- function (params, con, scriptVersion, save.data, scaling, cmetrics) {

  ## COLLECTING DATA FROM THE DATABASE
  message("Collecting data from the database\n")

  ## Collecting data for ALL probes (used to be: oligo's) JN09032011
  inclsamples <- paste("sampleID='", params$samples$sampleID, "'", sep = "", collapse = " OR ")

  rs <- dbSendQuery(con,"DROP TABLE IF EXISTS tmp1")
  query <- paste("CREATE TEMPORARY TABLE tmp1 ",
                 "SELECT s.sampleID, s.projectID, h.hybridisationID, fe.extractionID, h.dye, p.oligoID, p.probeID, af.featureID, ",
                 " fm.isOutlier, spatNormSignal AS featureSignal ",
                 "FROM sample s USE INDEX (PRIMARY) ",
                 "JOIN hybridisation h USING (sampleID) ",
                 "JOIN featureextraction fe USING (hybridisationID) ",
                 "JOIN featuremeasurement fm USING (extractionID) ",
                 "JOIN arrayfeature af USE INDEX (PRIMARY) USING (featureID) ",
                 "JOIN probe p USING (probeID) ",
                 "WHERE (",inclsamples,") ",
                 "AND normalisationFinished ",
                 "AND NOT (oligoID IS NULL) ",
                 "AND NOT isDiscarded ",
                 "AND NOT noSampleNormalisation ",
                 "ORDER BY s.sampleID, h.hybridisationID,  p.probeID, af.featureID",
                 sep ="")

  rs <- dbSendQuery(con, query)
  rs <- dbSendQuery(con, "ALTER TABLE tmp1 ADD INDEX (featureID)")
  rs <- dbSendQuery(con, 'SELECT * FROM tmp1')

  rawdata <- fetch(rs, n = -1)
  rawdataDim <- dim(rawdata)

  ## Check if there is any data
  if(rawdataDim[1]==0)
    stop("No data found for these samples (perhaps they are not normalized yet?).\n\n")

  ## Create the data matrix (featuretab) for clustering based on all array features, 
  ## each hybridisation having one column in this table and each feature having one row. 
  ## The first column contains the oligoID's
  message("Create the FULL data matrix (featuretab) for clustering.\n")

  ## Change outlier values to NAs
  rawdata$featureSignal[as.logical(rawdata$isOutlier)] <- NA

  ## Get the hybIDs and initialize featuretab with the first hyb
  hybIDs <- unique(rawdata$hybridisationID)
  featuretab <- rawdata[rawdata$hybridisationID==hybIDs[1], c("oligoID","probeID","featureID", "featureSignal")]

  ## Name the column with sampleID, hybID, and dye
  samplename <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"sampleID"])
  dye <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"dye"])
  projectname <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"projectID"])
  columnnames <- paste(samplename, hybIDs[1], dye, projectname, sep=".")

  message(str(featuretab))
  if (length(hybIDs)>1) {
    for (i in hybIDs[2:length(hybIDs)]) {
      addcol <- rawdata[rawdata$hybridisationID==i,"featureSignal"]
      featuretab <- cbind(featuretab, addcol)

      ## Name the columns with sampleID, hybID, and dye
      samplename <- unique(rawdata[rawdata$hybridisationID==i,"sampleID"])
      dye <- unique(rawdata[rawdata$hybridisationID==i,"dye"])
      projectname <- unique(rawdata[rawdata$hybridisationID==i,"projectID"])
      columnnames <- c(columnnames, paste(samplename, i, dye, projectname, sep="."))
    }
  }

  colnames(featuretab) = c("oligoID","probeID","featureID", columnnames)
  rownames(featuretab) <- featuretab$featureID 

  ## Discard the hybs contains only NAs
  onlyNA <- colSums(is.na(featuretab))==dim(featuretab)[1]
  naHybs <- names(onlyNA)[onlyNA]
  if(sum(onlyNA)>0){
    cat("Removing the following hybs, because they contain only NAs:\n")
    cat(naHybs,"\n\n")
    featuretab <- featuretab[,!onlyNA]
  }

  # Remove rmoligos
  featuretab <- featuretab[!featuretab$oligoID %in% params$rm.phylotypes$oligos, ]

  ## Write log of parameters used in profiling in to the file
  tmpTime <- strsplit(as.character(Sys.time()), split=" ")[[1]]
  tmpDate <- tmpTime[1]
  tmpTime <- paste(strsplit(tmpTime[2], split=":")[[1]], collapse=".")
  profTime <- paste(tmpDate,tmpTime,sep="_")
  logfilename <- paste(params$wdir,"/",profTime,"_profiling_log.txt", sep="")

  cat("Log of profiling script\n", "\n", file=logfilename)
  cat("profiling date: ",profTime, "\n", file=logfilename, append=T)
  cat("script version: ",scriptVersion,  "\n",file=logfilename, append=T)
  cat("data retrieved from db: ",params$useDB,  "\n", file=logfilename, append=T)
  cat("project IDs: ",params$prj$projectID,  "\n", file=logfilename, append=T)
  cat("sample IDs: ",params$samples$sampleID,  "\n", file=logfilename, append=T)
  cat("excluded oligos: ",params$rmoligos,  "\n", file=logfilename, append=T)
  cat("excluded hybridisations: ",naHybs,  "\n", file=logfilename, append=T)
  cat("phylogeny: ",params$phylogeny,  "\n", file=logfilename, append=T)
  cat("scaling: ",params$scal,  "\n", file=logfilename, append=T)
  cat("clustering tree in: ",params$clusterGraphFile,  "\n", file=logfilename, append=T)
  cat("tree ratio: ",params$figureratio, "\n",file=logfilename, append=T)
  cat("clustering metric: ",params$clmet, "\n",file=logfilename, append=T)
  cat("phylogeny level in figure: ",params$lev, "\n",file=logfilename, append=T)
  cat("figure coloring: ", params$pal, "\n",file=logfilename, append=T)
  cat("figure fontsize: ", params$fontsize, "\n",file=logfilename, append=T)
  cat("data saved: ", save.data, "\n",file=logfilename, append=T)
  cat("data in directory: ",params$wdir, "\n",file=logfilename, append=T)

  ## Write parameters used in profiling to the file
  paramfilename <- paste(params$wdir,"/",profTime,"_profiling_params.Rdata", sep="")

  save(logfilename, profTime, scriptVersion, params, naHybs, scaling, cmetrics, save.data, file=paramfilename)
  
  ## Collect the full phylogenetic information for oligos  
  cat("Collect the full 16S phylogeny\n")

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
                        "WHERE ph.name='",params$phylogeny,"' ",
                        "AND l1.taxonLevel='level 1' ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        "ORDER BY l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")
  rs <- dbSendQuery(con, full16Squery)
  full16S <- fetch(rs, n = -1)
  full16S <- full16S[full16S$oligoID %in% params$rm.phylotypes$oligos,]          
	  
  message("FINISHED COLLECTING THE DATA\n")

  list(full16S = full16S, featuretab = featuretab, logfilename = logfilename)

}



#' Description: List probes for each probeset
#'
#' Arguments:
#'   @param phylo data.frame with oligo - phylotype mapping information
#'   @param level phylotype level for probesets
#' Returns:
#'   @return A list. Probes for each phylotype.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

retrieve.probesets <- function (phylo, level = "species") {

  # phylo <- pruned16S
  phylo.list <- split(phylo, phylo[[level]])
  probesets <- lapply(phylo.list, function(x) {unique(x$oligoID)})	     
  names(probesets) <- names(phylo.list)

  probesets
}



#' Description: Calculate species summaries and possibly update d.oligo2
#'
#' Arguments:
#'   @param d.oligo2 d.oligo2
#'   @param bgc.method background correction method
#' Returns:
#'   @return Background-corrected data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

oligo.bg.correction <- function (d.oligo2, bgc.method) {

  if ( bgc.method == "2*sd bkg intensity" ){ bgth <- 2 }
  if ( bgc.method == "6*sd bkg intensity" ){ bgth <- 6 }

  d.oligo2 <- threshold.data(d.oligo2, bgth)
  d.oligo2 <- apply(d.oligo2, c(1,2), function(x) max(0, x))
  
  d.oligo2

}


#' Description: Calculate d.oligo2
#'
#' Arguments:
#'   @param featuretab featuretab
#'   @param d.scaled d.scaled
#'   @param oligo.ids oligo.ids
#' Returns:
#'   @return d.oligo2
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.doligo2 <- function (featuretab, d.scaled, oligo.ids) {

  # FIXME: replace into profiling script with d.oligo3 which is used in atlas  

  d.oSplit <- split(cbind(featuretab[,1:3],d.scaled), featuretab$oligoID)
  d.oSplit.pruned <- d.oSplit[oligo.ids]
  d.oligo <- t(sapply(d.oSplit, function(x) apply((x[,4:dim(x)[2]]), 2, mean, na.rm=T))) # hybs separate
  sampleID <- get.sampleid(d.oligo)

  d.oligo2 <- t(sapply(d.oSplit,
                     function(x){
                       temp <- apply((x[,4:dim(x)[2]]), 2, mean, na.rm=T)
                       temp2 <- sapply(split(temp, sampleID), mean, na.rm=T)
                       return(temp2)
                     }
                     ))# hybs averaged
  
  d.oligo2
}


#' Description: Compress normalized raw data matrix into final probe-level matrix:
#'              summarize oligos into probes and hybridisations into samples
#'
#' Arguments:
#'   @param fdat normalized raw data matrix oligos x hybridisations
#'   @param fdat.hybinfo hybridization info table
#'   @param fdat.oligoinfo oligo info table
#'   @param oligo.ids oligo.ids
#' Returns:
#'   @return probes x samples matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

summarize.rawdata <- function (fdat, fdat.hybinfo, fdat.oligoinfo, oligo.ids) {

  # List rows for each oligo (each oligo has multiple features which will be averaged)
  d.oSplit <- split(1:nrow(fdat), fdat.oligoinfo$oligoID)[oligo.ids] 

  # probes x hybs: oligo summary as means of log feature signals per oligo, hybs separate
  message("probe summary as means of log feature signals per oligo, hybs separate")
  oligo.data  <- t(sapply(d.oSplit, function(x) colMeans(fdat[x,], na.rm = TRUE)))

  ## Average over all hybridisations/extractions associated with this sample
  # List hybridisations associated with the same sample
  indlist <- split(1:ncol(oligo.data), fdat.hybinfo$sampleID)

  # Keep only cases with multiple (typically 2) hybs per sample (see table(sapply(indlist, length)))
  message("Removing samples with only one hybridisation")
  indlist <- indlist[sapply(indlist, length) > 1]

  message("Forming probes x samples matrix (average over hybridisations for each sample)")
  oligo.data <- sapply(indlist, function(inds) { rowMeans(oligo.data[, inds]) } )

  oligo.data
}


#' Description: Pick sampleIDs from d.oligo column names
#'
#' Arguments:
#'   @param d.oligo d.oligo matrix
#' Returns:
#'   @return sampleID vector
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.sampleid <- function (d.oligo) {

  sampleID <- sapply(colnames(d.oligo), function(z) { # edit by S.Tims
    # allows samples with "." in the samplename
    s <- strsplit(z, split="\\.")[[1]] 
    if(length(s)>4){
       s <- head(s,-3)
       s <- paste(s, collapse = ".")
    } else { s <- s[1] }
   return(s)})

   sampleID
}

#' Description: Load/install necessary packages and check OS
#'
#' Arguments:
#'
#' Returns:
#'   @return operating system string
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities


check.dependencies <- function () {

  ## determine if using windows or mac/linux, based on the path style
  if(strsplit(Sys.getenv()["R_HOME"],split="")[[1]][1]=="/"){
    os <- "unix"
  } else {
    os <- "win"
  }

  os
}

#' Description: Between-arrays normalization 
#'
#' Arguments:
#'   @param r.feature data matrix in logarithmic scale
#'   @param method normalization method
#'   @param bg.adjust background adjustment 
#'   @param minmax.quantiles quantiles for minmax
#' Returns:
#'   @return Normalized data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

ScaleProfile <- function (r.feature, method = 'minmax', bg.adjust = NULL, minmax.quantiles = c(0.005, 0.995)) {

  method <- list.scaling.methods()[[method]]

  message(paste("Normalizing with", method))
            
  ## Table r.feature is a copy of featuretab containing 
  ## logarithms of the values in 
  ## featuretab

  if (sd(na.omit(r.feature[,1])) > 100) {
    warning("Please check that the input matrix to ScaleProfile is in logarithmic scale")
  }

  r <- r.feature

    if (method=='minmax') {
      r <- scaling.minmax(r.feature, quantile.points = minmax.quantiles, robust = FALSE)
    } else if (method=='minmax.robust') {
      r <- scaling.minmax(r.feature, quantile.points = minmax.quantiles, robust = TRUE)
    } else if (method=='quant') {
      dn <- dimnames(r)
      r <- normalize.quantiles(r)
      dimnames(r) <- dn
    } else if (method=='normExpQuant') {
            ## Impute NA's with sample medians
            na.inds <- which(is.na(r), arr.ind=T)
            r <- apply(r,2,function(x){x[is.na(x)] <- median(x, na.rm=T); return(x)})
            rc <- apply(10^(r), 2, bg.adjust)
            dn <- dimnames(r)
            r <- normalize.quantiles(log10(rc+1))
            dimnames(r) <- dn
            r[na.inds] <- NA
     } else {
       stop("No between-array normalization recognized!!")
     }
 
  return(r)
}


#' Description: Minmax scaling. 
#'
#' Arguments:
#'   @param r data matrix in log10 scale
#'   @param quantile.points quantiles for minmax
#'   @param robust Select minmax version. 
#' 
#' Returns:
#'   @return normalized data matrix
#'
#' @note With robust = FALSE, the standard minmax is carried out. This
#'  shifts and scales each array such that their min and max values are
#'  identical across arrays. The robust = TRUE will perform the scaling
#'  such that the upper quantiles and minimum values of the data match (instead of
#'  maximum values).
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

scaling.minmax <- function (r, quantile.points, robust = FALSE) {

  # return original scale 
  rc <- 10^r 

  maxabs <- mean(apply(rc, 2, quantile, max(quantile.points), na.rm = TRUE))
  minabs <- mean(apply(rc, 2, quantile, min(quantile.points), na.rm = TRUE))

  if (!robust) {

    r <- apply(rc, 2, function (x) { 
                 x = (((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))*(maxabs - minabs)) + minabs;
       return(x)})

  } else {
    
    r <- apply(rc, 2, function (x) {

    # Shift data to start from zero
    xz <- x - min(x, na.rm = T);
    # Check the quantile points
    maxq <- quantile(xz, max(quantile.points), na.rm = TRUE);   
    # Determine the scaling factor such that the max quantiles will match between arrays
    k <- maxabs/maxq;
    # Scale the data to match max quantiles
    xs <- k * xz + min(rc, na.rm = TRUE);
    xs})
  }

  # Return to the input scale (log10)
  log10(r)
}






#' Description: determine detection threshold for the data
#'
#' Arguments:
#'   @param dat data
#'   @param sd.times standard deviation threshold
#'
#' Returns:
#'   @return thresholded data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

threshold.data <- function(dat, sd.times = 6){

  thr <- apply(dat, 2, function(x){
      DD <- density(as.numeric(x),adjust=1.2,na.rm=T);
      noise_mode <- DD$x[which(DD$y==max(DD$y))[1]];
      noise_sd   <- sd(x[x < noise_mode],na.rm=T);
      low.thresh <- noise_mode + sd.times*noise_sd;
      low.thresh 
      })

  # Subtract background from signal intensities in each sample
  data.mat<-t(apply(dat, 1, function(Tr){ Tr-thr })) 
  return(data.mat)
}



#' Description: get background parameters
#'
#' Arguments:
#'
#' Returns:
#'   @return TBA
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

get.bkg.params <- function(){
  tt <- tktoplevel()
  tkwm.title(tt,"Background Subtraction")

  frm <- populate.radiobuttons(tt,title="Background Subtraction Method",var.names=c("min. 500 oligos","2*sd bkg intensity","6*sd bkg intensity","none"),var.values=c("min. 500 oligos","2*sd bkg intensity","6*sd bkg intensity","none"),var.init=tclVar("2*sd bkg intensity"))

  frm2 <- tkframe(tt)

  frm2.up <- populate.radiobuttons(frm2,title="Handling of negatives in ave data",var.names=c("Keep negative values","Set negatives to zero"),var.values=c("TRUE","FALSE"),var.init=tclVar("FALSE"))

  frm2.down <- tkframe(frm2)
  button.OK <- tkbutton(frm2.down, text="OK", command=function(){
    tkdestroy(tt)
  })

  tkpack(frm2.down, button.OK)
  tkpack(frm2.up$frame,frm2.down,side="top",pady=5) 
  tkpack(frm2,frm$frame,side="left",padx=5) 
  tkpack(frm2)
  tkwait.window(tt) 
  return(list(method=tclvalue(frm$var),keep.neg=as.logical(tclvalue(frm2.up$var))))
}





#' Description: Default list of removed phylotypes and oligos
#'
#' Arguments:
#'  @param chip Chip name (HIT/MIT/PIT/Chick)Chip
#' Returns:
#'   @return List of removed oligos and phylotypes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

phylotype.rm.list <- function (chip) {

  rm.phylotypes <- list()

  if (chip == "HITChip") {
    
    rm.phylotypes[["oligos"]] <- c("UNI 515", "HIT 5658", "HIT 1503", "HIT 1505", "HIT 1506")
    rm.phylotypes[["species"]] <- c("Victivallis vadensis")
    rm.phylotypes[["level 1"]] <- c("Lentisphaerae")
    rm.phylotypes[["level 2"]] <- c("Victivallis")

  } else if (chip == "MITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["level 1"]] <- c()
    rm.phylotypes[["level 2"]] <- c()

  } else if (chip == "PITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["level 1"]] <- c()
    rm.phylotypes[["level 2"]] <- c()

  } else if (chip == "ChickChip") {
    warning("No universal probes excluded from ChichChip yet!")
  }

  rm.phylotypes

}




#' Description: Stability analysis. Calculates average Pearson '
#  correlation between samples in the input data and picks the lower '
#  triangular matrix to avoid duplicating the correlations. Returns 
#  correlations and stability estimate (average of the correlations).
#'
#' Arguments:
#'   @param dat data matrix phylotypes vs. samples
#'
#' Returns:
#'   @return List with correlations and astability estimate
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

calculate.stability <- function (dat) {
  cors <- lower.triangle(cor(dat))
  list(correlations = cors, stability = mean(cors))
}

