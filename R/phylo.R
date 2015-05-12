#' GetPhylogeny
#' 
#' Get Chip phylogeny
#'
#'   @param chip chip type (e.g. 'HITChip')
#'   @param phylogeny.version 'full' or 'filtered' 
#'           (latter is the basis for species/L1/L2 summarization)
#'   @param data.dir Data directory path
#'
#'   @return phylogeny mapping table
#'
#' @export
#'
#' @examples 
#'   phylogeny.info <- GetPhylogeny('HITChip', 'full')
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GetPhylogeny <- function(chip, phylogeny.version = "full", data.dir = NULL) {

    if (is.null(data.dir)) {
      data.dir <- system.file("extdata", package = "microbiome")
    }
    
    if (chip == "HITChip") {
        
      # Phylogeny
      f <- paste0(data.dir, "/phylogeny.", phylogeny.version, ".tab")
      tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
      phylogeny.info <- polish.phylogeny.info(tab)
      
      # Get the phylogeny from Github url <-
      # 'raw.github.com/microbiome/data/master/example-datasets/phylogeny' fnam
      # <- paste(url, '.', phylogeny.version, '.tab', sep = '') 
      # phylogeny.info <- read.csv(text = RCurl::getURL(fnam), sep = '\t')
        
    } else {

        message(paste("GetPhylogeny not implemented for", chip))
        phylogeny.info <- NULL  

    }
    
    phylogeny.info
    
}



#' levelmap
#' 
#' Map phylotypes between hierarchy levels
#'
#' @param phylotypes phylotypes to convert; 
#' 	  if NULL then considering all phylotypes in the phylogeny.info
#' @param level.from convert from 
#' 	  Options: 'L0', 'L1', 'L2', 'species', 'oligo'
#' @param level.to conver to Options: 'L0', 'L1', 'L2', 'species', 'oligo'
#' @param phylogeny.info phylogeny.info
#'
#' @return mappings
#'
#' @examples 
#'   phylogeny.info <- GetPhylogeny('HITChip', 'filtered')
#'   levelmap(phylotypes = 'Akkermansia', 'L2', 'L1', phylogeny.info)
#'                   
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

levelmap <- function(phylotypes = NULL, level.from, level.to, phylogeny.info) {
    
    if (level.from == level.to) {
        df <- list()
        df[[level.to]] <- factor(phylotypes)
        df <- as.data.frame(df)
        return(df)
    }
    
    if (level.from == "level 0") {
        level.from <- "L0"
    }
    if (level.from == "level 1") {
        level.from <- "L1"
    }
    if (level.from == "level 2") {
        level.from <- "L2"
    }
    if (level.from == "oligo") {
        level.from <- "oligoID"
    }
    
    if (level.to == "level 0") {
        level.to <- "L0"
    }
    if (level.to == "level 1") {
        level.to <- "L1"
    }
    if (level.to == "level 2") {
        level.to <- "L2"
    }
    if (level.to == "oligo") {
        level.to <- "oligoID"
    }
    
    phylogeny.info <- polish.phylogeny.info(phylogeny.info)
    
    if (is.null(phylotypes)) {
        phylotypes <- as.character(unique(phylogeny.info[[level.from]]))
    }
    
    if (level.from == "species" && level.to %in% c("L0", "L1", "L2")) {
        sl <- species2levels(phylotypes, phylogeny.info)[, level.to]
    }
    
    if (level.from == "oligoID" && 
        level.to %in% c("L0", "L1", "L2", "species")) {
        sl <- oligoTOhigher(phylotypes, phylogeny.info, level.to = level.to)
    }
    
    if (level.from == "L2" && level.to == "L1") {
        sl <- level2TOlevel1(phylotypes, phylogeny.info)[, 2]
    }
    
    if (level.from == "L2" && level.to == "L0") {
        sl <- level2TOlevel0(phylotypes, phylogeny.info)[, 2]
    }
    
    if (level.from == "L1" && level.to == "L0") {
        sl <- level1TOlevel0(phylotypes, phylogeny.info)[, 2]
    }
    
    if (level.from == "L1" && level.to == "L2") {
        sl <- list()
        for (pt in phylotypes) {
            sl[[pt]] <- as.character(unique(
                phylogeny.info[phylogeny.info[["L1"]] == 
                pt, "L2"]))
        }
    }
    
    if (level.from == "L1" && level.to == "L0") {
        sl <- list()
        for (pt in phylotypes) {
            sl[[pt]] <- as.character(unique(
                phylogeny.info[phylogeny.info[["L1"]] == 
                pt, "L0"]))
        }
    }
    
    if (level.from == "L0" && level.to %in% c("L1", "L2")) {
        sl <- list()
        for (pt in phylotypes) {
            sl[[pt]] <- as.character(unique(
                phylogeny.info[phylogeny.info[[level.from]] == 
                pt, level.to]))
        }
    }
    
    if (level.from %in% c("L0", "L1", "L2") && level.to == "species") {
        sl <- list()
        for (pt in phylotypes) {
            sl[[pt]] <- as.character(unique(
                phylogeny.info[phylogeny.info[[level.from]] == 
                pt, level.to]))
        }
    }
    
    
    if (level.from %in% c("L0", "L1", "L2", "species") && 
        level.to == "oligoID") {
        sl <- list()
        for (pt in phylotypes) {
            sl[[pt]] <- as.character(unique( 
                phylogeny.info[phylogeny.info[[level.from]] == 
                pt, level.to]))
        }
    }
    
    sl
    
}

#' retrieve.probesets
#' 
#' List probes for each probeset
#'
#' @param phylogeny.info data.frame with oligo - phylotype 
#' 	  		 mapping information
#' @param level phylotype level for probesets
#' @param name specify phylotypes to check (optional)
#'
#' @return A list. Probes for each phylotype.
#'
#' @examples 
#'   phylogeny.info <- GetPhylogeny('HITChip')
#'   sets <- retrieve.probesets(phylogeny.info, 'species', 'Weissella confusa')
#'                         
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

retrieve.probesets <- function(phylogeny.info, level = "species", name = NULL) {

    # If name not given, pick all
    if (is.null(name)) {
        name <- unique(as.character(phylogeny.info[[level]]))
    }
    
    phylo <- phylogeny.info[phylogeny.info[[level]] %in% name, ]
    
    if (is.factor(phylo[[level]])) {
        phylo[[level]] <- droplevels(phylo[[level]])
    }
    
    phylo.list <- split(phylo, phylo[[level]])
    probesets <- lapply(phylo.list, function(x) {
        as.character(unique(x$oligoID))
    })
    names(probesets) <- names(phylo.list)
    
    probesets
    
} 
