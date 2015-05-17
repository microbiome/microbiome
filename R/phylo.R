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
#'   tax.table <- GetPhylogeny('HITChip', 'full')
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
      tax.table <- polish.tax.table(tab)
      
      # Get the phylogeny from Github url <-
      # 'raw.github.com/microbiome/data/master/example-datasets/phylogeny' fnam
      # <- paste(url, '.', phylogeny.version, '.tab', sep = '') 
      # tax.table <- read.csv(text = RCurl::getURL(fnam), sep = '\t')
        
    } else {

        message(paste("GetPhylogeny not implemented for", chip))
        tax.table <- NULL  

    }
    
    tax.table

}



#' levelmap
#' 
#' Map phylotypes between hierarchy levels
#'
#' @param phylotypes phylotypes to convert; 
#' 	  if NULL then considering all phylotypes in the tax.table
#' @param from convert from taxonomic level 
#' @param to convert to taxonomic level
#' @param tax.table tax.table
#'
#' @return mappings
#'
#' @examples 
#'   tax.table <- GetPhylogeny('HITChip', 'filtered')
#'   levelmap(phylotypes = 'Akkermansia', 'L2', 'L1', tax.table)
#'
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
levelmap <- function(phylotypes = NULL, from, to, tax.table) {

  # If taxonomy table is from phyloseq, pick the data matrix separately	 
  if (class(tax.table) == "taxonomyTable") {
    tax.table <- tax_table(tax.table)@.Data
  }

  if (from == to) {
    df <- list()
    df[[to]] <- factor(phylotypes)
    df <- as.data.frame(df)
    return(df)
  }

  tax.table <- polish.tax.table(tax.table)

  if (is.null(phylotypes)) {
    phylotypes <- as.character(unique(tax.table[, from]))
  }

    # From higher to lower level
    if (length(unique(df[[from]])) <= length(unique(df[[to]]))) {

        sl <- list()
        for (pt in phylotypes) {

	    inds <- which(as.vector(as.character(tax.table[, from])) == pt)
	    pi <- tax.table[inds, to]
            sl[[pt]] <- as.character(unique(na.omit(pi)))
        }
    
    } else {

      # From lower to higher level
      inds <- match(as.character(phylotypes), tax.table[[from]])
      omap <- tax.table[inds, ]
      sl <- omap[[to]]

    }
     
    sl
    
}

#' polish.tax.table
#'
#' Ensure tax.table is in correct format
#' 
#' @param tax.table tax.table data frame
#'
#' @return polished tax.table
#' @references See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @examples 
#'   #tax.table <- GetPhylogeny('HITChip', 'filtered')
#'   #tax.table <- polish.tax.table(tax.table)
#' @keywords internal
polish.tax.table <- function(tax.table) {
    
    colnames(tax.table)[
          which(colnames(tax.table) == "level.0")] <- "L0"
    colnames(tax.table)[
          which(colnames(tax.table) == "level.1")] <- "L1"
    colnames(tax.table)[
          which(colnames(tax.table) == "level.2")] <- "L2"
    
    colnames(tax.table)[
          which(colnames(tax.table) == "level 0")] <- "L0"
    colnames(tax.table)[
          which(colnames(tax.table) == "level 1")] <- "L1"
    colnames(tax.table)[
          which(colnames(tax.table) == "level 2")] <- "L2"

    # Fix some names	  
    if ("L2" %in% colnames(tax.table)) {
      tax.table$L2 <- gsub("^Clostridiales$", "Clostridium \\(sensu stricto\\)", tax.table$L2)
    }

    # Convert into phyloseq taxonomyTable format
    tax.table <- tax_table(as.matrix(tax.table))    
    
    as.data.frame(tax.table)

}


#' retrieve.probesets
#' 
#' List probes for each probeset
#'
#' @param tax.table data.frame with oligo - phylotype 
#' 	  		 mapping information
#' @param level phylotype level for probesets
#' @param name specify phylotypes to check (optional)
#'
#' @return A list. Probes for each phylotype.
#'
#' @examples 
#'   tax.table <- GetPhylogeny('HITChip')
#'   sets <- retrieve.probesets(tax.table, 'species', 'Weissella confusa')
#'                         
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

retrieve.probesets <- function(tax.table, level = "species", name = NULL) {

    # If name not given, pick all
    if (is.null(name)) {
        name <- unique(as.character(tax.table[[level]]))
    }
    
    phylo <- tax.table[tax.table[[level]] %in% name, ]
    
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
