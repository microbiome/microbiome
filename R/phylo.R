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
      #tax.table <- as.data.frame(tax_table(tab))
      
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


#' polish.tax.table
#'
#' Ensure tax.table is in correct format
#' 
#' @param tax.table tax.table data frame
#'
#' @return polished tax.table
#' @references See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
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
      tax.table[, "L2"] <- gsub("^Clostridiales$", 
      		  	        "Clostridium \\(sensu stricto\\)", 
				tax.table[, "L2"])
    }

    # Convert into phyloseq taxonomyTable format
    tax.table <- tax_table(as.matrix(tax.table))    
    
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
    tax.table <- as.data.frame(tax_table(tax.table))
  }

  if (from == to) {
    df <- list()
    df[[to]] <- factor(phylotypes)
    df <- as.data.frame(df)
    return(df)
  }

  #df <- polish.tax.table(tax.table)
  df <- tax_table(as.matrix(tax.table))    
  
  if (is.null(phylotypes)) {
    phylotypes <- as.character(unique(df[, from]))
  }

  # From higher to lower level
  if (length(unique(df[, from])) <= length(unique(df[, to]))) {

    sl <- list()
    for (pt in phylotypes) {

      inds <- which(as.vector(as.character(df[, from])) == pt)
      pi <- df[inds, to]
      sl[[pt]] <- as.character(unique(na.omit(pi)))

    }
    
  } else {

    # From lower to higher level
    inds <- match(as.character(phylotypes), df[, from])
    omap <- df[inds, ]
    sl <- omap[,to]

  }
     
  sl
    
}

