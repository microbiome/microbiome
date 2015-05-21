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
      f <- paste0(data.dir, "/taxonomy.", phylogeny.version, ".tab")
      tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
      tax.table <- tax_table(as.matrix(tab))    
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



#' levelmap
#' 
#' Map taxa between hierarchy levels
#'
#' @param taxa taxa to convert; 
#' 	  if NULL then considering all taxa in the tax.table
#' @param from convert from taxonomic level 
#' @param to convert to taxonomic level
#' @param tax.table tax.table
#'
#' @return mappings
#'
#' @examples 
#'   tax.table <- GetPhylogeny('HITChip', 'filtered')
#'   levelmap('Akkermansia', 'L2', 'L1', tax.table)
#'
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
levelmap <- function(taxa = NULL, from, to, tax.table) {

  # If taxonomy table is from phyloseq, pick the data matrix separately	 
  if (class(tax.table) == "taxonomyTable") {
    tax.table <- as.data.frame(tax_table(tax.table))
  }

  if (from == to) {
    df <- list()
    df[[to]] <- factor(taxa)
    df <- as.data.frame(df)
    return(df)
  }

  df <- tax_table(as.matrix(tax.table))    
  
  if (is.null(taxa)) {
    taxa <- as.character(unique(df[, from]))
  }

  # From higher to lower level
  if (length(unique(df[, from])) <= length(unique(df[, to]))) {

    sl <- list()
    for (pt in taxa) {

      inds <- which(as.vector(as.character(df[, from])) == pt)
      pi <- df[inds, to]
      sl[[pt]] <- as.character(unique(na.omit(pi)))

    }
    
  } else {

    # From lower to higher level
    inds <- match(as.character(taxa), df[, from])
    omap <- df[inds, ]
    sl <- omap[,to]

  }
     
  as.vector(sl@.Data)
    
}

