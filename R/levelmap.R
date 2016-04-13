#' @title Map taxonomic levels
#' @description Map taxa between hierarchy levels.
#' @param taxa taxa to convert; 
#' 	  if NULL then considering all taxa in the tax.table
#' @param from convert from taxonomic level 
#' @param to convert to taxonomic level
#' @param tax.table tax.table
#' @return mappings
#' @examples 
#'   tax.table <- GetPhylogeny('HITChip', 'filtered')
#'   levelmap('Akkermansia', 'L2', 'L1', tax.table)
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

