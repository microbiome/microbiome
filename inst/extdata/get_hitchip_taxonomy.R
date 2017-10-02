#' @title HITChip Taxonomy Table
#' @description Get taxonomy table for phylogenetic microarrays (mainly HITChip).
#' @param chip chip type (e.g. 'HITChip')
#' @param phylogeny.version 'full' or 'filtered' 
#'           (latter is the basis for species/L1/L2 summarization)
#' @param data.dir Data directory path
#' @return phylogeny mapping table
#' @export
#' @examples # tax.table <- get_hitchip_taxonomy('HITChip', 'full')
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get_hitchip_taxonomy <- function(chip, phylogeny.version = "full", data.dir = NULL) {

  # Was: GetPhylogeny
  
  hitchip.taxonomy <- NULL

    if (is.null(data.dir)) {
      data.dir <- system.file("extdata", package = "microbiome")
    }
    
    if (chip == "HITChip") {

      load(system.file("data/hitchip.taxonomy.rda", package = "microbiome"))
      tax.table <- hitchip.taxonomy[[phylogeny.version]]

      ## Phylogeny
      #f <- paste0(data.dir, "/taxonomy.", phylogeny.version, ".tab")
      #tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
      #tax.table <- tax_table(as.matrix(tab))    
      #tax.table <- as.data.frame(tax_table(tab))
      
      ## Get the phylogeny from Github url <-
      ## 'raw.github.com/microbiome/data/master/example-datasets/phylogeny' fnam
      ## <- paste(url, '.', phylogeny.version, '.tab', sep = '') 
      ## tax.table <- read.csv(text = RCurl::getURL(fnam), sep = '\t')
        
    } else {

      message(paste("get_hitchip_taxonomy not implemented for", chip))
      tax.table <- NULL  

    }

    df <- as.data.frame(tax.table)

    df
}


