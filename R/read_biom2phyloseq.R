#' @title Read BIOM File into a Phyloseq Object
#' @description Read biom and mapping files into a \code{\link{phyloseq-class}} object.
#' @details Biom file and mapping files will be converted to \code{\link{phyloseq-class}}.
#' @param otu.file A biom file with ".biom" extension 
#' @param taxonomy.file NULL the latest version has taxonomic information within the biom 
#' @param metadata.file A simple metadata/mapping file with .csv extension
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#'   \dontrun{
#'     # Example data
#'     library(microbiome)
#'     p0 <- read_biom2phyloseq(otu.file = "otu_table.biom",
#'                              metadata.file = "mapping.csv",
#'				taxonomy.file = NULL)
#'   }
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_biom2phyloseq <- function(otu.file = NULL, taxonomy.file = NULL, metadata.file = NULL){
  otu_biom <- import_biom(otu.file, parseFunction=parse_taxonomy_default)
  map <-read.csv(metadata.file,check.names=FALSE, row.names = 1)
  s.map <- sample_data(map)
  phyobj <-merge_phyloseq(otu_biom, s.map)
  if (ncol(tax_table(phyobj)) == 6){
    colnames(tax_table(phyobj)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
    return(phyobj)
  } else if (ncol(tax_table(phyobj)) == 7){
    colnames(tax_table(phyobj)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  return(phyobj)
}
