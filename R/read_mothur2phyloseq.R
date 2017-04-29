#' @title Read and create a phyloseq object from mothur shared and consensus taxonomy files
#' @description Read and convert mothur output to Phyloseq object.
#' @details Mothur shared and consensus taxonomy files will be converted to \code{\link{phyloseq-class}}.
#' @param SharedFile =  Shared file ".shared file" extension. A \href{http://www.mothur.org/wiki/Shared_file}{shared file} produced by \emph{mothur}.
#' @param ConTaxonomy = consensus taxonomy file ".taxonomy" extension. Details \href{http://www.mothur.org/wiki/Constaxonomy_file}{consensus taxonomy file} produced by \emph{mothur}.
#' @param MappingFile = Metadata/mapping file in ".csv" extension
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     p0 <- read_mothur2phyloseq(SharedFile = "mothur_shared.shared", ConTaxonomy = "mothur_taxonomy.taxonomy", MappingFile = "mothur_mapping.csv")
#'           }
#' @keywords utilities
read_mothur2phyloseq <- function(SharedFile, ConTaxonomy, MappingFile){
  require(phyloseq);
  require(dplyr);
  ############### Mothur Shared file to otu table #############################
    m.otu <- read.table(SharedFile, check.names=FALSE, header = T, sep = "\t", stringsAsFactors = F)
    m.otu$label <- NULL
    m.otu$numOtus <- NULL
    x <- m.otu$Group
    m.otu$Group <- NULL
    rownames(m.otu) <- NULL
    rownames(m.otu) <- x
    m.otu[, c(2:(ncol(m.otu)-1))] <- sapply(m.otu[, c(2:(ncol(m.otu)-1))], as.numeric)
    mothur_otu_table <- (otu_table(t(as.matrix(m.otu)), taxa_are_rows=TRUE))
    print("Converted Shared to OTU table")

    ############### Consensus Taxonomy file to taxa table #############################
    # the file extension has to be .taxonomy!
    TaxonomyFile<-read.table(ConTaxonomy,check.names=FALSE, header = TRUE, sep = "\t",stringsAsFactors = FALSE)
    require(tidyr)
    head(TaxonomyFile)
    TaxonomyFile$Taxonomy <- gsub("[\"]", "", TaxonomyFile$Taxonomy)
    TaxonomyFile$Taxonomy <- gsub("[(1-100)]", "", TaxonomyFile$Taxonomy)
    mothur_tax <- separate(TaxonomyFile, "Taxonomy", into = c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus"), sep = ";", extra = "merge")
    mothur_tax$Genus <- gsub(";", "", mothur_tax$Genus)
    mothur_tax$Size <- NULL
    head(mothur_tax)
    rownames(mothur_tax) <- mothur_tax$OTU
    mothur_tax_mat<-as.matrix(mothur_tax)
    mothur_taxonomy = tax_table(mothur_tax_mat)
    print("Converted Contaxonomy to Taxa table")
    print("reading mapping file and creating a phyloseq object")
    ############### Mapping file to sampledata table #############################
    mothur_map <- read.csv(MappingFile,row.names=1,check.names=FALSE)
    mothur_map_data=sample_data(mothur_map)
    mothur_pseq <- merge_phyloseq(mot,mothur_taxonomy,mothur_map_data)
    print(mothur_pseq)
    return(mothur_pseq)
}
