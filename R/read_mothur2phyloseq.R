#' @title Read Mothur Output into a Phyloseq Object
#' @description Read mothur shared and consensus taxonomy files into a
#' \code{\link{phyloseq-class}} object.
#' @details Mothur shared and consensus taxonomy files will be converted to
#' \code{\link{phyloseq-class}}.
#' @param shared.file A
#' \href{http://www.mothur.org/wiki/Shared_file}{shared file}
#' produced by \emph{mothur}. Identified from the .shared extension
#' @param consensus.taxonomy.file Consensus taxonomy file
#' produced by \emph{mothur}. Identified from with the .taxonomy extension.
#' See \url{http://www.mothur.org/wiki/ConTaxonomy_file}.
#' @param mapping.file Metadata/mapping file with .csv extension
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#' #otu.file <- system.file(
#' #"extdata/Baxter_FITs_Microbiome_2016_fit.final.tx.1.subsample.shared",
#' #   package='microbiome')
#'
#' #tax.file <- system.file(
#' #"extdata/Baxter_FITs_Microbiome_2016_fit.final.tx.1.cons.taxonomy",
#' #   package='microbiome')
#'
#' #meta.file <- system.file(
#' #"extdata/Baxter_FITs_Microbiome_2016_mapping.csv",
#' #   package='microbiome')
#' 
#' #p0 <- read_mothur2phyloseq(
#' #       shared.file=otu.file,
#' #       consensus.taxonomy.file=tax.file,
#' #       mapping.file=meta.file)
#'
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
#'
read_mothur2phyloseq <- function(shared.file, consensus.taxonomy.file, 
    mapping.file=NULL) {
    
    ### Mothur Shared file to otu table ###
    
    m.otu <- read.table(shared.file, check.names=FALSE, header=TRUE, 
        sep="\t", stringsAsFactors=FALSE)
    m.otu$label <- NULL
    m.otu$numOtus <- NULL
    x <- m.otu$Group
    m.otu$Group <- NULL
    rownames(m.otu) <- NULL
    rownames(m.otu) <- x
    #m.otu[, c(2:(ncol(m.otu) - 1))] <- sapply(m.otu[, c(2:(ncol(m.otu) - 
    #    1))], as.numeric, 1)
    m.otu[, c(2:(ncol(m.otu) - 1))] <- vapply(m.otu[, c(2:(ncol(m.otu) - 
        1))], as.numeric, as.numeric(m.otu[,1]))    
    mothur_otu_table <- (otu_table(t(as.matrix(m.otu)), taxa_are_rows=TRUE))
    
    # message('Converted Shared to OTU table')
    
    ### Consensus Taxonomy file to taxa table ###
    
    # the file extension has to be .taxonomy!
    TaxonomyFile <- read.table(consensus.taxonomy.file, check.names=FALSE, 
        header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    TaxonomyFile$Taxonomy <- gsub("[\"]", "", TaxonomyFile$Taxonomy)
    TaxonomyFile$Taxonomy <- gsub("[(1-100)]", "", TaxonomyFile$Taxonomy)
    mothur_tax <- separate(TaxonomyFile, "Taxonomy", into=c("Kingdom", 
        "Phylum", "Order", "Class", "Family", "Genus"), sep=";",
        extra="merge")
    mothur_tax$Genus <- gsub(";", "", mothur_tax$Genus)
    mothur_tax$Size <- NULL
    
    rownames(mothur_tax) <- mothur_tax$OTU
    mothur_tax_mat <- as.matrix(mothur_tax)
    mothur_taxonomy=tax_table(mothur_tax_mat)
    
    # message('Converted Contaxonomy to Taxa table') message('reading
    # mapping file and creating a phyloseq object')
    
    ### Mapping file to sampledata table ###
    
    mothur_map_data <- NULL
    if (!is.null(mapping.file)) {
        mothur_map <- read.csv(mapping.file, row.names=1, check.names=FALSE)
        mothur_map_data <- sample_data(mothur_map)
    }
    
    ### Final phyloseq object ###
    pseq <- merge_phyloseq(mothur_otu_table, mothur_taxonomy, mothur_map_data)
    
    return(pseq)
    
}
