#' @title Read BIOM File into a Phyloseq Object
#' @description Read biom and mapping files into a \code{\link{phyloseq-class}}
#' object.
#' @details Biom file and mapping files will be converted to
#' \code{\link{phyloseq-class}}.
#' @param biom.file A biom file with '.biom' extension 
#' @param taxonomy.file NULL the latest version has taxonomic information
#' within the biom 
#' @param metadata.file A simple metadata/mapping file with .csv extension
#' @param ... Arguments to pass for import_biom
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#' p0 <- read_biom2phyloseq() 
#' #biom.file <- qiita1629.biom"
#' #meta.file <- qiita1629_mapping.csv"
#' #p0 <- read_biom2phyloseq(biom.file = biom.file, 
#' #                       metadata.file = meta.file, 
#' #                       taxonomy.file = NULL)
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_biom2phyloseq <- function(biom.file = NULL, 
                            taxonomy.file = NULL, metadata.file = NULL, ...)
{

    if (is.null(biom.file)) {return(NULL)}

    levels <- c("Domain", "Phylum", 
            "Class", "Order", "Family", 
            "Genus")

    otu_biom <- import_biom(biom.file, ...)
    phyobj <- otu_biom

    if (!is.null(metadata.file)) {
        map <- read.csv(metadata.file, 
                check.names = FALSE, row.names = 1)
        s.map <- sample_data(map)
        phyobj <- merge_phyloseq(otu_biom, 
                        s.map)
    }

    taxtab <- tax_table(phyobj)    
    if (!is.null(taxonomy.file) && is.null(tax_table(phyobj))) {
        taxtab <- read_taxtable(taxonomy.file)
        tax_table(phyobj) <- tax_table(taxtab)    
    } else if (!is.null(taxonomy.file) && !is.null(tax_table(phyobj))) {
        warning("Taxonomy is available both in the biom file and in 
        the separate taxonomy.file. Using the biom file version here. 
        To change this original taxonomy, do it explicitly in your code 
        by modifying tax_table(physeq).")
    }

    if (ncol(tax_table(phyobj)) == 6) {
        colnames(tax_table(phyobj)) <- levels
    
        return(phyobj)
    } else if (ncol(taxtab) == 7) {
        colnames(tax_table(phyobj)) <- c(levels, 
                                    "Species")
    }
    
    return(phyobj)

}


