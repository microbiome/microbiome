#' @title Import phyloseq Data 
#' @description Read the otu, taxonomy and metadata from various formats.
#' @param otu.file File containing the OTU table (for mothur this is the file
#' with the .shared extension)
#' @param taxonomy.file (for mothur this is typically the consensus taxonomy
#' file with the .taxonomy extension)
#' @param metadata.file File containing samples x variables metadata.
#' @param type Input data type: 'mothur' or 'simple' or 'biom' type.
#' @param sep CSV file separator
#' @return \code{\link{phyloseq-class}} object
#' @export
#' @details See help(read_mothur2phyloseq) for details on the Mothur input
#' format; and help(read_biom2phyloseq) for details on the biom format.
#' The simple format refers to the set of CSV files.
#' @examples 
#' #  pseq <- read_phyloseq(otu.file, 
#' #                     taxonomy.file, 
#' #                     metadata.file, 
#' #                     type=c('mothur', 'simple', 'biom'))
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_phyloseq <- function(otu.file=NULL, taxonomy.file=NULL,
    metadata.file=NULL, type=c("simple", "mothur", "biom"), sep = ",") {
    
    message("Time to complete depends on OTU file size")
    # TODO add automated recognition of the type?
    
    if (type == "mothur") {
        pseq <- read_mothur2phyloseq(otu.file, 
                                    taxonomy.file, 
                                    metadata.file)
    } else if (type == "simple") {
        pseq <- read_csv2phyloseq(otu.file, 
                                    taxonomy.file, 
                                    metadata.file, sep = sep)
    } else if (type == "biom") {
        pseq <- read_biom2phyloseq(biom.file = otu.file, 
                                    taxonomy.file, 
                                    metadata.file)
    } else {
        stop("Unrecognized type in read_phyloseq input. Exiting.")
    }
    
    pseq
    
}


