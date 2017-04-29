#' @title Importing phyloseq Data 
#' @description Reads the otu, taxonomy and metadata from various formats.
#' @param otu.file File containing the OTU table (for mothur this is the file with the .shared extension)
#' @param taxonomy.file (for mothur this is typically the consensus taxonomy file with the .taxonomy extension)
#' @param metadata.file File containing samples x variables metadata.
#' @param type Input data type: "mothur" implemented now.
#' @return \code{\link{phyloseq-class}} object
#' @export
#' @details See help(read_mothur2phyloseq) for details on the Mothur input format.
#' @examples \dontrun{
#'   pseq <- read_phyloseq(...)
#'  }
#' @author Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
read_phyloseq <- function(otu.file = NULL, taxonomy.file = NULL, metadata.file = NULL, type = "mothur"){

  # TODO make read_phyloseq for CSV files outputted by write_phyloseq function
  # TODO add automated recognition of the type?

  if (type == "mothur") {
    pseq <- read_mothur2phyloseq(otu.file, taxonomy.file, metadata.file)
  }

  pseq
  
}
