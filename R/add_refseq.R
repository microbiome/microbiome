#' @title Add \code{refseq} Slot for \code{dada2} based \code{phyloseq} Object
#' @description Utility to add refseq slot for \code{dada2} based 
#'              \code{phyloseq} Object. Here, the taxa_names which are unique 
#'              sequences, are stored in \code{refseq} slot of \code{phyloseq}.
#'              Sequence ids are converted to ids using tag option.  
#' @param x \code{\link{phyloseq-class}} object with sequences as rownames.
#' @param tag Provide name for Ids, Default="ASV".
#' @return \code{\link{phyloseq-class}} object 
#' @examples
#' 
#' # ps <- add_refseq(p0,tag="ASV")
#' # ps
#' 
#' @export
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
#' @importFrom Biostrings DNAStringSet
add_refseq <- function(x, tag="ASV"){
    
    if (class(x)!="phyloseq"){
        stop("Input is not an object of phyloseq class")
    }
    
    nucl <- Biostrings::DNAStringSet(taxa_names(x))
    names(nucl) <- taxa_names(x)
    x <- merge_phyloseq(x, nucl)
    
    rm(nucl)
    
    if(is.na(tag) || is.null(tag)){
        taxa_names(x) <- paste0("taxa", seq(ntaxa(x)))
        return(x)
    } else{
        taxa_names(x) <- paste0(tag, seq(ntaxa(x)))
        return(x)
    }
    
}