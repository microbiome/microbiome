#' @title OTU Prevalence 
#' @description Simple prevalence measure.
#' @param x A vector, data matrix or \code{\link{phyloseq}} object
#' @param sort Sort the groups by prevalence
#' @param count Logical. Indicate prevalence as fraction of samples
#' (in percentage [0, 1]; default); or in absolute counts indicating
#' the number of samples where the OTU is detected (strictly) above the given
#' abundance threshold.
#' @inheritParams core
#' @details For vectors, calculates the fraction (count=FALSE) or
#' number (count=TRUE) of samples that exceed the
#' detection. For matrices, calculates this for each matrix
#' column. For phyloseq objects, calculates this for each OTU. The
#' relative prevalence (count=FALSE) is simply the absolute
#' prevalence (count=TRUE) divided by the number of samples.
#' @return For each OTU, the fraction of samples where a given OTU is
#' detected. The output is readily given as a percentage.
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' data(peerj32)     
#' pr <- prevalence(peerj32$phyloseq, detection=0, sort=TRUE, count=TRUE)
prevalence <- function(x, detection=0, sort=FALSE, count=FALSE,
    include.lowest=FALSE) {
    
    if (is.null(detection)) {
        detection <- (-Inf)
    }
    
    if (is.null(x)) {
        warning("x is NULL - returning NULL")
        return(NULL)
    }
    
    # Convert phyloseq to matrix
    if (is.phyloseq(x)) {
        x <- abundances(x)
    }
    
    if (is.vector(x)) {
        
        if (include.lowest) {
            prev <- sum(x >= detection)
        } else {
            prev <- sum(x > detection)
        }
        
    } else if (is.matrix(x) || is.data.frame(x)) {
        
        if (include.lowest) {
            prev <- rowSums(x >= detection)
        } else {
            prev <- rowSums(x > detection)
        }
    }
    
    if (!count) {
        prev <- prev/prevalence_nsamples(x)
    }
    
    if (sort) {
        prev <- rev(sort(prev))
    }
    
    # Return
    prev
    
}




# Internal auxiliary function
prevalence_nsamples <- function(x) {
    
    if (is.vector(x)) {
        n <- length(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
        n <- ncol(x)
    }
    
    n
    
}
