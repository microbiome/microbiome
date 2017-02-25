#' @title Core Matrix 
#' @description Creates the core matrix.
#' @param x \code{\link{phyloseq}} object or a taxa x samples abundance matrix
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range
#' @return Estimated core microbiota
#' @examples
#'   library(microbiome)
#'   data(peerj32)
#'   core <- core_matrix(peerj32$phyloseq)
#' @export 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_matrix <- function(x,  
          prevalence.intervals = seq(5, 100, 5), 
          detection.thresholds = NULL) {

    data <- x	  
    if (class(x) == "phyloseq") {
      # Pick abundance matrix
      data <- abundances(x)
    }

    # Convert prevalences from percentages to sample counts
    p.seq <- 0.01 * prevalence.intervals * ncol(data)

    ## Intensity vector
    if (is.null(detection.thresholds)) {
      detection.thresholds <- seq(min(data), max(data), length = 10)
    }
    i.seq <- detection.thresholds

    coreMat <- matrix(NA, nrow = length(i.seq), ncol = length(p.seq), 
                      	  dimnames = list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
      for (p in p.seq) { 
        # Number of OTUs above a given prevalence threshold     
        coreMat[as.character(i), as.character(p)] <- sum(rowSums(data > i)>= p)
      }
    }
    
    # Convert Prevalences to percentages
    colnames(coreMat) <- as.numeric(colnames(coreMat))/ncol(data)
    rownames(coreMat) <- as.character(as.numeric(rownames(coreMat)))
    
    coreMat

}

