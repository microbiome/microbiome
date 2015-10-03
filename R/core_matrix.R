#' core_matrix
#'
#' create core matrix 
#'
#' @param x \code{\link{phyloseq}} object
#' @param verbose verbose
#' @param prevalence.intervals a vector of prevalence percentages in [0,100]
#' @param detection.thresholds a vector of intensities around the data range
#'
#' @return Estimated core microbiota
#'
#' @examples 
#'   #pseq <- download_microbiome("peerj32")$physeq
#'   #core <- core_matrix(pseq)
#'
#' @export 
#' 
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_matrix <- function(x,  
          prevalence.intervals = seq(5, 100, 5), 
          detection.thresholds = NULL, verbose = FALSE) {
    
    # Convert into OTU matrix
    data <- otu_table(x)@.Data

    # Use log10 
    data <- log10(data)

    # Convert prevalences from percentages to numerics
    p.seq <- 0.01 * prevalence.intervals * ncol(data)

    ## Intensity vector
    if (is.null(detection.thresholds)) {
      i.seq <- seq(min(data), max(data), length = 10)
    } else {   
      # Use log10
      i.seq <- log10(detection.thresholds)
    }
    
    coreMat <- matrix(NA, nrow = length(i.seq), ncol = length(p.seq), 
                      	  dimnames = list(i.seq, p.seq))
    
    n <- length(i.seq) * length(p.seq)
    cnt <- 0
    for (i in i.seq) {
      for (p in p.seq) {
        if (verbose) {
          cnt <- cnt + 1
        }
          coreMat[as.character(i), as.character(p)] <- core.sum(data, i, p)
        }
    }
    
    # Convert Prevalences to percentages
    colnames(coreMat) <- 100 * as.numeric(colnames(coreMat))/ncol(data)
    rownames(coreMat) <- as.character(10^as.numeric(rownames(coreMat)))
    
    coreMat

}

