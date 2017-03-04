#' @title Prevalence for Phyloseq OTUs
#' @description Simple prevalence measure.
#' @param x A vector, data matrix or phyloseq object
#' @param detection Detection threshold for absence/presence.
#' @param sort Sort the groups by prevalence
#' @param count Logical. Indicate prevalence as fraction of samples
#' (in percentage [0, 1]; default); or in absolute counts indicating
#' the number of samples where the OTU is detected above the given
#' abundance threshold.
#' @details For vectors, calculates the fraction (count = FALSE) or
#' number (count = TRUE) of samples that exceed the
#' detection. For matrices, calculates this for each matrix
#' column. For phyloseq objects, calculates this for each OTU. The
#' relative prevalence (count = FALSE) is simply the absolute
#' prevalence (count = TRUE) divided by the number of samples.
#' @return For each OTU, the fraction of samples where a given OTU is
#' detected. The output is readily given as a percentage.
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   data(peerj32)     
#'   ## With matrix
#'   prevalence(peerj32$data$microbes, detection = 200, sort = TRUE)
#'   ## With phyloseq
#'   prevalence(peerj32$phyloseq, detection = 200, sort = TRUE)
#'   prevalence(peerj32$phyloseq, detection = 200, sort = TRUE, count = TRUE)
prevalence <- function (x, detection = 0, sort = FALSE, count = FALSE) {

  if (is.null(x)) {
    warning("x is NULL - returning NULL")
    return(NULL)
  }

  if (is.vector(x)) {
    prev <- sum(x > detection)
  } else if (is.matrix(x) || is.data.frame(x)) {
    prev <- rowSums(x > detection)
  } else if (class(x) == "phyloseq") {
    # At this point necessary to have count = TRUE
    prev <- prevalence(abundances(x),
    	      detection = detection, count = TRUE)
  } 

  if (!count) {
    prev <- 1 * prev/prevalence_nsamples(x)
  }

  if (sort) {
    prev <- rev(sort(prev))
  }
  
  # Return
  prev

}




# Internal auxiliary function
prevalence_nsamples = function (x) {

  if (is.vector(x)) {
    n <- length(x)
  } else if (is.matrix(x) || is.data.frame(x)) {
    n <- ncol(x)    
  } else if (class(x) == "phyloseq") {
    n <- nsamples(x)
  }

  n

}
