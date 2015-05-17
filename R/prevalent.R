#' Simple prevalence measure
#'
#' 
#'   @param x A vector, data matrix or phyloseq object
#'   @param detection.threshold Detection threshold for absence/presence.
#'   @param sort Sort the groups by prevalence
#'
#'
#' @details For vectors, calculates the fraction of samples that exceed the
#' detection threshold. For matrices, calculates for each matrix column
#' the fraction of entries that #' exceed the detection threshold. For
#' phyloseq object, calculates for each OTU the fraction of samples that
#' exceed the detection threshold
#'
#' @return For each OTU, the fraction of samples where a given OTU is detected
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   peerj32 <- download_microbiome("peerj32")
#'   # With matrix
#'   prevalence(peerj32$data$microbes, detection.threshold = 200, sort = TRUE)
#'   # With phyloseq
#'   prevalence(peerj32$physeq, detection.threshold = 200, sort = TRUE)
#' 
prevalence <- function (x, detection.threshold, sort = FALSE) {

  if (is.vector(x)) {
    prev <- mean(x > detection.threshold)
  } else if (is.matrix(x)) {
    prev <- colMeans(x > detection.threshold)
  } else if (class(x) == "phyloseq") {
    x <- t(otu_table(x)@.Data)
    prev <- prevalence(x, detection.threshold = detection.threshold)
  }

  if (sort) {
    prev <- rev(sort(prev))
  }

  prev

}



#' List prevalent groups
#'
#'   @param x A matrix or a x \code{\link{phyloseq}} object
#'   @param detection.threshold Detection threshold for absence/presence.
#'   @param prevalence.threshold Detection threshold for prevalence
#'
#' @details For phyloseq object, lists taxa that are more prevalent with the given detection
#'          threshold. For matrix, lists columns that satisfy these criteria.
#'
#'   @return Vector of prevalent taxa names
#'
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   peerj32 <- download_microbiome("peerj32")
#'   prevalent_taxa(peerj32$data$microbes, 10^1.8 + 100, 0.2) # matrix
#'   prevalent_taxa(peerj32$physeq, 100, 0.2) # phyloseq object
#' 
prevalent_taxa <- function (x, detection.threshold, prevalence.threshold) {

  if (class(x) == "phyloseq") {
    # Convert into OTU matrix
    x <- t(otu_table(x)@.Data)    
  } 

  sort(names(which(prevalence(x, detection.threshold) > prevalence.threshold)))

}


