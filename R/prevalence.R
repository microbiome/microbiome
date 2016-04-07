#' @title Prevalence for phyloseq OTUs
#' @description Simple prevalence measure
#' @param x A vector, data matrix or phyloseq object
#' @param detection.threshold Detection threshold for absence/presence.
#' @param sort Sort the groups by prevalence
#' @details For vectors, calculates the fraction of samples that exceed the
#' detection threshold. For matrices, calculates for each matrix column
#' the fraction of entries that #' exceed the detection threshold. For
#' phyloseq object, calculates for each OTU the fraction of samples that
#' exceed the detection threshold
#' @return For each OTU, the fraction of samples where a given OTU is detected. The output is readily given as a percentage.
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   #peerj32 <- download_microbiome("peerj32")
#'   ## With matrix
#'   #prevalence(peerj32$data$microbes, detection.threshold = 200, sort = TRUE)
#'   ## With phyloseq
#'   #prevalence(peerj32$physeq, detection.threshold = 200, sort = TRUE)
prevalence <- function (x, detection.threshold, sort = FALSE) {

  if (is.vector(x)) {
    prev <- 100 * mean(x > detection.threshold)
  } else if (is.matrix(x)) {
    prev <- 100 * colMeans(x > detection.threshold)
  } else if (class(x) == "phyloseq") {
    x <- otu_table(x)@.Data
    prev <- prevalence(x, detection.threshold = detection.threshold)
  }

  if (sort) {
    prev <- rev(sort(prev))
  }

  prev

}



