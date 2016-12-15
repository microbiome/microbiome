#' @title Prevalent taxa
#' @description List prevalent taxa.
#' @param x A matrix or a x \code{\link{phyloseq}} object
#' @param detection.threshold Detection threshold for absence/presence.
#' @param prevalence.threshold Detection threshold for prevalence,
#'        provided as percentages [0, 100]
#' @param sort Logical. Sort the taxa.
#' @details For phyloseq object, lists taxa that are more prevalent with the
#'   given detection threshold. For matrix, lists columns that satisfy
#'   these criteria.
#' @return Vector of prevalent taxa names
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
#'   prevalent_taxa(peerj32$data$microbes, 10^1.8 + 100, 0.2) # for matrix
#'   prevalent_taxa(peerj32$physeq, 100, 0.2) # for phyloseq object
prevalent_taxa <- function (x, detection.threshold, prevalence.threshold, sort = TRUE) {

  if (class(x) == "phyloseq") {
    x <- taxa_abundances(x)
  } 
  taxa <- names(which(prevalence(x, detection.threshold) > prevalence.threshold))

  if (sort) {
    taxa <- sort(taxa)
  }

  taxa

}


