#' Description: Simple prevalence measure
#'
#' Arguments:
#'   @param x Abundance data matrix: samples x features (microbes) 
#'   @param detection.threshold Detection threshold for absence/presence.
#'
#' @details Calculates for each sample the fraction of samples that
#'          exceed the detection threshold
#'
#'
#' Returns:
#'   @return Number of OTUs.
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
#' data(peerj32)
#' prevalence(peerj32$microbes, detection.threshold, sort = T)
prevalence <- function (x, detection.threshold, sort = F) {
  prev <- colMeans(x > detection.threshold)
  if (sort) {prev <- rev(sort(prev))}
  prev
}



#' Description: List prevalent groups
#'
#' Arguments:
#'   @param x Abundance data matrix: samples x features (microbes) 
#'   @param detection.threshold Detection threshold for absence/presence.
#'   @param prevalence.threshold Detection threshold for prevalence
#'
#' @details Lists groups that are more prevalent above the detection
#'          threshold as specified by the detection and prevalence threshold
#' 	    arguments
#'
#' Returns:
#'   @return Vector of prevalent groups
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
#' data(peerj32)
#' list_prevalent_groups(peerj32$microbes, 2, 0.2)
#' 
list_prevalent_groups <- function (x, detection.threshold, prevalence.threshold) {
  sort(names(which(prevalence(x, detection.threshold) > prevalence.threshold)))
}


