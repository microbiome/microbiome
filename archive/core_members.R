#' @title Core Taxa
#' @description Determine members of the core microbiota with given abundance
#''             and prevalences.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold (non-negative real)
#' @param prevalence Prevalence threshold (in [0, 100])
#' @param sort Logical. Sort the taxa.
#' @return Vector of core members
#' @details For phyloseq object, lists taxa that are more prevalent with the
#'   given detection. For matrix, lists columns that satisfy
#'   these criteria.
#' @examples
#'   data(dietswap)
#'   a <- core_members(dietswap, 1, 95)
#' @export
#' @aliases prevalent_taxa
#' @references 
#'   A Salonen et al. The adult intestinal core microbiota is determined by 
#'   analysis depth and health status. Clinical Microbiology and Infection 
#'   18(S4):16 20, 2012. 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_members <- function(x, detection = 1, prevalence = 95,
                    sort = TRUE)  {

  if (class(x) == "phyloseq") {
    x <- abundances(x)
  } 
  taxa <- names(which(prevalence(x, detection) > prevalence))

  if (sort) {
    taxa <- sort(taxa)
  }

  taxa

}
