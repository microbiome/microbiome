#' @title Rare Taxa
#' @description Determine members of the rare microbiota with given abundance
#''             and prevalences.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold (non-negative real)
#' @param prevalence Prevalence threshold (in [0, 100])
#' @param sort Logical. Sort the taxa.
#' @return Vector of rare members
#' @details For phyloseq object, lists taxa that are more prevalent with the
#'   given detection. For matrix, lists columns that satisfy
#'   these criteria.
#' @examples
#'   data(dietswap)
#'   a <- rare_members(dietswap, 1, 95)
#' @export
#' @references 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rare_members <- function(x, detection = 1, prevalence = 95,
                    sort = TRUE)  {

  # First identify core based on these criteria
  cm <- core_members(x, detection, prevalence, sort)

  # Then pick rare taxa as those that do not belong to core
  # ie have prevalence below the given threshold at the given
  # detection limit
  taxa <- setdiff(taxa(x), cm)

  if (sort) {
    taxa <- sort(taxa)
  }

  taxa

}
