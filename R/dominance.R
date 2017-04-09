#' @title Dominance Index
#' @description Calculates the community dominance index.
#' @param rank Optional. The rank of the dominant taxa to consider.
#' @param aggregate Aggregate (TRUE; default) the top members or not. If aggregate=TRUE, then the sum of relative abundances is returned. Otherwise the relative abundance is returned for the single taxa with the indicated rank. 
#' @inheritParams diversity
#' @return A vector of dominance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- dominance(dietswap)
#' @details The dominance index gives the relative abundance of the most abundant species in [0,1]. This simple diversity index is occasionally used in ecological literature. With rank = 2, the sum of abundances for the two most abundant taxa are returned etc. However, if aggregate=FALSE, the abundance for the single n'th most dominant taxa (n = rank) is returned instead the sum of abundances up to that rank.
#' @references
#'   Dominance has been used in this sense for instance in
#'   Kenneth J. Locey and Jay T. Lennon. Scaling laws predict global microbial diversity. PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso coverage, core_abundance, rarity, global
#' @keywords utilities
dominance <- function(x, rank = 1, aggregate = TRUE, split = TRUE) {

  # Pick the OTU data
  otu <- abundances(x)

  if (!split) {
    otu <- as.matrix(rowSums(otu), nrow = nrow(otu))
  }

  if (!aggregate) {
    do <- apply(otu, 2, function (x) {rev(sort(x/sum(x, na.rm = TRUE)))[[rank]]})
  } else {
    do <- apply(otu, 2, function (x) {sum(rev(sort(x/sum(x, na.rm = TRUE)))[1:rank])})
  }
  names(do) <- sample_names(x)
  
  do

}


