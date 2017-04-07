#' @title Dominance Index
#' @description Calculates the community dominance index.
#' @param threshold Indicates the fraction of the ecosystem that has to be covered by the dominant species.
#' @inheritParams diversity
#' @return A vector of dominance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- dominance(dietswap, threshold = 0.5)
#' @details The dominance index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso top_abundance, diversity
#' @keywords utilities
dominance <- function(x, threshold = 0.5, split = TRUE) {

  # Pick the OTU data
  otu <- abundances(x)

  if (!split) {
    otu <- as.matrix(rowSums(otu)/sum(otu), nrow = nrow(otu))
  }

  # Number of groups needed to have 50% of the ecosystem occupied
  do <- apply(otu, 2, function (x) {min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))})

  names(do) <- sample_names(x)
  
  do

}


