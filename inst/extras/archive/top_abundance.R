#' @title Top Abundance Index
#' @description Calculates the community top_abundance index.
#' @param rank Optional. The rank of the dominant taxa to consider.
#' @param aggregate Aggregate (TRUE; default) the top members or not. If aggregate=TRUE, then the sum of relative abundances is returned. Otherwise the relative abundance is returned for the single taxa with the indicated rank. 
#' @inheritParams diversity
#' @return A vector of top_abundance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- top_abundance(dietswap)
#' @details The top_abundance index gives the abundance of the most abundant species in [0,1]. This simple diversity index is occasionally used in ecological literature, and sometimes also called dominance. However, note that the microbiome::dominance function uses a different definition. With rank = 2, the sum of abundances for the two most abundant taxa are returned etc. However, if aggregate=FALSE, the abundance for the single n'th most dominant taxa (n = rank) is returnde instead the sum of abundances up to that rank.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso dominance, diversity
#' @keywords utilities
top_abundance <- function(x, rank = 1, aggregate = TRUE, split = TRUE) {

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


