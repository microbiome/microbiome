#' @title Dominance Index
#' @description Calculates the community dominance index.
#' @param index If the index is given, it will override the other parameters. See the details below for description and references of the standard dominance indices. By default, this function returns the Berger-Parker index, ie relative dominance at rank 1.
#' @param rank Optional. The rank of the dominant taxa to consider.
#' @param relative Use relative abundances (default: TRUE)
#' @param aggregate Aggregate (TRUE; default) the top members or not. If aggregate=TRUE, then the sum of relative abundances is returned. Otherwise the relative abundance is returned for the single taxa with the indicated rank. 
#' @inheritParams diversity
#' @return A vector of dominance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- dominance(dietswap, rank = 1, relative = TRUE)
#' @details The dominance index gives the abundance of the most abundant species and has been used in microbiomics context for instance in Locey & Lennon (2016). The following indices are provided: 1) "absolute_dominance" is the most simple variant, giving the absolute abundance of the most abundant species (Magurran & McGill 2011). By default, this refers to the single most dominant species (rank = 1) but it is possible to calculate the absolute dominance with rank n based on the abundances of top-n species by tuning the rank argument. 2) "relative_dominance" gives the relative abundance of the most abundant species. This is with rank = 1 by default but can be calculated for other ranks. 3) "DBP" is the Berger–Parker index, a special case of relative dominance with rank 1; 4) "DMN" is the McNaughton’s dominance. This is the sum of the relative abundance of the two most abundant taxa, or a special case of relative dominance with rank 2; 5) "Simpson" Simpson's index has also an interpretation as a dominance measure. Finally, it is also possible to calculated dominances up to an arbitrary rank by setting the rank argument. Finally, by setting aggregate=FALSE, the abundance for the single n'th most dominant taxa (n = rank) is returned instead the sum of abundances up to that rank (the default). 
#' @references
#'
#'   Kenneth J. Locey and Jay T. Lennon. Scaling laws predict global microbial diversity. PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#'   Magurran AE, McGill BJ, eds (2011) Biological Diversity: Frontiers in Measurement and Assessment (Oxford Univ Press, Oxford), Vol 12
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso coverage, core_abundance, rarity, global
#' @keywords utilities
dominance <- function(x, index = NULL, rank = 1, relative = TRUE, aggregate = TRUE, split = TRUE) {

  if (is.null(index)) {
    rank <- rank
  } else if (index == "absolute_dominance") {
    relative <- FALSE # Rank = 1 by default but can be tuned
  } else if (index %in% c("relative_dominance")) {
    relative <- TRUE # Rank = 1 by default but can be tuned
  } else if (index %in% c("DBP")) {
    # Berger-Parker
    rank <- 1
    relative <- TRUE
  } else if (index %in% c("DMN")) {
    # McNaughton's dominance
    rank <- 2
    relative <- TRUE
    aggregate <- TRUE
  } else if (index %in% c("Simpson")) {
    ret <- estimate_richness(x, measures = "Simpson")$Simpson
    return(ret)
  } 

  if (relative) {
    x <- transform(x, "compositional")
  }

  # Pick the OTU data
  otu <- abundances(x)

  if (!split) {
    otu <- as.matrix(rowSums(otu), nrow = nrow(otu))
  }

  if (!aggregate) {
    do <- apply(otu, 2, function (x) {rev(sort(x))[[rank]]})
  } else {
    do <- apply(otu, 2, function (x) {sum(rev(sort(x))[1:rank])})
  }
  names(do) <- sample_names(x)
  
  do

}


