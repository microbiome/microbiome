#' @title Coverage Index
#' @description Community coverage index.
#' @details The coverage index gives the number of groups needed to
#'   have a given proportion of the ecosystem occupied (by default 0.5 ie 50%). 
#' @param threshold Indicates the fraction  of the ecosystem to be occupied by the N most abundant species (N is returned by this function).
#' @inheritParams global
#' @return A vector of coverage indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- coverage(dietswap, threshold = 0.5)
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso dominance, global
#' @keywords utilities
coverage <- function(x, threshold = 0.5) {

  otu <- pick_data(x, compositional = FALSE)

  # Number of groups needed to have 50% of the ecosystem occupied
  do <- apply(otu, 2, function (x) {min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))})

  names(do) <- sample_names(x)
  
  do

}
