#' @title Gini Index
#' @description Calculates the Gini index for phyloseq object
#' @inheritParams diversity
#' @return A vector of Gini indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- gini(dietswap)
#' @references
#'     Relative Distribution Methods in the Social Sciences. Mark S.
#'     Handcock and Martina Morris, Springer-Verlag, Inc., New York,
#'     1999. ISBN 0387987789.
#' @seealso diversity, reldist::gini (inspired by that implementation but independently written here to avoid external depedencies)
#' @details Gini index is a common measure for inequality in economical income, but can also be used as a community diversity measure.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
gini <- function(x, detection, prevalence, split = TRUE) {


  # Pick the OTU data
  otu <- abundances(x)

  if (!split) {
    otu <- as.matrix(rowSums(otu), nrow = nrow(otu))
  }

  # Gini index for each sample
  do <- apply(otu, 2, function (x) {gini_index(x)})
  names(do) <- sample_names(x)
  
  do

}



gini_index <- function (x, w = rep(1, length(x))) {
  # See also reldist::gini for an independent implementation
    o  <- order(x)
    x  <- x[o]
    w  <- w[o]/sum(w)
    p  <- cumsum(w)
    nu <- cumsum(w*x)
    n  <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}




