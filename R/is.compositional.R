#' @title Check Compositionality 
#' @description Check if a phyloseq object is compositional.
#' @param x \code{\link{phyloseq-class}} object
#' @return Logical TRUE/FALSE
#' @examples
#' data(dietswap)
#' d <- is.compositional(dietswap)
#' d <- is.compositional(transform(dietswap, "compositional"))
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
is.compositional <- function(x) {

    # Test that all sample abundances sum to one with a tolerance of 1e-10
    max(abs(colSums(abundances(x)) - 1)) < 1e-10

}





