#' @title Test Compositionality
#' @description Test if phyloseq object is compositional.
#' @param x \code{\link{phyloseq-class}} object
#' @param tolerance Tolerance for detecting compositionality.
#' @return Logical TRUE/FALSE
#' @keywords utilities
#' @details This function tests that the sum of abundances within each sample
#' is almost zero, within the tolerance of 1e-6 by default.
#' @export
#' @seealso transform
#' @examples
#' data(dietswap)
#' a <- is_compositional(dietswap)
#' b <- is_compositional(transform(dietswap, "identity"))
#' c <- is_compositional(transform(dietswap, "compositional"))
is_compositional <- function(x, tolerance = 1e-6) {
    # Test that the data is sufficiently close to compositional        
    all(abs(colSums(abundances(x)) - 1) < tolerance)
}






