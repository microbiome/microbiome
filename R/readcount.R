#' @title Total Read Count
#' @description Total Read Count
#' @param x \code{\link{phyloseq-class}} object
#' @return Vector of read counts.
#' @examples
#' data(dietswap)
#' d <- readcount(dietswap)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
readcount <- function(x) {

    colSums(abundances(x))
    
}





