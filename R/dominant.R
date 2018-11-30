#' @title Dominant taxa
#' @description Returns the dominant taxonomic group for each sample.
#' @param x A \code{\link{phyloseq-class}} object
#' @return A vector of dominance indices
#' @export
#' @examples
#' data(dietswap)
#' # vector
#' d <- dominant(dietswap)
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
dominant <- function(x) {

    # TODO add support to other taxonomic levels
    taxa(x)[apply(abundances(x), 2, which.max)]

}
