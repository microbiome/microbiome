#' @title Dominant taxa
#' @description Returns the dominant taxonomic group for each sample.
#' @param x A \code{\link{phyloseq-class}} object
#' @param level Optional. Taxonomic level.
#' @return A vector of dominance indices
#' @export
#' @examples
#' data(dietswap)
#' # vector
#' d <- dominant(dietswap)
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
dominant <- function(x, level = NULL) {

    if (!is.null(level)) {
        x <- aggregate_taxa(x, level = level)
    }

    # TODO add support to other taxonomic levels
    taxa(x)[apply(abundances(x), 2, which.max)]

}
