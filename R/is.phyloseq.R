#' @title Identify Phyloseq Objects
#' @description Identifies whether a given object is from the phyloseq class
#' @param x object to test
#' @return Logical
#' @export
#' @examples
#' data(dietswap)
#' is.phyloseq(dietswap)
is.phyloseq <- function(x) {
    class(x) == "phyloseq"
}
