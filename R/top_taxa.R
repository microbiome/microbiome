#' @title Top Taxa
#' @description Return n most abundant taxa (based on total abundance over
#' all samples), sorted by abundance
#' @param x phyloseq object
#' @param n Number of top taxa to return (default: all)
#' @return Character vector listing the top taxa
#' @export
#' @examples
#' data(dietswap)
#' topx <- top_taxa(dietswap, n=10)
top_taxa <- function(x, n=ntaxa(x)) {
    
    names(sort(rowSums(abundances(x)), decreasing=TRUE)[seq_len(n)])
    
}
