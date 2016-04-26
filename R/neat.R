#' @title Neatmap sorting for matrices
#' @description Order matrix based on the neatmap approach.
#' @param x A matrix.
#' @param arrange Order "rows" or "cols" or "both".
#' @param method Ordination method. Only NMDS implemented for now.
#' @param distance Distance method. See \code{\link{vegdist}} function from the \pkg{vegan} package.
#' @param first.row Optionally provide the name of the first row to start the ordering (the ordering is cyclic so we can start at any point). The choice of the first sample may somewhat affect the overall ordering.
#' @param first.col Optionally provide the name of the first col to start the ordering (the ordering is cyclic so we can start at any point). The choice of the first sample may somewhat affect the overall ordering.
#' @param ... Arguments to pass.
#' @return Sorted matrix
#' @export
#' @examples 
#'   data(peerj32)
#'   x <- peerj32$microbes
#'   xo <- neat(x, "both", method = "NMDS", distance = "bray") 
#' @references This function is partially based on code derived from the \pkg{phyloseq} package. However for the original
#'   neatmap approach for heatmap sorting, see (and cite):
#'   Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details Borrows elements from the heatmap implementation in the \pkg{phyloseq} package. The row/column sorting is not available there as a separate function. Therefore I implemented this function to provide an independent method for easy sample/taxon reordering for phyloseq objects.
#' @keywords utilities
neat <- function (x, arrange = "both", method = "NMDS", distance = "bray", first.row = NULL, first.col = NULL, ...) {

  if (arrange %in% c("rows", "both")) {
    if (nrow(x) > 2) {
      sr <- neatsort(x, "rows", method = method, distance = distance, first = first.row, ...)
      x <- x[sr,]
    } 

  }
  if (arrange %in% c("cols", "both")) {
    if (ncol(x) > 2) {
      sc <- neatsort(x, "cols", method = method, distance = distance, first = first.col, ...)
      x <- x[,sc]
    }
  }

  x
}
