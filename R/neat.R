#' @title Neatmap Sorting
#' @description Order matrix or phyloseq OTU table based on the neatmap
#' approach.
#' @param x A matrix or phyloseq object.
#' @param arrange Order 'features', 'samples' or 'both' (for matrices).
#' For matrices, it is assumed that the samples are on the columns and
#' features are on the rows. For phyloseq objects, features are the taxa of
#' the OTU table. 
#' @param method Ordination method. Only NMDS implemented for now.
#' @param distance Distance method. See \code{\link{vegdist}} function from
#' the \pkg{vegan} package.
#' @param first.feature Optionally provide the name of the first feature to
#' start the ordering 
#' @param first.sample Optionally provide the name of the first sample to
#' start the ordering 
#' @param ... Arguments to pass.
#' @return Sorted matrix
#' @export
#' @examples
#' data(peerj32)
#' # Take subset to speed up example
#' x <- peerj32$microbes[1:10,1:10]
#' xo <- neat(x, 'both', method='NMDS', distance='bray')
#' 
#' @references This function is partially based on code derived from the
#' \pkg{phyloseq} package. However for the original
#' neatmap approach for heatmap sorting, see (and cite):
#' Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map
#' alternatives in R. BMC Bioinformatics, 11, 45.
#' @details Borrows elements from the heatmap implementation in the
#' \pkg{phyloseq} package. The row/column sorting is not available there
#' as a separate function. Therefore I implemented this function to
#' provide an independent method for easy sample/taxon reordering for
#' phyloseq objects. The ordering is cyclic so we can start at any
#' point. The choice of the first sample may somewhat affect the overall
#' ordering
#' @keywords utilities
neat <- function(x, arrange="both", method="NMDS", distance="bray",
    first.feature=NULL, first.sample=NULL, ...) {
    
    # Ensure data is in matrix form
    x <- abundances(x)
    
    if (is.null(rownames(x))) {
        rownames(x) <- as.character(seq_len(nrow(x)))
    }
    
    if (is.null(colnames(x))) {
        colnames(x) <- as.character(seq_len(ncol(x)))
    }
    
    if (arrange %in% c("features", "both")) {
        if (nrow(x) > 2) {
            sr <- neatsort(x, "features", method=method,
                    distance=distance, first=first.feature, 
                ...)
            x <- x[sr, ]
        }
        
    }
    if (arrange %in% c("samples", "both")) {
        if (ncol(x) > 2) {
            sc <- neatsort(x, "samples", method=method,
            distance=distance, first=first.sample, 
                ...)
            x <- x[, sc]
        }
    }
    
    x
}
