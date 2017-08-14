#' @title Neatmap Ordering for Matrix Features or Samples
#' @description Order features or samples based on the neatmap approach.
#' @param x A matrix.
#' @param arrange Order "features" or "samples". For matrices, it is assumed that the samples are on the columns and features are on the rows. For phyloseq objects, features are the taxa of the OTU table.
#' @param method Ordination method. Only NMDS implemented for now.
#' @param distance Distance method. See \code{\link{vegdist}} function from the \pkg{vegan} package.
#' @param first Optionally provide the name of the first sample (or feature) to start the ordering (the ordering is cyclic so we can start at any point). The choice of the first sample/feature may somewhat affect the overall ordering.
#' @param ... Arguments to pass.
#' @return Vector of ordered elements
#' @export
#' @examples \dontrun{
#'   data(peerj32)
#'   x <- peerj32$microbes
#'   x <- neatsort(x, "features", method = "NMDS", distance = "bray")
#'                   }
#' @references This function is partially based on code derived from the \pkg{phyloseq} package. However for the original
#'   neatmap approach for heatmap sorting, see (and cite):
#'   Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details This function borrows elements from the heatmap implementation in the \pkg{phyloseq} package. The row/column sorting is there
#' not available as a separate function at present, however, hindering reuse in other tools. This function provides an independent
#' method for easy row/column reordering for matrices. This a quick hack and using the ordination could be expanded further
#' (now only NMDS is available, and the sorting is done independently for features and samples).
#' @keywords utilities
neatsort <- function (x, arrange, method = "NMDS", distance = "bray", first = NULL, ...) {

  if (arrange == "samples") {
    x <- t(x)
  }
  
  # Neatmap sorting for matrices with NMDS 
  d <- vegdist(x, distance)

  # Order		     
  # Capture the output to keep the screen clean
  junk <- capture.output(
    ord <- metaMDS(d, wascores = FALSE, autotransform = FALSE, noshare = FALSE), file=NULL)


  # Order items with the  NeatMap procedure
  # Reorder by the angle in radial coordinates on the 2-axis plane.
  DF <- NULL

  # Define new sample ordering based on the ordination
  tmp <- try({
    DF <- scores(ord, choices = c(1, 2), display = "sites")}, silent = TRUE)

  if(inherits(tmp, "try-error")){
    warning(paste("Order failed. Using default ordering.", sep = ""))
  }
  
  if(!is.null(DF)){
    # If the score accession worked, replace order
    ordering <- rownames(x)[order(radial_theta(DF))] 
  }
  # Determine the starting item (OTU or sample)
  if( !is.null(first) ){
    ordering <- chunk_reorder(ordering, first)
  }

  ordering
  
}



