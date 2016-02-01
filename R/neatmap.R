#' @title Neatmap ordering
#' @description Order samples or taxa based on the neatmap approach.
#' @param x \code{\link{phyloseq-class}} object
#' @param target Order "sites" (samples) or "species" (taxa/OTUs)
#' @param method Ordination method. See \code{\link{ordinate}} from \pkg{phyloseq} package.
#' @param distance Distance method. See \code{\link{ordinate}} from \pkg{phyloseq} package.
#' @param first Optionally provide the name of the first sample/taxon to start the ordering (the ordering is cyclic so we can start at any point). The choice of the first sample may somewhat affect the overall ordering.
#' @param ... Arguments to be passed.
#' @return Vector of ordered elements
#' @export
#' @examples \dontrun{
#'    data(peerj32)
#'    pseq <- peerj32$phyloseq
#'    order.sample <- order_neatmap(pseq, target = "sites", method = "NMDS", distance = "bray", first = NULL) 
#'    order.otu <- order_neatmap(pseq, target = "species", method = "NMDS", distance = "bray", first = NULL)
#'                   }
#' @references This function is partially based on code derived from the \pkg{phyloseq} package. However for the original
#'   neatmap approach for heatmap sorting, see (and cite):
#'   Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details This function borrows elements from the heatmap implementation in the \pkg{phyloseq} package. The row/column sorting is there
#' not available as a separate function at present, however, hindering reuse in other tools. Therefore I implemented this function to
#' provide an independent method for easy sample/taxon reordering for phyloseq objects.
#' @importFrom vegan scores
#' @importFrom phyloseq ordinate
#' @keywords utilities
order_neatmap <- function (x, target, method = "NMDS", distance = "bray", first = NULL, ...) {

  # Capture the output to keep the screen clean
  junk <- capture.output(
    ord <- ordinate(x, method, distance, ...), file=NULL)

  # Order items with the  NeatMap procedure
  # Reorder by the angle in radial coordinates on the 2-axis plane.
  DF <- NULL

  # Define new sample ordering based on the ordination
  trash <- try({
    DF <- scores(ord, choices = c(1, 2), display = target)}, silent = TRUE)
    #Was : DF <- scores(ord, choices = c(1, 2), display = target)}, silent = TRUE, x = x)

  if(inherits(trash, "try-error")){
    warning(paste("Ordering failed for ", target, ". Using default ordering.", sep = ""))
  }
  if(!is.null(DF)){
    # If the score accession worked, replace order
    if (target == "sites") {
      ordering <- sample_names(x)[order(RadialTheta(DF))] 
    } else if (target == "species") {
      ordering <- taxa_names(x)[order(RadialTheta(DF))] 
    }
  }
  # Determine the starting item (OTU or sample)
  if( !is.null(first) ){
    ordering <- chunkReOrder(ordering, first)
  }

  ordering

}




#' @title Neatmap ordering for matrices
#' @description Order rows or columns based on the neatmap approach.
#' @param x A matrix.
#' @param arrange Order "rows" or "cols" or "both".
#' @param method Ordination method. Only NMDS implemented for now.
#' @param distance Distance method. See \code{\link{vegdist}} function from the \pkg{vegan} package.
#' @param first Optionally provide the name of the first sample/taxon to start the ordering (the ordering is cyclic so we can start at any point). The choice of the first sample may somewhat affect the overall ordering.
#' @param ... Arguments to pass.
#' @return Vector of ordered elements
#' @export
#' @examples \dontrun{
#'    data(peerj32)
#'    x <- peerj32$microbes
#'    x <- neat(x, "both", method = "NMDS", distance = "bray") 
#'                   }
#' @references This function is partially based on code derived from the \pkg{phyloseq} package. However for the original
#'   neatmap approach for heatmap sorting, see (and cite):
#'   Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details This function borrows elements from the heatmap implementation in the \pkg{phyloseq} package. The row/column sorting is there
#' not available as a separate function at present, however, hindering reuse in other tools. Therefore I implemented this function to
#' provide an independent method for easy sample/taxon reordering for phyloseq objects.
#' @importFrom vegan scores
#' @importFrom vegan vegdist
#' @importFrom vegan metaMDS
#' @keywords utilities
neat <- function (x, arrange = "both", method = "NMDS", distance = "bray", first = NULL, ...) {

  if (arrange == "both") {
    x <- neat(x, "rows", method = method, distance = distance, first = first, ...)
    x <- neat(x, "cols", method = method, distance = distance, first = first, ...)    
    return(x)
  }


  if (arrange == "cols") {
    x <- t(x)
  }

  # neatmap sorting for matrices with NMDS  
  # Distance
  d <- vegdist(x, distance)

  # Order		     
  # Capture the output to keep the screen clean
  junk <- capture.output(
    ord <- metaMDS(d, wascores = FALSE, autotransform = FALSE, noshare = FALSE), file=NULL)

  # Order items with the  NeatMap procedure
  # Reorder by the angle in radial coordinates on the 2-axis plane.
  DF <- NULL

  # Define new sample ordering based on the ordination
  trash <- try({
    DF <- scores(ord, choices = c(1, 2), display = "sites")}, silent = TRUE)

  if(inherits(trash, "try-error")){
    warning(paste("Ordering failed. Using default ordering.", sep = ""))
  }
  if(!is.null(DF)){
    # If the score accession worked, replace order
    ordering <- rownames(x)[order(RadialTheta(DF))] 
  }
  # Determine the starting item (OTU or sample)
  if( !is.null(first) ){
    ordering <- chunkReOrder(ordering, first)
  }

  x <- x[ordering,]

  # Return to the original shape
  if (arrange == "cols") {
    x <- t(x)
  }

  x
}






#' @title Chunk reorder
#' @description Chunk re-order a vector so that specified newstart is first. Different than relevel.
#' @keywords internal
#' @details Borrowed from \pkg{phyloseq} package as needed here and not exported there.
#' @examples 
#' # Typical use-case
#' # chunkReOrder(1:10, 5)
#' # # Default is to not modify the vector
#' # chunkReOrder(1:10)
#' # # Another example not starting at 1
#' # chunkReOrder(10:25, 22)
#' # # Should silently ignore the second element of `newstart`
#' # chunkReOrder(10:25, c(22, 11))
#' # # Should be able to handle `newstart` being the first argument already
#' # # without duplicating the first element at the end of `x`
#' # chunkReOrder(10:25, 10)
#' # all(chunkReOrder(10:25, 10) == 10:25)
#' # # This is also the default
#' # all(chunkReOrder(10:25) == 10:25)
#' # # An example with characters
#' # chunkReOrder(LETTERS, "G") 
#' # chunkReOrder(LETTERS, "B") 
#' # chunkReOrder(LETTERS, "Z") 
#' # # What about when `newstart` is not in `x`? Return x as-is, throw warning.
#' # chunkReOrder(LETTERS, "g") 
chunkReOrder <- function(x, newstart = x[[1]]){
  pivot <- match(newstart[1], x, nomatch = NA)
  # If pivot `is.na`, throw warning, return x as-is
  if(is.na(pivot)){
    warning("The `newstart` argument was not in `x`. Returning `x` without reordering.")
    newx <- x
  } else {
    newx <- c(tail(x, {length(x) - pivot + 1}), head(x, pivot - 1L))
  }
  return(newx)
}
