#' @title Neatmap Ordering
#' @description Order samples or taxa based on the neatmap approach.
#' @param x \code{\link{phyloseq-class}} object
#' @param target Order "sites" (samples) or "species" (taxa/OTUs)
#' @param method Ordination method. See \code{\link{ordinate}}
#'                 from \pkg{phyloseq} package.
#' @param distance Distance method. See \code{\link{ordinate}}
#'                 from \pkg{phyloseq} package.
#' @param first Optionally provide the name of the first sample/taxon to
#'        start the ordering (the ordering is cyclic so we can start at any
#'        point). The choice of the first sample may somewhat affect the
#'        overall ordering.
#' @param ... Arguments to be passed.
#' @return Vector of ordered elements
#' @export
#' @examples
#'   \dontrun{
#'     library(microbiome)
#'     data(peerj32)
#'     pseq <- peerj32$phyloseq
#'     order.sample <- order_neatmap(pseq, target = "sites",
#'                                 method = "NMDS", distance = "bray",
#'				   first = NULL) 
#'     order.otu <- order_neatmap(pseq, target = "species", method = "NMDS",
#'     	       	                distance = "bray", first = NULL)
#'   }  
#' @references This function is partially based on code derived from the
#'             \pkg{phyloseq} package. For the original
#'             neatmap approach for heatmap sorting, see (and cite):
#'   Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering heat map
#'   alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' @details This function borrows elements from the heatmap implementation in
#' the \pkg{phyloseq} package. The row/column sorting is there not available
#' as a separate function at present, however, hindering reuse in other tools.
#' Implemented in the microbiome package to provide an independent method for
#' easy sample/taxon reordering for phyloseq objects.
#' @keywords utilities
order_neatmap <- function (x, target, method = "NMDS", distance = "bray",
	      	           first = NULL, ...) {

  # Capture the output to keep the screen clean
  junk <- capture.output(
    ord <- ordinate(x, method, distance, ...), file=NULL)

  # Order items with the  NeatMap procedure
  # Reorder by the angle in radial coordinates on the 2-axis plane.
  DF <- NULL

  # Define new sample ordering based on the ordination
  tmp <- try({
    DF <- scores(ord, choices = c(1, 2), display = target)}, silent = TRUE)

  if(inherits(tmp, "try-error")){
    warning(paste("Order failed with ", target, ". 
    			 Using default ordering.", sep = ""))
  }


  if(!is.null(DF)){
    # If the score accession worked, replace order
    if (target == "sites") {
      ordering <- sample_names(x)[order(radial_theta(DF))] 
    } else if (target == "species") {
      ordering <- taxa(x)[order(radial_theta(DF))] 
    } else {
      stop("Target should be either sites or species")
    }
  } else if (length(DF) > 1) {

    if (target == "sites") {
      ordering <- sample_names(x)[order(DF)] # 1:nsamples(x)
    } else if (target == "species") {
      ordering <- taxa(x)[order(DF)] # 1:ntaxa(x)
    } else {
      stop("Target should be either sites or species")
    }
  } else {
    if (target == "sites") {
      ordering <- sample_names(x)
    } else if (target == "species") {
      ordering <- taxa(x)
    } else {
      stop("Target should be either sites or species")
    }
  }

  # Determine the starting item (OTU or sample)
  if( !is.null(first) ){
    ordering <- chunk_reorder(ordering, first)
  }

  ordering

}


