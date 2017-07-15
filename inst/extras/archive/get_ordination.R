#' @title Get Ordination
#' @description Ordinate phyloseq data and merge it with sample metadata
#' @param x \code{\link{phyloseq-class}} object or a data matrix 
#' (features x samples; eg. HITChip taxa vs. samples)
#' @param method Ordination method, see phyloseq::plot_ordination
#' @param distance Ordination distance, see phyloseq::plot_ordination
#' @return data.frame with ordination coordinates and metadata
#' @examples
#' data(dietswap)
#' fc <- get_ordination(dietswap)
#' @seealso phyloseq::plot_ordination
#' @export
#' @details This is a wrapper for phyloseq ordination functions, providing
#' smooth access to ordinated data.frame with full info on the projection
#' and metadata necessary for further visualizations.
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get_ordination <- function(x, method="NMDS", distance="bray") {
    
    x.ord <- ordinate(x, method, distance)
    
    # Pick the projected data (first two columns + metadata)
    proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
    
    # Rename the projection axes
    names(proj)[1:2] <- paste("Comp", 1:2, sep=".")
    
    proj
    
}

