#' @title Filter Prevalent
#' @description Filter the phyloseq object to include only prevalent taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection.threshold Detection threshold for absence/presence.
#' @param prevalence.threshold Prevalence threshold 
#' @return Filtered phyloseq object including only prevalent taxa
#' @references 
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#'   data("peerj32")
#'   filter_prevalent(peerj32$physeq, 200, 0.2)
filter_prevalent <- function (x, detection.threshold, prevalence.threshold) {
  taxa <- core(x, detection.threshold, prevalence.threshold)
  prune_taxa(taxa, x)
}



