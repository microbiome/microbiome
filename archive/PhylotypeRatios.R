#' @title Phylotype ratios
#' @description Calculate phylotype ratios (eg. Bacteroides vs.
#'          Prevotella etc.) for a given phylotypes vs. samples matrix
#' @param dat phylotypes vs. samples data matrix in log10 scale
#' @return phylotype pairs x samples matrix indicating the ratio 
#'                 (in log10 domain) between each unique pair
#' @export 
#' @examples 
#'   data(peerj32)
#'   ratios <- PhylotypeRatios(peerj32$microbes)
#' @references See citation('microbiome')
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
PhylotypeRatios <- function(dat) {
    
    phylogroups <- rownames(dat)
    Nratios <- (length(phylogroups)^2 - length(phylogroups))/2
    Nsamples <- ncol(dat)
    ratios <- list()
    for (i in 1:(length(phylogroups) - 1)) {
        for (j in (i + 1):length(phylogroups)) {
            pt1 <- phylogroups[[i]]
            pt2 <- phylogroups[[j]]
            ratios[[paste(pt1, pt2, sep = "-")]] <- dat[pt1, ] - dat[pt2, ]
        }
    }
    ratios <- do.call(cbind, ratios)
    
    t(ratios)
}

