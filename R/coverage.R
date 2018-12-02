#' @title Coverage Index
#' @description Community coverage index.
#' @details The coverage index gives the number of groups needed to
#' have a given proportion of the ecosystem occupied (by default 0.5 ie 50%). 
#' @param threshold Indicates the fraction of the ecosystem to be occupied by
#' the N most abundant species (N is returned by this function). If the
#' detection argument is a vector, then a data.frame is returned, one
#' column for each detection threshold.
#' @inheritParams alpha
#' @return A vector of coverage indices
#' @export
#' @examples
#' data(dietswap)
#' d <- coverage(dietswap, threshold=0.5)
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso dominance, alpha
#' @keywords utilities
coverage <- function(x, threshold=0.5) {
    
    if (length(threshold) > 1) {
        tab <- vapply(threshold, function(th) {
            coverage(x, threshold=th)
        }, 1)
        colnames(tab) <- as.character(threshold)
        rownames(tab) <- colnames(abundances(x))
        return(tab)
    }
    
    otu <- abundances(x, transform="compositional")
    
    # Number of groups needed to have 50% of the ecosystem occupied
    do <- apply(otu, 2, function(x) {
        min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))
    })
    names(do) <- colnames(otu)
    
    do
    
}
