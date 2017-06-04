#' @title Rare Taxa
#' @description Determine members of the rare microbiota with given abundance
#' and prevalence threshold.
#' @inheritParams core
#' @return Vector of rare taxa
#' @details For phyloseq object, lists taxa that are less prevalent than the
#' given prevalence threshold. Optionally, never exceeds the given abundance
#' threshold (by default, all abundanecs accepted). For matrix, lists
#' columns that satisfy these criteria.
#' @examples
#' data(dietswap)
#' # Detection threshold: the taxa never exceed the given detection threshold
#' # Prevalence threshold 20 percent (strictly greater by default)
#' a <- rare_members(dietswap, detection = 100/100, prevalence = 20/100)
#' @seealso core_members, noncore_members
#' @export
#' @references
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rare_members <- function(x, detection = 100/100, prevalence = 20/100) {
    
    # Pick taxa x samples matrix
    x <- abundances(x)
    
    taxa1 <- names(which(
    
        # Detected in at most given prevalence
        prevalence(x, detection = 0,
                include.lowest = FALSE) <= prevalence
        ))
        
    # .. and never exceeds detection threshold
    taxa2 <- names(which(apply(x <= detection, 1, all)))

    taxa <- intersect(taxa1, taxa2)

    taxa
    
}


