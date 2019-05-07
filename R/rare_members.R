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
#' a <- rare_members(dietswap, detection=100/100, prevalence=20/100)
#' @seealso core_members, noncore_members
#' @export
#' @references
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

#' @inheritParams core_members
#' @return Filtered phyloseq object including only rare taxa
#' @references 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' data(peerj32)
#' pseq <- noncore_members(peerj32$phyloseq, 200, 20/100)
noncore_members <- function(x, detection, prevalence, include.lowest=FALSE) {
    
    # TODO: add optional renormalization such that the abundances
    # would sum up to 1 ?
    
    # Core taxa
    cm <- core_members(x, detection, prevalence,
                    include.lowest=include.lowest)
    
    # Non-core taxa as complement of core taxa
    rt <- setdiff(taxa(x), cm)
    
    # Pick the subset
    ret <- NULL
    if (length(rt) > 0) {
        ret <- prune_taxa(rt, x)
    } else {
        warning("No rare taxa with the given thresholds. Returning NULL.")
    }
    
    ret
    
}




