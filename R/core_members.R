#' @title Core Taxa
#' @description Determine members of the core microbiota with given abundance
#' and prevalences
#' @inheritParams core
#' @return Vector of core members
#' @details For phyloseq object, lists taxa that are more prevalent with the
#' given detection threshold. For matrix, lists columns that satisfy
#' these criteria.
#' @examples
#' data(dietswap)
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (dietswap) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' a <- core_members(dietswap, 1, 50/100)
#' @aliases prevalent_taxa
#' @export
#' @references 
#' A Salonen et al. The adult intestinal core microbiota is determined by 
#' analysis depth and health status. Clinical Microbiology and Infection 
#' 18(S4):16 20, 2012. 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_members <- function(x, detection=1/100, prevalence=50/100,
    include.lowest=FALSE) {

    if ((prevalence < 0) | (prevalence > 1)) {
        stop("The prevalence argument should be in [0, 1].")
    }

    if (is_compositional(x)) {
        if ((detection < 0) | (detection > 1)) {
            stop("The detection arguments should be in [0, 1] for 
            compositional data.")
        }      
    }

    # Pick taxa x samples matrix
    x <- abundances(x)
    
    if (include.lowest) {
        taxa <- names(which(prevalence(x, detection,
            include.lowest=include.lowest) >= prevalence))
    } else {
        taxa <- names(which(prevalence(x, detection,
        include.lowest=include.lowest) > prevalence))
    }
    
    taxa
    
}


