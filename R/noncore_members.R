#' @title Select Non-Core Taxa
#' @description Filter the phyloseq object to include only taxa that do not
#' belong to the core set defined by the detection and prevalence thresholds.
#' @inheritParams core_members
#' @return Filtered phyloseq object including only rare taxa
#' @references 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' #data(peerj32)
#' #pseq <- noncore_members(peerj32$phyloseq, 200, 20/100)
noncore_members <- function(x, detection, prevalence, include.lowest=FALSE) {

    .Deprecated("rare_members", "The microbiome::noncore_members has been deprecated.")

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




