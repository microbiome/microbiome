#' @title Variable Microbiota
#' @description Filter the phyloseq object to include only variable taxa.
#' @inheritParams core
#' @return Filtered phyloseq object including only variable taxa.
#' @details Variable taxa are members of the microbiota that are not rare, and
#' not members of the core microbiota. Such taxa reach notable abundance in a
#' notable fraction of the population. The variable members are calculated by
#' using the indicated prevalence threshold to define rare taxa; and its
#' complement 1-prevalence to define core taxa. The remaining taxa are
#' considered variable. The same detection threshold is used in all cases.
#' @references
#'   To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_members, core, rare
#' @export
#' @examples
#'   data(dietswap)
#'   # Detection threshold 0 (strictly greater by default);
#'   # Prevalence threshold 50 percent (strictly greater by default)
#'   pseq.rel <- transform(dietswap, 'compositional')
#'   pseq <- variable_members(pseq.rel, 0.5/100, 20/100)
#'
variable_members <- function(x, detection, prevalence, include.lowest = FALSE) {
    
    xorig <- x
    
    # Rare taxa (e.g. 20 percent prevalence)
    rm <- taxa(rare(x, detection, prevalence, include.lowest = include.lowest))
    
    # Core taxa (e.g. 80 percent prevalence - complement)
    cm <- core_members(x, detection, 1 - prevalence, include.lowest = include.lowest)
    
    taxa <- setdiff(taxa(x), c(rm, cm))
    if (length(taxa) == 0) {
        taxa <- NULL
    }
    # prune_taxa(taxa, xorig)
    
    if (is.null(taxa)) {
        warning("No variable taxa with the given thresholds. Returning NULL.")
    }
    
    taxa
}






