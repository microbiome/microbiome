#' @title Rare Microbiota
#' @description Filter the phyloseq object to include only rare (non-core) taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold for absence/presence
#' (strictly greater by default).
#' @param prevalence Prevalence threshold (in [0, 1];
#' strictly greater by default)
#' @param include.lowest Include the lower boundary of the detection and
#' prevalence cutoffs in core calculation. FALSE by default.
#' @param ... Arguments to pass.
#' @return Filtered phyloseq object including only rare taxa
#' @references
#' Salonen A, Salojarvi J, Lahti L, de Vos WM. The adult intestinal
#' core microbiota is determined by analysis depth and health
#' status. Clinical Microbiology and Infection 18(S4):16-20, 2012
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @seealso core_members
#' @examples
#' data(dietswap)
#' # Detection threshold 0 (strictly greater by default);
#' # Prevalence threshold 50 percent (strictly greater by default)
#' pseq <- rare(dietswap, 0, 50/100)
rare <- function(x, detection, prevalence, include.lowest=FALSE, ...) {
    
    xorig <- x
    
    # TODO: add optional renormalization such that the core member
    # abundances would
    # sum up to 1 ?
    taxa <- core_members(x, detection, prevalence,
        include.lowest=include.lowest)
    
    remove_taxa(taxa, xorig)
    
}






