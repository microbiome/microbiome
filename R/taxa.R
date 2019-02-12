#' @title Taxa Names
#' @description List the names of taxonomic groups in a phyloseq object.
#' @inheritParams core_members
#' @return A vector with taxon names.
#' @details A handy shortcut for phyloseq::taxa_names, with a potential to add
#' to add some extra tweaks later.
#' @references 
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @examples
#' data(dietswap)
#' x <- taxa(dietswap)
taxa <- function(x) {

    if (taxa_are_rows(x)) {
        rownames(otu_table(x))    
    } else {
        colnames(otu_table(x))    
    }  
}




