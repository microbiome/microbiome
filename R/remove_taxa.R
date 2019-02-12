#' @title Exclude Taxa
#' @description Filter out selected taxa from a phyloseq object.
#' @param taxa Names of taxa to be removed.
#' @param x \code{\link{phyloseq-class}} object
#' @return Filtered phyloseq object.
#' @references
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @seealso phyloseq::prune_taxa, phyloseq::subset_taxa
#' @details This complements the phyloseq function prune_taxa by providing
#' a way to exclude given groups from a phyloseq object.
#' @examples
#' data(dietswap)
#' pseq <- remove_taxa(c("Akkermansia", "Dialister"), dietswap)
remove_taxa <- function(taxa = NULL, x) {

    if (is.null(taxa)) {return(x)}
    i <- taxa %in% taxa(x)

    if (!all(i)) {
        warning(paste("Of the given OTU removal list, ",
        round(100*mean(i)), "% (n=", sum(i), ") match the data. 
        Removing these.", sep = ""))
    }    

    keep <- setdiff(taxa(x), taxa)

    prune_taxa(keep, x)

}



