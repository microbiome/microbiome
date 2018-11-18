#' @title Exclude Samples
#' @description Filter out selected samples from a phyloseq object.
#' @param samples Names of samples to be removed.
#' @param x \code{\link{phyloseq-class}} object
#' @return Filtered phyloseq object.
#' @references
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @export
#' @seealso phyloseq::prune_samples, phyloseq::subset_samples
#' @details This complements the phyloseq function prune_samples by providing
#' a way to exclude given groups from a phyloseq object.
#' @examples
#' data(dietswap)
#' pseq <- remove_samples(c("Sample-182", "Sample-222"), dietswap)
remove_samples <- function(samples = NULL, x) {

    if (is.null(samples)) {return(x)}
    
    i <- samples %in% sample_names(x)

    if (!all(i)) {
        warning(paste("Of the given sample removal list, ",
        round(100*mean(i)), "% (n=", sum(i), ") match the data. 
        Removing these.", sep = ""))
    }    

    keep <- setdiff(sample_names(x), samples)
    prune_samples(keep, x)
    
}






