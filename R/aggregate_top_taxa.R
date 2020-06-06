#' @title Aggregate Top Taxa
#' @description Summarize phyloseq: combine other than the most abundant taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param top Keep the top-n taxa, and merge the rest under the category
#' 'Other'. Instead of top-n numeric this can also be a character vector
#' listing the groups to combine.
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @return \code{\link{phyloseq-class}} object
#' @examples
#' data(dietswap)
#' s <- aggregate_top_taxa(dietswap, top = 3, 'Phylum')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_top_taxa <- function (x, top, level) {

   .Deprecated("aggregate_rare", "The microbiome::aggregate_top_taxa function is deprecated.")

        x <- aggregate_taxa(x, level)

        tops <- top_taxa(x, top)
        tax <- tax_table(x)

        inds <- which(!rownames(tax) %in% tops)

        tax[inds, level] <- "Other"

        tax_table(x) <- tax

        tt <- tax_table(x)[, level]
        tax_table(x) <- tax_table(tt)

        aggregate_taxa(x, level)

}


