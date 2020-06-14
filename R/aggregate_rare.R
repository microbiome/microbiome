#' @title Aggregate Rare Groups
#' @description Combining rare taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @inheritParams core
#' @return \code{\link{phyloseq-class}} object
#' @examples
#' data(dietswap)
#' s <- aggregate_rare(dietswap, level = 'Phylum',
#'     detection = 0.1/100, prevalence = 5/100)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_rare <- function (x, level, detection, prevalence, include.lowest=FALSE, ...) {

    x <- aggregate_taxa(x, level)

    rare <- rare_members(x, detection, prevalence, include.lowest)
    
    tax <- tax_table(x)

    inds <- which(rownames(tax) %in% rare)

    tax[inds, level] <- "Other"

    tax_table(x) <- tax

    tt <- tax_table(x)[, level]

    tax_table(x) <- tax_table(tt)

    aggregate_taxa(x, level)

}


