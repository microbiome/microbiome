#' @title aggregate_taxa
#' @description Aggregate phyloseq data into higher phylogenetic level
#' when the phylogenetic tree is missing. If the tree is available, 
#' uses the \code{\link{tax_glom}} function from the \pkg{phyloseq} package
#' @param pseq \code{\link{phyloseq-class}} object
#' @param level Aggregation level (from \code{colnames(tax_table(pseq))})
#' @return Aggregated phyloseq object
#' @examples # aggregate_taxa(pseq, "L1")
#' @importFrom phyloseq tax_glom
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function (pseq, level) {

  # Agglomerate taxa
  tg <- tax_glom(pseq, level) 

  if (is.null(pseq@phy_tree)) {

    # If the tree is not available, need to rename the abundance table 
    # rows such that they correspond to the higher level aggregate variable
    # Pick the agglomerated data
    rownames(tg@otu_table) <- as.character(as.data.frame(tax_table(tg))[[level]]) 

  }

  tg

}