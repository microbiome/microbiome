#' @title Summarize Taxa
#' @description Summarize phyloseq data into a higher phylogenetic level.
#' @details This provides a convenient way to aggregate phyloseq OTUs
#'   (or other taxa) when the phylogenetic tree is missing. Calculates the
#'   sum of OTU abundances over all OTUs that map to the same higher-level
#'   group. Removes ambiguous levels from the taxonomy table. Returns a
#'   phyloseq object with the summarized abundances.
#' @param pseq \code{\link{phyloseq-class}} object
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @return Summarized phyloseq object
#' @examples
#'   data(dietswap)
#'   s <- aggregate_taxa(dietswap, "Phylum")
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function (pseq, level) {

  if (!is.null(pseq@phy_tree)) {
  
    # Agglomerate taxa
    pseq2 <- tax_glom(pseq, level)
    
  } else {

    # Split the OTUs in tax_table by the given taxonomic level	       
    #otus <- split(rownames(tax_table(pseq)), tax_table(pseq)[, level])
    current.level <- names(which(apply(tax_table(pseq), 2,
                       function (x) {length(unique(x))}) == ntaxa(pseq)))
    otus <- map_levels(data = pseq, to = current.level, from = level)

    ab <- matrix(NA, nrow = length(otus), ncol = nsamples(pseq))
    rownames(ab) <- names(otus)
    colnames(ab) <- sample_names(pseq)

    d <- abundances(pseq)

    for (nam in names(otus)) {
      taxa <- otus[[nam]]
      ab[nam,] <- colSums(matrix(d[taxa,], ncol = nsamples(pseq)))
    }

    # Create phyloseq object
    OTU <- otu_table(ab, taxa_are_rows = TRUE)
    pseq2 <- phyloseq(OTU)

    # Remove all ambiguous levels 
    keep <- colnames(tax_table(pseq))[which(sapply(1:ncol(tax_table(pseq)),
    	      function (k) sum(sapply(split(as.character(tax_table(pseq)[,k]),
	      as.character(tax_table(pseq)[,level])), function (x)
	      {length(unique(x))}) > 1)) == 0)]
    tax <- unique(tax_table(pseq)[, keep])

    # Rename the lowest level
    tax <- as.data.frame(tax)
    rownames(tax) <- tax[, level]
    tax$OTU <- rownames(tax)
    tax <- as.matrix(tax)

    # Convert to taxonomy table
    TAX <- tax_table(tax)

    # Combine OTU and Taxon matrix into Phyloseq object
    pseq2 <- merge_phyloseq(pseq2, TAX)

    # Add the metadata as is
    pseq2 <- merge_phyloseq(pseq2, sample_data(pseq))

  }

  pseq2

}

