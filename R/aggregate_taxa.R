#' @title Aggregate Taxa
#' @description Summarize phyloseq data into a higher phylogenetic level.
#' @details This provides a convenient way to aggregate phyloseq OTUs
#' (or other taxa) when the phylogenetic tree is missing. Calculates the
#' sum of OTU abundances over all OTUs that map to the same higher-level
#' group. Removes ambiguous levels from the taxonomy table. Returns a
#' phyloseq object with the summarized abundances.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @param top Keep the top-n taxa, and merge the rest under the category
#' 'Other'. Instead of top-n numeric this can also be a character vector
#' listing the groups to combine.
#' @return Summarized phyloseq object
#' @examples
#'   data(dietswap)
#'   s <- aggregate_taxa(dietswap, 'Phylum')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function(x, level, top = NULL) {

    # FIXME: this function contains quick hacks to circumvent
    # missing tax_table and sample_data. Those would be better handled
    # in the original reading functions.

    mypseq <- x
    
    if (!is.null(mypseq@phy_tree)) {
        
        if (!is.null(top)) {
            warning("The top parameter to be implemented when phy_tree 
                is available.")
        }
        
        # Agglomerate taxa
        mypseq2 <- tax_glom(mypseq, level)
        
    } else {
        
        tt <- tax_table(mypseq)
        if (!is.null(top)) {

            # Merge the remaining taxa into a single group named "Other"
            if (is.numeric(top)) {
                top <- top_taxa(aggregate_taxa(mypseq, level), top)
            }
            
            tt[which(!tt[, level] %in% top), level] <- "Other"
            tax_table(mypseq) <- tt
        }
        
        # Split the OTUs in tax_table by the given taxonomic level otus <-
        # split(rownames(tax_table(mypseq)), tax_table(mypseq)[, level])
        current.level <- names(which.max(apply(tt, 2, function(x) {
            mean(taxa(mypseq) %in% unique(x))
        })))
    if (length(current.level) == 0) {
            current.level <- "unique"
        tax_table(mypseq) <- tax_table(cbind(tax_table(mypseq),
	    unique = rownames(tax_table(mypseq))))
        }

        otus <- map_levels(data=mypseq, to=current.level, from=level)
        
        ab <- matrix(NA, nrow=length(otus), ncol=nsamples(mypseq))
        rownames(ab) <- names(otus)
        colnames(ab) <- sample_names(mypseq)
        
        d <- abundances(mypseq)
        
        for (nam in names(otus)) {
            taxa <- otus[[nam]]
            ab[nam, ] <- colSums(matrix(d[taxa, ], ncol=nsamples(mypseq)))
        }
        
        # Create phyloseq object
        OTU <- otu_table(ab, taxa_are_rows=TRUE)
        mypseq2 <- phyloseq(OTU)
        
        # Remove ambiguous levels
    ## First remove NA entries from the target level    
    tax_table(mypseq) <- tax_table(mypseq)[!is.na(tax_table(mypseq)[, level]),]
        keep <- colnames(
        tax_table(mypseq))[which(vapply(1:ncol(tax_table(mypseq)),
            function(k)
            sum(vapply(split(as.character(tax_table(mypseq)[, k]),
            as.character(tax_table(mypseq)[, level])), function(x) {
            length(unique(x))
        }) > 1)) == 0)]
        tax <- unique(tax_table(mypseq)[, keep])

        # Rename the lowest level
        tax <- as.data.frame(tax)
        rownames(tax) <- tax[, level]
        tax$OTU <- rownames(tax)
        tax <- as.matrix(tax)
        
        # Convert to taxonomy table
        TAX <- tax_table(tax)
        
        # Combine OTU and Taxon matrix into Phyloseq object
        mypseq2 <- merge_phyloseq(mypseq2, TAX)

        # Add the metadata as is
        if (!is.null(mypseq@sam_data)) {
            mypseq2 <- merge_phyloseq(mypseq2, sample_data(mypseq))
        }      
    }
    
    mypseq2
    
}

