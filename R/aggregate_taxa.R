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
#' @param fill_na_taxa If TRUE, the NA entries in tax_table(x) will be
#'    replaced by "Unknown"
#' @return Summarized phyloseq object
#' @examples
#' data(dietswap)
#' s <- aggregate_taxa(dietswap, 'Phylum')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function(x, level, top = NULL, fill_na_taxa = FALSE) {

    # Check if the object is already at the given level	       
    if (all(tax_table(x)[, level] == tax_table(x)[, ncol(tax_table(x))])) {
        return(x)
    }

    # Sanity checks for a phyloseq object. Required with some methods.
    if (!taxa_are_rows(x)) {
        x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)
        taxa_are_rows(x) <- TRUE
    }

    if (fill_na_taxa == TRUE && !is.character(fill_na_taxa)) {
        fill_na_taxa <- "Unknown"
    }

    if (!is.logical(fill_na_taxa) && is.character(fill_na_taxa)) {    
        M <- as.matrix(tax_table(x))

        # Fill in missing entries down to the given level
        # (do not fill lower levels)
        for (i in seq_len(match(level, colnames(M)))) {
            if (any(is.na(M[,i]))) {
                M[which(is.na(M[, i])), i] <- fill_na_taxa
            }
        }

        # Ensure that the filled entries are unique
        inds <- which(M[, level] == fill_na_taxa)
        inds2 <- match(level, colnames(M))
        M[inds, inds2] <- apply(M[inds, seq_len(inds2)], 1,
            function (xx) {paste(xx, collapse = "_")})
        x@tax_table <- tax_table(M)
        
    }
    tt <- tax_table(x)

    if (!is.null(top)) {
        # Merge the remaining taxa into a single group named "Other"
        if (is.numeric(top)) {
            top <- top_taxa(aggregate_taxa(x, level), top)
        }    
        tt[which(!tt[, level] %in% top), level] <- "Other"
        tax_table(x) <- tt
    }

    mytaxa <- taxa(x)
    # Split the OTUs in tax_table by the given taxonomic level 
    v <- apply(tt, 2, function(i) {mean(mytaxa %in% unique(i))})

    if (max(v) > 0) {
        current.level <- names(which.max(v))
    } else {
        stop("The taxa are not found in tax_table in aggregate_taxa") 
    }

    if (length(current.level) == 0) {
        current.level <- "unique"
            tax_table(x) <- tax_table(cbind(tax_table(x),
                unique = rownames(tax_table(x))))
    }

    otus <- map_levels(data=x, to=current.level, from=level)
    ab <- matrix(NA, nrow=length(otus), ncol=nsamples(x))
    rownames(ab) <- names(otus)
    colnames(ab) <- sample_names(x)

    d <- abundances(x)
    for (nam in names(otus)) {
        taxa <- otus[[nam]]
        ab[nam, ] <- colSums(matrix(d[taxa, ], ncol=nsamples(x)), na.rm = TRUE)
    }

    # Create phyloseq object
    OTU <- otu_table(ab, taxa_are_rows=TRUE)
    x2 <- phyloseq(OTU)

    # Remove ambiguous levels
    ## First remove NA entries from the target level
    keep <- !is.na(tax_table(x)[, level])
    tax_table(x) <- tax_table(x)[keep,]

    keep <- colnames(
        tax_table(x))[
            which(
                vapply(seq(ncol(tax_table(x))),
                    function(k)
                        sum(
                vapply(split(as.character(tax_table(x)[, k]),
                    as.character(tax_table(x)[, level])),
                    function(x) {length(unique(x))
                }, 1) > 1), 1) == 0)]

    #tax <- unique(tax_table(x)[, 1:match(keep, colnames(tax_table(x)))])
    tax <- unique(tax_table(x)[, keep])    

    # Rename the lowest level
    tax <- as.data.frame(tax)

    rownames(tax) <- tax[, level]

    tax$OTU <- rownames(tax)
    tax <- as.matrix(tax)
    
    # Convert to taxonomy table
    TAX <- tax_table(tax)

    # Combine OTU and Taxon matrix into Phyloseq object
    x2 <- merge_phyloseq(x2, TAX)

    # Add the metadata as is
    if (!is.null(x@sam_data)) {
        x2 <- merge_phyloseq(x2, sample_data(x))
    }      

    
    x2
    
}

