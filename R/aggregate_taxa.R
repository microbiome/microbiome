#' @title Aggregate Taxa
#' @description Summarize phyloseq data into a higher phylogenetic level.
#' @details This provides a convenient way to aggregate phyloseq OTUs
#' (or other taxa) when the phylogenetic tree is missing. Calculates the
#' sum of OTU abundances over all OTUs that map to the same higher-level
#' group. Removes ambiguous levels from the taxonomy table. Returns a
#' phyloseq object with the summarized abundances.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @param verbose verbose
#' @return Summarized phyloseq object
#' @examples
#' data(dietswap)
#' s <- aggregate_taxa(dietswap, 'Phylum')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function(x, level, verbose = FALSE) {

    if (!level %in% rank_names(x)) {
        stop("The level argument should be one of the options 
            given by rank_names(x): ",
            paste(rank_names(x), collapse = " / "))
    }

    # Check if the object is already at the given level
    inds <- all(tax_table(x)[, level] == tax_table(x)[, ncol(tax_table(x))])
    inds <- which(inds)
    check1 <- length(inds) > 0
    check2 <- !any(duplicated(tax_table(x)[, level]))
    if (check1 && check2) {
        return(x)
    }

    # Sanity checks for a phyloseq object. Required with some methods.
    if (!taxa_are_rows(x)) {
        x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)
        taxa_are_rows(x) <- TRUE
    }

    fill_na_taxa <- "Unknown"

    if (verbose) {print("Remove taxonomic information below the target level")}
    M <- as.matrix(tax_table(x))
    inds2 <- match(level, colnames(M))    
    M <- M[, seq_len(inds2)]
    M[is.na(M)] <- fill_na_taxa
    # Ensure that the filled entries are unique
    inds <- which(M[, level] == fill_na_taxa)
    M[inds, seq_len(inds2)] <- fill_na_taxa

    unique <- apply(M, 1, function (x) {paste(x, collapse = "_")})
    M <- cbind(M, unique = unique)
    x@tax_table <- tax_table(M)

    if (!nrow(tax_table(x)) == nrow(otu_table(x))) {
        stop("Taxonomic table and OTU table dimensions do not match.")
    }

    if (verbose) {print("Mark the potentially ambiguous taxa")}
    # Some genera for instance belong to multiple Phyla and perhaps these
    # are different
    # genera. For instance there is genus Clostridium in Tenericutes
    # and Firmicutes.
    # (GlobalPatterns data set) and even more families.
    
    tt <- tax_table(x)
    if (verbose) {print("-- split")}
    otus <- split(rownames(tt), as.character(tt[, "unique"]))

    ab <- matrix(NA, nrow=length(otus), ncol=nsamples(x))
    rownames(ab) <- names(otus)
    colnames(ab) <- sample_names(x)

    if (verbose) {print("-- sum")}
    d <- abundances(x)

    ab <- t(vapply(otus, function (taxa) {
        as.numeric(colSums(matrix(d[taxa, ], ncol=nsamples(x)), na.rm = TRUE))
    }, FUN.VALUE = unname(as.numeric(d[1,]))))
    
    colnames(ab) <- colnames(d)
    rownames(ab) <- names(otus)

    if (verbose) {print("Create phyloseq object")}
    OTU <- otu_table(ab, taxa_are_rows=TRUE)
    x2 <- phyloseq(OTU)

    if (verbose) {print("Remove ambiguous levels")}
    ## First remove NA entries from the target level
    inds3 <- match(level, colnames(tt@.Data))
    inds4 <- match("unique", colnames(tt@.Data))    
    taxtab <- tt@.Data[which(!is.na(tt@.Data[, level])),
        c(seq_len(inds3), inds4)]

    if (verbose) {print("-- unique")}
    tax <- unique(taxtab)
    if (verbose) {print("-- Rename the lowest level")}
    tax <- as.data.frame(tax)
    if (verbose) {print("-- rownames")}
    rownames(tax) <- tax[, "unique"]

    if (verbose) {print("-- taxa")}
    tax <- as.matrix(tax)

    if (verbose) {print("Convert to taxonomy table")}
    TAX <- tax_table(tax)

    if (verbose) {print("Combine OTU and Taxon matrix into Phyloseq object")}
    x2 <- merge_phyloseq(x2, TAX)

    # Then keep short names for those taxa where short names are unique
    tt <- tax_table(x2)

    uni <- names(which(table(as.vector(tt[, level])) == 1))
    inds <- which(tt[, level] %in% uni)
    taxa <- tt[inds, level]
    tt[inds, "unique"] <- taxa
    rownames(tt)[inds] <- taxa
    ab <- abundances(x2)
    rnams <- rownames(ab)
    rnams[inds] <- taxa
    rownames(ab) <- rnams

    x2 <- phyloseq(otu_table(ab, taxa_are_rows=TRUE),
            tax_table(tt)) 

    if (verbose) {print("Add the metadata as is")}
    if (!is.null(x@sam_data)) {
        x2 <- phyloseq(otu_table(ab, taxa_are_rows=TRUE),
                    tax_table(tt),
                    sample_data(x))
    }

    x2
    
}




