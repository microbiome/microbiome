#' @title Aggregate Taxa
#' @description Summarize phyloseq data into a higher phylogenetic level.
#' @details This provides a convenient way to aggregate phyloseq OTUs
#' (or other taxa) when the phylogenetic tree is missing. Calculates the
#' sum of OTU abundances over all OTUs that map to the same higher-level
#' group. Removes ambiguous levels from the taxonomy table. Returns a
#' phyloseq object with the summarized abundances.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Summarization level (from \code{rank_names(pseq)})
#' @param rm.na Remove NA taxa
#' @param verbose verbose
#' @return Summarized phyloseq object
#' @examples
#' data(dietswap)
#' s <- aggregate_taxa(dietswap, 'Phylum')
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
aggregate_taxa <- function(x, level, rm.na = TRUE, verbose = FALSE) {

<<<<<<< HEAD
    # FIXME: this function contains quick hacks to circumvent
    # missing tax_table and sample_data. Those would be better handled
    # in the original reading functions.

    mypseq <- check_phyloseq(x)

    if (!is.null(mypseq@phy_tree)) {
        
        if (!is.null(top)) {
            warning("The top parameter to be implemented when phy_tree 
                is available.")
        }
        
        # Agglomerate taxa
        mypseq2 <- tax_glom(mypseq, level)
        mypseq2@phy_tree <- NULL # Remove tree 
        a <- abundances(mypseq2)
        nams <- as.character(tax_table(mypseq2)[, level])
        rownames(a) <- nams
        tt <- tax_table(mypseq2)[, seq_len(match(level,
        colnames(tax_table(mypseq2))))]
        rownames(tt) <- nams

        mypseq2 <- phyloseq(otu_table(a, taxa_are_rows=TRUE), 
                sample_data(mypseq2), 
                tax_table(tt))

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

        # Split the OTUs in tax_table by the given taxonomic level
        if (all(taxa(mypseq) == rownames(tax_table(mypseq)))) {
            level.from <- taxa(mypseq)
        } else {     
            v <- apply(tt, 2, function(x) {
                mean(taxa(mypseq) %in% unique(x))})            
            if (max(v) > 0) {
                current.level <- names(which.max(v))
                level.from <-  as.vector(tt[, 4])
            } else {
                stop("Taxa not found in aggregate_taxa. Halting.")
            }
        }

        level.to <- as.vector(tax_table(mypseq)[, level])        
        otus <- split(level.from, level.to)

        ab <- matrix(NA, nrow=length(otus), ncol=nsamples(mypseq))
        rownames(ab) <- names(otus)
        colnames(ab) <- sample_names(mypseq)

        d <- abundances(mypseq)

        for (nam in names(otus)) {
            taxa <- otus[[nam]]
            ab[nam, ] <- colSums(matrix(d[taxa, ], ncol=nsamples(mypseq)),
            na.rm = TRUE)
        }

        # Create phyloseq object
        OTU <- otu_table(ab, taxa_are_rows=TRUE)
        mypseq2 <- phyloseq(OTU)

        # Remove ambiguous levels
        ## First remove NA entries from the target level
        keep <- !is.na(tax_table(mypseq)[, level])
        tax_table(mypseq) <- tax_table(mypseq)[keep,]
        keep <- colnames(
        tax_table(mypseq))[
        which(
            vapply(seq(ncol(tax_table(mypseq))),
            function(k)
            sum(
            vapply(split(as.character(tax_table(mypseq)[, k]),
            as.character(tax_table(mypseq)[, level])), function(x) {
            length(unique(x))
        }, 1) > 1), 1) == 0)]
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
||||||| merged common ancestors
    # FIXME: this function contains quick hacks to circumvent
    # missing tax_table and sample_data. Those would be better handled
    # in the original reading functions.
    mypseq <- check_phyloseq(x)
    
    if (!is.null(mypseq@phy_tree)) {
        
        if (!is.null(top)) {
            warning("The top parameter to be implemented when phy_tree 
                is available.")
        }
        
        # Agglomerate taxa
        mypseq2 <- tax_glom(mypseq, level)
        mypseq2@phy_tree <- NULL # Remove tree 
        a <- abundances(mypseq2)
        nams <- as.character(tax_table(mypseq2)[, level])
        rownames(a) <- nams
        tt <- tax_table(mypseq2)[, seq_len(match(level,
        colnames(tax_table(mypseq2))))]
        rownames(tt) <- nams

    mypseq2 <- phyloseq(otu_table(a, taxa_are_rows=TRUE), 
                sample_data(mypseq2), 
                tax_table(tt))

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
        v <- apply(tt, 2, function(x) {mean(taxa(mypseq) %in% unique(x))})
        if (max(v) > 0) {
            current.level <- names(which.max(v))
        } else {
            stop("The taxa are not found in tax_table in aggregate_taxa") 
        }
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
            ab[nam, ] <- colSums(matrix(d[taxa, ], ncol=nsamples(mypseq)),
            na.rm = TRUE)
        }


        # Create phyloseq object
        OTU <- otu_table(ab, taxa_are_rows=TRUE)
        mypseq2 <- phyloseq(OTU)

        # Remove ambiguous levels
        ## First remove NA entries from the target level
        keep <- !is.na(tax_table(mypseq)[, level])
        tax_table(mypseq) <- tax_table(mypseq)[keep,]
        keep <- colnames(
        tax_table(mypseq))[
        which(
            vapply(seq(ncol(tax_table(mypseq))),
            function(k)
            sum(
            vapply(split(as.character(tax_table(mypseq)[, k]),
            as.character(tax_table(mypseq)[, level])), function(x) {
            length(unique(x))
        }, 1) > 1), 1) == 0)]
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
=======
    # x <- GlobalPatterns; level <- "Phylum"; top <- NULL; fill_na_taxa <- FALSE;

    # Check if the object is already at the given level

    check1 <- length(which(all(tax_table(x)[, level] == tax_table(x)[, ncol(tax_table(x))]))) > 0
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
>>>>>>> db93d7a2dda19314c25956730051e0f5b79b68e7
    }

    if (verbose) {print("Mark the potentially ambiguous taxa")}
    # Some genera for instance belong to multiple Phyla and perhaps these are different
    # genera. For instance there is genus Clostridium in Tenericutes and Firmicutes.
    # (GlobalPatterns data set) and even more families.
    
    tt <- tax_table(x)
    if (verbose) {print("-- split")}
    otus <- split(rownames(tt), as.character(tt[, "unique"]))

    ab <- matrix(NA, nrow=length(otus), ncol=nsamples(x))
    rownames(ab) <- names(otus)
    colnames(ab) <- sample_names(x)

    if (verbose) {print("-- sum")}
    d <- abundances(x)
    ab <- t(sapply(otus, function (taxa) {colSums(matrix(d[taxa, ], ncol=nsamples(x)), na.rm = TRUE)}))
    colnames(ab) <- colnames(d)
    rownames(ab) <- names(otus)

    if (verbose) {print("Create phyloseq object")}
    OTU <- otu_table(ab, taxa_are_rows=TRUE)
    x2 <- phyloseq(OTU)

    if (verbose) {print("Remove ambiguous levels")}
    ## First remove NA entries from the target level
    inds3 <- match(level, colnames(tt@.Data))
    inds4 <- match("unique", colnames(tt@.Data))    
    taxtab <- tt@.Data[which(!is.na(tt@.Data[, level])), c(seq_len(inds3), inds4)]

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




