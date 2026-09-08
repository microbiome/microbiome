#' @title Abundance Matrix from Phyloseq
#' @description Retrieves the taxon abundance table from a
#' phyloseq-class or SummarizedExperiment-derived object (including
#' TreeSummarizedExperiment) and ensures it is systematically returned as
#' taxa x samples matrix.
#' @inheritParams transform
#' @param assay.type Name of the assay to pick when \code{x} is a
#' SummarizedExperiment-derived object. Defaults to \code{NULL}, which
#' selects the \code{counts} assay when present and otherwise the first
#' assay. Ignored for phyloseq objects.
#' @return Abundance matrix (OTU x samples).
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @aliases ab, otu
#' @examples
#' data(dietswap)
#' a <- abundances(dietswap)
#' # b <- abundances(dietswap, transform='compositional')
#' @keywords utilities
abundances <- function(x, transform="identity", assay.type=NULL) {

    # Pick the OTU data
    if (.is_se(x)) {

        # Features are always rows in SummarizedExperiment, so no
        # transposition is needed here.
        otu <- .se_assay(x, assay.type)

    } else if (any(c("phyloseq", "otu_table") %in% is(x))) {

        .deprecate_phyloseq(x)

        # Pick OTU matrix
        otu <- as(otu_table(x), "matrix") # get_taxa(x)

        # Ensure that taxa are on the rows
        if (!taxa_are_rows(x) && ntaxa(x) > 1 && nsamples(x) > 1) {
            otu <- t(otu)
        }

        if (ntaxa(x) == 1) {
            otu <- matrix(otu, nrow=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }

        if (nsamples(x) == 1) {
            otu <- matrix(otu, ncol=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }

    } else if (is.vector(x)) {

        otu <- as.matrix(x, ncol=1)

    } else {

        # If x is not vector or phyloseq object then let us assume it is a
        # taxa x samples
        # count matrix
        otu <- x

    }

    # Apply the indicated transformation
    if (!transform == "identity") {
        otu <- transform(otu, transform)
    }
    otu

}





