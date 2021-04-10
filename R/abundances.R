#' @title Abundance Matrix from Phyloseq
#' @description Retrieves the taxon abundance table from
#' phyloseq-class object and ensures it is systematically returned as
#' taxa x samples matrix.
#' @inheritParams transform
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
abundances <- function(x, transform="identity") {

    # Pick the OTU data
    if (any(c("phyloseq", "otu_table") %in% is(x))) {

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





