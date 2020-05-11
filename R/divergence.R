#' @title Divergence within a Sample Group
#' @description Quantify microbiota divergence (heterogeneity) within a
#' given sample set with respect to a reference.
#'
#' @details
#' Microbiota divergence (heterogeneity / spread) within a given sample
#' set can be quantified by the average sample dissimilarity or beta
#' diversity with respect to a given reference sample.
#'
#' This measure is sensitive to sample size.
#' Subsampling or bootstrapping can be applied to equalize sample sizes
#' between comparisons.
#' 
#' @param x phyloseq object or a vector
#' @param y Reference sample. A vector. 
#' @param method dissimilarity method: any method available via
#' phyloseq::distance function. Note that some methods
#' ("jsd" and 'unifrac' for instance) do not work with the group divergence.
#' @return Vector with dissimilarities; one for each sample, quantifying the
#' dissimilarity of the sample from the reference sample.
#' @export
#'
#' @examples
#' # Assess beta diversity among the African samples
#' # in a diet swap study (see \code{help(dietswap)} for references)
#' data(dietswap)
#' pseq <- subset_samples(dietswap, nationality == 'AFR')
#' reference <- apply(abundances(pseq), 1, median)
#' b <- divergence(pseq, reference, method = "bray")
#'
#' @references
#'
#' To cite this R package, see citation('microbiome')
#' 
#' @seealso the vegdist function from the \pkg{vegan} package provides many
#' standard beta diversity measures
#'
#' @author Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
divergence <- function(x, y, method="bray") {

    # Abundance matrix (taxa x samples)
    if (is.phyloseq(x)) {
        x <- abundances(x)
    }

    if (!is.matrix(x)) {
        x <- matrix(x, nrow = length(x))
    }

    y <- as.vector(y)

    # Divergence against the reference
    b <- c()
    for (i in seq_len(ncol(x))) {
        xx <- rbind(x[, i], y)

        xxx <- distance(otu_table(t(xx), taxa_are_rows = TRUE), method=method)

        b[[i]] <- as.matrix(xxx)[1, 2]
    }
    
    # Add sample names
    names(b) <- colnames(x)

    unlist(b)
    
}


    

