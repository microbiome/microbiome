#' @title Diversity within a Sample Group
#' @description Quantify microbiota divergence (heterogeneity) within a
#' given sample set.
#'
#' @details
#' Microbiota divergence (heterogeneity / spread) within a given sample
#' set can be quantified by the average sample dissimilarity or beta
#' diversity. Taking average over
#' all pairwise dissimilarities is sensitive to sample size and heavily biased
#' as the similarity values are not independent. To reduce this bias, the
#' dissimilarity of each sample against the group mean is calculated. This
#' generates one value per sample. These can be compared between groups in
#' order to compare differences in group homogeneity. 
#'
#' Note that this measure is still affected by sample size.
#' Subsampling or bootstrapping can be applied to equalize sample sizes
#' between comparisons.
#' 
#' The spearman mode is a simple indicator that returns
#' average spearman correlation between samples of the input data and
#' the overall group-wise average. The inverse of this measure
#' (ie rho instead of 1-rho as in here) was used in Salonen et al. (2014)
#' to quantify group homogeneity.
#' 
#' @param x phyloseq object 
#' @param method dissimilarity method: any method available via
#' stats::cor or phyloseq::distance function. Note that some methods
#' ("jsd" and 'unifrac' for instance) do not work with the group divergence.
#' @param coreset phyloseq object; the samples to be used to define the centroid
#' @return Vector with dissimilarities; one for each sample, quantifying the
#' dissimilarity of the sample from the group-level mean.
#' @export
#' @examples
#' # Assess beta diversity among the African samples
#' # in a diet swap study (see \code{help(dietswap)} for references)
#' data(dietswap)
#' b <- divergence(subset_samples(dietswap, nationality == 'AFR'))
#' @references
#'
#' The inter- and intra-individual homogeneity measures used in
#' Salonen et al. ISME J. 8:2218-30, 2014 were obtained as
#' 1 - beta where beta is the group diversity as quantified by the
#' spearman method.
#' 
#' To cite this R package, see citation('microbiome')
#' 
#' @seealso the vegdist function from the \pkg{vegan} package provides many
#' standard beta diversity measures
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
divergence <- function(x, method="bray", coreset = NULL) {

    if (is.null(coreset)) {
        coreset <- x
    }

    # Abundance matrix (taxa x samples)
    x <- abundances(x)
    coreset <- abundances(coreset)    

    if (method %in% c("spearman", "pearson", "kendall")) {
        b <- correlation_divergence(x, method = method, coreset = coreset)
    } else {
        b <- beta.mean(x, method = method, coreset = coreset)
    }

    # Add sample names
    names(b) <- colnames(x)

    b
    
}

correlation_divergence <- function(x, method="spearman", coreset) {

    # Correlations calculated against the mean of the sample set
    m <- matrix(rowMeans(coreset))
    cors <- as.vector(cor(x, m,
        method=method, use="pairwise.complete.obs"))
    
    1 - cors
    
}


beta.mean <- function(x, method="bray", coreset) {
    
    # Divergence calculated against the mean of the sample set
    b <- c()
    m <- rowMeans(coreset)

    for (i in seq_len(ncol(x))) {
        xx <- rbind(x[, i], m)
        xxx <- distance(otu_table(t(xx), taxa_are_rows = TRUE), method=method)
        b[[i]] <- as.matrix(xxx)[1, 2]
    }
    
    b
    
}





beta.pairs <- function(x, method="bray", n = ncol(x)) {
    
    # Divergence calculated against the mean of the sample set
    b <- c()

    # All index pairs
    # pairs <- combn(1:ncol(x), 2)

    # Pick n pairs without replacement
    # inds <- sample(ncol(pairs), n)

    # Distance between each sample and a random pair
    for (i in seq_len(ncol(x))) {
        i2 <- sample(setdiff(seq_len(ncol(x)), i), 1)
        xx <- t(x[, c(i, i2)])
        b[[i]] <- as.vector(vegdist(xx, method=method)) 
    }
    
    b
    
}


