#' @title Beta Diversity
#' @description Quantification of beta diversity.
#' @details
#' 
#' Beta diversity quantifies community similarity between
#' two samples. At the group level, it can be used to quantify group
#' variability (or spread). This function provides a wrapper that is handy
#' for group-level comparisons of beta diversity.
#'
#' When the sample size n=2, this function calculates the beta diversity
#' between the two samples.
#'
#' When the sample size n>2, it calculates beta
#' diversity of each sample against the group mean. This can be
#' compared between groups in order to compare differences in group
#' homogeneity. 
#'
#' Note that this homogeneity measure is affected by sample size.
#' Subsampling or bootstrapping can be applied to equalize sample sizes
#' between comparisons.
#' 
#' More advanced beta correlation measures are available but not implemented
#' here yet.
#' 
#' The anticorrelation mode is a simple educational indicator that returns
#' average spearman correlation between samples of the input data and
#' the overall group-wise average.
#' 
#' @param x phyloseq object 
#' @param mode Beta diversity method
#' @return Vector with beta diversities; one for each sample, quantifying the
#' dissimilarity of the sample from the group-level mean
#' @export
#' @examples 
#'   # Example data
#'   library(microbiome)
#'   data(peerj32)
#'   # Assess beta diversity among the African samples
#'   # in a diet swap study
#'   b <- group_diversity(subset_samples(dietswap, group == "AFR"))
#' @references 
#' The inter- and intra-individual homogeneity in
#' Salonen et al. ISME J. 8:2218-30, 2014 are obtained as
#' 1 - beta where beta is the "anticorrelation" beta diversity.
#' To cite this R package, see citation('microbiome')
#' @seealso the vegdist function from the \pkg{vegan} package provides many
#' standard beta diversity measures
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
group_diversity <- function(x, mode = "anticorrelation") {

  # Abundance matrix (taxa x samples)
  if (is.phyloseq(x)) {
    x <- abundances(x)
  }

  if (ncol(x) == 2) {
    if (mode == "anticorrelation") {
      ret <- 1 - cor(x[,1], x[,2], method = "spearman", use = "pairwise.complete.obs")
    }
  } else if (ncol(x) > 2 & mode == "anticorrelation") {
    ret <- anticorrelation(x, "spearman")
  }

  ret
  
}


anticorrelation <- function(x, method = "spearman") {

  # Correlations calculated against the mean of the sample set
  cors <- as.vector(cor(
       	    x, matrix(rowMeans(x)),
       	    method = method,
	    use = "pairwise.complete.obs"))

  1 - cors    

}
