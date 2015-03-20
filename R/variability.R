#' Stability analysis. Calculates average Pearson '
#' correlation between samples in the input data and picks the lower '
#' triangular matrix to avoid duplicating the correlations. Returns 
#' correlations and stability estimate (average of the correlations). 
#' Can also be used to calculate stability between two data sets. 
#' Then provide two data sets as inputs.
#'
#' @param dat1 data matrix phylotypes vs. samples (in log10 scale)
#' @param dat2 Optional. Second data matrix phylotypes vs. samples. 
#'          Provide this to calculate stability between two (paired) 
#'          data sets.
#' @param method Correlation method (see ?cor)
#'
#' @return List with correlations and astability estimate
#'
#' @export
#' @examples 
#'   data(peerj32)
#'   s <- estimate_variability(t(peerj32$microbes)[, 1:5])
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate_variability <- function(dat1, dat2 = NULL, method = "pearson") {
    
    if (is.null(dat2)) {
        # Within-matrix stability NOTE: earlier this was calculated as
        # the average of upper triangular correlation matrix This is
        # heavily biased since the values are dependent Now replaced
        # by calculating correlations against the mean of the whole
        # sample set cors <- lower.triangle(cor(dat1))
        cors <- as.vector(cor(dat1, matrix(rowMeans(dat1)), method = method))
        names(cors) <- colnames(dat1)
        stab <- list(correlations = cors, stability = mean(cors))
    } else {
        # Between-matrices stability
        cors <- diag(cor(dat1, dat2, method = method))
        stab <- list(correlations = cors, stability = mean(cors))
    }
    
    stab
    
}
