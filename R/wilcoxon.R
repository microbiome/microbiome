#' Calculate Wilcoxon test (unpaired; BH correction) for the 
#' specified sample groups. 
#'             
#'   @param dat data matrix (features x samples)
#'   @param G1 Sample group 1 (for comparison) 
#'   @param G2 Sample group 2 (for comparison)
#'   @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). If NULL, no correction will be performed.
#'   @param sort sort the results
#'   @param paired paired Wilcoxon test
#'
#'   @return (Corrected) p-values for two-group comparison.
#'
#' @aliases check.wilcoxon
#' @examples 
#'  data(peerj32)
#'  pval <- check_wilcoxon(t(peerj32$microbes), G1 = 1:22, G2 = 23:44)
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_wilcoxon <- function(dat, G1, G2, 
	       	  	   p.adjust.method = "BH", 
			   sort = FALSE, paired = FALSE) {
    
    samples <- colnames(dat)
    levels <- rownames(dat)
    
    M <- matrix(data = NA, length(levels), 1)
    rownames(M) <- levels
    
    for (i in 1:length(levels)) {
        
        lvl <- levels[i]
        l.g1 <- dat[lvl, G1]
        l.g2 <- dat[lvl, G2]
        
        p <- wilcox.test(as.numeric(l.g1), as.numeric(l.g2), 
                paired = paired)$p.value
        
        # message(lvl, ' p-value: ', p, '\n')        
        M[i, 1] <- p
        
    }
    
    ## To Adjust P-values for Multiple Comparisons with 
    ## Benjamini & Hochberg (1995)
    ## ('BH' or its alias 'fdr')
    if (!is.null(p.adjust.method)) {
        cor.p <- p.adjust(M, method = p.adjust.method)
        names(cor.p) <- rownames(M)
        
    } else {
        
        # Skip p-value correction
        cor.p <- as.vector(M)
        names(cor.p) <- rownames(M)
        
    }
    
    # Sort the values
    if (sort) {
        cor.p <- sort(cor.p)
    }
    
    cor.p
    
}

