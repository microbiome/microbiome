

#' Description: Calculate ANOVA test (BH correction) for multi-group comparison 
#' Arguments:
#'   @param dat data matrix (features x samples; eg. HITChip taxa vs. samples)
#'   @param group Vector with specifying the groups
#'   @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). If NULL, no correction will be performed.
#'   @param sort sort the results
#'
#' Returns:
#'   @return (Corrected) p-values for multi-group comparison.
#'
#' @examples data(peerj32); 
#'          pval <- check.anova(t(peerj32$microbes), 
#'                             peerj32$meta$time)
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check.anova <- function(dat, group, p.adjust.method = "BH", sort = FALSE) {
    
    if (is.null(p.adjust.method)) {
        p.adjust.method <- "none"
    }
    
    pvals <- c()
    for (tax in rownames(dat)) {
        pvals[[tax]] <- anova(lm(dat[tax, ] ~ group))["group", "Pr(>F)"]
    }
    
    pvals <- p.adjust(pvals, method = p.adjust.method)
    
    if (sort) {
        pvals <- sort(pvals)
    }
    
    pvals
}
