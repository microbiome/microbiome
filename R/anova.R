#' Calculate ANOVA test (BH correction) for multi-group comparison 
#'
#' @param x \code{\link{phyloseq-class}} object or a data matrix 
#'            (features x samples; eg. HITChip taxa vs. samples)
#' @param group Vector with specifying the groups
#' @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). For other options, see ?p.adjust
#' @param sort sort the results
#'
#' @return Corrected p-values for multi-group comparison.
#'
#' @examples 
#'   data(peerj32)
#'   pval <- check_anova(t(peerj32$microbes), peerj32$meta$time)
#'
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_anova <- function(x, group, p.adjust.method = "BH", sort = FALSE) {

  if (class(x) == "phyloseq") {    
    x <- log10(otu_table(x)@.Data)
  }

    pvals <- c()
    for (tax in rownames(x)) {
        pvals[[tax]] <- anova(lm(x[tax, ] ~ group))["group", "Pr(>F)"]
    }
    
    pvals <- p.adjust(pvals, method = p.adjust.method)
    
    if (sort) {
        pvals <- sort(pvals)
    }
    
    pvals
}



