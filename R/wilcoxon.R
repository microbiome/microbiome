#' Calculate Wilcoxon test (with BH correction) for two-group comparison 
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
#'   pseq <- download_microbiome("peerj32")$physeq
#'   pval <- check_wilcoxon(pseq, "time")
#'
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_wilcoxon <- function (x, group, p.adjust.method = "BH", sort = FALSE) {

  # Pick the grouping variable from sample metadata
  g <- group
  if (length(g) == 1) {
    g <- sample_data(x)[[g]]
    if (length(unique(g)) == 2) {
      g <- factor(g)
    }
    if (!is.factor(g)) {
      stop(paste("The grouping variable", g, "must be a factor"))
    }    
    g <- droplevels(g)
    if (!length(levels(g)) == 2) {
      stop("check_wilcoxon is valid only for two-group comparisons")
    }
  }

  if (class(x) == "phyloseq") {    
    x <- log10(otu_table(x)@.Data)
  }

  # Calculate Wilcoxon test with BH p-value correction for gender
  pval <- suppressWarnings(apply(x, 1, function (xi) {wilcox.test(xi ~ g)$p.value}))

  # Multiple testing correction
  pval <- p.adjust(pval, method = p.adjust.method)

  # Sort p-values
  if (sort) {
    pval <- sort(pval)
  }

  pval

}

