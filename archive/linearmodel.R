#' @title Linear model with random subject effects for phyloseq
#' @description Linear model test for paired comparisons for phyloseq objects with random effect subject term.
#' @param x \code{\link{phyloseq-class}} object 
#' @param group Metadata field specifying the groups.
#' @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). For other options, see ?p.adjust
#' @return Corrected p-values. 
#' @examples 
#'   data(peerj32)
#'   pval <- check_lmer(peerj32$phyloseq, "group")
#' @export
#' @importFrom lme4 lmer
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_lmer <- function (x, group, p.adjust.method = "BH") {

  # We need taxa x samples matrix
  mydata <- otu_table(x)@.Data
  if (!taxa_are_rows(x)) {mydata <- t(mydata)}

  metadata <- sample_data(x)
  metadata$group <- droplevels(factor(metadata[[group]]))

  # Pvalues
  pv <- c()
  for (tax in rownames(mydata)) {

    dfs <- metadata
    dfs$signal <- mydata[tax, rownames(dfs)]

    # Paired comparison 
    out <- lmer(signal ~ group + (1|subject), data = dfs)
    out0 <- lmer(signal ~ (1|subject), data = dfs)
    comp <- anova(out0, out)
    pv[[tax]] <- comp[["Pr(>Chisq)"]][[2]]
    
  }

  # Adjust ANOVA group-level p-values
  pv <- p.adjust(pv[, "p.anova"], method = p.adjust.method)  
  
  # Order by p 
  pv <- sort(pv)

}



