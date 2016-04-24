#' @title Limma test for phyloseq
#' @description Limma for a single continuous variable vs. OTUs with a phyloseq object.
#' @param x \code{\link{phyloseq-class}} object 
#' @param varname Metadata field specifying the investigated metadata variable.
#' @param transformation Transformation for the original phyloseq values. By default a log10 transformation is done to better approximate the requirements of a linear model.
#' @param p.adj.method p-value correction method for p.adjust function 
#'               (default 'BH'). For other options, see ?p.adjust
#' @return Limma output table.
#' @examples 
#'  data("atlas1006")
#'  tab <- lm_phyloseq(atlas1006, "age")
#' @export
#' @importFrom phyloseq get_variable
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
lm_phyloseq <- function (x, varname, transformation = "log10", p.adj.method = "BH") {

  # Transformation useful for linear models  
  x <- transform_phyloseq(x, transformation)

  # Limma significance analysis
  # Prepare the design matrix which states the variable values for each sample
  design <- cbind(intercept = 1, var = get_variable(x, varname))
  rownames(design) <- rownames(sample_data(x)) 

  # OTU matrix
  otu <- otu_table(x)@.Data
  if (!taxa_are_rows(x)) {otu <- t(otu)}

  # Remove missing vals
  keep = which(rowSums(is.na(design)) == 0)

  if (length(keep) > 0) {
    otu = otu[,keep]
    design = design[keep,]
  } else {
    warning("All samples have missing values in the metadata design matrix. Halting the computation.")
    return(NULL)
  }
  
  # Fit the limma model
  fit <- lmFit(otu, design)
  fit2 <- eBayes(fit)

  # NOTE: results and p-values are given for all groupings in the design matrix
  # Now focus on the second grouping ie. pairwise comparison
  topTable(fit2, coef = 2)

}