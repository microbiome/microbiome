#' @title Standard data transformations for phyloseq objects
#' @description Provides phyloseq transformations with log10(x), log10(1+x), z transformation, and relative abundance.
#' @param x \code{\link{phyloseq-class}} object
#' @param transformation Transformation to apply: 'relative.abundance', 'Z', 'log10', 'hellinger', 'identity', or any method from the vegan::decostand function.
#' @param target Apply the transformation for 'sample' or 'OTU'. Does not affect the log transformation.
#' @return Transformed \code{\link{phyloseq}} object
#' @details The relative abunance are returned as percentages in [0, 100]. The Hellinger transformation is square root of the relative abundance but instead given at the scale [0,1].
#' @export
#' @examples
#' \dontrun{
#'
#'   # OTU relative abundances
#'   xt <- transform_phyloseq(x, "relative.abundance", "OTU")
#' 
#'   # Z-transformation for OTUs
#'   xt <- transform_phyloseq(x, "Z", "OTU")
#'
#'   # Z-transformation for samples
#'   xt <- transform_phyloseq(x, "Z", "sample")
#'
#'   # Log10 transformation (log(1+x) if the data contains zeroes)
#'   xt <- transform_phyloseq(x, "log10")
#'
#' }
#' @keywords utilities
transform_phyloseq <- function (x, transformation = "relative.abundance", target = "OTU") {

  y <- NULL

  if (!all(sample(round(prod(dim(otu_table(x)))/10) ))%%1 == 0) {
    warning("The OTU abundances are not integers. Check that the OTU input data is given as original counts to avoid transformation errors!")
  }

  if (transformation == "relative.abundance") {
    if (target == "OTU") {
      xt <- transform_sample_counts(x, function (x) {100 * x/sum(x)})
    } else {
      stop(paste("transform_phyloseq not implemented for transformation", transformation, "with target", target))
    }
  } else if (transformation == "Z") {

    # Z transformation for sample or OTU
    xt <- ztransform_phyloseq(x, target)

  } else if (transformation == "log10") {
  
    # Log transformation:
    if (min(otu_table(x)) == 0) {
      warning("OTU table contains zeroes. Using log10(1 + x) transformation.")
      # target does not affect the log transformation 
      xt <- transform_sample_counts(x, function(x) log10(1 + x))      
    } else {
      xt <- transform_sample_counts(x, function(x) log10(x))      
    }

  } else if (transformation == "identity") {

    # No transformation
    xt <- x
    
  } else {
  
    if (target == "OTU") {
    
      xt <- x
      a <- try(xx <- decostand(taxa_abundances(xt), method = transformation, MARGIN = 2))
            
      if (class(a) == "try-error") {
        xt <- NULL
        stop(paste("Transformation", transformation, "not defined."))
      }

      if (!taxa_are_rows(xt)) {xx <- t(xx)}
      otu_table(xt)@.Data <- xx



     } else {
    
      stop(paste("Transformation", transformation, "not defined for", target))
      
    }    
  }

  xt
}
