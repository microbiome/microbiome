#' @title Standard phyloseq transformations
#' @description Transform phyloseq objects with standard transformations including log10(x), log10(1+x), z, relative abundance.
#' @param x \code{\link{phyloseq-class}} object
#' @param transformation Transformation to apply: 'relative.abundance', 'Z', or 'log10'.
#' @param target Apply the transformation for 'sample' or 'OTU'. Does not affect the log transformation.
#' @return Transformed \code{\link{phyloseq}} object
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
    }
    # target does not affect the log transformation 
    xt <- transform_sample_counts(x, function(x) log10(1 + x))      
  }
  
  xt
}
