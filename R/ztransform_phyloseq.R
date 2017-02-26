#' @title Phyloseq Z Transformation
#' @description Z transform phyloseq objects.
#' @details Performs centering (to zero) and scaling (to unit
#'   variance) across samples for each taxa.
#' @param x \code{\link{phyloseq-class}} object 
#' @param which Specify Z transform for "sample" or "OTU"
#' @return Z-transformed phyloseq object
#' @examples \dontrun{
#'   data(peerj32)
#'   pseqz <- ztransform_phyloseq(peerj32$phyloseq)
#' }
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
ztransform_phyloseq <- function (x, which) {

  taxa_are_rows <- y <- NULL

  if (!all(sample(abundances(x), 100)%%1 == 0)) {
    warning("phyloseq object may already have been log transformed - the 
             abundances are not counts - log10 omitted in Z transform. 
	     Perform manually if needed.")
  } else {
    # Start with log10 transform
    x <- transform_phyloseq(x, "log10")
  }
  
  if (which == "OTU") {

    # taxa x samples
    ddd <- abundances(x)

    # Z transform OTUs
    trans <- as.matrix(scale(t(ddd)))

    nullinds <- which(rowMeans(is.na(trans)) == 1)
    if (length(nullinds) > 0 & min(ddd) == 1) {
      warning("Setting undetected OTUs to zero in ztransform_phyloseq")
      # Some OTUs have minimum signal in all samples and scaling gives NA.
      # In these cases just give 0 signal for these OTUs in all samples
      trans[names(which(rowMeans(is.na(trans)) == 1)),] <- 0
    }

    # Use the same matrix format than in original data
    # (taxa x samples or samples x taca)
    xz <- x
    if (taxa_are_rows(x)) { trans = t(trans) }
    otu_table(xz) <- otu_table(trans, taxa_are_rows = taxa_are_rows(x))  

  } else if (which == "sample") {

    # Z transform samples
    xz <- transform_sample_counts(x, function(x) {(y - mean(y))/sd(y) })

  }
  
  xz

}




