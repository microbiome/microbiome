#' @title ztransform_phyloseq
#' @description Z transform phyloseq objects
#' @param x \code{\link{phyloseq-class}} object 
#' @param which Specify Z transformation for "sample" or "OTU"
#' @return Z-transformed phyloseq object
#' @examples
#'   pseq <- download_microbiome("peerj32")$physeq
#'   pseqz <- ztransform_phyloseq(pseq, "OTU")
#' @importFrom phyloseq transform_sample_counts
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq otu_table<-
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
ztransform_phyloseq <- function (x, which) {

  taxa_are_rows <- y <- NULL

  # Start with log10 transformation
  x <- transform_phyloseq(x, "log10")

  if (which == "OTU") {

    ddd <- otu_table(x)@.Data
    if (taxa_are_rows(x)) {
      ddd <- t(ddd)
    }

    # Z transform OTUs
    trans <- as.matrix(t(scale(ddd)))

    nullinds <- which(rowMeans(is.na(trans)) == 1)
    if (length(nullinds) > 0 & min(ddd) == 1) {
      warning("Setting undetected OTUs to zero in ztransform_phyloseq")
      # Some OTUs have minimum signal in all samples and scaling gives NA.
      # In these cases just give 0 signal for these OTUs in all samples
      trans[names(which(rowMeans(is.na(trans)) == 1)),] <- 0
    }
    
    xz <- x
    if (!taxa_are_rows(xz)) {trans <- t(trans)}
    otu_table(xz) <- otu_table(trans, taxa_are_rows=taxa_are_rows(xz))  

  } else if (which == "sample") {

    # Z transform samples
    xz <- transform_sample_counts(x, function(x) {(y - mean(y))/sd(y) })

  }
  
  xz

}




