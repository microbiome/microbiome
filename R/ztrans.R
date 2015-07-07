#' Z transform phyloseq objects
#'
#' @param x \code{\link{phyloseq-class}} object 
#' @param which Specify Z transformation for "sample" or "OTU"
#'
#' @return Z-transformed phyloseq object
#'
#' @examples
#'   pseq <- download_microbiome("peerj32")$physeq
#'   pseqz <- ztransform_phyloseq(pseq, "OTU")
#'
#' @importFrom phyloseq transform_sample_counts
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq otu_table<-
#' @export
#'
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
ztransform_phyloseq <- function (x, which) {

  taxa_are_rows <- NULL

  if (which == "OTU") {
  
    # Z transform OTUs
    trans <- as.matrix(t(scale(t(log10(1 + otu_table(x)@.Data)))))

    nullinds <- which(rowMeans(is.na(trans)) == 1)
    if (length(nullinds) > 0) {
      warning("Setting undetected OTUs to zero in ztransform_phyloseq")
      # Some OTUs have minimum signal in all samples and scaling gives NA.
      # In these cases just give 0 signal for these OTUs in all samples
      trans[names(which(rowMeans(is.na(trans)) == 1)),] <- 0
    }
    
    xz <- x
    otu_table(xz) <- otu_table(trans, taxa_are_rows=taxa_are_rows(xz))  

  } else if (which == "sample") {

    # Z transform samples
    xz <- transform_sample_counts(x, function(x) {
       	  	y <- log10(1 + x); (y - mean(y))/sd(y) })

  } else {
    stop("Specify the target for Z transformation: 'sample' or 'OTU'")
  }
  
  xz

}




