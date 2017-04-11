#' @title Abundance matrix from phyloseq object
#' @description Retrieves the taxon abundance table from
#'    \code{\link{phyloseq-class}} object and ensures it is returned as
#'    taxa x samples matrix.
#' @inheritParams transform
#' @return Abundance matrix (OTU x samples).
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @aliases ab, otu
#' @examples
#'   data(dietswap)
#'   a <- abundances(dietswap)
#'   # b <- abundances(dietswap, transform = "compositional")
#' @keywords utilities
abundances <- function (x, transform = "identity") {

  otu <- pick_abundances(x)

  # Apply the indicated transformation
  otu <- transform(otu, transform)

  # Ensure the output is a matrix
  #otu <- pick_abundances(otu)
  otu <- otu@.Data

  otu
}


pick_abundances <- function (x) {

  # Pick the OTU data
  if (is.phyloseq(x)) {
    # Do not apply transformation at this point yet
    otu <- abundances_help(x, transform = "identity")
    
  } else if (is.vector(x)) {
    otu <- as.matrix(x, ncol = 1)    
  } else {
    # If x is not vector or phyloseq object
    # then let us assume it is a taxa x samples count matrix
    otu <- x
  }

  otu
}

abundances_help <- function (x, transform = "identity") {

  # Pick OTU matrix
  otu <- get_taxa(x)

  if (ntaxa(x) == 1) {
    otu <- matrix(otu, nrow = 1)
    rownames(otu) <- taxa(x)
    colnames(otu) <- sample_names(x)
  }

  if (nsamples(x) == 1) {
    otu <- matrix(otu, ncol = 1)
    rownames(otu) <- taxa(x)
    colnames(otu) <- sample_names(x)
  }

  # Ensure that taxa are on the rows
  if (!taxa_are_rows(x) && nrow(otu) > 1 && nsamples(x) > 1) {
    otu <- t(otu)
  }

  if (nrow(otu) == 1) {
    rownames(otu) <- taxa(x)
    colnames(otu) <- sample_names(x)    
  }

  otu@.Data

}




