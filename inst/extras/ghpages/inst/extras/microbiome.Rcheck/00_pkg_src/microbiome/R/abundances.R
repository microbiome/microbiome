#' @title Abundance matrix from phyloseq object
#' @description Retrieves the taxon abundance table from
#'    \code{\link{phyloseq-class}} object and ensures it is returned as
#'    taxa x samples matrix.
#' @inheritParams transform
#' @return Abundance matrix.
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @aliases ab
#' @examples
#'   data(dietswap)
#'   abundances(dietswap)
#'   abundances(dietswap, transform = "identity")
#'   abundances(dietswap, transform = "compositional")
#'   abundances(dietswap, transform = "clr")
#'   abundances(dietswap, transform = "Z")
#'   abundances(dietswap, transform = "log10")
#' @keywords utilities
abundances <- function (x, transform = "identity") {

  # Was taxa_abundances

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

  otu

}