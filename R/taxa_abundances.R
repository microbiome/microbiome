#' @title Abundance matrix from phyloseq object
#' @description Retrieves the taxon abundance table from \code{\link{phyloseq-class}} object and ensures it is returned as taxa x samples matrix.
#' @param x \code{\link{phyloseq}} object
#' @return Abundance matrix.
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @examples
#'   data(dietswap)
#'   taxa_abundances(dietswap)
#' @keywords utilities
taxa_abundances <- function (x) {

  # Pick OTU matrix
  otu <- get_taxa(x)

  if (ntaxa(x) == 1) {
    otu = matrix(otu, nrow = 1)
    rownames(otu) = taxa_names(x)
    colnames(otu) = sample_names(x)
  }

  # Ensure that taxa are on the rows
  if (!taxa_are_rows(x) && ntaxa > 1) {
    otu <- t(otu)
  }

  otu

}