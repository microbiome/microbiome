#' @title Bacteroidetes-Firmicutes Ratio
#' @description Estimates Bacteroidetes-Firmicutes ratio.
#' @inheritParams transform
#' @return Numeric vector (B/F-ratio)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @examples
#' data(dietswap)
#' bf <- bfratio(dietswap)
#' @keywords utilities
bfratio <- function (x) {
    a <- transform(aggregate_taxa(x, level = "Phylum"), "compositional")
    b <- abundances(a)["Bacteroidetes",]
    f <- abundances(a)["Firmicutes",]  
    b/f
}
