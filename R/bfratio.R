#' @title Bacteroidetes-Firmicutes Ratio
#' @description Estimates Bacteroidetes-Firmicutes ratio.
#' @inheritParams transform
#' @param b Term used for Bacteroidetes
#' @param f Term used for Firmicutes
#' @return Numeric vector (B/F-ratio)
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @examples
#' data(dietswap)
#' bf <- bfratio(dietswap, "Bacteroidetes", "Firmicutes")
#' @keywords utilities
bfratio <- function (x, b = "Bacteroidetes", f = "Firmicutes") {

    .Deprecated("The bfratio function will be removed in the future releases of microbiome R package. We recommend to calculate this manually.")

    a <- transform(aggregate_taxa(x, level = "Phylum"), "compositional")
    if (all(c(b, f) %in% taxa(a))) {
        b <- abundances(a)[b,]
        f <- abundances(a)[f,]
        return(b/f)
    } else {
        stop(paste("Check that the names for the two phyla are same 
        as in the data (see arguments b and f in the function call)."))
    }
    
}
